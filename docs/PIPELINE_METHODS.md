# VAPOR — Pipeline Methods & Formulas Reference

Complete technical reference for the VAPOR metagenomics/virome pipeline: architecture,
every processing stage, and **all quantitative formulas** (normalization, abundance,
diversity, thresholds, quality tiers). Kept in sync with the code — every formula below
is transcribed from the actual implementation (file:function noted).

> Scope note: this reflects the pipeline after the **tracks + co-assembly/co-binning**
> restructure (Plans 1–6). Co-assembly group analyses mirror the per-sample formulas
> exactly unless stated otherwise.

---

## 1. Architecture

### 1.1 Tracks (selectable via `config.yaml → tracks`)
- **reads** — reference-based read classification (sylph), no assembly.
- **viral** — assembly-based viral detection → vOTUs (per-sample).
- **prok** — assembly-based prokaryotic MAGs (per-sample).

Resolution + validation is a pure module: `pipeline_config.py::resolve_pipeline_config` /
`validate_pipeline_config`. Derived flags: `INTEGRATION_ENABLED = viral AND prok AND
use_host_defense` (host prediction + defense arms-race only run when both tracks + flag on).

`rule all` is composed from per-track target builders (`Snakefile::_t_foundation/_t_reads/
_t_viral/_t_prok/_t_integration/_t_coassembly/_t_report`). The report/finalize/diversity
aggregators are **track-aware** (conditional inputs) so selecting a subset of tracks
actually prunes the DAG.

### 1.2 Central hub
Per sample, the deduplicated representative set `{sample}_rep_seq.fasta` (MMseqs2, 95% id)
is the single reference for viral detection, mapping, binning, taxonomy, host.

### 1.3 Read modes
- **SR/PE** (paired), **SR/SE** (single-end) — MEGAHIT + metaSPAdes.
- **LR** — metaFlye + hifiasm + (ONT) Medaka + metaMDBG.
Mapping: bwa-mem2 (SR) / minimap2 (LR). Mapping is against the **full** assembly and is
**shared** by viral (vRhyme) and prok (binners); viral/prok separation happens at the
**binning input** stage, not the mapping stage.

### 1.4 Co-assembly / co-binning (opt-in, `config.yaml → coassembly`)
Samples are grouped (metadata TSV `group` column | `all` | `none`). Per group:
- **Co-assembly** — MEGAHIT (SR) / metaFlye (LR) on the group's pooled reads →
  `coassembly/{group}/contigs.fa` (canonical, mode-agnostic).
- **Co-binning** (short reads only) — each sample mapped back → multi-sample abundance
  matrix → **VAMB** on the *viral-filtered* contigs → group MAGs → CheckM2 + GTDB.
- **Multi-split** (`cobinning_multisplit`) — VAMB over the concatenated *individual*
  assemblies (preserves strains), independent of co-assembly.
- **Viral consumer** (SR+LR) — full viral pipeline on the co-assembly → group vOTUs.
- Group MAGs/vOTUs get the same functional layer (AMR, defense, annotation, host, vRhyme).

VAMB group bins are written as **`*.fna`** (per-sample Binette bins are `*.fa`).

---

## 2. Quality control & assembly
- **SR QC**: fastp (adapter/quality trim). **LR QC**: NanoPlot + Porechop + Filtlong
  (`lr_min_len`, `lr_min_mean_q`).
- **Optional host removal**: bwa-mem2/minimap2 vs `host_genome`, before assembly.
- **Assembly (SR)**: MEGAHIT (`-m` bytes, `--min-contig-len MIN_CONTIG`, preset) +
  metaSPAdes (PE) / SPAdes `-s` (SE) + metaviralSPAdes (PE). Merged with tool prefixes,
  filtered `< MIN_CONTIG`, deduplicated by **MMseqs2 at `MIN_SEQ_ID` (95%)** → `rep_seq`.
- **Assembly (LR)**: metaFlye `--meta` (`--nano-raw`/`--nano-hq`/`--pacbio-hifi` by
  `lr_tech`/`lr_ont_chem`, `--min-overlap LR_FLYE_OVERLAP`) + hifiasm-meta + Medaka (ONT).

---

## 3. Viral detection & consensus

Three detectors run on `rep_seq` (or, for co-assembly, `coassembly/{group}/contigs.fa`):
- **VirSorter2** (`--min-score SCORE_VS2_MIN` = 0.5), **GeNomad** (`--min-score
  SCORE_GENOMAD_MIN` = 0.5, score calibration), **VIBRANT** (binary call, no numeric score).

**Consensus** (`viral_consensus_mode`, `viral_detection.smk::viral_consensus`):
- **count**: keep a contig called viral by **≥ `MIN_VIRAL_TOOLS`** tools (default 2 of 3).
- **score**: keep if any tool score ≥ its threshold (VS2 ≥ 0.5 OR GeNomad ≥ 0.5).
- **hybrid** (default): union — pass if it meets the count rule **OR** the score rule.

Contig IDs are normalized across tools (VS2 appends `||...`, GeNomad `|provirus_...`) before
intersection so the same contig from different tools is matched.

---

## 4. CheckV & provirus trimming

- **CheckV** `end_to_end` on the consensus set → `quality_summary.tsv` + `viruses.fna`
  (full viral contigs) + `proviruses.fna` (proviral regions **excised** from host, header
  `orig_id|start_end`).
- **CheckV quality tiers**: Complete / High-quality / Medium-quality / Low-quality /
  Not-determined (CheckV's own `completeness` + AAI/HMM evidence).
- **Provirus trim** (`viral_binning.smk::viral_nonredundant`; co-assembly:
  `coassembly_viral_trimmed`): for every contig, if CheckV trimmed it (provirus) use the
  trimmed sequence (`viruses.fna`/`proviruses.fna`), else keep the original — so host DNA
  flanking proviruses is removed before vOTU clustering. One contig can yield >1 trimmed
  entry (multiple provirus regions).

---

## 5. vOTU clustering (species-level viral OTUs)

`viral_binning.smk::skani_votu` → `skani_cluster`, ICTV / **Roux et al. 2019** definition:

> Two viral genomes are in the same vOTU iff **ANI ≥ `VOTU_ANI` (95%) AND
> max(af_query, af_ref) ≥ `VOTU_AF` (85%)**.

- ANI/aligned-fraction matrix from **skani triangle**.
- Clustering: connected components (union-find) over the edge set defined by the rule above.
- **Representative** per cluster = the member with the highest CheckV completeness
  (tie-break: earliest in the member list).
- Representative sets: `votu_all_reps.fasta` (all), `votu_mq_reps.fasta` (the
  **annotation subset** feeding taxonomy / PHIST / pharokka / genome maps),
  `votu_hq_10kb_reps.fasta` (HQ+ and ≥ 10 kb, for vConTACT3).
- **`votu_mq_reps` quality gate is configurable** via `viral_min_quality`
  (`complete | high | medium | low | not_determined`, default **medium** =
  Complete/HQ/MQ or completeness ≥ 50%). Lower it for high-novelty /
  short-fragment data (e.g. IonTorrent viromes) where most contigs are
  Low-quality/Not-determined and the default would leave taxonomy/host/annotation
  empty. Implemented in `pipeline_config.py::viral_keep_tiers`; the completeness
  ≥ 50% fallback applies only at the `medium` threshold. **vConTACT3 keeps its own
  HQ+/≥10 kb gate regardless** (genome-network clustering needs near-complete
  genomes). Same gate mirrored in `coassembly_viral_votu_reps`.

---

## 6. vRhyme (viral MAGs / vMAGs)

`viral_binning.smk::vrhyme` (co-assembly: `coassembly_vrhyme`, multi-sample). vRhyme bins
viral contigs into **vMAGs** using coverage (per-sample BAM, or **all group BAMs** for
differential coverage) + protein homology, minimum contig length `MIN_CONTIG`. CheckV is
re-run on the vMAGs (`checkv_vrhyme`).

---

## 7. Prokaryotic binning & MAG QC

- **Per-sample binners** (on the viral-filtered contig set, `filter_viral_for_prok`):
  MetaBAT2 + SemiBin2 → **Binette** (best-bin consolidation).
  VAMB was **removed** from per-sample (it is used only for co-binning, where its
  differential-coverage design applies).
- **Viral filter for prok** (`filter_viral_for_prok`): remove = viral_consensus **MINUS**
  provirus set (CheckV `provirus=Yes` ∪ GeNomad `|provirus_` suffix). Provirus-bearing
  contigs **stay** in prok input (the provirus is inside a host contig). Co-assembly mirror:
  `coassembly_filter_viral_for_prok`; VAMB then bins the nonviral subset and its abundance
  matrix is filtered to match.
- **CheckM2** — completeness/contamination (universal ML model).
- **MIMAG quality tiers** (`prok_binning.smk::checkm2`):
  - **High-quality (HQ)**: completeness ≥ 90% AND contamination ≤ 5%.
  - **Medium-quality (MQ)**: completeness ≥ 50% AND contamination ≤ 10%.
- **GUNC** — chimerism/contamination (progenomes DB), flagged in report.
- **MAG dereplication** (galah, `mag_derep`): cluster MAGs at **ANI ≥ `MAG_DEREP_ANI`
  (95%)**, keep the CheckM2-best representative per cluster.
- **GTDB-Tk** — taxonomy (bac120 + ar53).
- **Prok annotation quality gate (Bakta)**: MAGs are annotated by Bakta only if
  completeness ≥ `prok_min_completeness` (default **50%**, MIMAG-MQ) AND
  contamination ≤ `prok_max_contamination` (default **10%**) — the prokaryotic
  analogue of `viral_min_quality`. Lower `prok_min_completeness` for fragmented /
  high-novelty MAGs (e.g. IonTorrent) that rarely reach MQ. **AMR, defense and
  GTDB-Tk are NOT gated** — they run on every Binette/VAMB bin (or the whole
  contig set as a pseudo-genome under `low_depth_mode`). Legacy
  `bakta_min_completeness`/`bakta_max_contamination` keys still override the gate.

---

## 8. Co-binning abundance matrix (VAMB input)

`coassembly.smk::coassembly_abundance`. For a group, each sample's reads are mapped to the
co-assembly; per-contig depth from `jgi_summarize_bam_contig_depths` (`totalAvgDepth`
column). The multi-sample matrix is:

```
contigname   <s1>            <s2>            ...
contig_k     depth_k,s1      depth_k,s2      ...
```

fed to VAMB as `--abundance_tsv` (rows filtered to exactly match the `--fasta` contig set,
required by VAMB v5). VAMB uses composition (tetranucleotide) + this multi-sample
differential coverage in its variational autoencoder.

---

## 9. Taxonomy

### 9.1 Viral taxonomy merge (`taxonomy.smk::viral_taxonomy`)
Sources, combined by **deepest-recognized-rank-wins** over the 8-level ICTV scheme
(realm→kingdom→phylum→class→order→family→subfamily→genus):
- **vConTACT3** (genome clustering; family/genus/order predictions; Novel/Shared/Assigned).
- **MMseqs2/INPHARED** — real per-query LCA against an INPHARED seqTaxDB (`_mmseqs_lca_rollup`:
  per-protein LCA rolled up to per-contig by taking the longest/most-specific).
- **MMseqs2/custom** (optional, e.g. IMG/VR).
- **GeNomad** taxonomy.
- Source priority (tie-break at equal depth): **vConTACT3 > mmseqs_inphared >
  mmseqs_custom > genomad**. Output: `viral_taxonomy_merged.tsv` (`final_family` etc.).
- **Co-assembly core** uses **GeNomad + MMseqs2/INPHARED only** (vConTACT3/custom deferred);
  priority collapses to mmseqs_inphared > genomad.

### 9.2 Prokaryotic taxonomy
GTDB-Tk (primary) with an MMseqs2-LCA fallback (`mmseqs_taxonomy_prok`, `_mmseqs_lca_rollup`).

---

## 10. Abundance & normalization

### 10.1 Per-contig CoverM (`abundance.smk::coverm_viral`, `coverm_prok`)
CoverM on the sorted BAM. `coverm_method ∈ {rpkm, tpm, mean, covered_fraction}` (default
**rpkm**), plus `mean covered_fraction count` always emitted. Each metric is of the form:

```
value_i = (count_i / length_i_kb) * scale        # rpkm, tpm, mean depth
covered_fraction_i ∈ [0, 1]                       # bounded, NOT of that form
```

- **RPKM** = reads per kilobase per million mapped reads: `scale = 1e6 / total_mapped_reads`.
- **TPM** = transcripts per million (length-normalize first, then scale to sum 1e6).
- **mean** = mean read depth. **covered_fraction** = fraction of the contig with ≥1× coverage.
- `coverm_prok` uses **genome mode** (each MAG = one genome unit).

### 10.2 vOTU-level aggregation (`abundance.smk::votu_abundance`)
Members of a vOTU (95% ANI/85% AF) represent the same viral population; summing their
length-normalized values (rpkm/tpm/mean) would double-count. So:

1. **Sum raw read COUNTS** across all cluster members (counts are additive):
   `count_vOTU = Σ_{m ∈ cluster} count_m`.
2. Re-apply the per-sample normalization to the summed count. The per-sample **scale
   factor is backed out empirically** from the existing per-contig CoverM values (since
   `value_i = count_i/length_i_kb * scale`, `scale` is recovered as
   `median_i(value_i * length_i_kb / count_i)`), then
   `value_vOTU = (count_vOTU / length_rep_kb) * scale`.
3. **covered_fraction** is NOT summable (bounded) → the representative's own value is used.
- Degrades to a 1:1 identity mapping when vOTU clustering is off or found no clusters.

---

## 11. Diversity (`scripts/compute_diversity.py`)

Computed per domain (viral, prok) and combined, across all samples. The feature×sample
matrix is built from the per-sample abundance tables (`build_matrix`, values in the chosen
CoverM metric).

### 11.1 Alpha diversity (per sample)
For a sample with feature abundances `a = (a_1,…,a_S)`, total `N = Σ a_i`:

- **Richness** = number of features with `a_i > 0`.
- **Shannon**: `H = − Σ_{a_i>0} p_i · ln(p_i)`, with `p_i = a_i / N`. (H = 0 if N = 0.)
- **Simpson** (Gini-Simpson, 1−D): `1 − Σ a_i(a_i−1) / (N(N−1))`. (0 if N ≤ 1.)
- **Chao1** (integer counts; floats floored): with `S_obs` observed, `f1` singletons,
  `f2` doubletons:
  - `f2 > 0`:  `Chao1 = S_obs + f1² / (2·f2)`
  - `f2 = 0`:  `Chao1 = S_obs + f1·(f1−1) / 2`

### 11.2 Beta diversity (between samples)
- **Bray-Curtis** dissimilarity between samples `x`, `y`:
  `BC(x,y) = Σ_i |x_i − y_i| / (Σ_i x_i + Σ_i y_i)`  (0 if denominator 0).
- **PCoA** (classical MDS, `pcoa`): from the Bray-Curtis distance matrix `D`:
  1. Square: `D² `. 2. **Double-center**: `G = −½ (D² − rowmean − colmean + grandmean)`.
  3. Eigendecompose symmetric `G` → eigenvalues `λ`, eigenvectors `E` (descending).
  4. Keep positive `λ`: **coordinates** `X = E[:,pos] · √(λ_pos)`.
  5. Axis variance explained: `PC_k_var = λ_k / Σλ · 100 %` (first `n_axes = 5` reported).
- **Procrustes** (viral vs prok ordinations, `run_procrustes` / `_procrustes`):
  on the common samples, both PCoA coordinate matrices are centered and Frobenius-normalized
  (`mtx/‖mtx‖`), then the optimal **orthogonal rotation** `R` is found by SVD of `BᵀA`
  (`R = U·Vᵀ`), and **disparity** `M² = Σ (mtx1 − mtx2·R)²` measures viral↔prok
  ordination congruence (lower = more congruent).

---

## 12. Reads-based classification (sylph track)

`reads_classify.smk` (sylph): reads sketched (k-mer) and profiled against pre-built
databases (IMG/VR, UHGV, GTDB); taxonomy via **sylph-tax**. Outputs `sylphmpa`
(`|`-separated lineages, viral realm prefix `r__`). **No RPKM/TPM here** — sylph reports its
own two abundance types (both already **percentages, 0–100**, merged cross-sample with
`sylph-tax merge --column ...`):

- **`relative_abundance`** (taxonomic abundance): community proportion by **genome/organism**
  (length-corrected, like a taxonomic profile) — sums ≈ 100% across taxa.
- **`sequence_abundance`**: proportion of **reads/sequence** assigned (longer genomes get
  proportionally more) — sums ≈ 100%.

Key property: sylph parent-node abundances are **aggregate estimates, not sums of children**
(a class can read 99% while its individual genomes are ~0%).

Derived tables/metrics (`scripts/reads_classify/`):
- **Prevalence filter** (`filter_by_prevalence.py`): `prevalence_taxon = (#samples with
  value > 0) / (#samples)`; keep taxa with `prevalence > min_prevalence`
  (`reads_classify_min_prevalence`).
- **OTU table** (`make_otu.py`): reshape only (clade → `#OTU_ID`), no re-normalization
  (QIIME2/phyloseq-compatible).
- **Host collapse** (`collapse_by_host.py`): sum per-sample viral abundances grouped by
  predicted host genus (`Σ abundance over taxa sharing a host genus`) + taxon count.
- **BACPHLIP lifestyle** (`bacphlip_lifestyle.py`): per-genome `Virulent` if BACPHLIP score
  > `reads_classify_virulence_threshold` (0.5), else `Temperate`; sample-level
  `virulent_ratio = n_virulent / n_total`.

---

## 13b. Report-computed metrics (`scripts/report/`)

Metrics computed at report time (not stored by rules):

- **Novelty** (`renderer.py`): per sample, `unclassified = max(0, total_viral − classified)`
  (viral consensus contigs minus taxonomy-classified); `pct_novel = 100 × unclassified /
  total_viral`.
- **MIMAG bin counts** (`renderer.py`, from CheckM2): per sample count HQ (comp ≥ 90 AND
  cont ≤ 5), MQ (comp ≥ 50 AND cont ≤ 10), else LQ.
- **Host-collapse abundance** (`data_loaders.py::build_host_collapse`): viral abundance
  weighted by the CoverM metric (rpkm), summed over vOTU reps sharing a PHIST-predicted host.
- **QC percentages** (`data_loaders.py`): `gc_pct = gc_content × 100`; `adapter_pct =
  adapter_reads / reads_in × 100`; `bp_removed_pct = (bases_in − bases_out) / bases_in × 100`.

### Genome-map tracks (`scripts/genome_map.py`, sliding windows 500 bp / 100 bp step)
- **GC content**: `GC(w) = (n_G + n_C) / |w|` per window `w`.
- **GC skew**: `skew(w) = (n_G − n_C) / (n_G + n_C)` (0 if `n_G + n_C = 0`).
- Genome scale: 1 kb minor / 5 kb major ticks.

### vOTU membership table (`scripts/make_votu_table.py`)
One row per cluster member (rep included); representative annotations (CheckV quality,
taxonomy, lifestyle, host) are **propagated to all members** — no new metric, just the join
that makes the vOTU the reporting unit.

---

## 13. AMR / defense / annotation

- **AMR**: AMRFinderPlus (NCBI), RGI/CARD (`CARD_DB`), DeepARG (`DEEPARG_DB`), normalized to
  ARO by **argNorm**, combined by **amr_consensus** (`scripts/consolidate_amr.py`): hits are
  joined by CDS locus / normalized ARO, and **`consensus_score = n_tools_that_detected /
  N_TOOLS`** (N_TOOLS = 3). Curated (AMRFinder/RGI) vs exploratory (DeepARG) tiers are kept
  separate (never-merge rule). ABRicate → VFDB (virulence) + PlasmidFinder.
- **Defense**: DefenseFinder / AntiDefenseFinder (`DEFENSE_FINDER_MODELS_DB`) on prok ORFs
  and viral ORFs; dbAPIS (Diamond, `APIS_DB`) for viral anti-defense; defense islands computed.
- **Annotation**: Bakta (prok MAGs, `BAKTA_DB`, gated by completeness/contamination
  thresholds), eggNOG (COG/KEGG, `EGGNOG_DB`) → KEGG-Decoder; pharokka + phold (phage
  annotation, `PHAROKKA_DB`/`PHOLD_DB`, `pharokka_min_completeness`).
- Every AMR/defense/annotation tool consumes a shared **per-genome Prodigal** protein set
  (`prok_bin_proteins` / `coassembly_prok_bin_proteins`), predicted once per MAG.

---

## 14. Host prediction
**PHIST** — k-mer phage↔host prediction; viral genomes (vOTU MQ+ reps) vs the recovered
prokaryotic MAGs as candidate hosts (per-sample Binette bins; co-assembly: group VAMB MAGs).
iPHoP was removed. Co-assembly `coassembly_phist` links **group vOTUs → group MAGs** (needs
both viral + prok tracks + short reads).

---

## 15. Config → formula parameters (quick index)
| Param | Meaning | Default |
|---|---|---|
| `min_contig` | min contig length (bp) | 1000+ |
| `min_seq_id` | MMseqs2 dedup identity | 0.95 |
| `viral_consensus_mode` / `min_viral_tools` | consensus rule / N tools | hybrid / 2 |
| `score_vs2_min` / `score_genomad_min` | detector score thresholds | 0.5 / 0.5 |
| `votu_ani` / `votu_af` | vOTU ANI / aligned fraction (%) | 95 / 85 |
| `viral_min_quality` | CheckV tier gate for viral annotation subset (`complete`…`not_determined`) | medium |
| `prok_min_completeness` / `prok_max_contamination` | MAG quality gate for Bakta annotation (%) | 50 / 10 |
| `mag_derep_ani` | MAG dereplication ANI (%) | 95 |
| `coverm_method` | abundance metric | rpkm |
| CheckM2 MIMAG | HQ ≥90%/≤5%, MQ ≥50%/≤10% | — |
| `reads_classify_virulence_threshold` | BACPHLIP virulent cutoff | 0.5 |

---

*Maintenance: when a formula/threshold changes in code, update the matching section here
(each cites its source file/function). Co-assembly group rules (`coassembly_*`) mirror the
per-sample formulas verbatim unless a section says otherwise.*
