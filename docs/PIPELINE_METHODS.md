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
- **viral** — assembly-based viral detection, feeding a single global vOTU catalog
  built once across every sample and co-assembly group (§5).
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
- **SR/PE** (paired), **SR/SE** (single-end) — MEGAHIT.
- **LR** — one assembler by `lr_tech`: metaFlye + Medaka (ONT) or metaMDBG (HiFi).
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
- **Viral consumer** (SR+LR) — full viral detection pipeline runs on the co-assembly;
  the group's trimmed viral contigs are pooled into the same global vOTU catalog as every
  sample (§5), namespaced by group ID instead of sample ID.
- Group MAGs get the same functional layer (AMR, defense, annotation, host, vRhyme); group
  contigs feed vOTU clustering exactly as sample contigs do.

VAMB group bins are written as **`*.fna`** (per-sample Binette bins are `*.fa`).

---

## 2. Quality control & assembly
- **SR QC**: fastp (adapter/quality trim). **LR QC**: NanoPlot + Porechop + Filtlong
  (`lr_min_len`, `lr_min_mean_q`).
- **Optional host removal**: bwa-mem2/minimap2 vs `host_genome`, before assembly.
- **Assembly (SR)**: MEGAHIT (`-m` bytes, `--min-contig-len MIN_CONTIG`, preset), the
  single short-read assembler. Headers prefixed `MEGAHIT_`, filtered `< MIN_CONTIG`,
  deduplicated by **MMseqs2 at `MIN_SEQ_ID` (95%)** → `rep_seq`.
- **Assembly (LR)**: a single assembler chosen by `lr_tech`. **ONT** — metaFlye `--meta`
  (`--nano-raw`/`--nano-hq` by `lr_ont_chem`, `--min-overlap LR_FLYE_OVERLAP`), polished
  with Medaka. **HiFi** — metaMDBG `--in-hifi`, no polishing. There is no merge step:
  the assembler output feeds MMseqs2 directly.

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

## 5. Global vOTU catalog (species-level viral OTUs)

A vOTU is defined **once**, over the pooled viral sets of every sample and every
co-assembly group — not per sample. Clustering the same virus separately in each sample
it happens to appear in would assign it a different vOTU ID per sample, making richness
and presence/absence numbers incomparable across samples and inflated relative to the true
population count. `rules/votu_catalog.smk` (`votu_catalog_pool` → `votu_catalog_skani` →
`votu_catalog_cluster` → `votu_catalog_reps`) builds the catalog once; per-sample
presence is then derived from read recruitment against it (§5.2), not from re-clustering.
Pure logic lives in `scripts/votu_catalog.py`.

Measured on a 6-sample slice of the reference dataset: **9653 viral contigs collapse
into 5524 vOTUs (42.8% redundancy removed)**; summing each sample's own viral contig
count instead of using the catalog inflates apparent richness **1.75×** over the same
slice, because the same virus recovered in several samples is counted once per sample
rather than once per vOTU.

### 5.1 Pooling and clustering
- **`votu_catalog_pool`** concatenates every sample's `{sample}_viral_nonredundant.fasta`
  and every co-assembly group's trimmed viral set into one FASTA. Contig IDs are only
  unique within their own assembly, so every sequence is renamed
  `{source_id}|{contig_id}` (`source_id` = sample name or group name) before pooling —
  otherwise contigs from different assemblies that happen to share an ID would silently
  merge. A `provenance.tsv` records `source_type`/`source_id`/`member_id` for every
  pooled sequence.
- **`votu_catalog_skani`** runs `skani triangle --sparse` once over the whole pool.
  `--sparse` is required, not an optimization: skani's default dense matrix reports ANI
  only, with no aligned fraction, so the AF criterion below cannot be evaluated from it.
- **`votu_catalog_cluster`** applies the ICTV / **Roux et al. 2019** species definition
  (MIUViG standard):

  > Two viral genomes belong to the same vOTU iff **ANI ≥ `votu_ani` (95%) AND
  > max(af_query, af_ref) ≥ `votu_af` (85%)**.

  Clustering is single-linkage (connected components) over the edge set defined by that
  rule. **Representative** per cluster = the member with the highest CheckV completeness
  (tie-break: earliest in a deterministic traversal order), so vOTU IDs stay stable
  across reruns.
- **`votu_catalog_reps`** extracts two representative tiers, same quality gates as
  before but applied once over the global catalog instead of per sample:
  `catalog_all_reps.fasta` (all vOTUs), `catalog_mq_reps.fasta` (the **annotation
  subset** feeding taxonomy / PHIST / pharokka / genome maps).
- **`catalog_mq_reps` quality gate is configurable** via `viral_min_quality`
  (`complete | high | medium | low | not_determined`, default **medium** =
  Complete/HQ/MQ or completeness ≥ 50%). Lower it for high-novelty / short-fragment data
  (e.g. IonTorrent viromes) where most contigs are Low-quality/Not-determined and the
  default would leave taxonomy/host/annotation empty. Implemented in
  `pipeline_config.py::viral_keep_tiers`; the completeness ≥ 50% fallback applies only
  at the `medium` threshold.
  (genome-network clustering needs near-complete genomes).

### 5.2 Presence: two independent signals
Per-sample presence in the catalog is reported as two signals, side by side, **never
merged into one**, because they answer different questions:

- **`assembled`** — a contig from that sample's own assembly is a member of the vOTU
  (read straight off `provenance.tsv` + the cluster membership, no extra computation).
- **`recruited`** — that sample's reads, mapped competitively against the catalog's
  `all` representative set (`votu_catalog_map` → `votu_catalog_sort` →
  `votu_catalog_coverm`, bwa-mem2 for short reads / minimap2 for ONT and HiFi), cover
  ≥ `votu_presence_min_coverage` (default **75%**) of the representative's length —
  the breadth-of-coverage presence cutoff from **Roux et al. 2017**. This is what makes
  presence assembly-independent: a virus too low-abundance to assemble in a sample can
  still be detected here.

A vOTU can be `assembled` only, `recruited` only, `both`, or `absent` in a given sample.
Read identity for recruitment (`votu_recruit_min_identity`) defaults to 95% for
short/HiFi reads or 85% for ONT (configurable; `null` picks the technology-appropriate
default).

### 5.3 Output matrices (`votu_catalog_matrices`)
- **`presence_matrix.tsv`** — one row per vOTU, one column per sample, cell ∈
  `{assembled, recruited, both, absent}`.
- **`votu_abundance_matrix.tsv`** — one row per vOTU, one column per sample, cell =
  the configured CoverM metric (default rpkm) for that sample against the
  representative, when recruitment clears the presence cutoff.

Co-assembly groups have no reads of their own — recruitment is always per sample — so a
vOTU recovered only in a co-assembly group gets `assembled` presence attributed to the
group (visible in `provenance.tsv`) while its `recruited` presence in the matrix comes
from whichever samples' reads actually cover it; `presence_matrix.tsv` itself only has
sample columns.

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
- **MMseqs2/INPHARED** — real per-query LCA against an INPHARED seqTaxDB (`_mmseqs_lca_rollup`:
  per-protein LCA rolled up to per-contig by taking the longest/most-specific).
- **MMseqs2/custom** (optional, e.g. IMG/VR).
- **GeNomad** taxonomy.
- Source priority (tie-break at equal depth): **mmseqs_inphared >
  mmseqs_custom > genomad**. Output: `viral_taxonomy_merged.tsv` (`final_family` etc.).
- **Co-assembly core** uses **GeNomad + MMseqs2/INPHARED only** (MMseqs2/custom deferred);
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

### 10.2 Per-sample vOTU re-projection (`abundance.smk::votu_abundance`)
This is a **per-sample view** of the global catalog, kept alongside the catalog-level
matrices (§5.3) for per-sample reporting: it takes that sample's own per-contig CoverM
table and re-keys it onto the global vOTU clusters (`votu_catalog_cluster` output,
namespaced IDs stripped of this sample's own prefix). Members of a vOTU (95% ANI/85% AF,
§5.1) represent the same viral population; summing their length-normalized values
(rpkm/tpm/mean) would double-count. So:

1. **Sum raw read COUNTS** across all cluster members present in this sample (counts are
   additive): `count_vOTU = Σ_{m ∈ cluster ∩ sample} count_m`.
2. Re-apply the per-sample normalization to the summed count. The per-sample **scale
   factor is backed out empirically** from the existing per-contig CoverM values (since
   `value_i = count_i/length_i_kb * scale`, `scale` is recovered as
   `median_i(value_i * length_i_kb / count_i)`), then
   `value_vOTU = (count_vOTU / length_rep_kb) * scale`.
3. **covered_fraction** is NOT summable (bounded) → the max across this sample's members
   is used.
- Degrades to a 1:1 identity mapping for any contig missing from the global clusters file.
- This differs from `votu_abundance_matrix.tsv` (§5.3), which is recruitment-based (reads
  mapped against the catalog representative, whether or not this sample assembled it);
  `votu_abundance` here only sums contigs this sample actually assembled.

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
One row per cluster member (rep included), joined against the global catalog's
`vOTU_clusters.tsv` by member ID; representative annotations (CheckV quality, taxonomy,
lifestyle, host) are **propagated to all members** — no new metric, just the join that
makes the vOTU the reporting unit.

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
**PHIST** — k-mer phage↔host prediction; the global catalog's MQ+ representatives
(`catalog_mq_reps.fasta`, §5.1) vs the recovered prokaryotic MAGs as candidate hosts
(per-sample Binette bins; co-assembly: group VAMB MAGs). iPHoP was removed. Co-assembly
`coassembly_phist` links the same global vOTU representatives → group MAGs (needs both
viral + prok tracks + short reads).

---

## 15. Config → formula parameters (quick index)
| Param | Meaning | Default |
|---|---|---|
| `min_contig` | min contig length (bp) | 1000+ |
| `min_seq_id` | MMseqs2 dedup identity | 0.95 |
| `viral_consensus_mode` / `min_viral_tools` | consensus rule / N tools | hybrid / 2 |
| `score_vs2_min` / `score_genomad_min` | detector score thresholds | 0.5 / 0.5 |
| `votu_ani` / `votu_af` | vOTU ANI / aligned fraction (%), Roux et al. 2019 | 95 / 85 |
| `votu_catalog_enabled` | build the global vOTU catalog (§5) | true |
| `votu_presence_min_coverage` | % of representative covered by reads to count as `recruited` presence, Roux et al. 2017 | 75.0 |
| `votu_recruit_min_identity` | min. read identity (%) for recruitment; `null` = 95 (SR/HiFi) or 85 (ONT) | null |
| `viral_min_quality` | CheckV tier gate for viral annotation subset (`complete`…`not_determined`) | medium |
| `prok_min_completeness` / `prok_max_contamination` | MAG quality gate for Bakta annotation (%) | 50 / 10 |
| `mag_derep_ani` | MAG dereplication ANI (%) | 95 |
| `coverm_method` | abundance metric | rpkm |
| CheckM2 MIMAG | HQ ≥90%/≤5%, MQ ≥50%/≤10% | — |
| `reads_classify_virulence_threshold` | BACPHLIP virulent cutoff | 0.5 |

---

## Referências

- Roux, S., Adriaenssens, E. M., Dutilh, B. E., Koonin, E. V., Kropinski, A. M.,
  Krupovic, M., Kuhn, J. H., Lavigne, R., Brister, J. R., Varsani, A., Amid, C.,
  Aziz, R. K., Bordenstein, S. R., Bork, P., Breitbart, M., Cochrane, G. R.,
  Daly, R. A., Desnues, C., Duhaime, M. B., … Eloe-Fadrosh, E. A. (2019). Minimum
  Information about an Uncultivated Virus Genome (MIUViG). *Nature Biotechnology*,
  37(1), 29–37. https://doi.org/10.1038/nbt.4306
- Roux, S., Emerson, J. B., Eloe-Fadrosh, E. A., & Sullivan, M. B. (2017). Benchmarking
  viromics: an in silico evaluation of metagenome-enabled estimates of viral community
  composition and diversity. *PeerJ*, 5, e3817. https://doi.org/10.7717/peerj.3817

---

*Maintenance: when a formula/threshold changes in code, update the matching section here
(each cites its source file/function). Co-assembly group rules (`coassembly_*`) mirror the
per-sample formulas verbatim unless a section says otherwise.*
