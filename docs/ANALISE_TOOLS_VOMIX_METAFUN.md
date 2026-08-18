# Análise ferramenta a ferramenta — vOMIX-MEGA e metaFun

Complementa `docs/BENCHMARK_VOMIX_METAFUN.md` (arquitetura e performance) com o inventário completo e o *porquê* de cada escolha. Levantado do código dos dois repositórios, não das descrições.

---

## 0. O achado que atravessa os dois papers

Os dois grupos chegaram, por caminhos independentes, à mesma conclusão contraintuitiva: **mais ferramentas não melhora o resultado, e às vezes piora.**

| Etapa | Quem testou | Resultado |
|---|---|---|
| Identificação viral | vOMIX | geNomad sozinho: BA **0,9785**. Consenso de 6 por voto: **0,7010** |
| Montagem | metaFun | 3 configurações testadas, ficou com **uma** (MEGAHIT) |
| Binning | metaFun | 5 combinações testadas, refinamento ajuda **mas não de graça** (ver §4.2) |
| Perfil taxonômico | metaFun | manteve **duas** ferramentas — e é a exceção que confirma a regra (§4.3) |

O mecanismo é o mesmo nos dois casos: ferramentas com alta sensibilidade e baixa especificidade contaminam o consenso. Um voto por contagem **herda os falsos positivos de todos os votantes**.

---

## 1. vOMIX-MEGA — inventário por módulo

Módulo default: `viral-end-to-end`.

### 1.1 `preprocess`
| Ferramenta | Papel |
|---|---|
| **sra-tools + pigz** | download direto de acessos SRA |
| **fastp** | trimming/QC |
| **hostile** | descontaminação de hospedeiro, índice `human-t2t-hla` (opcional, `decontam-host`) |
| **MultiQC** | agregação |

### 1.2 `assembly`
**MEGAHIT** (default, `--prune-level 3`) ou **metaSPAdes**. Um ou outro, nunca ambos — `assembler:` no config, com validação que aborta se não for um dos dois.

### 1.3 `viral-identify` — **uma ferramenta só**

```
filter_contigs (>= 1000 bp)
   └─ geNomad end-to-end  --enable-score-calibration --relaxed --min-score 0.7
        └─ genomad_filter → viral.contigs.fa  (por amostra)
             └─ cat_contigs → pool GLOBAL
                  └─ cluster-fast (árvore binária)
                       └─ CheckV-PyHMMER
                            └─ combine_classifications → consensus_filtering → vOTUs
```

**Só o geNomad identifica.** O "consensus_filtering" do nome não é consenso entre ferramentas — é o cruzamento do score do geNomad com a qualidade do CheckV.

Justificativa no paper, em duas frentes:
- **acurácia**: geNomad tem sensibilidade 0,9743 **com** especificidade 0,9826; as alternativas trocam uma pela outra
- **memória**: *"omitting alternatives like DeepVirFinder that exceed 300 GB RAM allocations at high thread counts"*

### 1.4 `viral-benchmark` — **respondendo à pergunta das 6 ferramentas**

**Sim, é só benchmark.** Confirmado no gating do `workflow/Snakefile`:

```python
# viral-identify — entra no default
if config["module"] in ["viral-identify", "viral-end-to-end", "run-all"]:

# viral-benchmark — NÃO entra no viral-end-to-end
if config["module"] in ["viral-benchmark", "run-all"]:
```

Como o default é `viral-end-to-end`, as 6 ferramentas **nunca rodam** a menos que você peça `--module viral-benchmark` explicitamente.

As 6: geNomad, DeepVirFinder, PhaMer (PhaBOX2), VirSorter2, VirFinder, VIBRANT. O paper avalia ainda PPR-Meta e Seeker fora do módulo.

O propósito é declarado: *"to provide users with direct comparative controls"* — é uma ferramenta de auditoria para o usuário reproduzir o benchmark no próprio dado, não um caminho de produção. É honesto: eles não escondem as alternativas, entregam o meio de contestar a escolha deles.

**Isso é o que a vapor poderia copiar sem custo:** manter VirSorter2 e VIBRANT como um módulo de benchmark opcional em vez de no caminho default.

### 1.5 `viral-taxonomy` — **duas fontes, 183 linhas**

```
geNomad classify (reaproveitado) ──┐
                                   ├─ merge_taxonomy (ete3) → taxonomia final
PhaGCN (PhaBOX2)  ─────────────────┘
```

Nenhuma busca por homologia (sem MMseqs2, sem Diamond, sem INPHARED). O `ete3` resolve a hierarquia NCBI.

Comparação direta com a vapor:

| | vOMIX | vapor |
|---|---|---|
| Fontes | 2 (geNomad + PhaGCN) | 3 (geNomad + MMseqs2/INPHARED + MMseqs2/custom) |
| Linhas | 183 | ~650 |
| Estratégia | merge por hierarquia NCBI | rank mais profundo vence, desempate por confiança |

A vapor é mais rica em evidência; o vOMIX é mais barato. Nenhum dos dois usa vConTACT — a vapor removeu em 2026-08-17, o vOMIX tem `vcontact2.yml` no `envs/` mas **nenhuma regra o invoca**.

### 1.6 `viral-host`
| Ferramenta | Papel | Default |
|---|---|---|
| **CHERRY** (PhaBOX2) | predição de hospedeiro | sim |
| **PhaTYP** (PhaBOX2) | estilo de vida (virulento/temperado) | sim |
| **iPHoP** | hospedeiro, mais sensível | opcional (`iphop-host`) |

Sem PHIST. O CHERRY é rede neural sobre grafo; o PHIST é k-mers compartilhados — abordagens diferentes, complementares.

### 1.7 `viral-annotate` — **e o que está desligado**
| Ferramenta | Papel | No alvo do módulo? |
|---|---|---|
| **prodigal-gv** | ORFs (paralelizado por chunks) | sim |
| **eggNOG-mapper v2** | ortologia/COG/KEGG | sim |
| **PhaVIP** (PhaBOX2) | anotação de proteína viral | sim |
| **MetaCerberus** | HMM contra PVOG + VOG | sim |
| **Pharokka** | anotação de fago (PHROGS) | sim |
| **VirSorter2 + DRAM-v** | AMGs com curadoria | **comentado no `done_log`** |

O DRAM-v está implementado (`dramv_annotate`, `dramv_distill`) mas **fora do alvo** — as linhas estão comentadas em `viral-annotate.smk`. Provavelmente por custo: DRAM-v exige rodar o VirSorter2 de novo só para gerar o input.

Nota: **não há `phold`**. A vapor tem pharokka→phold (anotação estrutural por ProstT5), que é mais moderno que o par pharokka+PhaVIP.

### 1.8 `viral-binning` e `viral-community`
**vRhyme** (+ prodigal-gv) → **dRep**. Abundância por **CoverM** com `--mapper minimap2-sr --min-read-percent-identity 95 --min-read-aligned-percent 75 --trim-min 10 --trim-max 90`.

### 1.9 `prok-binning`, `prok-community`, `prok-annotate`
| Etapa | Ferramentas |
|---|---|
| Mapeamento | **strobealign** |
| Cobertura | `jgi_summarize_bam_contig_depths` |
| Binning | **VAMB + MetaBAT2 + MaxBin2 + CONCOCT** (quatro!) |
| Consolidação | **DAS Tool** |
| Qualidade | **CheckM2** |
| Dereplicação | **galah** |
| Taxonomia | **GTDB-Tk** |
| Comunidade | **MetaPhlAn 4** |
| Função | **HUMAnN 3** (com 7 mapeamentos: EC, eggNOG, GO, KO, Pfam, MetaCyc) |

**Contradição interna do vOMIX:** eles pregam "menos ferramentas, melhor escolhida" no lado viral e usam **quatro binners** no lado procariótico, sem benchmark que justifique. O rigor deles é assimétrico — o paper é sobre viroma, e o lado procarioto parece herdado sem a mesma auditoria.

---

## 2. metaFun — inventário por módulo

Sete módulos independentes, executados um a um pelo usuário.

### 2.1 `RAWREAD_QC`
**FastQC** (bruto) → **fastp** → **bowtie2** (remoção de hospedeiro) → **FastQC** (filtrado) → **MultiQC**.

### 2.2 `ASSEMBLY_BINNING`
**MEGAHIT** → **bowtie2** (map-back) → **MetaBAT2** + **SemiBin2** → **DAS Tool**.

### 2.3 `BIN_ASSESSMENT`
```
CheckM2 → GUNC → filtro de qualidade combinado → GTDB-Tk → metadados
```
Exatamente a cadeia da vapor (`checkm2 → gunc → gtdbtk`). Convergência independente.

### 2.4 `WMS_TAXONOMY`
**Kraken2 + Bracken** (confidence 0,1 default) **ou** **sylph** (c200, GTDB r220), à escolha por `--profiler`. Default: **sylph**.

### 2.5 `WMS_FUNCTION`
**HUMAnN 3** + parsing + análise estatística em R.

### 2.6 `COMPARATIVE_ANNOTATION` — o módulo mais denso
```
prepare_genomes → prokka (por genoma)
                     └─ PPanGGOLiN  ←── HUB
                          ├─ KofamScan  → matriz KO → completude de módulo KEGG
                          ├─ VFDB       → fatores de virulência
                          ├─ RGI/CARD   → AMR
                          ├─ dbCAN      → CAZymes
                          ├─ eggNOG     → ortologia
                          ├─ DefenseFinder → sistemas de defesa
                          ├─ Scoary2    → GWAS microbiano
                          └─ genePA cluster → PCoA interativo
skani (ANI) e dRep correm em paralelo
```

> ⚠️ **`panaroo` está definido no `.nf` mas nunca é chamado no bloco `workflow {}`.** É código morto, como as linhas `//conda` comentadas. Quem roda é o **PPanGGOLiN**.

### 2.7 `WMS_STRAIN`
**inStrain** (`profile` + `compare`), com bowtie2, prodigal e eggNOG-mapper como apoio. Microdiversidade e SNVs — capacidade que nem a vapor nem o vOMIX têm.

---

## 3. Respostas diretas

### 3.1 Por que o metaFun só usa MEGAHIT?

Testaram metaSPAdes (default) contra duas configurações de MEGAHIT, em três ambientes CAMI simulados (n=5 cada). Resultado:

- MEGAHIT venceu em **comprimento total de contigs** e **fração de bp alinhados à referência**, com vantagem maior na rizosfera (alta complexidade)
- metaSPAdes só ganhou em **comprimento do contig mal-montado** no gut humano e no marinho

Decidiram *"based on the overall metrics and resource efficiency"* — ou seja, o desempate foi **custo computacional**, não só acurácia.

**Para a vapor:** hoje rodamos MEGAHIT **e** metaSPAdes **e** metaviralSPAdes, depois deduplicamos com MMseqs2. Isso pode ser defensável (montadores diferentes recuperam contigs diferentes; a dedup resolve a redundância), mas **nunca foi medido**. É o mesmo padrão do consenso viral: acumular ferramentas sem evidência de que somam.

### 3.2 Por que só GTDB-Tk na classificação?

Não é "só GTDB-Tk" — são **duas camadas com propósitos diferentes**, e é aí que está a lição:

| Camada | Ferramenta | Objeto | Pergunta que responde |
|---|---|---|---|
| MAG | GTDB-Tk | genoma montado | "que organismo é este MAG?" |
| Reads | sylph / Kraken2+Bracken | leituras brutas | "o que existe na amostra e em que proporção?" |

São perguntas distintas — não há redundância. A vapor faz o mesmo (GTDB-Tk nos bins, sylph nos reads).

E o achado metodológico que vale independente de ferramenta: **Kraken2 é violentamente sensível a parâmetro.** Variando confidence e filtro de abundância, as métricas oscilaram **3× em distância composicional e 45× em F1**. O sylph mostrou variação mínima. Por isso o metaFun mantém os dois: sylph venceu em acurácia *e* eficiência, mas **Kraken2 classifica reads individuais**, coisa que o sylph não faz. Mantiveram por **capacidade**, não por acurácia — a exceção que confirma a regra do §0.

### 3.3 PPanGGOLiN serve para a vapor?

**Depende de uma pergunta que a vapor hoje não faz.**

O PPanGGOLiN constrói o pangenoma de um conjunto de genomas da **mesma espécie**, particionando genes em persistente/shell/cloud. No metaFun ele é o hub: KO, VFDB, CARD, dbCAN e Scoary2 todos consomem a saída dele.

Para servir, a vapor precisaria:
1. agrupar MAGs por espécie entre amostras (o `galah` já derreplica, mas não agrupa por espécie para pangenoma)
2. ter **vários MAGs da mesma espécie** — o que exige muitas amostras do mesmo ambiente

**Meu parecer:** não é prioridade. A vapor é uma pipeline de *viroma* com trilha procariótica de apoio; o pangenoma é genômica comparativa, um projeto diferente. O metaFun é uma ferramenta de genômica comparativa com trilha metagenômica de apoio — o eixo é invertido.

**Mas há um alerta do paper que se aplica à vapor hoje**, independente de pangenoma:

> Fluidez genômica em MAG incompleto é armadilha. *Bordetella holmesii* saltou de **0,002 → 0,148** no conjunto "Fragmented Incomplete".

Qualquer interpretação comparativa a partir de MAG fragmentado precisa dessa ressalva. E o threshold de core gene a 99% praticamente **zera** o core genome com MAGs — eles usam 90%.

### 3.4 Como cada um faz anotação?

| | vOMIX (viral) | metaFun (procarioto) | vapor |
|---|---|---|---|
| ORFs | prodigal-gv | prokka | prodigal-gv / bakta |
| Ortologia | eggNOG-mapper | eggNOG-mapper | eggNOG-mapper |
| Fago | Pharokka + PhaVIP | — | Pharokka + **phold** |
| HMM | MetaCerberus (PVOG/VOG) | KofamScan (KEGG) | — |
| AMR | — | RGI/CARD | AMRFinderPlus + RGI + DeepARG + ABRicate + argNorm |
| Virulência | — | VFDB | ABRicate (VFDB) |
| Defesa | — | DefenseFinder | DefenseFinder + AntiDefenseFinder + dbAPIS |
| CAZymes | — | dbCAN | — |
| AMG | DRAM-v (**desligado**) | — | — |
| Pangenoma | — | PPanGGOLiN + Scoary2 | — |

**A vapor é a mais completa em AMR/defesa dos três, com folga.** Nenhum dos dois tem anti-defesa (dbAPIS) nem consenso de AMR normalizado (argNorm).

Lacunas reais da vapor: **CAZymes (dbCAN)** e **AMGs curados (DRAM-v)** — e o vOMIX desligou o DRAM-v, o que sugere que o custo não compensou nem para eles.

### 3.5 Como o vOMIX fica rápido e com pouca memória?

Quatro mecanismos, detalhados em `BENCHMARK_VOMIX_METAFUN.md` §3. Resumo com os números:

| Mecanismo | Ganho medido |
|---|---|
| **Escolha por estabilidade de memória** — geNomad em vez de DVF | evita picos de 300+ GB com muitas threads |
| **CheckV-PyHMMER** — pré-popula o `tmp/` do CheckV com prodigal-gv + 80 shards de PyHMMER, depois roda `checkv end_to_end -t 1` | 300k contigs: 72 h/939 GB → 19,2 h/34 GB (**3,7× tempo, 27× memória**) |
| **Árvore binária no clustering** — divide-and-conquer com merge estilo torneio, contra o O(N²) do CD-HIT/megablast | **−63% tempo, −64% RSS**, com 93% de identidade nos representantes |
| **Data splitting** (`--contig-splits`, `--checkv-splits`) | com `--checkv-splits 7`, pico cai para 18 GB |

O princípio comum: **trocar um job grande por muitos jobs pequenos e independentes**, deixando o escalonador do Snakemake paralelizar. O CheckV-PyHMMER é o exemplo mais claro — com 8 splits vira 640 jobs de 1 thread e 4 GB.

Ressalva honesta, declarada por eles: em n=1.000 e n=10.000 o CheckV **nativo é mais rápido** (58 s vs 192 s), por overhead de inicialização dos bindings Python. O ganho só aparece acima de ~10⁵.

---

## 4. O que fazer com isso na vapor

### 4.1 Detecção viral — o item de maior impacto
Rodar geNomad-sozinho contra os três modos de consenso num mock com verdade conhecida, medindo **especificidade** e não só sensibilidade. Se o padrão do vOMIX se confirmar, caem VirSorter2, VIBRANT e três modos de config.

Meio-termo que o vOMIX oferece de graça: mover as ferramentas extras para um **módulo de benchmark opcional**, como eles fazem. Preserva a capacidade de auditar sem pagar o custo em todo run.

### 4.2 Binning — cuidado, aqui a lição é mais sutil

O metaFun **não** concluiu que menos binners é melhor. Concluiu:

- MetaBAT2 sozinho gerou **menos MAGs** e teve o **pior ARI**
- SemiBin2 em modo self-supervised melhorou ARI em metagenomas complexos
- **DAS Tool melhorou a acurácia de atribuição de bp, mas com leve queda de ARI** no marinho e na rizosfera — no gut humano foi o oposto

Ou seja: **refinamento não é de graça**, e o efeito depende do ambiente. Eles escolheram `DAS Tool(MetaBAT2 + SemiBin2 self-supervised)`.

A vapor usa `Binette(MetaBAT2 + SemiBin2)` — a mesma família de escolha, com o Binette no lugar do DAS Tool. **Está alinhada com a evidência do metaFun.** Não vejo motivo para mexer, mas vale saber que o refinamento tem custo em ARI e que isso varia por ambiente.

### 4.3 Montagem — a lacuna não medida
A vapor roda três montadores. O metaFun mediu e escolheu um. Não estou dizendo que a vapor está errada — a dedup por MMseqs2 muda o cálculo. Estou dizendo que **nunca foi medido**, e é o mesmo padrão do consenso viral.

### 4.4 Adoções concretas, por razão custo/benefício

| Item | De onde | Por quê |
|---|---|---|
| Módulo de benchmark viral opcional | vOMIX | preserva auditoria, tira custo do caminho default |
| **PhaTYP** | vOMIX | estilo de vida sobre vOTUs montados (o BACPHLIP da vapor só cobre genomas de referência do sylph) |
| **CHERRY** | vOMIX | complementa o PHIST por abordagem diferente (grafo vs. k-mer) |
| **hostile** | vOMIX | índices de hospedeiro curados |
| **CheckV-PyHMMER** | vOMIX | só acima de ~10⁵ contigs; opcional via config |
| **dbCAN** | metaFun | CAZymes, lacuna real |
| **inStrain** | metaFun | microdiversidade — capacidade nova, custo alto |
| ~~PPanGGOLiN~~ | metaFun | escopo diferente; ver §3.3 |
| ~~DRAM-v~~ | vOMIX | eles próprios desligaram |

---

## Referências

Lee, H. G., Song, J. Y., Yoon, J., Chung, Y., Kwon, S.-K., & Kim, J. F. (2026). metaFun: An analysis pipeline for metagenomic big data with fast and unified functional searches. *Gut Microbes*, 18(1), 2611544. https://doi.org/10.1080/19490976.2025.2611544

Shekarriz, E., Vijendran, E., & Ho, J. W. K. (2026). vOMIX-MEGA: An ultra-fast end-to-end pipeline for terabyte-scale viral metagenomics analysis. *bioRxiv*. https://doi.org/10.64898/2026.07.28.741255

Camargo, A. P., Roux, S., Schulz, F., Babinski, M., Xu, Y., Hu, B., Chain, P. S. G., Nayfach, S., & Kyrpides, N. C. (2024). Identification of mobile genetic elements with geNomad. *Nature Biotechnology*, 42(8), 1303–1312. https://doi.org/10.1038/s41587-023-01953-y

Gauthier, J., Vincent, A. T., Charette, S. J., & Derome, N. (2019). A brief history of bioinformatics. *Briefings in Bioinformatics*, 20(6), 1981–1996. https://doi.org/10.1093/bib/bby063

Gautreau, G., Bazin, A., Gachet, M., Planel, R., Burlot, L., Dubois, M., Perrin, A., Médigue, C., Calteau, A., Cruveiller, S., Matias, C., Ambroise, C., Rocha, E. P. C., & Vallenet, D. (2020). PPanGGOLiN: Depicting microbial diversity via a partitioned pangenome graph. *PLOS Computational Biology*, 16(3), e1007732. https://doi.org/10.1371/journal.pcbi.1007732

Olm, M. R., Crits-Christoph, A., Bouma-Gregson, K., Firek, B. A., Morowitz, M. J., & Banfield, J. F. (2021). inStrain profiles population microdiversity from metagenomic data and sensitively detects shared microbial strains. *Nature Biotechnology*, 39(6), 727–736. https://doi.org/10.1038/s41587-020-00797-0

Pan, S., Zhao, X.-M., & Coelho, L. P. (2023). SemiBin2: self-supervised contrastive learning leads to better MAGs for short- and long-read sequencing. *Bioinformatics*, 39(Supplement_1), i21–i29. https://doi.org/10.1093/bioinformatics/btad209

Shaw, J., & Yu, Y. W. (2024). Rapid species-level metagenome profiling and containment estimation with sylph. *Nature Biotechnology*, 43(7), 1122–1131. https://doi.org/10.1038/s41587-024-02412-y

Shang, J., Peng, C., Tang, X., & Sun, Y. (2024). PhaBOX: a web server for identifying and characterizing phage contigs in metagenomic data. *Bioinformatics Advances*, 4(1), vbae052. https://doi.org/10.1093/bioadv/vbae052

Sieber, C. M. K., Probst, A. J., Sharrar, A., Thomas, B. C., Hess, M., Tringe, S. G., & Banfield, J. F. (2018). Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. *Nature Microbiology*, 3(7), 836–843. https://doi.org/10.1038/s41564-018-0171-1
