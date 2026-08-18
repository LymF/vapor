# vapor × vOMIX-MEGA × metaFun — análise comparativa de arquitetura e ferramentas

Documento de referência interno. Base:

- **vOMIX-MEGA** — repo `vOMIX-MEGA/` (commit shallow, 2026-08), paper: Shekarriz E, Vijendran E, Ho JWK. *vOMIX-MEGA: An ultra-fast end-to-end pipeline for terabyte-scale viral metagenomics analysis.* bioRxiv, 2026. https://doi.org/10.64898/2026.07.28.741255 — texto em `docs/refs/vOMIX-MEGA_biorxiv_2026.pdf`
- **metaFun** — repo `metaFun/`, paper: Lee HG, Song JY, Yoon J, Chung Y, Kwon S-K, Kim JF. *metaFun: An analysis pipeline for metagenomic big data with fast and unified functional searches.* *Gut Microbes*, 18(1), 2611544, 2026. https://doi.org/10.1080/19490976.2025.2611544 — texto em `docs/refs/metaFun_GutMicrobes_2026.md`
- **vapor** — este repo, `Snakefile` + `rules/*.smk` + `envs/*.yaml`

---

## 1. Comparação de arquitetura

| | **vapor** | **vOMIX-MEGA** | **metaFun** |
|---|---|---|---|
| Motor | Snakemake | Snakemake 9.x | Nextflow 24.04 |
| Interface | `snakemake` direto (+ `vapor.py`) | CLI `vomix` (click) que monta a chamada Snakemake | CLI `metafun` (bash) que chama `nextflow run` |
| Isolamento | 26 envs conda + `container:` por regra (Docker) | ~40 envs conda por regra + 1 `.sif` Apptainer global | 5 imagens `.sif` Apptainer, **zero** conda para as ferramentas |
| Unidade de trabalho | por amostra, com catálogo global de vOTU no final | por amostra até `viral.contigs.fa`, depois **pool global** | por amostra, módulos independentes encadeados manualmente |
| Modularidade | alvos Snakemake (`results/.../done.txt`) | `--module` seleciona quais `include:` entram no DAG | módulos = scripts `.nf` separados, executados um a um |
| Escopo | viral + procarioto + defesa/AMR + relatório HTML | viral (foco) + procarioto básico | procarioto (foco), genômica comparativa, **sem viral** |
| Long reads | sim (Flye, hifiasm, Medaka, minimap2) | não | não |

**Leitura rápida:** metaFun não é concorrente da vapor no viroma — é complementar (genômica comparativa/pangenoma/strain-level, que a vapor não tem). vOMIX-MEGA é o concorrente direto, e o interesse dele é quase todo de *engenharia*, não de escopo biológico: o pipeline dele faz **menos** que a vapor, mas roda em outra ordem de grandeza.

---

## 2. Ferramentas — o que é diferente

### 2.1 Onde as três coincidem
MEGAHIT/metaSPAdes, CheckV, geNomad, MetaBAT2, CheckM2, GTDB-Tk, DAS Tool/Binette, CoverM, MultiQC, eggNOG-mapper, prodigal-gv.

### 2.2 Ferramentas que **vOMIX-MEGA tem e a vapor não**

| Ferramenta | Etapa | Vale importar? |
|---|---|---|
| **hostile** | remoção de hospedeiro (índices human-t2t-hla prontos) | **Sim.** Substitui a `host_removal.smk` artesanal por índice curado e versionado |
| **strobealign** | mapeamento SR para binning | **Sim, avaliar.** Substituto direto do bwa-mem2, ~3× mais rápido em leitura curta |
| **PhaBOX2** (PhaMer, PhaGCN, CHERRY, PhaTYP, PhaVIP) | ident./taxonomia/hospedeiro/estilo de vida/anotação proteica | **Sim, em parte.** CHERRY complementa o PHIST. PhaTYP cobriria estilo de vida **de vOTUs montados** — a vapor tem BACPHLIP, mas só na trilha de reads e só sobre genomas de referência do sylph (`reads_classify_genome_fasta`) |
| **DRAM-v** (+ VirSorter2 `--prep-for-dramv`) | AMGs virais com curadoria | Talvez. Sobrepõe pharokka/phold, mas o *distill* de AMG é melhor |
| **MetaCerberus** | anotação HMM multi-DB | Baixa prioridade |
| **CONCOCT** | binning procarioto | Não — Binette já consolida MetaBAT2+SemiBin2 |
| **galah** | dereplicação de MAGs | Já temos (`env_derep.yaml`) |
| **MetaPhlAn 4 + HUMAnN 3** | perfil taxonômico/funcional por reads | Concorre com nosso sylph; ver §4 |
| **iPHoP** | predição de hospedeiro (opcional) | **Sim.** Mais sensível que o PHIST sozinho |

### 2.3 Ferramentas que **metaFun tem e a vapor não**

| Ferramenta | Etapa | Comentário |
|---|---|---|
| **sylph** (c200, GTDB r220) | perfil taxonômico de abundância | **Já temos** (`reads_classify.smk`) — e metaFun valida a escolha (ver §5) |
| **Kraken2 + Bracken** | perfil de abundância de sequência | Complementar ao sylph — Kraken2 classifica *reads individuais*, sylph não |
| **GUNC** | quimerismo em MAGs | **Já temos** (`env_gunc.yaml`) |
| **prokka / panaroo / PPanGGOLiN** | pangenoma | Gap real, se quisermos genômica comparativa |
| **skani** | ANI entre genomas | **Já temos** — usamos no `votu_catalog.smk` |
| **inStrain** (profile + compare) | microdiversidade / SNVs / strain-level | **Gap interessante.** Nada equivalente na vapor |
| **KofamScan, dbCAN, Scoary2** | anotação KO / CAZy / GWAS microbiano | Parcial (temos eggNOG/KEGG) |
| **DefenseFinder, RGI/CARD, VFDB** | defesa e AMR | **Já temos os três** — a vapor é mais completa aqui |

### 2.4 Ferramentas que **só a vapor tem**

vRhyme (binning viral), vConTACT3, MMseqs2 taxonomy, VIBRANT/VirSorter2 em consenso, COBRA (extensão de contigs), phold, bakta, dbAPIS (anti-defesa), AMRFinderPlus, DeepARG, ABRicate, argNorm, ilhas de defesa, relatório HTML interativo, suporte a long reads. **Em escopo biológico a vapor é a mais ampla das três.**

---

## 3. A pergunta central: como o vOMIX-MEGA fica rápido

O paper reivindica 54 min / 20,1 GB de pico em 100.000 contigs contra 464 min / 384 GB (ViroProfiler) e 998 min / 223 GB (Nayfach et al.); VIRify e ViWrap não terminaram em 14 dias. São **quatro** mecanismos independentes, e todos são portáveis para a vapor.

### 3.1 Clustering por redução em árvore binária (`--cluster-iter`)
`workflow/rules/cluster-fast.smk`

O gargalo: CD-HIT e o MEGABLAST all-vs-all do CheckV têm RSS que escala com **O(N²)** no número de sequências. Em escala de terabyte isso trava.

A solução é divide-and-conquer com merge estilo torneio:

1. `seqkit split2` parte o FASTA em `2^(iter-1)` chunks.
2. Cada chunk é clusterizado **em paralelo e isoladamente** (camada 1).
3. Camada *L* pega os representantes de **dois** chunks da camada *L−1*, concatena e re-clusteriza.
4. Repete até sobrar `layer_N/chunk_0.fa`.

O `get_iter_inputs(wildcards)` é o coração: para `layer > 1` devolve `chunk_{2c}` e `chunk_{2c+1}` da camada anterior. Snakemake resolve a árvore inteira sozinho a partir do alvo final.

Resultado medido pelo paper: **−63% de tempo e −64% de RSS de pico** vs. o MEGABLAST original, com **93% de identidade** nos representantes finais entre 1, 2, 4 e 8 iterações (Supplementary Figure 6). Ou seja: o custo em acurácia é pequeno e mensurável.

> **Nota para a vapor:** nosso `votu_catalog.smk` usa `skani triangle --sparse`, que já é sub-quadrático em memória e muito melhor que CD-HIT/megablast. O problema aqui é *outro*: `--slow` com `-t {THREADS}` num único job. A árvore binária ainda ajudaria acima de ~10⁵ genomas, mas **não é a nossa prioridade** — o skani já resolve o pior do problema que eles atacaram.

### 3.2 CheckV-PyHMMER
`workflow/rules/checkv-pyhmmer.smk` + `workflow/scripts/pyhmmer_wrapper.py`

Este é o ganho maior, e o truque é elegante: eles **não modificam o CheckV**. Eles pré-populam o diretório `tmp/` que o CheckV usaria e depois rodam `checkv end_to_end -t 1`, que detecta os arquivos e pula as etapas caras.

Sequência exata das regras:
1. `split_contigs` — `seqkit split2` em `checkv-splits + 1` partes.
2. `checkv_prodigalgv` — gera `split-{part}/checkv/tmp/proteins.faa`.
3. `checkv_pyhmmer` — **um job por (part × shard de HMM)**. O DB do CheckV já vem em 80 arquivos `checkv_hmms/{1..80}.hmm`, então com 8 splits isso vira **640 jobs independentes de 1 thread e 4 GB** cada. É aí que o paralelismo real acontece.
4. `checkv_hmm_merge` — `cat` dos 80 `.hmmout` em `tmp/hmmsearch.txt`.
5. `checkv` — `checkv end_to_end -t 1` sobre o `tmp/` já preenchido.
6. `checkv_merge` — concatena FASTAs e faz row-bind dos `quality_summary.tsv`.

O `pyhmmer_wrapper.py` tem uma decisão de memória explícita: se `proteins.faa + hmmdb < 10%` da RAM disponível, faz `seq_file.read_block()` (pré-carrega tudo na memória); caso contrário, faz streaming do `SequenceFile`. Também expõe `--z_flip`, que ajusta o valor `Z` para que um `hmmsearch` produza resultado idêntico a um `hmmscan` — hmmsearch é muito mais rápido.

Números do paper, 64 CPUs, n=300.000 contigs:

| | CheckV nativo | CheckV-PyHMMER |
|---|---|---|
| Walltime | 72,0 h | 19,2 h (**3,7×**) |
| RSS de pico | 939,1 GB | 34,2 GB (**27,4×**) |
| Memória virtual | 1,31 TB | 36,1 GB |

Com `--checkv-splits 7`, o pico cai para **18 GB**. Contrapartida honesta e declarada: em n=1.000 e n=10.000 o CheckV nativo é *mais rápido* (58 s vs. 192 s) por causa do overhead de inicialização dos bindings Python — o ganho só aparece a partir de ~10⁵.

> ⚠️ **Discrepância encontrada entre paper e código.** O paper afirma que substituíram o Prodigal por **Pyrodigal-gv** com "paralelização OpenMP". O código (`workflow/scripts/parallel_prodigal_gv.py`, linha ~68) na verdade faz `subprocess.run(["prodigal-gv", ...])` sobre chunks num `ThreadPoolExecutor`, e `workflow/envs/prodigal-gv.yml` contém apenas `prodigal-gv=2.11.0` — **nenhum pyrodigal**. É paralelismo por fork de processo, não a biblioteca. Isso importa para nós: a vapor já tem `pyrodigal-gv=0.3.2` em `env_annotation.yaml`, `env_genomad.yaml` e `env_vcontact3.yaml`, então podemos implementar a versão que o paper *descreve* — que é melhor que a que eles entregaram.

### 3.3 Data splitting genérico (`--contig-splits`, `--checkv-splits`)
O mesmo padrão split→processa em paralelo→concatena sem perda é aplicado também ao geNomad (`--splits`) e às ferramentas de identificação do módulo de benchmark. Baixa complexidade, ganho direto em RSS de pico.

### 3.4 Escolha deliberada de ferramenta por estabilidade de memória

Este é o ponto mais contraintuitivo e o mais relevante para a vapor. Citando o paper:

> "`vomix viral-identify` uses geNomad as its core identification framework, omitting alternatives like DeepVirFinder that exceed 300 GB RAM allocations at high thread counts."

E, sobre consenso multi-ferramenta:

> "consensus algorithms reduced taxonomic accuracy; the vOMIX-MEGA (multi-tool) voting scheme dropped BA to 0.7010 (mock) and 0.6189 (experimental), showing that a single optimized machine-learning classifier provides superior resolution while avoiding latency introduced by using multiple software."

Acurácia balanceada (BA) medida:

| Abordagem | BA mock | BA experimental |
|---|---|---|
| **geNomad sozinho (default vOMIX)** | **0,9785** | **0,8702** |
| VIRify | 0,8985 | 0,7845 |
| PPR-Meta | 0,8784 | 0,8259 |
| ViWrap | 0,8380 | 0,8296 |
| ViroProfiler | 0,8380 | 0,8251 |
| **vOMIX multi-tool (consenso por voto)** | **0,7010** | **0,6189** |
| Nayfach et al. | 0,5041 | 0,5377 |

Especificidade no mock: geNomad 0,9826; VirSorter2 0,5728; consenso multi-tool 0,4101; VirFinder 0,2631; **VIBRANT 0,0843**.

A leitura: as ferramentas alternativas têm sensibilidade alta (multi-tool chega a 0,9918) mas especificidade péssima, e um voto por contagem **herda os falsos positivos de todo mundo**. O geNomad sozinho entrega sensibilidade 0,9743 *com* especificidade 0,9826.

> **Isso é um desafio direto ao `VIRAL_CONSENSUS_MODE` da vapor.** Nosso `viral_detection.smk` roda VirSorter2 + geNomad + VIBRANT e combina por `"count"` / `"score"` / `"hybrid"`. Se o resultado do vOMIX se sustentar no nosso tipo de dado, o modo `count` está *degradando* a especificidade em relação a usar só o geNomad. **Não vamos aceitar isso de cara** — o benchmark deles é sobre um mock próprio e um dataset do Wu et al., e o `hybrid` ponderado por score não é o mesmo que o voto simples que eles testaram. Mas é a pergunta empírica mais importante que este material levanta para a vapor, e dá para responder com dados próprios (ver §6).

---

## 4. Como cada um consegue "um ambiente só"

A pergunta era como eles chegam a um único ambiente conda. As duas respostas são diferentes, e nenhuma das duas é "colocaram tudo num env".

### metaFun — conda só para o *launcher*
`meta.yaml` (receita bioconda) declara como dependências apenas: `python`, `apptainer=1.3.0`, `nextflow=24.04.2`, `dash`, `plotly`, `pandas`, `sylph=0.6.1`, `squashfuse`. **Nenhuma ferramenta de bioinformática.**

Todo o resto vive em 5 imagens Apptainer pré-construídas, atribuídas por processo em `config/nextflow.config`:

```groovy
process {
    container = "${projectDir}/../sif_images/metafun_v0.1.sif"   // default
    withName: run_prokka          { container = ".../prokka_latest.sif" }
    withName: sylph_sketch_all    { container = ".../interactive_wms_taxonomy_v02.sif" }
    withName: 'prevalence_filter_phyloseq|concat_genomes|bowtie2_build|prodigal|...' {
        container = ".../instrain_wms_strain_v03.sif" }
}
```

Os `.sif` são baixados sob demanda (`$PREFIX/share/metafun/sif_images/`, criado com `chmod 777` no `build.sh`). Dentro dos containers eles usam `micromamba activate <env>` — as linhas `//conda "$HOME/miniforge3/envs/..."` nos `.nf` são resquício comentado da fase de desenvolvimento.

**Custo:** o usuário não resolve nenhuma dependência, mas depende de imagens que só o autor consegue reconstruir. Reprodutibilidade alta, auditabilidade baixa.

### vOMIX-MEGA — conda por regra, empacotado num `.sif` só
`environment.yml` do env `vomix` tem 7 pacotes: python, pip, click, rich-click, `snakemake>=9.23.1`, biopython, plugin de cluster. Também é só launcher.

A diferença é o que acontece depois. O `workflow/Snakefile` abre com:

```python
containerized: "oras://ghcr.io/erfanshekarriz/vomix:v0.1.0-beta.1"
```

Essa diretiva do Snakemake diz: "quando rodar com `--sdm apptainer`, todas as regras com `conda:` devem usar os envs **já materializados dentro desta imagem**". A imagem é gerada por `snakemake --containerize`, que percorre todos os `envs/*.yml`, resolve cada um e os assa em camadas sobre um SO base. Da doc deles (`docs/install.md`):

> "The container image generated contains explicitly each conda environment mounted on top of a base operating system."

E o aviso operacional, que é a parte fácil de errar:

> "If you only use `--sdm apptainer`, Snakemake will not launch any conda environments and hence all jobs will fail. If you use `--sdm apptainer --use-conda` it will try and re-install conda environments in your local `.snakemake/conda` folder, which counteracts the purpose of containers."

Ou seja: **`--sdm apptainer` sozinho** usa a imagem; `--sdm conda --use-conda` usa conda local. São dois caminhos, o mesmo `envs/*.yml` alimenta os dois. Há ainda `conda-lock.yml` como fallback de resolução, e `workflow/envs/CentOS8/` com variantes de env para um SO específico.

### O que isso significa para a vapor

A vapor **já está no modelo do vOMIX** — 26 `envs/*.yaml` + `container:` por regra. O que falta é uma linha:

```python
containerized: "oras://ghcr.io/<org>/vapor:<tag>"
```

no `Snakefile`, mais um `snakemake --containerize > Dockerfile` no CI. Isso daria a instalação "um comando só" **sem** abandonar os envs conda e **sem** exigir imagens artesanais como as do metaFun. É a mudança de melhor razão esforço/benefício de todo este documento.

> ⚠️ **Correção (2026-08-17).** A frase original aqui dizia que os `container:` por regra da vapor eram "um modelo mais frágil, cada imagem tem que existir e ser mantida". **Isso está errado**, e a recomendação abaixo foi feita sobre essa premissa falsa.
>
> A vapor tem uma estratégia deliberada: `containers.yaml` (47 ferramentas → pacote bioconda + versão) → `pin_containers.py` → `containers.lock.yaml` (URIs `quay.io/biocontainers` pinadas). **As imagens são construídas automaticamente pelo projeto bioconda** — ninguém na vapor mantém 47 imagens, mantém uma lista de versões. Só duas imagens são próprias (`genome-map`, `medaka-gpu`), publicadas pelo CI.
>
> O problema real era que o `containers.lock.yaml` **nunca havia sido gerado**, então `CONTAINERS = {}` e os 47 `container:` resolviam para `None` — a estratégia inteira estava inerte. Resolvido gerando e commitando o lock, mais um job de CI que roda `pin_containers.py --check`.
>
> E os dois mecanismos coexistem mal: `container:` por regra **tem precedência** sobre `containerized:`, então com o lock presente a imagem global seria quase toda ignorada. Adotar `containerized:` exigiria **remover** o esquema atual, não somar a ele.

---

## 5. Os benchmarks do metaFun — e o que eles decidiram

O metaFun é o oposto do vOMIX: quase nenhuma engenharia de performance, mas **escolha de parâmetro por benchmark sistemático** em dados CAMI simulados de três ambientes (gut humano, marinho, rizosfera; n=5 cada).

| Decisão | Alternativas testadas | Escolha | Justificativa medida |
|---|---|---|---|
| Assembler | metaSPAdes vs. MEGAHIT (2 configs) | **MEGAHIT** | Vence em comprimento total e fração de bp alinhados a referência, sobretudo em rizosfera (alta complexidade); metaSPAdes só ganha em misassembly em gut/marinho. Decidiu no conjunto + eficiência de recursos |
| Binning | MetaBAT2, SemiBin2 (pré-treinado), SemiBin2 (self-supervised), DAS Tool×2 | **DAS Tool(MetaBAT2 + SemiBin2 self-supervised)** | Refinamento melhora acurácia de atribuição de bp; modo self-learning do SemiBin2 melhora ARI em metagenomas complexos. Métricas via AMBER |
| Clustering de família gênica | 90% id/50% cov vs. 80% id/80% cov | **80/80** | 90/50 inflou artificialmente o nº de famílias e piorou a acurácia de anotação (Manhattan normalizada + F1) |
| Threshold de core gene | 90%, 95%, 99% | **90%** | Com MAGs fragmentados, 99% praticamente zera o core genome |
| Perfil taxonômico | Kraken2+Bracken (confidence 0–0,5 × filtro 0–0,01%) vs. sylph (c50–c1000) | **os dois** — Kraken2 conf 0,1/0,25 + filtro 0,01%; sylph c200 | sylph venceu em L1, Bray-Curtis e F1 *com* melhor eficiência; mas Kraken2 classifica reads individuais, o que sylph não faz. Mantiveram ambos por capacidade, não por acurácia |

Dois achados metodológicos que valem para nós independentemente de ferramenta:

- **Kraken2 é violentamente sensível a parâmetro:** variando confidence e filtro, as métricas oscilaram **3× em distância composicional e 45× em F1**. sylph mostrou variação mínima entre parâmetros de compressão. Argumento forte a favor da nossa escolha de sylph em `reads_classify.smk`.
- **Fluidez genômica em MAGs incompletos é uma armadilha:** *Bordetella holmesii* saltou de 0,002 para 0,148 no conjunto "Fragmented Incomplete". Interpretação de pangenoma a partir de MAG precisa de ressalva explícita.

---

## 6. Recomendações para a vapor

Ordenadas por razão benefício/esforço. Nada aqui foi implementado — é uma proposta.

**Alto retorno, baixo risco**

1. ~~**`containerized:` no Snakefile**~~ — **descartado**, ver a correção em §4. O equivalente de valor era ativar a estratégia de containers que já existia: `containers.lock.yaml` gerado e commitado (48 URIs BioContainers pinadas) + `pin_containers.py --check` no CI. Feito em 2026-08-17. O `containerized:` continua sendo uma opção válida para quem precise de um artefato único offline, mas exigiria remover `containers.yaml`, `pin_containers.py` e as 47 diretivas `container:`.
2. **CheckV-PyHMMER na `viral_binning.smk`.** Padrão de priming do `tmp/`: `pyrodigal-gv` → 80 jobs `pyhmmer` de 1 thread → `checkv end_to_end -t 1`. Já temos `pyhmmer=0.12.0` e `pyrodigal-gv=0.3.2` em `env_annotation.yaml`. Fazer a versão *com* pyrodigal-gv de verdade, que é o que o paper descreve e o código deles não faz. Ganho esperado só acima de ~10⁵ contigs — abaixo disso, o CheckV nativo é mais rápido, então precisa ser opcional via config. (§3.2)
3. **`hostile` na `host_removal.smk`.** Índices curados e versionados no lugar do nosso mapeamento manual. (§2.2)

**Alto retorno, precisa de decisão nossa**

4. **Reavaliar o `VIRAL_CONSENSUS_MODE`.** Rodar geNomad-sozinho contra os nossos três modos de consenso num mock com verdade conhecida, medindo especificidade — não só sensibilidade. Se o padrão do vOMIX se confirmar, o default deveria mudar. Este é o item mais importante do documento. (§3.4)
5. **PhaTYP e CHERRY** via PhaBOX2. PhaTYP daria estilo de vida sobre os vOTUs montados (hoje só temos BACPHLIP na trilha de reads, sobre genomas de referência); CHERRY complementa o PHIST. (§2.2)

**Considerar depois**

6. **strobealign** como alternativa ao bwa-mem2 para binning de leitura curta.
7. **Árvore binária de clustering** no `votu_catalog.smk` — só se passarmos de ~10⁵ genomas no pool. O `skani triangle --sparse` já é bem melhor que o CD-HIT que eles substituíram, então isto é otimização de segunda ordem para nós.
8. **inStrain** para microdiversidade/strain-level — capacidade nova, não otimização.
9. **`--splits` no geNomad** (`contig-splits`) para reduzir RSS de pico. Barato.

---

## Referências

Lee, H. G., Song, J. Y., Yoon, J., Chung, Y., Kwon, S.-K., & Kim, J. F. (2026). metaFun: An analysis pipeline for metagenomic big data with fast and unified functional searches. *Gut Microbes*, 18(1), 2611544. https://doi.org/10.1080/19490976.2025.2611544

Shekarriz, E., Vijendran, E., & Ho, J. W. K. (2026). vOMIX-MEGA: An ultra-fast end-to-end pipeline for terabyte-scale viral metagenomics analysis. *bioRxiv*. https://doi.org/10.64898/2026.07.28.741255

Camargo, A. P., Roux, S., Schulz, F., Babinski, M., Xu, Y., Hu, B., Chain, P. S. G., Nayfach, S., & Kyrpides, N. C. (2024). Identification of mobile genetic elements with geNomad. *Nature Biotechnology*, 42(8), 1303–1312. https://doi.org/10.1038/s41587-023-01953-y

Larralde, M., & Zeller, G. (2023). PyHMMER: a Python library binding to HMMER for efficient sequence analysis. *Bioinformatics*, 39(5), btad214. https://doi.org/10.1093/bioinformatics/btad214

Nayfach, S., Camargo, A. P., Schulz, F., Eloe-Fadrosh, E., Roux, S., & Kyrpides, N. C. (2021). CheckV assesses the quality and completeness of metagenome-assembled viral genomes. *Nature Biotechnology*, 39(5), 578–585. https://doi.org/10.1038/s41587-020-00774-7

Shaw, J., & Yu, Y. W. (2023). Fast and robust metagenomic sequence comparison through sparse chaining with skani. *Nature Methods*, 20(11), 1661–1665. https://doi.org/10.1038/s41592-023-02018-3

Shaw, J., & Yu, Y. W. (2024). Rapid species-level metagenome profiling and containment estimation with sylph. *Nature Biotechnology*, 43(7), 1122–1131. https://doi.org/10.1038/s41587-024-02412-y
