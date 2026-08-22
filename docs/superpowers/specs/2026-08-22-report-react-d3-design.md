# Report HTML em React + D3 — desenho

**Criado:** 2026-08-22
**Estado:** desenho aprovado, implementação não iniciada
**Motivação:** o report atual conta a história que a pipeline tinha até julho —
centrada na amostra — enquanto a pipeline de agosto passou a girar em torno de
dois catálogos globais. Ao mesmo tempo, dados que já existem em disco há semanas
(pangenoma, catálogo de MAGs, GUNC, bakta, descartes virais, plasmídeos) não
aparecem em lugar nenhum.

---

## 1. O problema

Três problemas distintos, que este desenho resolve de uma vez porque separá-los
custaria dois rewrites.

**(a) O eixo narrativo está errado.** O catálogo global de vOTUs
(`votu_catalog.smk`, 2026-08-18) e o catálogo global de MAGs (`mag_catalog.smk`,
2026-08-19) tornaram-se a espinha dorsal: nada procariótico a jusante do binning
é per-sample, e o vOTU passou a ter um espaço de IDs único para todas as amostras
e grupos. O report ainda organiza tudo por amostra, então ele não mostra a
herança representante→membro — e como toda tabela procariótica a jusante é
*herdada*, quem lê o report não tem como interpretar o que está vendo.

**(b) Há dados invisíveis.** Sem loader e sem gráfico:
`pangenome/` inteiro (`gene_by_member.tsv`, `candidates.tsv`), o catálogo de MAGs
(`provenance.tsv`, `mag_membership.tsv`, clusters do galah), GUNC, bakta, o
sidecar de descarte viral do item (e) (`viral_discarded.tsv`) e a camada de
plasmídeos (`plasmidfinder_results.tsv`) — esta última nunca cruzada com AMR nem
com defesa, que é justamente onde ela informa alguma coisa.

**(c) A camada de rede não existe.** O vConTACT3 saiu em 2026-08-17 e nada
ocupou o lugar: não há hoje nenhum agrupamento de vOTUs por conteúdo gênico, e a
predição de hospedeiro do PHIST é mostrada como tabela, quando ela é uma relação
bipartida.

O que **não** é problema, e foi verificado: o tamanho do HTML. O `report.html`
de 428 MB no diretório do projeto é de 2026-08-12 e antecede a projeção de campos
em `load_gtdbtk` (`data_loaders.py:544`, hoje 12 campos em vez da linha crua do
`classify_wf`). O código atual gera report de tamanho normal. Fica só o
guarda-corpo da §8.

---

## 2. Decisões tomadas (não reabrir sem motivo novo)

| decisão | valor | por quê |
|---|---|---|
| Toolchain | **bundle pré-compilado, versionado no git** | Snakemake continua sem Node em runtime; rodada offline/HPC não quebra. Node só para quem edita o report |
| Framework | **React + D3** | componentes com estado (filtro global de amostra, seletor de rank, drill-down de rede) são o que o vanilla atual não sustenta |
| Biblioteca de charts | **D3 puro, ECharts removida** | um só modelo mental e estética unificada; o custo é uma camada de primitivas, escrita uma vez |
| Navegação | **por objeto biológico, amostra como filtro** | reflete a pipeline pós-agosto; amostra deixa de ser aba |
| Taxonomia | **seletor de rank, nunca rank fixo** | filo/classe/ordem/família/gênero trocáveis no lugar; clique no sunburst desce e filtra |
| Treemap | **proibido** | decisão do usuário; o sunburst com drill cobre o caso de "uma unidade, hierarquia inteira" |
| Redes | **graph-tool, SBM aninhado, layout no Python** | layout determinístico com semente fixa: o mesmo dado desenha igual, o que a figura de paper exige |
| Entrada da rede de genes | **grafo de blocos, não a rede inteira** | acima de ~150 nós node-link não se lê (`REPORT_VIZ_GUIDE.md` §4) |
| Fago–hospedeiro | **matriz ordenada por bloco como forma primária** | node-link fica só na vista ego |
| Plasmídeo–AMR–defesa | **sem grafo** | conjuntos pequenos; UpSet + trilha genômica respondem melhor |
| `REPORT_VIZ_GUIDE.md` | **continua sendo lei** | é ativo, não legado; os gatilhos migram para dentro dos componentes |

---

## 3. Arquitetura

Três camadas com contratos explícitos.

### (a) Produção de dados — Python, no Snakemake

`scripts/report/` segue sendo quem lê o disco, mas muda de saída. Cada domínio
escreve um bloco JSON cujo formato é declarado em `scripts/report/schemas/`. O
schema é a fronteira: **um gráfico só consome campo que o schema declara**, e o
Python projeta apenas esses campos. É a regra que impede a reincidência do
padrão `base = dict(gtdb_bins[key])` — copiar a linha crua da ferramenta e
embarcar trinta campos para consumir sete.

`renderer.py` encolhe para quatro passos: montar o JSON, ler o bundle, ler o CSS,
escrever o HTML. Some o `str.replace` de blocos de JS.

### (b) Pré-cálculo de redes — Python, env novo

`envs/env_network.yaml` (graph-tool + MMseqs2) e `rules/report_network.smk`.
Detalhado na §6.

### (c) Interface — React + D3, `src/report-ui/`

```
src/report-ui/
  viz/        escalas, eixos, legenda, tooltip, useResize, useTheme, paleta
  charts/     um componente por FORMA (não por dataset)
  panels/     um por aba: escolhe formas, passa dados
  state/      filtro global de amostra, rank taxonômico, seleção
```

Compilado por esbuild em `scripts/report/assets/report-ui.js`, versionado.

**Os gatilhos numéricos do guia migram para dentro de `charts/`.** Hoje `VIZ`
(n>12 vira horizontal, n<20 vira strip, >150 nós vira matriz) depende de o autor
do gráfico lembrar de chamar o helper certo. Dentro da forma, o padrão fica
correto sem depender de disciplina.

---

## 4. Contrato de dados

| bloco | fonte em disco | chave |
|---|---|---|
| `run` | `done.txt` de cada regra, `benchmarks/` | regra × amostra |
| `sequencing` | fastp, NanoPlot, QUAST, mapping, depth | amostra |
| `viral_catalog` | `vOTU_clusters.tsv`, `provenance.tsv`, `presence_matrix.tsv`, `votu_abundance_matrix.tsv`, `viral_taxonomy_merged.tsv`, `votu_lifestyle.tsv`, CheckV | vOTU representante (namespaced `{source}\|{contig}`) |
| `viral_discarded` | `final/viral/viral_discarded.tsv` | contig + motivo |
| `mag_catalog` | `mag_membership.tsv`, `provenance.tsv`, CheckM2, GUNC, GTDB-Tk global | representante + membros |
| `annotation` | `module_completeness.tsv`, `cazy_per_mag.tsv`, `ko_per_mag.tsv`, bakta, pharokka/phold | genoma |
| `defense_amr` | `defensefinder_systems.tsv`, `amr_consensus.tsv`, `plasmidfinder_results.tsv`, `vfdb_results.tsv`, ilhas | genoma + contig + coordenada |
| `pangenome` | `gene_by_member.tsv`, `candidates.tsv` | cluster × membro × gene |
| `networks` | `nodes.tsv` / `edges.tsv` das duas redes | nó |
| `diversity` | alfa, PCoA, procrustes | amostra |
| `reads` | sylph merged, host map | amostra |

**Duas armadilhas de ID que o schema tem de tornar explícitas**, porque já
causaram bug silencioso na pipeline:

1. Um ID de proteína do catálogo é `S1__binette_bin1__k141_1_5`. **Cortar no
   primeiro `__` atribui o achado à AMOSTRA em vez do MAG.** A view já reescreve
   o prefixo para o nome original antes de o report ler — o report herda essa
   garantia e não pode refazer o corte por conta própria.
2. O espaço de IDs do sylph (`t__IMGVR_UViG_…`) **não conversa** com o dos
   contigs montados (`k141_…`). Nenhum join entre as duas trilhas é válido, e a
   aba de leituras diz isso na cara do usuário.

---

## 5. Catálogo de gráficos por tipo de dado

Cada item: fonte → forma → a regra do `REPORT_VIZ_GUIDE.md` que a justifica.

### 5.1 Visão geral

- **Fita de KPIs** (amostras, contigs, vOTUs, MAGs, representantes) como stat
  tiles. §2 do guia: valor único não vira gráfico de uma barra.
- **Funil de atrição**: reads → trimados → contigs ≥ `MIN_CONTIG` → candidatos
  virais → pós-CheckV → pós-vRhyme → vOTUs retidos, com o ramo de descarte
  saindo lateral em cada etapa. Clicar na perda abre a repartição por motivo
  vinda de `viral_discarded.tsv` (tier CheckV, comprimento, sem bin). Sem isto,
  quem perdeu 80% do viroma no portão do item (e) não descobre.
- **Matriz de status** (amostra × regra) dos `done.txt`, com `ok`/`skipped`/
  `failed` visualmente distintos. **Ferramenta que falhou é lacuna, nunca zero
  biológico** — o bug do AMRFinderPlus com disco cheio.

### 5.2 Sequenciamento e montagem

fastp antes/depois em barra empilhada 100% + Q30; QUAST em heatmap amostra ×
métrica normalizado **por métrica**; comprimento de contig em ridgeline log;
taxa de mapeamento em lollipop; cobertura em hexbin comprimento × profundidade.
Boxplot de comprimento é proibido (§7): montagem metagenômica é bimodal e a
caixa esconde exatamente isso.

### 5.3 Catálogo viral

- **Composição taxonômica** com **seletor de rank** (filo→classe→ordem→família→
  gênero). Barra horizontal ordenada + cauda em "Other"; nunca pizza, nunca uma
  nona cor. Sunburst para *uma* unidade com clique que desce o rank e filtra a
  aba; barra empilhada no rank selecionado para *comparar* amostras — hierarquia
  não se lê lado a lado.
- **Qualidade CheckV** como escada ordinal Complete→HQ→MQ→LQ→ND, rampa boa→ruim.
- **vOTU × amostra**: heatmap clusterizado com dendrograma; bolha quando esparso
  (<20% das células). Em vOTU o esparso é o caso comum — o componente decide.
- **Curva de acumulação global.** O guia dizia que acumulação só valia na
  co-montagem porque os vOTUs per-sample eram clusterizados independentemente.
  **O catálogo global de 2026-08-18 tornou isso obsoleto**: há espaço de features
  compartilhado entre todas as amostras. Média sobre ordens aleatórias, banda
  10–90. *Esta linha do guia precisa ser corrigida junto com a implementação.*
- **Estilo de vida** (BACPHLIP) × qualidade, com aviso de que BACPHLIP só é
  confiável em genoma completo.
- **Concordância de detectores** (VirSorter2 × geNomad × consenso): **UpSet**.
- **Explorador de vOTU**: trilha genômica em bp real — ORFs do prodigal-gv com
  fita, PHROGs coloridos por categoria, AMGs, anti-defesa e hits dbAPIS na
  coordenada. Ordem de gene desenhada como mapa é anti-padrão; a coordenada vem
  de graça no cabeçalho do Prodigal.

### 5.4 Catálogo de MAGs

- **Qualidade**: scatter completude × contaminação com as zonas MIMAG
  desenhadas. As linhas de corte são a razão de o gráfico existir.
- **GUNC** (inédito): CSS por MAG, e o cruzamento que importa — scatter
  contaminação × CSS destacando o quadrante "passa no CheckM2, reprova no GUNC".
- **Estrutura do catálogo** (inédito): tamanho dos clusters do galah em barra,
  mais proveniência por amostra/grupo. Como toda anotação procariótica é herdada
  do representante, a herança precisa estar visível para qualquer tabela a
  jusante ser interpretável.
- **Taxonomia GTDB** com o mesmo seletor de rank da §5.3.
- **Metabolismo**: heatmap MAG × módulo KEGG com `missing_ko` no tooltip — é o
  que torna um módulo incompleto interpretável. CAZy por classe (GH/GT/PL/CE/
  AA/CBM) em barra empilhada.

### 5.5 Defesa, AMR e plasmídeos

- **Sistemas de defesa**: heatmap MAG × tipo, linhas ordenadas pela taxonomia
  GTDB, de modo que "este clado carrega CBASS" seja legível.
- **ARGs**: barra por classe de droga a partir do consenso (`n_tools ≥ 2`), com
  o número de ferramentas concordantes como encoding secundário — nunca só cor.
- **A ligação plasmídeo–AMR–defesa.** O ABRicate dá `SEQUENCE/START/END` do
  replicon; AMRFinder e DefenseFinder dão IDs de proteína cujas coordenadas saem
  do cabeçalho do Prodigal. Isso permite **colocalização no mesmo contig**, que é
  uma afirmação biológica de verdade: *este ARG está num contig que carrega
  replicon plasmidial*. Três formas: (1) UpSet de MAGs por {replicon, ARG de
  consenso, ilha de defesa}; (2) trilha genômica do contig com replicon, ARG e
  sistema nas coordenadas reais; (3) stat tile da fração de ARGs em contig com
  replicon. **Ressalva escrita na aba:** replicon em contig de MAG é evidência de
  origem plasmidial, não prova de plasmídeo intacto — a montagem pode ter
  quebrado o elemento.
- **Ilhas de defesa**: trilha com a janela destacada + ranking por densidade.

### 5.6 Pangenoma

- **Matriz gene × membro em TRÊS estados.** `x` presente, `.` ausente, `?` não
  avaliável (membro abaixo de 70% de completude). **O `?` não pode parecer
  ausência** — é a distinção que a regra do Snakemake se deu ao trabalho de
  codificar, e uma escala sequencial de duas cores a destruiria. Categórica de
  três estados, `?` hachurado, denominador de frequência excluindo `?` declarado
  no tooltip.
- **Core/shell/cloud** por cluster; barra empilhada comparando clusters.
- **Candidatos**: tabela de `candidates.tsv` com o critério que qualificou cada
  cluster. O sinal do PlasmidFinder aparece como sinal de mobilidade, **não como
  critério** — na regra ele não elege cluster sozinho.

### 5.7 Diversidade e leituras

Alfa em strip plot com poucas amostras, ridgeline com muitas (§4). PCoA com **%
de variância nos dois eixos** e elipses quando há grupos. Procrustes viral ×
procariótico como setas ligando pares. **Simpson e Chao1 só aparecem quando
vieram de contagens**; escritos vazios por falta de contagem, a aba mostra a
lacuna em vez de plotar zero. Trilha sylph: OTU por prevalência, hospedeiro
colapsado com `host_source` (`db`/`none`) visível, e o aviso de espaço de IDs da
§4.

---

## 6. A camada de redes

### 6.1 Rede de compartilhamento gênico entre vOTUs

**Construção.** MMseqs2 `easy-cluster` sobre o `.faa` global do `votu_prodigal`
→ clusters de proteína. Cada vOTU vira um vetor de clusters; a aresta entre dois
vOTUs é a significância do compartilhamento por teste hipergeométrico, no
espírito do gene-sharing network do vConTACT (Bin Jang et al., 2019). O corte
fica em `config.yaml`, não escondido no script. Sobre a rede roda o **SBM
aninhado** (Peixoto, 2014), que devolve hierarquia de blocos.

**Como se entra na rede.** Não pela rede inteira — milhares de nós, e acima de
~150 node-link não se lê. A vista padrão é o **grafo de blocos**: nó = bloco,
tamanho = número de vOTUs, cor = família majoritária do MMseqs2, espessura =
compartilhamento entre blocos. Dezenas de nós, legível. Clicar desce um nível da
hierarquia; no fundo chega-se ao node-link dos vOTUs do bloco, aí sim <150. A
rede completa fica atrás de um botão, em canvas, com o aviso de que ali se vê
estrutura global e não identidade individual.

**Por que isso é ciência e não enfeite.** O bloco do SBM é um agrupamento tipo-
gênero derivado de conteúdo gênico, independente da atribuição por homologia.
Pintar os blocos com as famílias do MMseqs2 mostra **onde os dois discordam** — e
essa discordância é onde mora vírus novo.

### 6.2 Rede bipartida fago–hospedeiro

SBM bipartido sobre o resultado do PHIST — o graph-tool infere a bipartição como
estrutura em vez de precisar ser informado dela. Forma primária: **matriz de
adjacência ordenada pelos blocos**, vOTUs nas linhas, MAGs nas colunas, filos do
GTDB anotados na margem. Blocos diagonais densos = faixa estreita de hospedeiro;
linhas espalhadas = generalista. Node-link só como **vista ego** (um vOTU e seus
hospedeiros, ou um MAG e seus fagos, com pivô entre os lados). O `adj-pvalue`
precisa estar codificado na aresta: PHIST dá evidência, não certeza.

### 6.3 Contrato Python → browser

`nodes.tsv`: `id, tipo, bloco_l0..bloco_lN, x, y, grau, rotulo`.
`edges.tsv`: `origem, destino, peso`.

**O layout SFDP sai do Python com semente fixa.** O D3 não simula força: ele
desenha, filtra e interage. Rede desenhada por simulação muda a cada abertura, o
que é inaceitável numa figura que vai para o paper.

Trilha desligada ou graph-tool ausente → os arquivos não existem → a aba mostra
estado vazio. Mesma disciplina que o resto do report já segue.

---

## 7. Regras Snakemake novas

| regra | entrada | saída |
|---|---|---|
| `report_gene_network` | `.faa` global do `votu_prodigal` | `report/networks/gene_sharing/{nodes,edges}.tsv` |
| `report_host_network` | `phist_results.csv` do catálogo + GTDB global | `report/networks/phage_host/{nodes,edges}.tsv` |

Ambas em `rules/report_network.smk`, env `env_network.yaml`, contador de jobs
**fixo** — duas, para a rodada inteira. Nenhum `checkpoint`: DAG dinâmico
quebraria o invariante de dry-run de que a verificação do
`ROADMAP_SIMPLIFICACAO.md` depende.

`rule report` ganha as duas saídas como entrada opcional.

---

## 8. Verificação

1. **Testes de loader/schema** — cada bloco JSON valida contra seu schema; campo
   não declarado é erro, não aviso.
2. **Orçamento de payload** — teste falha se o bundle de dados passar de 25 MB
   num dataset de referência. Guarda-corpo contra a reincidência do padrão
   descrito no fim da §1 — embarcar a linha crua da ferramenta.
3. **Determinismo da rede** — rodar duas vezes o pré-cálculo com a mesma semente
   produz `nodes.tsv` byte a byte idêntico.
4. **Degradação** — com cada trilha desligada em `config.yaml`, o report gera e
   abre; nenhuma aba com eixo quebrado.
5. **Dry-run** — `snakemake -n` continua com contagem de jobs previsível.
6. **Paleta** — `validate_palette.js` em modo claro e escuro antes de qualquer
   mudança de cor. ΔE não se avalia a olho.

---

## 9. Migração

O report atual continua funcionando enquanto o novo é construído: `rule report`
ganha um alvo paralelo, e a troca é uma linha no `Snakefile` quando a paridade
for atingida. `scripts/report/components/*.js` e `assets/echarts.min.js` só saem
do repositório depois disso, num commit separado, para que o diff da remoção não
se misture ao da construção.

---

## 10. Fora de escopo

- **PPanGGOLiN** (fase 2 do pangenoma) — o report lê o que a fase 1 produz.
- **Servir o report** — segue arquivo único, sem servidor.
- **Rede plasmídeo–AMR–defesa como grafo** — decidido contra: conjuntos pequenos,
  pergunta de coocorrência e coordenada, respondida melhor por UpSet e trilha.
- **Reprocessar dados** — o report não computa biologia; toda análise nova vira
  regra, e as duas da §7 são as únicas deste desenho.

---

## Referências

Bin Jang, H., Bolduc, B., Zablocki, O., Kuhn, J. H., Roux, S., Adriaenssens, E. M., Brister, J. R., Kropinski, A. M., Krupovic, M., Lavigne, R., Turner, D., & Sullivan, M. B. (2019). Taxonomic assignment of uncultivated prokaryotic virus genomes is enabled by gene-sharing networks. *Nature Biotechnology*, 37(6), 632–639. https://doi.org/10.1038/s41587-019-0100-8

Bowers, R. M., Kyrpides, N. C., Stepanauskas, R., Harmon-Smith, M., Doud, D., Reddy, T. B. K., Schulz, F., Jarett, J., Rivers, A. R., Eloe-Fadrosh, E. A., Tringe, S. G., Ivanova, N. N., Copeland, A., Clum, A., Becraft, E. D., Malmstrom, R. R., Birren, B., Podar, M., Bork, P., … Woyke, T. (2017). Minimum information about a single amplified genome (MISAG) and a metagenome-assembled genome (MIMAG) of bacteria and archaea. *Nature Biotechnology*, 35(8), 725–731. https://doi.org/10.1038/nbt.3893

Hockenberry, A. J., & Wilke, C. O. (2021). BACPHLIP: predicting bacteriophage lifestyle from conserved protein domains. *PeerJ*, 9, e11396. https://doi.org/10.7717/peerj.11396

Lex, A., Gehlenborg, N., Strobelt, H., Vuillemot, R., & Pfister, H. (2014). UpSet: Visualization of intersecting sets. *IEEE Transactions on Visualization and Computer Graphics*, 20(12), 1983–1992. https://doi.org/10.1109/TVCG.2014.2346248

Nayfach, S., Camargo, A. P., Schulz, F., Eloe-Fadrosh, E., Roux, S., & Kyrpides, N. C. (2021). CheckV assesses the quality and completeness of metagenome-assembled viral genomes. *Nature Biotechnology*, 39(5), 578–585. https://doi.org/10.1038/s41587-020-00774-7

Orakov, A., Fullam, A., Coelho, L. P., Khedkar, S., Szklarczyk, D., Mende, D. R., Schmidt, T. S. B., & Bork, P. (2021). GUNC: detection of chimerism and contamination in prokaryotic genomes. *Genome Biology*, 22, 178. https://doi.org/10.1186/s13059-021-02393-0

Peixoto, T. P. (2014). Hierarchical block structures and high-resolution model selection in large networks. *Physical Review X*, 4(1), 011047. https://doi.org/10.1103/PhysRevX.4.011047

Steinegger, M., & Söding, J. (2017). MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. *Nature Biotechnology*, 35(11), 1026–1028. https://doi.org/10.1038/nbt.3988

Zielezinski, A., Deorowicz, S., & Gudyś, A. (2022). PHIST: fast and accurate prediction of prokaryotic hosts from metagenomic viral sequences. *Bioinformatics*, 38(5), 1447–1449. https://doi.org/10.1093/bioinformatics/btab837
