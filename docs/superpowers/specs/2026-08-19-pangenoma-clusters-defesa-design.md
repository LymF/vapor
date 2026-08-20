# Pangenoma dos clusters com ilha de defesa/AMR — desenho

**Criado:** 2026-08-19
**Estado:** desenho aprovado, implementação não iniciada
**Motivação:** responder *"este sistema de defesa (ou ARG) é fixo na espécie ou
foi adquirido por esta linhagem?"* — a pergunta natural para uma ferramenta cujo
foco é defesa/antidefesa/AMR, e que a vapor hoje não consegue nem formular.

---

## 1. O problema

Sistemas de defesa e ARGs vivem no genoma **acessório**: ilhas de defesa,
elementos móveis, plasmídeos. Saber que um MAG tem CBASS é menos informativo que
saber se o CBASS está em todas as linhagens daquela espécie ou só nesta.

O metaFun responde isso com PPanGGOLiN como hub (`ANALISE_TOOLS_VOMIX_METAFUN.md`
§2.6). A vapor não, e o §3.3 daquele documento descartou o pangenoma por escopo.
**Este desenho reabre a decisão em escopo reduzido**, porque a pergunta é
central ao foco da ferramenta — não periférica como o §3.3 supunha.

### O obstáculo que o princípio (h) cria

`docs/ROADMAP_SIMPLIFICACAO.md` item (h): computa-se na representante, herda-se
no membro. Uma análise de core/acessório precisa comparar o conteúdo gênico dos
membros **entre si**, e os membros não têm anotação própria — têm uma cópia da
anotação da representante. Sob o desenho atual, todo cluster do galah pareceria
ter genoma acessório **vazio**.

**A saída, decidida com o usuário:** reintroduzir anotação por membro apenas nos
clusters que passam um portão de interesse. Não desfaz a economia do (h) — que
eliminou ~300 jobs — porque incide sobre dezenas de genomas, não sobre todos.

---

## 2. Decisões tomadas (não reabrir sem motivo novo)

| decisão | valor | por quê |
|---|---|---|
| Mínimo de membros (fase 1) | **3** | com 2 não há frequência que signifique nada |
| Escopo dos genes | **defesa + AMR** | ARG é o exemplo canônico de gene acessório; a máquina é a mesma |
| Faseamento | **fase 1 mede, fase 2 depende** | o nº de clusters elegíveis é desconhecido; construir a fase 2 antes de medir é construir no escuro |
| Ativação | **sem flag de config** | roda quando houver cluster elegível, não roda quando não houver — decisão dos dados, não do usuário |
| Limiar de core gene | **90%, não 99%** | com MAG o padrão zera o core (metaFun) |
| `K` do PPanGGOLiN | **3, fixo** | recomendação explícita da doc do PPanGGOLiN para MAGs |

### Sobre "sem flag": por que NÃO é um checkpoint

O usuário pediu que rode "se tiver MAGs com ilha — nem sempre é, mas geralmente
é". O reflexo seria um `checkpoint` do Snakemake (DAG dependente de dados). **Não
é o caminho**, por duas razões:

1. **O repo não tem nenhum checkpoint** (`grep -rn "^checkpoint " rules/` →
   vazio). Introduzir o padrão para este caso é custo desproporcional.
2. **Quebraria o método de verificação do roadmap.** A invariante usada em todo
   o trabalho de 2026-08 é o diff do conjunto de outputs do dry-run
   (`ROADMAP_SIMPLIFICACAO.md` "Como verificar cada passo"). Um DAG que só se
   resolve depois de um job executar torna essa contagem não-determinística
   antes da execução — perderíamos a ferramenta que pegou os últimos dez bugs.

**O padrão correto já existe no repo: `mag_bakta`.** Ela lê
`qualifying_bins.txt`, itera dentro de UM job e grava um sumário com status por
genoma. Seleção dependente de dados, DAG estático. As regras abaixo seguem essa
forma: número fixo de jobs, cada um decidindo internamente sobre quantos
clusters agir. Zero clusters elegíveis → `skipped: no eligible clusters` no
`done.txt`, que é exatamente a convenção de status honesto que `load_tool_status`
consome.

---

## 3. Fase 1 — matriz gene × membro

### 3.1 Regras

| regra | entrada | saída |
|---|---|---|
| `mag_pangenome_select` | `mag_membership.tsv`, tabelas de defesa/AMR das representantes, `checkm2_quality_report.tsv` | `pangenome/candidates.tsv`, `pangenome/members.txt` |
| `mag_pangenome_proteins` | `members.txt`, `mag_catalog/genomes/*.fa` | `pangenome/proteins/manifest.txt` |
| `mag_pangenome_defensefinder` | manifesto | `pangenome/defensefinder/*.tsv` |
| `mag_pangenome_amr` | manifesto | `pangenome/amr/*.tsv` |
| `mag_pangenome_matrix` | as duas acima + `candidates.tsv` | `pangenome/gene_by_member.tsv`, `pangenome/cluster_summary.tsv` |

**`mag_pangenome_proteins` é uma cópia paramétrica de `mag_catalog_proteins`**
(`rules/mag_catalog.smk:504`), que já é genérica: prodigal sobre um diretório de
genomas, emitindo o manifesto `name / mode / fna / faa / gff`. Muda o diretório
de entrada (`genomes/` filtrado por `members.txt`, em vez de `representatives/`).

**`mag_pangenome_defensefinder` e `mag_pangenome_amr` são
`use rule mag_* as ... with:`** apontando para o manifesto novo. Os corpos
consomem manifesto e não conhecem `{sample}` — foi exatamente para isso que o
formato do manifesto foi preservado na migração de 2026-08-19.

### 3.2 O portão (`mag_pangenome_select`)

```
elegível(cluster) :=
      n_membros >= 3
  AND ( tem_ilha_de_defesa(representante)
     OR n_sistemas_defesa(representante) >= 3
     OR n_args_consenso(representante) >= 1        # n_tools >= 2
     OR tem_hit_plasmidfinder(representante) )
```

PlasmidFinder entra como **sinal de mobilidade que reforça**, nunca como
critério isolado — plasmídeo sem defesa nem ARG não motiva um pangenoma.

`candidates.tsv` registra a evidência por cluster (`n_membros`, `n_ilhas`,
`n_sistemas`, `n_args`, `n_plasmid`, `criterio_disparado`), para que a seleção
seja auditável e não um número que apareceu.

### 3.3 Dívida técnica que isto paga: `scripts/defense_islands.py`

`compute_defense_islands` vive hoje em `scripts/report/data_loaders.py:959` — na
camada de **relatório**. A pipeline não alcança. Para gatear por ilha:

- extrair o núcleo para `scripts/defense_islands.py` (função pura, testável);
- `data_loaders.py` passa a importar de lá;
- `mag_pangenome_select` importa a mesma função.

Hoje há uma definição de ilha só, mas em lugar que a pipeline não vê. Sem a
extração passariam a existir **duas**, e divergiriam em silêncio — a família de
defeito que o roadmap inteiro persegue.

### 3.4 Ausência em MAG não é ausência no organismo

**O ponto científico central deste desenho.** Um `.` na matriz pode significar
"o organismo não tem o gene" ou "o MAG está 74% completo e a região não montou".
Tratar os dois como o mesmo zero produz conclusão errada com aparência de dado —
mesma família do `done.txt` vazio lido como zero biológico.

Duas defesas, **ambas obrigatórias**:

1. **Três estados, nunca dois.** `x` presente; `.` ausente; **`?` não
   avaliável** — membro com completude CheckM2 < 70% (o mesmo piso do
   `mag_bakta`). Um membro `?` não entra no denominador da frequência.
2. **A completude viaja junto.** `gene_by_member.tsv` traz a completude de cada
   membro no cabeçalho, e `cluster_summary.tsv` traz a mediana do cluster. `4/6`
   só é interpretável com o denominador à vista.

O metaFun mediu este fenômeno: fluidez genômica de *Bordetella holmesii* saltou
de **0,002 → 0,148** no conjunto "Fragmented Incomplete". É o mesmo efeito.

### 3.5 O que a fase 1 mede (e que decide a fase 2)

`cluster_summary.tsv` traz, por cluster: `n_membros`, `n_membros_avaliaveis`,
distribuição de tamanho, `n_genes_core` (em ≥ 90% dos avaliáveis),
`n_genes_variaveis`, `n_genes_singleton`, e a taxonomia GTDB da representante.

**Critério de decisão da fase 2, definido agora para não ser racionalizado
depois:** a fase 2 só se justifica se existir **pelo menos um cluster com ≥ 5
membros avaliáveis** (piso do paradigma core/acessório) **e com variação real**
(`n_genes_variaveis > 0`). Sem isso, a matriz da fase 1 é a resposta honesta e a
fase 2 não deve ser construída.

---

## 4. Fase 2 — PPanGGOLiN (condicional)

**Não implementar antes de a fase 1 rodar nos dados reais.**

- Ferramenta: `ppanggolin` 2.3.1 (bioconda, verificado).
- Entrada: bakta por membro (o PPanGGOLiN quer GFF **com sequência**), depois
  `ppanggolin annotate --anno <lista.gff>` → `cluster` → `graph` → `partition`.
- **`-K 3` fixo** e core a **90%**.
- Saída: `pangenome/<cluster>/partition_map.tsv` — cada sistema de defesa e cada
  ARG com sua partição (persistente/shell/cloud) e a frequência entre membros.

### Ressalva de validade, a ser repetida em qualquer figura

Com 3–6 membros o PPanGGOLiN **roda**, mas a separação shell/cloud não tem
sustentação estatística: a recomendação para resultado robusto é **≥ 15 genomas
com variação genômica real** (não só SNPs), e **5** é o piso do paradigma
core/acessório. Ruído com nome de resultado é pior que ausência de resultado.

Estimativa do porte: 106 bins em 5 amostras de `~/global/results` (~21/amostra)
extrapolam para ~700 MAGs nas 32 amostras da Amazônia. Depois de agrupar a 95%
ANI, **é plausível que nenhum cluster alcance 15 membros**. Daí o faseamento.

---

## 5. Resolução de ID — onde isto pode quebrar em silêncio

O repo tem um histórico longo de bugs de namespace (`ROADMAP_SIMPLIFICACAO.md`,
oito manifestações). Os pontos de risco aqui:

1. **Nome de membro.** `mag_catalog/genomes/` usa `{source}__{bin}.fa`, e o bin
   pode ser `binette_bin1` (Binette) ou um inteiro nu (VAMB). **Nunca cortar no
   primeiro `__`**: usar `resolve_prefixed_id` contra os membros conhecidos, como
   as vistas fazem.
2. **Proteína de membro.** O prodigal emite `{contig}_{n}`; o contig do membro
   não é o do representante. A matriz casa gene por **nome de sistema/ARG**, não
   por ID de proteína entre membros diferentes.
3. **AMR é por gene, não por sistema.** DefenseFinder agrupa proteínas num
   sistema nomeado; AMRFinderPlus/RGI/DeepARG devolvem hits por proteína com
   nomenclatura divergente. A matriz usa o nome **normalizado pelo argNorm** e o
   consenso `n_tools >= 2` — senão a coluna "variável" enche de falso positivo
   de uma ferramenta só.

---

## 6. Verificação

Pelo método do roadmap, sem exceção:

- **DAG:** dry-run `-n --forceall` com `config_amazon_18-08-26.yaml` e outdir
  temporário, antes e depois. Delta esperado: **+9 jobs globais fixos**, zero
  por amostra. Qualquer job por amostra é bug de desenho. O AMR sozinho conta
  cinco regras encadeadas (`mag_pangenome_amrfinderplus` →
  `mag_pangenome_rgi_card` → `mag_pangenome_deeparg` → `mag_pangenome_argnorm`
  → `mag_pangenome_amr_consensus`), não uma — cada ferramenta e cada etapa de
  normalização/consenso é um job próprio, como no catálogo global de MAGs.
  Somadas a `mag_pangenome_select`, `mag_pangenome_proteins`,
  `mag_pangenome_defensefinder` e `mag_pangenome_matrix`, o total medido foi
  1404 → 1413.
- **Testes:** funções puras de `scripts/defense_islands.py` (extração) e do
  seletor/matriz. Piso: os 171 atuais continuam passando.
- **Formato:** toda suposição conferida contra arquivo real em
  `~/global/results`, nunca contra docstring — inclusive o `?` de completude,
  que exige o `checkm2_quality_report.tsv` do catálogo.
- **Status:** `mag_pangenome_*` em `STATUS_TRACKED_GLOBAL_TOOLS`. "Nenhum
  cluster elegível" é resultado legítimo e tem de aparecer como `skipped:`,
  jamais como painel vazio.

---

## 7. Fora de escopo

- **Scoary2 / GWAS microbiano** — exige fenótipo por genoma, que a vapor não tem.
- **A trilha `multisplit`** — fora do catálogo por decisão anterior.
- **Pangenoma de todos os clusters** — é o custo que o (h) eliminou.
- **Anotação funcional completa dos membros** (eggNOG, KEGG, CAZy) — os membros
  recebem só o necessário ao portão: proteínas, defesa e AMR.

---

## Referências

Gautreau, G., Bazin, A., Gachet, M., Planel, R., Burlot, L., Dubois, M., Perrin, A., Médigue, C., Calteau, A., Cruveiller, S., Matias, C., Ambroise, C., Rocha, E. P. C., & Vallenet, D. (2020). PPanGGOLiN: Depicting microbial diversity via a partitioned pangenome graph. *PLOS Computational Biology*, 16(3), e1007732. https://doi.org/10.1371/journal.pcbi.1007732

Lee, H. G., Song, J. Y., Yoon, J., Chung, Y., Kwon, S.-K., & Kim, J. F. (2026). metaFun: An analysis pipeline for metagenomic big data with fast and unified functional searches. *Gut Microbes*, 18(1), 2611544. https://doi.org/10.1080/19490976.2025.2611544

### Documentos irmãos neste repo

- `docs/ROADMAP_SIMPLIFICACAO.md` — item (h) e as oito manifestações de bug de namespace
- `docs/ANALISE_TOOLS_VOMIX_METAFUN.md` — §2.6 (módulo comparativo do metaFun) e §3.3 (a decisão que este desenho reabre)
- `CLAUDE.md` — contrato das vistas e do catálogo de MAGs
