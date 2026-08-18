# Roadmap — simplificação da vapor

**Criado:** 2026-08-18
**Motivação:** a vapor acumulou ferramentas redundantes. O objetivo desta linha de
trabalho é **remover complexidade**, não adicionar, e só adicionar onde houver
lacuna funcional real.

Este documento existe para uma sessão nova conseguir continuar sem redescobrir as
decisões. Leia as seções "Decisões já tomadas" e "Achados que não devem ser
redescobertos" ANTES de propor alternativas — várias já foram avaliadas e
descartadas com motivo.

---

## Estado atual

| item | estado |
|---|---|
| status real no `done.txt` (15 regras) | **feito** — commit `4ff4507`, já no master |
| Remoção do VIBRANT | **feito** — commit `d8dbf7d` |
| (f) `bacphlip_votu` + (g) `eggnog_viral` | **feito** — commit `bf2d358`, regras globais |
| Relatório religado (Lifestyle + Putative AMGs) | **feito** — commit `1c66c6c` |
| (b) metaSPAdes + metaviralSPAdes | **feito** — só o MEGAHIT nas short reads |
| (c) um montador só para long reads | **feito** — Flye+Medaka (ONT) / metaMDBG (HiFi) |
| Restante | **próximo passo: (h)** — ver "Ordem de execução" |

Branch de trabalho: `master` (o `refactor/unit-wildcard` já foi mergeado por
fast-forward e pode ser apagado).

---

## Sequência acordada

**As letras são identidade, não ordem.** Elas foram atribuídas na ordem em que
os itens surgiram na discussão; a ordem de execução foi revisada depois e está
abaixo. Não execute em ordem alfabética.

### Ordem de execução

```
(a) VIBRANT ........................ FEITO  d8dbf7d
(f) bacphlip_votu .................. FEITO  bf2d358
(g) eggnog_viral ................... FEITO  bf2d358
        ↓
(b) remover metaSPAdes + metaviralSPAdes ... FEITO
(c) um montador só para long reads ....... FEITO
        ↓
(h) migrar regras para o representante   ← ANTES de (d)
        ↓
(d) remover merge_contigs / merge_lr / mmseqs2
(e) mover o limiar de comprimento
```

**Por que (h) vem antes de (d).** As duas mexem em *quem consome o quê*, e
incidem sobre o mesmo conjunto de pontos de referência — as 22 referências ao
`rep_seq.fasta` e as 7 regras que consomem `votu_catalog_reps` se sobrepõem.
Fazer as duas em sequência sobre os mesmos consumidores dobra a chance de erro
silencioso, que é exatamente a classe de falha que já custou caro aqui (ver
"Lição registrada sobre método").

Fazendo (h) primeiro, as regras que rodam no representante global **saem do
caminho** do hub `rep_seq.fasta` antes de ele ser desmontado. O (d) fica com um
conjunto menor e mais homogêneo de consumidores para reapontar.

(b) e (c) vêm antes das duas porque são pré-requisito de (d): o dedup só pode
sair depois que sobrar **um montador por trilha**, que é o que remove a razão de
existir dele.

Cada passo é verificado antes do seguinte (ver "Como verificar").

### (a) Remover VIBRANT — FEITO (commit `d8dbf7d`)

Motivo: parte relevante do falso positivo viral vem dele.

Superfície: `rules/viral_detection.smk` (regra + parsing no `viral_consensus`),
`rules/coassembly.smk` (gêmea), `rules/viral_binning.smk` (`make_votu_table`),
`scripts/make_votu_table.py`, `scripts/report/{data_loaders,renderer}.py`,
`scripts/report/components/{overview,viral,annotation,about}.js`, `Snakefile`
(`_VIBRANT_BASE` + alvo no `rule all`), `config.yaml`, `vapor.py`,
`envs/phage_vibrant.yaml` (apagar), `envs/env_viral.yaml`, `containers.yaml`,
`containers.lock.yaml`, docs.

**Perdas funcionais aceitas** (o VIBRANT era fonte única das duas):
- detecção de AMG → reposta em (g)
- lifestyle lítico/lisogênico dos vOTUs montados → reposta em (f)

Decisão de schema: as colunas de lifestyle e AMG do `make_votu_table.py`
**continuam existindo, preenchidas como vazio**, para não quebrar quem consome o
arquivo e para religar a fonte depois sem migração.

### (b) Remover metaSPAdes e metaviralSPAdes — FEITO

Fica só o MEGAHIT nas short reads.

Consequência: `rule cobra_spades` morre junto (o COBRA tem uma regra por
montador). `rules/quast.smk` e o relatório referenciam os montadores.

Precedente: o metaFun testou MEGAHIT vs metaSPAdes em 3 ambientes CAMI e ficou
só com o MEGAHIT — venceu em comprimento total e fração de bp alinhados, com
desempate por eficiência de recursos.

**O que saiu:** regras `metaspades`, `metaviral_spades` e `cobra_spades`; chaves
`use_spades`, `spades_mem`, `spades_kmers`, `spades_kmer_list`,
`cobra_spades_mink/maxk`; `spades=4.2.0` do `env_assembly`; entrada `spades` do
`containers.yaml`/`.lock.yaml`; rótulos do QUAST (5 → 3 conjuntos) e do relatório.

**O que ficou de propósito:** o `merge_contigs` continua existindo e continua
prefixando `MEGAHIT_` — ele só morre em (d), e manter o prefixo evita renomear
contig em todo o resultado já produzido. O `cobra_merge` também ficou, agora com
um diretório de corrida só.

**Verificação:** dry-run PE 195 → 189 jobs, 331 → 323 outputs. Delta = 2 amostras
× 4 arquivos (metaspades `contigs.fasta` + `.gfa`, metaviral `contigs.fasta`,
`cobra/spades/done.txt`) e 2 × 3 regras. Nada mais mudou. Dry-run também limpo
em SE (189 jobs) e long reads (167 jobs).

### (c) Um montador só para long reads — FEITO

Escolha por `lr_tech`: **Flye para ONT, metaMDBG para HiFi**. O `hifiasm_lr` sai
sempre; o `merge_lr` morre.

**O que saiu:** regras `hifiasm_lr` e `merge_lr`; o ramo hifiasm do `medaka_lr`
(output `hifiasm_pol`); o guard `USE_METAMDBG` e o ramo `--in-ont` do
`metaMDBG_lr`; constantes `LR_HIFIASM_HOM`, `LR_METAMDBG` e `USE_METAMDBG` do
`Snakefile`; chaves `lr_hifiasm_hom` e `lr_metaMDBG` do `config.yaml`;
`hifiasm`/`hifiasm_meta` do `env_flye`; entrada `hifiasm` dos containers;
`scripts/merge_lr_assemblies.py` (apagado, sem consumidor).

**Estrutura nova:** dentro de `if LONG_READS:` as regras passam a ser
instanciadas por `if LR_TECH == "ont":` / `else:`, então nunca existe montador
no DAG que não vai ser usado. O `mmseqs2` consome:
- ONT → `{sample}/assembly/lr/flye_polished/assembly.fasta`
- HiFi → `{sample}/assembly/lr/metaMDBG/assembly.fasta`

**Prefixo de header:** o `merge_lr` prefixava `FLYE_`/`HIFIASM_`/`MDBG_`. Com um
montador só o prefixo some para long reads. Verificado por grep que **ninguém
faz parsing dele** em `scripts/`, `rules/` ou `scripts/report/` — e o
`votu_catalog_pool` já prefixa por origem no catálogo global.

**Alvo do `rule all` reapontado:** o `Snakefile` ainda pedia
`assembly/lr/merged/done.txt` (o sentinel do `merge_lr`). Passou a pedir o
sentinel do montador único — `medaka_done.txt` (ONT) ou `metaMDBG/done.txt`
(HiFi). Esse era o único ponto que o grep por "merge_lr" não pegava, porque a
referência era um caminho literal.

**Verificação:** dry-run limpo nos quatro modos.

| modo | jobs antes | jobs depois |
|---|---|---|
| PE  | 189 | 189 (diff de outputs **vazio**) |
| SE  | 189 | 189 |
| LR ONT | 167 | 161 |
| LR HiFi | — | 159 |

Delta LR-ONT = 6 = 2 amostras × 3 regras (`hifiasm_lr`, `merge_lr`,
`metaMDBG_lr`); o conjunto de regras do DAG não mudou em mais nada. Short reads
intactas, como exigido.

**Bug pré-existente encontrado e corrigido (fora de (c)):** a trilha **HiFi não
rodava desde antes deste roadmap**. Em `rules/qc.smk` o `filtlong_lr` recebia,
quando `lr_tech != "ont"`, o caminho
`{sample}/lr_filtered/{sample}_placeholder.fastq.gz` — **arquivo que nenhuma
regra produz** → `MissingInputException` antes mesmo da montagem. Como o
`porechop_lr` já trata HiFi internamente (faz `cat` das reads cruas), o
`filtlong_lr` passou a consumir sempre a saída dele: uma trilha só para as duas
tecnologias, sem condicional. O condicional antigo tinha ainda um segundo
problema — o `and LONG_READS` era morto, já que a regra só existe dentro do
`if LONG_READS:`.

Depois do conserto o HiFi monta o DAG: **159 jobs**, com
`porechop_lr` (pass-through) → `filtlong_lr` → `metaMDBG_lr`, sem Flye nem
Medaka. ONT e short reads não mudaram nada.

### (d) Remover `merge_contigs`, `merge_lr` e `mmseqs2` (dedup)

> **Pré-requisitos: (b), (c) e (h).** (b)/(c) porque o dedup só pode sair
> depois de sobrar um montador por trilha. (h) porque as regras que rodam
> no representante global precisam sair do caminho do hub `rep_seq.fasta`
> antes dele ser desmontado — ver "Ordem de execução".

**Esta é a mudança de maior risco.** O hub deixa de ser
`{sample}_rep_seq.fasta` e passa a ser o contig do MEGAHIT — 22 referências a
reapontar (`viral_detection`, `mapping`, `prok_binning`, `cobra`, `quast`,
`coassembly`, `Snakefile`).

Justificativa: o `easy-linclust` existia para colapsar redundância **entre
montadores**. Com um montador por trilha essa fonte deixa de existir.

**Precedente decisivo:** a trilha de co-assembly **já roda assim** — MEGAHIT →
`contigs.fa` → tudo, sem merge e sem dedup — e está em produção
(`results/coassembly/RNG/` tem corrida completa com bins, checkm2, gtdbtk e
viral). O padrão de destino já existe e está testado; é reapontar, não inventar.

Ganho colateral: sem dedup, o COBRA passa a receber a mesma montagem contra a
qual o BAM foi mapeado, que é o que ele exige.

### (e) Mover o limiar de comprimento em vez de apagá-lo

| onde | hoje | destino |
|---|---|---|
| MEGAHIT | `--min-contig-len 3000` | **1000** |
| `merge_contigs` | filtra 3000 de novo | removido em (d) |
| detecção viral | **sem piso** | **`viral_min_contig`, 3000** (parâmetro novo) |
| binners | clampam sozinhos | inalterado |
| COBRA | recebe filtrado | recebe as curtas |

Racional: o filtro do `merge_contigs` é **peso morto** — o MEGAHIT já filtra na
montagem (`assembly.smk:50,59`), então ele roda sobre contigs que já passaram
pelo mesmo corte. E os binners se protegem: MetaBAT2 tem piso rígido de 1500
(já clampado), vRhyme recusa abaixo de 2000 (já clampado), SemiBin2 tem o dele.

**Mas o limiar não pode ir a zero:** a detecção viral não tem piso nenhum, e
geNomad/VirSorter2 em contig curta é a maior fonte isolada de FP viral. Soltar
o limiar trabalharia contra o motivo de ter tirado o VIBRANT.

Custo a antecipar: baixar para 1000 bp multiplica a contagem de contigs — o
mapeamento fica bem mais pesado.

### (f) `bacphlip_votu` — repor lifestyle — FEITO

O BACPHLIP **já está instalado** (`envs/env_reads_classify.yaml`, via pip) e o
script `scripts/reads_classify/bacphlip_lifestyle.py` já existe — mas hoje recebe
genomas de referência do sylph, não contigs montados.

Regra nova: entra o FASTA viral filtrado por qualidade CheckV, sai lifestyle;
religa a coluna do `make_votu_table`.

**Parâmetro `bacphlip_min_quality`, default `HQ+Complete`.** Reusar o vocabulário
que já existe em `VIRAL_KEEP_TIERS` (`rules/votu_catalog.smk:202` — o tier `mq`
já está definido como "Complete/HQ/MQ or completeness >= 50%").

**Por que o default é conservador:** o BACPHLIP decide pela presença de domínios
associados à lisogenia (integrase, recombinase, repressor). Em MQ falta até
metade do genoma; se o pedaço ausente carregava a integrase, um fago lisogênico
é chamado de virulento. **O erro não é simétrico — a fragmentação enviesa na
direção de "virulento"** e infla artificialmente o percentual de líticos. Rodar
em MQ é permitido, mas é escolha consciente do usuário. Registrar em coluna qual
tier foi usado.

### (g) `eggnog_viral` — repor AMG como *putative* — FEITO

**Nada a instalar.** A vapor já tem as duas peças:
- `rule prodigal_viral` (`taxonomy.smk:76`) já produz os ORFs virais (`.faa`)
- eggNOG-mapper 2.1.15 já instalado e `eggnog_db` já baixado (`config.yaml:179`),
  com KEGG, CAZy e GO

AMG, na definição operacional do VIBRANT e do DRAM-v, é **gene metabólico KEGG
em genoma viral**. Rodar eggNOG sobre o `.faa` viral e filtrar KOs metabólicos
entrega o mesmo tipo de candidato com zero DB novo.

**Nomear como "putative AMGs" no relatório e nas colunas, nunca como "AMGs"
cravado** — decisão explícita do usuário. Justificativa: existe literatura de
2025 alertando que chamada de AMG sofre muita má-anotação e exige inspeção
manual do contexto genômico. Isso vale igual para DRAM-v e VIBRANT: não estamos
trocando método curado por bruto, estamos trocando um não-curado por outro mais
barato.

### (h) Princípio: computar no representante, herdar no membro

> **Executar ANTES de (d).** As duas mexem nos mesmos consumidores; fazer
> (h) primeiro reduz o conjunto que (d) precisa reapontar — ver
> "Ordem de execução".

**Esta é a mudança estruturalmente mais importante do roadmap.** Não é um
conserto pontual — é o modelo que o resto das regras deveria seguir.

#### O princípio

Um vOTU é um agrupamento a 95% ANI / 85% AF — nível de espécie pelo padrão ICTV.
Todo membro de um vOTU é, por construção, a mesma entidade biológica que o
representante. Logo:

> **Toda análise que depende só da sequência roda UMA vez, no representante
> global. O membro herda o resultado por join na provenance.**

Só permanece por amostra o que depende de **dados daquela amostra** — não da
sequência.

| natureza da análise | onde roda | exemplos |
|---|---|---|
| depende só da sequência | **uma vez, no representante global** | predição de ORF, anotação funcional, taxonomia, lifestyle, AMG, qualidade |
| depende de reads/bins da amostra | **por amostra** | recrutamento, cobertura, abundância, presença/ausência, predição de hospedeiro |

A tabela por amostra deixa de ser um recálculo e passa a ser uma **visão**: junta
a anotação global (via `provenance.tsv`) com a abundância local.

#### Por que isso é possível hoje

A parte difícil já está pronta. O `votu_catalog_pool` prefixa cada contig com sua
origem (*"Contig IDs are only unique within an assembly, so they are prefixed
with their source"*) e emite `provenance.tsv`. **Pool global com rastreio de
origem já existe** — é exatamente o que essa arquitetura precisa.

As regras `bacphlip_votu` e `eggnog_viral` (itens f e g) **já nasceram assim** e
servem de molde: global, no representante, com o membro herdando no
`make_votu_table.py`.

#### O ganho

Com 32 amostras, cada regra migrada deixa de rodar 32 vezes e passa a rodar 1.
Para `pharokka`/`phold` isso é a diferença entre horas e minutos. E elimina uma
classe inteira de inconsistência: hoje nada garante que duas amostras anotem o
mesmo representante do mesmo jeito.

#### Estado verificado (2026-08-18)

Sete regras por amostra consomem `votu_catalog_reps`. Verificado input a input.
**Cuidado com o método:** referências `rules.X.output` são por amostra sem conter
a string `{sample}` — uma checagem por texto literal dá falso positivo. Foi assim
que a primeira contagem desta análise saiu errada.

O `CLAUDE.md` registra que o catálogo global de vOTU *"Replaces the former
per-sample clustering"*. A parte difícil está pronta: o `votu_catalog_pool`
prefixa cada contig com sua origem (*"Contig IDs are only unique within an
assembly, so they are prefixed with their source"*) e emite `provenance.tsv` —
ou seja, **pool global com rastreio de origem já existe**.

Mas nem toda regra a jusante foi movida junto. Sete regras por amostra consomem
`votu_catalog_reps`. Verificado input a input (cuidado: referências
`rules.X.output` são por amostra sem conter a string `{sample}` — uma checagem
por texto literal dá falso positivo):

| regra | inputs | veredito |
|---|---|---|
| `prodigal_viral` | **só o FASTA global** | **N× idêntico — desperdício confirmado** |
| `pharokka` | global + CheckV **por amostra** | **BUG confirmado — não anota nada** |
| `viral_taxonomy` | + geNomad/MMseqs por amostra | **BUG confirmado — geNomad sai do merge** |
| `genome_map_virus` | + CheckV por amostra | **BUG confirmado — zero mapas** |
| `genome_map_phage` | + phold/pharokka por amostra | quebrado em cascata pelo `pharokka` |
| `phist` | + `bins_dir` por amostra | legítimo |
| `make_votu_table` | dados por amostra | **correto** — já namespaceia (ver abaixo) |

**`prodigal_viral`**: entrada exclusivamente global, saída por amostra. Com 32
amostras são 32 execuções produzindo bytes idênticos. Deve virar regra global.

#### O bug de namespace de ID (analisado em 2026-08-18)

O `pharokka` **não estava só suspeito — está quebrado**, e não sozinho.

**Causa raiz.** O commit `8ef8bb4` (2026-08-13, *"remove a cadeia per-amostra e
reponta os 9 consumidores no catálogo global"*) trocou o input dos consumidores
de `viral_votu_reps` (por amostra) para `votu_catalog_reps` (global). Os dois
produzem **espaços de nome de ID diferentes**:

| fonte | ID |
|---|---|
| `viral_votu_reps` (antigo, por amostra) | `MEGAHIT_k141_10006` |
| `votu_catalog_reps` (global) | `S1\|MEGAHIT_k141_10006` |

O `build_pool` prefixa por origem de propósito. Mas os *joins por ID* dos
consumidores não foram movidos junto: continuam cruzando o catálogo global com
tabelas por amostra que têm ID nu. Nenhum deles falha — todos devolvem conjunto
vazio, que passa por resultado biológico.

**Três regras afetadas, verificadas uma a uma:**

1. **`pharokka`** (`annotation.smk:82-107`) monta `hq_set` do CheckV da amostra
   (IDs nus) e filtra o FASTA global (namespaced). Zero match → `hq_fa` vazio →
   cai no `if not hq_set or getsize == 0` → grava *"No HQ phages found —
   skipping"* e toca outputs vazios. Sai como *"essa amostra não tem fago de alta
   completude"*. Em cascata: `phold` recebe GBK vazio, `genome_map_phage`
   também.
2. **`genome_map_virus`** (`genome_map_universal.py:714-733`), independente do
   pharokka: `comp_map` vem do CheckV (nu), `seqs` do FASTA global (namespaced),
   então `comp_map.get(gid, 0.0)` sempre devolve 0.0 e nenhum genoma passa do
   `min_comp`. Zero mapas.
3. **`viral_taxonomy`** (`taxonomy.smk:376-434`) itera os contigs do FASTA global
   (namespaced) e busca `genomad_tax.get(contig)`, cujas chaves são o `seq_name`
   do geNomad **daquela amostra** (nu). O geNomad roda sobre o `rep_seq.fasta` da
   amostra, então nunca casa — **o geNomad sai do merge de taxonomia**. O MMseqs
   escapa porque roda sobre as proteínas do `prodigal_viral`, que já é global.

**Reprodução** (cópia literal da lógica do `pharokka`):

```
hq_set (do CheckV da amostra) : ['MEGAHIT_k141_10006', 'MEGAHIT_k141_777']
headers do catálogo global    : ['S1|MEGAHIT_k141_10006', 'S2|...', 'S1|...']
sequências extraídas          : 0
```

**Por que ninguém viu.** O `results/` da corrida atual **não tem
`votu_catalog/`** — os outputs de pharokka em disco são da arquitetura anterior
(headers `k141_2550`, sem prefixo). O código atual nunca rodou até o fim; a
regressão está latente desde 13/08.

**A convenção era conhecida e só não foi propagada:** `votu_catalog.smk:79` monta
a chave como `f"{source_id}|{contig}"`, e `make_votu_table.py:309` compara com
`f"{sample}|{mem}" == rep` enquanto busca o CheckV com o ID nu. Esses dois estão
certos e servem de molde.

**O conserto não é adicionar o prefixo no join.** Mesmo com o namespace certo, a
semântica continua errada: cada amostra filtraria o catálogo global pelo CheckV
dela, e um representante originário da S2 seria aceito ou descartado pela tabela
da S1. O conserto é o próprio (h) — as regras viram globais e usam a completude
do próprio catálogo (`_load_catalog_completeness()`, que já namespaceia certo). O
namespace se resolve como consequência.

**Isso reclassifica o (h): não é só otimização de custo, é correção de
resultado.**

As duas regras adicionadas em (f) e (g) já nasceram globais e não aumentam esta
dívida.

---

## Achados que não devem ser redescobertos

### O "consenso" viral hoje é uma união

`viral_consensus_mode: "hybrid"` com `min_viral_tools: 2`, mas a lógica
(`viral_detection.smk:252-287`) faz `len(tools) >= 2 **OR** n in high_conf`, e:
- geNomad roda com `--min-score 0.7` contra um portão `score_genomad_min: 0.5`
  → **toda** chamada do geNomad passa
- VIBRANT entrava no `high_conf` se `type` fosse lytic/lysogenic → **toda**
  chamada passava

Ou seja o consenso é **≈ união** das ferramentas, e o `min_viral_tools` é
praticamente decorativo nesse modo. Isso explica o excesso de FP. Depois de (a),
revisar se `hybrid` ainda faz sentido ou se vira `score` puro.

### Descartado: strobealign como mapeador — por ora

Avaliado a fundo em 2026-08-18. Ver `docs/ANALISE_TOOLS_VOMIX_METAFUN.md` e as
referências abaixo.

- Nem o vOMIX nem o metaFun usam strobealign (o vOMIX usa **minimap2 + CoverM**).
- Existe paper de metagenômica — **AEMB** (Pan et al. 2025) — mas ele é sobre
  viabilizar binning **multi-against-multi** (N² alinhamentos), não sobre trocar
  o alinhador.
- A vantagem só aparece com **≥30 amostras**. Nos dados atuais: 32 amostras em 7
  grupos de co-assembly (~4-5 por grupo) — **nenhum dos dois modos da vapor é o
  modo em que o AEMB ganha**.
- No oceano (n=109) o AEMB deu **−11,3% bins HQ** vs Bowtie2(concat). Efeito
  dependente de ambiente.
- **Limitação estrutural:** Methods 4.2.5 do paper diz textualmente que o AEMB
  *"cannot give the per-base coverage"*. O MetaBAT2 da vapor depende do
  `jgi_summarize_bam_contig_depths`, que usa profundidade **e variância** por
  base. O AEMB não pode alimentar o MetaBAT2.
- O BAM continua necessário de todo jeito (vRhyme, CoverM, COBRA, MetaBAT2).

**Reabrir só se** a vapor for reestruturada para binning multi-against-multi.

### Descartado: `Binette(MetaBAT2 + SemiBin2)` está correto

O metaFun mediu: MetaBAT2 sozinho dá menos MAGs e **pior ARI**; o refinamento do
DAS Tool melhora acurácia de bp **mas reduz ARI** no marinho e na rizosfera (no
gut humano, o oposto). A lição **não** é "menos binners" — é que refinamento não
é de graça e o efeito depende do ambiente. A escolha atual da vapor está alinhada
com a evidência. **Não mexer.**

### Descartado: PPanGGOLiN

Exige vários MAGs da mesma espécie e responde pergunta de genômica comparativa.
O eixo do metaFun é invertido ao nosso (ele é comparativa com metagenômica de
apoio; a vapor é viroma com procarioto de apoio).

### Assimetria do vOMIX

Ele prega "uma ferramenta bem escolhida" no viral e usa **quatro binners**
(VAMB + MetaBAT2 + MaxBin2 + CONCOCT → DAS Tool) no procariótico, sem benchmark
que sustente. Não importar o lado procarioto deles por autoridade.

### As 6 ferramentas virais do vOMIX são só benchmark

Gating verificado no `Snakefile` deles: `viral-benchmark` não está no
`viral-end-to-end` (o módulo default). Padrão que vale copiar — mover VirSorter2
para um módulo opcional de benchmark em vez do caminho default.

---

## Pendências menores

- 28 `touch {output.done}` cosméticos restantes (regras que falham alto; o
  relatório vê "existe" em vez de "ok"). Baixa prioridade.
- Apagar as branches `backup/master-pre-rewrite` e
  `backup/unit-wildcard-pre-rewrite` — **só depois** de confirmar que o pull deu
  certo em todas as máquinas. São a rede de segurança da reescrita de autoria.
- Apagar a branch `refactor/unit-wildcard` (já mergeada).
- `config_amazon_12-08-26.yaml` (untracked) tem chaves órfãs `vcontact3_db` e
  `vcontact3_ver`; ganhará `vibrant_base` órfã depois de (a). Inofensivas — o
  Snakemake ignora chaves desconhecidas.
- Se algum dia entrar strobealign ou inStrain: a imagem BioContainers está bem
  atrás do bioconda (strobealign 0.9.0 vs 0.17.0), então pinar hoje faria o
  `scripts/check_env_container_sync.py` acusar divergência. Resolver com imagem
  custom no GHCR (a infra já existe: `ghcr.io/lymf/vapor-genome-map` com
  `custom: true`, mais `.github/workflows/build-containers.yml`).

---

## Ideias ainda não avaliadas

Nenhuma foi analisada a fundo — não tratar como aprovadas.

- **dbCAN** (CAZymes) — lacuna real, o metaFun tem e a vapor não. Disponível:
  bioconda 5.2.9, quay.io `5.2.9--pyhdfd78af_0`.
- **hostile** para remoção de hospedeiro (bioconda 2.0.2) — substituiria
  bwa-mem2 + filtro samtools.
- **PhaBOX2** (PhaTYP, CHERRY) — bioconda 2.1.13.
- **CheckV-PyHMMER** — só compensa acima de ~10⁵ contigs.
- **inStrain** — nível de cepa.

---

## Como verificar cada passo

Invariante do DAG: o conjunto ordenado de arquivos de output do dry-run.

```bash
conda activate snakemake
SP=/tmp/.../scratchpad   # helpers: snap.sh, ruleedit.py, config_dagtest.yaml
snakemake -n -p --configfile $SP/config_dagtest.yaml > $SP/dag_X.txt
bash $SP/snap.sh $SP/dag_X.txt $SP/outputs_X.txt
diff $SP/outputs_ANTERIOR.txt $SP/outputs_X.txt
```

Baseline em 2026-08-18 (pós status, pré-VIBRANT): **1616 jobs, 2625 outputs**.

**Regra de ouro:** mudança que só mexe em corpo de `shell:` deve dar diff
**vazio**. Mudança que remove ferramenta deve dar um delta **explicável item a
item** — nunca aceitar "mudou mais ou menos o esperado".

### Lição registrada sobre método

Já houve um erro caro neste projeto: uma comparação automatizada que colapsava
diffs de corpo executável **não viu** um `*.fa` vs `*.fna` no `gunc`, o que teria
feito a trilha de grupo achar zero bins **em silêncio**. Ver §0 de
`docs/AUDITORIA_COASSEMBLY_PARES.md`.

> Contagem de diff normalizado serve para **priorizar**, nunca para **decidir**.

Vale igual para trabalho delegado a subagente: a spec pode ser mecânica, a
revisão não.

---

## Referências

Pan, S., Tolstoganov, I., Sahlin, K., Martin, M., Zhao, X.-M., & Coelho, L. P. (2025). AEMB: efficient abundance estimation for metagenomic binning. *bioRxiv*. https://doi.org/10.1101/2025.07.30.667338

Sahlin, K. (2022). Strobealign: flexible seed size enables ultra-fast and accurate read alignment. *Genome Biology*, 23, 260. https://doi.org/10.1186/s13059-022-02831-7

Aroney, S. T. N., Newell, R. J. P., Nissen, J. N., Camargo, A. P., Tyson, G. W., & Woodcroft, B. J. (2025). CoverM: read alignment statistics for metagenomics. *Bioinformatics*, 41(4), btaf147. https://doi.org/10.1093/bioinformatics/btaf147

Pratama, A. A., Bolduc, B., Zayed, A. A., Zhong, Z.-P., Guo, J., Vik, D. R., Gazitúa, M. C., Wainaina, J. M., Roux, S., & Sullivan, M. B. (2021). Expanding standards in viromics: in silico evaluation of dsDNA viral genome identification, classification, and auxiliary metabolic gene curation. *PeerJ*, 9, e11447. https://doi.org/10.7717/peerj.11447

Shaffer, M., Borton, M. A., McGivern, B. B., Zayed, A. A., La Rosa, S. L., Solden, L. M., Liu, P., Narrowe, A. B., Rodríguez-Ramos, J., Bolduc, B., Gazitúa, M. C., Daly, R. A., Smith, G. J., Vik, D. R., Pope, P. B., Sullivan, M. B., Roux, S., & Wrighton, K. C. (2020). DRAM for distilling microbial metabolism to automate the curation of microbiome function. *Nucleic Acids Research*, 48(16), 8883–8900. https://doi.org/10.1093/nar/gkaa621

Hockenberry, A. J., & Wilke, C. O. (2021). BACPHLIP: predicting bacteriophage lifestyle from conserved protein domains. *PeerJ*, 9, e11396. https://doi.org/10.7717/peerj.11396

### Documentos irmãos neste repo

- `docs/ANALISE_TOOLS_VOMIX_METAFUN.md` — análise ferramenta a ferramenta das duas pipelines de referência
- `docs/BENCHMARK_VOMIX_METAFUN.md` — arquitetura e performance
- `docs/VAPOR_TOOLS_MAP.md` — inventário por estágio das ferramentas da vapor
- `docs/AUDITORIA_COASSEMBLY_PARES.md` — laudo par a par da unificação amostra/grupo
