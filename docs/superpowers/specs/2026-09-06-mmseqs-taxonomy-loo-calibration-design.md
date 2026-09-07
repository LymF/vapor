# Calibração leave-one-out de cutoffs por rank/taxon para os rollups de taxonomia viral/prok (mmseqs2)

Status: **roadmap — design aprovado, implementação não iniciada.**

## Motivação

Investigando o VITAP (https://github.com/DrKaiyangZheng/VITAP, Zheng, K. et
al. Nat Commun 16, 2226, 2025, https://doi.org/10.1038/s41467-025-57500-7)
como candidato a ferramenta de classificação taxonômica viral, rodamos o
pipeline completo (ICTV VMR + fallback UniRef90) contra o catálogo real de
vOTUs do run "Amazon" (`results-18-08-26_bin_coassembly`, 1.852 vOTUs) e
comparamos com a abordagem já usada pelo vapor (`votu_mmseqs_taxonomy` via
INPHARED + `votu_mmseqs_taxonomy_custom` + geNomad, mergeadas em
`votu_taxonomy`).

Resultado: a abordagem atual do vapor resolveu ~17x mais vOTUs até Family e
~26x mais até Genus do que o VITAP nesse catálogo (794 vs 46 em Family, 514
vs 20 em Genus) — principalmente porque o INPHARED é uma base muito mais
densa para fagos ambientais do que o ICTV VMR (curado, só espécies
oficialmente ratificadas). O VITAP também expôs um bug de dado conhecido: a
família `Rhodogtaviriformidae` no VMR MSL41 tem "genomas" de referência que
na verdade são cromossomos bacterianos quase inteiros (600kb–2.5Mb, quando um
GTA real tem ~14kb) — confirmado como corrupção de range por Excel já citada
no próprio README do VITAP — e isso gerou 6 falsos positivos de família no
nosso catálogo de teste (verificados: a abordagem atual do vapor classifica
essas 6 sequências como `Unclassified`, corroborando que o call do VITAP era
espúrio).

**Conclusão da investigação**: não adotar o VITAP como está — ele perde para
o que já existe e carrega um bug de referência. Mas o *algoritmo* dele tem
uma ideia genuinamente melhor que o que o vapor faz hoje na etapa de
consenso: um **cutoff de score calibrado empiricamente por
leave-one-out**, em vez de um limiar fixo ou de exigir unanimidade.

### O problema concreto no vapor hoje

Três pontos do código fazem hoje a mesma coisa — agregam múltiplas
chamadas de taxonomia por-PROTEÍNA (já uma LCA conservadora do próprio
mmseqs2) em uma única chamada por-GENOMA/CONTIG, exigindo que **100% dos
proteins concordem** em cada rank, senão trunca para o rank anterior:

- `_mmseqs_lca_rollup()` em `rules/votu_catalog.smk:224` — compartilhada
  por `votu_mmseqs_taxonomy` (INPHARED) e `votu_mmseqs_taxonomy_custom`
  (ex: IMG/VR).
- `load_mmseqs_taxonomy_prok()` em `scripts/report/data_loaders.py:708` —
  usada pelo `custom_prok_mmseqs_db` (ex: IMG_NR).

Essa unanimidade é uma decisão deliberada e documentada — ambos os locais
citam von Meijenfeldt et al. 2019 (CAT/BAT) e o problema de "spurious
specificity" de abordagens best-hit — e o loader do prok diz explicitamente
que "a vote would just be a softer version of the same 'spurious
specificity' this path exists to avoid". **Essa objeção é válida contra um
limiar arbitrário tipo "maioria simples"**, mas não se aplica a um cutoff
calibrado empiricamente contra a própria base de referência: a diferença é
que o segundo é validado com dado real (quanto de concordância entre
proteínas é preciso, POR RANK e POR TAXON, para que a chamada continue
batendo com a verdade conhecida), não um número escolhido a dedo.

O efeito prático da unanimidade estrita é sub-classificação (falso
negativo), não falso positivo: um contig com 50 proteínas em que 49
concordam em "Family X" e 1 diverge (HGT, prófago, erro de anotação) perde
o rank inteiro e vira "só Class" — descartando sinal real. Foi exatamente
esse tipo de perda, olhando o VITAP, que motivou essa investigação (embora
lá o mecanismo fosse outro: cutoffs por taxon já calibrados, mas contra uma
referência pequena e com bug).

## Ideia calibrada (o que reaproveitar do VITAP)

VITAP calibra, por taxon e por rank, um cutoff de **score** via leave-one-
out: `same_taxon_score.min()` vs `cross_taxon_score.max()`, cutoff = ponto
médio entre os dois (ou `same.min() * 0.75` sem dado cross-taxon
suficiente). Ele faz isso de forma barata: **um único self-alignment
all-vs-all** da base de referência contra ela mesma (visto rodando no nosso
teste real: "[INFO] Self-aligning of ICTV reference proteins" seguido de
"[INFO] Calculating best-fit taxonomic threshold for Species/Genus/.../
Realm") — o "leave-one-out" é só filtrar, da matriz de scores pareados já
computada, os hits de uma sequência contra si mesma.

**Correção feita durante o `writing-plans` (2026-09-06)**: a versão inicial
deste spec propunha calibrar um cutoff de **bitscore**, igual ao VITAP.
Verificado contra um output real do `votu_mmseqs_taxonomy` do run Amazon:
o `mmseqs createtsv` que os 3 rollups consomem hoje só tem as colunas
`qseqid, taxid, rank, name, lineage` — **sem nenhum score por proteína**.
Isso porque `mmseqs taxonomy` já resolve cada proteína internamente via
voto majoritário ponderado por score (`--vote-mode 2`, `--majority`) e só
expõe o resultado final da LCA, não os scores brutos por hit. Um cutoff de
bitscore não tem onde encaixar nos rollups sem mudar as 3 rules de
taxonomia para também emitirem score por proteína (via
`--tax-output-mode`) — investigação não feita, comportamento dessa flag
não verificado contra uma DB real.

**Decisão**: calibrar uma **fração de concordância mínima** por rank/taxon,
não um score — plugável direto nos dados que os rollups já produzem hoje,
sem mudar nenhuma das 3 rules de taxonomia existentes. O all-vs-all
all-vs-all continua sendo a fonte de dado (nenhuma mudança de escopo/custo
computacional), só a métrica final calibrada muda de "bitscore" para
"fração de proteínas concordando".

### Mecânica exata

Reaproveita dois artefatos que `prepare_mmseqs_taxdb.py` já produz para
qualquer formato (`img`/`ncbi`/`inphared`/`imgvr`), sem parsing
formato-específico:
- `mapping.tsv` (`{protein_id}\t{taxid}`) — já escrito por `write_mapping()`.
- `taxdump/nodes.dmp` + `taxdump/names.dmp` — já escritos por
  `build_taxdump()`; dão `parent_of[taxid]`, `rank_of[taxid]`,
  `name_of[taxid]`, de onde sai `lineage_by_rank(taxid) -> {rank: nome}`
  andando até a raiz.

Como o taxid já é o nó folha do lineage (`build_taxdump` reusa o mesmo
taxid pra duas entradas com lineage idêntico), agrupar proteínas por taxid
agrupa por espécie/unidade taxonômica — não precisa reconstruir "o mesmo
genoma" a partir do protein_id (que exigiria parsing específico por
formato, incluindo o caso `img`, cujo genoma só é recuperável do HEADER
completo, não do protein_id isolado).

Passo a passo, para cada rank R (ex: Family):
1. Para cada proteína Q, achar seu melhor hit (`max(score)`) entre as
   linhas do self-align **cuja unidade-alvo (taxid do target) é diferente
   da unidade de Q** — isso já filtra o self-hit trivial (`query==target`)
   e qualquer hit contra outra proteína da mesma espécie, sem precisar de
   um segundo arquivo de mapeamento proteína→genoma.
2. `agree(Q, R)` = `lineage(unidade de Q)[R] == lineage(unidade do melhor
   hit)[R]`.
3. Por unidade U com taxon verdadeiro `T = lineage(U)[R]`:
   `f_same(U, R) = proporção das proteínas de U com agree(Q,R) == True`.
4. Por unidade U' com taxon verdadeiro `T' != T`: quando `agree == False`
   e o melhor hit apontou (erroneamente) para `T`, registra a fração de
   proteínas de U' "atribuídas" a T por engano — `f_cross(U', T, R)`.
5. `cutoff[R][T] = (min(f_same(U,R) para toda unidade U de taxon T) +
   max(f_cross(U',T,R) para toda unidade U' de outro taxon)) / 2`, ou
   `min(f_same) * 0.75` quando não há dado cross-taxon.

Isso é a mesma fórmula de ponto-médio do VITAP, só trocando "score de
alinhamento" por "fração de proteínas concordando" como unidade de medida
— e o resultado (`cutoff[R][T]`, um número entre 0 e 1) é diretamente
comparável à fração de proteínas de um contig/genoma novo que concordam
num taxon T naquele rank, que é exatamente o dado que `_mmseqs_lca_rollup`
e `load_mmseqs_taxonomy_prok` já calculam ao decidir se mantêm ou truncam
um rank.

## Escopo desta calibração

Cobre as **3 DBs de taxonomia viral/prok baseadas em mmseqs2** que o vapor
já mantém:

1. INPHARED (`votu_mmseqs_taxonomy`)
2. Custom viral, ex. IMG/VR (`votu_mmseqs_taxonomy_custom`)
3. Custom prok, ex. IMG_NR (`mag_mmseqs_taxonomy_prok`, consumido via
   `load_mmseqs_taxonomy_prok`)

As três já passam por `scripts/prepare_mmseqs_taxdb.py` na etapa de setup
(INSTALL.md) — calibrar as três ali, com o mesmo módulo, tem custo
incremental baixo.

**Fora de escopo** (não avaliado nesta investigação): GTDB-Tk,
`mag_mmseqs_taxonomy_prok`'s outros consumidores fora do rollup de
report, e qualquer coisa fora dos 3 rollups acima.

## Arquitetura

### Componente novo: `scripts/mmseqs_taxonomy_calibration.py`

Funções puras, sem I/O de rede nem chamada a `mmseqs` diretamente (isso
fica no script que orquestra):

- `load_taxdump(taxdump_dir) -> (parent_of, rank_of, name_of)` — parseia
  `nodes.dmp`/`names.dmp` (formato já escrito por `build_taxdump()`).
- `lineage_by_rank(taxid, parent_of, rank_of, name_of) -> {rank: nome}` —
  anda de `taxid` até a raiz (id `1`) coletando `{rank_of[t]: name_of[t]}`
  por ancestral.
- `load_protein_taxid_map(mapping_path) -> {protein_id: taxid}` — lê o
  `mapping.tsv` que `write_mapping()` já escreve.
- `iter_self_align_rows(self_align_path) -> iterator[(query, target,
  score)]` — lê o TSV sem header produzido por `mmseqs convertalis
  --format-output query,target,bits`.
- `compute_rank_cutoffs(self_align_rows, protein_taxid, taxid_lineage,
  ranks) -> dict` — implementa a mecânica descrita acima (melhor hit
  cross-unidade por proteína, `f_same`/`f_cross` por rank/taxon, cutoff =
  ponto médio ou `min(f_same) * 0.75`). Retorna `{rank: {taxon: cutoff}}`,
  só com ranks/taxa que tiveram dado suficiente (>= 1 unidade em
  `f_same`) — ranks/taxa ausentes do dict sinalizam "sem calibração para
  isso", tratado pelo consumo como não-calibrado.
- `write_rank_cutoffs(cutoffs, out_path)` / `load_rank_cutoffs(out_path) ->
  dict | None` — serialização JSON; `load_rank_cutoffs` retorna `None` se
  o arquivo não existe.

### Geração: novo modo em `scripts/prepare_mmseqs_taxdb.py`

Flag nova (`--calibrate`), rodada manualmente pelo usuário depois de montar
a seqTaxDB — reaproveita o padrão idempotente do resto do arquivo: rodar o
comando de novo com `--calibrate` sobre uma DB já construída pula o
`createdb`/`createtaxdb` (já existem) e só executa a calibração. Não
autoexecuta dentro de uma rule Snakemake, mesma razão já documentada para o
build da própria seqTaxDB (evitar corrida entre samples sob `--cores >1`).

Fluxo:
```
seqTaxDB (já existe)
   │  mmseqs search seqTaxDB seqTaxDB result tmp     (all-vs-all, 1x manual)
   │  mmseqs convertalis seqTaxDB seqTaxDB result
   │      self_align.tsv --format-output query,target,bits
   ▼
self_align.tsv  (query, target, score — sem header)
   │  compute_rank_cutoffs(), usando mapping.tsv + taxdump/ já existentes
   ▼
rank_cutoffs.json   {realm: {taxon: cutoff}, ..., genus: {taxon: cutoff}}
```

`rank_cutoffs.json` é escrito ao lado da própria seqTaxDB (mesmo
diretório), não em `CATALOG_DIR` — é uma propriedade da base de
referência, não do catálogo de um run específico, e várias
rodadas/catálogos podem reusar a mesma DB calibrada uma vez.

### Consumo: mudança pontual nos 2 pontos de rollup existentes

Em `_mmseqs_lca_rollup()` (`rules/votu_catalog.smk`) e
`load_mmseqs_taxonomy_prok()` (`scripts/report/data_loaders.py`): em vez de
`if len(set(level)) == 1` (unanimidade), calcular a fração do taxon
majoritário entre as proteínas que alcançam aquele rank e comparar contra
`rank_cutoffs[rank][taxon]`. Essa fração substitui `n_proteins` como sinal
de confiança reportado (hoje um contador bruto, não normalizado pelo total
de proteínas do contig/genoma).

Ambos os pontos continuam com assinatura de entrada/saída compatível com
quem os consome hoje (`votu_taxonomy`, o report) — a mudança é só a regra
de decisão interna de quando truncar.

**Fallback sem calibração (decisão fechada, escopo desta fase)**: os dois
call sites recebem os cutoffs via um parâmetro (`cutoffs: dict | None`).
Quando `cutoffs` é `None` (arquivo `rank_cutoffs.json` ausente) **ou**
quando `rank`/`taxon` não aparece no dict calibrado (dado insuficiente
naquele ponto específico), a decisão cai no comportamento atual de
unanimidade (`fração == 1.0` é o único caso aceito) — zero regressão para
quem não rodou `--calibrate` ainda, e degradação graciosa por
rank/taxon dentro de uma DB parcialmente calibrada. Em
`rules/votu_catalog.smk`, `cutoffs` é carregado uma vez via
`load_rank_cutoffs(seqtaxdb_path + ".rank_cutoffs.json")` dentro da regra
`votu_taxonomy` (que já lê os dois `hits` de INPHARED/custom) e passado
para `_mmseqs_lca_rollup`; em `data_loaders.py`, `load_rank_cutoffs` é
chamado com o path derivado de `CUSTOM_PROK_MMSEQS_DB` (já disponível em
`paths_d`/config no momento em que o report é montado).

## Testes

`compute_rank_cutoffs`, `lineage_by_rank`, `load_protein_taxid_map` e
`iter_self_align_rows` são puras (recebem/retornam dicts, listas ou
iteráveis, sem tocar disco além do parsing de um path recebido) —
testáveis com um `nodes.dmp`/`names.dmp`/`mapping.tsv` sintéticos pequenos
em `tmp_path`, seguindo o padrão já estabelecido em `scripts/mag_catalog.py`
(21 testes unitários puros, citado no CLAUDE.md do projeto).
`load_rank_cutoffs` testa o caso de arquivo ausente e o de JSON válido. A
mudança em `_mmseqs_lca_rollup`/`load_mmseqs_taxonomy_prok` testa os 3
casos: unanimidade (sem cutoffs), fração acima do cutoff calibrado (mantém
o rank), fração abaixo (trunca).

## Referências

Zheng, K., Zhang, W., Su, F. et al. VITAP: a high precision tool for DNA
and RNA viral classification based on meta-omic data. *Nature
Communications*, 16, 2226 (2025).
https://doi.org/10.1038/s41467-025-57500-7

von Meijenfeldt, F. A. B., Arkhipova, K., Cambuy, D. D., Coutinho, F. H., &
Dutilh, B. E. Robust taxonomic classification of uncharted microbial
sequences and bins with CAT and BAT. *Genome Biology*, 20, 217 (2019).
https://doi.org/10.1186/s13059-019-1817-x
