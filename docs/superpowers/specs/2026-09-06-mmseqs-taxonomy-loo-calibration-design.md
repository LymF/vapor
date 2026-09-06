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

VITAP calibra, por taxon e por rank, um cutoff de score via leave-one-out:
`same_taxon_score.min()` vs `cross_taxon_score.max()`, cutoff = ponto médio
entre os dois (ou `same.min() * 0.75` quando não há dado cross-taxon
suficiente). Ele faz isso de forma barata: **um único self-alignment
all-vs-all** da base de referência contra ela mesma (visto rodando no nosso
teste real: "[INFO] Self-aligning of ICTV reference proteins" seguido de
"[INFO] Calculating best-fit taxonomic threshold for Species/Genus/.../
Realm"), não um re-alinhamento por genoma excluído — o "leave-one-out" é só
filtrar, da matriz de scores pareados já computada, os hits de uma sequência
contra si mesma.

A mesma técnica se aplica ao mmseqs2: `mmseqs search seqTaxDB seqTaxDB
result tmp` roda a DB contra ela mesma uma única vez; dali sai uma tabela de
scores pareados com taxonomia de cada lado, e o cutoff por rank/taxon sai de
um groupby simples em cima dessa tabela.

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

Duas funções puras, sem I/O de rede nem chamada a `mmseqs` diretamente
(isso fica no script que orquestra):

- `compute_rank_cutoffs(self_align_rows, taxonomy_map) -> dict` — recebe a
  tabela de scores pareados (query, target, score, taxid_query,
  taxid_target) já carregada e o mapa taxid→lineage; filtra hits de uma
  sequência contra si mesma (mesmo ID de origem, não apenas mesmo taxid —
  duas cópias do mesmo genoma não devem se validar mutuamente); agrupa por
  rank (realm..genus, mesmo esquema ICTV das outras rules) e por taxon;
  calcula `same_taxon_score.min()` vs `cross_taxon_score.max()`; cutoff =
  ponto médio, ou `same.min() * 0.75` sem dado cross-taxon. Retorna
  `{rank: {taxon_name: cutoff}}`.
- `load_rank_cutoffs(cutoffs_path) -> dict | None` — lê
  `rank_cutoffs.json`; retorna `None` se o arquivo não existe (sinaliza "DB
  não calibrada ainda").

### Geração: novo modo em `scripts/prepare_mmseqs_taxdb.py`

Flag nova (ex: `--calibrate`), rodada manualmente pelo usuário depois de
montar a seqTaxDB (mesmo padrão manual de setup do resto do arquivo — não
autoexecuta dentro de uma rule Snakemake, mesma razão já documentada para o
build da própria seqTaxDB: evitar corrida entre samples sob `--cores >1`).

Fluxo:
```
seqTaxDB (já existe)
   │  mmseqs search seqTaxDB seqTaxDB result tmp   (all-vs-all, 1x manual)
   ▼
self_align.tsv  (query, target, score, taxid_query, taxid_target)
   │  compute_rank_cutoffs()
   ▼
rank_cutoffs.json   {realm: {taxon: cutoff}, ..., genus: {taxon: cutoff}}
```

`rank_cutoffs.json` é escrito ao lado da própria seqTaxDB (mesmo diretório),
não em `CATALOG_DIR` — é uma propriedade da base de referência, não do
catálogo de um run específico, e várias rodadas/catálogos podem reusar a
mesma DB calibrada uma vez.

### Consumo: mudança pontual nos 2 pontos de rollup existentes

Em `_mmseqs_lca_rollup()` (`rules/votu_catalog.smk`) e
`load_mmseqs_taxonomy_prok()` (`scripts/report/data_loaders.py`): em vez de
`if len(set(level)) == 1` (unanimidade), calcular a fração do taxon
majoritário entre as proteínas que alcançam aquele rank e comparar contra
`rank_cutoffs[rank][taxon]` (quando existir) antes de decidir manter ou
truncar o rank. Essa fração substitui `n_proteins` como sinal de confiança
reportado (hoje um contador bruto, não normalizado pelo total de proteínas
do contig/genoma).

Ambos os pontos continuam com assinatura de entrada/saída compatível com
quem os consome hoje (`votu_taxonomy`, o report) — a mudança é só a regra de
decisão interna de quando truncar.

## Gap conhecido / decisão adiada

**O que fazer quando `rank_cutoffs.json` não existe** (DB instalada antes
desta feature, ou usuário que nunca rodou `--calibrate`) ainda não foi
implementado — fica como TODO explícito para a fase de implementação. A
decisão já tomada é: cair no comportamento atual de unanimidade (zero
regressão para quem não recalibrar), mas o encaixe fino disso nos 3 call
sites (onde checar, como logar que está em modo não-calibrado) precisa ser
resolvido durante o `writing-plans`, não nesta fase de design.

## Testes (quando for implementado)

`compute_rank_cutoffs` é pura (recebe/retorna dicts, sem tocar disco nem
subprocess) — testável com uma tabela sintética pequena de scores/taxa,
seguindo o padrão já estabelecido em `scripts/mag_catalog.py` (21 testes
unitários puros, citado no CLAUDE.md do projeto). `load_rank_cutoffs` testa
o caso de arquivo ausente e o de JSON válido.

## Referências

Zheng, K., Zhang, W., Su, F. et al. VITAP: a high precision tool for DNA
and RNA viral classification based on meta-omic data. *Nature
Communications*, 16, 2226 (2025).
https://doi.org/10.1038/s41467-025-57500-7

von Meijenfeldt, F. A. B., Arkhipova, K., Cambuy, D. D., Coutinho, F. H., &
Dutilh, B. E. Robust taxonomic classification of uncharted microbial
sequences and bins with CAT and BAT. *Genome Biology*, 20, 217 (2019).
https://doi.org/10.1186/s13059-019-1817-x
