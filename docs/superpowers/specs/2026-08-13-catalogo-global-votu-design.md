# Design: Catálogo global de vOTU + abundância por recrutamento

**Data:** 2026-08-13
**Status:** Aprovado (estrutura) — pronto para plano de implementação
**Escopo:** Substituir o clustering de vOTU por amostra por um catálogo global único,
construído em um passo de skani sobre o pool completo, com abundância por amostra
derivada de recrutamento competitivo de leituras, e migração da anotação viral
(taxonomia, host prediction, anotação funcional) para os representantes do catálogo.

**Antecedente:** Bloco B do inventário em
`~/analyses/amazon/PIPELINE_erros_e_correcoes.md` §10. O Bloco A (falhas silenciosas,
clamp do vRhyme, `*` do Prodigal no RGI, idempotência do galah, remoção do COMEBin)
foi entregue em `6f4aa0e`.

---

## 1. Motivação

### 1.1 O defeito que reordena o problema

O documento de auditoria descreve a §10 como "o clustering de vOTU é por amostra".
A investigação mostrou algo pior: **o clustering nunca rodou, em nenhum nível.**

São três defeitos empilhados em `rules/viral_binning.smk:248-357` e no espelho
`rules/coassembly.smk:1707+`:

1. **Modo de saída errado.** `skani triangle` é chamado sem `--sparse`, e portanto
   escreve uma matriz triangular estilo PHYLIP:

   ```
   1871
   MEGAHIT_k141_27250
   MEGAHIT_k141_3145	0.00
   MEGAHIT_k141_68803	0.00	0.00
   ```

   O parser em `skani_cluster` espera lista de arestas (`qname rname ani af_q af_r`).
   Lendo a matriz densa, ele atribui `q` = nome do genoma e `r` = um valor numérico
   como `"0.00"`; a guarda `if q in neigh and r in neigh` falha em toda linha e
   **nenhuma aresta é criada**.

2. **Aligned fraction ausente.** O modo denso emite apenas ANI. O critério
   ICTV/Roux 2019 exige ANI ≥ 95 % **e** AF ≥ 85 %. Mesmo com o parser corrigido,
   a coluna de AF não existiria. `--sparse` é requisito, não otimização.

3. **Índices de coluna errados no modo esparso.** Verificado empiricamente com
   skani 0.3.2, a saída de `--sparse` é:

   | col | nome |
   |---|---|
   | 0 | `Ref_file` |
   | 1 | `Query_file` |
   | 2 | `ANI` |
   | 3 | `Align_fraction_ref` |
   | 4 | `Align_fraction_query` |
   | 5 | `Ref_name` |
   | 6 | `Query_name` |

   Os **nomes dos genomas estão nas colunas 5 e 6**, não 0 e 1 — as duas primeiras
   são caminhos de arquivo. Acrescentar `--sparse` sem corrigir os índices trocaria
   um bug silencioso por outro.

**Evidência nos dados da corrida de junho/2026:** `genomes=1871 clusters=1871` em
toda amostra, `genomes=857 clusters=857` em todo grupo de co-montagem, e cada linha
de `vOTU_clusters.tsv` é `X → X`. Os 44.123 "vOTUs" são os contigs não-redundantes
sem nenhum agrupamento.

### 1.2 O que muda com o clustering correto

Medido rodando skani 0.3.2 sobre os dados reais de junho:

| escopo | contigs | vOTUs reais | redução |
|---|---|---|---|
| 1 amostra (P06_TAP_3_958) | 772 | 674 | 12,7 % |
| pool de 6 amostras (4×RNG + 2×TAP) | 9.653 | 5.524 | 42,8 % |

No pool de 6 amostras a soma ingênua infla a riqueza em **1,75×**, e 10,4 % dos
vOTUs aparecem em mais de uma amostra. Com 32 amostras o fator é maior.

**Custo:** 17 s para 9.653 genomas (8 threads, `--slow --sparse`). Extrapolando por
n², ~6 min para as 44.123 do pool completo. O skani global não é o gargalo.

### 1.3 Por que um passo único, e não dois estágios

Clustering em dois estágios (local, depois global sobre os representantes) é prática
comum na literatura por tratabilidade, mas **não reproduz a partição de um passo
único**: a etapa global só compara representantes, então se o rep da amostra A tem
< 95 % ANI com o rep da amostra B mas ≥ 95 % com um *membro* do cluster de B, o elo
se perde.

Como o pool completo (44.123 sequências) roda em ~6 min, a aproximação não se
justifica aqui. **Decisão: passo único sobre o pool completo.**

Consequência: a cadeia por amostra (`skani_votu` → `skani_cluster` →
`viral_votu_reps`) deixa de alimentar qualquer coisa e é **removida**. Mantê-la
consertada recriaria os dois conceitos concorrentes de vOTU que este trabalho existe
para eliminar.

---

## 2. Arquitetura

Novo estágio global em `{OUTDIR}/votu_catalog/`, sem wildcard de amostra:

```
{sample}/viral/consensus/*_viral_nonredundant.fasta        (32 amostras)
coassembly/{group}/viral/consensus/*_viral_nonredundant.fasta  (7 grupos)
        │
        ├─► votu_catalog_pool     → pool.fasta
        │                           provenance.tsv
        ├─► votu_catalog_skani    → skani_ani.tsv        (--sparse --slow)
        ├─► votu_catalog_cluster  → vOTU_clusters.tsv    (IDs vOTU_00001…)
        ├─► votu_catalog_reps     → catalog_all_reps.fasta
        │                           catalog_mq_reps.fasta
        │                           catalog_hq_10kb_reps.fasta
        ├─► votu_catalog_map      → {sample}.catalog.sorted.bam   (por amostra)
        └─► votu_catalog_abundance → presence_matrix.tsv
                                     votu_abundance_matrix.tsv
```

### 2.1 `votu_catalog_pool`

Concatena os conjuntos não-redundantes de todas as amostras e grupos, **prefixando
cada ID** com a origem: `{origem}|{contig_id}`.

O prefixo é obrigatório, não estético: **3.668 IDs colidem entre amostras** nos dados
de junho (`MEGAHIT_k141_10006` existe em várias). Sem prefixo, o pool funde contigs
não relacionados e o catálogo fica corrompido de forma silenciosa.

Saídas:

- `pool.fasta` — todas as sequências, IDs prefixados.
- `provenance.tsv` — colunas `member_id`, `source_type` (`sample` \| `group`),
  `source_id`, `original_contig_id`.

Entrada de grupos condicionada a `coassembly.enabled` e `coassembly.viral`; quando
desligadas, o pool contém apenas amostras.

### 2.2 `votu_catalog_skani`

```
skani triangle -i pool.fasta -o skani_ani.tsv -t {threads} --slow --sparse
```

`--sparse` emite lista de arestas com aligned fraction, e apenas pares acima de
`--min-af`, o que mantém o arquivo tratável em 44 k genomas.

### 2.3 `votu_catalog_cluster`

Single-linkage guloso sobre as arestas, com o critério ICTV/Roux 2019:
`ANI ≥ votu_ani` **e** `max(af_q, af_r) ≥ votu_af`.

Parser com os índices corretos da §1.1: `ani = parts[2]`, `af_ref = parts[3]`,
`af_query = parts[4]`, `ref_name = parts[5]`, `query_name = parts[6]`. Auto-comparações
(`ref_name == query_name`) são descartadas.

Representante do cluster: maior completeness do CheckV entre os membros, empate
resolvido por ordem no pool — mesma regra do `skani_cluster` atual, agora aplicada
sobre o conjunto global. Requer os `checkv/quality_summary.tsv` de todas as amostras
e grupos como entrada, com as chaves de contig re-prefixadas para casar com o pool.

Saída `vOTU_clusters.tsv`: colunas `votu_id`, `representative`, `member`. Os IDs são
estáveis e sequenciais (`vOTU_00001`…), ordenados por tamanho de cluster decrescente
e, em empate, pelo ID do representante — para que uma re-execução sobre o mesmo pool
produza os mesmos rótulos.

**Validação obrigatória:** a regra falha alto se o número de clusters for igual ao
número de sequências de entrada, porque isso é exatamente a assinatura do defeito da
§1.1. Um catálogo sem nenhum agrupamento em 44 k genomas virais não é um resultado
biológico plausível.

### 2.4 `votu_catalog_reps`

Extrai as três camadas de representantes que a cadeia por amostra produzia, agora
uma vez só:

- `catalog_all_reps.fasta` — um por vOTU; base do recrutamento e do relatório.
- `catalog_mq_reps.fasta` — MQ+ (Complete/HQ/MQ ou completeness ≥ 50 %); alimenta
  taxonomia, PHIST e anotação.
- `catalog_hq_10kb_reps.fasta` — HQ+/Complete e ≥ 10 kb; alimenta o vConTACT3.

Mesmos limiares do `viral_votu_reps` atual (`VIRAL_MIN_QUALITY`), sem mudança de
critério — apenas de escopo.

---

## 3. Presença e abundância por amostra

Requisito explícito: o relatório precisa mostrar **em qual amostra está cada vírus**,
e a diferença entre a contagem por amostra e o total.

Dois sinais independentes, reportados lado a lado e **nunca fundidos**:

### 3.1 Presença por montagem

Deriva de `provenance.tsv` × `vOTU_clusters.tsv`, sem custo computacional: um vOTU
está "montado" na amostra X se algum de seus membros veio de X.

Responde "este vírus montou nesta amostra". É enviesado — um vírus presente mas de
cobertura baixa demais para montar fica invisível.

### 3.2 Presença por recrutamento

`votu_catalog_map` mapeia as leituras limpas de cada amostra contra
`catalog_all_reps.fasta` (bwa-mem2 para leituras curtas, minimap2 para longas —
mesma escolha que `mapping.smk` já faz), e `votu_catalog_abundance` roda CoverM
sobre os BAMs.

Um vOTU conta como presente na amostra quando **≥ 75 % do comprimento do
representante está coberto** (Roux et al., 2017), com os mesmos filtros de
identidade e comprimento alinhado já usados em `coverm_viral`
(`--min-read-percent-identity 95`, `--min-read-aligned-length 45`,
`--contig-end-exclusion 75`).

Responde "este vírus está nesta amostra", tenha montado ou não. É o sinal
comparável com a literatura.

**Custo:** `bwa_mem` leva no máximo 37 s por amostra na corrida de junho, então a
passada extra sai por ~20 min nas 32. O mapeamento não é o gargalo.

### 3.3 Saídas

- `presence_matrix.tsv` — vOTU × amostra, valores `assembled` / `recruited` /
  `both` / `absent`.
- `votu_abundance_matrix.tsv` — vOTU × amostra, métrica configurada em
  `COVERM_METHOD`, restrita às células que passam o corte de cobertura.

---

## 4. Migração da anotação para o catálogo

As oito dependências mapeadas do `viral_votu_reps` perdem o wildcard `{sample}` e
passam a rodar uma vez sobre os representantes do catálogo:

| regra | arquivo | entrada antes | entrada depois |
|---|---|---|---|
| `prodigal_viral` | `taxonomy.smk:80` | `mq_fasta` da amostra | `catalog_mq_reps.fasta` |
| `viral_taxonomy` | `taxonomy.smk:438` | `mq_fasta` da amostra | `catalog_mq_reps.fasta` |
| `vcontact3` | `taxonomy.smk:346` | `hq_10kb_fasta` da amostra | `catalog_hq_10kb_reps.fasta` |
| `phist` | `host_prediction.smk:23` | `mq_fasta` da amostra | `catalog_mq_reps.fasta` |
| anotação (×3) | `annotation.smk:35,461,507` | `mq_fasta` da amostra | `catalog_mq_reps.fasta` |
| `make_votu_table` | `viral_binning.smk:522` | `all_fasta` da amostra | `catalog_all_reps.fasta` |

Ganho colateral: cada população viral é anotada uma vez, não uma vez por amostra em
que aparece. Nos dados de junho isso elimina a anotação redundante dos 10 % de vOTUs
compartilhados (fração que cresce com o número de amostras).

**Consequência para o relatório:** deixa de existir "taxonomia da amostra X" como
dado direto. A visão por amostra é reconstruída por junção de *anotação global ×
matriz de presença*. É o modelo correto e é o que torna os números por amostra
comparáveis entre si, mas é a maior parte do esforço de implementação — maior que
o skani global.

### 4.1 Regras removidas

- `skani_votu` (`viral_binning.smk:248`)
- `skani_cluster` (`viral_binning.smk:291`)
- `viral_votu_reps` (`viral_binning.smk:393`)

Consumidores repontados para o catálogo: `make_votu_table`, `votu_abundance`
(`abundance.smk:63`) e `coverm_viral` (`abundance.smk:11`).

### 4.2 Co-montagem preservada

A aba de co-montagem continua exibindo vOTUs **por grupo**, como decidido. As regras
`coassembly_skani_votu` e `coassembly_skani_cluster` permanecem, mas **recebem a
mesma correção da §1.1** (`--sparse` + índices certos + validação anti-singleton).
Sem isso a aba continua reportando `N contigs = N vOTUs`.

Os contigs dos grupos entram também no pool global (§2.1), de modo que um vírus
recuperado só na co-montagem existe no catálogo — sem que isso mexa na contagem por
grupo mostrada na aba.

---

## 5. Relatório

| item | hoje | depois |
|---|---|---|
| KPI de riqueza | soma de `viral_consensus` por amostra (`overview.js:46`) — infla | contagem de vOTUs do catálogo |
| vOTUs por amostra | contagem local, não comparável entre amostras | nº de vOTUs do catálogo presentes na amostra |
| em qual amostra está cada vírus | não existe | matriz de presença navegável, com os dois sinais da §3 |
| curva de acumulação | só na track de co-montagem (`data_loaders.py:1681`) | vale para o conjunto todo |
| aba de co-montagem | vOTUs por grupo | inalterada, com o parser corrigido |

Novos carregadores em `scripts/report/data_loaders.py`: `load_votu_catalog()`,
`load_votu_presence()`, `load_votu_abundance_matrix()`. O comentário em
`data_loaders.py:1684-1688`, que hoje explica por que a curva de acumulação só existe
na co-montagem, deixa de valer e deve ser reescrito.

---

## 6. Configuração

Chaves novas em `config.yaml`:

| chave | padrão | efeito |
|---|---|---|
| `votu_catalog_enabled` | `true` | liga o estágio global |
| `votu_presence_min_coverage` | `75.0` | % do representante coberto para contar presença (Roux 2017) |

Chaves reaproveitadas sem mudança de significado: `votu_ani` (95,0), `votu_af` (85,0),
`coverm_method`, `viral_min_quality`.

`votu_clustering_enabled` perde o sentido no escopo por amostra e passa a governar o
catálogo, ou é substituída por `votu_catalog_enabled` — decisão de implementação.

---

## 7. Verificação

1. **Anti-regressão do defeito original:** sobre `pool.fasta`, o número de clusters
   tem de ser estritamente menor que o de sequências. A regra falha alto caso
   contrário (§2.3).
2. **Prefixação:** nenhum ID duplicado em `pool.fasta`; contagem de sequências do
   pool igual à soma das entradas.
3. **Reprodutibilidade dos IDs:** duas execuções sobre o mesmo pool produzem os
   mesmos rótulos `vOTU_*`.
4. **Coerência das matrizes:** todo vOTU com presença por montagem em X tem membro
   com `source_id = X` no `provenance.tsv`.
5. **Baseline empírico:** o pool de 6 amostras já medido (9.653 → 5.524 vOTUs,
   42,8 %) serve de teste de regressão para o clustering.
6. **DAG:** `snakemake -n` com as tracks viral e prok ligadas, antes e depois, sem
   erro e sem referência às regras removidas.

---

## 8. Fora de escopo

- §4, §5 e §8 do inventário (congelar config, comparar config efetivo, vConTACT3).
- O molde de falha silenciosa nos ~115 pontos restantes; o Bloco A cobriu cinco
  ferramentas. As regras novas deste design gravam status real desde o início.
- Migração da anotação **procariótica** para qualquer catálogo — só o lado viral.

---

## Referências

Roux, S., Emerson, J. B., Eloe-Fadrosh, E. A., & Sullivan, M. B. (2017). Benchmarking
viromics: an in silico evaluation of metagenome-enabled estimates of viral community
composition and diversity. *PeerJ*, 5, e3817. https://doi.org/10.7717/peerj.3817

Roux, S., Adriaenssens, E. M., Dutilh, B. E., Koonin, E. V., Kropinski, A. M.,
Krupovic, M., … Eloe-Fadrosh, E. A. (2019). Minimum Information about an Uncultivated
Virus Genome (MIUViG): a community consensus on standards and best practices for
describing genome sequences from uncultivated viruses. *Nature Biotechnology*, 37(1),
29–37. https://doi.org/10.1038/s41587-018-0100-8

Nayfach, S., Páez-Espino, D., Call, L., Low, S. J., Sberro, H., Ivanova, N. N.,
Proal, A. D., Fischbach, M. A., Bhatt, A. S., Hugenholtz, P., & Kyrpides, N. C.
(2021). Metagenomic compendium of 189,680 DNA viruses from the human gut microbiome.
*Nature Microbiology*, 6(7), 960–970. https://doi.org/10.1038/s41564-021-00928-6

Gregory, A. C., Zablocki, O., Zayed, A. A., Howell, A., Bolduc, B., & Sullivan, M. B.
(2020). The Gut Virome Database Reveals Age-Dependent Patterns of Virome Diversity in
the Human Gut. *Cell Host & Microbe*, 28(5), 724–740.e8.
https://doi.org/10.1016/j.chom.2020.08.003

Shaw, J., & Yu, Y. W. (2023). Fast and robust metagenomic sequence comparison through
sparse chaining with skani. *Nature Methods*, 20(11), 1661–1665.
https://doi.org/10.1038/s41592-023-02018-3
