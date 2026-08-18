# Laudo — auditoria par a par: regras per-sample × gêmeas `coassembly_*`

Levantado em 2026-08-17 na branch `refactor/unit-wildcard`. Método: extração de cada par, remoção de docstrings e comentários, normalização de caminhos e nomes de wildcard, diff do **código restante**. Script em `scratchpad/audit.py`.

Objetivo: separar duplicação acidental (unificável) de divergência intencional (tem que ficar) e de *drift* (bug).

---

## Resumo

| Categoria | Pares | Ação |
|---|---|---|
| **A** — zero/quase-zero diferença de código | 9 | Unificar por herança. Risco nulo. |
| **B** — divergência legítima de input/params | 6 | Unificar com override no `use rule`. |
| **C** — **drift real (regressão)** | 2 | **Não unificar antes de corrigir.** Ver §3. |
| **D** — divergência intencional de algoritmo | 3 | Manter separado, documentar o porquê. |

Ou seja: dos 20 pares, **15 são unificáveis** — bem mais do que eu estimei quando vi só o diff bruto (o docstring inflava a contagem). E 2 escondem um bug ativo.

---

## 1. Categoria A — só docstring (unificar já)

Diferença de código medida = 0, exceto onde anotado.

| Par | Difs de código | Observação |
|---|---|---|
| `abricate` | 0 | |
| `mmseqs_taxonomy_viral` | 0 | |
| `defensefinder_viral` | 0 | |
| `defensefinder` | 0 | |
| `argnorm_normalize` | 0 | |
| `phold` | 0 | |
| `dbapis_viral` | 1 | só o guard `if` do bloco |
| `extract_kegg_kos` | 1 | só o guard `if` do bloco |
| `pharokka` | 2 | só a origem do input (`votu_catalog_reps` vs `viral_votu_reps`) |

São ~700 linhas de código idêntico mantidas em duplicata.

## 2. Categoria B — divergência legítima (unificar com override)

Todas têm a mesma causa: **a trilha de grupo usa VAMB, a per-sample usa Binette**, então mudam os diretórios de bin e a regra de origem. É uma diferença de *entrada*, não de lógica — exatamente o que o `use rule ... with:` sobrescreve.

| Par | Difs | Natureza |
|---|---|---|
| `gunc` | 4 | `bins/binette/final_bins` → `vamb/run/bins`; `rules.binette` → `rules.vamb_cobinning` |
| `amr_consensus` | 5 | idem |
| `bakta` | 6 | idem |
| `galah_derep` | 6 | idem + `rules.checkm2` → `rules.checkm2_group` |
| `deeparg` | 7 | idem |
| `filter_viral_for_prok` | 14 | idem + nome do output (`_rep_seq_nonviral` vs `_contigs_nonviral`) |

## 3. Categoria C — DRIFT REAL: `amrfinderplus` e `rgi_card`

**Este é o achado principal do laudo.**

As duas versões per-sample foram corrigidas para registrar o desfecho real da execução; as gêmeas de coassembly **ficaram para trás**:

```python
# rules/defense_amr.smk (per-sample) — CORRIGIDO
def write_empty(msg, status):
    write_status(str(output.done), status)      # 'ok' | 'skipped: ...' | 'failed: ...'
...
try:
    ...">> {log} 2>&1"
except sp.CalledProcessError as exc:
    amr_err = exc.returncode
    write_status(str(output.done), f"failed: amrfinder exit {amr_err}")

# rules/coassembly.smk (gêmea) — AINDA COM O PADRÃO ANTIGO
def write_empty(msg):
    Path(str(output.done)).touch()              # done.txt VAZIO
...
"...>> {log} 2>&1 || echo '[amrfinderplus] WARNING: amrfinder failed' >> {log}"
```

A gêmea engole o código de saída com `|| echo` e escreve um `done.txt` vazio. É **exatamente o bug que o próprio `write_status()` documenta** no `Snakefile`:

> *"An empty done.txt makes 'the tool crashed' and 'the tool found nothing' indistinguishable downstream, which is how a disk-full AMRFinderPlus run was read as a biological zero."*

O padrão é sistemático, não pontual:

| Arquivo | `write_status` | `touch` vazio |
|---|---|---|
| `rules/coassembly.smk` | **0** | **46** |
| `rules/defense_amr.smk` | 6 | 14 |
| `rules/annotation.smk` | 0 | 13 |
| `rules/prok_binning.smk` | 0 | 7 |
| `rules/taxonomy.smk` | 0 | 7 |

**Impacto real, delimitado:** `load_tool_status()` (`scripts/report/data_loaders.py:61`) itera sobre `samples` e sobre as regras globais. Grupos de co-assembly **não estão em `samples`**, então esses `done.txt` nunca chegam a ser lidos pelo relatório. Consequência: uma ferramenta que falha na trilha de co-assembly não aparece como lacuna nem como erro — **ela simplesmente não aparece**. O relatório está correto no que lê; o problema é que não lê nada de coassembly.

Vale notar que o `_read_status_file()` já trata `done.txt` vazio como `unknown` (não como sucesso), e `tool_failed()` já converte `unknown` em lacuna. A infraestrutura da correção existe — só não foi ligada na trilha de grupo.

**Recomendação:** propagar `write_status` para as 46 chamadas de `coassembly.smk` **antes** de unificar esses dois pares, e estender `STATUS_TRACKED_TOOLS` para cobrir grupos. É correção de bug, não refatoração — e muda o que o relatório mostra.

## 4. Categoria D — divergência intencional (manter separado)

### `vrhyme`
```python
rule coassembly_vrhyme:
    input:
        bams = lambda wc: expand(..., sample=GROUPS[wc.group])   # TODOS os BAMs do grupo
```
O per-sample usa **um** BAM; o de grupo usa **todos**, para cobertura diferencial — análogo ao co-binning do VAMB no lado procariótico. São algoritmos diferentes. Unificar seria bug.

### `phist`
Mesma causa da categoria B (`bins_dir` VAMB vs Binette), mas o bloco vizinho `_coas_final_inputs` é lógica exclusiva de grupo. Unificável em parte; o entorno fica.

### `viral_taxonomy` — o mais grave depois da categoria C

**Idêntico:** o schema de saída (mesmas 15 colunas, mesma ordem), o dicionário `_PRIORITY`, o algoritmo de merge (rank mais profundo vence, empate por confiança), o fallback para geNomad e o tratamento de `unclassified`.

**Diferente:** a gêmea zera duas das quatro fontes, em código:

```python
# rules/taxonomy.smk (per-sample) — 4 fontes
vc3_tax    = {...}                                     # vConTACT3, populado
custom_tax = _mmseqs_lca_rollup(input.custom_hits, _RANKS)   # MMseqs2/custom, populado

# rules/coassembly.smk (gêmea) — 2 fontes
vc3_tax    = {}      # hardcoded vazio
custom_tax = {}      # hardcoded vazio
```

Consequências:

1. **A trilha de grupo classifica com metade da evidência.** O mesmo contig, analisado per-sample ou via co-assembly, pode receber táxons diferentes — e o de grupo é estritamente mais pobre.
2. **As colunas `vc3_status`, `vc3_novel_anchor`, `custom_rank`, `custom_lineage`, `custom_n_proteins` saem sempre vazias** na saída de grupo, mas continuam no cabeçalho. Quem lê o TSV não distingue "não rodou" de "não teve hit" — a **mesma classe** de ambiguidade da categoria C, agora no nível de coluna.
3. O `_PRIORITY` da gêmea ainda lista `vcontact3: 0` e `mmseqs_custom: 2`, fontes que ela nunca produz. Código morto que sugere um comportamento que não existe.

**Nota de oportunidade:** com o vConTACT3 saindo da pipeline (decisão de 2026-08-17), a distância entre as duas versões cai de 2 fontes para 1 (só o MMseqs2/custom). Faz sentido **remover o vConTACT3 primeiro e reavaliar este par depois** — a unificação fica bem mais barata, e talvez trivial.

---

## 5. Ordem sugerida

1. **Categoria A** (9 pares) — herança direta, invariante do DAG preservado por construção. Ganho ~700 linhas.
2. **Categoria B** (6 pares) — herança com override de input/params. Invariante preservado.
3. **Remover vConTACT3** — decisão já tomada; encolhe o problema do `viral_taxonomy`.
4. **Corrigir o `write_status` na trilha de coassembly** (46 pontos) + estender `STATUS_TRACKED_TOOLS` para grupos. **Muda o relatório** — é correção de bug, precisa ser um commit próprio e anunciado.
5. **Categoria C** (2 pares) — unificar só depois de (4).
6. **`viral_taxonomy`** — reavaliar depois de (3).
7. **`vrhyme`** — não unificar. Adicionar um comentário explicando que a duplicação é deliberada, para ninguém "consertar" isso depois.

## Verificação usada

Baseline com 32 amostras e 7 grupos (`coassembly.viral` e `coassembly.binning` ligados): **1641 jobs, 2683 arquivos de output**. Toda unificação feita até agora manteve o diff do invariante vazio. Os passos (3) e (4) são os primeiros que **vão** mudar comportamento — por isso ficam separados e sinalizados.
