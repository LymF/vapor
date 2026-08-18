# Laudo — auditoria par a par: regras per-sample × gêmeas `coassembly_*`

Levantado em 2026-08-17 na branch `refactor/unit-wildcard`. **Revisado no mesmo dia** após a primeira versão se mostrar imprecisa — ver §0.

---

## 0. Correção da primeira versão deste laudo

A v1 classificou os pares com um comparador que normalizava caminhos e nomes de wildcard e diffava o resto. **Esse método subestimou as diferenças** e por pouco não causou um erro sério.

Duas falhas do comparador:

1. **Colapsava trechos do corpo executável.** Reportou `gunc` com "4 difs (input + bins_dir)", quando na verdade o `shell:` também diferia: `find -name "*.fa"` e `--file_suffix .fa` contra `*.fna`/`.fna`. Bins do Binette são `.fa`, os do VAMB são `.fna`. Converter no automático confiando nessa contagem teria feito a trilha de grupo **não encontrar bin nenhum** — e o resultado seria "0 bins", não um erro. Falha silenciosa, a pior categoria.
2. **Atribuía helpers de módulo à regra anterior.** Reportou que `deeparg` perdia `_has_data_rows()`. Falso: essa função é definida em nível de módulo em `defense_amr.smk`, logo **depois** da regra. O corpo do `deeparg` era idêntico. O mesmo artefato inflou `phist`, `bakta`, `viral_taxonomy` e outros.

**Método da v2, usado daqui em diante:** extrair o bloco `shell:`/`run:` de cada regra delimitando pela indentação, normalizar só os rótulos de log (`[nome]`), e diffar. Quando um par vai ser convertido, o diff é lido por inteiro antes — não por contagem.

**Lição operacional:** contagem de diff normalizado serve para *priorizar*, nunca para *decidir*. A decisão exige ler os dois corpos.

---

## 1. Resultado final

Dos 20 pares auditados, **19 foram unificados** por herança (`use rule ... as ... with:`). Só o `vrhyme` continua separado, de propósito.

`rules/coassembly.smk`: **3605 → 1948 linhas (−46%)**. 28 regras herdadas, 12 gêmeas próprias restantes — todas específicas de grupo, sem equivalente per-sample.

### Unificados — diferença só de docstring/rótulo (9)
`abricate`, `mmseqs_taxonomy_viral`, `defensefinder_viral`, `defensefinder`, `argnorm_normalize`, `phold`, `dbapis_viral`, `extract_kegg_kos`, `pharokka`.

Mais, nas levas seguintes: `amr_consensus`, `deeparg`, `prodigal_viral`, `checkv_vrhyme`, `phist` e `viral_taxonomy` — a v1 os listou como divergentes, mas os corpos eram idênticos byte a byte fora o rótulo do log.

### Unificados — divergência legítima de entrada, resolvida por override (4)
| Par | O que difere | Como foi resolvido |
|---|---|---|
| `gunc` | bins `.fa` vs `.fna` | `params.bin_ext` na regra base, sobrescrito na herdada |
| `bakta` | idem (`BIN_FA`) | idem |
| `galah_derep` | idem, em 5 pontos | idem — **mais o §2** |
| `filter_viral_for_prok` | entrada (rep_seq vs contigs), nome do output, `genomad_dir` | override de input/output; `genomad_dir` passou a derivar de `input.genomad` |

### Já unificados antes desta auditoria (8)
`virsorter2`, `genomad`, `vibrant`, `viral_consensus`, `checkv`, `eggnog_prok`, `amrfinderplus`, `rgi_card`.

---

## 2. Drift real encontrado — correções que existiam só do lado per-sample

Três casos onde a cópia de grupo estava **atrás** da per-sample. Unificar não foi limpeza, foi correção de bug, e cada um saiu em commit próprio marcado como tal.

### `amrfinderplus` e `rgi_card`
As versões per-sample ganharam `write_status()` (`ok` / `skipped: <motivo>` / `failed: <motivo>`); as gêmeas continuavam engolindo o exit code com `|| echo WARNING` e escrevendo `done.txt` vazio. É o caso que o próprio docstring de `write_status()` descreve:

> *"An empty done.txt makes 'the tool crashed' and 'the tool found nothing' indistinguishable downstream, which is how a disk-full AMRFinderPlus run was read as a biological zero."*

### `galah_derep` — três correções ausentes
1. `rm -rf {params.repdir}` antes de rodar. Sem isso o galah aborta quando o diretório de representantes existe e não está vazio, ou seja **um `snakemake` retomado sobre execução anterior falhava sempre**.
2. Captura do exit code (`GALAH_EXIT`) em vez de `|| echo`.
3. `done.txt` com status real em vez de `touch` vazio. Um fallback para os bins originais significa que a dereplicação **não aconteceu** — agora fica registrado em vez de passar por sucesso.

### O agravante comum: o relatório não lia nada disso
`load_tool_status()` iterava só sobre `samples`. Grupos vivem em `coassembly/<grupo>/` e **nunca eram lidos**, então uma ferramenta que falhasse na trilha de grupo não aparecia nem como lacuna nem como erro. Corrigido: a função aceita `groups` e indexa por `GROUP_STATUS_PREFIX`, com `ValueError` em caso de colisão com nome de amostra. Dois testes cobrem isso.

Escopo do que **não** foi feito: `coassembly.smk` tinha 46 `touch` vazios. Só 2 eram regressão (as gêmeas de rules já corrigidas per-sample). Os outros espelham regras per-sample que também usam `touch` simples — `annotation.smk` (13), `prok_binning.smk` (7), `taxonomy.smk` (7). Estender `write_status` ao pipeline inteiro é uma mudança de convenção separada, e agora mais barata: a leitura por grupo já existe.

---

## 3. Mantido separado de propósito

```python
rule coassembly_vrhyme:
    input:
        bams = lambda wc: expand(..., sample=GROUPS[wc.group])   # TODOS os BAMs do grupo
```
O per-sample usa **um** BAM; o de grupo usa **todos**, para cobertura diferencial — análogo ao co-binning do VAMB. São algoritmos diferentes. **Unificar seria bug.**

---

## 4. `viral_taxonomy` — resolvido ligando a fonte que faltava

O algoritmo sempre foi idêntico: mesmo schema, mesmo `_PRIORITY`, mesmo merge por rank mais profundo, mesmo fallback. O que diferia era a **base de evidência**: a gêmea zerava `custom_tax = {}` em código, então a trilha de grupo classificava com 2 fontes contra 3 da per-sample. O mesmo contig recebia taxonomia diferente conforme a trilha, e as colunas `custom_*` saíam sempre vazias mantendo-se no cabeçalho.

Resolvido (2026-08-17) ligando `coassembly_mmseqs_taxonomy_custom_viral` — herdado, com auto-skip se `custom_viral_mmseqs_db` não estiver configurado — e convertendo o merge para herança. As duas trilhas passam a usar as mesmas três fontes.

**Custo:** +1 execução de MMseqs2 por grupo. **Mudança de resultado:** sim, a taxonomia do co-assembly muda (para melhor, e agora comparável com a per-sample).

---

## 5. Gêmeas próprias restantes (12)

Nenhuma tem equivalente per-sample a herdar:

| Regras | Por quê |
|---|---|
| `viral_trimmed`, `skani_votu`, `skani_cluster`, `viral_votu_reps` | o per-sample usa o **catálogo global** de vOTU; o grupo faz clustering próprio |
| `index`, `map`, `sort`, `mapback`, `abundance` | mapeamento das amostras do grupo contra o co-assembly |
| `prok_bin_proteins`, `organize_outputs` | específicas de grupo |
| `vrhyme` | divergência intencional — ver §3 |

---

## 6. Verificação

Baseline: 32 amostras, 7 grupos, `coassembly.viral` e `coassembly.binning` ligados.

O invariante é o **conjunto de arquivos de output do DAG**, extraído do dry-run e diffado a cada leva. Nomes de regra mudam na refatoração; os outputs não podem.

| Etapa | Jobs | Outputs |
|---|---|---|
| Baseline | 1641 | 2683 |
| Após remoção do vConTACT3 | 1609 | 2619 |
| Após remoção do tier `hq_10kb` | 1609 | 2611 |
| Após ligar MMseqs2/custom no grupo | 1616 | 2625 |

Cada delta foi conferido item a item e continha **apenas** os caminhos esperados: −32 jobs/−64 outputs só de `.../vcontact3/`; −8 outputs só de `*hq_10kb*`; +7 jobs/+14 outputs só de `mmseqs_vs_custom.tsv` e `mmseqs_custom_viral_done.txt`, sem nada removido.

**Nenhuma das 19 unificações por herança alterou o invariante.** As três mudanças acima são de escopo, decididas explicitamente, não efeitos colaterais de refatoração.

As mudanças do §2 são de **conteúdo** do `done.txt` e de comportamento em resume — não alteram o conjunto de arquivos, por isso o invariante permanece válido como rede de segurança mas não as detecta. Elas são cobertas pelos testes (67, era 65).
