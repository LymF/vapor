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

Dos 20 pares auditados, **15 foram unificados** por herança (`use rule ... as ... with:`), 2 continuam separados de propósito e 3 seguem pendentes.

`rules/coassembly.smk`: **3605 → 2185 linhas (−39%)**. 23 regras herdadas, 16 gêmeas próprias restantes.

### Unificados — diferença só de docstring/rótulo (9)
`abricate`, `mmseqs_taxonomy_viral`, `defensefinder_viral`, `defensefinder`, `argnorm_normalize`, `phold`, `dbapis_viral`, `extract_kegg_kos`, `pharokka`.

Mais, na segunda leva: `amr_consensus` e `deeparg` — a v1 os listou como divergentes, mas os corpos eram idênticos byte a byte fora o rótulo do log.

### Unificados — divergência legítima de entrada, resolvida por override (4)
| Par | O que difere | Como foi resolvido |
|---|---|---|
| `gunc` | bins `.fa` vs `.fna` | `params.bin_ext` na regra base, sobrescrito na herdada |
| `bakta` | idem (`BIN_FA`) | idem |
| `galah_derep` | idem, em 5 pontos | idem — **mais o §2** |
| `filter_viral_for_prok` | entrada (rep_seq vs contigs), nome do output, `genomad_dir` | override de input/output; `genomad_dir` passou a derivar de `input.genomad` |

### Já unificados antes desta auditoria (6)
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

## 3. Mantidos separados de propósito

### `vrhyme`
```python
rule coassembly_vrhyme:
    input:
        bams = lambda wc: expand(..., sample=GROUPS[wc.group])   # TODOS os BAMs do grupo
```
O per-sample usa **um** BAM; o de grupo usa **todos**, para cobertura diferencial — análogo ao co-binning do VAMB. São algoritmos diferentes. **Unificar seria bug.**

### `viral_taxonomy`
O algoritmo é idêntico: mesmo schema de saída, mesmo `_PRIORITY`, mesmo merge por rank mais profundo, mesmo fallback. O que difere é a **base de evidência**: a gêmea zera `custom_tax = {}`, então a trilha de grupo classifica com 2 fontes (MMseqs2/INPHARED + geNomad) contra 3 da per-sample.

Com a saída do vConTACT3 (2026-08-17) a distância caiu de 2 fontes para 1. Unificar exigiria ligar o `mmseqs_taxonomy_custom_viral` na trilha de grupo — isso **acrescenta** uma execução de ferramenta, não é refatoração. Decisão pendente.

---

## 4. Pendentes

| Par | Situação |
|---|---|
| `prodigal_viral` | corpo aparenta ser idêntico; falta confirmar e converter |
| `checkv_vrhyme` | idem |
| `phist` | `bins_dir` VAMB vs Binette (como categoria B), mas o bloco vizinho `_coas_final_inputs` é lógica exclusiva de grupo |
| `viral_trimmed`, `skani_votu`, `skani_cluster`, `viral_votu_reps` | sem regra base per-sample equivalente (o per-sample usa o catálogo global de vOTU) |
| `index`, `map`, `sort`, `mapback`, `abundance`, `prok_bin_proteins`, `organize_outputs` | específicos de grupo |

---

## 5. Verificação

Baseline: 32 amostras, 7 grupos, `coassembly.viral` e `coassembly.binning` ligados.

O invariante é o **conjunto de arquivos de output do DAG**, extraído do dry-run e diffado a cada leva. Nomes de regra mudam na refatoração; os outputs não podem.

| Etapa | Jobs | Outputs |
|---|---|---|
| Baseline | 1641 | 2683 |
| Após todas as heranças | 1609 | 2611 |

As duas reduções são explicadas e foram conferidas item a item: −32 jobs e −64 outputs da remoção do vConTACT3 (o diff continha **só** caminhos `.../vcontact3/`), −8 outputs da remoção do tier `hq_10kb` (só caminhos `*hq_10kb*`). **Nenhuma unificação por herança alterou o invariante.**

As mudanças do §2 são de **conteúdo** do `done.txt` e de comportamento em resume — não alteram o conjunto de arquivos, por isso o invariante permanece válido como rede de segurança mas não as detecta. Elas são cobertas pelos testes (67, era 65).
