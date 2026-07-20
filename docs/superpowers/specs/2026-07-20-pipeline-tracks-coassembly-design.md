# Design: Reestruturação em tracks + co-assembly/co-binning

**Data:** 2026-07-20
**Status:** Aprovado (estrutura) — pronto para plano de implementação
**Escopo:** Restruturação do VAPOR em tracks selecionáveis por configuração, com um
novo track opcional de co-assembly/co-binning, e rearranjo dos binners procarióticos.

---

## 1. Motivação

Hoje o VAPOR roda como um bloco monolítico: o `rule all` (~100 linhas) sempre puxa o
pipeline inteiro (assembly → viral + prok → integração → report). Não há como rodar
"só reads" ou "só a track viral" sem gastar tempo/RAM com o resto.

Objetivos desta reestruturação (confirmados com o usuário):

1. **Tracks independentes** — poder rodar reads-based, viral ou prokariótico
   isoladamente, por parâmetro.
2. **Entry/exit points** — poder parar em estágios (`--until`).
3. **Co-assembly + co-binning** — nova análise paralela opcional (grupo de amostras).
4. **Rearranjo de binners** — VAMB sai dos binners per-sample e passa a ser usado
   exclusivamente no co-binning (joga a favor da natureza multi-amostra do VAMB).

**Decisão explícita:** NÃO migrar para o sistema de `module` do Snakemake 8. As tracks
são fortemente acopladas pelo hub central (`rep_seq.fasta`, mapping/depth, sample
discovery, report) e a integração cruza as fronteiras dos módulos — o custo de reescrita
supera o ganho. ATLAS e aviary (os dois Snakemake metagenômicos mais citados com essa
estrutura) usam `include:` + composição de alvos + CLI, não módulos. Mantemos o mesmo
padrão. Migração para módulos fica como decisão futura, caso o objetivo passe a ser
publicar tracks como ferramentas independentes reutilizáveis.

---

## 2. Diferencial preservado: cruzamento vírus↔bactéria

O valor central do VAPOR é o cruzamento vírus↔bactéria (host prediction, defesa↔anti-defesa,
report integrado). A separação em tracks **não** perde isso — a integração vira uma
**camada condicional** que acende quando as tracks de que depende estão presentes.

- **Host prediction (PHIST)** — liga fagos aos MAGs procarióticos recuperados localmente.
  Requer track viral **e** prok. (iPHoP foi removido do pipeline; PHIST é a única via de
  host prediction, e é intrinsecamente cross-domain.)
- **Arms-race defesa↔anti-defesa** — sistemas de defesa nas bactérias vs anti-defesa nos
  fagos. Requer as duas tracks.
- **Report adaptativo** — já renderiza abas conforme a presença de dados; a reestruturação
  apenas formaliza isso.

| Track(s) rodada(s)   | Abas do report                                                        |
|----------------------|-----------------------------------------------------------------------|
| só reads             | Overview, Reads Survey, Diversity(reads), About                       |
| só viral             | + Viral, Annotation(viral)                                            |
| só prok              | + Prokaryotic, Host&Defense(só defesa bacteriana + AMR)              |
| viral + prok         | + cruzamento: PHIST, arms-race, abundância/diversidade cruzada        |
| + co-assembly        | + seção "Co-assembly" (vOTUs de grupo e/ou MAGs de grupo)             |
| full                 | tudo + Reads Survey                                                   |

---

## 3. Superfície de configuração (`config.yaml`)

```yaml
# ── Seleção de tracks ──────────────────────────
tracks:
  reads:  false     # classificação reads-based (sylph) — independe de assembly
  viral:  true      # track viral (assembly-based, per-sample)
  prok:   true      # track procariótica (assembly-based, per-sample)

# ── Integração cruzada (acende com viral+prok; flag pra desligar) ──
use_host_defense: true    # PHIST + arms-race defesa↔anti-defesa

# ── Co-assembly (track opcional, nível de grupo) — OFF por padrão ──
# A co-assembly produz contigs de grupo que alimentam DOIS consumidores
# (simétrico ao per-sample): detecção viral e/ou co-binning.
coassembly:
  enabled:  false         # OFF por padrão — feature opt-in
  grouping: metadata      # metadata (coluna 'group') | all | none
  viral:    true          # (quando enabled) detecção viral na co-assembly → vOTUs de grupo
  binning:  true          # (quando enabled) VAMB co-binning na co-assembly → MAGs de grupo

# VAMB multi-split (independe de co-assembly): binning multi-amostra sobre
# assemblies individuais concatenadas. Ligar junto com coassembly.binning = comparação.
cobinning_multisplit: false   # OFF por padrão
```

**Padrão de execução:** com `coassembly.enabled: false` e `cobinning_multisplit: false`
(ambos default), o pipeline roda exatamente como hoje (per-sample puro). As análises de
grupo são inteiramente opt-in. Os sub-toggles `viral`/`binning` só têm efeito quando
`coassembly.enabled: true`.

### Variáveis derivadas no Snakefile

```python
TRACK_READS   = config.get("tracks", {}).get("reads", False)
TRACK_VIRAL   = config.get("tracks", {}).get("viral", True)
TRACK_PROK    = config.get("tracks", {}).get("prok",  True)
USE_HOST_DEFENSE = config.get("use_host_defense", True)

COASSEMBLY_ENABLED  = config.get("coassembly", {}).get("enabled", False)
COASSEMBLY_GROUPING = config.get("coassembly", {}).get("grouping", "metadata")
COASSEMBLY_VIRAL    = config.get("coassembly", {}).get("viral", True)
COASSEMBLY_BINNING  = config.get("coassembly", {}).get("binning", True)
COBINNING_MULTISPLIT = config.get("cobinning_multisplit", False)

# integração só existe quando ambas as tracks assembly-based rodam
INTEGRATION_ENABLED = TRACK_VIRAL and TRACK_PROK and USE_HOST_DEFENSE
```

### Validações no startup (falha cedo, com mensagem clara)

- `coassembly.enabled and not (coassembly.viral or coassembly.binning)` ⇒ aviso
  (co-assembly ligada mas sem consumidor: só produz contigs, útil só p/ QUAST comparativo).
- `coassembly.grouping == "metadata"` exige o TSV de metadata com coluna `group`.
- `coassembly.grouping == "none"` ⇒ trata como co-assembly desligado.
- `cobinning_multisplit` roda independente de `coassembly.enabled` (usa assemblies indiv.).
- Nenhuma track ligada (`reads/viral/prok` todas false) e nenhuma análise de grupo ⇒ erro.

---

## 4. Estrutura de execução

```
FASTQs
  │
  ▼ QC ──────────────────────────────────────────► SEMPRE
  │
  ├─ [track reads] ──► sylph classify ──► report(Reads Survey)   [se TRACK_READS]
  │                    (não precisa de assembly)
  │
  ▼ assembly + dedup + QUAST + mapping/depth ─────► se TRACK_VIRAL ou TRACK_PROK
  │
  ├─ [track viral] (per-sample) ─────────────────► [se TRACK_VIRAL]
  │     detecção → CheckV → vRhyme → vOTU → tax → annot → defesa viral
  │
  ├─ [track prok]  (per-sample) ─────────────────► [se TRACK_PROK]
  │     MetaBAT2 + SemiBin2 + COMEBin(se GPU)      ← VAMB REMOVIDO daqui
  │        → Binette → derep(per-sample) → CheckM2 → GTDB → tax → AMR → annot
  │
  ├─ [track co-assembly] (nível de grupo) ───────► [se COASSEMBLY_ENABLED]
  │     grupos do metadata TSV
  │        MEGAHIT co-assembly(grupo) ──► contigs de grupo
  │           ├─ [viral]   ────────────────────► [se coassembly.viral]
  │           │     detecção → CheckV → vRhyme → vOTU → tax → vOTUs de grupo
  │           │     (map cada si de volta p/ abundância per-sample)
  │           └─ [binning] ────────────────────► [se coassembly.binning]
  │                 VAMB co-binning → CheckM2(grupo) → GTDB(grupo) → MAGs de grupo
  │     (análises PARALELAS, SEM derep global — lado a lado com per-sample)
  │
  ├─ [VAMB multi-split] ─────────────────────────► [se COBINNING_MULTISPLIT]
  │     concatena assemblies individuais → map todas → VAMB → split → MAGs
  │     (independe de co-assembly)
  │
  ├─ [integração] ───────────────────────────────► [se INTEGRATION_ENABLED]
  │     PHIST (fago↔seus MAGs) + arms-race defesa↔anti-defesa + abund/div cruzada
  │
  ▼ report.html ─────────────────────────────────► SEMPRE, adaptativo
```

### Fundação compartilhada

- **QC** → sempre.
- **assembly + dedup + QUAST** → só se `TRACK_VIRAL or TRACK_PROK`.
- **mapping/depth** → se viral (cobertura vRhyme) ou prok (depth de binning).

---

## 5. Track co-assembly (detalhe)

Decisão-chave: **análise paralela independente, SEM derep global.** A ideia inicial de
reconciliar MAGs per-sample + MAGs de grupo num derep global foi descartada — introduzia
acoplamento cross-sample novo (hoje o derep é per-sample) e misturava duas coisas
conceitualmente distintas. As saídas de grupo ficam lado a lado com as per-sample; ter o
mesmo genoma/vOTU nas duas é feature (permite comparar recuperação per-sample vs
co-assembly), não bug.

### A co-assembly é um "espelho de grupo" simétrico ao per-sample

A co-assembly produz um conjunto de contigs de grupo que alimenta **dois consumidores**,
espelhando o per-sample: detecção **viral** (→ vOTUs de grupo) e/ou **co-binning** VAMB
(→ MAGs de grupo). Cada um é opt-in via sub-toggle (`coassembly.viral`, `coassembly.binning`).

**Motivação do co-assembly→viral:** vírus de baixa abundância frequentemente não montam
per-sample (poucos reads → contig fragmentado); na co-assembly os reads do grupo somam e
o fago monta. Para ambientes de alta novelty (ex: Amazônia), recupera genomas virais raros
que a montagem per-sample perde.

### Fluxo

```
GRUPO g = {s1..sn}   (definido pela coluna 'group' do metadata TSV)

  reads(g) ──► MEGAHIT co-assembly(g) ──► coassembly(g).fa
                 │
                 ├─► [viral]   (se coassembly.viral)
                 │      detecção (VS2+GeNomad+VIBRANT) → CheckV → vRhyme → vOTU → tax
                 │      map cada si de volta → abundância per-sample
                 │      → vOTUs de grupo → seção "Co-assembly" no report
                 │
                 └─► [binning] (se coassembly.binning)
                        map cada si → coassembly(g) → depth multi-amostra → VAMB(g)
                        → CheckM2(g) → GTDB(g) → MAGs de grupo → seção "Co-assembly"

VAMB multi-split (se cobinning_multisplit; independe de co-assembly):
  concatena assemblies individuais (contigs renomeados por amostra)
  → map todas as amostras ao catálogo → VAMB → split bins de volta → MAGs
  (ligar junto com coassembly.binning = comparação co-assembly vs multi-split)
```

### Custos

- **Detecção viral na co-assembly dobra o custo** de VS2+GeNomad+VIBRANT para os grupos
  (opt-in via `coassembly.viral`).
- Co-assembly pode gerar **contigs quiméricos** entre cepas de fago similares de amostras
  diferentes — mitigado por CheckV (filtragem de qualidade) e pelo consenso multi-tool.

### Grouping (metadata)

- Coluna `group` num TSV de metadata (estilo nf-core/mag). A biologia dirige o
  agrupamento (mesmo sítio de coleta, tempo, tratamento).
- Reaproveita `sample_metadata.tsv` quando presente.
- Valores especiais: `all` (um grupo único), `none` (desliga).
- Agrupamento automático por similaridade (sourmash/mash) fica como opção avançada
  **fora deste corte** — anotado como trabalho futuro.

Exemplo:
```
sample          group
P01_RNG_08_947  rio_negro
P01_RNG_08_948  rio_negro
P01_RNG_3_947   solo_igapo
P01_RNG_3_948   solo_igapo
```

---

## 6. Rearranjo dos binners procarióticos

- **VAMB removido** dos binners per-sample. Justificativa: VAMB em single-sample é
  subótimo (o VAE precisa da variação de cobertura entre amostras). Dedicá-lo ao
  co-binning joga a favor da ferramenta.
- Binners per-sample passam a ser **MetaBAT2 + SemiBin2 + COMEBin**.
  - MetaBAT2 e SemiBin2: **obrigatórios**.
  - COMEBin: **obrigatório-se-GPU-disponível** (mantém o escape hatch atual para
    servidores sem GPU; não travar em ambiente sem GPU).
- **Binette** passa a consolidar 3 binners (antes 4).

---

## 7. Composição do `rule all` (funções-alvo)

Substitui o bloco chapado por funções que retornam listas de alvos:

```python
def _t_foundation():
    t = qc_targets()
    if TRACK_VIRAL or TRACK_PROK:
        t += [rep_seq_targets, quast_targets, mapping_depth_targets]
    return t

def _t_reads():        return reads_targets()      if TRACK_READS         else []
def _t_viral():        return viral_targets()      if TRACK_VIRAL         else []
def _t_prok():         return prok_targets()       if TRACK_PROK          else []
def _t_coassembly():   return coassembly_targets() if COASSEMBLY_ENABLED  else []
def _t_integration():  return integration_targets()if INTEGRATION_ENABLED else []
def _t_report():       return report_targets()     # sempre

rule all:
    input:
        *_t_foundation(), *_t_reads(), *_t_viral(),
        *_t_prok(), *_t_coassembly(), *_t_integration(), *_t_report(),
```

Includes (`rules/*.smk`) permanecem — regra cuja saída ninguém pede simplesmente não
dispara. Só o `rule all` decide o que efetivamente roda.

---

## 8. Organização de arquivos

| Arquivo | Mudança |
|---|---|
| `Snakefile` | funções-alvo `_t_*()`; parse de `tracks`/`coassembly`/`cobinning`; validações |
| `rules/prok_binning.smk` | remove VAMB dos binners per-sample; Binette consolida 3 |
| `rules/coassembly.smk` | **NOVO** — grouping do metadata, MEGAHIT co-assembly; consumidores de grupo: detecção viral (→ vOTUs) e/ou VAMB co-binning (→ MAGs, CheckM2/GTDB); VAMB multi-split |
| `rules/reads_classify.smk` | inalterado (já opcional/independente) |
| `rules/report.smk` + `scripts/report/data_loaders.py` | **aba única** "Co-assembly" (vOTUs de grupo e MAGs de grupo na mesma aba, sem sub-abas); demais abas condicionais às tracks |
| `vapor.py` (CLI) | atalhos `--track`; aliases `--until <stage>` |

---

## 9. Entry/exit points (CLI)

Aproveita `--until`/`--omit-from` nativos do Snakemake, com aliases amigáveis:

```
vapor --until qc          # para no QC
vapor --until assembly    # para no rep_seq.fasta
vapor --until viral       # para no fim da track viral
vapor --track reads       # roda só reads (override do config via --set)
vapor --track viral,prok  # roda as duas assembly-based
```

Convenção de configuração: `config.yaml` é canônico; `vapor --set KEY=VALUE` sobrepõe
(já existe). Os atalhos `--track` traduzem para overrides de `--config` do Snakemake.

---

## 10. Fora de escopo (trabalho futuro)

- Agrupamento automático de amostras por similaridade (sourmash/mash).
- Migração para o sistema de `module` do Snakemake 8.
- Derep global reconciliando MAGs per-sample + de grupo.
- Reconciliação/cross-comparação automática de vOTUs de grupo vs per-sample (ficam lado a lado).
- Kraken2/Bracken/Phanta — decidido rodar **fora** do pipeline, como análise standalone
  (ver `scripts/install_kraken_phanta.sh`).
