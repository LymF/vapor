# VAPOR — Report Visualization Guide

The chart standard for the VAPOR HTML report (`scripts/report/`). It is the visual
counterpart to `docs/PIPELINE_METHODS.md`: that file fixes the *formulas*, this one
fixes the *forms*.

**Scope.** VAPOR is a general-purpose pipeline, so every rule here is keyed to the
**shape of the data** (how many samples, how many categories, how many points) and
never to the biology of one dataset. A rule that needs to know "these samples are
rivers" does not belong in this guide — grouping semantics live in the co-assembly
track (`GROUPS`), not in chart logic.

**Maintenance.** When a chart is added or changed, check it against §3 (form), §4
(triggers) and §7 (anti-patterns). If a chart matches an anti-pattern, it is wrong.

---

## 1. The procedure

Always in this order. **Color comes last** — most bad charts pick color first.

1. **Pick the form** from the data's job (§3). Sometimes the answer is not a chart.
2. **Assign color by the job it does** (§5): categorical, sequential, diverging, status.
3. **Validate the palette** — run the checker, never eyeball CVD (§6).
4. **Apply mark specs**: thin marks, 2px gap between adjacent fills, ≥8px markers,
   recessive grid and axes, selective direct labels (never a number on every point).
5. **Add the hover layer** — a tooltip is the default, not a bonus. Only a bare stat
   tile skips it.
6. **Accessibility pass**: legend present for ≥2 series; identity never by color
   alone; a table view exists; dark mode is *designed*, not an inverted flip.
7. **Render it and look at it.** The validator checks color, not layout.

---

## 2. Is it even a chart?

| The data is… | Use | Not |
|---|---|---|
| A single current value | **stat tile** / hero number | a one-bar bar chart |
| A handful of headline numbers | **KPI row** of stat tiles | a grouped bar chart |
| A single ratio against a limit | **meter** | a 2-slice pie |
| More than ~7 classes that all carry meaning | a **table** (or table + chart) | more colors |

---

## 3. Data shape → form

The core table. "n" is the number of items actually plotted.

| Data shape | Form | Notes |
|---|---|---|
| **One count per sample** (contigs, vOTUs, MAGs) | sorted bar; **lollipop** when labels are long | past the sample trigger (§4) switch to a heatmap row or fold to top-N |
| **Several metrics × samples** | **heatmap**, sequential, normalized *per metric* | one column per metric; never normalize across metrics of different units |
| **Part-to-whole** (quality tiers, domain, composition at a rank) | **stacked bar** — 100% when proportion is the point, absolute when magnitude is | horizontal when category names are long; ≤8 categories, rest to "Other" |
| **Distribution of one continuous variable** (completeness, GC, length, depth) | **histogram / density**; **dot / strip** at small n; **ridgeline** across many samples | see the n triggers in §4 — this is where boxplots quietly lie |
| **Two continuous variables with thresholds** (completeness × contamination) | **scatter + threshold zones**; **hexbin / 2-D density** once it saturates | always draw the cutoff lines; they are the reason the chart exists |
| **Set membership / tool agreement** | **UpSet** | for ≥3 sets. Venn is unreadable past 3 and worse even at 3 |
| **Sample × feature matrix** (abundance, gene presence) | **clustered heatmap** with dendrogram; **bubble** when sparse | sparse ≈ <20% of cells filled |
| **Hierarchy** (full taxonomy) | **sunburst / treemap** for *one* unit; **stacked bar at a fixed rank** to *compare* units | a sunburst cannot compare samples side by side |
| **Sequential attrition** (reads → trimmed → assembled → vOTUs) | **funnel** or waterfall | only for genuinely ordered stages with loss |
| **Genomic coordinates** (defense islands, genome maps) | **track / browser**: bp ruler, strand-aware arrows, true feature widths | never plot gene *order* when coordinates exist |
| **Pairwise relations** (gene sharing, phage–host) | node-link when sparse; **adjacency matrix** when dense; bipartite for two entity types | past ~150 nodes a node-link is a hairball — go to the matrix |
| **Ordination** (PCoA / NMDS) | scatter, **% variance on both axes**, ellipses when groups exist | axes without variance % are incomplete |
| **Effort vs discovery** | **accumulation / rarefaction curve** | needs a feature space shared across samples. **Desde o catálogo global de vOTUs (2026-08-18) esse espaço existe para toda a rodada**, não só na co-assembly: o catálogo clusteriza uma vez sobre todas as amostras e grupos, então a curva global é legítima. Antes disso os vOTUs per-sample eram clusterizados independentemente e acumulá-los contava o mesmo vírus várias vezes. Sempre promediar sobre ordens aleatórias — uma ordenação é arbitrária |
| **Ranked list with long names** | **horizontal bar**, top-N + "Other" | horizontal whenever labels exceed ~12 characters |

---

## 4. Numeric triggers

The rules that make the report adapt to any dataset. Implement them as thresholds,
never as per-dataset special cases.

All thresholds live in one object, `window.VIZ` (`app.js`) — change one there and
every chart built on the helpers follows.

| Trigger | Threshold | Behavior |
|---|---|---|
| **Samples on a category axis** | n > 12 | the axis flips **horizontal** and rows sort by total, so names read straight and every sample stays visible (`samplesBar`) |
| **Distribution, few points** | n < 20 | **dot / strip plot** — a KDE over few points invents structure that is not in the data |
| **Distribution, enough points** | n ≥ 20 | histogram or density |
| **Distributions across samples** | samples > 8 | **ridgeline** instead of a row of boxplots |
| **Scatter overplotting** | n > ~500 | **hexbin** or 2-D density instead of raw points |
| **Categorical series** | > 8 | fold the tail into "Other" (`window.foldOther`), facet, or use composite encoding — **never generate a 9th hue** |
| **Network nodes** | > ~150 | adjacency matrix instead of node-link |
| **Table rows** | > ~200 | paginate + search (already the `makeTable` default) |

**Why the n<20 rule matters.** A density curve asserts a continuous underlying
distribution. Drawn over 3 MAGs it is a fabrication — the reader sees a smooth,
confident shape produced by the bandwidth, not by the data. Show the points.

---

## 5. Color

Four jobs, one rule each:

| Job | Palette | Rule |
|---|---|---|
| **Identity** — the series *are* the subject | categorical | fixed order, never cycled |
| **Magnitude** — more vs less | sequential | one hue, light → dark. Never a rainbow |
| **Polarity** — above/below a reference | diverging | two hues + a **neutral gray** midpoint |
| **State** — good / warning / serious / critical | status | reserved; never reused as "series 4"; always with an icon or label |

**Series-count ladder**

| Series | Treatment |
|---|---|
| 1–3 | color alone is comfortable; direct-label |
| 4 | direct labels become mandatory |
| 5–6 | legend, or small multiples |
| 7–8 | ceiling — beyond it, fold to "Other" or facet |

**Non-negotiables**

- One y-axis. Never dual-axis — two scales on one frame is the single most common
  chart error. Two measures → two charts or index them to a common base.
- Color follows the entity, not its rank: filtering must not repaint the survivors.
- Text wears text tokens, never the series color. A colored mark beside the label
  carries the identity.
- Quality tiers (Complete / HQ / MQ / LQ / Not-determined) are an **ordinal ladder**,
  so they take an ordered good→bad ramp, not eight categorical hues.

---

## 6. The validated palette

`window.PAL` in `scripts/report/components/app.js`, checked with the palette
validator in **both** modes:

```
#0d9488  #d97706  #7c3aed  #0891b2  #16a34a  #db2777  #9333ea  #ef4444
```

`window.PAL_MUTED = #64748b` is the neutral for "Other"/"Unknown" — an escape hatch
for the tail, never a categorical identity slot.

**Validation record**

- Light: PASS. Dark: PASS.
- One standing warning: `#db2777` ↔ `#16a34a` at CVD ΔE 6.1 (deutan), inside the
  6–8 floor band. This is legal **only** because every categorical chart also
  carries a legend or direct labels plus 2px gaps between fills. If a future chart
  drops that secondary encoding, this pair must be re-stepped.
- Slot 6 was `#f59e0b` and was replaced: it FAILED the dark-mode lightness band
  (L 0.769) and sat at 2.09:1 contrast against the light surface.
- A 9–15 slot extension was attempted and FAILED (two slots below the chroma floor,
  one adjacent pair at ΔE 2.9). Eight is the ceiling; use `foldOther`.

**Re-run before changing any color:**

```bash
node scripts/validate_palette.js "<hex,hex,…>" --mode light
node scripts/validate_palette.js "<hex,hex,…>" --mode dark
```

Never reason about ΔE by eye — the check is computable, so compute it.

---

## 7. Anti-patterns

If a chart matches a row here, it is wrong.

| Anti-pattern | Why it fails | Use instead |
|---|---|---|
| Pie chart of taxonomy | angles do not compare, and taxonomy always exceeds ~7 classes | sorted horizontal bar |
| Boxplot of completeness | **hides bimodality** — the classic boxplot failure; two very different distributions can share a box | histogram / density |
| A mean plotted as a bar (e.g. mean GC per sample) | discards the entire distribution, which is often the diagnostic signal | density / ridgeline |
| Sunburst used to compare samples | a hierarchy cannot be read side by side | stacked bar at a fixed rank |
| Saturated scatter | overplotting hides where the mass is | hexbin / 2-D density |
| Rainbow on a continuous scale | not perceptually ordered; invents boundaries | single-hue sequential |
| Dual y-axis | any correlation it shows is an artifact of scaling | two charts |
| Colour as the only encoding | fails for CVD readers, print and forced-colors | add label, shape or texture |
| A 9th generated hue | indistinguishable under CVD; breaks every check | "Other", facet, composite encoding |
| Gene *order* drawn as a genome map | fake spacing, no strand, uniform widths | real bp track (coordinates come free in the Prodigal `.faa` header) |
| Zero-baseline violation on bars | exaggerates differences | bars always start at zero (lines need not) |

---

## 8. Implementation notes

Shared helpers in `scripts/report/components/app.js` — use them rather than
re-implementing:

| Helper | Purpose |
|---|---|
| `mkChart(id, option)` | ECharts wrapper: registers the instance, handles resize and theme. Honours `option.__height` for forms whose height depends on the data |
| `VIZ` | the numeric triggers of §4, in one object |
| `samplesBar({samples, series, stack, valueName, sort})` | one value per sample, N series. Vertical bars while small; **horizontal + sorted** past `VIZ.manySamples` |
| `distPlot({groups, xName, log, cutoffs, xMin, xMax, colors})` | distribution per unit. Picks **strip plot** (median n < `VIZ.densityMinN`), **density**, or **ridgeline** (groups > `VIZ.manyGroups`); `cutoffs` draws threshold lines |
| `upsetPlot({sets, combos, valueName})` | set-intersection sizes + membership matrix — the readable form for tool/detector agreement |
| `hexbin(pts, {threshold, cols})` | hex-grid density binning for saturated scatters; returns `null` below `VIZ.denseScatter` so the caller keeps drawing points |
| `load_votu_accumulation()` (Python) | per-group collector curve from the co-binning abundance matrix + vOTU clusters; permutation mean + 10–90 band |
| `makeTable(id, rows, cols, opts)` | paginated + searchable table (the "table view" the a11y pass requires) |
| `boxStats(values)` | quartiles + outliers |
| `foldOther(map, max)` | folds a `{name: count}` tail into "Other" — the ≤8 series rule |
| `qualBadge(q)` | quality tier badge with the ordinal ramp |
| `PAL`, `PAL_MUTED` | the validated categorical palette |
| `TRIGGERS` (`src/report-ui/src/viz/triggers.js`) | os gatilhos do §4 no report v2, aplicados DENTRO de cada forma em vez de dependerem de o autor chamar o helper certo |
| `foldOther` (`src/report-ui/src/viz/palette.js`) | versão v2 do dobrador de cauda; "Other" sempre em último |
| `useResize` (`src/report-ui/src/viz/useResize.js`) | dimensões do container para gráficos SVG responsivos |

Every chart lives in a component under `scripts/report/components/` and is assembled
by `renderer.py`; data loading is `data_loaders.py`. Charts must degrade to an empty
state (never a broken axis) when their input is missing, because any track can be
switched off in `config.yaml`.
