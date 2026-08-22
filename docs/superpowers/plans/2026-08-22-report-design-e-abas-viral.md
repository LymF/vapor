# Report React + D3 — Design system e abas Sequenciamento/Viral (plano 2 de 4)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** dar ao report v2 um design system de verdade e as duas primeiras abas completas — Sequenciamento e Catálogo viral, incluindo o explorador de vOTU — com as formas de gráfico que elas exigem.

**Architecture:** fatia vertical. Cada forma nova nasce com o painel que a consome, nunca antes. O design system vem primeiro porque tudo depois se apoia nele. Os gatilhos numéricos do guia moram dentro das formas (`viz/triggers.js`), não no chamador.

**Tech Stack:** React 18, D3 v7, esbuild, vitest + @testing-library/react, pytest.

**Spec:** `docs/superpowers/specs/2026-08-22-report-react-d3-design.md`

## Global Constraints

- **O report é UM arquivo HTML standalone.** Sem CDN, sem fetch, sem servidor.
- **Node nunca é dependência de runtime.** O bundle compilado é versionado; recompile e commite junto com qualquer mudança em `src/report-ui/`.
- **Orçamento de payload: 25 MB**, checado por `scripts/report/schema.py` antes de escrever.
- **`docs/REPORT_VIZ_GUIDE.md` é lei.** Forma pelo shape do dado; a paleta de 8 cores é fixa e validada sob daltonismo — nunca uma nona cor, nunca reordenar.
- **Nunca só cor.** Toda identidade carrega rótulo, glifo ou posição além da cor.
- **Ferramenta que falhou é lacuna, nunca zero.** Os quatro estados (`ok`/`skipped`/`failed`/`unknown`) permanecem distintos em qualquer componente que os mostre.
- **Todo painel degrada para estado vazio**, jamais para eixo quebrado: qualquer trilha pode ser desligada em `config.yaml`.
- **Zero real ≠ fonte ausente.** Contagem sem fonte é `None` e some; zero medido é desenhado como zero.
- **Nenhum gráfico taxonômico preso a um rank.** Rank é estado, trocável no lugar.
- **Proibido treemap** (decisão do usuário).
- **Commits sem rodapé de coautoria.**
- Testes: `python -m pytest tests/ -q` (env `snakemake`) e `npx vitest run` em `src/report-ui/` (env `./.envs/reportui`).

---

### Task 1: Design system

**Files:**
- Modify: `src/report-ui/src/styles.css`
- Create: `src/report-ui/src/viz/theme.js`
- Create: `src/report-ui/test/theme.test.jsx`
- Modify: `src/report-ui/src/App.jsx`

**Interfaces:**
- Consumes: nada.
- Produces:
  - `styles.css` com escala tipográfica, espaçamento, superfícies, sombras, raios, e **modo escuro** por `@media (prefers-color-scheme: dark)` sobre os mesmos tokens.
  - `useTheme()` → `{ theme: 'light'|'dark', toggle() }`, persistido em `localStorage` com try/catch (pode lançar em contexto sem storage).
  - `<ThemeToggle />`.

- [ ] **Step 1: Escrever o teste que falha**

`src/report-ui/test/theme.test.jsx`:

```jsx
import { render, screen, fireEvent } from '@testing-library/react';
import { useTheme, ThemeToggle } from '../src/viz/theme.js';

function Sonda() {
  const { theme } = useTheme();
  return <span data-testid="t">{theme}</span>;
}

beforeEach(() => { document.documentElement.removeAttribute('data-theme'); });

test('comeca em light quando nada foi escolhido', () => {
  render(<Sonda />);
  expect(screen.getByTestId('t').textContent).toBe('light');
});

test('o toggle troca o tema e marca a raiz', () => {
  render(<><Sonda /><ThemeToggle /></>);
  fireEvent.click(screen.getByRole('button', { name: /tema/i }));
  expect(document.documentElement.getAttribute('data-theme')).toBe('dark');
});

test('nao quebra quando localStorage lanca', () => {
  const orig = Object.getOwnPropertyDescriptor(window, 'localStorage');
  Object.defineProperty(window, 'localStorage', {
    configurable: true,
    get() { throw new Error('bloqueado'); },
  });
  expect(() => render(<Sonda />)).not.toThrow();
  if (orig) Object.defineProperty(window, 'localStorage', orig);
});
```

- [ ] **Step 2: Rodar e confirmar que falha**

Run: `npx vitest run test/theme.test.jsx` — FAIL: `Failed to resolve import "../src/viz/theme.js"`

- [ ] **Step 3: Implementar o tema**

`src/report-ui/src/viz/theme.js`:

```jsx
import { useCallback, useEffect, useState } from 'react';

const CHAVE = 'vapor-report-theme';

function lido() {
  // localStorage lanca em contexto sem storage (captura de thumbnail, aba
  // privada, navegador com dados de site bloqueados) -- ler tem de ser seguro.
  try {
    const v = window.localStorage.getItem(CHAVE);
    return v === 'dark' || v === 'light' ? v : null;
  } catch { return null; }
}

function grava(v) {
  try { window.localStorage.setItem(CHAVE, v); } catch { /* sem storage */ }
}

export function useTheme() {
  const [theme, setTheme] = useState(() => lido() ?? 'light');

  useEffect(() => {
    document.documentElement.setAttribute('data-theme', theme);
  }, [theme]);

  const toggle = useCallback(() => {
    setTheme((t) => {
      const novo = t === 'dark' ? 'light' : 'dark';
      grava(novo);
      return novo;
    });
  }, []);

  return { theme, toggle };
}

export function ThemeToggle() {
  const { theme, toggle } = useTheme();
  return (
    <button className="theme-toggle" onClick={toggle}
            aria-label={`Alternar tema (atual: ${theme === 'dark' ? 'escuro' : 'claro'})`}>
      {theme === 'dark' ? '☾' : '☀'}
    </button>
  );
}
```

- [ ] **Step 4: Reescrever o CSS como design system**

Substitua `src/report-ui/src/styles.css` inteiro. Requisitos, não decoração:

- **Tokens em `:root`** para a paleta clara COMPLETA: `--bg`, `--surface`, `--surface-2`, `--border`, `--border-strong`, `--text`, `--text-2`, `--text-3`, mais os quatro de estado (`--ok`, `--skipped`, `--failed`, `--unknown`) e os oito categóricos (`--cat-1`..`--cat-8`) com os valores de `PAL`.
- **Escala tipográfica** (`--fs-xs` 0.75rem, `--fs-sm` 0.85rem, `--fs-base` 1rem, `--fs-lg` 1.15rem, `--fs-xl` 1.4rem, `--fs-2xl` 1.9rem) e **escala de espaçamento** (`--sp-1` .25rem … `--sp-8` 2rem). Nada de valor mágico fora dos tokens.
- **Modo escuro** redefinindo SÓ os tokens, em dois seletores: `@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) { … } }` e `:root[data-theme="dark"] { … }`. Nenhuma cor pode ter sua única definição dentro de um bloco de media.
- **`body` com background explícito** do token.
- Componentes: `.nav` (fixa, com sombra ao rolar), `.card` (superfície, borda, raio, título e subtítulo), `.kpi-row`/`.stat-tile`, `.status-matrix`, `.empty`, `.chart` (container responsivo), `.legend`, `.tooltip` (posicionado, `pointer-events:none`), `.table` (zebra, cabeçalho fixo, scroll horizontal próprio).
- **Conteúdo largo rola dentro do próprio container** (`overflow-x:auto`); o `body` nunca rola na horizontal.
- `@media (max-width: 820px)`: nav vira coluna, `.kpi-row` empilha.

- [ ] **Step 5: Ligar o toggle ao shell**

Em `src/report-ui/src/App.jsx`, importe `ThemeToggle` e monte-o dentro de `.nav`, depois do filtro de amostra.

- [ ] **Step 6: Rodar a suíte inteira, compilar, commitar**

```bash
npx vitest run                     # 23 anteriores + 3 novos = 26
npm run build
git add src/report-ui scripts/report/assets/report-ui.js scripts/report/assets/report-ui.css
git commit -m "feat(report-v2): design system com tokens, escala e modo escuro"
```

---

### Task 2: Primitivas de gráfico

**Files:**
- Create: `src/report-ui/src/viz/Chart.jsx`
- Create: `src/report-ui/src/viz/Axis.jsx`
- Create: `src/report-ui/src/viz/Tooltip.jsx`
- Create: `src/report-ui/src/viz/Legend.jsx`
- Create: `src/report-ui/test/primitives.test.jsx`

**Interfaces:**
- Consumes: `useResize`, `PAL`.
- Produces:
  - `<Chart title sub height children/>` — mede a largura, passa `{width, height}` via render-prop, desenha `<figure>` + `<figcaption>` e o estado vazio quando `empty` é verdadeiro.
  - `<AxisBottom scale width height tickFormat/>` e `<AxisLeft scale …/>` — grade recessiva, ticks legíveis, sem linha de domínio grossa.
  - `useTooltip()` → `{ show(evt, conteudo), hide(), node }`.
  - `<Legend items={[{label, color}]} />`.

- [ ] **Step 1: Escrever o teste que falha**

`src/report-ui/test/primitives.test.jsx`:

```jsx
import { render, screen } from '@testing-library/react';
import { Chart } from '../src/viz/Chart.jsx';
import { Legend } from '../src/viz/Legend.jsx';

test('Chart mostra titulo e passa dimensoes ao filho', () => {
  render(<Chart title="Contigs" sub="por amostra">{({ width, height }) =>
    <text data-testid="dim">{`${width}x${height}`}</text>}</Chart>);
  expect(screen.getByText('Contigs')).toBeTruthy();
  expect(screen.getByTestId('dim').textContent).toMatch(/^\d+x\d+$/);
});

test('Chart em estado vazio nao chama o filho', () => {
  let chamou = false;
  render(<Chart title="X" empty emptyLabel="Sem dados nesta rodada">{() => { chamou = true; return null; }}</Chart>);
  expect(chamou).toBe(false);
  expect(screen.getByText('Sem dados nesta rodada')).toBeTruthy();
});

test('Legend carrega rotulo alem da cor', () => {
  render(<Legend items={[{ label: 'Bacteria', color: '#0d9488' }]} />);
  expect(screen.getByText('Bacteria')).toBeTruthy();
});
```

- [ ] **Step 2: Rodar e confirmar que falha**

Run: `npx vitest run test/primitives.test.jsx` — FAIL na resolução de `../src/viz/Chart.jsx`.

- [ ] **Step 3: Implementar**

`Chart.jsx` mede com `useResize`, usa `height` (default 280) e chama `children({width, height})` dentro de um `<svg role="img">`; com `empty` verdadeiro renderiza `<p className="empty">{emptyLabel}</p>` e **não** chama `children`. `Axis.jsx` usa `scale.ticks()`/`scale.domain()` conforme o tipo e desenha linhas de grade com `stroke: var(--border)`. `Tooltip.jsx` mantém um nó absoluto único, posicionado por `clientX/clientY`, com `pointer-events:none`. `Legend.jsx` renderiza `<ul>` com um quadrado colorido e o texto do rótulo.

- [ ] **Step 4: Passar, compilar, commitar**

```bash
npx vitest run
npm run build
git add -A && git commit -m "feat(report-v2): primitivas de grafico (Chart, Axis, Tooltip, Legend)"
```

---

### Task 3: Barras e distribuições

**Files:**
- Create: `src/report-ui/src/charts/BarChart.jsx`
- Create: `src/report-ui/src/charts/DistPlot.jsx`
- Create: `src/report-ui/test/bar-dist.test.jsx`

**Interfaces:**
- Consumes: `Chart`, `Axis`, `useTooltip`, `PAL`, `PAL_MUTED`, `TRIGGERS`, `pickAxisOrientation`, `pickDistributionForm`, `foldOther`.
- Produces:
  - `<BarChart data={[{name, value}]} orientation="auto"|"vertical"|"horizontal" sort="desc"|"none" valueName maxSeries/>` — em `auto`, vira horizontal acima de `TRIGGERS.manySamples` **ou** quando algum rótulo passa de 12 caracteres. Barras sempre a partir do zero.
  - `<StackedBar data={[{name, parts: {categoria: valor}}]} order={[...]} normalize={false} colors/>` — 100% quando `normalize`.
  - `<DistPlot groups={[{name, values: []}]} xName log cutoffs/>` — escolhe strip / histograma / ridgeline pelos gatilhos; `cutoffs: [{value, label}]` desenha linha tracejada com rótulo.

- [ ] **Step 1: Escrever os testes que falham**

```jsx
import { render, screen } from '@testing-library/react';
import { BarChart, StackedBar } from '../src/charts/BarChart.jsx';
import { DistPlot } from '../src/charts/DistPlot.jsx';

const muitas = Array.from({ length: 15 }, (_, i) => ({ name: `S${i}`, value: i + 1 }));

test('poucas amostras ficam na vertical', () => {
  render(<BarChart data={[{ name: 'a', value: 1 }]} valueName="n" />);
  expect(screen.getByTestId('barchart').getAttribute('data-orientation')).toBe('vertical');
});

test('muitas amostras viram horizontal', () => {
  render(<BarChart data={muitas} valueName="n" />);
  expect(screen.getByTestId('barchart').getAttribute('data-orientation')).toBe('horizontal');
});

test('rotulo longo tambem vira horizontal', () => {
  render(<BarChart data={[{ name: 'Pseudomonadota_muito_longo', value: 3 }]} valueName="n" />);
  expect(screen.getByTestId('barchart').getAttribute('data-orientation')).toBe('horizontal');
});

test('a barra parte do zero', () => {
  render(<BarChart data={[{ name: 'a', value: 5 }, { name: 'b', value: 10 }]} valueName="n" />);
  const barras = [...document.querySelectorAll('[data-bar]')];
  const razao = Number(barras[1].getAttribute('data-len')) / Number(barras[0].getAttribute('data-len'));
  expect(razao).toBeCloseTo(2, 1);   // zero-baseline: 10 e o dobro de 5
});

test('mais de 8 categorias dobram em Other', () => {
  const partes = Object.fromEntries(Array.from({ length: 12 }, (_, i) => [`c${i}`, 1]));
  render(<StackedBar data={[{ name: 'S1', parts: partes }]} />);
  expect(screen.getByTestId('stacked').getAttribute('data-series')).toBe('8');
});

test('poucos pontos viram strip, muitos viram densidade', () => {
  const poucos = { name: 'g', values: [1, 2, 3] };
  const muitos = { name: 'g', values: Array.from({ length: 40 }, (_, i) => i) };
  const { rerender } = render(<DistPlot groups={[poucos]} xName="x" />);
  expect(screen.getByTestId('distplot').getAttribute('data-form')).toBe('strip');
  rerender(<DistPlot groups={[muitos]} xName="x" />);
  expect(screen.getByTestId('distplot').getAttribute('data-form')).toBe('density');
});

test('muitos grupos viram ridgeline', () => {
  const g = Array.from({ length: 10 }, (_, i) =>
    ({ name: `S${i}`, values: Array.from({ length: 40 }, (_, j) => j) }));
  render(<DistPlot groups={g} xName="x" />);
  expect(screen.getByTestId('distplot').getAttribute('data-form')).toBe('ridgeline');
});

test('linhas de corte aparecem com rotulo', () => {
  render(<DistPlot groups={[{ name: 'g', values: [1, 2, 3] }]} xName="x"
                   cutoffs={[{ value: 2, label: 'MQ' }]} />);
  expect(screen.getByText('MQ')).toBeTruthy();
});
```

- [ ] **Step 2: Rodar e confirmar que falham**
- [ ] **Step 3: Implementar as duas formas**

`BarChart` marca `data-testid="barchart"` e `data-orientation`; cada barra leva `data-bar` e `data-len` (o comprimento em px, para o teste do zero-baseline). `StackedBar` aplica `foldOther` às categorias e marca `data-series` com a contagem final. `DistPlot` marca `data-form` com `strip`/`density`/`ridgeline`, escolhendo pela **mediana** de `values.length` entre os grupos e por `groups.length > TRIGGERS.manyGroups`; densidade por KDE gaussiano com regra de Silverman.

- [ ] **Step 4: Passar, compilar, commitar**

---

### Task 4: Heatmap e UpSet

**Files:**
- Create: `src/report-ui/src/charts/Heatmap.jsx`
- Create: `src/report-ui/src/charts/UpSet.jsx`
- Create: `src/report-ui/test/heatmap-upset.test.jsx`

**Interfaces:**
- Produces:
  - `<Heatmap rows cols values={{[row]: {[col]: n}}} normalize="per-col"|"none" sparseAsBubble/>` — sequencial de um matiz; **normaliza por coluna** quando as colunas têm unidades diferentes; abaixo de 20% de células preenchidas vira bolha (área ∝ valor).
  - `<UpSet sets={{nome: tamanho}} combos={[{tools: [...], count}]} valueName/>` — barras de interseção ordenadas + matriz de pertencimento.

- [ ] **Step 1: Testes**

```jsx
test('matriz esparsa vira bolha', () => {
  render(<Heatmap rows={['a','b','c','d','e']} cols={['x','y','z','w']}
                  values={{ a: { x: 1 }, b: { y: 2 } }} sparseAsBubble />);
  expect(screen.getByTestId('heatmap').getAttribute('data-mode')).toBe('bubble');
});

test('matriz densa fica heatmap', () => {
  const vals = Object.fromEntries(['a','b'].map(r => [r, { x: 1, y: 2 }]));
  render(<Heatmap rows={['a','b']} cols={['x','y']} values={vals} sparseAsBubble />);
  expect(screen.getByTestId('heatmap').getAttribute('data-mode')).toBe('heat');
});

test('normalizacao por coluna nao mistura unidades', () => {
  render(<Heatmap rows={['a','b']} cols={['n50','contigs']}
                  values={{ a: { n50: 1000, contigs: 10 }, b: { n50: 2000, contigs: 5 } }}
                  normalize="per-col" />);
  const cel = (r, c) => document.querySelector(`[data-cell="${r}|${c}"]`).getAttribute('data-norm');
  expect(Number(cel('b','n50'))).toBeCloseTo(1, 2);
  expect(Number(cel('a','contigs'))).toBeCloseTo(1, 2);
});

test('UpSet ordena as intersecoes por tamanho', () => {
  render(<UpSet sets={{ geNomad: 10, VirSorter2: 8 }}
                combos={[{ tools: ['geNomad'], count: 3 },
                         { tools: ['geNomad','VirSorter2'], count: 7 }]} valueName="contigs" />);
  const barras = [...document.querySelectorAll('[data-combo]')].map(b => b.getAttribute('data-combo'));
  expect(barras[0]).toBe('geNomad+VirSorter2');
});
```

- [ ] **Step 2: Falhar** → **Step 3: Implementar** → **Step 4: Passar, compilar, commitar**

---

### Task 5: Sunburst e seletor de rank

**Files:**
- Create: `src/report-ui/src/charts/Sunburst.jsx`
- Create: `src/report-ui/src/viz/RankSelector.jsx`
- Modify: `src/report-ui/src/state/store.jsx`
- Create: `src/report-ui/test/taxonomy.test.jsx`

**Interfaces:**
- Produces:
  - `RANKS = ['Phylum','Class','Order','Family','Genus']` e `RANK_LABEL = {Phylum:'Filo', Class:'Classe', Order:'Ordem', Family:'Família', Genus:'Gênero'}` — as chaves são as colunas do dado (inglês, como saem do MMseqs2 e do GTDB-Tk), os rótulos são o que o usuário lê. Nunca traduza a chave.
  - No store: `rank`/`setRank` (default `'Family'`) e `taxonFilter`/`setTaxonFilter` (`{rank, name}` ou `null`).
  - `<RankSelector />` — controle segmentado; muda o rank global.
  - `<Sunburst rows={[{Phylum, Class, ..., count}]} onDrill={(rank, name) => void} />` — clique num arco **desce o rank global e fixa o filtro de táxon**.

- [ ] **Step 1: Escrever os testes que falham**

`src/report-ui/test/taxonomy.test.jsx`. Note os dois auxiliares no topo: `Prov`
monta o provider (o rank é estado global, não prop) e `Espia` expõe esse estado
para o teste ler.

```jsx
import { render, screen, fireEvent } from '@testing-library/react';
import { ReportProvider, useReport } from '../src/state/store.jsx';
import { RankSelector, RANKS, RANK_LABEL } from '../src/viz/RankSelector.jsx';
import { Sunburst } from '../src/charts/Sunburst.jsx';

function Prov({ children }) {
  return <ReportProvider data={{ run: { title: 'VAPOR', samples: [] } }}>{children}</ReportProvider>;
}

function Espia() {
  const { rank, taxonFilter } = useReport();
  return (<>
    <span data-testid="rank">{rank}</span>
    <span data-testid="taxon">{taxonFilter?.name ?? ''}</span>
  </>);
}

test('os rotulos sao portugues sobre as chaves em ingles', () => {
  expect(RANKS).toEqual(['Phylum', 'Class', 'Order', 'Family', 'Genus']);
  expect(RANK_LABEL.Family).toBe('Família');
  expect(RANK_LABEL.Phylum).toBe('Filo');
});

test('o seletor troca o rank global', () => {
  render(<Prov><RankSelector /><Espia /></Prov>);
  fireEvent.click(screen.getByRole('button', { name: 'Família' }));
  expect(screen.getByTestId('rank').textContent).toBe('Family');
});

test('clicar num arco desce o rank e fixa o taxon', () => {
  const linhas = [{ Phylum: 'Uroviricota', Class: 'Caudoviricetes', count: 5 }];
  render(<Prov><Sunburst rows={linhas} /><Espia /></Prov>);
  fireEvent.click(screen.getByTestId('arc-Uroviricota'));
  expect(screen.getByTestId('rank').textContent).toBe('Class');
  expect(screen.getByTestId('taxon').textContent).toBe('Uroviricota');
});

test('rank mais fundo que o dado nao inventa nivel', () => {
  const linhas = [{ Phylum: 'Uroviricota', count: 2 }];
  render(<Prov><Sunburst rows={linhas} /></Prov>);
  expect(screen.queryByTestId('arc-undefined')).toBeNull();
});
```

- [ ] **Step 2: Rodar e confirmar que falham**
- [ ] **Step 3: Implementar** `RankSelector.jsx`, `Sunburst.jsx` e as duas fatias novas do store (`rank`, `taxonFilter`), mantendo `sample`/`tab` intactos.
- [ ] **Step 4: Rodar a suíte inteira, compilar e commitar.**

**Nenhum treemap** — decisão do usuário, registrada no spec.

---

### Task 6: Trilha genômica

**Files:**
- Create: `src/report-ui/src/charts/GenomeTrack.jsx`
- Create: `src/report-ui/test/genometrack.test.jsx`

**Interfaces:**
- Produces: `<GenomeTrack length features={[{start, end, strand, label, kind}]} kinds={{kind: {color, label}}} />` — régua em bp, seta por fita, **largura real da feature**, colisão resolvida em faixas.

- [ ] **Step 1: Testes**

```jsx
const feats = [
  { start: 1, end: 900, strand: '+', label: 'terminase', kind: 'phrog' },
  { start: 1200, end: 1500, strand: '-', label: 'anti-CBASS', kind: 'antidefense' },
];

test('a largura da feature e proporcional ao seu tamanho em bp', () => {
  render(<GenomeTrack length={3000} features={feats} kinds={{}} />);
  const a = Number(document.querySelector('[data-feature="terminase"]').getAttribute('data-w'));
  const b = Number(document.querySelector('[data-feature="anti-CBASS"]').getAttribute('data-w'));
  expect(a / b).toBeCloseTo(900 / 300, 1);
});

test('a fita determina o sentido da seta', () => {
  render(<GenomeTrack length={3000} features={feats} kinds={{}} />);
  expect(document.querySelector('[data-feature="anti-CBASS"]').getAttribute('data-strand')).toBe('-');
});

test('features sobrepostas vao para faixas diferentes', () => {
  const sobre = [{ start: 1, end: 1000, strand: '+', label: 'a', kind: 'x' },
                 { start: 500, end: 1500, strand: '+', label: 'b', kind: 'x' }];
  render(<GenomeTrack length={2000} features={sobre} kinds={{}} />);
  const fa = document.querySelector('[data-feature="a"]').getAttribute('data-lane');
  const fb = document.querySelector('[data-feature="b"]').getAttribute('data-lane');
  expect(fa).not.toBe(fb);
});

test('sem features desenha so a regua', () => {
  render(<GenomeTrack length={1000} features={[]} kinds={{}} />);
  expect(document.querySelectorAll('[data-feature]').length).toBe(0);
});
```

- [ ] Steps 2-4 como nas anteriores.

---

### Task 7: Dados das duas abas

**Files:**
- Modify: `scripts/report/renderer_v2.py`
- Modify: `tests/test_renderer_v2.py`

**Interfaces:**
- Consumes de `data_loaders.py` (já existem, não reescreva): `parse_fastp_json`, `parse_quast_all`, `parse_mapping_rate`, `collect_depth_data`, `parse_fasta_lengths`, `parse_tsv` (CheckV), `parse_support_combos`, `load_viral_taxonomy`, `load_votu_catalog`, `load_votu_presence`, `load_votu_lifestyle`, `load_phrogs`, `load_putative_amgs`.
- Produces, dentro de `build_data`:
  - `data.sequencing = {qc: [...], quast: {sample: {metrica: valor}}, mapping: [...], lengths: {sample: [int]}, depth: {sample: [[len, depth]]}}`
  - `data.viral = {taxonomy: [{sample, rank fields..., count}], checkv_tiers: {sample: {tier: n}}, detectors: {sets, combos}, catalog: {n_votus, n_pool, reduction_pct}, presence: [...], lifestyle: {...}, explorer: [{votu_id, length, features: [...]}]}`

**Regras que os testes travam:**
1. Toda tabela embarcada passa por `schema.project` com os campos declarados — nada de linha crua de ferramenta.
2. Fonte ausente → chave ausente (não `[]` nem `0`), para o painel distinguir "desligado" de "vazio".
3. `explorer` limita-se aos **50 vOTUs mais longos** com anotação, e o `check_budget` continua abaixo de 25 MB.

- [ ] **Step 1: Escrever os testes que falham**

Acrescentar a `tests/test_renderer_v2.py`:

```python
from report.renderer_v2 import build_sequencing, build_viral


class _P:
    def __init__(self, outdir, samples):
        self.outdir = outdir
        self.samples = samples
        self.coassembly_groups = []


class _S:
    def __init__(self, outdir, samples):
        self.params = _P(outdir, samples)


def test_sequencing_ausente_omite_a_chave_em_vez_de_lista_vazia(tmp_path):
    # Fonte ausente e "trilha desligada", nao "medi zero" -- o painel precisa
    # distinguir os dois casos, entao a chave nao pode existir vazia.
    saida = build_sequencing(str(tmp_path), ["S1"])
    assert "quast" not in saida


def test_sequencing_le_quast_no_nivel_certo(tmp_path):
    d = tmp_path / "S1" / "quast"
    d.mkdir(parents=True)
    (d / "report.tsv").write_text(
        "Assembly\tassembly\n# contigs\t1234\nN50\t5678\n", encoding='utf-8')
    saida = build_sequencing(str(tmp_path), ["S1"])
    assert saida["quast"]["S1"]["# contigs"] == "1234"


def test_viral_projeta_somente_campos_declarados(tmp_path):
    d = tmp_path / "S1" / "viral" / "checkv"
    d.mkdir(parents=True)
    (d / "quality_summary.tsv").write_text(
        "contig_id\tcheckv_quality\tcompleteness\tcampo_gigante\n"
        "k141_1\tHigh-quality\t95.0\t" + "x" * 5000 + "\n", encoding='utf-8')
    saida = build_viral(str(tmp_path), ["S1"])
    tiers = saida["checkv_tiers"]["S1"]
    assert tiers == {"High-quality": 1}
    assert "campo_gigante" not in json.dumps(saida)


def test_explorer_limita_a_50_votus(tmp_path):
    feats = [{"votu_id": f"v{i}", "length": i, "features": []} for i in range(200)]
    assert len(_limita_explorer(feats)) == 50
    assert _limita_explorer(feats)[0]["length"] == 199   # os mais longos primeiro
```

O import de `_limita_explorer` entra junto com os outros no topo do arquivo de teste.

- [ ] **Step 2: Rodar e confirmar que falham**

Run: `python -m pytest tests/test_renderer_v2.py -v` — FAIL: `cannot import name 'build_sequencing'`

- [ ] **Step 3: Implementar `build_sequencing`, `build_viral` e `_limita_explorer`**

Cada um monta seu bloco a partir dos loaders já existentes, projetando com `schema.project` e **omitindo a chave** quando a fonte não existe. `_limita_explorer(linhas)` ordena por `length` decrescente e corta em 50.

- [ ] **Step 4: Ligar os dois blocos em `build_data`** (`data["sequencing"]`, `data["viral"]`), só quando não vierem vazios.

- [ ] **Step 5: Rodar, conferir o orçamento com a rodada real e commitar**

```bash
python -m pytest tests/ -q
git add scripts/report/renderer_v2.py tests/test_renderer_v2.py
git commit -m "feat(report-v2): blocos de dados das abas sequenciamento e viral"
```

---

### Task 8: Aba Sequenciamento

**Files:**
- Create: `src/report-ui/src/panels/Sequencing.jsx`
- Modify: `src/report-ui/src/App.jsx`
- Create: `src/report-ui/test/sequencing.test.jsx`

**Conteúdo** (formas ditadas pelo guia):
- fastp antes/depois: `StackedBar` normalizada (retido/descartado) + Q30 em `BarChart`.
- QUAST: `Heatmap` amostra × métrica, `normalize="per-col"` — nunca normalizar entre métricas de unidades diferentes.
- Comprimento de contig: `DistPlot` em log (vira ridgeline acima de 8 amostras).
- Mapeamento: `BarChart` ordenado.
- Cobertura: `Heatmap` de densidade (hexbin fica para quando houver dado denso; use bolha até lá).

- [ ] **Step 1: Escrever os testes que falham**

`src/report-ui/test/sequencing.test.jsx`:

```jsx
import { render, screen, fireEvent } from '@testing-library/react';
import { App } from '../src/App.jsx';

const base = { run: { title: 'VAPOR', samples: ['S1', 'S2'] }, overview: { kpis: [] } };
const comDados = {
  ...base,
  sequencing: {
    qc: [{ sample: 'S1', reads_before: 100, reads_after: 90, q30: 0.93 },
         { sample: 'S2', reads_before: 200, reads_after: 150, q30: 0.88 }],
    quast: { S1: { '# contigs': '1200', N50: '3400' }, S2: { '# contigs': '900', N50: '2900' } },
    mapping: [{ sample: 'S1', rate: 0.82 }, { sample: 'S2', rate: 0.61 }],
    lengths: { S1: [1000, 2000, 3000], S2: [1500, 2500] },
  },
};

test('a aba Sequenciamento aparece quando ha dado', () => {
  render(<App data={comDados} />);
  expect(screen.getByRole('tab', { name: 'Sequenciamento' })).toBeTruthy();
});

test('sem o bloco sequencing a aba nao e listada', () => {
  render(<App data={base} />);
  expect(screen.queryByRole('tab', { name: 'Sequenciamento' })).toBeNull();
});

test('o heatmap do QUAST normaliza por metrica, nao entre metricas', () => {
  render(<App data={comDados} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Sequenciamento' }));
  expect(screen.getByTestId('heatmap').getAttribute('data-normalize')).toBe('per-col');
});

test('o filtro de amostra restringe o painel', () => {
  render(<App data={comDados} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Sequenciamento' }));
  fireEvent.change(screen.getByLabelText('Amostra'), { target: { value: 'S1' } });
  expect(screen.getByTestId('seq-scope').textContent).toBe('S1');
});
```

- [ ] **Step 2: Rodar e confirmar que falham**
- [ ] **Step 3: Implementar `Sequencing.jsx`** com os cinco cards descritos acima, cada um dentro de `<Chart>` e com estado vazio próprio quando a sua fonte não veio.
- [ ] **Step 4: Registrar a aba em `ABAS`** (`App.jsx`), condicionada a `data.sequencing` existir.
- [ ] **Step 5: Rodar a suíte inteira, compilar e commitar.**

---

### Task 9: Aba Catálogo viral e explorador de vOTU

**Files:**
- Create: `src/report-ui/src/panels/ViralCatalog.jsx`
- Create: `src/report-ui/src/panels/VotuExplorer.jsx`
- Modify: `src/report-ui/src/App.jsx`
- Create: `src/report-ui/test/viral.test.jsx`

**Conteúdo:**
- Composição taxonômica: `RankSelector` + `BarChart` horizontal no rank ativo (cauda em "Other") + `Sunburst` para o catálogo inteiro, com drill.
- Qualidade CheckV: `StackedBar` com a escada ordinal Complete→HQ→MQ→LQ→ND (rampa boa→ruim, **não** oito matizes).
- vOTU × amostra: `Heatmap` com `sparseAsBubble`.
- Concordância de detectores: `UpSet`.
- Estilo de vida: `StackedBar` + nota visível de que BACPHLIP só é confiável em genoma completo.
- Explorador: seletor de vOTU + `GenomeTrack` com PHROGs, AMGs e anti-defesa nas coordenadas reais.

- [ ] **Step 1: Escrever os testes que falham**

`src/report-ui/test/viral.test.jsx`:

```jsx
import { render, screen, fireEvent } from '@testing-library/react';
import { App } from '../src/App.jsx';

const dados = {
  run: { title: 'VAPOR', samples: ['S1'] },
  overview: { kpis: [] },
  viral: {
    taxonomy: [
      { sample: 'S1', Phylum: 'Uroviricota', Class: 'Caudoviricetes', Family: 'Straboviridae', count: 12 },
      { sample: 'S1', Phylum: 'Uroviricota', Class: 'Caudoviricetes', Family: 'Drexlerviridae', count: 7 },
      { sample: 'S1', Phylum: '', Class: '', Family: '', count: 40 },
    ],
    checkv_tiers: { S1: { 'Complete': 2, 'High-quality': 5, 'Low-quality': 30, '': 11 } },
    detectors: { sets: { geNomad: 40, VirSorter2: 30 },
                 combos: [{ tools: ['geNomad'], count: 15 },
                          { tools: ['geNomad', 'VirSorter2'], count: 25 }] },
    catalog: { n_votus: 800, n_pool: 2000, reduction_pct: 60 },
    explorer: [{ votu_id: 'S1|k141_1', length: 30000, features: [
      { start: 100, end: 900, strand: '+', label: 'terminase', kind: 'phrog' }] }],
  },
};

test('a composicao taxonomica respeita o rank escolhido', () => {
  render(<App data={dados} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Catálogo viral' }));
  fireEvent.click(screen.getByRole('button', { name: 'Família' }));
  expect(screen.getByText('Straboviridae')).toBeTruthy();
  fireEvent.click(screen.getByRole('button', { name: 'Filo' }));
  expect(screen.queryByText('Straboviridae')).toBeNull();
  expect(screen.getByText('Uroviricota')).toBeTruthy();
});

test('tier vazio do CheckV nao vira Complete nem some', () => {
  render(<App data={dados} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Catálogo viral' }));
  const barra = screen.getByTestId('checkv-tiers');
  expect(barra.getAttribute('data-series')).toContain('Not-determined');
  expect(barra.getAttribute('data-complete')).toBe('2');
});

test('linhagem vazia vira Unclassified, nao um taxon de nome vazio', () => {
  render(<App data={dados} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Catálogo viral' }));
  expect(screen.getByText(/Unclassified/)).toBeTruthy();
});

test('o explorador desenha a trilha do vOTU escolhido', () => {
  render(<App data={dados} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Catálogo viral' }));
  expect(document.querySelector('[data-feature="terminase"]')).toBeTruthy();
});

test('sem anotacao o explorador mostra estado vazio, nao trilha vazia', () => {
  const semAnot = { ...dados, viral: { ...dados.viral, explorer: [] } };
  render(<App data={semAnot} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Catálogo viral' }));
  expect(screen.getByText(/sem anotação/i)).toBeTruthy();
});
```

- [ ] **Step 2: Rodar e confirmar que falham**
- [ ] **Step 3: Implementar `ViralCatalog.jsx` e `VotuExplorer.jsx`** com os seis blocos descritos acima.
- [ ] **Step 4: Registrar a aba em `ABAS`**, condicionada a `data.viral`.
- [ ] **Step 5: Rodar a suíte inteira, compilar e commitar.**

---

## Verificação final do plano 2

- [ ] `python -m pytest tests/ -q` e `npx vitest run` passam inteiros
- [ ] `npm run build` sem chamada de rede no bundle
- [ ] O report gerado da rodada real (`/home/alumnos/lmelo/amazon/results`, 32 amostras) abre, e as abas Sequenciamento e Catálogo viral mostram dado real
- [ ] Modo escuro conferido nas duas abas
- [ ] Nenhum gráfico taxonômico preso a um rank; nenhum treemap
- [ ] Trilha desligada em `config.yaml` → aba com estado vazio, nunca eixo quebrado

## O que fica para os planos 3 e 4

**Plano 3:** abas Catálogo de MAGs, Defesa/AMR/Plasmídeos, Pangenoma, Diversidade e Leituras; a vista por amostra nos dois modos com a marca de procedência; o inventário de regras da matriz de status (hoje só três regras per-sample são rastreadas); a remoção do report antigo.

**Plano 4:** `envs/env_network.yaml`, `rules/report_network.smk`, o SBM do graph-tool e as abas de rede.
