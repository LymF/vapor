# Report React + D3 — Fundação (plano 1 de 3)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** entregar um report novo em React + D3 que abre, navega, filtra por amostra e tem a aba Visão geral completa, gerado pelo Snakemake em paralelo ao report atual.

**Architecture:** três camadas com contratos explícitos. O Python lê o disco e escreve um JSON validado por schema (nenhum campo não declarado atravessa); o esbuild compila `src/report-ui/` num bundle único versionado no git; o `renderer_v2.py` costura JSON + bundle + CSS num HTML standalone. O report atual continua intacto e funcionando durante toda a construção.

**Tech Stack:** Python 3.10+, Snakemake, React 18, D3 v7, esbuild, vitest + @testing-library/react, pytest.

**Spec:** `docs/superpowers/specs/2026-08-22-report-react-d3-design.md`

## Global Constraints

- **O report é UM arquivo HTML standalone.** Sem CDN, sem fetch, sem servidor. Tudo inlined.
- **Node NUNCA é dependência de runtime da pipeline.** `envs/env_reportui.yaml` é ambiente de desenvolvimento; nenhuma regra do Snakemake o declara. O bundle compilado é versionado no git.
- **Nenhum campo não declarado no schema chega ao browser.** Projeção explícita, sempre.
- **Orçamento de payload: 25 MB.** Acima disso o build falha.
- **`docs/REPORT_VIZ_GUIDE.md` é lei.** Forma escolhida pelo shape do dado; paleta validada com `node scripts/validate_palette.js` nos modos claro e escuro antes de qualquer mudança de cor.
- **Ferramenta que falhou é lacuna, nunca zero.** Todo componente distingue os quatro estados de `load_tool_status` — `ok` / `skipped` / `failed` / `unknown` — visualmente. `unknown` (sem `done.txt`, ou vazio) nunca é sucesso.
- **Commits sem rodapé de coautoria** (convenção deste repositório).
- **O report atual não é tocado.** `scripts/report/renderer.py` e `components/*.js` só saem do repositório no plano 2, em commit separado.

---

### Task 1: Toolchain de build

**Files:**
- Create: `envs/env_reportui.yaml`
- Create: `src/report-ui/package.json`
- Create: `src/report-ui/build.mjs`
- Create: `src/report-ui/src/index.jsx`
- Create: `src/report-ui/vitest.config.js`
- Create: `src/report-ui/test/smoke.test.jsx`
- Create: `src/report-ui/.gitignore`
- Create: `scripts/report/assets/report-ui.js` (gerado, versionado)

**Interfaces:**
- Consumes: nada.
- Produces: `npm run build` em `src/report-ui/` escreve `scripts/report/assets/report-ui.js` (IIFE) e `scripts/report/assets/report-ui.css`. O bundle monta em `#vapor-root` lendo `window.VAPOR_DATA`.

- [ ] **Step 1: Criar o ambiente conda de desenvolvimento**

`envs/env_reportui.yaml`:

```yaml
# Ambiente de DESENVOLVIMENTO do report — nenhuma regra do Snakemake o usa.
# A pipeline consome apenas o bundle compilado, versionado em
# scripts/report/assets/report-ui.js. Node nao e dependencia de runtime.
name: env_reportui
channels:
  - conda-forge
dependencies:
  - nodejs>=20
  - esbuild>=0.21
```

Rodar:

```bash
mamba env create -f envs/env_reportui.yaml -p ./.envs/reportui
conda activate ./.envs/reportui
node -v   # esperado: v20 ou superior
npm -v    # esperado: 10 ou superior
```

- [ ] **Step 2: Declarar as dependências do front**

`src/report-ui/package.json`:

```json
{
  "name": "vapor-report-ui",
  "private": true,
  "type": "module",
  "scripts": {
    "build": "node build.mjs",
    "test": "vitest run"
  },
  "dependencies": {
    "d3": "^7.9.0",
    "react": "^18.3.1",
    "react-dom": "^18.3.1"
  },
  "devDependencies": {
    "@testing-library/jest-dom": "^6.4.8",
    "@testing-library/react": "^16.0.1",
    "esbuild": "^0.23.0",
    "jsdom": "^24.1.1",
    "vitest": "^2.0.5"
  }
}
```

`src/report-ui/.gitignore`:

```
node_modules/
```

Rodar: `cd src/report-ui && npm install`

- [ ] **Step 3: Escrever o script de build**

`src/report-ui/build.mjs`:

```js
import * as esbuild from 'esbuild';
import { fileURLToPath } from 'node:url';
import path from 'node:path';

const here = path.dirname(fileURLToPath(import.meta.url));
const out = path.resolve(here, '../../scripts/report/assets/report-ui.js');

await esbuild.build({
  entryPoints: [path.join(here, 'src/index.jsx')],
  bundle: true,
  format: 'iife',
  minify: true,
  jsx: 'automatic',
  target: ['es2020'],
  define: { 'process.env.NODE_ENV': '"production"' },
  loader: { '.jsx': 'jsx' },
  outfile: out,
});

console.log(`[report-ui] bundle escrito em ${out}`);
```

- [ ] **Step 4: Escrever o teste que falha**

`src/report-ui/vitest.config.js`:

```js
export default { test: { environment: 'jsdom', globals: true } };
```

`src/report-ui/test/smoke.test.jsx`:

```jsx
import { render, screen } from '@testing-library/react';
import { App } from '../src/index.jsx';

test('monta e mostra o titulo do report', () => {
  render(<App data={{ run: { title: 'VAPOR' } }} />);
  expect(screen.getByText('VAPOR')).toBeTruthy();
});
```

- [ ] **Step 5: Rodar o teste e confirmar que falha**

Run: `cd src/report-ui && npx vitest run`
Expected: FAIL — `Failed to resolve import "../src/index.jsx"`

- [ ] **Step 6: Implementar o mínimo**

`src/report-ui/src/index.jsx`:

```jsx
import { createRoot } from 'react-dom/client';

export function App({ data }) {
  return <h1>{data?.run?.title ?? 'VAPOR'}</h1>;
}

const el = typeof document !== 'undefined' && document.getElementById('vapor-root');
if (el) createRoot(el).render(<App data={window.VAPOR_DATA ?? {}} />);
```

- [ ] **Step 7: Rodar o teste e confirmar que passa**

Run: `cd src/report-ui && npx vitest run`
Expected: PASS, 1 teste

- [ ] **Step 8: Compilar e verificar que o bundle é autossuficiente**

```bash
cd src/report-ui && npm run build
grep -c "https\?://" ../../scripts/report/assets/report-ui.js || echo "sem URL externa: ok"
```

Expected: nenhuma URL de CDN. Se `grep` achar algo, é violação da restrição global de arquivo único.

- [ ] **Step 9: Commit**

```bash
git add envs/env_reportui.yaml src/report-ui scripts/report/assets/report-ui.js
git commit -m "build(report): toolchain esbuild + React para o report novo

Ambiente conda de desenvolvimento (nodejs + esbuild); a pipeline segue sem
Node em runtime e consome apenas o bundle versionado."
```

---

### Task 2: Contrato de dados — projeção por schema

**Files:**
- Create: `scripts/report/schema.py`
- Create: `tests/test_report_schema.py`

**Interfaces:**
- Consumes: nada.
- Produces:
  - `Block(name: str, fields: tuple[str, ...], key: str | None)`
  - `project(block: Block, rows: list[dict]) -> list[dict]` — devolve só os campos declarados, na ordem declarada, campo ausente vira `''`.
  - `UndeclaredField` (exceção) — levantada por `project_strict`.
  - `project_strict(block, rows)` — igual a `project`, mas levanta `UndeclaredField` se a linha trouxer campo não declarado. Usada onde a tabela de origem é nossa (e portanto um campo novo é bug), não onde é saída crua de ferramenta.

- [ ] **Step 1: Escrever o teste que falha**

`tests/test_report_schema.py`:

```python
import pytest
from scripts.report.schema import Block, project, project_strict, UndeclaredField


def test_project_mantem_apenas_campos_declarados():
    b = Block(name="mag", fields=("Bin", "Completeness"), key="Bin")
    rows = [{"Bin": "bin1", "Completeness": "95.1", "other_related_references": "x" * 5000}]
    assert project(b, rows) == [{"Bin": "bin1", "Completeness": "95.1"}]


def test_project_preenche_campo_ausente_com_vazio():
    b = Block(name="mag", fields=("Bin", "Contamination"), key="Bin")
    assert project(b, [{"Bin": "bin1"}]) == [{"Bin": "bin1", "Contamination": ""}]


def test_project_preserva_a_ordem_declarada():
    b = Block(name="mag", fields=("Bin", "Completeness"), key="Bin")
    out = project(b, [{"Completeness": "9", "Bin": "b"}])
    assert list(out[0].keys()) == ["Bin", "Completeness"]


def test_project_strict_recusa_campo_nao_declarado():
    b = Block(name="mag", fields=("Bin",), key="Bin")
    with pytest.raises(UndeclaredField) as e:
        project_strict(b, [{"Bin": "b", "surpresa": 1}])
    assert "surpresa" in str(e.value)
```

- [ ] **Step 2: Rodar e confirmar que falha**

Run: `pytest tests/test_report_schema.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'scripts.report.schema'`

- [ ] **Step 3: Implementar**

`scripts/report/schema.py`:

```python
"""Contrato de dados do report: um grafico so consome campo declarado aqui.

Existe por causa de um bug real: um loader copiava a linha crua do
classify_wf do GTDB-Tk (`base = dict(gtdb_bins[key])`), embarcando trinta
campos para consumir sete -- entre eles other_related_references, que e uma
string enorme por genoma. A projecao explicita mata essa classe de bug.
"""
from dataclasses import dataclass


class UndeclaredField(Exception):
    """Linha traz campo que o bloco nao declara."""


@dataclass(frozen=True)
class Block:
    name: str
    fields: tuple
    key: str = None


def project(block, rows):
    return [{f: row.get(f, '') for f in block.fields} for row in rows]


def project_strict(block, rows):
    declared = set(block.fields)
    for row in rows:
        extra = set(row) - declared
        if extra:
            raise UndeclaredField(
                f"bloco '{block.name}' recebeu campo(s) nao declarado(s): "
                f"{sorted(extra)}"
            )
    return project(block, rows)
```

- [ ] **Step 4: Rodar e confirmar que passa**

Run: `pytest tests/test_report_schema.py -v`
Expected: PASS, 4 testes

- [ ] **Step 5: Commit**

```bash
git add scripts/report/schema.py tests/test_report_schema.py
git commit -m "feat(report): contrato de dados por projecao de campos declarados"
```

---

### Task 3: Orçamento de payload

**Files:**
- Modify: `scripts/report/schema.py`
- Modify: `tests/test_report_schema.py`

**Interfaces:**
- Consumes: `scripts/report/schema.py` da Task 2.
- Produces:
  - `PayloadOverBudget` (exceção)
  - `payload_report(data: dict) -> list[tuple[str, int]]` — bytes por bloco, do maior para o menor.
  - `check_budget(data: dict, limit_mb: float = 25.0) -> list[tuple[str, int]]` — levanta `PayloadOverBudget` acima do limite, com os três maiores blocos na mensagem.

- [ ] **Step 1: Escrever o teste que falha**

Acrescentar a `tests/test_report_schema.py`:

```python
from scripts.report.schema import payload_report, check_budget, PayloadOverBudget


def test_payload_report_ordena_do_maior_para_o_menor():
    data = {"pequeno": [1], "grande": ["x" * 1000]}
    nomes = [nome for nome, _ in payload_report(data)]
    assert nomes == ["grande", "pequeno"]


def test_check_budget_passa_abaixo_do_limite():
    assert check_budget({"a": [1, 2, 3]}, limit_mb=1.0)


def test_check_budget_falha_acima_do_limite_e_nomeia_o_culpado():
    data = {"culpado": ["x" * 200_000], "inocente": [1]}
    with pytest.raises(PayloadOverBudget) as e:
        check_budget(data, limit_mb=0.1)
    assert "culpado" in str(e.value)
```

- [ ] **Step 2: Rodar e confirmar que falha**

Run: `pytest tests/test_report_schema.py -v`
Expected: FAIL — `ImportError: cannot import name 'payload_report'`

- [ ] **Step 3: Implementar**

Acrescentar a `scripts/report/schema.py`:

```python
import json


class PayloadOverBudget(Exception):
    """O JSON embarcado passou do orcamento."""


def payload_report(data):
    sizes = [(name, len(json.dumps(obj, ensure_ascii=False).encode('utf-8')))
             for name, obj in data.items()]
    return sorted(sizes, key=lambda item: item[1], reverse=True)


def check_budget(data, limit_mb=25.0):
    sizes = payload_report(data)
    total = sum(size for _, size in sizes)
    if total > limit_mb * 1024 * 1024:
        piores = ", ".join(f"{n} ({s / 1024 / 1024:.1f} MB)" for n, s in sizes[:3])
        raise PayloadOverBudget(
            f"payload de {total / 1024 / 1024:.1f} MB excede o orcamento de "
            f"{limit_mb} MB. Maiores blocos: {piores}. "
            f"Projete os campos com schema.project antes de embarcar."
        )
    return sizes
```

- [ ] **Step 4: Rodar e confirmar que passa**

Run: `pytest tests/test_report_schema.py -v`
Expected: PASS, 7 testes

- [ ] **Step 5: Commit**

```bash
git add scripts/report/schema.py tests/test_report_schema.py
git commit -m "feat(report): orcamento de payload de 25 MB com o bloco culpado nomeado"
```

---

### Task 4: Primitivas de visualização

**Files:**
- Create: `src/report-ui/src/viz/palette.js`
- Create: `src/report-ui/src/viz/triggers.js`
- Create: `src/report-ui/src/viz/useResize.js`
- Create: `src/report-ui/test/viz.test.jsx`

**Interfaces:**
- Consumes: nada.
- Produces:
  - `PAL: string[]` (8 cores), `PAL_MUTED: string`
  - `foldOther(counts: Record<string, number>, max = 8): Array<[string, number]>` — dobra a cauda em `"Other"`, sempre em último lugar.
  - `TRIGGERS: { manySamples: 12, densityMinN: 20, manyGroups: 8, denseScatter: 500, networkNodes: 150, tableRows: 200 }`
  - `pickDistributionForm(n: number): 'strip' | 'density'`
  - `pickAxisOrientation(nSamples: number): 'vertical' | 'horizontal'`
  - `useResize(ref): { width: number, height: number }`

- [ ] **Step 1: Escrever o teste que falha**

`src/report-ui/test/viz.test.jsx`:

```jsx
import { PAL, PAL_MUTED, foldOther } from '../src/viz/palette.js';
import { pickDistributionForm, pickAxisOrientation, TRIGGERS } from '../src/viz/triggers.js';

test('a paleta validada tem exatamente 8 cores', () => {
  expect(PAL).toHaveLength(8);
  expect(PAL).toContain('#0d9488');
});

test('foldOther dobra a cauda e poe Other em ultimo', () => {
  const counts = { a: 10, b: 9, c: 8, d: 7, e: 6, f: 5, g: 4, h: 3, i: 2, j: 1 };
  const out = foldOther(counts, 8);
  expect(out).toHaveLength(8);
  expect(out[7][0]).toBe('Other');
  expect(out[7][1]).toBe(2 + 1 + 3);
});

test('foldOther nao inventa Other quando cabe', () => {
  expect(foldOther({ a: 1, b: 2 }, 8).map((r) => r[0])).toEqual(['b', 'a']);
});

test('poucos pontos viram strip plot, muitos viram densidade', () => {
  expect(pickDistributionForm(TRIGGERS.densityMinN - 1)).toBe('strip');
  expect(pickDistributionForm(TRIGGERS.densityMinN)).toBe('density');
});

test('muitas amostras viram eixo horizontal', () => {
  expect(pickAxisOrientation(TRIGGERS.manySamples)).toBe('vertical');
  expect(pickAxisOrientation(TRIGGERS.manySamples + 1)).toBe('horizontal');
});
```

- [ ] **Step 2: Rodar e confirmar que falha**

Run: `cd src/report-ui && npx vitest run test/viz.test.jsx`
Expected: FAIL — `Failed to resolve import "../src/viz/palette.js"`

- [ ] **Step 3: Implementar a paleta**

`src/report-ui/src/viz/palette.js`:

```js
// Paleta validada com scripts/validate_palette.js nos modos claro e escuro.
// Oito e o teto: uma nona cor reprova sob CVD (REPORT_VIZ_GUIDE.md §6).
export const PAL = [
  '#0d9488', '#d97706', '#7c3aed', '#0891b2',
  '#16a34a', '#db2777', '#9333ea', '#ef4444',
];

// Neutro para "Other"/"Unknown". Nunca e um slot de identidade categorica.
export const PAL_MUTED = '#64748b';

export function foldOther(counts, max = 8) {
  const ordenado = Object.entries(counts).sort((a, b) => b[1] - a[1]);
  if (ordenado.length <= max) return ordenado;
  const cabeca = ordenado.slice(0, max - 1);
  const cauda = ordenado.slice(max - 1).reduce((soma, [, v]) => soma + v, 0);
  return [...cabeca, ['Other', cauda]];
}
```

- [ ] **Step 4: Implementar os gatilhos**

`src/report-ui/src/viz/triggers.js`:

```js
// Os gatilhos numericos do REPORT_VIZ_GUIDE.md §4. Vivem DENTRO das formas:
// no report antigo dependiam de o autor do grafico lembrar de chamar o helper
// certo, e o padrao so ficava correto por disciplina.
export const TRIGGERS = {
  manySamples: 12,
  densityMinN: 20,
  manyGroups: 8,
  denseScatter: 500,
  networkNodes: 150,
  tableRows: 200,
};

export function pickDistributionForm(n) {
  // Uma curva de densidade sobre poucos pontos afirma uma distribuicao
  // continua que o dado nao tem: a forma vem da largura de banda.
  return n < TRIGGERS.densityMinN ? 'strip' : 'density';
}

export function pickAxisOrientation(nSamples) {
  return nSamples > TRIGGERS.manySamples ? 'horizontal' : 'vertical';
}
```

- [ ] **Step 5: Rodar e confirmar que passa**

Run: `cd src/report-ui && npx vitest run test/viz.test.jsx`
Expected: PASS, 5 testes

- [ ] **Step 6: Escrever o teste do hook de resize**

Acrescentar a `src/report-ui/test/viz.test.jsx`:

```jsx
import { render } from '@testing-library/react';
import { useRef } from 'react';
import { useResize } from '../src/viz/useResize.js';

function Sonda() {
  const ref = useRef(null);
  const { width } = useResize(ref);
  return <div ref={ref} data-testid="alvo" data-w={width} />;
}

test('useResize comeca em zero e nao quebra sem ResizeObserver real', () => {
  const { getByTestId } = render(<Sonda />);
  expect(getByTestId('alvo').getAttribute('data-w')).toBe('0');
});
```

- [ ] **Step 7: Rodar e confirmar que falha**

Run: `cd src/report-ui && npx vitest run test/viz.test.jsx`
Expected: FAIL — `Failed to resolve import "../src/viz/useResize.js"`

- [ ] **Step 8: Implementar o hook**

`src/report-ui/src/viz/useResize.js`:

```js
import { useEffect, useState } from 'react';

export function useResize(ref) {
  const [tamanho, setTamanho] = useState({ width: 0, height: 0 });

  useEffect(() => {
    const alvo = ref.current;
    if (!alvo || typeof ResizeObserver === 'undefined') return undefined;
    const obs = new ResizeObserver(([entrada]) => {
      const { width, height } = entrada.contentRect;
      setTamanho({ width, height });
    });
    obs.observe(alvo);
    return () => obs.disconnect();
  }, [ref]);

  return tamanho;
}
```

- [ ] **Step 9: Rodar e confirmar que passa**

Run: `cd src/report-ui && npx vitest run test/viz.test.jsx`
Expected: PASS, 6 testes

- [ ] **Step 10: Commit**

```bash
git add src/report-ui/src/viz src/report-ui/test/viz.test.jsx
git commit -m "feat(report): primitivas de viz com os gatilhos do guia dentro das formas"
```

---

### Task 5: Formas da Visão geral — stat tile e matriz de status

**Files:**
- Create: `src/report-ui/src/charts/StatTile.jsx`
- Create: `src/report-ui/src/charts/StatusMatrix.jsx`
- Create: `src/report-ui/test/charts.test.jsx`

**Interfaces:**
- Consumes: nada de `viz/palette.js` — as cores dos quatro estados vêm de CSS custom properties definidas em `styles.css` (Task 7). Não importe `PAL` aqui: seria import morto.
- Produces:
  - `<StatTile label value sub />`
  - `<StatusMatrix rows={[{ rule, sample, status, reason }]} />` — `status` ∈ `ok` | `skipped` | `failed` | `unknown`; cada estado recebe classe CSS própria **e** um rótulo textual. Os quatro estados são os de `load_tool_status` (`scripts/report/data_loaders.py:96`), e `unknown` — `done.txt` ausente ou vazio — **não pode ser apresentado como sucesso**: foi o que fez uma rodada de AMRFinderPlus com disco cheio passar por zero biológico.

- [ ] **Step 1: Escrever o teste que falha**

`src/report-ui/test/charts.test.jsx`:

```jsx
import { render, screen } from '@testing-library/react';
import { StatTile } from '../src/charts/StatTile.jsx';
import { StatusMatrix } from '../src/charts/StatusMatrix.jsx';

test('StatTile mostra rotulo e valor', () => {
  render(<StatTile label="vOTUs" value={1234} sub="catalogo global" />);
  expect(screen.getByText('vOTUs')).toBeTruthy();
  expect(screen.getByText('1.234')).toBeTruthy();
});

test('StatusMatrix distingue falha de zero biologico', () => {
  render(<StatusMatrix rows={[
    { rule: 'defensefinder', sample: 'S1', status: 'failed', reason: 'disco cheio' },
    { rule: 'defensefinder', sample: 'S2', status: 'ok', reason: '' },
  ]} />);
  const falhou = screen.getByTestId('cell-defensefinder-S1');
  const passou = screen.getByTestId('cell-defensefinder-S2');
  expect(falhou.getAttribute('data-status')).toBe('failed');
  expect(passou.getAttribute('data-status')).toBe('ok');
  expect(falhou.getAttribute('aria-label')).toContain('disco cheio');
});

test('StatusMatrix nunca representa failed e skipped com o mesmo estado', () => {
  render(<StatusMatrix rows={[
    { rule: 'gunc', sample: 'S1', status: 'failed', reason: 'x' },
    { rule: 'gunc', sample: 'S2', status: 'skipped', reason: 'desligado' },
  ]} />);
  expect(screen.getByTestId('cell-gunc-S1').getAttribute('data-status'))
    .not.toBe(screen.getByTestId('cell-gunc-S2').getAttribute('data-status'));
});

test('unknown nao e apresentado como ok', () => {
  render(<StatusMatrix rows={[
    { rule: 'bakta', sample: 'S1', status: 'unknown', reason: 'no status recorded' },
    { rule: 'bakta', sample: 'S2', status: 'ok', reason: '' },
  ]} />);
  const desconhecido = screen.getByTestId('cell-bakta-S1');
  expect(desconhecido.getAttribute('data-status')).toBe('unknown');
  expect(desconhecido.getAttribute('aria-label')).toContain('sem status');
});
```

- [ ] **Step 2: Rodar e confirmar que falha**

Run: `cd src/report-ui && npx vitest run test/charts.test.jsx`
Expected: FAIL — `Failed to resolve import "../src/charts/StatTile.jsx"`

- [ ] **Step 3: Implementar o StatTile**

`src/report-ui/src/charts/StatTile.jsx`:

```jsx
const fmt = new Intl.NumberFormat('pt-BR');

export function StatTile({ label, value, sub }) {
  const texto = typeof value === 'number' ? fmt.format(value) : value;
  return (
    <div className="stat-tile">
      <span className="stat-tile__label">{label}</span>
      <span className="stat-tile__value">{texto}</span>
      {sub ? <span className="stat-tile__sub">{sub}</span> : null}
    </div>
  );
}
```

- [ ] **Step 4: Implementar a StatusMatrix**

`src/report-ui/src/charts/StatusMatrix.jsx`:

```jsx
// Os QUATRO estados de load_tool_status, visualmente distintos. Ferramenta que
// falhou tem de ler como LACUNA, jamais como contagem zero: um done.txt vazio
// ja fez uma rodada de AMRFinderPlus com disco cheio passar por zero
// biologico. 'unknown' (done.txt ausente ou vazio) entra no mesmo grupo: nao e
// sucesso.
const ROTULO = {
  ok: 'ok',
  skipped: 'pulado',
  failed: 'falhou',
  unknown: 'sem status registrado',
};
const GLIFO = { ok: '●', skipped: '–', failed: '✕', unknown: '?' };

export function StatusMatrix({ rows }) {
  const regras = [...new Set(rows.map((r) => r.rule))];
  const amostras = [...new Set(rows.map((r) => r.sample))];
  const porChave = new Map(rows.map((r) => [`${r.rule}||${r.sample}`, r]));

  return (
    <table className="status-matrix">
      <thead>
        <tr>
          <th scope="col">regra</th>
          {amostras.map((s) => <th scope="col" key={s}>{s}</th>)}
        </tr>
      </thead>
      <tbody>
        {regras.map((regra) => (
          <tr key={regra}>
            <th scope="row">{regra}</th>
            {amostras.map((s) => {
              const cel = porChave.get(`${regra}||${s}`);
              const status = cel?.status ?? 'unknown';
              const motivo = cel?.reason ? ` — ${cel.reason}` : '';
              return (
                <td
                  key={s}
                  data-testid={`cell-${regra}-${s}`}
                  data-status={status}
                  className={`status-matrix__cell status-matrix__cell--${status}`}
                  aria-label={`${regra} em ${s}: ${ROTULO[status]}${motivo}`}
                  title={`${ROTULO[status]}${motivo}`}
                >
                  <span className="status-matrix__glyph" aria-hidden="true">
                    {GLIFO[status] ?? '?'}
                  </span>
                </td>
              );
            })}
          </tr>
        ))}
      </tbody>
    </table>
  );
}
```

- [ ] **Step 5: Rodar e confirmar que passa**

Run: `cd src/report-ui && npx vitest run test/charts.test.jsx`
Expected: PASS, 4 testes

- [ ] **Step 6: Commit**

```bash
git add src/report-ui/src/charts src/report-ui/test/charts.test.jsx
git commit -m "feat(report): stat tile e matriz de status com falha distinta de zero"
```

---

### Task 6: Funil de atrição

**Files:**
- Create: `src/report-ui/src/charts/AttritionFunnel.jsx`
- Create: `src/report-ui/test/funnel.test.jsx`

**Interfaces:**
- Consumes: `PAL`, `PAL_MUTED` de `src/viz/palette.js`; `useResize` de `src/viz/useResize.js`.
- Produces: `<AttritionFunnel stages={[{ name, value }]} losses={{ [stageName]: [{ reason, count }] }} onSelectLoss={(stageName) => void} />` — desenha em SVG com D3 `scaleLinear`; cada perda é clicável e chama `onSelectLoss`.

- [ ] **Step 1: Escrever o teste que falha**

`src/report-ui/test/funnel.test.jsx`:

```jsx
import { render, screen, fireEvent } from '@testing-library/react';
import { AttritionFunnel } from '../src/charts/AttritionFunnel.jsx';

const stages = [
  { name: 'contigs', value: 1000 },
  { name: 'candidatos virais', value: 400 },
  { name: 'vOTUs retidos', value: 120 },
];
const losses = {
  'candidatos virais': [{ reason: 'sem sinal viral', count: 600 }],
  'vOTUs retidos': [
    { reason: 'curto e sem bin', count: 200 },
    { reason: 'CheckV baixo', count: 80 },
  ],
};

test('desenha uma barra por etapa', () => {
  render(<AttritionFunnel stages={stages} losses={losses} />);
  stages.forEach((s) => expect(screen.getByTestId(`stage-${s.name}`)).toBeTruthy());
});

test('a perda de cada etapa e a diferenca para a etapa anterior', () => {
  render(<AttritionFunnel stages={stages} losses={losses} />);
  expect(screen.getByTestId('loss-vOTUs retidos').getAttribute('data-loss')).toBe('280');
});

test('clicar na perda emite a etapa', () => {
  const visto = [];
  render(<AttritionFunnel stages={stages} losses={losses} onSelectLoss={(n) => visto.push(n)} />);
  fireEvent.click(screen.getByTestId('loss-vOTUs retidos'));
  expect(visto).toEqual(['vOTUs retidos']);
});

test('a primeira etapa nao tem perda', () => {
  render(<AttritionFunnel stages={stages} losses={losses} />);
  expect(screen.queryByTestId('loss-contigs')).toBeNull();
});
```

- [ ] **Step 2: Rodar e confirmar que falha**

Run: `cd src/report-ui && npx vitest run test/funnel.test.jsx`
Expected: FAIL — `Failed to resolve import "../src/charts/AttritionFunnel.jsx"`

- [ ] **Step 3: Implementar**

`src/report-ui/src/charts/AttritionFunnel.jsx`:

```jsx
import { useRef } from 'react';
// scaleLinear vem do pacote guarda-chuva 'd3', a unica dependencia declarada;
// o esbuild faz tree-shaking de ESM, entao o bundle leva so o que se usa.
import { scaleLinear } from 'd3';
import { PAL, PAL_MUTED } from '../viz/palette.js';
import { useResize } from '../viz/useResize.js';

const ALTURA_LINHA = 46;
const ROTULO_W = 190;

export function AttritionFunnel({ stages, losses = {}, onSelectLoss }) {
  const ref = useRef(null);
  const { width } = useResize(ref);
  const largura = width || 720;
  const maximo = Math.max(...stages.map((s) => s.value), 1);
  const x = scaleLinear().domain([0, maximo]).range([0, Math.max(largura - ROTULO_W - 24, 80)]);

  return (
    <div ref={ref} className="funnel">
      <svg width="100%" height={stages.length * ALTURA_LINHA + 8} role="img"
           aria-label="funil de atricao da rodada">
        {stages.map((etapa, i) => {
          const anterior = i > 0 ? stages[i - 1].value : null;
          const perda = anterior === null ? null : anterior - etapa.value;
          const y = i * ALTURA_LINHA;
          const motivos = losses[etapa.name] ?? [];
          return (
            <g key={etapa.name} transform={`translate(0,${y})`}>
              <text x={0} y={ALTURA_LINHA / 2} dominantBaseline="middle"
                    className="funnel__label">{etapa.name}</text>
              <rect data-testid={`stage-${etapa.name}`} x={ROTULO_W} y={8}
                    width={x(etapa.value)} height={ALTURA_LINHA - 20}
                    fill={PAL[0]} rx={3} />
              <text x={ROTULO_W + x(etapa.value) + 8} y={ALTURA_LINHA / 2}
                    dominantBaseline="middle" className="funnel__value">
                {etapa.value.toLocaleString('pt-BR')}
              </text>
              {perda !== null && perda > 0 ? (
                <rect
                  data-testid={`loss-${etapa.name}`}
                  data-loss={String(perda)}
                  x={ROTULO_W + x(etapa.value)} y={8}
                  width={x(perda)} height={ALTURA_LINHA - 20}
                  fill={PAL_MUTED} fillOpacity={0.35} rx={3}
                  style={{ cursor: motivos.length ? 'pointer' : 'default' }}
                  onClick={() => onSelectLoss?.(etapa.name)}
                >
                  <title>
                    {`perdidos: ${perda.toLocaleString('pt-BR')}`}
                    {motivos.length ? ` — ${motivos.map((m) => `${m.reason}: ${m.count}`).join('; ')}` : ''}
                  </title>
                </rect>
              ) : null}
            </g>
          );
        })}
      </svg>
    </div>
  );
}
```

- [ ] **Step 4: Rodar e confirmar que passa**

Run: `cd src/report-ui && npx vitest run test/funnel.test.jsx`
Expected: PASS, 4 testes

- [ ] **Step 5: Commit**

```bash
git add src/report-ui/src/charts/AttritionFunnel.jsx src/report-ui/test/funnel.test.jsx
git commit -m "feat(report): funil de atricao com a perda por etapa clicavel"
```

---

### Task 7: Shell, navegação e filtro global de amostra

**Files:**
- Create: `src/report-ui/src/state/store.jsx`
- Create: `src/report-ui/src/panels/Overview.jsx`
- Create: `src/report-ui/src/App.jsx`
- Create: `src/report-ui/src/styles.css`
- Modify: `src/report-ui/src/index.jsx`
- Create: `src/report-ui/test/app.test.jsx`
- Modify: `src/report-ui/test/smoke.test.jsx`

**Interfaces:**
- Consumes: `StatTile`, `StatusMatrix`, `AttritionFunnel`.
- Produces:
  - `<ReportProvider data>{children}</ReportProvider>` e `useReport()` → `{ data, sample, setSample, samples, tab, setTab }`. `sample === '__all__'` significa todas.
  - `<App data />` — nav + painel ativo.
  - o build passa a emitir também `scripts/report/assets/report-ui.css`, porque o bundle importa `styles.css`.

- [ ] **Step 1: Escrever o teste que falha**

`src/report-ui/test/app.test.jsx`:

```jsx
import { render, screen, fireEvent } from '@testing-library/react';
import { App } from '../src/App.jsx';

const data = {
  run: { title: 'VAPOR', samples: ['S1', 'S2'] },
  overview: {
    kpis: [{ label: 'Amostras', value: 2 }],
    status: [{ rule: 'megahit', sample: 'S1', status: 'ok', reason: '' }],
    funnel: {
      __all__: { stages: [{ name: 'contigs', value: 10 }], losses: {} },
      S1: { stages: [{ name: 'contigs', value: 4 }], losses: {} },
    },
  },
};

test('abre na Visao geral', () => {
  render(<App data={data} />);
  expect(screen.getByRole('tab', { selected: true }).textContent).toBe('Visão geral');
});

test('o filtro de amostra lista todas as amostras mais a opcao agregada', () => {
  render(<App data={data} />);
  const seletor = screen.getByLabelText('Amostra');
  expect([...seletor.options].map((o) => o.value)).toEqual(['__all__', 'S1', 'S2']);
});

test('trocar de amostra troca os dados do painel', () => {
  render(<App data={data} />);
  expect(screen.getByTestId('stage-contigs').getAttribute('width')).not.toBe('0');
  fireEvent.change(screen.getByLabelText('Amostra'), { target: { value: 'S1' } });
  expect(screen.getByTestId('funnel-scope').textContent).toBe('S1');
});

test('aba sem dados mostra estado vazio, nao eixo quebrado', () => {
  render(<App data={{ run: { title: 'VAPOR', samples: [] }, overview: {} }} />);
  expect(screen.getByText(/sem dados/i)).toBeTruthy();
});
```

- [ ] **Step 2: Rodar e confirmar que falha**

Run: `cd src/report-ui && npx vitest run test/app.test.jsx`
Expected: FAIL — `Failed to resolve import "../src/App.jsx"`

- [ ] **Step 3: Implementar o estado**

`src/report-ui/src/state/store.jsx`:

```jsx
import { createContext, useContext, useMemo, useState } from 'react';

const Ctx = createContext(null);

export const TODAS = '__all__';

export function ReportProvider({ data, children }) {
  const [sample, setSample] = useState(TODAS);
  const [tab, setTab] = useState('overview');
  const samples = data?.run?.samples ?? [];
  const valor = useMemo(
    () => ({ data, sample, setSample, samples, tab, setTab }),
    [data, sample, samples, tab],
  );
  return <Ctx.Provider value={valor}>{children}</Ctx.Provider>;
}

export function useReport() {
  const v = useContext(Ctx);
  if (!v) throw new Error('useReport fora de ReportProvider');
  return v;
}
```

- [ ] **Step 4: Implementar o painel Visão geral**

`src/report-ui/src/panels/Overview.jsx`:

```jsx
import { StatTile } from '../charts/StatTile.jsx';
import { StatusMatrix } from '../charts/StatusMatrix.jsx';
import { AttritionFunnel } from '../charts/AttritionFunnel.jsx';
import { useReport, TODAS } from '../state/store.jsx';

export function Overview() {
  const { data, sample } = useReport();
  const ov = data?.overview ?? {};
  const kpis = ov.kpis ?? [];
  const status = ov.status ?? [];
  const funil = (ov.funnel ?? {})[sample] ?? (ov.funnel ?? {})[TODAS] ?? null;

  if (!kpis.length && !status.length && !funil) {
    return <p className="empty">Sem dados para esta aba nesta rodada.</p>;
  }

  return (
    <div className="panel">
      <div className="kpi-row">
        {kpis.map((k) => <StatTile key={k.label} {...k} />)}
      </div>

      {funil ? (
        <section className="card">
          <h2>Atrição da rodada</h2>
          <p className="card__scope">
            escopo: <span data-testid="funnel-scope">{sample === TODAS ? 'todas as amostras' : sample}</span>
          </p>
          <AttritionFunnel stages={funil.stages} losses={funil.losses} />
        </section>
      ) : null}

      {status.length ? (
        <section className="card">
          <h2>Status das ferramentas</h2>
          <StatusMatrix rows={sample === TODAS ? status : status.filter((r) => r.sample === sample)} />
        </section>
      ) : null}
    </div>
  );
}
```

- [ ] **Step 5: Implementar o App**

`src/report-ui/src/App.jsx`:

```jsx
import { ReportProvider, useReport, TODAS } from './state/store.jsx';
import { Overview } from './panels/Overview.jsx';

// As abas do desenho. Os paineis dos planos 2 e 3 entram aqui; enquanto nao
// existem, a aba nao e listada -- nunca uma aba vazia sem explicacao.
const ABAS = [
  { id: 'overview', label: 'Visão geral', Painel: Overview },
];

function Shell() {
  const { data, sample, setSample, samples, tab, setTab } = useReport();
  const ativa = ABAS.find((a) => a.id === tab) ?? ABAS[0];
  const { Painel } = ativa;

  return (
    <div className="app">
      <nav className="nav">
        <span className="nav__brand">{data?.run?.title ?? 'VAPOR'}</span>
        <div role="tablist" className="nav__tabs">
          {ABAS.map((a) => (
            <button key={a.id} role="tab" aria-selected={a.id === ativa.id}
                    className={`nav__tab${a.id === ativa.id ? ' is-active' : ''}`}
                    onClick={() => setTab(a.id)}>
              {a.label}
            </button>
          ))}
        </div>
        <label className="nav__filter">
          <span>Amostra</span>
          <select value={sample} onChange={(e) => setSample(e.target.value)}>
            <option value={TODAS}>todas</option>
            {samples.map((s) => <option key={s} value={s}>{s}</option>)}
          </select>
        </label>
      </nav>
      <main className="main"><Painel /></main>
    </div>
  );
}

export function App({ data }) {
  return (
    <ReportProvider data={data}>
      <Shell />
    </ReportProvider>
  );
}
```

- [ ] **Step 6: Escrever o CSS mínimo**

`src/report-ui/src/styles.css`:

```css
:root {
  --bg: #f3f6fa;
  --surface: #ffffff;
  --border: #e2e8f0;
  --text: #1e293b;
  --text-2: #475569;
  --ok: #16a34a;
  --skipped: #94a3b8;
  --failed: #ef4444;
  --unknown: #d97706;
  --radius: 14px;
  --font: 'Inter', system-ui, -apple-system, sans-serif;
}

* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: var(--font); color: var(--text); background: var(--bg); }

.nav { display: flex; align-items: center; gap: 1.5rem; padding: 0 1.5rem;
       height: 58px; background: var(--surface); border-bottom: 1px solid var(--border); }
.nav__brand { font-weight: 700; letter-spacing: .02em; }
.nav__tabs { display: flex; gap: .25rem; margin-right: auto; }
.nav__tab { border: 0; background: none; padding: .5rem .9rem; border-radius: 8px;
            font: inherit; color: var(--text-2); cursor: pointer; }
.nav__tab.is-active { background: var(--bg); color: var(--text); font-weight: 600; }
.nav__filter { display: flex; align-items: center; gap: .5rem; font-size: .85rem; }

.main { padding: 1.5rem; }
.kpi-row { display: flex; flex-wrap: wrap; gap: 1rem; margin-bottom: 1.5rem; }
.stat-tile { background: var(--surface); border: 1px solid var(--border);
             border-radius: var(--radius); padding: 1rem 1.25rem; min-width: 150px;
             display: flex; flex-direction: column; gap: .15rem; }
.stat-tile__label { font-size: .78rem; color: var(--text-2); text-transform: uppercase; }
.stat-tile__value { font-size: 1.7rem; font-weight: 700; }
.stat-tile__sub { font-size: .78rem; color: var(--text-2); }

.card { background: var(--surface); border: 1px solid var(--border);
        border-radius: var(--radius); padding: 1.25rem; margin-bottom: 1.5rem; }
.card h2 { font-size: 1rem; margin-bottom: .25rem; }
.card__scope { font-size: .78rem; color: var(--text-2); margin-bottom: 1rem; }
.empty { color: var(--text-2); font-style: italic; }

.funnel__label { font-size: 12px; fill: var(--text-2); }
.funnel__value { font-size: 12px; fill: var(--text); }

.status-matrix { border-collapse: collapse; font-size: .8rem; }
.status-matrix th { text-align: left; font-weight: 600; padding: .3rem .5rem; color: var(--text-2); }
.status-matrix__cell { text-align: center; padding: .3rem .5rem; }
.status-matrix__cell--ok .status-matrix__glyph { color: var(--ok); }
.status-matrix__cell--skipped .status-matrix__glyph { color: var(--skipped); }
.status-matrix__cell--failed .status-matrix__glyph { color: var(--failed); font-weight: 700; }
.status-matrix__cell--unknown .status-matrix__glyph { color: var(--unknown); font-weight: 700; }
```

- [ ] **Step 7: Ligar CSS e App ao ponto de entrada**

Substituir `src/report-ui/src/index.jsx` por:

```jsx
import { createRoot } from 'react-dom/client';
import { App } from './App.jsx';
import './styles.css';

export { App };

const el = typeof document !== 'undefined' && document.getElementById('vapor-root');
if (el) createRoot(el).render(<App data={window.VAPOR_DATA ?? {}} />);
```

`build.mjs` não muda: quando o bundle importa um `.css`, o esbuild emite
`report-ui.css` ao lado do `outfile`, com o mesmo nome base. Confirme depois do
build que os dois arquivos existem.

Ajustar `src/report-ui/test/smoke.test.jsx` para importar de `../src/App.jsx`:

```jsx
import { render, screen } from '@testing-library/react';
import { App } from '../src/App.jsx';

test('monta e mostra o titulo do report', () => {
  render(<App data={{ run: { title: 'VAPOR', samples: [] } }} />);
  expect(screen.getByText('VAPOR')).toBeTruthy();
});
```

- [ ] **Step 8: Rodar toda a suíte JS**

Run: `cd src/report-ui && npx vitest run`
Expected: PASS — smoke 1, viz 6, charts 4, funnel 4, app 4 (19 testes)

- [ ] **Step 9: Compilar e commitar**

```bash
cd src/report-ui && npm run build && cd ../..
ls -la scripts/report/assets/report-ui.js scripts/report/assets/report-ui.css
git add src/report-ui scripts/report/assets/report-ui.js scripts/report/assets/report-ui.css
git commit -m "feat(report): shell React com navegacao e filtro global de amostra"
```

---

### Task 8: Geração do HTML pelo Snakemake

**Files:**
- Create: `scripts/report/renderer_v2.py`
- Create: `scripts/report/components/shell_v2.html`
- Create: `scripts/generate_report_v2.py`
- Create: `tests/test_renderer_v2.py`
- Modify: `rules/report.smk` (acrescentar `rule generate_report_v2` ao fim, antes de `rule multiqc`)

**Interfaces:**
- Consumes: `check_budget` de `scripts/report/schema.py`; o bundle da Task 7.
- Produces:
  - `build_data(snakemake) -> dict` — monta `{run, overview}`; `overview.funnel` é `{ "__all__": {...}, "<amostra>": {...} }`.
  - `_quebra_por_tier(tsv_path) -> list[dict]` — descartes por tier do CheckV (`reason`, `count`).
  - `_etapas(contagens: dict) -> list[dict]` — etapas do funil na ordem, omitindo a sem fonte.
  - `_funil_da_amostra(outdir, sample) -> dict`, `_funil_agregado(outdir, samples) -> dict`.
  - `render_html(data: dict, assets_dir: str, comp_dir: str) -> str`
  - `write_report(data, out_html, assets_dir, comp_dir) -> str` — checa o orçamento, escreve, devolve o caminho.

- [ ] **Step 1: Escrever o teste que falha**

`tests/test_renderer_v2.py`:

```python
import os
import pytest
from scripts.report.renderer_v2 import render_html, write_report
from scripts.report.schema import PayloadOverBudget

DADOS = {"run": {"title": "VAPOR", "samples": ["S1"]}, "overview": {"kpis": []}}


def _assets(tmp_path):
    a = tmp_path / "assets"; a.mkdir()
    (a / "report-ui.js").write_text("console.log('bundle');", encoding='utf-8')
    (a / "report-ui.css").write_text(".app{color:red}", encoding='utf-8')
    c = tmp_path / "components"; c.mkdir()
    (c / "shell_v2.html").write_text(
        "<html><head><style>{{CSS}}</style></head>"
        "<body><div id=\"vapor-root\"></div>{{DATA_JSON}}<script>{{APP_JS}}</script></body></html>",
        encoding='utf-8')
    return str(a), str(c)


def test_render_inlineia_bundle_css_e_dados(tmp_path):
    a, c = _assets(tmp_path)
    html = render_html(DADOS, a, c)
    assert "console.log('bundle');" in html
    assert ".app{color:red}" in html
    assert '"title": "VAPOR"' in html or '"title":"VAPOR"' in html
    assert "{{" not in html


def test_render_escapa_fechamento_de_script(tmp_path):
    a, c = _assets(tmp_path)
    html = render_html({"run": {"title": "</script><script>alert(1)"}}, a, c)
    assert "</script><script>alert(1)" not in html
    assert "<\\/script>" in html


def test_write_report_recusa_payload_acima_do_orcamento(tmp_path):
    a, c = _assets(tmp_path)
    gordo = {"run": {"title": "VAPOR"}, "lixo": ["x" * 300_000]}
    with pytest.raises(PayloadOverBudget):
        write_report(gordo, str(tmp_path / "r.html"), a, c, limit_mb=0.1)
    assert not os.path.exists(tmp_path / "r.html")


def test_write_report_escreve_o_arquivo(tmp_path):
    a, c = _assets(tmp_path)
    destino = str(tmp_path / "sub" / "r.html")
    assert write_report(DADOS, destino, a, c) == destino
    assert os.path.getsize(destino) > 0
```

- [ ] **Step 2: Rodar e confirmar que falha**

Run: `pytest tests/test_renderer_v2.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'scripts.report.renderer_v2'`

- [ ] **Step 3: Escrever o shell HTML**

`scripts/report/components/shell_v2.html`:

```html
<!DOCTYPE html>
<html lang="pt-BR">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>VAPOR — Report</title>
<style>{{CSS}}</style>
</head>
<body>
<div id="vapor-root"></div>
{{DATA_JSON}}
<script>{{APP_JS}}</script>
</body>
</html>
```

- [ ] **Step 4: Implementar o renderer**

`scripts/report/renderer_v2.py`:

```python
"""renderer_v2.py — costura JSON + bundle React + CSS num HTML standalone.

Quatro passos e nada mais: montar o JSON, ler o bundle, ler o CSS, escrever.
Toda a logica de grafico vive em src/report-ui/, compilada em
scripts/report/assets/report-ui.js.
"""
import json
import os

from .schema import check_budget

_HERE = os.path.dirname(__file__)
_ASSETS = os.path.join(_HERE, "assets")
_COMP = os.path.join(_HERE, "components")


def _read(path):
    with open(path, encoding='utf-8') as f:
        return f.read()


def _data_script(data):
    # O escape de "</" impede que uma string do dado feche o <script> que a
    # carrega -- e o mesmo cuidado do _jsstr do renderer antigo.
    payload = json.dumps(data, ensure_ascii=False).replace("</", "<\\/")
    return f"<script>window.VAPOR_DATA = {payload};</script>"


def render_html(data, assets_dir=_ASSETS, comp_dir=_COMP):
    shell = _read(os.path.join(comp_dir, "shell_v2.html"))
    return (shell
            .replace("{{CSS}}", _read(os.path.join(assets_dir, "report-ui.css")))
            .replace("{{DATA_JSON}}", _data_script(data))
            .replace("{{APP_JS}}", _read(os.path.join(assets_dir, "report-ui.js"))))


def write_report(data, out_html, assets_dir=_ASSETS, comp_dir=_COMP, limit_mb=25.0):
    check_budget(data, limit_mb=limit_mb)
    html = render_html(data, assets_dir, comp_dir)
    os.makedirs(os.path.dirname(out_html) or '.', exist_ok=True)
    with open(out_html, 'w', encoding='utf-8') as f:
        f.write(html)
    print(f"[VAPOR] Report (v2) escrito em {out_html}")
    return out_html
```

- [ ] **Step 5: Rodar e confirmar que passa**

Run: `pytest tests/test_renderer_v2.py -v`
Expected: PASS, 4 testes

- [ ] **Step 6: Escrever o montador de dados e o ponto de entrada**

Acrescentar a `scripts/report/renderer_v2.py`:

```python
from .data_loaders import load_tool_status, parse_fasta_lengths  # ampliado no Step 7


def build_data(snakemake):
    """Monta o dicionario do report a partir do que ja existe em disco.

    So o bloco 'overview' nesta fase; as demais abas entram no plano 2.
    """
    outdir = snakemake.params.outdir
    samples = list(snakemake.params.samples)
    grupos = list(getattr(snakemake.params, 'coassembly_groups', []) or [])

    # load_tool_status devolve {sample: {tool: {state, reason, raw}}}, com as
    # pseudo-amostras "(global)" e "(coassembly) <grupo>". A StatusMatrix quer
    # linhas, entao achatamos aqui -- e 'raw' NAO atravessa: e o conteudo bruto
    # do done.txt, que nenhum componente le.
    status = [
        {"rule": tool, "sample": unidade,
         "status": entrada.get("state", "unknown"),
         "reason": entrada.get("reason", "")}
        for unidade, ferramentas in load_tool_status(outdir, samples, grupos).items()
        for tool, entrada in sorted(ferramentas.items())
    ]

    contigs = {}
    for s in samples:
        caminho = os.path.join(outdir, s, "final", "viral", "viral_nonredundant.fasta")
        contigs[s] = len(parse_fasta_lengths(caminho)) if os.path.exists(caminho) else 0

    kpis = [
        {"label": "Amostras", "value": len(samples)},
        {"label": "Grupos", "value": len(grupos)},
        {"label": "vOTUs retidos", "value": sum(contigs.values())},
    ]

    funil = {TODAS: _funil_agregado(outdir, samples)}
    for s in samples:
        funil[s] = _funil_da_amostra(outdir, s)

    return {
        "run": {"title": "VAPOR", "samples": samples, "groups": grupos},
        "overview": {"kpis": kpis, "status": status, "funnel": funil},
    }
```

`scripts/generate_report_v2.py`:

```python
"""Ponto de entrada Snakemake do report v2. Toda a logica esta no pacote."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from report.renderer_v2 import build_data, write_report  # noqa: E402

write_report(build_data(snakemake), snakemake.output.html)  # noqa: F821
```

- [ ] **Step 7: Escrever o teste do funil e implementá-lo**

O funil é a peça central da aba, e a repartição do descarte tem uma sutileza que
o teste trava: **toda linha do `viral_discarded.tsv` falhou as TRÊS armas do
portão composto ao mesmo tempo** — não foi binada pelo vRhyme, não atingiu tier
MQ+ no CheckV, e ficou abaixo de `VIRAL_MIN_CONTIG`. Não existe "o motivo" de um
descarte, então a quebra é **por tier do CheckV**, e `checkv_quality` vazio
("nunca avaliado") continua distinto de tier presente-porém-baixo — essa
distinção é a razão de o sidecar existir.

Acrescentar a `tests/test_renderer_v2.py`:

```python
from scripts.report.renderer_v2 import _quebra_por_tier, _etapas


def test_quebra_por_tier_separa_nunca_avaliado_de_tier_baixo(tmp_path):
    tsv = tmp_path / "viral_discarded.tsv"
    tsv.write_text(
        "contig_id\tlength\tcheckv_quality\tcheckv_completeness\tin_vrhyme_bin\tsource_id\n"
        "k141_1\t1200\tLow-quality\t12.0\tFalse\tS1\n"
        "k141_2\t900\t\t\tFalse\tS1\n"
        "k141_3\t800\tLow-quality\t8.0\tFalse\tS1\n",
        encoding='utf-8')
    quebra = {d["reason"]: d["count"] for d in _quebra_por_tier(str(tsv))}
    assert quebra == {"Low-quality": 2, "sem avaliação CheckV": 1}


def test_quebra_por_tier_sem_arquivo_e_lista_vazia(tmp_path):
    assert _quebra_por_tier(str(tmp_path / "nao_existe.tsv")) == []


def test_etapas_omite_a_etapa_cuja_fonte_nao_existe():
    etapas = _etapas({"contigs": 10, "candidatos virais": 0, "vOTUs retidos": 3})
    assert [e["name"] for e in etapas] == ["contigs", "vOTUs retidos"]
```

Run: `pytest tests/test_renderer_v2.py -v`
Expected: FAIL — `ImportError: cannot import name '_quebra_por_tier'`

Acrescentar a `scripts/report/renderer_v2.py`:

```python
import csv

TODAS = "__all__"

# checkv_quality vazio significa "CheckV nunca avaliou este contig", que NAO e o
# mesmo que "avaliado e ruim". viral_length_gate.format_discard_row preserva
# essa diferenca de proposito; aqui ela vira um rotulo proprio.
SEM_AVALIACAO = "sem avaliação CheckV"


def _quebra_por_tier(tsv_path):
    # Descartes do portao composto (item (e)) agrupados por tier do CheckV.
    #
    # NAO e uma quebra por "motivo": toda linha desse arquivo falhou as tres
    # armas do portao ao mesmo tempo (sem bin do vRhyme, tier abaixo de MQ, e
    # comprimento abaixo de VIRAL_MIN_CONTIG). Eleger uma das armas como causa
    # seria inventar informacao que o dado nao tem.
    if not os.path.exists(tsv_path):
        return []
    contagem = {}
    with open(tsv_path, encoding='utf-8', newline='') as fh:
        for linha in csv.DictReader(fh, delimiter='\t'):
            tier = (linha.get("checkv_quality") or "").strip() or SEM_AVALIACAO
            contagem[tier] = contagem.get(tier, 0) + 1
    return [{"reason": t, "count": n}
            for t, n in sorted(contagem.items(), key=lambda kv: -kv[1])]


def _etapas(contagens):
    # Etapas na ordem do funil, omitindo aquela cuja fonte nao existe. Zero aqui
    # significa "nao consegui ler a fonte", nao "zero biologico" -- desenhar uma
    # barra zerada afirmaria o segundo. A etapa some.
    ordem = ["reads", "contigs", "candidatos virais", "vOTUs retidos"]
    return [{"name": nome, "value": contagens[nome]}
            for nome in ordem
            if contagens.get(nome)]


def _conta_fasta(caminho):
    return len(parse_fasta_lengths(caminho)) if os.path.exists(caminho) else 0


def _funil_da_amostra(outdir, sample):
    quast = parse_quast_all(os.path.join(outdir, sample, "quast", "report.tsv"))
    n_contigs = safe_int((quast or {}).get("# contigs", 0))
    descartado = os.path.join(outdir, sample, "final", "viral", "viral_discarded.tsv")
    return {
        "stages": _etapas({
            "reads": parse_total_reads(outdir, sample),
            "contigs": n_contigs,
            "candidatos virais": _conta_fasta(os.path.join(
                outdir, sample, "viral", "consensus",
                f"{sample}_viral_consensus.fasta")),
            "vOTUs retidos": _conta_fasta(os.path.join(
                outdir, sample, "final", "viral", "viral_nonredundant.fasta")),
        }),
        "losses": {"vOTUs retidos": _quebra_por_tier(descartado)},
    }


def _funil_agregado(outdir, samples):
    por_amostra = [_funil_da_amostra(outdir, s) for s in samples]
    soma_etapas, perdas = {}, {}
    for f in por_amostra:
        for etapa in f["stages"]:
            soma_etapas[etapa["name"]] = soma_etapas.get(etapa["name"], 0) + etapa["value"]
        for motivo in f["losses"].get("vOTUs retidos", []):
            perdas[motivo["reason"]] = perdas.get(motivo["reason"], 0) + motivo["count"]
    return {
        "stages": _etapas(soma_etapas),
        "losses": {"vOTUs retidos": [
            {"reason": r, "count": n}
            for r, n in sorted(perdas.items(), key=lambda kv: -kv[1])]},
    }
```

Ajustar o import de `data_loaders` no topo de `scripts/report/renderer_v2.py` para:

```python
from .data_loaders import (
    load_tool_status, parse_fasta_lengths, parse_quast_all,
    parse_total_reads, safe_int,
)
```

Run: `pytest tests/test_renderer_v2.py -v`
Expected: PASS, 7 testes

- [ ] **Step 8: Acrescentar a regra ao Snakemake**

Em `rules/report.smk`, imediatamente antes de `rule multiqc`:

```python
rule generate_report_v2:
    """Report novo (React + D3), em paralelo ao rule generate_report enquanto
    a paridade nao e atingida. Consome o bundle versionado em
    scripts/report/assets/report-ui.js -- sem Node em runtime."""
    input:
        rules.generate_report.input,
    output:
        html = f"{OUTDIR}/report_v2.html",
    params:
        samples           = list(SAMPLES.keys()),
        outdir            = OUTDIR,
        coassembly_groups = list(GROUPS.keys()) + (
            ["multisplit"] if (COBINNING_MULTISPLIT and not LONG_READS) else []),
    benchmark:
        f"{OUTDIR}/benchmarks/generate_report_v2.tsv"
    script:
        "../scripts/generate_report_v2.py"
```

- [ ] **Step 9: Verificar que o DAG continua previsível**

```bash
conda activate snakemake
snakemake -n --use-conda --cores 1 2>&1 | tail -20
```

Expected: o dry-run conclui sem erro; `generate_report_v2` aparece com **1** job.

- [ ] **Step 10: Rodar a suíte inteira e commitar**

```bash
pytest tests/ -q
cd src/report-ui && npx vitest run && cd ../..
git add scripts/report/renderer_v2.py scripts/report/components/shell_v2.html \
        scripts/generate_report_v2.py tests/test_renderer_v2.py rules/report.smk
git commit -m "feat(report): regra generate_report_v2 gerando o HTML standalone novo"
```

---

### Task 9: Documentação

**Files:**
- Modify: `CLAUDE.md` (seção "HTML Report — `scripts/report/` package")
- Modify: `docs/REPORT_VIZ_GUIDE.md` (§4 e §8)

**Interfaces:**
- Consumes: tudo acima.
- Produces: nada de código.

- [ ] **Step 1: Atualizar o CLAUDE.md**

Acrescentar ao fim da seção do report:

```markdown
### Report v2 (React + D3) — em construção

Vive em `src/report-ui/` (React 18 + D3 v7), compilado por esbuild num bundle
único versionado em `scripts/report/assets/report-ui.js`. **Node nunca é
dependência de runtime**: `envs/env_reportui.yaml` é ambiente de
desenvolvimento e nenhuma regra o declara. `scripts/report/renderer_v2.py` faz
quatro coisas — monta o JSON, lê o bundle, lê o CSS, escreve o HTML — e checa o
orçamento de 25 MB (`scripts/report/schema.py`) antes de escrever qualquer
coisa. A regra é `generate_report_v2` (`{OUTDIR}/report_v2.html`), paralela ao
`generate_report` até a paridade. Desenho:
`docs/superpowers/specs/2026-08-22-report-react-d3-design.md`.

Para editar: `mamba env create -f envs/env_reportui.yaml -p ./.envs/reportui`,
`cd src/report-ui && npm install`, `npx vitest run`, `npm run build`, **e
commitar o bundle recompilado junto com a mudança** — o bundle no git é o que a
pipeline consome.
```

- [ ] **Step 2: Corrigir a linha obsoleta do guia de viz**

Em `docs/REPORT_VIZ_GUIDE.md` §3, na linha "Effort vs discovery", substituir o
texto da coluna Notes por:

```
needs a feature space shared across samples. **Desde o catálogo global de vOTUs
(2026-08-18) esse espaço existe para toda a rodada**, não só na co-assembly: o
catálogo clusteriza uma vez sobre todas as amostras e grupos, então a curva
global é legítima. Antes disso os vOTUs per-sample eram clusterizados
independentemente e acumulá-los contava o mesmo vírus várias vezes. Sempre
promediar sobre ordens aleatórias — uma ordenação é arbitrária
```

- [ ] **Step 3: Registrar as primitivas novas no §8 do guia**

Acrescentar à tabela de helpers do §8:

```
| `TRIGGERS` (`src/report-ui/src/viz/triggers.js`) | os gatilhos do §4 no report v2, aplicados DENTRO de cada forma em vez de dependerem de o autor chamar o helper certo |
| `foldOther` (`src/report-ui/src/viz/palette.js`) | versão v2 do dobrador de cauda; "Other" sempre em último |
| `useResize` (`src/report-ui/src/viz/useResize.js`) | dimensões do container para gráficos SVG responsivos |
```

- [ ] **Step 4: Commit**

```bash
git add CLAUDE.md docs/REPORT_VIZ_GUIDE.md
git commit -m "docs(report): registra o report v2 e corrige a regra da curva de acumulacao"
```

---

## Verificação final do plano 1

- [ ] `pytest tests/ -q` passa inteiro
- [ ] `cd src/report-ui && npx vitest run` passa inteiro (19 testes)
- [ ] o funil da aba Visão geral renderiza com dados reais da rodada, e a quebra do descarte aparece por tier do CheckV
- [ ] `npm run build` produz bundle sem nenhuma URL externa
- [ ] `snakemake -n --use-conda --cores 1` conclui com `generate_report_v2` em 1 job
- [ ] `{OUTDIR}/report_v2.html` abre no navegador via `file://`, navega e o filtro de amostra funciona
- [ ] `{OUTDIR}/report.html` (o antigo) continua sendo gerado e funcionando

---

## O que fica para os planos 2 e 3

**Plano 2 — paridade e lacunas:** abas Sequenciamento, Catálogo viral, Catálogo
de MAGs, Defesa/AMR/Plasmídeos, Pangenoma, Diversidade, Leituras; a vista por
amostra com os dois modos e a marca de procedência (§5.8 do spec); os loaders
dos dados hoje invisíveis; a remoção do report antigo em commit separado.

**Plano 3 — redes:** `envs/env_network.yaml`, `rules/report_network.smk` com as
duas regras, o contrato `nodes.tsv`/`edges.tsv`, o teste de determinismo de
layout e as abas de rede com drill-down por bloco do SBM.
