// Aba Sequenciamento (task 8 do plano 2 do report). Cinco cards, cada um com
// seu proprio estado vazio -- as cinco fontes (qc, quast, mapping, lengths,
// depth) sao independentes: uma trilha desligada em config.yaml derruba uma
// so, nunca todas.
import { StackedBar, BarChart } from '../charts/BarChart.jsx';
import { Heatmap } from '../charts/Heatmap.jsx';
import { DistPlot } from '../charts/DistPlot.jsx';
import { useReport, TODAS } from '../state/store.jsx';

const METRICAS_QUAST_ORDEM = [
  '# contigs', 'Largest contig', 'Total length', 'GC (%)', 'N50', 'N75', 'L50', 'L75',
];

function escopoTexto(sample) {
  return sample === TODAS ? 'todas as amostras' : sample;
}

function filtraPorAmostra(linhas, sample) {
  return sample === TODAS ? linhas : linhas.filter((l) => l.sample === sample);
}

function filtraEntradas(obj, sample) {
  return Object.entries(obj ?? {}).filter(([s]) => sample === TODAS || s === sample);
}

// O QUAST manda numero como string ("1200") -- converte, descartando o que
// nao e numerico em vez de propagar NaN para a escala.
function numerico(v) {
  if (v === undefined || v === null) return null;
  const n = Number(String(v).replace(/,/g, ''));
  return Number.isFinite(n) ? n : null;
}

// q30/rate chegam como fracao (0-1) no exemplo do brief, e como percentual
// ja calculado na rodada real (renderer_v2.py:_q30_pct). Os dois convivem:
// uma fracao <=1 vira percentual, um percentual so passa direto.
function paraPercentual(v) {
  if (v === undefined || v === null) return null;
  return v <= 1 ? v * 100 : v;
}

// `lengths[amostra]` chega em tres formas: array cru (brief), {bins,...}
// (rodada real, binada em log), ou {values,...} (amostra com <20 contigs --
// ver renderer_v2.py LIMIAR_BRUTO). DistPlot ja aceita bins/values; so
// precisa do nome do grupo anexado.
function normalizaGrupoComprimento(nome, entrada) {
  if (Array.isArray(entrada)) return { name: nome, values: entrada };
  if (entrada?.bins) return { name: nome, bins: entrada.bins, n: entrada.n, min: entrada.min, max: entrada.max };
  if (entrada?.values) return { name: nome, values: entrada.values, n: entrada.n };
  return { name: nome, values: [] };
}

// Borda real do bin `i` dado o par de extremos do grid e a escala que o
// binador python usou (a mesma formula do comentario do brief: linear e
// `x0 + i*(x1-x0)/n`; log10/log1p espacam uniforme no espaco transformado).
function bordaBin(i, lo, hi, nBins, escala) {
  const frac = i / nBins;
  if (escala === 'log10') {
    const lLo = Math.log10(lo || 1e-9);
    const lHi = Math.log10(hi || 1);
    return 10 ** (lLo + frac * (lHi - lLo));
  }
  if (escala === 'log1p') {
    const lLo = Math.log1p(lo);
    const lHi = Math.log1p(hi);
    return Math.expm1(lLo + frac * (lHi - lLo));
  }
  return lo + frac * (hi - lo);
}

// Amostras com <20 contigs chegam com pares crus em vez de grade -- binamos
// em poucas celulas aqui mesmo, so para o Heatmap ter algo categorico para
// desenhar (a mesma logica de binagem 2D do python, so que menor).
function normalizaGradeCobertura(entrada) {
  if (entrada?.bins2d && entrada?.grid) return entrada;
  const pares = entrada?.values ?? [];
  if (!pares.length) return null;
  const nBins = 6;
  const xs = pares.map((p) => p[0]);
  const ys = pares.map((p) => p[1]);
  const x0 = Math.min(...xs);
  const x1 = Math.max(...xs);
  const y0 = Math.min(...ys);
  const y1 = Math.max(...ys);
  const grade = {};
  pares.forEach(([l, d]) => {
    const ix = Math.max(0, Math.min(nBins - 1, Math.floor(((l - x0) / ((x1 - x0) || 1)) * nBins)));
    const iy = Math.max(0, Math.min(nBins - 1, Math.floor(((d - y0) / ((y1 - y0) || 1)) * nBins)));
    const chave = `${ix}|${iy}`;
    grade[chave] = (grade[chave] || 0) + 1;
  });
  const bins2d = Object.entries(grade).map(([chave, count]) => {
    const [ix, iy] = chave.split('|').map(Number);
    return { ix, iy, count };
  });
  return { bins2d, grid: { n_bins: nBins, x0, x1, y0, y1, x_scale: 'linear', y_scale: 'linear' } };
}

function gradeParaHeatmap(entrada) {
  const norm = normalizaGradeCobertura(entrada);
  if (!norm || !norm.bins2d.length) return { rows: [], cols: [], values: {} };
  const { grid, bins2d } = norm;
  const rotuloX = (ix) => {
    const a = Math.round(bordaBin(ix, grid.x0, grid.x1, grid.n_bins, grid.x_scale));
    const b = Math.round(bordaBin(ix + 1, grid.x0, grid.x1, grid.n_bins, grid.x_scale));
    return `${a.toLocaleString('pt-BR')}–${b.toLocaleString('pt-BR')} bp`;
  };
  const rotuloY = (iy) => {
    const a = bordaBin(iy, grid.y0, grid.y1, grid.n_bins, grid.y_scale);
    const b = bordaBin(iy + 1, grid.y0, grid.y1, grid.n_bins, grid.y_scale);
    return `${a.toFixed(1)}–${b.toFixed(1)}x`;
  };
  const rows = [];
  const cols = [];
  const values = {};
  bins2d.forEach(({ ix, iy, count }) => {
    const r = rotuloX(ix);
    const c = rotuloY(iy);
    if (!rows.includes(r)) rows.push(r);
    if (!cols.includes(c)) cols.push(c);
    values[r] = values[r] || {};
    values[r][c] = (values[r][c] || 0) + count;
  });
  return { rows, cols, values };
}

export function Sequencing() {
  const { data, sample } = useReport();
  const seq = data?.sequencing;

  if (!seq) {
    return <p className="empty">Sem dados para esta aba nesta rodada.</p>;
  }

  const qc = filtraPorAmostra(seq.qc ?? [], sample);
  const mapping = filtraPorAmostra(seq.mapping ?? [], sample);
  const quastEntradas = filtraEntradas(seq.quast, sample);
  const lengthsEntradas = filtraEntradas(seq.lengths, sample);
  const depthEntradas = filtraEntradas(seq.depth, sample);
  const amostraCobertura = depthEntradas[0];
  const gradeCobertura = amostraCobertura ? gradeParaHeatmap(amostraCobertura[1]) : null;

  return (
    <div className="panel">
      <p className="card__scope">escopo: <span data-testid="seq-scope">{escopoTexto(sample)}</span></p>

      <section className="card">
        <h2>Qualidade da leitura (fastp)</h2>
        {qc.length ? (
          <>
            <StackedBar
              normalize
              order={['retido', 'descartado']}
              colors={{ retido: '#0d9488', descartado: '#64748b' }}
              data={qc.map((q) => ({
                name: q.sample,
                parts: { retido: q.reads_after, descartado: Math.max(q.reads_before - q.reads_after, 0) },
              }))}
            />
            <BarChart
              valueName="Q30 (%)"
              data={qc.map((q) => ({ name: q.sample, value: paraPercentual(q.q30) ?? 0 }))}
            />
          </>
        ) : <p className="empty">Sem dados de fastp nesta rodada.</p>}
      </section>

      <section className="card">
        <h2>Montagem (QUAST)</h2>
        {quastEntradas.length ? (
          <Heatmap
            normalize="per-col"
            rows={quastEntradas.map(([s]) => s)}
            cols={METRICAS_QUAST_ORDEM.filter((m) => quastEntradas.some(([, v]) => numerico(v[m]) !== null))}
            values={Object.fromEntries(quastEntradas.map(([s, v]) => [
              s,
              Object.fromEntries(
                Object.entries(v)
                  .map(([k, val]) => [k, numerico(val)])
                  .filter(([, n]) => n !== null),
              ),
            ]))}
          />
        ) : <p className="empty">Sem relatório QUAST nesta rodada.</p>}
      </section>

      <section className="card">
        <h2>Comprimento de contigs</h2>
        {lengthsEntradas.length ? (
          <DistPlot
            log
            xName="comprimento (bp)"
            groups={lengthsEntradas.map(([s, e]) => normalizaGrupoComprimento(s, e))}
          />
        ) : <p className="empty">Sem dados de comprimento nesta rodada.</p>}
      </section>

      <section className="card">
        <h2>Mapeamento de leituras</h2>
        {mapping.length ? (
          <BarChart
            sort="desc"
            valueName="% mapeado"
            data={mapping.map((m) => ({ name: m.sample, value: paraPercentual(m.rate) ?? 0 }))}
          />
        ) : <p className="empty">Sem dados de mapeamento nesta rodada.</p>}
      </section>

      <section className="card">
        <h2>Cobertura (profundidade × comprimento)</h2>
        {gradeCobertura && gradeCobertura.rows.length ? (
          <>
            <p className="card__scope">amostra: {amostraCobertura[0]}</p>
            <Heatmap
              sparseAsBubble
              rows={gradeCobertura.rows}
              cols={gradeCobertura.cols}
              values={gradeCobertura.values}
            />
          </>
        ) : <p className="empty">Sem dados de cobertura nesta rodada.</p>}
      </section>
    </div>
  );
}
