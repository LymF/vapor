// Matriz linha x coluna. Sequencial de um matiz (magnitude), nunca arco-iris.
// Normaliza POR COLUNA quando as colunas carregam unidades diferentes (ex.:
// N50 em bp e numero de contigs do QUAST) -- normalizar entre colunas
// misturaria "N50 grande" com "muitos contigs" na mesma escala, o que nao faz
// sentido (REPORT_VIZ_GUIDE.md secao 3). Abaixo de 20% de celulas preenchidas
// a matriz vira quase toda vazia com pontinhos: a bolha, area proporcional ao
// valor, le melhor (secao 3, "sparse ~<20%").
import { scaleLinear, scaleBand, scaleSqrt, max as d3max } from 'd3';
import { Chart } from '../viz/Chart.jsx';
import { useTooltip } from '../viz/Tooltip.jsx';
import { PAL } from '../viz/palette.js';

const MARGEM = { top: 30, right: 16, bottom: 10, left: 110 };
const LIMIAR_ESPARSO = 0.2;
const CELULA_MIN = 4;

function densidade(rows, cols, values) {
  let preenchidas = 0;
  rows.forEach((r) => {
    cols.forEach((c) => {
      const v = values?.[r]?.[c];
      if (v !== undefined && v !== null) preenchidas += 1;
    });
  });
  const total = rows.length * cols.length || 1;
  return preenchidas / total;
}

// Normaliza cada coluna [0,1] com base no seu proprio min/max -- nunca
// compara colunas entre si.
function normalizaPorColuna(rows, cols, values) {
  const norm = {};
  cols.forEach((c) => {
    const vals = rows.map((r) => values?.[r]?.[c]).filter((v) => v !== undefined && v !== null);
    const min = vals.length ? Math.min(...vals) : 0;
    const max = vals.length ? Math.max(...vals) : 1;
    const amplitude = max - min || 1;
    rows.forEach((r) => {
      const v = values?.[r]?.[c];
      norm[r] = norm[r] || {};
      norm[r][c] = v === undefined || v === null ? null : (v - min) / amplitude;
    });
  });
  return norm;
}

export function Heatmap({ rows = [], cols = [], values = {}, normalize = 'none', sparseAsBubble = false }) {
  const { show, hide, node } = useTooltip();
  const vazio = rows.length === 0 || cols.length === 0;
  const dens = densidade(rows, cols, values);
  const modo = sparseAsBubble && dens < LIMIAR_ESPARSO ? 'bubble' : 'heat';
  const norm = normalize === 'per-col' ? normalizaPorColuna(rows, cols, values) : null;

  const todosValores = rows.flatMap((r) => cols.map((c) => values?.[r]?.[c]))
    .filter((v) => v !== undefined && v !== null);
  const maxGlobal = d3max(todosValores) || 1;

  return (
    <div data-testid="heatmap" data-mode={modo} className="chart-wrap">
      <Chart height={Math.max(rows.length * 28 + MARGEM.top, 120)} empty={vazio}>
        {({ width, height }) => {
          const innerW = Math.max(width - MARGEM.left - MARGEM.right, 10);
          const innerH = Math.max(height - MARGEM.top - MARGEM.bottom, 10);
          const x = scaleBand().domain(cols).range([0, innerW]).padding(0.08);
          const y = scaleBand().domain(rows).range([0, innerH]).padding(0.08);
          const corSeq = scaleLinear().domain([0, 1]).range(['#e6f4f1', PAL[0]]);
          const raio = scaleSqrt().domain([0, maxGlobal]).range([0, Math.min(x.bandwidth(), y.bandwidth()) / 2]);

          return (
            <g transform={`translate(${MARGEM.left},${MARGEM.top})`}>
              {cols.map((c) => (
                <text key={c} x={x(c) + x.bandwidth() / 2} y={-8} textAnchor="middle"
                      className="heatmap__axis-label">{c}</text>
              ))}
              {rows.map((r) => (
                <text key={r} x={-8} y={y(r) + y.bandwidth() / 2} dy="0.32em" textAnchor="end"
                      className="heatmap__axis-label">{r}</text>
              ))}
              {rows.map((r) => cols.map((c) => {
                const bruto = values?.[r]?.[c];
                if (bruto === undefined || bruto === null) return null;
                const normalizado = norm ? norm[r][c] : bruto / maxGlobal;
                const cx = x(c) + x.bandwidth() / 2;
                const cy = y(r) + y.bandwidth() / 2;
                const chave = `${r}|${c}`;
                const comum = {
                  'data-cell': chave,
                  'data-norm': normalizado,
                  onMouseMove: (e) => show(e, `${r} x ${c}: ${bruto}`),
                  onMouseLeave: hide,
                };
                if (modo === 'bubble') {
                  return (
                    <circle key={chave} {...comum} cx={cx} cy={cy}
                            r={Math.max(raio(bruto), CELULA_MIN)} fill={PAL[0]} fillOpacity={0.75} />
                  );
                }
                return (
                  <rect key={chave} {...comum} x={x(c)} y={y(r)} width={x.bandwidth()} height={y.bandwidth()}
                        fill={corSeq(normalizado ?? 0)} rx={2} />
                );
              }))}
            </g>
          );
        }}
      </Chart>
      {node}
    </div>
  );
}
