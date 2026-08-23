// Matriz categorica de TRES estados. Nao e um heatmap: heatmap e escala
// sequencial de magnitude, e aqui os valores nao tem ordem de grandeza --
// presente, ausente e NAO AVALIAVEL sao categorias.
//
// A distincao que este componente existe para preservar: '?' (membro abaixo
// do piso de completude, ou com falha de anotacao) NAO e ausencia. Uma regiao
// nao montada num MAG incompleto nao e evidencia de que o organismo nao tem o
// gene, e por isso ela sai do denominador da frequencia. Pintar '?' com um
// tom mais claro da mesma rampa de '.' destruiria exatamente essa diferenca
// -- daí a hachura: textura, nao matiz, e ela sobrevive ao P&B e ao daltonismo.
import { scaleBand } from 'd3';
import { Chart } from '../viz/Chart.jsx';
import { useTooltip } from '../viz/Tooltip.jsx';
import { Legend } from '../viz/Legend.jsx';
import { PAL, PAL_MUTED } from '../viz/palette.js';

const MARGEM = { top: 70, right: 16, bottom: 8, left: 160 };
const ALTURA_LINHA = 22;

export const CORES_ESTADO = {
  x: PAL[0],
  '.': 'var(--surface-2, #e2e8f0)',
  '?': PAL_MUTED,
};

const ROTULO_ESTADO = {
  x: 'presente',
  '.': 'ausente',
  '?': 'não avaliável (fora do denominador)',
};

export function StateMatrix({ rows = [], cols = [], states = {}, rowLabel, tooltipOf }) {
  const { show, hide, node } = useTooltip();
  const vazio = rows.length === 0 || cols.length === 0;
  const rotulo = rowLabel ?? ((r) => r);
  const estadosPresentes = Array.from(new Set(
    rows.flatMap((r) => cols.map((c) => states?.[r]?.[c])).filter(Boolean),
  ));

  return (
    <div data-testid="state-matrix" data-cols={cols.join(',')} className="chart-wrap">
      <Chart height={Math.max(rows.length * ALTURA_LINHA + MARGEM.top, 140)} empty={vazio}>
        {({ width, height }) => {
          const innerW = Math.max(width - MARGEM.left - MARGEM.right, 10);
          const innerH = Math.max(height - MARGEM.top - MARGEM.bottom, 10);
          const x = scaleBand().domain(cols).range([0, innerW]).padding(0.06);
          const y = scaleBand().domain(rows).range([0, innerH]).padding(0.12);

          return (
            <g transform={`translate(${MARGEM.left},${MARGEM.top})`}>
              <defs>
                {/* Hachura: a unica marca que nao pode ser confundida com um
                    tom mais fraco de presenca/ausencia. */}
                <pattern id="hatch-nao-avaliavel" width="6" height="6"
                         patternUnits="userSpaceOnUse" patternTransform="rotate(45)">
                  <rect width="6" height="6" fill="var(--surface, #fff)" />
                  <line x1="0" y1="0" x2="0" y2="6" stroke={PAL_MUTED} strokeWidth="2.5" />
                </pattern>
              </defs>

              {cols.map((c) => (
                <text key={c} className="axis__name"
                      transform={`translate(${x(c) + x.bandwidth() / 2},-8) rotate(-40)`}
                      textAnchor="start">{c}</text>
              ))}

              {rows.map((r) => (
                <text key={`lbl-${r}`} className="matrix__row-label"
                      x={-8} y={y(r) + y.bandwidth() / 2} dy="0.32em"
                      textAnchor="end">{rotulo(r)}</text>
              ))}

              {rows.map((r) => cols.map((c) => {
                const estado = states?.[r]?.[c];
                if (!estado) return null;
                return (
                  <rect
                    key={`${r}-${c}`}
                    data-state={estado}
                    data-cell={`${r}|${c}`}
                    x={x(c)} y={y(r)}
                    width={x.bandwidth()} height={y.bandwidth()}
                    fill={estado === '?' ? 'url(#hatch-nao-avaliavel)' : CORES_ESTADO[estado]}
                    stroke="var(--border)"
                    onMouseMove={(e) => show(e, tooltipOf
                      ? tooltipOf(r, c, estado)
                      : `${rotulo(r)} × ${c}: ${ROTULO_ESTADO[estado] ?? estado}`)}
                    onMouseLeave={hide}
                  />
                );
              }))}
            </g>
          );
        }}
      </Chart>
      <Legend items={estadosPresentes.map((e) => ({
        label: ROTULO_ESTADO[e] ?? e, color: CORES_ESTADO[e],
      }))} />
      {node}
    </div>
  );
}
