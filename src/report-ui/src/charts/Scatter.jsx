// Dispersao com zonas de corte desenhadas ao fundo. As zonas nao sao
// decoracao: num scatter de completude x contaminacao elas SAO o motivo de o
// grafico existir, porque a leitura util nao e "onde caiu o ponto" e sim "de
// que lado da linha do MIMAG caiu o ponto" (Bowers et al. 2017). Sem elas o
// usuario teria de carregar os cortes de cabeca.
import { scaleLinear, extent } from 'd3';
import { Chart } from '../viz/Chart.jsx';
import { useTooltip } from '../viz/Tooltip.jsx';
import { PAL, PAL_MUTED } from '../viz/palette.js';

const MARGEM = { top: 16, right: 18, bottom: 40, left: 52 };
const RAIO = 4;

// Cortes do MIMAG. `contamination < 5` e estrito de proposito: 5% exatos NAO
// e alta qualidade na norma, e arredondar aqui reclassificaria MAG.
export const ZONAS_MIMAG = [
  { key: 'medium-quality', label: 'Qualidade média (≥50% compl., <10% cont.)',
    x0: 50, x1: 100, y0: 0, y1: 10, color: PAL[1] },
  { key: 'high-quality', label: 'Alta qualidade (≥90% compl., <5% cont.)',
    x0: 90, x1: 100, y0: 0, y1: 5, color: PAL[0] },
];

export function zonaMIMAG(completeness, contamination) {
  if (completeness === null || completeness === undefined
      || contamination === null || contamination === undefined) {
    // Sem CheckM2 nao e "baixa qualidade": e ausencia de medida, e as duas
    // precisam continuar distinguiveis (ferramenta que falhou e lacuna).
    return 'unknown';
  }
  if (completeness >= 90 && contamination < 5) return 'high-quality';
  if (completeness >= 50 && contamination < 10) return 'medium-quality';
  return 'low-quality';
}

export const COR_ZONA = {
  'high-quality': PAL[0],
  'medium-quality': PAL[1],
  'low-quality': PAL[7],
  unknown: PAL_MUTED,
};

export function Scatter({
  points = [], xName = 'x', yName = 'y', xDomain, yDomain, zones = [],
  colorOf, flagAttr, flagOf, tooltipOf, testid = 'scatter', dataAttrs = {},
}) {
  const { show, hide, node } = useTooltip();
  const validos = points.filter((p) => p.x !== null && p.x !== undefined
                                    && p.y !== null && p.y !== undefined);
  const vazio = validos.length === 0;

  return (
    <div data-testid={testid} className="chart-wrap" {...dataAttrs}>
      <Chart height={320} empty={vazio}>
        {({ width, height }) => {
          const innerW = Math.max(width - MARGEM.left - MARGEM.right, 10);
          const innerH = Math.max(height - MARGEM.top - MARGEM.bottom, 10);
          const dx = xDomain ?? extent(validos, (p) => p.x);
          const dy = yDomain ?? extent(validos, (p) => p.y);
          const x = scaleLinear().domain(dx).nice().range([0, innerW]);
          const y = scaleLinear().domain(dy).nice().range([innerH, 0]);

          return (
            <g transform={`translate(${MARGEM.left},${MARGEM.top})`}>
              {zones.map((z) => (
                <rect
                  key={z.key}
                  data-mimag={z.key}
                  x={x(Math.max(z.x0, x.domain()[0]))}
                  y={y(Math.min(z.y1, y.domain()[1]))}
                  width={Math.max(x(Math.min(z.x1, x.domain()[1]))
                                  - x(Math.max(z.x0, x.domain()[0])), 0)}
                  height={Math.max(y(Math.max(z.y0, y.domain()[0]))
                                   - y(Math.min(z.y1, y.domain()[1])), 0)}
                  fill={z.color}
                  fillOpacity={0.08}
                  stroke={z.color}
                  strokeDasharray="4 3"
                />
              ))}

              {validos.map((p) => {
                const marcado = flagOf ? flagOf(p) : null;
                const attrs = marcado && flagAttr ? { [flagAttr]: marcado } : {};
                return (
                  <circle
                    key={p.id}
                    {...attrs}
                    cx={x(p.x)}
                    cy={y(p.y)}
                    r={marcado ? RAIO + 2 : RAIO}
                    fill={colorOf ? colorOf(p) : PAL[0]}
                    fillOpacity={marcado ? 0.25 : 0.75}
                    // Glifo, nao so cor: o ponto marcado ganha anel grosso,
                    // legivel em P&B e sob daltonismo.
                    stroke={marcado ? 'var(--text-1, #111)' : 'none'}
                    strokeWidth={marcado ? 2 : 0}
                    onMouseMove={(e) => show(e, tooltipOf ? tooltipOf(p) : p.id)}
                    onMouseLeave={hide}
                  />
                );
              })}

              <text className="axis__name" x={innerW / 2} y={innerH + 32}
                    textAnchor="middle">{xName}</text>
              <text className="axis__name" transform={`translate(-38,${innerH / 2}) rotate(-90)`}
                    textAnchor="middle">{yName}</text>
              <g className="axis axis--bottom" transform={`translate(0,${innerH})`}>
                {x.ticks(5).map((t) => (
                  <text key={t} x={x(t)} y={16} textAnchor="middle">{t}</text>
                ))}
              </g>
              <g className="axis axis--left">
                {y.ticks(5).map((t) => (
                  <text key={t} x={-8} y={y(t)} dy="0.32em" textAnchor="end">{t}</text>
                ))}
              </g>
            </g>
          );
        }}
      </Chart>
      {node}
    </div>
  );
}
