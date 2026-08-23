// Procrustes: duas ordenações alinhadas, uma seta por AMOSTRA ligando o ponto
// viral ao procariótico. A unidade desenhada é o PAR, não o ponto — o que a
// figura responde é "as duas comunidades contam a mesma história nesta
// amostra?", e o comprimento da seta é a resposta.
import { scaleLinear, extent } from 'd3';
import { Chart } from '../viz/Chart.jsx';
import { useTooltip } from '../viz/Tooltip.jsx';
import { Legend } from '../viz/Legend.jsx';
import { PAL } from '../viz/palette.js';

const MARGEM = { top: 16, right: 18, bottom: 34, left: 46 };

export function ProcrustesPlot({ pairs = [], disparity }) {
  const { show, hide, node } = useTooltip();
  const vazio = pairs.length === 0;

  return (
    <div data-testid="procrustes" className="chart-wrap">
      <Chart height={340} empty={vazio}>
        {({ width, height }) => {
          const innerW = Math.max(width - MARGEM.left - MARGEM.right, 10);
          const innerH = Math.max(height - MARGEM.top - MARGEM.bottom, 10);
          const xs = pairs.flatMap((p) => [p.viral_pc1, p.prok_pc1]);
          const ys = pairs.flatMap((p) => [p.viral_pc2, p.prok_pc2]);
          const x = scaleLinear().domain(extent(xs)).nice().range([0, innerW]);
          const y = scaleLinear().domain(extent(ys)).nice().range([innerH, 0]);

          return (
            <g transform={`translate(${MARGEM.left},${MARGEM.top})`}>
              <defs>
                <marker id="proc-arrow" viewBox="0 0 10 10" refX="9" refY="5"
                        markerWidth="6" markerHeight="6" orient="auto-start-reverse">
                  <path d="M 0 0 L 10 5 L 0 10 z" fill={PAL[2]} />
                </marker>
              </defs>
              {pairs.map((p) => (
                <g key={p.sample} data-procrustes-pair={p.sample}
                   onMouseMove={(e) => show(e, `${p.sample}: viral (${p.viral_pc1}, ${p.viral_pc2}) → procariótico (${p.prok_pc1}, ${p.prok_pc2})`)}
                   onMouseLeave={hide}>
                  <line x1={x(p.viral_pc1)} y1={y(p.viral_pc2)}
                        x2={x(p.prok_pc1)} y2={y(p.prok_pc2)}
                        stroke={PAL[2]} strokeWidth={1.5}
                        markerEnd="url(#proc-arrow)" />
                  <circle cx={x(p.viral_pc1)} cy={y(p.viral_pc2)} r={4}
                          fill={PAL[0]} />
                  <circle cx={x(p.prok_pc1)} cy={y(p.prok_pc2)} r={4}
                          fill={PAL[1]} />
                </g>
              ))}
            </g>
          );
        }}
      </Chart>
      <Legend items={[
        { label: 'viral (origem da seta)', color: PAL[0] },
        { label: 'procariótico (ponta da seta)', color: PAL[1] },
      ]} />
      {disparity !== null && disparity !== undefined ? (
        <p className="chart__sub">
          Disparidade de Procrustes: {disparity} — quanto menor, mais as duas
          ordenações concordam.
        </p>
      ) : null}
      {node}
    </div>
  );
}
