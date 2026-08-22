// Eixos recessivos: sem <line> de dominio grossa, grade fina em var(--border),
// rotulos em var(--text-3) (REPORT_VIZ_GUIDE.md secao 3 -- "recessive grid and
// axes"). scale.ticks()/scale.domain() cobrem tanto escalas continuas
// (scaleLinear/scaleTime, tem .ticks) quanto escalas de banda (scaleBand, so
// tem .domain) sem depender de d3-axis, que nao esta no package.json.

function ticksDe(scale, n) {
  if (typeof scale.ticks === 'function') return scale.ticks(n);
  if (typeof scale.domain === 'function') return scale.domain();
  return [];
}

function posDe(scale, valor) {
  const p = scale(valor);
  const largura = typeof scale.bandwidth === 'function' ? scale.bandwidth() : 0;
  return p + largura / 2;
}

export function AxisBottom({ scale, width, height, tickFormat, tickCount = 5, gridLines = true }) {
  if (!scale || !width || !height) return null;
  const ticks = ticksDe(scale, tickCount);
  const fmt = tickFormat ?? ((v) => v);

  return (
    <g className="axis axis--bottom" transform={`translate(0,${height})`}>
      {ticks.map((t, i) => {
        const x = posDe(scale, t);
        return (
          <g key={`${String(t)}-${i}`} transform={`translate(${x},0)`}>
            {gridLines ? (
              <line className="axis__grid" y1={-height} y2={0} />
            ) : null}
            <text y={16} textAnchor="middle">{fmt(t)}</text>
          </g>
        );
      })}
    </g>
  );
}

export function AxisLeft({ scale, width, height, tickFormat, tickCount = 5, gridLines = true }) {
  if (!scale || !width || !height) return null;
  const ticks = ticksDe(scale, tickCount);
  const fmt = tickFormat ?? ((v) => v);

  return (
    <g className="axis axis--left">
      {ticks.map((t, i) => {
        const y = posDe(scale, t);
        return (
          <g key={`${String(t)}-${i}`} transform={`translate(0,${y})`}>
            {gridLines ? (
              <line className="axis__grid" x1={0} x2={width} />
            ) : null}
            <text x={-8} dy="0.32em" textAnchor="end">{fmt(t)}</text>
          </g>
        );
      })}
    </g>
  );
}
