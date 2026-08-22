// Barras: sempre partem de zero (REPORT_VIZ_GUIDE.md secao 7 -- violar a
// linha de base exagera diferencas). A orientacao e escolha do proprio
// componente (secao 4 do brief), nao do chamador: acima de TRIGGERS.manySamples
// OU com algum rotulo > 12 caracteres, vira horizontal.
import { scaleLinear, scaleBand, max as d3max } from 'd3';
import { Chart } from '../viz/Chart.jsx';
import { AxisBottom, AxisLeft } from '../viz/Axis.jsx';
import { useTooltip } from '../viz/Tooltip.jsx';
import { Legend } from '../viz/Legend.jsx';
import { PAL, PAL_MUTED, foldOther } from '../viz/palette.js';
import { TRIGGERS, pickAxisOrientation } from '../viz/triggers.js';

const ROTULO_LONGO = 12;
const MARGEM = { top: 10, right: 16, bottom: 36, left: 100 };

function decideOrientacao(data, orientation) {
  if (orientation === 'vertical' || orientation === 'horizontal') return orientation;
  const porTrigger = pickAxisOrientation(data.length);
  const temRotuloLongo = data.some((d) => String(d.name).length > ROTULO_LONGO);
  return temRotuloLongo ? 'horizontal' : porTrigger;
}

export function BarChart({
  data = [], orientation = 'auto', sort = 'none', valueName = 'valor', maxSeries,
}) {
  const { show, hide, node } = useTooltip();
  const ordenado = sort === 'desc' ? [...data].sort((a, b) => b.value - a.value) : [...data];
  const limitado = maxSeries ? ordenado.slice(0, maxSeries) : ordenado;
  const orient = decideOrientacao(limitado, orientation);
  const vazio = limitado.length === 0;

  return (
    <div data-testid="barchart" data-orientation={orient} className="chart-wrap">
      <Chart height={orient === 'horizontal' ? Math.max(limitado.length * 30, 80) : 280} empty={vazio}>
        {({ width, height }) => {
          const innerW = Math.max(width - MARGEM.left - MARGEM.right, 10);
          const innerH = Math.max(height - MARGEM.top - MARGEM.bottom, 10);
          const maximo = d3max(limitado, (d) => d.value) || 1;

          if (orient === 'horizontal') {
            const y = scaleBand().domain(limitado.map((d) => d.name)).range([0, innerH]).padding(0.25);
            const x = scaleLinear().domain([0, maximo]).range([0, innerW]);
            return (
              <g transform={`translate(${MARGEM.left},${MARGEM.top})`}>
                <AxisBottom scale={x} width={innerW} height={innerH} />
                {limitado.map((d) => {
                  const len = x(d.value);
                  return (
                    <g key={d.name} transform={`translate(0,${y(d.name)})`}>
                      <text x={-8} y={y.bandwidth() / 2} dy="0.32em" textAnchor="end"
                            className="bar__label">{d.name}</text>
                      <rect data-bar data-len={len} x={0} y={0} width={len} height={y.bandwidth()}
                            fill={PAL[0]} rx={2}
                            onMouseMove={(e) => show(e, `${d.name}: ${d.value.toLocaleString('pt-BR')} ${valueName}`)}
                            onMouseLeave={hide} />
                    </g>
                  );
                })}
              </g>
            );
          }

          const x = scaleBand().domain(limitado.map((d) => d.name)).range([0, innerW]).padding(0.25);
          const y = scaleLinear().domain([0, maximo]).range([innerH, 0]);
          return (
            <g transform={`translate(${MARGEM.left},${MARGEM.top})`}>
              <AxisLeft scale={y} width={innerW} height={innerH} />
              <AxisBottom scale={x} width={innerW} height={innerH} gridLines={false} />
              {limitado.map((d) => {
                const len = innerH - y(d.value);
                return (
                  <rect key={d.name} data-bar data-len={len}
                        x={x(d.name)} y={y(d.value)} width={x.bandwidth()} height={len}
                        fill={PAL[0]} rx={2}
                        onMouseMove={(e) => show(e, `${d.name}: ${d.value.toLocaleString('pt-BR')} ${valueName}`)}
                        onMouseLeave={hide} />
                );
              })}
            </g>
          );
        }}
      </Chart>
      {node}
    </div>
  );
}

export function StackedBar({ data = [], order, normalize = false, colors }) {
  const { show, hide, node } = useTooltip();
  const totalPorCategoria = {};
  data.forEach((d) => {
    Object.entries(d.parts || {}).forEach(([k, v]) => {
      totalPorCategoria[k] = (totalPorCategoria[k] || 0) + v;
    });
  });
  const dobrado = foldOther(totalPorCategoria, 8);
  const categorias = order && order.length ? order.filter((c) => dobrado.some(([k]) => k === c)) : dobrado.map(([k]) => k);
  const cores = colors || Object.fromEntries(categorias.map((c, i) => [c, c === 'Other' ? PAL_MUTED : PAL[i % PAL.length]]));
  const vazio = data.length === 0;

  // Refolda cada linha individualmente para "Other" bater com as categorias globais.
  function linhaDobrada(parts) {
    const saida = {};
    categorias.forEach((c) => { saida[c] = 0; });
    Object.entries(parts || {}).forEach(([k, v]) => {
      saida[categorias.includes(k) ? k : 'Other'] = (saida[categorias.includes(k) ? k : 'Other'] || 0) + v;
    });
    return saida;
  }

  return (
    <div data-testid="stacked" data-series={String(categorias.length)} className="chart-wrap">
      <Chart height={Math.max(data.length * 34, 80)} empty={vazio}>
        {({ width, height }) => {
          const innerW = Math.max(width - MARGEM.left - MARGEM.right, 10);
          const innerH = Math.max(height - MARGEM.top - MARGEM.bottom, 10);
          const y = scaleBand().domain(data.map((d) => d.name)).range([0, innerH]).padding(0.25);
          const maxTotal = normalize ? 1 : d3max(data, (d) => Object.values(linhaDobrada(d.parts)).reduce((a, b) => a + b, 0)) || 1;
          const x = scaleLinear().domain([0, maxTotal]).range([0, innerW]);

          return (
            <g transform={`translate(${MARGEM.left},${MARGEM.top})`}>
              <AxisBottom scale={x} width={innerW} height={innerH} />
              {data.map((d) => {
                const linha = linhaDobrada(d.parts);
                const totalLinha = Object.values(linha).reduce((a, b) => a + b, 0) || 1;
                let acumulado = 0;
                return (
                  <g key={d.name} transform={`translate(0,${y(d.name)})`}>
                    <text x={-8} y={y.bandwidth() / 2} dy="0.32em" textAnchor="end"
                          className="bar__label">{d.name}</text>
                    {categorias.map((cat) => {
                      const valor = linha[cat] || 0;
                      const proporcao = normalize ? valor / totalLinha : valor;
                      const largura = x(proporcao) - x(0);
                      const xPos = x(acumulado);
                      acumulado += proporcao;
                      if (valor === 0) return null;
                      return (
                        <rect key={cat} data-part={cat} x={xPos} y={0} width={largura} height={y.bandwidth()}
                              fill={cores[cat]}
                              onMouseMove={(e) => show(e, `${cat}: ${valor.toLocaleString('pt-BR')}`)}
                              onMouseLeave={hide} />
                      );
                    })}
                  </g>
                );
              })}
            </g>
          );
        }}
      </Chart>
      <Legend items={categorias.map((c) => ({ label: c, color: cores[c] }))} />
      {node}
    </div>
  );
}
