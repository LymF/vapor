// Distribuicao de uma variavel continua. A FORMA e escolhida aqui dentro,
// nunca pelo chamador (REPORT_VIZ_GUIDE.md secao 4): poucos pontos (mediana de
// values.length < TRIGGERS.densityMinN) vira strip -- uma curva de densidade
// sobre 3 pontos afirma uma distribuicao continua que o dado nao tem, porque a
// forma viria da largura de banda, nao da medida. Muitos grupos
// (> TRIGGERS.manyGroups) vira ridgeline em vez de uma fileira de boxplots
// (boxplot e proibido neste report: esconde bimodalidade).
import { scaleLinear, max as d3max, min as d3min, deviation, mean as d3mean } from 'd3';
import { Chart } from '../viz/Chart.jsx';
import { AxisBottom } from '../viz/Axis.jsx';
import { useTooltip } from '../viz/Tooltip.jsx';
import { PAL } from '../viz/palette.js';
import { TRIGGERS, pickDistributionForm } from '../viz/triggers.js';

const MARGEM = { top: 10, right: 16, bottom: 36, left: 60 };

function mediana(xs) {
  const s = [...xs].sort((a, b) => a - b);
  const n = s.length;
  if (n === 0) return 0;
  const meio = Math.floor(n / 2);
  return n % 2 ? s[meio] : (s[meio - 1] + s[meio]) / 2;
}

// KDE gaussiano com largura de banda de Silverman: 1.06 * sigma * n^(-1/5).
// sigma=0 (todos os valores iguais) usa um fallback pequeno para nao dividir
// por zero e ainda desenhar uma curva estreita e honesta.
function kde(values, xVals) {
  const n = values.length;
  const sigma = deviation(values) || 1;
  const bw = Math.max(1.06 * sigma * Math.pow(n, -1 / 5), 1e-6);
  const norm = 1 / (n * bw * Math.sqrt(2 * Math.PI));
  return xVals.map((xv) => {
    let soma = 0;
    for (const v of values) {
      const u = (xv - v) / bw;
      soma += Math.exp(-0.5 * u * u);
    }
    return soma * norm;
  });
}

function decideForma(groups) {
  if (groups.length > TRIGGERS.manyGroups) return 'ridgeline';
  const contagens = groups.map((g) => g.values.length);
  return pickDistributionForm(mediana(contagens));
}

export function DistPlot({ groups = [], xName = 'valor', log = false, cutoffs = [] }) {
  const { show, hide, node } = useTooltip();
  const vazio = groups.length === 0 || groups.every((g) => (g.values || []).length === 0);
  const forma = decideForma(groups);

  const todosValores = groups.flatMap((g) => g.values || []);
  const dominioMin = d3min(todosValores) ?? 0;
  const dominioMax = d3max(todosValores) ?? 1;

  const altura = forma === 'ridgeline' ? Math.max(groups.length * 46, 120) : 240;

  return (
    <div data-testid="distplot" data-form={forma} className="chart-wrap">
      <Chart height={altura} empty={vazio}>
        {({ width, height }) => {
          const innerW = Math.max(width - MARGEM.left - MARGEM.right, 10);
          const innerH = Math.max(height - MARGEM.top - MARGEM.bottom, 10);
          const x = scaleLinear().domain([dominioMin, dominioMax]).nice().range([0, innerW]);

          const linhasCorte = (alturaFaixa) => cutoffs.map((c) => (
            <g key={c.label} transform={`translate(${x(c.value)},0)`}>
              <line className="distplot__cutoff" y1={0} y2={alturaFaixa} strokeDasharray="4,3" />
              <text y={-4} textAnchor="middle" className="distplot__cutoff-label">{c.label}</text>
            </g>
          ));

          if (forma === 'strip') {
            const y0 = innerH / 2;
            return (
              <g transform={`translate(${MARGEM.left},${MARGEM.top})`}>
                <AxisBottom scale={x} width={innerW} height={innerH} />
                {linhasCorte(innerH)}
                {groups.map((g, gi) => (g.values || []).map((v, i) => (
                  <circle key={`${g.name}-${i}`} cx={x(v)} cy={y0} r={5}
                          fill={PAL[gi % PAL.length]} fillOpacity={0.75}
                          onMouseMove={(e) => show(e, `${g.name}: ${v}`)}
                          onMouseLeave={hide} />
                )))}
              </g>
            );
          }

          if (forma === 'density') {
            const xVals = Array.from({ length: 100 }, (_, i) => dominioMin + (i / 99) * (dominioMax - dominioMin || 1));
            return (
              <g transform={`translate(${MARGEM.left},${MARGEM.top})`}>
                <AxisBottom scale={x} width={innerW} height={innerH} />
                {linhasCorte(innerH)}
                {groups.map((g, gi) => {
                  const densidades = kde(g.values, xVals);
                  const maxD = d3max(densidades) || 1;
                  const y = scaleLinear().domain([0, maxD]).range([innerH, 0]);
                  const linha = xVals.map((xv, i) => `${i === 0 ? 'M' : 'L'}${x(xv)},${y(densidades[i])}`).join(' ');
                  return (
                    <path key={g.name} d={linha} fill="none" stroke={PAL[gi % PAL.length]} strokeWidth={2}
                          onMouseMove={(e) => show(e, g.name)} onMouseLeave={hide} />
                  );
                })}
              </g>
            );
          }

          // ridgeline: uma faixa por grupo. O gatilho de muitos grupos decide
          // SO o empilhamento (caber na tela) -- cada faixa ainda escolhe seu
          // proprio conteudo pelo n daquele grupo (TRIGGERS.densityMinN):
          // faixas com poucos pontos desenham strip, nunca uma curva de
          // densidade que a largura de banda inventaria sozinha.
          const faixaH = innerH / groups.length;
          const xVals = Array.from({ length: 80 }, (_, i) => dominioMin + (i / 79) * (dominioMax - dominioMin || 1));
          return (
            <g transform={`translate(${MARGEM.left},${MARGEM.top})`}>
              <AxisBottom scale={x} width={innerW} height={innerH} />
              {linhasCorte(innerH)}
              {groups.map((g, gi) => {
                const formaFaixa = pickDistributionForm((g.values || []).length);
                const baseY = gi * faixaH + faixaH * 0.9;
                const cor = PAL[gi % PAL.length];

                if (formaFaixa === 'strip') {
                  return (
                    <g key={g.name} data-lane={g.name} data-lane-form="strip">
                      {(g.values || []).map((v, i) => (
                        <circle key={i} cx={x(v)} cy={baseY - (faixaH * 0.9) / 2} r={4}
                                fill={cor} fillOpacity={0.75}
                                onMouseMove={(e) => show(e, `${g.name}: ${v}`)}
                                onMouseLeave={hide} />
                      ))}
                      <text x={4} y={gi * faixaH + 12} className="distplot__ridge-label">{g.name}</text>
                    </g>
                  );
                }

                const densidades = kde(g.values, xVals);
                const maxD = d3max(densidades) || 1;
                const y = scaleLinear().domain([0, maxD]).range([faixaH * 0.9, 0]);
                const area = [
                  `M${x(xVals[0])},${baseY}`,
                  ...xVals.map((xv, i) => `L${x(xv)},${gi * faixaH + y(densidades[i])}`),
                  `L${x(xVals[xVals.length - 1])},${baseY}`,
                  'Z',
                ].join(' ');
                return (
                  <g key={g.name} data-lane={g.name} data-lane-form="density">
                    <path d={area} fill={cor} fillOpacity={0.55}
                          stroke={cor} strokeWidth={1}
                          onMouseMove={(e) => show(e, g.name)} onMouseLeave={hide} />
                    <text x={4} y={gi * faixaH + 12} className="distplot__ridge-label">{g.name}</text>
                  </g>
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
