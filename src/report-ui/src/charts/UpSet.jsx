// Concordancia entre conjuntos (ex.: ferramentas de deteccao viral). UpSet,
// nunca Venn: a partir de tres conjuntos o Venn fica ilegivel. As barras de
// interseccao vem ordenadas por tamanho, e a matriz de pertencimento mostra
// quais conjuntos compoem cada barra.
import { scaleLinear, max as d3max } from 'd3';
import { Chart } from '../viz/Chart.jsx';
import { AxisLeft } from '../viz/Axis.jsx';
import { useTooltip } from '../viz/Tooltip.jsx';

const MARGEM = { top: 10, right: 16, bottom: 10, left: 60 };
const ALTURA_BARRAS = 140;
const ALTURA_LINHA_MATRIZ = 24;
const LARGURA_LABEL_SET = 110;

export function UpSet({ sets = {}, combos = [], valueName = 'itens' }) {
  const { show, hide, node } = useTooltip();
  const nomesSets = Object.keys(sets);
  const vazio = nomesSets.length === 0 || combos.length === 0;
  const ordenados = [...combos].sort((a, b) => b.count - a.count);

  const altura = ALTURA_BARRAS + MARGEM.top + nomesSets.length * ALTURA_LINHA_MATRIZ + 20;

  return (
    <div data-testid="upset" className="chart-wrap">
      <Chart height={altura} empty={vazio}>
        {({ width }) => {
          const innerW = Math.max(width - MARGEM.left - MARGEM.right - LARGURA_LABEL_SET, 10);
          const passo = innerW / Math.max(ordenados.length, 1);
          const larguraBarra = Math.max(passo * 0.55, 6);
          const maximo = d3max(ordenados, (d) => d.count) || 1;
          const yBar = scaleLinear().domain([0, maximo]).range([ALTURA_BARRAS, 0]);
          const yEixoBarras = scaleLinear().domain([0, maximo]).range([ALTURA_BARRAS, 0]);
          const matrizTop = ALTURA_BARRAS + 24;

          return (
            <g transform={`translate(${MARGEM.left + LARGURA_LABEL_SET},${MARGEM.top})`}>
              <AxisLeft scale={yEixoBarras} width={innerW} height={ALTURA_BARRAS} tickCount={4} />
              {ordenados.map((combo, i) => {
                const cx = i * passo + passo / 2;
                const rotulo = combo.tools.join('+');
                const alturaBarra = ALTURA_BARRAS - yBar(combo.count);
                return (
                  <g key={rotulo}>
                    <rect data-combo={rotulo} data-count={combo.count}
                          x={cx - larguraBarra / 2} y={yBar(combo.count)}
                          width={larguraBarra} height={alturaBarra}
                          fill="var(--cat-1, #0d9488)"
                          onMouseMove={(e) => show(e, `${rotulo}: ${combo.count.toLocaleString('pt-BR')} ${valueName}`)}
                          onMouseLeave={hide} />
                    {nomesSets.map((s, si) => {
                      const on = combo.tools.includes(s);
                      const cy = matrizTop + si * ALTURA_LINHA_MATRIZ + ALTURA_LINHA_MATRIZ / 2;
                      return (
                        <circle key={s} cx={cx} cy={cy} r={on ? 6 : 3}
                                className={on ? 'upset__dot--on' : 'upset__dot--off'} />
                      );
                    })}
                    {(() => {
                      const onIdx = nomesSets.map((s, si) => (combo.tools.includes(s) ? si : null)).filter((v) => v !== null);
                      if (onIdx.length < 2) return null;
                      const y1 = matrizTop + onIdx[0] * ALTURA_LINHA_MATRIZ + ALTURA_LINHA_MATRIZ / 2;
                      const y2 = matrizTop + onIdx[onIdx.length - 1] * ALTURA_LINHA_MATRIZ + ALTURA_LINHA_MATRIZ / 2;
                      return <line className="upset__link" x1={cx} x2={cx} y1={y1} y2={y2} />;
                    })()}
                  </g>
                );
              })}
              {nomesSets.map((s, si) => (
                <text key={s} x={-LARGURA_LABEL_SET + 8}
                      y={matrizTop + si * ALTURA_LINHA_MATRIZ + ALTURA_LINHA_MATRIZ / 2}
                      dy="0.32em" className="upset__set-label">
                  {s} ({(sets[s] ?? 0).toLocaleString('pt-BR')})
                </text>
              ))}
            </g>
          );
        }}
      </Chart>
      {node}
    </div>
  );
}
