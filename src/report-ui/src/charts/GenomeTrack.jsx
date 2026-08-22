// Trilha genomica em pares de base: coordenada real, largura real -- desenhar
// ordem de gene quando a coordenada existe eh anti-padrao declarado no
// REPORT_VIZ_GUIDE.md (espacamento falso, larguras uniformes). As coordenadas
// vem de graca no cabecalho do Prodigal. Fita determina o sentido da seta;
// features sobrepostas empilham em faixas em vez de se esconder.
import { scaleLinear } from 'd3';
import { Chart } from '../viz/Chart.jsx';
import { AxisBottom } from '../viz/Axis.jsx';
import { useTooltip } from '../viz/Tooltip.jsx';
import { Legend } from '../viz/Legend.jsx';
import { PAL_MUTED } from '../viz/palette.js';

const MARGEM = { top: 10, right: 16, bottom: 8, left: 16 };
const ALTURA_FAIXA = 22;
const ALTURA_REGUA = 24;
const LARGURA_PONTA = 8;

// Aloca cada feature na primeira faixa cujo fim ja passou do inicio desta --
// resolve colisao empilhando, nunca escondendo (varredura tipo interval
// scheduling, greedy por inicio crescente).
function alocarFaixas(features) {
  const ordenado = [...features]
    .map((f, i) => ({ ...f, _i: i }))
    .sort((a, b) => a.start - b.start);
  const fimPorFaixa = [];
  const comFaixa = ordenado.map((f) => {
    let faixa = fimPorFaixa.findIndex((fim) => fim <= f.start);
    if (faixa === -1) {
      faixa = fimPorFaixa.length;
      fimPorFaixa.push(f.end);
    } else {
      fimPorFaixa[faixa] = f.end;
    }
    return { ...f, lane: faixa };
  });
  // Devolve na ordem original -- a alocacao de faixa nao deve reordenar o
  // desenho para quem consome data-lane por posicao.
  return comFaixa.sort((a, b) => a._i - b._i);
}

export function GenomeTrack({ length, features = [], kinds = {} }) {
  const { show, hide, node } = useTooltip();
  const comFaixa = alocarFaixas(features);
  const numFaixas = comFaixa.length ? Math.max(...comFaixa.map((f) => f.lane + 1)) : 1;
  const alturaFaixas = numFaixas * ALTURA_FAIXA;
  const alturaTotal = MARGEM.top + ALTURA_REGUA + alturaFaixas + MARGEM.bottom;

  const corDoTipo = (kind) => kinds?.[kind]?.color || PAL_MUTED;
  const rotuloDoTipo = (kind) => kinds?.[kind]?.label || kind;
  const kindsPresentes = Array.from(new Set(features.map((f) => f.kind)));

  return (
    <div data-testid="genometrack" className="chart-wrap">
      <Chart height={alturaTotal} empty={false}>
        {({ width }) => {
          const innerW = Math.max(width - MARGEM.left - MARGEM.right, 10);
          const x = scaleLinear().domain([0, length]).range([0, innerW]);
          return (
            <g transform={`translate(${MARGEM.left},${MARGEM.top})`}>
              <AxisBottom scale={x} width={innerW} height={ALTURA_REGUA} tickFormat={(v) => `${v.toLocaleString('pt-BR')} bp`} />
              <g transform={`translate(0,${ALTURA_REGUA})`}>
                {comFaixa.map((f) => {
                  const x0 = x(f.start);
                  const x1 = x(f.end);
                  const w = Math.max(x1 - x0, 1);
                  const y = f.lane * ALTURA_FAIXA + 2;
                  const alturaCaixa = ALTURA_FAIXA - 4;
                  const ponta = Math.min(w / 3, LARGURA_PONTA);
                  // Seta como codificacao da fita, nao a cor: '+' aponta para
                  // a direita (corte na ponta direita), '-' para a esquerda.
                  const pontos = f.strand === '-'
                    ? [
                        [x0 + ponta, y],
                        [x1, y],
                        [x1, y + alturaCaixa],
                        [x0 + ponta, y + alturaCaixa],
                        [x0, y + alturaCaixa / 2],
                      ]
                    : [
                        [x0, y],
                        [x1 - ponta, y],
                        [x1, y + alturaCaixa / 2],
                        [x1 - ponta, y + alturaCaixa],
                        [x0, y + alturaCaixa],
                      ];
                  return (
                    <polygon
                      key={f.label}
                      data-feature={f.label}
                      data-w={w}
                      data-strand={f.strand}
                      data-lane={f.lane}
                      points={pontos.map((p) => p.join(',')).join(' ')}
                      fill={corDoTipo(f.kind)}
                      onMouseMove={(e) => show(e, `${f.label} (${f.kind}) ${f.start.toLocaleString('pt-BR')}-${f.end.toLocaleString('pt-BR')} bp, fita ${f.strand}`)}
                      onMouseLeave={hide}
                    />
                  );
                })}
              </g>
            </g>
          );
        }}
      </Chart>
      <Legend items={kindsPresentes.map((k) => ({ label: rotuloDoTipo(k), color: corDoTipo(k) }))} />
      {node}
    </div>
  );
}
