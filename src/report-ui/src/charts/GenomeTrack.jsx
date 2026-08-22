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
  // Normaliza start/end aqui tambem -- a colisao e sobre a posicao real na
  // regua, nao sobre a ordem crua dos numeros (que start > end pode inverter).
  const ordenado = [...features]
    .map((f, i) => ({ ...f, _i: i, _inicio: Math.min(f.start, f.end), _fim: Math.max(f.start, f.end) }))
    .sort((a, b) => a._inicio - b._inicio);
  const fimPorFaixa = [];
  const comFaixa = ordenado.map((f) => {
    let faixa = fimPorFaixa.findIndex((fim) => fim <= f._inicio);
    if (faixa === -1) {
      faixa = fimPorFaixa.length;
      fimPorFaixa.push(f._fim);
    } else {
      fimPorFaixa[faixa] = f._fim;
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

  // Um genoma sem comprimento conhecido nao tem trilha para mostrar --
  // regua com dominio [0, NaN]/[0, 0] mentiria a escala inteira (NaN se
  // propaga para data-w e para os pontos do poligono). Estado vazio
  // explicito e preferivel a inventar uma regua.
  const comprimentoValido = Number.isFinite(length) && length > 0;

  return (
    <div data-testid="genometrack" className="chart-wrap">
      <Chart
        height={alturaTotal}
        empty={!comprimentoValido}
        emptyLabel="Comprimento do genoma desconhecido -- sem trilha para desenhar"
      >
        {({ width }) => {
          const innerW = Math.max(width - MARGEM.left - MARGEM.right, 10);
          const x = scaleLinear().domain([0, length]).range([0, innerW]);
          return (
            <g transform={`translate(${MARGEM.left},${MARGEM.top})`}>
              <AxisBottom scale={x} width={innerW} height={ALTURA_REGUA} tickFormat={(v) => `${v.toLocaleString('pt-BR')} bp`} />
              <g transform={`translate(0,${ALTURA_REGUA})`}>
                {comFaixa.map((f) => {
                  // start > end acontece de verdade (formato de origem
                  // grava fita reversa assim) -- normaliza a caixa pela
                  // posicao, mas a fita continua sendo a UNICA fonte do
                  // sentido da seta, nunca a ordem crua dos numeros.
                  const inicio = Math.min(f.start, f.end);
                  const fim = Math.max(f.start, f.end);
                  const x0 = x(inicio);
                  const x1 = x(fim);
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
                      onMouseMove={(e) => show(e, `${f.label} (${f.kind}) ${inicio.toLocaleString('pt-BR')}-${fim.toLocaleString('pt-BR')} bp, fita ${f.strand}`)}
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
