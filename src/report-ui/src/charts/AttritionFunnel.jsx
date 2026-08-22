import { useRef } from 'react';
// scaleLinear vem do pacote guarda-chuva 'd3', a unica dependencia declarada;
// o esbuild faz tree-shaking de ESM, entao o bundle leva so o que se usa.
import { scaleLinear } from 'd3';
import { PAL, PAL_MUTED } from '../viz/palette.js';
import { useResize } from '../viz/useResize.js';

const ALTURA_LINHA = 46;
const ROTULO_W = 190;

export function AttritionFunnel({ stages, losses = {}, onSelectLoss }) {
  const ref = useRef(null);
  const { width } = useResize(ref);
  const largura = width || 720;
  const maximo = Math.max(...stages.map((s) => s.value), 1);
  const x = scaleLinear().domain([0, maximo]).range([0, Math.max(largura - ROTULO_W - 24, 80)]);

  return (
    <div ref={ref} className="funnel">
      <svg width="100%" height={stages.length * ALTURA_LINHA + 8} role="img"
           aria-label="funil de atricao da rodada">
        {stages.map((etapa, i) => {
          const anterior = i > 0 ? stages[i - 1].value : null;
          const perda = anterior === null ? null : anterior - etapa.value;
          const y = i * ALTURA_LINHA;
          const motivos = losses[etapa.name] ?? [];
          return (
            <g key={etapa.name} transform={`translate(0,${y})`}>
              <text x={0} y={ALTURA_LINHA / 2} dominantBaseline="middle"
                    className="funnel__label">{etapa.name}</text>
              <rect data-testid={`stage-${etapa.name}`} x={ROTULO_W} y={8}
                    width={x(etapa.value)} height={ALTURA_LINHA - 20}
                    fill={PAL[0]} rx={3} />
              <text x={ROTULO_W + x(etapa.value) + 8} y={ALTURA_LINHA / 2}
                    dominantBaseline="middle" className="funnel__value">
                {etapa.value.toLocaleString('pt-BR')}
              </text>
              {perda !== null && perda > 0 ? (
                <rect
                  data-testid={`loss-${etapa.name}`}
                  data-loss={String(perda)}
                  x={ROTULO_W + x(etapa.value)} y={8}
                  width={x(perda)} height={ALTURA_LINHA - 20}
                  fill={PAL_MUTED} fillOpacity={0.35} rx={3}
                  style={{ cursor: motivos.length ? 'pointer' : 'default' }}
                  onClick={() => onSelectLoss?.(etapa.name)}
                >
                  <title>
                    {`perdidos: ${perda.toLocaleString('pt-BR')}`}
                    {motivos.length ? ` — ${motivos.map((m) => `${m.reason}: ${m.count}`).join('; ')}` : ''}
                  </title>
                </rect>
              ) : null}
            </g>
          );
        })}
      </svg>
    </div>
  );
}
