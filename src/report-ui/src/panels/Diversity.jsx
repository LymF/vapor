// Aba Diversidade (§5.7 do desenho).
//
// A regra que governa esta aba: índice que a rodada não calculou vira lacuna
// NOMEADA, nunca barra em zero. Simpson e Chao1 são estimadores de contagem
// (f1/f2 e a*(a-1)) e `compute_diversity.py` se recusa a calculá-los sobre
// RPKM — plotar 0 aqui afirmaria "uma única espécie domina", que é uma
// afirmação biológica forte que ninguém fez.
import { DistPlot } from '../charts/DistPlot.jsx';
import { Scatter } from '../charts/Scatter.jsx';
import { ProcrustesPlot } from '../charts/ProcrustesPlot.jsx';
import { useReport, TODAS } from '../state/store.jsx';
import { PAL } from '../viz/palette.js';

const ROTULO_INDICE = {
  observed: 'Riqueza observada',
  shannon: 'Shannon',
  simpson: 'Simpson',
  chao1: 'Chao1',
};

const ROTULO_DOMINIO = { viral: 'viral', prok: 'procariótico', combined: 'combinado' };

function pct(v) {
  return v === null || v === undefined ? '?' : (v * 100).toFixed(1);
}

export function Diversity() {
  const { data, sample } = useReport();
  const bloco = data?.diversity;
  if (!bloco) {
    return <p className="empty">Sem dados para esta aba nesta rodada.</p>;
  }

  const alpha = bloco.alpha ?? [];
  const indices = Array.from(new Set(alpha.map((r) => r.index)));
  const ausentes = bloco.alpha_missing ?? [];

  return (
    <div className="panel">
      <section className="card">
        <h2>Diversidade alfa</h2>
        {ausentes.length ? (
          <p className="notice" data-testid="alpha-missing">
            Não calculados nesta rodada: {ausentes.map((i) => ROTULO_INDICE[i] ?? i).join(', ')}.
            Ambos são estimadores de <strong>contagens</strong> de reads, e esta
            rodada não as tem — a aba mostra a lacuna em vez de plotar zero.
          </p>
        ) : null}
        {indices.map((idx) => {
          const dominios = Array.from(new Set(
            alpha.filter((r) => r.index === idx).map((r) => r.domain)));
          const grupos = dominios.map((d) => ({
            name: ROTULO_DOMINIO[d] ?? d,
            values: alpha.filter((r) => r.index === idx && r.domain === d)
              .map((r) => r.value),
          }));
          // A posição relativa é o que dá sentido a um valor de alfa isolado:
          // a amostra filtrada é marcada SOBRE a distribuição das demais, não
          // mostrada sozinha.
          const marcada = sample === TODAS ? [] : alpha
            .filter((r) => r.index === idx && r.sample === sample)
            .map((r) => ({ label: sample, value: r.value }));
          return (
            <div key={idx} data-testid={`alpha-${idx}`}>
              <h3>{ROTULO_INDICE[idx] ?? idx}</h3>
              <DistPlot groups={grupos} xName={ROTULO_INDICE[idx] ?? idx}
                        cutoffs={marcada} />
            </div>
          );
        })}
        {!alpha.length ? <p className="empty">Sem diversidade alfa nesta rodada.</p> : null}
      </section>

      <section className="card">
        <h2>Diversidade beta (PCoA, Bray-Curtis)</h2>
        {bloco.pcoa && Object.keys(bloco.pcoa).length ? (
          Object.entries(bloco.pcoa).map(([trilha, pontos]) => (
            <div key={trilha} data-testid={`pcoa-${trilha}`}>
              <h3>{ROTULO_DOMINIO[trilha] ?? trilha}</h3>
              <Scatter
                testid={`pcoa-scatter-${trilha}`}
                points={pontos.map((p) => ({ id: p.sample, x: p.pc1, y: p.pc2 }))}
                xName={`PC1 (${pct(pontos[0]?.var_pc1)}%)`}
                yName={`PC2 (${pct(pontos[0]?.var_pc2)}%)`}
                colorOf={() => PAL[0]}
                flagAttr="data-selected"
                flagOf={(p) => (sample !== TODAS && p.id === sample ? p.id : null)}
                tooltipOf={(p) => `${p.id}: PC1 ${p.x}, PC2 ${p.y}`}
              />
            </div>
          ))
        ) : <p className="empty">Sem ordenação beta nesta rodada.</p>}
      </section>

      <section className="card">
        <h2>Procrustes viral × procariótico</h2>
        {bloco.procrustes?.pairs?.length ? (
          <ProcrustesPlot pairs={bloco.procrustes.pairs}
                          disparity={bloco.procrustes.disparity} />
        ) : <p className="empty">Sem alinhamento de Procrustes nesta rodada — ele exige as duas ordenações.</p>}
      </section>
    </div>
  );
}
