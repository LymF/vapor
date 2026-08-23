// Aba Pangenoma (§5.6 do desenho) — fase 1: os clusters do catálogo de MAGs
// que ganharam anotação por MEMBRO, em vez de só no representante.
//
// Todo o resto do report mostra valor computado no representante. Aqui é o
// contrário, e é o único lugar onde é: a pergunta do pangenoma é justamente
// se os membros do cluster diferem entre si, e ela não pode ser respondida
// por um genoma só.
import { useState } from 'react';
import { StateMatrix } from '../charts/StateMatrix.jsx';
import { StackedBar } from '../charts/BarChart.jsx';
import { StatTile } from '../charts/StatTile.jsx';
import { useReport } from '../state/store.jsx';
import { PAL } from '../viz/palette.js';

const CORES_CATEGORIA = { core: PAL[0], variáveis: PAL[1] };

function chaveDaLinha(r) {
  // tipo + gene: um sistema de defesa e um ARG podem ter o mesmo nome, e
  // fundir os dois daria um perfil de presença que nenhum dos dois tem.
  return `${r.tipo}|${r.gene}`;
}

export function Pangenome() {
  const { data } = useReport();
  const bloco = data?.pangenome;
  const [clusterAtivo, setClusterAtivo] = useState(null);

  if (!bloco) {
    return <p className="empty">Sem dados para esta aba nesta rodada.</p>;
  }

  const matriz = bloco.matrix ?? {};
  const repsComMatriz = Object.keys(matriz);
  const rep = clusterAtivo && matriz[clusterAtivo] ? clusterAtivo : repsComMatriz[0];
  const cluster = rep ? matriz[rep] : null;
  const resumo = (bloco.clusters ?? []).find((c) => c.representative === rep);
  const candidatos = bloco.candidates ?? [];
  const elegiveis = candidatos.filter((c) => c.eligible);

  const linhas = cluster?.rows ?? [];
  const chaves = linhas.map(chaveDaLinha);
  const estados = {};
  linhas.forEach((r) => { estados[chaveDaLinha(r)] = r.states; });
  const freqPorChave = {};
  linhas.forEach((r) => { freqPorChave[chaveDaLinha(r)] = r; });

  return (
    <div className="panel">
      <div className="kpi-row">
        <StatTile label="Clusters avaliados" value={candidatos.length} />
        <StatTile label="Clusters elegíveis" value={elegiveis.length}
                  sub="≥ 3 membros e evidência de defesa/AMR" />
        <StatTile label="Membros no cluster" value={cluster?.members?.length ?? 0} />
        <StatTile label="Membros no denominador" value={resumo?.n_evaluable ?? 0}
                  sub="≥ 70% de completude" />
      </div>

      <section className="card">
        <h2>Gene × membro</h2>
        {cluster ? (
          <>
            {repsComMatriz.length > 1 ? (
              <div className="rank-selector" role="group" aria-label="Cluster">
                {repsComMatriz.map((r) => (
                  <button key={r} type="button" className="rank-selector__btn"
                          aria-pressed={r === rep}
                          data-active={r === rep ? 'true' : 'false'}
                          onClick={() => setClusterAtivo(r)}>{r}</button>
                ))}
              </div>
            ) : null}
            <p className="notice">
              Membros abaixo de 70% de completude (ou com falha de anotação)
              aparecem hachurados e são <strong>excluídos do denominador</strong>{' '}
              da frequência: uma região não montada num MAG incompleto não é
              evidência de que o organismo não tem o gene.
            </p>
            <StateMatrix
              rows={chaves}
              cols={cluster.members}
              states={estados}
              rowLabel={(k) => {
                const r = freqPorChave[k];
                return `${r.gene} (${r.tipo}) ${r.freq}`;
              }}
              tooltipOf={(k, membro, estado) => {
                const r = freqPorChave[k];
                const compl = bloco.completeness?.[membro];
                const compTxt = compl === undefined ? 'sem CheckM2' : `${compl}% completo`;
                if (estado === '?') {
                  return `${membro}: não avaliável (${compTxt}) — fora do denominador ${r.freq}`;
                }
                return `${r.gene} (${r.tipo}) em ${membro}: `
                  + `${estado === 'x' ? 'presente' : 'ausente'} · frequência ${r.freq}`
                  + ` · ${compTxt}`;
              }}
            />
            <ul className="plain-list">
              {cluster.members.map((m) => (
                // Num nó de texto só: o ID e a completude são a mesma
                // informação (é ela que explica por que o membro é '?'), e
                // quebrá-los em dois elementos separa o número da causa.
                <li key={m}>
                  <code>{`${m} — ${bloco.completeness?.[m] ?? '—'}% de completude`}</code>
                </li>
              ))}
            </ul>
          </>
        ) : <p className="empty">Nenhum cluster com matriz gene × membro nesta rodada.</p>}
      </section>

      <section className="card">
        <h2>Core e variáveis</h2>
        {resumo ? (
          <div data-testid="pangenome-core"
               data-core={String(resumo.n_core)}
               data-variable={String(resumo.n_variable)}>
            <p className="chart__sub">
              {resumo.representative} · {resumo.n_members} membros,{' '}
              {resumo.n_evaluable} no denominador · {resumo.n_singleton} gene(s)
              em um único membro · {resumo.taxonomy || 'sem taxonomia GTDB'}
            </p>
            <StackedBar
              order={['core', 'variáveis']}
              colors={CORES_CATEGORIA}
              data={[{ name: resumo.representative,
                       parts: { core: resumo.n_core, variáveis: resumo.n_variable } }]}
            />
          </div>
        ) : <p className="empty">Sem sumário por cluster nesta rodada.</p>}
      </section>

      <section className="card">
        <h2>Candidatos</h2>
        <p className="chart__sub">
          A seleção é auditável: o cluster recusado aparece com o motivo. O
          PlasmidFinder entra como <strong>sinal de mobilidade</strong> e nunca
          elege um cluster sozinho — plasmídeo sem defesa nem ARG não motiva um
          pangenoma.
        </p>
        <div className="table-wrap">
        <table className="table">
          <thead>
            <tr>
              <th>Cluster</th><th>Membros</th><th>Ilhas</th><th>Sistemas</th>
              <th>ARGs</th><th>Replicons</th><th>Critério</th>
            </tr>
          </thead>
          <tbody>
            {candidatos.map((c) => (
              <tr key={c.representative}
                  data-testid={`candidate-${c.representative}`}
                  data-eligible={String(c.eligible)}
                  className={c.eligible ? 'is-eligible' : undefined}>
                <td><code>{c.representative}</code></td>
                <td>{c.n_members}</td>
                <td>{c.n_islands}</td>
                <td>{c.n_systems}</td>
                <td>{c.n_args}</td>
                <td>{c.n_plasmid}</td>
                <td>{c.criterio}</td>
              </tr>
            ))}
          </tbody>
        </table>
        </div>
      </section>
    </div>
  );
}
