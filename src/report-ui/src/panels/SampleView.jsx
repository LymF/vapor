// Aba Por amostra (§5.8 do desenho).
//
// As outras abas respondem "como está o catálogo?". Esta responde "como foi a
// amostra P01?", e nenhuma das outras a substitui: o filtro global repinta os
// painéis, mas não junta, num lugar só, tudo o que se sabe de uma amostra na
// ordem em que se investiga.
//
// Dois modos, porque são duas perguntas. Individual: uma amostra inteira.
// Comparação: pequenos múltiplos alinhados no MESMO eixo — nunca um sunburst
// por amostra lado a lado, porque hierarquia não se compara visualmente.
import { useState } from 'react';
import { StatTile } from '../charts/StatTile.jsx';
import { StatusMatrix } from '../charts/StatusMatrix.jsx';
import { AttritionFunnel } from '../charts/AttritionFunnel.jsx';
import { DistPlot } from '../charts/DistPlot.jsx';
import { Scatter, ZONAS_MIMAG, zonaMIMAG, COR_ZONA } from '../charts/Scatter.jsx';
import { Herdado, Medido } from '../viz/Provenance.jsx';
import { useReport, TODAS } from '../state/store.jsx';

function grupoDeComprimento(bloco, nome) {
  if (!bloco) return null;
  return { name: nome, values: bloco.values, bins: bloco.bins, n: bloco.n };
}

function Individual({ amostra, assumida }) {
  const { data } = useReport();
  const seq = data?.sequencing;
  const prok = data?.prokaryotic;

  const status = (data?.overview?.status ?? []).filter((r) => r.sample === amostra);
  const funil = data?.overview?.funnel?.[amostra];
  const qc = (seq?.qc ?? []).find((r) => r.sample === amostra);
  const mapping = (seq?.mapping ?? []).find((r) => r.sample === amostra);
  const comprimentos = grupoDeComprimento(seq?.lengths?.[amostra], amostra);

  const mags = (prok?.quality ?? []).filter((r) => r.source === amostra);
  const taxPorGenoma = {};
  (prok?.taxonomy ?? []).forEach((r) => { taxPorGenoma[r.genome] = r; });

  const alpha = data?.diversity?.alpha ?? [];
  const alfaDaAmostra = alpha.find(
    (r) => r.sample === amostra && r.index === 'shannon');
  const outrasAlfa = alpha.filter(
    (r) => r.index === 'shannon' && r.sample !== amostra);

  return (
    <div data-testid="sample-view" data-sample={amostra}>
      {assumida ? (
        <p className="notice">
          Nenhuma amostra selecionada no filtro do topo — mostrando{' '}
          <strong>{amostra}</strong>, a primeira da rodada.
        </p>
      ) : null}

      <section className="card" data-testid="sample-status">
        <h2>Execução</h2>
        {status.length ? (
          <>
            <StatusMatrix rows={status} />
            {/* O motivo escrito, e nao so o glifo da matriz: numa vista de UMA
                amostra, "gunc falhou" sem o porquê deixa o leitor concluir que
                a amostra não tem quimeras -- ferramenta que falhou é lacuna. */}
            {status.filter((r) => r.status === 'failed' || r.status === 'skipped')
              .map((r) => (
                <p key={r.rule} className="notice">
                  <code>{r.rule}</code> — {r.status === 'failed' ? 'falhou' : 'pulada'}
                  {r.reason ? `: ${r.reason}` : ''}. Os números desta regra são
                  {' '}lacuna, não zero biológico.
                </p>
              ))}
          </>
        ) : <p className="empty">Sem registro de execução para esta amostra.</p>}
      </section>

      <section className="card" data-testid="sample-sequencing">
        <h2>Sequenciamento e montagem <Medido o="QC, montagem e mapeamento" /></h2>
        <div className="kpi-row">
          <StatTile label="Reads antes do fastp" value={qc?.reads_before ?? '—'} />
          <StatTile label="Reads depois" value={qc?.reads_after ?? '—'} />
          <StatTile label="Q30" value={qc?.q30 !== undefined ? `${qc.q30}%` : '—'} />
          <StatTile label="Taxa de mapeamento"
                    value={mapping?.rate !== undefined ? `${mapping.rate}%` : '—'} />
        </div>
        {comprimentos ? (
          <DistPlot groups={[comprimentos]} xName="comprimento do contig (bp)" log />
        ) : <p className="empty">Sem distribuição de comprimento para esta amostra.</p>}
      </section>

      <section className="card">
        <h2>Viral <Medido o="detecção e retenção" /></h2>
        {funil?.stages?.length ? (
          <AttritionFunnel stages={funil.stages} losses={funil.losses} />
        ) : <p className="empty">Sem funil viral para esta amostra.</p>}
      </section>

      <section className="card">
        <h2>MAGs binados aqui</h2>
        <p className="chart__sub">
          CheckM2 e GUNC são medidos nesta amostra. Taxonomia GTDB e tudo o
          que vem depois do binning saem do representante do cluster.
        </p>
        {mags.length ? (
          <>
            <Scatter
              testid="sample-mag-quality"
              points={mags.map((m) => ({ id: m.genome, x: m.completeness,
                                         y: m.contamination }))}
              zones={ZONAS_MIMAG}
              xName="Completude (%)" yName="Contaminação (%)" xDomain={[0, 100]}
              colorOf={(p) => COR_ZONA[zonaMIMAG(p.x, p.y)]}
              tooltipOf={(p) => `${p.id}: ${p.x}% completo, ${p.y}% contaminado`}
            />
            <div className="table-wrap">
              <table className="table">
                <thead>
                  <tr><th>MAG</th><th>Cluster</th><th>Taxonomia GTDB</th></tr>
                </thead>
                <tbody>
                  {mags.map((m) => {
                    const tax = taxPorGenoma[m.genome];
                    const linhagem = tax
                      ? [tax.Phylum, tax.Class, tax.Order, tax.Family, tax.Genus]
                          .filter(Boolean).join(' · ')
                      : '—';
                    return (
                      <tr key={m.genome} data-testid={`mag-${m.genome}`}>
                        <td><code>{m.genome}</code></td>
                        <td>
                          {m.is_representative
                            ? <>representante <Medido o="binning" /></>
                            : <code>{m.representative}</code>}
                        </td>
                        <td>
                          {linhagem}{' '}
                          {tax?.inherited
                            ? <Herdado de={tax.representative} o="taxonomia GTDB" />
                            : null}
                        </td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            </div>
          </>
        ) : <p className="empty">Nenhum MAG foi binado nesta amostra.</p>}
      </section>

      <section className="card">
        <h2>Diversidade <Medido o="diversidade alfa" /></h2>
        {alfaDaAmostra ? (
          <div data-testid="sample-alpha"
               data-value={String(alfaDaAmostra.value)}
               data-n-outras={String(outrasAlfa.length)}>
            <p className="chart__sub">
              Shannon desta amostra contra a distribuição das demais — um valor
              de alfa isolado não significa nada; a posição relativa é o que
              lhe dá sentido.
            </p>
            <DistPlot
              groups={[{ name: 'demais amostras',
                         values: outrasAlfa.map((r) => r.value) }]}
              xName="Shannon"
              cutoffs={[{ label: amostra, value: alfaDaAmostra.value }]}
            />
          </div>
        ) : <p className="empty">Sem diversidade alfa para esta amostra.</p>}
      </section>
    </div>
  );
}

function Comparacao({ selecionadas, alternar }) {
  const { data, samples } = useReport();
  const seq = data?.sequencing;
  const grupos = selecionadas
    .map((s) => grupoDeComprimento(seq?.lengths?.[s], s))
    .filter(Boolean);

  return (
    <div data-testid="compare-view" data-samples={selecionadas.join(',')}>
      <section className="card">
        <h2>Amostras na comparação</h2>
        <div className="rank-selector" role="group" aria-label="Amostras">
          {samples.map((s) => (
            <button key={s} type="button" className="rank-selector__btn"
                    aria-pressed={selecionadas.includes(s)}
                    data-active={selecionadas.includes(s) ? 'true' : 'false'}
                    onClick={() => alternar(s)}>{s}</button>
          ))}
        </div>
      </section>

      <section className="card">
        <h2>Comprimento de contig</h2>
        <p className="chart__sub">
          Pequenos múltiplos no <strong>mesmo eixo</strong>: é o alinhamento
          que torna a comparação legível. Hierarquia taxonômica não entra aqui
          — ela não se compara lado a lado.
        </p>
        {grupos.length ? (
          <div data-testid="compare-lengths">
            <DistPlot groups={grupos} xName="comprimento do contig (bp)" log />
          </div>
        ) : <p className="empty">Sem distribuição de comprimento nas amostras escolhidas.</p>}
      </section>
    </div>
  );
}

export function SampleView() {
  const { data, sample, samples } = useReport();
  const [modo, setModo] = useState('individual');
  const [extras, setExtras] = useState([]);

  if (!samples.length) {
    return <p className="empty">Esta rodada não tem amostras.</p>;
  }

  // O filtro global manda; TODAS não é uma amostra, então a vista assume a
  // primeira e DIZ que assumiu -- ficar vazia faria o usuário achar que a
  // aba quebrou.
  const assumida = sample === TODAS;
  const amostra = assumida ? samples[0] : sample;

  const selecionadas = [amostra, ...extras.filter((s) => s !== amostra)];
  const alternar = (s) => setExtras((atuais) => (
    atuais.includes(s) ? atuais.filter((x) => x !== s) : [...atuais, s]));

  return (
    <div className="panel">
      <div className="rank-selector" role="group" aria-label="Modo">
        <button type="button" className="rank-selector__btn"
                aria-pressed={modo === 'individual'}
                data-active={modo === 'individual' ? 'true' : 'false'}
                onClick={() => setModo('individual')}>Individual</button>
        <button type="button" className="rank-selector__btn"
                aria-pressed={modo === 'comparacao'}
                data-active={modo === 'comparacao' ? 'true' : 'false'}
                onClick={() => setModo('comparacao')}>Comparar</button>
      </div>

      {modo === 'individual'
        ? <Individual amostra={amostra} assumida={assumida} />
        : <Comparacao selecionadas={selecionadas} alternar={alternar} />}
    </div>
  );
}
