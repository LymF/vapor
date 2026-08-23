// Aba Leituras — perfilagem independente de montagem (sylph + sylph-tax).
//
// O aviso do topo não é rodapé: o sylph identifica GENOMAS DE REFERÊNCIA do
// IMG/VR (`t__IMGVR_UViG_…`) e o resto do report identifica contigs montados
// (`k141_…`). São dois espaços de ID sem interseção, e qualquer cruzamento
// entre as duas trilhas seria inválido — dizê-lo na cara do usuário é mais
// barato que descobrir depois que alguém juntou as tabelas.
import { StackedBar } from '../charts/BarChart.jsx';
import { StatTile } from '../charts/StatTile.jsx';
import { useReport } from '../state/store.jsx';

// Procedência do hospedeiro. 'db' = anotação do banco; 'none' = nenhuma.
// Sem essa distinção, um gênero sem hospedeiro conhecido e um cujo
// hospedeiro foi anotado ficam idênticos no gráfico.
const ROTULO_FONTE = {
  db: 'anotado no banco',
  none: 'sem hospedeiro conhecido',
  '': 'procedência não registrada',
};

function ultimoNivel(clade) {
  const partes = String(clade || '').split('|');
  return partes[partes.length - 1] || clade;
}

function composicaoPorAmostra(linhas, samples) {
  return samples.map((s) => ({
    name: s,
    parts: Object.fromEntries(
      linhas.map((r) => [ultimoNivel(r.clade), r[s] || 0]).filter(([, v]) => v > 0),
    ),
  }));
}

export function Reads() {
  const { data } = useReport();
  const bloco = data?.reads;
  const samples = data?.run?.samples ?? [];
  if (!bloco) {
    return <p className="empty">Sem dados para esta aba nesta rodada.</p>;
  }

  const viral = bloco.viral ?? [];
  const prok = bloco.prok ?? [];
  const host = bloco.host ?? [];
  const truncados = bloco.truncated ?? {};
  const totalTruncado = Object.values(truncados).reduce((a, b) => a + b, 0);

  return (
    <div className="panel">
      <p className="notice">{bloco.id_space_warning}</p>

      <div className="kpi-row">
        <StatTile label="Táxons virais" value={viral.length} />
        <StatTile label="Táxons procarióticos" value={prok.length} />
        <StatTile label="Gêneros hospedeiros" value={host.length} />
      </div>

      {totalTruncado ? (
        <p className="chart__sub">
          Mostrando os táxons mais abundantes; {totalTruncado} ficaram de fora
          do report para caber no orçamento de payload — a tabela completa está
          em <code>reads_classify/otu_table.tsv</code>.
        </p>
      ) : null}

      <section className="card">
        <h2>Composição viral</h2>
        {viral.length ? (
          <StackedBar normalize data={composicaoPorAmostra(viral, samples)} />
        ) : <p className="empty">Nenhum táxon viral sobreviveu ao filtro de prevalência.</p>}
      </section>

      <section className="card">
        <h2>Composição procariótica</h2>
        {prok.length ? (
          <StackedBar normalize data={composicaoPorAmostra(prok, samples)} />
        ) : <p className="empty">Nenhum táxon procariótico nesta rodada.</p>}
      </section>

      <section className="card">
        <h2>Hospedeiro dos vírus detectados</h2>
        {host.length ? (
          <>
            <p className="chart__sub">
              O hospedeiro vem da anotação do <strong>banco</strong> (coluna
              <code> Virus_host</code> do IMG/VR), não de predição feita nesta
              rodada — nem BACPHLIP nem PHIST rodam nesta trilha, porque só
              existem os esboços k-mer do sylph, e um esboço não devolve
              sequência.
            </p>
            <div className="table-wrap">
              <table className="table">
                <thead>
                  <tr><th>Gênero hospedeiro</th><th>Procedência</th>
                    <th>Táxons virais</th></tr>
                </thead>
                <tbody>
                  {host.map((h) => (
                    <tr key={h.host_genus}
                        data-testid={`host-${h.host_genus}`}
                        data-source={h.host_source}>
                      <td>{h.host_genus}</td>
                      <td>{ROTULO_FONTE[h.host_source] ?? h.host_source}</td>
                      <td>{h.n_viral_taxa}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </>
        ) : <p className="empty">Sem colapso por hospedeiro nesta rodada.</p>}
      </section>
    </div>
  );
}
