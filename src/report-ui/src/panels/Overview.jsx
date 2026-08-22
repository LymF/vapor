import { StatTile } from '../charts/StatTile.jsx';
import { StatusMatrix } from '../charts/StatusMatrix.jsx';
import { AttritionFunnel } from '../charts/AttritionFunnel.jsx';
import { useReport, TODAS } from '../state/store.jsx';

export function Overview() {
  const { data, sample } = useReport();
  const ov = data?.overview ?? {};
  const kpis = ov.kpis ?? [];
  const status = ov.status ?? [];
  const funil = (ov.funnel ?? {})[sample] ?? (ov.funnel ?? {})[TODAS] ?? null;

  const temFunil = Boolean(funil?.stages?.length);

  if (!kpis.length && !status.length && !temFunil) {
    return <p className="empty">Sem dados para esta aba nesta rodada.</p>;
  }

  return (
    <div className="panel">
      <div className="kpi-row">
        {kpis.map((k) => <StatTile key={k.label} {...k} />)}
      </div>

      {temFunil ? (
        <section className="card">
          <h2>Atrição da rodada</h2>
          <p className="card__scope">
            escopo: <span data-testid="funnel-scope">{sample === TODAS ? 'todas as amostras' : sample}</span>
          </p>
          <AttritionFunnel stages={funil.stages} losses={funil.losses} />
        </section>
      ) : null}

      {status.length ? (
        <section className="card">
          <h2>Status das ferramentas</h2>
          <StatusMatrix rows={sample === TODAS ? status : status.filter((r) => r.sample === sample)} />
        </section>
      ) : null}
    </div>
  );
}
