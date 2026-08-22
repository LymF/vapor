// Os QUATRO estados de load_tool_status, visualmente distintos. Ferramenta que
// falhou tem de ler como LACUNA, jamais como contagem zero: um done.txt vazio
// ja fez uma rodada de AMRFinderPlus com disco cheio passar por zero
// biologico. 'unknown' (done.txt ausente ou vazio) entra no mesmo grupo: nao e
// sucesso.
const ROTULO = {
  ok: 'ok',
  skipped: 'pulado',
  failed: 'falhou',
  unknown: 'sem status registrado',
};
const GLIFO = { ok: '●', skipped: '–', failed: '✕', unknown: '?' };

export function StatusMatrix({ rows }) {
  const regras = [...new Set(rows.map((r) => r.rule))];
  const amostras = [...new Set(rows.map((r) => r.sample))];
  const porChave = new Map(rows.map((r) => [`${r.rule}||${r.sample}`, r]));

  return (
    <table className="status-matrix">
      <thead>
        <tr>
          <th scope="col">regra</th>
          {amostras.map((s) => <th scope="col" key={s}>{s}</th>)}
        </tr>
      </thead>
      <tbody>
        {regras.map((regra) => (
          <tr key={regra}>
            <th scope="row">{regra}</th>
            {amostras.map((s) => {
              const cel = porChave.get(`${regra}||${s}`);
              const status = cel?.status ?? 'unknown';
              const motivo = cel?.reason ? ` — ${cel.reason}` : '';
              return (
                <td
                  key={s}
                  data-testid={`cell-${regra}-${s}`}
                  data-status={status}
                  className={`status-matrix__cell status-matrix__cell--${status}`}
                  aria-label={`${regra} em ${s}: ${ROTULO[status]}${motivo}`}
                  title={`${ROTULO[status]}${motivo}`}
                >
                  <span className="status-matrix__glyph" aria-hidden="true">
                    {GLIFO[status] ?? '?'}
                  </span>
                </td>
              );
            })}
          </tr>
        ))}
      </tbody>
    </table>
  );
}
