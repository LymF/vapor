const fmt = new Intl.NumberFormat('pt-BR');

export function StatTile({ label, value, sub }) {
  const texto = typeof value === 'number' ? fmt.format(value) : value;
  return (
    <div className="stat-tile">
      <span className="stat-tile__label">{label}</span>
      <span className="stat-tile__value">{texto}</span>
      {sub ? <span className="stat-tile__sub">{sub}</span> : null}
    </div>
  );
}
