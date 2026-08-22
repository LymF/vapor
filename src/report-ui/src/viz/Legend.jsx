// Cor sozinha nao identifica nada (daltonismo, impressao P&B, alto contraste):
// o rotulo textual e obrigatorio, a cor e so o acompanhamento.
export function Legend({ items = [] }) {
  if (!items.length) return null;
  return (
    <ul className="legend">
      {items.map((item) => (
        <li key={item.label} className="legend__item">
          <span className="legend__swatch" style={{ backgroundColor: item.color }} aria-hidden="true" />
          <span>{item.label}</span>
        </li>
      ))}
    </ul>
  );
}
