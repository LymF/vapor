import { useCallback, useState } from 'react';

const MARGEM = 12;
// Chute de largura antes do primeiro layout -- soh usado para o clamp da
// borda direita na primeira renderizacao; apos montar, o proprio ref mede.
const LARGURA_ESTIMADA = 220;

// Nó único e absoluto: tooltip nunca é criado/destruído por gráfico, e nunca
// entra na árvore de eventos (pointer-events:none) para não disparar
// mouseenter/mouseleave sobre si mesmo (loop de entra/sai descrito no brief).
export function Tooltip({ content, x = 0, y = 0 }) {
  if (content == null) return null;
  const larguraJanela = typeof window !== 'undefined' ? window.innerWidth : 1200;
  const alturaJanela = typeof window !== 'undefined' ? window.innerHeight : 800;
  const left = Math.min(x + MARGEM, Math.max(larguraJanela - LARGURA_ESTIMADA - MARGEM, MARGEM));
  const top = Math.min(y + MARGEM, alturaJanela - MARGEM);

  return (
    <div className="tooltip" style={{ left, top }}>
      {content}
    </div>
  );
}

export function useTooltip() {
  const [estado, setEstado] = useState(null);

  const show = useCallback((evt, content) => {
    const x = evt?.clientX ?? 0;
    const y = evt?.clientY ?? 0;
    setEstado({ x, y, content });
  }, []);

  const hide = useCallback(() => setEstado(null), []);

  const node = estado
    ? <Tooltip content={estado.content} x={estado.x} y={estado.y} />
    : null;

  return { show, hide, node };
}
