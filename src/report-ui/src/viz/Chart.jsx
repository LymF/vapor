import { useRef } from 'react';
import { useResize } from './useResize.js';

// Largura de fallback usada antes do ResizeObserver disparar (jsdom nunca
// dispara; no navegador real o primeiro paint tambem chega com width:0). Sem
// isto o filho criaria escalas com range [0,0] e dividiria por zero.
const LARGURA_FALLBACK = 640;

export function Chart({ title, sub, height = 280, empty = false, emptyLabel = 'Sem dados nesta rodada', children }) {
  const ref = useRef(null);
  const { width } = useResize(ref);
  const larguraEfetiva = width > 0 ? width : LARGURA_FALLBACK;

  return (
    <figure className="chart" ref={ref}>
      {title ? <figcaption className="chart__title">{title}</figcaption> : null}
      {sub ? <p className="chart__sub">{sub}</p> : null}
      {empty ? (
        <p className="empty">{emptyLabel}</p>
      ) : (
        <svg role="img" aria-label={title || 'grafico'} width={larguraEfetiva} height={height}>
          {children({ width: larguraEfetiva, height })}
        </svg>
      )}
    </figure>
  );
}
