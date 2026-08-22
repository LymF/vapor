import { useEffect, useState } from 'react';

export function useResize(ref) {
  const [tamanho, setTamanho] = useState({ width: 0, height: 0 });

  useEffect(() => {
    const alvo = ref.current;
    if (!alvo || typeof ResizeObserver === 'undefined') return undefined;
    const obs = new ResizeObserver(([entrada]) => {
      const { width, height } = entrada.contentRect;
      setTamanho({ width, height });
    });
    obs.observe(alvo);
    return () => obs.disconnect();
  }, [ref]);

  return tamanho;
}
