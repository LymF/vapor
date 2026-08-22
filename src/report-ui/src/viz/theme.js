import { createElement, useCallback, useEffect, useState } from 'react';

const CHAVE = 'vapor-report-theme';

function lido() {
  // localStorage lanca em contexto sem storage (captura de thumbnail, aba
  // privada, navegador com dados de site bloqueados) -- ler tem de ser seguro.
  try {
    const v = window.localStorage.getItem(CHAVE);
    return v === 'dark' || v === 'light' ? v : null;
  } catch { return null; }
}

function grava(v) {
  try { window.localStorage.setItem(CHAVE, v); } catch { /* sem storage */ }
}

export function useTheme() {
  const [theme, setTheme] = useState(() => lido() ?? 'light');

  useEffect(() => {
    document.documentElement.setAttribute('data-theme', theme);
  }, [theme]);

  const toggle = useCallback(() => {
    setTheme((t) => {
      const novo = t === 'dark' ? 'light' : 'dark';
      grava(novo);
      return novo;
    });
  }, []);

  return { theme, toggle };
}

export function ThemeToggle() {
  const { theme, toggle } = useTheme();
  return createElement(
    'button',
    {
      className: 'theme-toggle',
      onClick: toggle,
      'aria-label': `Alternar tema (atual: ${theme === 'dark' ? 'escuro' : 'claro'})`,
    },
    theme === 'dark' ? '☾' : '☀',
  );
}
