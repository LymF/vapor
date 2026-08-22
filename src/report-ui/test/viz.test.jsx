import { test, expect } from 'vitest';
import { PAL, PAL_MUTED, foldOther } from '../src/viz/palette.js';
import { pickDistributionForm, pickAxisOrientation, TRIGGERS } from '../src/viz/triggers.js';
import { render } from '@testing-library/react';
import { useRef } from 'react';
import { useResize } from '../src/viz/useResize.js';

test('a paleta validada tem exatamente 8 cores', () => {
  expect(PAL).toHaveLength(8);
  expect(PAL).toContain('#0d9488');
});

test('foldOther dobra a cauda e poe Other em ultimo', () => {
  const counts = { a: 10, b: 9, c: 8, d: 7, e: 6, f: 5, g: 4, h: 3, i: 2, j: 1 };
  const out = foldOther(counts, 8);
  expect(out).toHaveLength(8);
  expect(out[7][0]).toBe('Other');
  expect(out[7][1]).toBe(2 + 1 + 3);
});

test('foldOther nao inventa Other quando cabe', () => {
  expect(foldOther({ a: 1, b: 2 }, 8).map((r) => r[0])).toEqual(['b', 'a']);
});

test('poucos pontos viram strip plot, muitos viram densidade', () => {
  expect(pickDistributionForm(TRIGGERS.densityMinN - 1)).toBe('strip');
  expect(pickDistributionForm(TRIGGERS.densityMinN)).toBe('density');
});

test('muitas amostras viram eixo horizontal', () => {
  expect(pickAxisOrientation(TRIGGERS.manySamples)).toBe('vertical');
  expect(pickAxisOrientation(TRIGGERS.manySamples + 1)).toBe('horizontal');
});

function Sonda() {
  const ref = useRef(null);
  const { width } = useResize(ref);
  return <div ref={ref} data-testid="alvo" data-w={width} />;
}

test('useResize comeca em zero e nao quebra sem ResizeObserver real', () => {
  const { getByTestId } = render(<Sonda />);
  expect(getByTestId('alvo').getAttribute('data-w')).toBe('0');
});
