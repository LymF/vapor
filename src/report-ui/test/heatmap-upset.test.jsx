import { render, screen } from '@testing-library/react';
import { Heatmap } from '../src/charts/Heatmap.jsx';
import { UpSet } from '../src/charts/UpSet.jsx';

test('matriz esparsa vira bolha', () => {
  render(<Heatmap rows={['a','b','c','d','e']} cols={['x','y','z','w']}
                  values={{ a: { x: 1 }, b: { y: 2 } }} sparseAsBubble />);
  expect(screen.getByTestId('heatmap').getAttribute('data-mode')).toBe('bubble');
});

test('matriz densa fica heatmap', () => {
  const vals = Object.fromEntries(['a','b'].map(r => [r, { x: 1, y: 2 }]));
  render(<Heatmap rows={['a','b']} cols={['x','y']} values={vals} sparseAsBubble />);
  expect(screen.getByTestId('heatmap').getAttribute('data-mode')).toBe('heat');
});

test('normalizacao por coluna nao mistura unidades', () => {
  render(<Heatmap rows={['a','b']} cols={['n50','contigs']}
                  values={{ a: { n50: 1000, contigs: 10 }, b: { n50: 2000, contigs: 5 } }}
                  normalize="per-col" />);
  const cel = (r, c) => document.querySelector(`[data-cell="${r}|${c}"]`).getAttribute('data-norm');
  expect(Number(cel('b','n50'))).toBeCloseTo(1, 2);
  expect(Number(cel('a','contigs'))).toBeCloseTo(1, 2);
});

test('UpSet ordena as intersecoes por tamanho', () => {
  render(<UpSet sets={{ geNomad: 10, VirSorter2: 8 }}
                combos={[{ tools: ['geNomad'], count: 3 },
                         { tools: ['geNomad','VirSorter2'], count: 7 }]} valueName="contigs" />);
  const barras = [...document.querySelectorAll('[data-combo]')].map(b => b.getAttribute('data-combo'));
  expect(barras[0]).toBe('geNomad+VirSorter2');
});
