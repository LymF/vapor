import { render } from '@testing-library/react';
import { GenomeTrack } from '../src/charts/GenomeTrack.jsx';

const feats = [
  { start: 1, end: 900, strand: '+', label: 'terminase', kind: 'phrog' },
  { start: 1200, end: 1500, strand: '-', label: 'anti-CBASS', kind: 'antidefense' },
];

test('a largura da feature e proporcional ao seu tamanho em bp', () => {
  render(<GenomeTrack length={3000} features={feats} kinds={{}} />);
  const a = Number(document.querySelector('[data-feature="terminase"]').getAttribute('data-w'));
  const b = Number(document.querySelector('[data-feature="anti-CBASS"]').getAttribute('data-w'));
  expect(a / b).toBeCloseTo(900 / 300, 1);
});

test('a fita determina o sentido da seta', () => {
  render(<GenomeTrack length={3000} features={feats} kinds={{}} />);
  expect(document.querySelector('[data-feature="anti-CBASS"]').getAttribute('data-strand')).toBe('-');
});

test('features sobrepostas vao para faixas diferentes', () => {
  const sobre = [{ start: 1, end: 1000, strand: '+', label: 'a', kind: 'x' },
                 { start: 500, end: 1500, strand: '+', label: 'b', kind: 'x' }];
  render(<GenomeTrack length={2000} features={sobre} kinds={{}} />);
  const fa = document.querySelector('[data-feature="a"]').getAttribute('data-lane');
  const fb = document.querySelector('[data-feature="b"]').getAttribute('data-lane');
  expect(fa).not.toBe(fb);
});

test('sem features desenha so a regua', () => {
  render(<GenomeTrack length={1000} features={[]} kinds={{}} />);
  expect(document.querySelectorAll('[data-feature]').length).toBe(0);
});
