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

function nenhumAtributoTemNaN() {
  const todos = document.querySelectorAll('svg, svg *');
  for (const el of todos) {
    for (const attr of el.attributes) {
      if (String(attr.value).includes('NaN')) return false;
    }
  }
  return true;
}

test('length ausente nao produz NaN em nenhum atributo do svg', () => {
  render(<GenomeTrack length={undefined} features={feats} kinds={{}} />);
  expect(nenhumAtributoTemNaN()).toBe(true);
  expect(document.querySelectorAll('[data-feature]').length).toBe(0);
});

test('length zero recebe o mesmo tratamento de length ausente', () => {
  render(<GenomeTrack length={0} features={feats} kinds={{}} />);
  expect(nenhumAtributoTemNaN()).toBe(true);
  expect(document.querySelectorAll('[data-feature]').length).toBe(0);
});

test('start > end com fita + ainda aponta para a direita, igual a uma feature bem formada', () => {
  const invertida = [{ start: 900, end: 1, strand: '+', label: 'invertida', kind: 'phrog' }];
  const bemFormada = [{ start: 1, end: 900, strand: '+', label: 'normal', kind: 'phrog' }];
  const { unmount } = render(<GenomeTrack length={3000} features={invertida} kinds={{}} />);
  const pontosInvertida = document.querySelector('[data-feature="invertida"]').getAttribute('points');
  unmount();
  render(<GenomeTrack length={3000} features={bemFormada} kinds={{}} />);
  const pontosNormal = document.querySelector('[data-feature="normal"]').getAttribute('points');
  expect(pontosInvertida).toBe(pontosNormal);
});

test('start > end com fita - aponta para a esquerda', () => {
  const invertida = [{ start: 900, end: 1, strand: '-', label: 'invertida-neg', kind: 'phrog' }];
  render(<GenomeTrack length={3000} features={invertida} kinds={{}} />);
  const el = document.querySelector('[data-feature="invertida-neg"]');
  expect(el.getAttribute('data-strand')).toBe('-');
  const pontos = el.getAttribute('points').split(' ').map((p) => p.split(',').map(Number));
  // Fita '-': primeiro ponto tem x maior que o ultimo ponto (a ponta da seta
  // fica no lado esquerdo/menor x) -- oposto de uma seta '+'.
  const primeiroX = pontos[0][0];
  const ultimoX = pontos[pontos.length - 1][0];
  expect(primeiroX).toBeGreaterThan(ultimoX);
});
