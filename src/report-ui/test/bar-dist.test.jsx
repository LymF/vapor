import { render, screen } from '@testing-library/react';
import { BarChart, StackedBar } from '../src/charts/BarChart.jsx';
import { DistPlot } from '../src/charts/DistPlot.jsx';

const muitas = Array.from({ length: 15 }, (_, i) => ({ name: `S${i}`, value: i + 1 }));

test('poucas amostras ficam na vertical', () => {
  render(<BarChart data={[{ name: 'a', value: 1 }]} valueName="n" />);
  expect(screen.getByTestId('barchart').getAttribute('data-orientation')).toBe('vertical');
});

test('muitas amostras viram horizontal', () => {
  render(<BarChart data={muitas} valueName="n" />);
  expect(screen.getByTestId('barchart').getAttribute('data-orientation')).toBe('horizontal');
});

test('rotulo longo tambem vira horizontal', () => {
  render(<BarChart data={[{ name: 'Pseudomonadota_muito_longo', value: 3 }]} valueName="n" />);
  expect(screen.getByTestId('barchart').getAttribute('data-orientation')).toBe('horizontal');
});

test('a barra parte do zero', () => {
  render(<BarChart data={[{ name: 'a', value: 5 }, { name: 'b', value: 10 }]} valueName="n" />);
  const barras = [...document.querySelectorAll('[data-bar]')];
  const razao = Number(barras[1].getAttribute('data-len')) / Number(barras[0].getAttribute('data-len'));
  expect(razao).toBeCloseTo(2, 1);   // zero-baseline: 10 e o dobro de 5
});

test('mais de 8 categorias dobram em Other', () => {
  const partes = Object.fromEntries(Array.from({ length: 12 }, (_, i) => [`c${i}`, 1]));
  render(<StackedBar data={[{ name: 'S1', parts: partes }]} />);
  expect(screen.getByTestId('stacked').getAttribute('data-series')).toBe('8');
});

test('poucos pontos viram strip, muitos viram densidade', () => {
  const poucos = { name: 'g', values: [1, 2, 3] };
  const muitos = { name: 'g', values: Array.from({ length: 40 }, (_, i) => i) };
  const { rerender } = render(<DistPlot groups={[poucos]} xName="x" />);
  expect(screen.getByTestId('distplot').getAttribute('data-form')).toBe('strip');
  rerender(<DistPlot groups={[muitos]} xName="x" />);
  expect(screen.getByTestId('distplot').getAttribute('data-form')).toBe('density');
});

test('muitos grupos viram ridgeline', () => {
  const g = Array.from({ length: 10 }, (_, i) =>
    ({ name: `S${i}`, values: Array.from({ length: 40 }, (_, j) => j) }));
  render(<DistPlot groups={g} xName="x" />);
  expect(screen.getByTestId('distplot').getAttribute('data-form')).toBe('ridgeline');
});

test('linhas de corte aparecem com rotulo', () => {
  render(<DistPlot groups={[{ name: 'g', values: [1, 2, 3] }]} xName="x"
                   cutoffs={[{ value: 2, label: 'MQ' }]} />);
  expect(screen.getByText('MQ')).toBeTruthy();
});
