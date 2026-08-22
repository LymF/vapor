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

test('ridgeline com poucos pontos por grupo: todas as faixas em strip', () => {
  const g = Array.from({ length: 10 }, (_, i) =>
    ({ name: `S${i}`, values: [1, 2, 3] }));
  render(<DistPlot groups={g} xName="x" />);
  expect(screen.getByTestId('distplot').getAttribute('data-form')).toBe('ridgeline');
  const faixas = [...document.querySelectorAll('[data-lane]')];
  expect(faixas).toHaveLength(10);
  faixas.forEach((f) => expect(f.getAttribute('data-lane-form')).toBe('strip'));
});

test('ridgeline com muitos pontos por grupo: todas as faixas em densidade', () => {
  const g = Array.from({ length: 10 }, (_, i) =>
    ({ name: `S${i}`, values: Array.from({ length: 40 }, (_, j) => j) }));
  render(<DistPlot groups={g} xName="x" />);
  expect(screen.getByTestId('distplot').getAttribute('data-form')).toBe('ridgeline');
  const faixas = [...document.querySelectorAll('[data-lane]')];
  expect(faixas).toHaveLength(10);
  faixas.forEach((f) => expect(f.getAttribute('data-lane-form')).toBe('density'));
});

test('ridgeline com grupos mistos: cada faixa escolhe pelo seu proprio n', () => {
  const poucos = Array.from({ length: 5 }, (_, i) => ({ name: `pouco${i}`, values: [1, 2, 3] }));
  const muitos = Array.from({ length: 5 }, (_, i) => ({ name: `muito${i}`, values: Array.from({ length: 40 }, (_, j) => j) }));
  render(<DistPlot groups={[...poucos, ...muitos]} xName="x" />);
  expect(screen.getByTestId('distplot').getAttribute('data-form')).toBe('ridgeline');
  poucos.forEach((g) => {
    expect(document.querySelector(`[data-lane="${g.name}"]`).getAttribute('data-lane-form')).toBe('strip');
  });
  muitos.forEach((g) => {
    expect(document.querySelector(`[data-lane="${g.name}"]`).getAttribute('data-lane-form')).toBe('density');
  });
});

test('DistPlot com bins desenha histograma sem chamar KDE e mostra n no rodape', () => {
  // kde() so produz <path> (a curva suavizada); o caminho de bins so produz
  // <rect> (barras a partir das contagens prontas). Group sem `values` faz
  // kde(g.values, ...) explodir se fosse chamado (values e undefined) --
  // entao o proprio render sem excecao, mais a ausencia de <path>, ja prova
  // que o KDE nunca roda para este grupo.
  const grupo = {
    name: 'S1',
    n: 24680,
    min: 1000,
    max: 400000,
    bins: [
      { x0: 1000, x1: 2000, count: 12000 },
      { x0: 2000, x1: 4000, count: 8000 },
      { x0: 4000, x1: 400000, count: 4680 },
    ],
  };
  render(<DistPlot groups={[grupo]} xName="comprimento" />);

  expect(screen.getByTestId('distplot').getAttribute('data-form')).toBe('density');
  expect(document.querySelector('[data-group-form="histogram"]')).toBeTruthy();
  expect(document.querySelectorAll('rect[fill]').length).toBeGreaterThan(0);
  expect(document.querySelector('path')).toBeNull();   // nenhuma curva de KDE desenhada
  expect(screen.getByTestId('distplot-n').textContent).toContain('24.680');
});

test('DistPlot com values continua identico (strip/densidade por pontos crus)', () => {
  const poucos = { name: 'g', values: [1, 2, 3] };
  const muitos = { name: 'g', values: Array.from({ length: 40 }, (_, i) => i) };
  const { rerender } = render(<DistPlot groups={[poucos]} xName="x" />);
  expect(screen.getByTestId('distplot').getAttribute('data-form')).toBe('strip');
  expect(screen.queryByTestId('distplot-n')).toBeNull(); // sem bins, sem rodape de n

  rerender(<DistPlot groups={[muitos]} xName="x" />);
  expect(screen.getByTestId('distplot').getAttribute('data-form')).toBe('density');
  expect(document.querySelector('path')).toBeTruthy();   // curva de densidade, nao barras
  expect(document.querySelector('[data-group-form="histogram"]')).toBeNull();
});

test('StackedBar: order sem Other nao faz a cauda dobrada sumir', () => {
  const partes = Object.fromEntries(Array.from({ length: 12 }, (_, i) => [`c${i}`, 12 - i]));
  render(<StackedBar data={[{ name: 'S1', parts: partes }]} order={['c0', 'c1', 'c2']} />);
  expect(screen.getByTestId('stacked').getAttribute('data-series')).toBe('4');
  expect(document.querySelector('[data-part="Other"]')).toBeTruthy();
});
