import { render, screen, fireEvent } from '@testing-library/react';
import { AttritionFunnel } from '../src/charts/AttritionFunnel.jsx';

const stages = [
  { name: 'contigs', value: 1000 },
  { name: 'candidatos virais', value: 400 },
  { name: 'vOTUs retidos', value: 120 },
];
const losses = {
  'candidatos virais': [{ reason: 'sem sinal viral', count: 600 }],
  'vOTUs retidos': [
    { reason: 'curto e sem bin', count: 200 },
    { reason: 'CheckV baixo', count: 80 },
  ],
};

test('desenha uma barra por etapa', () => {
  render(<AttritionFunnel stages={stages} losses={losses} />);
  stages.forEach((s) => expect(screen.getByTestId(`stage-${s.name}`)).toBeTruthy());
});

test('a perda de cada etapa e a diferenca para a etapa anterior', () => {
  render(<AttritionFunnel stages={stages} losses={losses} />);
  expect(screen.getByTestId('loss-vOTUs retidos').getAttribute('data-loss')).toBe('280');
});

test('clicar na perda emite a etapa', () => {
  const visto = [];
  render(<AttritionFunnel stages={stages} losses={losses} onSelectLoss={(n) => visto.push(n)} />);
  fireEvent.click(screen.getByTestId('loss-vOTUs retidos'));
  expect(visto).toEqual(['vOTUs retidos']);
});

test('a primeira etapa nao tem perda', () => {
  render(<AttritionFunnel stages={stages} losses={losses} />);
  expect(screen.queryByTestId('loss-contigs')).toBeNull();
});

test('unidades iguais desenham a barra de perda', () => {
  const comUnidade = [
    { name: 'contigs', value: 1000, unit: 'contig' },
    { name: 'candidatos virais', value: 400, unit: 'contig' },
  ];
  render(<AttritionFunnel stages={comUnidade} />);
  expect(screen.getByTestId('loss-candidatos virais')).toBeTruthy();
});

test('unidades diferentes nao desenham a barra de perda, mas a etapa aparece', () => {
  const unidadesDiferentes = [
    { name: 'candidatos virais', value: 400, unit: 'contig' },
    { name: 'sequências virais retidas', value: 120, unit: 'sequência' },
  ];
  render(<AttritionFunnel stages={unidadesDiferentes} />);
  expect(screen.getByTestId('stage-sequências virais retidas')).toBeTruthy();
  expect(screen.queryByTestId('loss-sequências virais retidas')).toBeNull();
});

test('ausencia de unit e tratada como mesma unidade (compatibilidade)', () => {
  const semUnit = [
    { name: 'contigs', value: 1000 },
    { name: 'candidatos virais', value: 400 },
  ];
  render(<AttritionFunnel stages={semUnit} />);
  expect(screen.getByTestId('loss-candidatos virais')).toBeTruthy();
});

test('sem stages nao quebra (default [])', () => {
  render(<AttritionFunnel />);
  expect(screen.queryByRole('img')).toBeTruthy();
});
