import { render, screen, fireEvent } from '@testing-library/react';
import { App } from '../src/App.jsx';

const dados = {
  run: { title: 'VAPOR', samples: ['S1'] },
  overview: { kpis: [] },
  viral: {
    taxonomy: [
      { sample: 'S1', Phylum: 'Uroviricota', Class: 'Caudoviricetes', Family: 'Straboviridae', count: 12 },
      { sample: 'S1', Phylum: 'Uroviricota', Class: 'Caudoviricetes', Family: 'Drexlerviridae', count: 7 },
      { sample: 'S1', Phylum: '', Class: '', Family: '', count: 40 },
    ],
    checkv_tiers: { S1: { 'Complete': 2, 'High-quality': 5, 'Low-quality': 30, '': 11 } },
    detectors: { sets: { geNomad: 40, VirSorter2: 30 },
                 combos: [{ tools: ['geNomad'], count: 15 },
                          { tools: ['geNomad', 'VirSorter2'], count: 25 }] },
    catalog: { n_votus: 800, n_pool: 2000, reduction_pct: 60 },
    explorer: [{ votu_id: 'S1|k141_1', length: 30000, features: [
      { start: 100, end: 900, strand: '+', label: 'terminase', kind: 'phrog' }] }],
  },
};

test('a composicao taxonomica respeita o rank escolhido', () => {
  render(<App data={dados} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Catálogo viral' }));
  fireEvent.click(screen.getByRole('button', { name: 'Família' }));
  expect(screen.getByText('Straboviridae')).toBeTruthy();
  fireEvent.click(screen.getByRole('button', { name: 'Filo' }));
  expect(screen.queryByText('Straboviridae')).toBeNull();
  expect(screen.getByText('Uroviricota')).toBeTruthy();
});

test('tier vazio do CheckV nao vira Complete nem some', () => {
  render(<App data={dados} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Catálogo viral' }));
  const barra = screen.getByTestId('checkv-tiers');
  expect(barra.getAttribute('data-series')).toContain('Not-determined');
  expect(barra.getAttribute('data-complete')).toBe('2');
});

test('linhagem vazia vira Unclassified, nao um taxon de nome vazio', () => {
  render(<App data={dados} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Catálogo viral' }));
  expect(screen.getByText(/Unclassified/)).toBeTruthy();
});

test('o explorador desenha a trilha do vOTU escolhido', () => {
  render(<App data={dados} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Catálogo viral' }));
  expect(document.querySelector('[data-feature="terminase"]')).toBeTruthy();
});

test('sem anotacao o explorador mostra estado vazio, nao trilha vazia', () => {
  const semAnot = { ...dados, viral: { ...dados.viral, explorer: [] } };
  render(<App data={semAnot} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Catálogo viral' }));
  expect(screen.getByText(/sem anotação/i)).toBeTruthy();
});
