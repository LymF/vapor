import { render, screen, fireEvent } from '@testing-library/react';
import { App } from '../src/App.jsx';
import { zonaMIMAG } from '../src/charts/Scatter.jsx';

const dados = {
  run: { title: 'VAPOR', samples: ['S1', 'S2'] },
  overview: { kpis: [] },
  prokaryotic: {
    quality: [
      { genome: 'S1__binette_bin1', source: 'S1', completeness: 95.2, contamination: 1.4,
        css: 0.02, gunc_pass: true, is_representative: true, representative: 'S1__binette_bin1' },
      { genome: 'S2__binette_bin1', source: 'S2', completeness: 72.0, contamination: 3.1,
        css: 0.11, gunc_pass: true, is_representative: false, representative: 'S1__binette_bin1' },
      { genome: 'S2__binette_bin2', source: 'S2', completeness: 91.0, contamination: 2.0,
        css: 0.87, gunc_pass: false, is_representative: true, representative: 'S2__binette_bin2' },
    ],
    clusters: {
      n_mags: 3,
      n_clusters: 2,
      sizes: [
        { representative: 'S1__binette_bin1', n_members: 2, n_sources: 2 },
        { representative: 'S2__binette_bin2', n_members: 1, n_sources: 1 },
      ],
    },
    taxonomy: [
      { genome: 'S1__binette_bin1', source: 'S1', Phylum: 'Pseudomonadota', Class: 'Gammaproteobacteria',
        Order: 'Enterobacterales', Family: 'Enterobacteriaceae', Genus: 'Escherichia',
        count: 1, inherited: false, representative: 'S1__binette_bin1' },
      { genome: 'S2__binette_bin1', source: 'S2', Phylum: 'Pseudomonadota', Class: 'Gammaproteobacteria',
        Order: 'Enterobacterales', Family: 'Enterobacteriaceae', Genus: 'Escherichia',
        count: 1, inherited: true, representative: 'S1__binette_bin1' },
      { genome: 'S2__binette_bin2', source: 'S2', Phylum: 'Bacillota', Class: 'Bacilli',
        Order: '', Family: '', Genus: '', count: 1, inherited: false,
        representative: 'S2__binette_bin2' },
    ],
    kegg: {
      genomes: ['S1__binette_bin1', 'S2__binette_bin2'],
      modules: [
        { module: 'M00001', name: 'Glycolysis', missing: { 'S2__binette_bin2': 'K01810' } },
        { module: 'M00002', name: 'Pentose phosphate', missing: { 'S1__binette_bin1': 'K00615' } },
      ],
      values: {
        'S1__binette_bin1': { M00001: 100.0, M00002: 66.7 },
        'S2__binette_bin2': { M00001: 80.0 },
      },
    },
    cazy: [
      { genome: 'S1__binette_bin1', parts: { GH: 2, GT: 1 } },
      { genome: 'S2__binette_bin2', parts: { PL: 1 } },
    ],
  },
};

function abreAba() {
  render(<App data={dados} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Catálogo de MAGs' }));
}

test('a aba some quando a rodada nao tem catalogo de MAGs', () => {
  const semProk = { run: dados.run, overview: { kpis: [] } };
  render(<App data={semProk} />);
  expect(screen.queryByRole('tab', { name: 'Catálogo de MAGs' })).toBeNull();
});

test('o scatter desenha as linhas de corte MIMAG', () => {
  abreAba();
  // As zonas sao a razao de o grafico existir: sem elas e um scatter
  // qualquer de dois numeros.
  expect(document.querySelector('[data-mimag="high-quality"]')).toBeTruthy();
  expect(document.querySelector('[data-mimag="medium-quality"]')).toBeTruthy();
});

test('zonaMIMAG segue os cortes de Bowers 2017, nao aproximacoes', () => {
  expect(zonaMIMAG(95, 1)).toBe('high-quality');
  // 90/5 sao os limites exatos: 90 completo entra, 5 de contaminacao NAO.
  expect(zonaMIMAG(90, 4.9)).toBe('high-quality');
  expect(zonaMIMAG(90, 5)).toBe('medium-quality');
  expect(zonaMIMAG(50, 9.9)).toBe('medium-quality');
  expect(zonaMIMAG(49.9, 1)).toBe('low-quality');
  expect(zonaMIMAG(null, 1)).toBe('unknown');
});

test('MAG que passa no CheckM2 e reprova no GUNC fica destacado', () => {
  abreAba();
  // O quadrante "boa completude, CSS alto" e o achado do painel: contaminacao
  // baixa no CheckM2 nao exclui quimera.
  const alerta = document.querySelector('[data-gunc-fail="S2__binette_bin2"]');
  expect(alerta).toBeTruthy();
  // O aviso escrito importa tanto quanto o glifo: "passou no CheckM2" nao
  // significa "nao e quimera", e essa frase precisa estar na aba.
  expect(screen.getByText(/não aparece na contaminação do CheckM2/i)).toBeTruthy();
});

test('a taxonomia GTDB respeita o seletor de rank global', () => {
  abreAba();
  fireEvent.click(screen.getByRole('button', { name: 'Gênero' }));
  expect(screen.getByText('Escherichia')).toBeTruthy();
  fireEvent.click(screen.getByRole('button', { name: 'Filo' }));
  expect(screen.queryByText('Escherichia')).toBeNull();
  expect(screen.getByText('Pseudomonadota')).toBeTruthy();
});

test('MAG com taxonomia herdada aparece marcado, com o representante', () => {
  abreAba();
  const marca = document.querySelector('[data-inherited-from="S1__binette_bin1"]');
  expect(marca).toBeTruthy();
  expect(marca.textContent).toMatch(/herdada/i);
  // O representante precisa estar nomeado: "herdado" sem dizer de quem nao
  // deixa o leitor conferir a afirmacao.
  expect(marca.textContent).toContain('S1__binette_bin1');
});

test('o filtro global de amostra repinta a estrutura do catalogo', () => {
  abreAba();
  expect(screen.getByTestId('mag-quality').getAttribute('data-n')).toBe('3');
  fireEvent.change(screen.getByRole('combobox'), { target: { value: 'S1' } });
  expect(screen.getByTestId('mag-quality').getAttribute('data-n')).toBe('1');
});

test('o heatmap KEGG carrega missing_ko por genoma', () => {
  abreAba();
  const heat = screen.getByTestId('mag-kegg');
  // Guardado no proprio DOM porque e o que o tooltip le: uma via a 66,7% sem
  // o passo que falta nao e interpretavel.
  expect(heat.getAttribute('data-missing-m00002')).toBe('S1__binette_bin1:K00615');
});

test('sem KEGG a secao mostra estado vazio, nao heatmap sem eixo', () => {
  const semKegg = { ...dados, prokaryotic: { ...dados.prokaryotic, kegg: undefined } };
  render(<App data={semKegg} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Catálogo de MAGs' }));
  expect(screen.queryByTestId('mag-kegg')).toBeNull();
  expect(screen.getByText(/sem completude de módulo/i)).toBeTruthy();
});
