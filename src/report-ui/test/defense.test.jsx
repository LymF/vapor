import { render, screen, fireEvent } from '@testing-library/react';
import { App } from '../src/App.jsx';

const dados = {
  run: { title: 'VAPOR', samples: ['S1', 'S2'] },
  overview: { kpis: [] },
  prokaryotic: {
    quality: [
      { genome: 'S1__binette_bin1', source: 'S1', completeness: 95, contamination: 1,
        css: 0.02, gunc_pass: true, is_representative: true, representative: 'S1__binette_bin1' },
      { genome: 'S2__binette_bin2', source: 'S2', completeness: 91, contamination: 2,
        css: 0.1, gunc_pass: true, is_representative: true, representative: 'S2__binette_bin2' },
    ],
    clusters: { n_mags: 2, n_clusters: 2, sizes: [] },
    taxonomy: [
      { genome: 'S2__binette_bin2', source: 'S2', Phylum: 'Bacillota', Class: 'Bacilli',
        Order: '', Family: '', Genus: '', count: 1, inherited: false,
        representative: 'S2__binette_bin2' },
      { genome: 'S1__binette_bin1', source: 'S1', Phylum: 'Pseudomonadota', Class: 'Gammaproteobacteria',
        Order: 'Enterobacterales', Family: 'Enterobacteriaceae', Genus: 'Escherichia',
        count: 1, inherited: false, representative: 'S1__binette_bin1' },
    ],
  },
  defense_amr: {
    defense: [
      { genome: 'S1__binette_bin1', system: 'CBASS', count: 1 },
      { genome: 'S1__binette_bin1', system: 'RM', count: 1 },
      { genome: 'S2__binette_bin2', system: 'Gabija', count: 1 },
    ],
    amr: [
      { genome: 'S1__binette_bin1', contig: 'k141_1', gene: 'tetA',
        drug_class: 'tetracycline', n_tools: 3, tools: 'amrfinder,rgi,deeparg' },
      { genome: 'S2__binette_bin2', contig: 'k141_1', gene: 'blaTEM',
        drug_class: 'beta-lactam', n_tools: 2, tools: 'amrfinder,rgi' },
    ],
    plasmids: [
      { genome: 'S1__binette_bin1', contig: 'k141_1', gene: 'IncFII', start: 100, end: 900 },
    ],
    islands: [
      { genome: 'S1__binette_bin1', contig: 'k141_1', start: 1000, end: 5800,
        n_genes: 5, n_systems: 3, systems: ['CBASS', 'RM', 'Wadjet'],
        genes: [
          { start: 1000, end: 1800, strand: '+', label: 'CBASS', kind: 'defense' },
          { start: 2000, end: 2800, strand: '-', label: '', kind: 'outro' },
          { start: 5000, end: 5800, strand: '+', label: 'Wadjet', kind: 'defense' },
        ] },
    ],
    colocalization: { n_args: 2, n_args_on_replicon: 1,
                      args_on_replicon: ['S1__binette_bin1|k141_1|tetA'] },
    upset: {
      sets: { replicon: 1, 'ARG de consenso': 2, 'sistema de defesa': 2, 'ilha de defesa': 1 },
      combos: [
        { tools: ['ARG de consenso', 'sistema de defesa'], count: 1 },
        { tools: ['ARG de consenso', 'ilha de defesa', 'replicon', 'sistema de defesa'], count: 1 },
      ],
    },
  },
};

function abreAba() {
  render(<App data={dados} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Defesa, AMR e plasmídeos' }));
}

test('a aba some sem a trilha de defesa/AMR', () => {
  const sem = { run: dados.run, overview: { kpis: [] }, prokaryotic: dados.prokaryotic };
  render(<App data={sem} />);
  expect(screen.queryByRole('tab', { name: 'Defesa, AMR e plasmídeos' })).toBeNull();
});

test('o heatmap de defesa ordena as linhas pela taxonomia GTDB', () => {
  abreAba();
  // Ordem alfabetica de MAG colocaria S1 antes de S2. Pela taxonomia,
  // Bacillota vem antes de Pseudomonadota -- e a ordem que torna legivel
  // "este clado carrega CBASS".
  const heat = screen.getByTestId('defense-heatmap');
  expect(heat.getAttribute('data-rows'))
    .toBe('S2__binette_bin2,S1__binette_bin1');
});

test('ARG traz o numero de ferramentas alem da cor', () => {
  abreAba();
  const barras = screen.getByTestId('amr-classes');
  // Cor sozinha nao identifica nada: o encoding secundario e o n de
  // ferramentas concordantes, e ele precisa estar legivel no DOM.
  expect(barras.getAttribute('data-tools-tetracycline')).toBe('3');
  expect(barras.getAttribute('data-tools-beta-lactam')).toBe('2');
});

test('a fracao de ARGs em contig com replicon vem com a ressalva escrita', () => {
  abreAba();
  expect(screen.getByText('1 de 2')).toBeTruthy();
  // A ressalva nao e opcional: replicon em contig de MAG e evidencia de
  // origem plasmidial, nao prova de plasmidio intacto.
  expect(screen.getByText(/não prova de plasmídeo intacto/i)).toBeTruthy();
});

test('a trilha da ilha desenha os genes de defesa distintos dos demais', () => {
  abreAba();
  expect(document.querySelector('[data-feature="CBASS"]')).toBeTruthy();
  const legenda = screen.getByTestId('island-track');
  expect(legenda.getAttribute('data-genome')).toBe('S1__binette_bin1');
});

test('o filtro de amostra recorta os MAGs mostrados', () => {
  abreAba();
  expect(screen.getByTestId('defense-heatmap').getAttribute('data-rows'))
    .toContain('S1__binette_bin1');
  fireEvent.change(screen.getByRole('combobox'), { target: { value: 'S2' } });
  expect(screen.getByTestId('defense-heatmap').getAttribute('data-rows'))
    .toBe('S2__binette_bin2');
});

test('sem plasmidio a secao de colocalizacao mostra lacuna, nunca 0%', () => {
  const sem = {
    ...dados,
    defense_amr: { ...dados.defense_amr, plasmids: undefined, colocalization: undefined },
  };
  render(<App data={sem} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Defesa, AMR e plasmídeos' }));
  expect(screen.queryByText('0 de 2')).toBeNull();
  expect(screen.getByText(/PlasmidFinder não rodou/i)).toBeTruthy();
});
