import { render, screen, fireEvent } from '@testing-library/react';
import { App } from '../src/App.jsx';

const dados = {
  run: { title: 'VAPOR', samples: ['S1', 'S2'], groups: [] },
  overview: {
    kpis: [],
    status: [
      { rule: 'megahit', sample: 'S1', status: 'ok', reason: '' },
      { rule: 'gunc', sample: 'S1', status: 'failed', reason: 'disco cheio', kind: 'rule' },
      { rule: 'gtdbtk', sample: 'S1', status: 'skipped', reason: 'sem MAG', kind: 'view' },
      { rule: 'megahit', sample: 'S2', status: 'ok', reason: '' },
    ],
    funnel: {
      S1: { stages: [{ name: 'contigs', value: 900, unit: 'contig' },
                     { name: 'sequências virais retidas', value: 12, unit: 'sequência' }],
            losses: { 'sequências virais retidas': [{ reason: 'Low-quality', count: 40 }] } },
      S2: { stages: [{ name: 'contigs', value: 700, unit: 'contig' }], losses: {} },
    },
  },
  sequencing: {
    qc: [{ sample: 'S1', reads_before: 1000, reads_after: 950, q30: 93.2 },
         { sample: 'S2', reads_before: 800, reads_after: 770, q30: 91.0 }],
    mapping: [{ sample: 'S1', rate: 88.5 }, { sample: 'S2', rate: 74.0 }],
    lengths: { S1: { values: [1000, 2000, 5000], n: 3 },
               S2: { values: [900, 1500], n: 2 } },
  },
  prokaryotic: {
    quality: [
      { genome: 'S1__binette_bin1', source: 'S1', completeness: 95, contamination: 1,
        css: 0.02, gunc_pass: true, is_representative: true, representative: 'S1__binette_bin1' },
      { genome: 'S2__binette_bin1', source: 'S2', completeness: 72, contamination: 3,
        css: 0.11, gunc_pass: true, is_representative: false, representative: 'S1__binette_bin1' },
    ],
    clusters: { n_mags: 2, n_clusters: 1, sizes: [{ representative: 'S1__binette_bin1', n_members: 2, n_sources: 2 }] },
    taxonomy: [
      { genome: 'S1__binette_bin1', source: 'S1', Phylum: 'Pseudomonadota', Class: 'Gammaproteobacteria',
        Order: 'Enterobacterales', Family: 'Enterobacteriaceae', Genus: 'Escherichia',
        count: 1, inherited: false, representative: 'S1__binette_bin1' },
      { genome: 'S2__binette_bin1', source: 'S2', Phylum: 'Pseudomonadota', Class: 'Gammaproteobacteria',
        Order: 'Enterobacterales', Family: 'Enterobacteriaceae', Genus: 'Escherichia',
        count: 1, inherited: true, representative: 'S1__binette_bin1' },
    ],
  },
  diversity: {
    alpha: [
      { sample: 'S1', domain: 'viral', index: 'shannon', value: 3.41 },
      { sample: 'S2', domain: 'viral', index: 'shannon', value: 3.02 },
    ],
    alpha_missing: [],
  },
};

function abre(amostra) {
  render(<App data={dados} />);
  if (amostra) {
    fireEvent.change(screen.getByRole('combobox'), { target: { value: amostra } });
  }
  fireEvent.click(screen.getByRole('tab', { name: 'Por amostra' }));
}

test('sem amostra escolhida a vista diz qual assumiu, em vez de ficar vazia', () => {
  abre(null);
  expect(screen.getByTestId('sample-view').getAttribute('data-sample')).toBe('S1');
  expect(screen.getByText(/nenhuma amostra selecionada/i)).toBeTruthy();
});

test('a vista individual mostra os blocos da amostra escolhida', () => {
  abre('S2');
  const vista = screen.getByTestId('sample-view');
  expect(vista.getAttribute('data-sample')).toBe('S2');
  // Execucao: a regra que falhou nesta amostra tem de aparecer como falha,
  // nunca como ausencia.
  expect(screen.getByTestId('sample-status').textContent).toContain('megahit');
});

test('regra que falhou aparece como falha, com o motivo', () => {
  abre('S1');
  const status = screen.getByTestId('sample-status');
  expect(status.textContent).toContain('gunc');
  expect(status.textContent).toContain('disco cheio');
});

test('o MAG binado aqui traz seu representante e o valor herdado marcado', () => {
  abre('S2');
  const linha = screen.getByTestId('mag-S2__binette_bin1');
  // Sob (h), a taxonomia deste MAG saiu do representante -- que veio de OUTRA
  // amostra. Sem a marca, "o MAG da S2 e Escherichia" seria uma afirmacao que
  // a pipeline nunca fez.
  expect(linha.querySelector('[data-inherited-from="S1__binette_bin1"]')).toBeTruthy();
});

test('o que foi medido na amostra e marcado como medido, nao deixado implicito', () => {
  abre('S1');
  const qc = screen.getByTestId('sample-sequencing');
  expect(qc.querySelector('[data-measured="true"]')).toBeTruthy();
  // O MAG de S1 E o representante: nada a herdar.
  expect(screen.getByTestId('mag-S1__binette_bin1')
    .querySelector('[data-inherited-from]')).toBeNull();
});

test('o alfa da amostra aparece contra a distribuicao das demais', () => {
  abre('S1');
  const alfa = screen.getByTestId('sample-alpha');
  // Um valor de alfa isolado nao significa nada: a posicao relativa e o que
  // lhe da sentido.
  expect(alfa.getAttribute('data-value')).toBe('3.41');
  expect(alfa.getAttribute('data-n-outras')).toBe('1');
});

test('no modo comparacao nao se desenha um sunburst por amostra', () => {
  abre('S1');
  fireEvent.click(screen.getByRole('button', { name: 'Comparar' }));
  fireEvent.click(screen.getByRole('button', { name: 'S2' }));
  // Hierarquia nao se compara lado a lado (§5.3): a comparacao usa pequenos
  // multiplos alinhados e barra empilhada, nunca N sunbursts.
  expect(document.querySelectorAll('[data-testid="sunburst"]').length).toBe(0);
  expect(screen.getByTestId('compare-view').getAttribute('data-samples')).toBe('S1,S2');
});

test('o modo comparacao alinha as distribuicoes no mesmo eixo', () => {
  abre('S1');
  fireEvent.click(screen.getByRole('button', { name: 'Comparar' }));
  fireEvent.click(screen.getByRole('button', { name: 'S2' }));
  // Um DistPlot com as duas amostras como grupos: mesmo eixo x, que e o que
  // torna a comparacao legivel.
  expect(screen.getByTestId('compare-lengths')).toBeTruthy();
});

test('uma vista que falhou diz que o calculo e do representante', () => {
  abre('S1');
  const status = screen.getByTestId('sample-status');
  // 'gtdbtk pulada' numa amostra nao significa que a amostra ficou sem
  // classificacao: o classify_wf roda uma vez, sobre as representantes.
  expect(status.textContent).toMatch(/vista.*representante do cluster/i);
});
