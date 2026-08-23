import { render, screen, fireEvent } from '@testing-library/react';
import { App } from '../src/App.jsx';

const dados = {
  run: { title: 'VAPOR', samples: ['S1', 'S2'] },
  overview: { kpis: [] },
  diversity: {
    alpha: [
      { sample: 'S1', domain: 'viral', index: 'observed', value: 120 },
      { sample: 'S2', domain: 'viral', index: 'observed', value: 95 },
      { sample: 'S1', domain: 'viral', index: 'shannon', value: 3.41 },
      { sample: 'S2', domain: 'viral', index: 'shannon', value: 3.02 },
      { sample: 'S1', domain: 'prok', index: 'shannon', value: 4.8 },
    ],
    alpha_missing: ['simpson', 'chao1'],
    pcoa: {
      viral: [
        { sample: 'S1', pc1: 0.31, pc2: -0.12, var_pc1: 0.425, var_pc2: 0.183 },
        { sample: 'S2', pc1: -0.29, pc2: 0.15, var_pc1: 0.425, var_pc2: 0.183 },
      ],
    },
    procrustes: {
      disparity: 0.214,
      pairs: [
        { sample: 'S1', viral_pc1: 0.30, viral_pc2: -0.10, prok_pc1: 0.28, prok_pc2: -0.14 },
        { sample: 'S2', viral_pc1: -0.28, viral_pc2: 0.13, prok_pc1: -0.26, prok_pc2: 0.16 },
      ],
    },
  },
  reads: {
    id_space_warning: 'Os identificadores desta aba são genomas de referência do banco do sylph (t__IMGVR_UViG_…) e não têm relação com os contigs montados (k141_…).',
    viral: [{ clade: 'r__Duplodnaviria|p__Uroviricota', _eff_rank: 'p__', S1: 12.5, S2: 8.0 }],
    prok: [{ clade: 'd__Bacteria|p__Pseudomonadota', _eff_rank: 'p__', S1: 40.0, S2: 35.5 }],
    host: [
      { host_genus: 'Escherichia', host_source: 'db', n_viral_taxa: 4, S1: 9.0, S2: 5.0 },
      { host_genus: 'desconhecido', host_source: 'none', n_viral_taxa: 2, S1: 3.5, S2: 3.0 },
    ],
    truncated: { prok: 340 },
  },
};

test('Simpson e Chao1 ausentes viram lacuna nomeada, nunca barra em zero', () => {
  render(<App data={dados} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Diversidade' }));
  const aviso = screen.getByTestId('alpha-missing');
  expect(aviso.textContent).toMatch(/simpson/i);
  expect(aviso.textContent).toMatch(/chao1/i);
  // A razao precisa estar escrita: sao estimadores de CONTAGEM, e esta
  // rodada nao tinha contagens de reads.
  expect(aviso.textContent).toMatch(/contagens/i);
});

test('o PCoA nomeia a variancia explicada nos dois eixos', () => {
  render(<App data={dados} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Diversidade' }));
  // Sem isso, a distancia entre dois pontos do PCoA nao tem interpretacao.
  expect(screen.getByText(/PC1 \(42[.,]5%\)/)).toBeTruthy();
  expect(screen.getByText(/PC2 \(18[.,]3%\)/)).toBeTruthy();
});

test('o Procrustes desenha uma seta por amostra, ligando os dois pontos', () => {
  render(<App data={dados} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Diversidade' }));
  expect(document.querySelectorAll('[data-procrustes-pair]').length).toBe(2);
  expect(screen.getByText(/0[.,]214/)).toBeTruthy();
});

test('sem alfa a aba de diversidade nao aparece', () => {
  render(<App data={{ run: dados.run, overview: { kpis: [] } }} />);
  expect(screen.queryByRole('tab', { name: 'Diversidade' })).toBeNull();
});

test('a aba de leituras avisa que seus IDs nao casam com os contigs', () => {
  render(<App data={dados} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Leituras' }));
  expect(screen.getByText(/não têm relação com os contigs montados/i)).toBeTruthy();
});

test('o hospedeiro mostra a procedencia da anotacao', () => {
  render(<App data={dados} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Leituras' }));
  const linha = screen.getByTestId('host-desconhecido');
  // 'none' e 'db' precisam ficar distinguiveis: hospedeiro desconhecido nao
  // e o mesmo que hospedeiro anotado.
  expect(linha.getAttribute('data-source')).toBe('none');
  expect(screen.getByTestId('host-Escherichia').getAttribute('data-source')).toBe('db');
});

test('o corte por abundancia e declarado, nao silencioso', () => {
  render(<App data={dados} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Leituras' }));
  expect(screen.getByText(/340 .*fora/i)).toBeTruthy();
});
