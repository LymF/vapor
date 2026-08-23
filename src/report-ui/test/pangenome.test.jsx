import { render, screen, fireEvent } from '@testing-library/react';
import { App } from '../src/App.jsx';

const dados = {
  run: { title: 'VAPOR', samples: ['S1', 'S2'] },
  overview: { kpis: [] },
  pangenome: {
    candidates: [
      { representative: 'S1__bin1', n_members: 3, n_islands: 1, n_systems: 4,
        n_args: 2, n_plasmid: 1, criterio: 'ilha de defesa', eligible: true },
      { representative: 'S2__bin9', n_members: 4, n_islands: 0, n_systems: 0,
        n_args: 0, n_plasmid: 3, criterio: 'sem evidencia de defesa/amr',
        eligible: false },
    ],
    clusters: [
      { representative: 'S1__bin1', n_members: 3, n_evaluable: 2, n_core: 1,
        n_variable: 2, n_singleton: 2, completeness_median: 84.8,
        size_median: 4200000,
        taxonomy: 'd__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria' },
    ],
    matrix: {
      S1__bin1: {
        members: ['S1__bin1', 'S1__bin7', 'S1__bin8'],
        rows: [
          { tipo: 'defesa', gene: 'CBASS', freq: '2/2', n_present: 2, n_evaluable: 2,
            states: { S1__bin1: 'x', S1__bin7: 'x', S1__bin8: '?' } },
          { tipo: 'amr', gene: 'tetA', freq: '1/2', n_present: 1, n_evaluable: 2,
            states: { S1__bin1: '.', S1__bin7: 'x', S1__bin8: '?' } },
        ],
      },
    },
    completeness: { S1__bin1: 98.0, S1__bin7: 71.5, S1__bin8: 42.0 },
  },
};

function abreAba() {
  render(<App data={dados} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Pangenoma' }));
}

test('a aba some sem a trilha de pangenoma', () => {
  render(<App data={{ run: dados.run, overview: { kpis: [] } }} />);
  expect(screen.queryByRole('tab', { name: 'Pangenoma' })).toBeNull();
});

test('o nao avaliavel e hachurado, nunca a cor da ausencia', () => {
  abreAba();
  const naoAvaliavel = document.querySelector('[data-cell="defesa|CBASS|S1__bin8"]');
  const ausente = document.querySelector('[data-cell="amr|tetA|S1__bin1"]');
  expect(naoAvaliavel.getAttribute('data-state')).toBe('?');
  // Textura, nao matiz: um tom mais claro da mesma cor de ausencia leria
  // como "meio ausente", que e precisamente a confusao a evitar.
  expect(naoAvaliavel.getAttribute('fill')).toContain('hatch');
  expect(ausente.getAttribute('fill')).not.toContain('hatch');
});

test('o denominador excluindo o nao avaliavel vai escrito na aba', () => {
  abreAba();
  expect(screen.getByText(/excluídos do denominador/i)).toBeTruthy();
  // A frequencia mostrada e a da regra (2/2), com tres membros no cluster --
  // o leitor precisa ver os dois numeros para nao ler 2/2 como "todos".
  expect(screen.getByText(/2\/2/)).toBeTruthy();
});

test('a matriz desenha so os membros do cluster', () => {
  abreAba();
  expect(screen.getByTestId('state-matrix').getAttribute('data-cols'))
    .toBe('S1__bin1,S1__bin7,S1__bin8');
});

test('defesa e amr com o mesmo gene ficam em linhas distintas', () => {
  const comColisao = {
    ...dados,
    pangenome: {
      ...dados.pangenome,
      matrix: {
        S1__bin1: {
          members: ['S1__bin1', 'S1__bin7'],
          rows: [
            { tipo: 'defesa', gene: 'tetA', freq: '1/2', n_present: 1, n_evaluable: 2,
              states: { S1__bin1: 'x', S1__bin7: '.' } },
            { tipo: 'amr', gene: 'tetA', freq: '1/2', n_present: 1, n_evaluable: 2,
              states: { S1__bin1: '.', S1__bin7: 'x' } },
          ],
        },
      },
    },
  };
  render(<App data={comColisao} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Pangenoma' }));
  expect(document.querySelector('[data-cell="defesa|tetA|S1__bin1"]')
    .getAttribute('data-state')).toBe('x');
  expect(document.querySelector('[data-cell="amr|tetA|S1__bin1"]')
    .getAttribute('data-state')).toBe('.');
});

test('a completude do membro fica visivel junto da matriz', () => {
  abreAba();
  // 42% e a razao de S1__bin8 ser '?'. Sem o numero a vista, o '?' vira
  // um simbolo arbitrario.
  expect(screen.getByText(/S1__bin8.*42/)).toBeTruthy();
});

test('plasmidio sozinho aparece como sinal de mobilidade, nao como criterio', () => {
  abreAba();
  const linha = screen.getByTestId('candidate-S2__bin9');
  expect(linha.textContent).toContain('sem evidencia de defesa/amr');
  expect(linha.getAttribute('data-eligible')).toBe('false');
  expect(screen.getByText(/sinal de mobilidade/i)).toBeTruthy();
});

test('core e variaveis aparecem por cluster', () => {
  abreAba();
  const barras = screen.getByTestId('pangenome-core');
  expect(barras.getAttribute('data-core')).toBe('1');
  expect(barras.getAttribute('data-variable')).toBe('2');
});
