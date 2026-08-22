import { render, screen, fireEvent } from '@testing-library/react';
import { App } from '../src/App.jsx';

const data = {
  run: { title: 'VAPOR', samples: ['S1', 'S2'] },
  overview: {
    kpis: [{ label: 'Amostras', value: 2 }],
    status: [{ rule: 'megahit', sample: 'S1', status: 'ok', reason: '' }],
    funnel: {
      __all__: { stages: [{ name: 'contigs', value: 10 }], losses: {} },
      S1: { stages: [{ name: 'contigs', value: 4 }], losses: {} },
    },
  },
};

test('abre na Visao geral', () => {
  render(<App data={data} />);
  expect(screen.getByRole('tab', { selected: true }).textContent).toBe('Visão geral');
});

test('o filtro de amostra lista todas as amostras mais a opcao agregada', () => {
  render(<App data={data} />);
  const seletor = screen.getByLabelText('Amostra');
  expect([...seletor.options].map((o) => o.value)).toEqual(['__all__', 'S1', 'S2']);
});

test('trocar de amostra troca os dados do painel', () => {
  render(<App data={data} />);
  expect(screen.getByTestId('stage-contigs').getAttribute('width')).not.toBe('0');
  fireEvent.change(screen.getByLabelText('Amostra'), { target: { value: 'S1' } });
  expect(screen.getByTestId('funnel-scope').textContent).toBe('S1');
});

test('aba sem dados mostra estado vazio, nao eixo quebrado', () => {
  render(<App data={{ run: { title: 'VAPOR', samples: [] }, overview: {} }} />);
  expect(screen.getByText(/sem dados/i)).toBeTruthy();
});
