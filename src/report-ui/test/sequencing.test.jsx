import { render, screen, fireEvent } from '@testing-library/react';
import { App } from '../src/App.jsx';

const base = { run: { title: 'VAPOR', samples: ['S1', 'S2'] }, overview: { kpis: [] } };
const comDados = {
  ...base,
  sequencing: {
    qc: [{ sample: 'S1', reads_before: 100, reads_after: 90, q30: 0.93 },
         { sample: 'S2', reads_before: 200, reads_after: 150, q30: 0.88 }],
    quast: { S1: { '# contigs': '1200', N50: '3400' }, S2: { '# contigs': '900', N50: '2900' } },
    mapping: [{ sample: 'S1', rate: 0.82 }, { sample: 'S2', rate: 0.61 }],
    lengths: { S1: [1000, 2000, 3000], S2: [1500, 2500] },
  },
};

test('a aba Sequenciamento aparece quando ha dado', () => {
  render(<App data={comDados} />);
  expect(screen.getByRole('tab', { name: 'Sequenciamento' })).toBeTruthy();
});

test('sem o bloco sequencing a aba nao e listada', () => {
  render(<App data={base} />);
  expect(screen.queryByRole('tab', { name: 'Sequenciamento' })).toBeNull();
});

test('o heatmap do QUAST normaliza por metrica, nao entre metricas', () => {
  render(<App data={comDados} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Sequenciamento' }));
  expect(screen.getByTestId('heatmap').getAttribute('data-normalize')).toBe('per-col');
});

test('o filtro de amostra restringe o painel', () => {
  render(<App data={comDados} />);
  fireEvent.click(screen.getByRole('tab', { name: 'Sequenciamento' }));
  fireEvent.change(screen.getByLabelText('Amostra'), { target: { value: 'S1' } });
  expect(screen.getByTestId('seq-scope').textContent).toBe('S1');
});
