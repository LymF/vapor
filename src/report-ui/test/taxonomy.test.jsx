import { render, screen, fireEvent } from '@testing-library/react';
import { ReportProvider, useReport } from '../src/state/store.jsx';
import { RankSelector, RANKS, RANK_LABEL } from '../src/viz/RankSelector.jsx';
import { Sunburst } from '../src/charts/Sunburst.jsx';

function Prov({ children }) {
  return <ReportProvider data={{ run: { title: 'VAPOR', samples: [] } }}>{children}</ReportProvider>;
}

function Espia() {
  const { rank, taxonFilter } = useReport();
  return (<>
    <span data-testid="rank">{rank}</span>
    <span data-testid="taxon">{taxonFilter?.name ?? ''}</span>
  </>);
}

test('os rotulos sao portugues sobre as chaves em ingles', () => {
  expect(RANKS).toEqual(['Phylum', 'Class', 'Order', 'Family', 'Genus']);
  expect(RANK_LABEL.Family).toBe('Família');
  expect(RANK_LABEL.Phylum).toBe('Filo');
});

test('o seletor troca o rank global', () => {
  render(<Prov><RankSelector /><Espia /></Prov>);
  fireEvent.click(screen.getByRole('button', { name: 'Família' }));
  expect(screen.getByTestId('rank').textContent).toBe('Family');
});

test('clicar num arco desce o rank e fixa o taxon', () => {
  const linhas = [{ Phylum: 'Uroviricota', Class: 'Caudoviricetes', count: 5 }];
  render(<Prov><Sunburst rows={linhas} /><Espia /></Prov>);
  fireEvent.click(screen.getByTestId('arc-Uroviricota'));
  expect(screen.getByTestId('rank').textContent).toBe('Class');
  expect(screen.getByTestId('taxon').textContent).toBe('Uroviricota');
});

test('rank mais fundo que o dado nao inventa nivel', () => {
  const linhas = [{ Phylum: 'Uroviricota', count: 2 }];
  render(<Prov><Sunburst rows={linhas} /></Prov>);
  expect(screen.queryByTestId('arc-undefined')).toBeNull();
});
