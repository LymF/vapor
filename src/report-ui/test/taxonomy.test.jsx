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

test('Phylum e Class vazios nao apagam a arvore -- desenha a partir de Family', () => {
  // Caso comum de taxonomia viral: classificado so ate familia. O rank
  // global comeca em 'Family' (default do ReportProvider), entao a arvore
  // deve desenhar o arco de familia em vez de ficar vazia por causa dos dois
  // ranks iniciais sem dado nenhum.
  const linhas = [
    { Phylum: '', Class: '', Order: '', Family: 'Straboviridae', count: 7 },
    { Phylum: '', Class: '', Order: '', Family: 'Straboviridae', count: 3 },
  ];
  render(<Prov><Sunburst rows={linhas} /></Prov>);
  expect(screen.getByTestId('arc-Straboviridae')).toBeTruthy();
  expect(screen.queryByText('Sem dados nesta rodada')).toBeNull();
});

test('sem nenhuma taxonomia atribuida mostra estado vazio especifico', () => {
  const linhas = [{ count: 5 }];
  render(<Prov><Sunburst rows={linhas} /></Prov>);
  expect(screen.getByText('Nenhuma taxonomia atribuída nesta rodada')).toBeTruthy();
});

test('sem nenhuma linha mostra o estado vazio generico, nao o especifico', () => {
  render(<Prov><Sunburst rows={[]} /></Prov>);
  expect(screen.getByText('Sem dados nesta rodada')).toBeTruthy();
});

test('ranks sem dado ficam desabilitados no seletor, com title explicando', () => {
  render(<Prov><RankSelector availableRanks={['Family', 'Genus']} /></Prov>);
  const filo = screen.getByRole('button', { name: 'Filo' });
  const classe = screen.getByRole('button', { name: 'Classe' });
  const familia = screen.getByRole('button', { name: 'Família' });
  const genero = screen.getByRole('button', { name: 'Gênero' });
  expect(filo.disabled).toBe(true);
  expect(classe.disabled).toBe(true);
  expect(familia.disabled).toBe(false);
  expect(genero.disabled).toBe(false);
  expect(filo.title).toMatch(/não tem taxonomia atribuída/);
});
