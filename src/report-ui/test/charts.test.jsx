import { render, screen } from '@testing-library/react';
import { StatTile } from '../src/charts/StatTile.jsx';
import { StatusMatrix } from '../src/charts/StatusMatrix.jsx';

test('StatTile mostra rotulo e valor', () => {
  render(<StatTile label="vOTUs" value={1234} sub="catalogo global" />);
  expect(screen.getByText('vOTUs')).toBeTruthy();
  expect(screen.getByText('1.234')).toBeTruthy();
});

test('StatusMatrix distingue falha de zero biologico', () => {
  render(<StatusMatrix rows={[
    { rule: 'defensefinder', sample: 'S1', status: 'failed', reason: 'disco cheio' },
    { rule: 'defensefinder', sample: 'S2', status: 'ok', reason: '' },
  ]} />);
  const falhou = screen.getByTestId('cell-defensefinder-S1');
  const passou = screen.getByTestId('cell-defensefinder-S2');
  expect(falhou.getAttribute('data-status')).toBe('failed');
  expect(passou.getAttribute('data-status')).toBe('ok');
  expect(falhou.getAttribute('aria-label')).toContain('disco cheio');
});

test('StatusMatrix nunca representa failed e skipped com o mesmo estado', () => {
  render(<StatusMatrix rows={[
    { rule: 'gunc', sample: 'S1', status: 'failed', reason: 'x' },
    { rule: 'gunc', sample: 'S2', status: 'skipped', reason: 'desligado' },
  ]} />);
  expect(screen.getByTestId('cell-gunc-S1').getAttribute('data-status'))
    .not.toBe(screen.getByTestId('cell-gunc-S2').getAttribute('data-status'));
});

test('unknown nao e apresentado como ok', () => {
  render(<StatusMatrix rows={[
    { rule: 'bakta', sample: 'S1', status: 'unknown', reason: 'no status recorded' },
    { rule: 'bakta', sample: 'S2', status: 'ok', reason: '' },
  ]} />);
  const desconhecido = screen.getByTestId('cell-bakta-S1');
  expect(desconhecido.getAttribute('data-status')).toBe('unknown');
  expect(desconhecido.getAttribute('aria-label')).toContain('sem status');
});
