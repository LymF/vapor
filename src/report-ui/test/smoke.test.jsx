import { render, screen } from '@testing-library/react';
import { App } from '../src/index.jsx';

test('monta e mostra o titulo do report', () => {
  render(<App data={{ run: { title: 'VAPOR' } }} />);
  expect(screen.getByText('VAPOR')).toBeTruthy();
});
