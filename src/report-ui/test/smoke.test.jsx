import { render, screen } from '@testing-library/react';
import { App } from '../src/App.jsx';

test('monta e mostra o titulo do report', () => {
  render(<App data={{ run: { title: 'VAPOR', samples: [] } }} />);
  expect(screen.getByText('VAPOR')).toBeTruthy();
});
