import { render, screen, fireEvent } from '@testing-library/react';
import { useTheme, ThemeToggle } from '../src/viz/theme.js';

function Sonda() {
  const { theme } = useTheme();
  return <span data-testid="t">{theme}</span>;
}

beforeEach(() => { document.documentElement.removeAttribute('data-theme'); });

test('comeca em light quando nada foi escolhido', () => {
  render(<Sonda />);
  expect(screen.getByTestId('t').textContent).toBe('light');
});

test('o toggle troca o tema e marca a raiz', () => {
  render(<><Sonda /><ThemeToggle /></>);
  fireEvent.click(screen.getByRole('button', { name: /tema/i }));
  expect(document.documentElement.getAttribute('data-theme')).toBe('dark');
});

test('nao quebra quando localStorage lanca', () => {
  const orig = Object.getOwnPropertyDescriptor(window, 'localStorage');
  Object.defineProperty(window, 'localStorage', {
    configurable: true,
    get() { throw new Error('bloqueado'); },
  });
  expect(() => render(<Sonda />)).not.toThrow();
  if (orig) Object.defineProperty(window, 'localStorage', orig);
});
