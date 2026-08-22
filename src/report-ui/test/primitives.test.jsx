import { render, screen } from '@testing-library/react';
import { scaleLinear } from 'd3';
import { Chart } from '../src/viz/Chart.jsx';
import { AxisBottom, AxisLeft } from '../src/viz/Axis.jsx';
import { Legend } from '../src/viz/Legend.jsx';
import { Tooltip, useTooltip } from '../src/viz/Tooltip.jsx';

test('Chart mostra titulo e passa dimensoes ao filho', () => {
  render(<Chart title="Contigs" sub="por amostra">{({ width, height }) =>
    <text data-testid="dim">{`${width}x${height}`}</text>}</Chart>);
  expect(screen.getByText('Contigs')).toBeTruthy();
  expect(screen.getByTestId('dim').textContent).toMatch(/^\d+x\d+$/);
});

test('Chart em estado vazio nao chama o filho', () => {
  let chamou = false;
  render(<Chart title="X" empty emptyLabel="Sem dados nesta rodada">{() => { chamou = true; return null; }}</Chart>);
  expect(chamou).toBe(false);
  expect(screen.getByText('Sem dados nesta rodada')).toBeTruthy();
});

test('Legend carrega rotulo alem da cor', () => {
  render(<Legend items={[{ label: 'Bacteria', color: '#0d9488' }]} />);
  expect(screen.getByText('Bacteria')).toBeTruthy();
});

test('Chart com container sem dimensao medida ainda assim renderiza sem lancar e passa largura > 0', () => {
  // jsdom nao dispara ResizeObserver: useResize fica em {width: 0}. Chart
  // precisa de um fallback de largura sensato em vez de desenhar SVG 0px
  // ou deixar o filho criar escalas com dominio vazio.
  let larguraRecebida = null;
  expect(() => {
    render(<Chart title="Y">{({ width }) => {
      larguraRecebida = width;
      return <text data-testid="w">{width}</text>;
    }}</Chart>);
  }).not.toThrow();
  expect(larguraRecebida).toBeGreaterThan(0);
});

test('AxisBottom e AxisLeft desenham sem lancar com escala valida', () => {
  const scale = scaleLinear().domain([0, 100]).range([0, 200]);
  expect(() => {
    render(
      <svg>
        <AxisBottom scale={scale} width={200} height={100} />
        <AxisLeft scale={scale} width={200} height={100} />
      </svg>,
    );
  }).not.toThrow();
});

test('useTooltip fornece show/hide/node e node comeca oculto', () => {
  function Sonda() {
    const { show, hide, node } = useTooltip();
    return (
      <div>
        <button type="button" onClick={(e) => show(e, 'conteudo')}>mostrar</button>
        <button type="button" onClick={hide}>esconder</button>
        {node}
      </div>
    );
  }
  render(<Sonda />);
  // sem interacao, nao deve haver tooltip visivel no DOM
  expect(screen.queryByText('conteudo')).toBeNull();
});

test('Tooltip standalone nao lanca sem conteudo', () => {
  expect(() => render(<Tooltip content={null} x={0} y={0} />)).not.toThrow();
});
