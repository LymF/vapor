import { ReportProvider, useReport, TODAS } from './state/store.jsx';
import { Overview } from './panels/Overview.jsx';
import { ThemeToggle } from './viz/theme.js';

// As abas do desenho. Os paineis dos planos 2 e 3 entram aqui; enquanto nao
// existem, a aba nao e listada -- nunca uma aba vazia sem explicacao.
const ABAS = [
  { id: 'overview', label: 'Visão geral', Painel: Overview },
];

function Shell() {
  const { data, sample, setSample, samples, tab, setTab } = useReport();
  const ativa = ABAS.find((a) => a.id === tab) ?? ABAS[0];
  const { Painel } = ativa;

  return (
    <div className="app">
      <nav className="nav">
        <span className="nav__brand">{data?.run?.title ?? 'VAPOR'}</span>
        <div role="tablist" className="nav__tabs">
          {ABAS.map((a) => (
            <button key={a.id} role="tab" aria-selected={a.id === ativa.id}
                    className={`nav__tab${a.id === ativa.id ? ' is-active' : ''}`}
                    onClick={() => setTab(a.id)}>
              {a.label}
            </button>
          ))}
        </div>
        <label className="nav__filter">
          <span>Amostra</span>
          <select value={sample} onChange={(e) => setSample(e.target.value)}>
            <option value={TODAS}>todas</option>
            {samples.map((s) => <option key={s} value={s}>{s}</option>)}
          </select>
        </label>
        <ThemeToggle />
      </nav>
      <main className="main"><Painel /></main>
    </div>
  );
}

export function App({ data }) {
  return (
    <ReportProvider data={data}>
      <Shell />
    </ReportProvider>
  );
}
