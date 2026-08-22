import { ReportProvider, useReport, TODAS } from './state/store.jsx';
import { Overview } from './panels/Overview.jsx';
import { Sequencing } from './panels/Sequencing.jsx';
import { ViralCatalog } from './panels/ViralCatalog.jsx';
import { ThemeToggle } from './viz/theme.js';

// As abas do desenho. Os paineis dos planos 3+ entram aqui; enquanto nao
// existem, a aba nao e listada -- nunca uma aba vazia sem explicacao.
// `temDado` decide isso por bloco do payload, nao por presenca do painel.
const ABAS = [
  { id: 'overview', label: 'Visão geral', Painel: Overview, temDado: () => true },
  { id: 'sequencing', label: 'Sequenciamento', Painel: Sequencing, temDado: (data) => Boolean(data?.sequencing) },
  { id: 'viral', label: 'Catálogo viral', Painel: ViralCatalog, temDado: (data) => Boolean(data?.viral) },
];

function Shell() {
  const { data, sample, setSample, samples, tab, setTab } = useReport();
  const abasVisiveis = ABAS.filter((a) => a.temDado(data));
  const ativa = abasVisiveis.find((a) => a.id === tab) ?? abasVisiveis[0];
  const { Painel } = ativa;

  return (
    <div className="app">
      <nav className="nav">
        <span className="nav__brand">{data?.run?.title ?? 'VAPOR'}</span>
        <div role="tablist" className="nav__tabs">
          {abasVisiveis.map((a) => (
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
