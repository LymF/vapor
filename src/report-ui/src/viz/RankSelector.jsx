// As chaves sao as colunas do dado (ingles, como saem do MMseqs2/GTDB-Tk) --
// nunca traduzidas. RANK_LABEL faz a ponte so na exibicao (REPORT_VIZ_GUIDE.md).
import { useReport } from '../state/store.jsx';

export const RANKS = ['Phylum', 'Class', 'Order', 'Family', 'Genus'];

export const RANK_LABEL = {
  Phylum: 'Filo',
  Class: 'Classe',
  Order: 'Ordem',
  Family: 'Família',
  Genus: 'Gênero',
};

// Controle segmentado: muda o rank global. Nenhum grafico de taxonomia nasce
// preso a um rank -- o rank e estado do ReportProvider, nao prop.
export function RankSelector() {
  const { rank, setRank } = useReport();
  return (
    <div className="rank-selector" role="group" aria-label="Nível taxonômico">
      {RANKS.map((r) => (
        <button
          key={r}
          type="button"
          className="rank-selector__btn"
          aria-pressed={r === rank}
          data-active={r === rank ? 'true' : 'false'}
          onClick={() => setRank(r)}
        >
          {RANK_LABEL[r]}
        </button>
      ))}
    </div>
  );
}
