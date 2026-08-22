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
//
// `availableRanks` (opcional, default = todos) vem de quem chama e ja tem as
// linhas de taxonomia em maos (ViralCatalog.jsx) -- evita duplicar aqui a
// logica de "quais linhas entram nesta aba" que o painel ja resolve (filtro
// de amostra etc). Um rank fora dessa lista fica desabilitado com `title`
// explicando o motivo: sem isso o usuario clica em "Filo", ve o sunburst
// vazio e conclui que a pipeline falhou, quando na verdade e so que a
// taxonomia viral desta rodada nao chega la.
export function RankSelector({ availableRanks = RANKS }) {
  const { rank, setRank } = useReport();
  return (
    <div className="rank-selector" role="group" aria-label="Nível taxonômico">
      {RANKS.map((r) => {
        const disponivel = availableRanks.includes(r);
        return (
          <button
            key={r}
            type="button"
            className="rank-selector__btn"
            aria-pressed={r === rank}
            data-active={r === rank ? 'true' : 'false'}
            disabled={!disponivel}
            title={disponivel ? undefined : 'Esta rodada não tem taxonomia atribuída neste nível'}
            onClick={() => setRank(r)}
          >
            {RANK_LABEL[r]}
          </button>
        );
      })}
    </div>
  );
}
