// Aba Catalogo viral (task 9 do plano 2). Seis blocos, cada um com sua
// propria fonte e seu proprio estado vazio -- nenhum grafico taxonomico fica
// preso a um rank fixo (RankSelector e Sunburst usam o rank global).
import { RankSelector, RANKS } from '../viz/RankSelector.jsx';
import { BarChart, StackedBar } from '../charts/BarChart.jsx';
import { Sunburst } from '../charts/Sunburst.jsx';
import { Heatmap } from '../charts/Heatmap.jsx';
import { UpSet } from '../charts/UpSet.jsx';
import { StatTile } from '../charts/StatTile.jsx';
import { VotuExplorer } from './VotuExplorer.jsx';
import { useReport, TODAS } from '../state/store.jsx';
import { foldOther, PAL_MUTED } from '../viz/palette.js';

// Linhagem incompleta e o caso comum em virus (item 4 do brief): a maioria
// dos contigs nao chega a genero. A ausencia vira uma categoria de verdade,
// nao um taxon de nome vazio nem um contig que desaparece do denominador.
const NAO_CLASSIFICADO = 'Unclassified';

function temValor(v) {
  return v !== undefined && v !== null && v !== '';
}

function agregaPorRank(rows, rankKey) {
  const totais = {};
  rows.forEach((r) => {
    const nome = temValor(r[rankKey]) ? r[rankKey] : NAO_CLASSIFICADO;
    totais[nome] = (totais[nome] || 0) + (r.count || 0);
  });
  return totais;
}

// Escada ordinal Complete -> ND: rampa boa->ruim de verdade (nao oito
// matizes categoricos -- item 5 do brief). 'Not-determined' e a traducao da
// chave vazia do CheckV, que significa "nunca avaliado", distinta de um tier
// baixo avaliado de fato.
const ORDEM_CHECKV = ['Complete', 'High-quality', 'Medium-quality', 'Low-quality', 'Not-determined'];
const CORES_CHECKV = {
  Complete: '#0d9488',
  'High-quality': '#5eead4',
  'Medium-quality': '#fbbf24',
  'Low-quality': '#f97316',
  'Not-determined': PAL_MUTED,
};

function normalizaTiers(tiers) {
  const saida = {};
  Object.entries(tiers ?? {}).forEach(([k, v]) => {
    const chave = k === '' ? 'Not-determined' : k;
    saida[chave] = (saida[chave] || 0) + v;
  });
  return saida;
}

export function ViralCatalog() {
  const { data, rank, sample } = useReport();
  const viral = data?.viral;

  if (!viral) {
    return <p className="empty">Sem dados para esta aba nesta rodada.</p>;
  }

  const taxonomyRows = viral.taxonomy ?? [];
  // So os ranks com pelo menos um valor real nesta rodada ficam clicaveis no
  // seletor -- ver RankSelector.jsx.
  const availableRanks = RANKS.filter((r) => taxonomyRows.some((row) => temValor(row[r])));
  const totaisRank = agregaPorRank(taxonomyRows, rank);
  const barrasTaxon = foldOther(totaisRank, 8).map(([name, value]) => ({ name, value }));

  const checkvTiers = viral.checkv_tiers ?? {};
  const amostrasCheckv = Object.keys(checkvTiers).filter((s) => sample === TODAS || s === sample);
  const linhasCheckv = amostrasCheckv.map((s) => ({ name: s, parts: normalizaTiers(checkvTiers[s]) }));
  const totalCheckv = {};
  linhasCheckv.forEach((l) => {
    Object.entries(l.parts).forEach(([k, v]) => { totalCheckv[k] = (totalCheckv[k] || 0) + v; });
  });
  const categoriasCheckv = ORDEM_CHECKV.filter((c) => totalCheckv[c] !== undefined);

  const presence = viral.presence;
  const votusPresenca = presence?.votus ?? [];
  const amostrasPresenca = data?.run?.samples ?? [];
  const valoresPresenca = {};
  votusPresenca.forEach((v) => {
    valoresPresenca[v.votu_id] = {};
    amostrasPresenca.forEach((s) => {
      const estado = v.samples?.[s];
      if (estado && estado !== 'absent') valoresPresenca[v.votu_id][s] = 1;
    });
  });

  const detectores = viral.detectors;
  const catalogo = viral.catalog;
  const lifestyle = viral.lifestyle;

  return (
    <div className="panel">
      {catalogo ? (
        <div className="kpi-row">
          <StatTile label="vOTUs no catálogo" value={catalogo.n_votus} />
          <StatTile label="Sequências no pool" value={catalogo.n_pool} />
          <StatTile label="Redução do pool" value={`${catalogo.reduction_pct}%`} />
        </div>
      ) : null}

      <section className="card">
        <h2>Composição taxonômica</h2>
        {taxonomyRows.length ? (
          <>
            <RankSelector availableRanks={availableRanks} />
            <BarChart orientation="horizontal" sort="desc" valueName="sequências" data={barrasTaxon} />
            <Sunburst rows={taxonomyRows} />
          </>
        ) : <p className="empty">Sem taxonomia viral nesta rodada.</p>}
      </section>

      <section className="card">
        <h2>Qualidade CheckV</h2>
        {linhasCheckv.length ? (
          <div
            data-testid="checkv-tiers"
            data-series={categoriasCheckv.join(',')}
            data-complete={String(totalCheckv.Complete ?? 0)}
          >
            <StackedBar normalize order={ORDEM_CHECKV} colors={CORES_CHECKV} data={linhasCheckv} />
          </div>
        ) : <p className="empty">Sem avaliação CheckV nesta rodada.</p>}
      </section>

      <section className="card">
        <h2>vOTU × amostra</h2>
        {votusPresenca.length ? (
          <Heatmap
            sparseAsBubble
            rows={votusPresenca.map((v) => v.votu_id)}
            cols={amostrasPresenca}
            values={valoresPresenca}
          />
        ) : <p className="empty">Sem matriz de presença/ausência nesta rodada.</p>}
      </section>

      <section className="card">
        <h2>Concordância de detectores</h2>
        {detectores ? (
          <UpSet sets={detectores.sets} combos={detectores.combos} valueName="contigs" />
        ) : <p className="empty">Sem consenso de detecção nesta rodada.</p>}
      </section>

      <section className="card">
        <h2>Estilo de vida</h2>
        {lifestyle?.rows?.length ? (
          <>
            <p className="notice">
              BACPHLIP só é confiável em genoma completo — resultados de vOTUs parciais são indicativos, não conclusivos.
            </p>
            <StackedBar
              normalize
              order={['lytic', 'lysogenic']}
              colors={{ lytic: '#0d9488', lysogenic: '#7c3aed' }}
              data={[{ name: 'catálogo', parts: lifestyle.counts }]}
            />
          </>
        ) : <p className="empty">Sem predição de estilo de vida nesta rodada.</p>}
      </section>

      <section className="card">
        <h2>Explorador de vOTU</h2>
        <VotuExplorer vOTUs={viral.explorer ?? []} />
      </section>
    </div>
  );
}
