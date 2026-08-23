// Aba Catalogo de MAGs (§5.4 do desenho). O eixo e o CATALOGO, nao a amostra:
// o galah desreplica uma vez para a rodada inteira e tudo a jusante do binning
// roda no representante. O filtro global de amostra continua funcionando, mas
// como recorte -- "quais destes MAGs foram binados aqui" -- e nao como o eixo
// da narrativa.
import { RankSelector, RANKS } from '../viz/RankSelector.jsx';
import { BarChart, StackedBar } from '../charts/BarChart.jsx';
import { Sunburst } from '../charts/Sunburst.jsx';
import { Heatmap } from '../charts/Heatmap.jsx';
import { StatTile } from '../charts/StatTile.jsx';
import { Scatter, ZONAS_MIMAG, zonaMIMAG, COR_ZONA } from '../charts/Scatter.jsx';
import { Legend } from '../viz/Legend.jsx';
import { useReport, TODAS } from '../state/store.jsx';
import { foldOther } from '../viz/palette.js';

const NAO_CLASSIFICADO = 'Unclassified';

// Classes CAZy na ordem em que a literatura as lista, nao por frequencia:
// ordem estavel entre rodadas e o que permite comparar dois reports.
const ORDEM_CAZY = ['GH', 'GT', 'PL', 'CE', 'AA', 'CBM', 'Other'];

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

export function MagCatalog() {
  const { data, rank, sample } = useReport();
  const prok = data?.prokaryotic;

  if (!prok) {
    return <p className="empty">Sem dados para esta aba nesta rodada.</p>;
  }

  const daAmostra = (r) => sample === TODAS || r.source === sample;
  const quality = (prok.quality ?? []).filter(daAmostra);
  const taxonomy = (prok.taxonomy ?? []).filter(daAmostra);

  const pontos = quality.map((r) => ({
    id: r.genome,
    x: r.completeness,
    y: r.contamination,
    zona: zonaMIMAG(r.completeness, r.contamination),
    row: r,
  }));

  // O quadrante que interessa: contaminacao baixa no CheckM2 nao exclui
  // quimera. Um MAG que passa num e reprova no outro precisa saltar aos
  // olhos, ou o usuario o trata como bom genoma.
  const quimeras = quality.filter((r) => r.gunc_pass === false);
  const idsQuimera = new Set(quimeras.map((r) => r.genome));

  const pontosGunc = quality
    .filter((r) => r.css !== null && r.css !== undefined)
    .map((r) => ({
      id: r.genome, x: r.contamination, y: r.css,
      zona: zonaMIMAG(r.completeness, r.contamination), row: r,
    }));

  const clusters = prok.clusters;
  const barrasCluster = (clusters?.sizes ?? [])
    .slice(0, 20)
    .map((c) => ({ name: c.representative, value: c.n_members }));

  const porFonte = {};
  (prok.quality ?? []).forEach((r) => {
    porFonte[r.source] = (porFonte[r.source] || 0) + 1;
  });
  const barrasFonte = Object.entries(porFonte)
    .map(([name, value]) => ({ name, value }));

  const ranksDisponiveis = RANKS.filter((r) => taxonomy.some((row) => temValor(row[r])));
  const barrasTaxon = foldOther(agregaPorRank(taxonomy, rank), 8)
    .map(([name, value]) => ({ name, value }));
  const herdados = taxonomy.filter((r) => r.inherited);

  const kegg = prok.kegg;
  const modulosKegg = kegg?.modules ?? [];
  const genomasKegg = (kegg?.genomes ?? []).filter(
    (g) => sample === TODAS || quality.some((r) => r.genome === g || r.representative === g),
  );
  const valoresKegg = {};
  genomasKegg.forEach((g) => { valoresKegg[g] = kegg?.values?.[g] ?? {}; });
  // O `missing_ko` viaja no proprio DOM porque e o que o tooltip le, e sem
  // ele um modulo a 66% e um numero sem interpretacao possivel.
  const attrsMissing = {};
  modulosKegg.forEach((m) => {
    const pares = Object.entries(m.missing ?? {})
      .filter(([g]) => genomasKegg.includes(g))
      .map(([g, ko]) => `${g}:${ko}`);
    if (pares.length) attrsMissing[`data-missing-${m.module.toLowerCase()}`] = pares.join(';');
  });

  const cazy = (prok.cazy ?? []).filter(
    (r) => sample === TODAS || quality.some((q) => q.genome === r.genome || q.representative === r.genome),
  );

  const nAlta = quality.filter(
    (r) => zonaMIMAG(r.completeness, r.contamination) === 'high-quality').length;

  return (
    <div className="panel">
      <div className="kpi-row">
        <StatTile label="MAGs" value={quality.length} />
        <StatTile label="Clusters do catálogo" value={clusters?.n_clusters ?? 0}
                  sub="95% ANI, galah" />
        <StatTile label="Alta qualidade (MIMAG)" value={nAlta} />
        <StatTile label="Reprovam no GUNC" value={quimeras.length} />
      </div>

      <section className="card">
        <h2>Qualidade dos MAGs</h2>
        <p className="chart__sub">
          CheckM2 e GUNC são medidos na própria amostra — a cobertura que separa
          os bins é dela. Nada aqui é herdado do representante.
        </p>
        <Scatter
          testid="mag-quality"
          dataAttrs={{ 'data-n': String(quality.length) }}
          points={pontos}
          zones={ZONAS_MIMAG}
          xName="Completude (%)"
          yName="Contaminação (%)"
          xDomain={[0, 100]}
          colorOf={(p) => COR_ZONA[p.zona]}
          flagAttr="data-gunc-fail"
          flagOf={(p) => (idsQuimera.has(p.id) ? p.id : null)}
          tooltipOf={(p) => `${p.id} — ${p.x}% completo, ${p.y}% contaminado`}
        />
        <Legend items={[
          { label: 'Alta qualidade', color: COR_ZONA['high-quality'] },
          { label: 'Qualidade média', color: COR_ZONA['medium-quality'] },
          { label: 'Baixa qualidade', color: COR_ZONA['low-quality'] },
          { label: 'Sem CheckM2', color: COR_ZONA.unknown },
        ]} />
      </section>

      <section className="card">
        <h2>Quimerismo (GUNC)</h2>
        {pontosGunc.length ? (
          <>
            <p className="notice">
              {quimeras.length
                ? `${quimeras.length} MAG(s) reprova(m) no GUNC — quimerismo não aparece na contaminação do CheckM2, então um genoma pode passar num e reprovar no outro.`
                : 'Nenhum MAG reprova no GUNC nesta rodada.'}
            </p>
            <Scatter
              testid="mag-gunc"
              points={pontosGunc}
              xName="Contaminação CheckM2 (%)"
              yName="Clade separation score (GUNC)"
              colorOf={(p) => COR_ZONA[p.zona]}
              flagAttr="data-gunc-fail"
              flagOf={(p) => (idsQuimera.has(p.id) ? p.id : null)}
              tooltipOf={(p) => `${p.id} — CSS ${p.y}, contaminação ${p.x}%`}
            />
          </>
        ) : <p className="empty">GUNC não rodou nesta rodada.</p>}
      </section>

      <section className="card">
        <h2>Estrutura do catálogo</h2>
        <p className="chart__sub">
          Toda anotação procariótica a jusante do binning é computada no
          representante do cluster e herdada pelos membros. O tamanho do
          cluster é, portanto, quantos MAGs dependem de uma única anotação.
        </p>
        {barrasCluster.length ? (
          <BarChart orientation="horizontal" sort="desc" valueName="MAGs no cluster"
                    data={barrasCluster} />
        ) : <p className="empty">Sem desreplicação nesta rodada.</p>}
        <h3>Proveniência</h3>
        <BarChart sort="desc" valueName="MAGs" data={barrasFonte} />
      </section>

      <section className="card">
        <h2>Taxonomia GTDB</h2>
        {taxonomy.length ? (
          <>
            <RankSelector availableRanks={ranksDisponiveis} />
            <BarChart orientation="horizontal" sort="desc" valueName="MAGs"
                      data={barrasTaxon} />
            <Sunburst rows={taxonomy} />
            {herdados.length ? (
              <p className="notice" data-inherited-from={herdados[0].representative}>
                {herdados.length} MAG(s) com taxonomia <strong>herdada</strong> do
                representante do cluster — o GTDB-Tk classificou o representante,
                não este genoma. Ex.: {herdados[0].genome} herda de {herdados[0].representative}.
              </p>
            ) : null}
          </>
        ) : <p className="empty">Sem classificação GTDB nesta rodada.</p>}
      </section>

      <section className="card">
        <h2>Metabolismo</h2>
        <h3>Completude de módulo KEGG</h3>
        {modulosKegg.length && genomasKegg.length ? (
          <div data-testid="mag-kegg" {...attrsMissing}>
            <p className="chart__sub">
              Calculado no representante do cluster (herdado pelos membros). O
              tooltip traz o <code>missing_ko</code>: uma via incompleta só é
              interpretável com o passo que falta.
            </p>
            <Heatmap
              rows={genomasKegg}
              cols={modulosKegg.map((m) => m.module)}
              values={valoresKegg}
            />
          </div>
        ) : <p className="empty">Sem completude de módulo KEGG nesta rodada.</p>}

        <h3>Famílias CAZy</h3>
        {cazy.length ? (
          <StackedBar order={ORDEM_CAZY}
                      data={cazy.map((r) => ({ name: r.genome, parts: r.parts }))} />
        ) : <p className="empty">Sem anotação CAZy nesta rodada.</p>}
      </section>
    </div>
  );
}
