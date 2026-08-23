// Aba Defesa, AMR e plasmideos (§5.5 do desenho).
//
// Tudo aqui foi computado no REPRESENTANTE do cluster do catalogo: o
// DefenseFinder roda uma vez por genoma representante, e as ferramentas de
// AMR uma vez sobre o .faa concatenado do catalogo. Um MAG membro nao tem
// resultado proprio -- ele herda o do representante. A aba diz isso em vez de
// deixar o leitor supor que cada MAG foi analisado.
import { useState } from 'react';
import { BarChart } from '../charts/BarChart.jsx';
import { Heatmap } from '../charts/Heatmap.jsx';
import { UpSet } from '../charts/UpSet.jsx';
import { StatTile } from '../charts/StatTile.jsx';
import { GenomeTrack } from '../charts/GenomeTrack.jsx';
import { useReport, TODAS } from '../state/store.jsx';
import { PAL, PAL_MUTED } from '../viz/palette.js';
import { RANKS } from '../viz/RankSelector.jsx';

const TIPOS_ILHA = {
  defense: { color: PAL[5], label: 'gene de sistema de defesa' },
  outro: { color: PAL_MUTED, label: 'outro gene na janela' },
};

// Ordem taxonomica: a linhagem inteira concatenada, do filo ao genero. Ordenar
// por nome de MAG e ordem de arquivo, nao de biologia -- e "este clado carrega
// CBASS" so se le quando clados vizinhos ficam vizinhos.
function chaveTaxonomica(taxonomy) {
  const porGenoma = {};
  (taxonomy ?? []).forEach((r) => {
    porGenoma[r.genome] = RANKS.map((k) => r[k] || '~').join('|');
  });
  // '~' para rank vazio empurra o nao-classificado para o fim, em vez de
  // ordena-lo junto do que comeca com letra.
  return (g) => porGenoma[g] ?? '~~';
}

export function DefenseAmr() {
  const { data, sample } = useReport();
  const bloco = data?.defense_amr;
  const prok = data?.prokaryotic;
  const [ilhaAtiva, setIlhaAtiva] = useState(0);

  if (!bloco) {
    return <p className="empty">Sem dados para esta aba nesta rodada.</p>;
  }

  // O genoma aqui e sempre um representante; a amostra vem da tabela de
  // qualidade, que e a unica que sabe de onde cada MAG foi binado.
  const fonteDoGenoma = {};
  (prok?.quality ?? []).forEach((r) => { fonteDoGenoma[r.genome] = r.source; });
  const daAmostra = (g) => sample === TODAS || fonteDoGenoma[g] === sample;

  const defense = (bloco.defense ?? []).filter((r) => daAmostra(r.genome));
  const amr = (bloco.amr ?? []).filter((r) => daAmostra(r.genome));
  const plasmids = (bloco.plasmids ?? []).filter((r) => daAmostra(r.genome));
  const islands = (bloco.islands ?? []).filter((r) => daAmostra(r.genome));

  const chave = chaveTaxonomica(prok?.taxonomy);
  const genomasDefesa = Array.from(new Set(defense.map((r) => r.genome)))
    .sort((a, b) => chave(a).localeCompare(chave(b)) || a.localeCompare(b));
  const sistemas = Array.from(new Set(defense.map((r) => r.system))).sort();
  const valoresDefesa = {};
  defense.forEach((r) => {
    valoresDefesa[r.genome] = valoresDefesa[r.genome] || {};
    valoresDefesa[r.genome][r.system] = r.count;
  });

  const porClasse = {};
  amr.forEach((r) => {
    const classe = r.drug_class || 'não classificada';
    const d = porClasse[classe] || { name: classe, value: 0, maxTools: 0 };
    d.value += 1;
    d.maxTools = Math.max(d.maxTools, r.n_tools || 0);
    porClasse[classe] = d;
  });
  const barrasAmr = Object.values(porClasse);
  const attrsTools = {};
  barrasAmr.forEach((d) => {
    attrsTools[`data-tools-${d.name.toLowerCase()}`] = String(d.maxTools);
  });

  const coloc = bloco.colocalization;
  const ilha = islands[Math.min(ilhaAtiva, islands.length - 1)];

  return (
    <div className="panel">
      <p className="card__scope">
        Computado nos <strong>representantes</strong> do catálogo de MAGs e
        herdado pelos membros do mesmo cluster (95% ANI). Um MAG membro não
        tem resultado próprio.
      </p>

      <div className="kpi-row">
        <StatTile label="MAGs com sistema de defesa" value={genomasDefesa.length} />
        <StatTile label="Tipos de sistema" value={sistemas.length} />
        <StatTile label="ARGs de consenso" value={amr.length}
                  sub="≥ 2 ferramentas concordantes" />
        <StatTile label="Replicons plasmidiais" value={plasmids.length} />
      </div>

      <section className="card">
        <h2>Sistemas de defesa por MAG</h2>
        {genomasDefesa.length ? (
          <div data-testid="defense-heatmap" data-rows={genomasDefesa.join(',')}>
            <p className="chart__sub">
              Linhas ordenadas pela taxonomia GTDB, não pelo nome do MAG.
            </p>
            <Heatmap rows={genomasDefesa} cols={sistemas} values={valoresDefesa}
                     sparseAsBubble />
          </div>
        ) : <p className="empty">Nenhum sistema de defesa detectado nesta seleção.</p>}
      </section>

      <section className="card">
        <h2>Genes de resistência</h2>
        {barrasAmr.length ? (
          <div data-testid="amr-classes" {...attrsTools}>
            <p className="chart__sub">
              Só o consenso de duas ou mais ferramentas entra aqui. O tooltip
              traz quantas concordaram — o número de ferramentas é encoding
              próprio, não uma nuance de cor.
            </p>
            <BarChart orientation="horizontal" sort="desc" valueName="ARGs"
                      data={barrasAmr} />
          </div>
        ) : <p className="empty">Nenhum ARG de consenso nesta seleção.</p>}
      </section>

      <section className="card">
        <h2>ARG em contig com replicon plasmidial</h2>
        {coloc ? (
          <>
            <div className="kpi-row">
              <StatTile label="ARGs em contig com replicon"
                        value={`${coloc.n_args_on_replicon} de ${coloc.n_args}`} />
            </div>
            <p className="notice">
              Replicon num contig de MAG é evidência de origem plasmidial,
              <strong> não prova de plasmídeo intacto</strong>: a montagem pode
              ter quebrado o elemento. O cruzamento exige o mesmo genoma e o
              mesmo contig — nomes de contig se repetem entre MAGs.
            </p>
            <ul className="plain-list">
              {coloc.args_on_replicon.map((a) => <li key={a}><code>{a}</code></li>)}
            </ul>
          </>
        ) : <p className="empty">PlasmidFinder não rodou nesta rodada — sem replicon não há pergunta de colocalização.</p>}
      </section>

      <section className="card">
        <h2>Coocorrência por MAG</h2>
        {bloco.upset ? (
          <UpSet sets={bloco.upset.sets} combos={bloco.upset.combos} valueName="MAGs" />
        ) : <p className="empty">Sem evidências para cruzar nesta rodada.</p>}
      </section>

      <section className="card">
        <h2>Ilhas de defesa</h2>
        {ilha ? (
          <div data-testid="island-track" data-genome={ilha.genome}>
            <p className="chart__sub">
              {ilha.genome} · {ilha.contig} · {ilha.n_genes} genes,{' '}
              {ilha.n_systems} sistemas ({ilha.systems.join(', ')})
            </p>
            {islands.length > 1 ? (
              <div className="rank-selector" role="group" aria-label="Ilha">
                {islands.map((i, idx) => (
                  <button key={`${i.genome}-${i.contig}-${i.start}`} type="button"
                          className="rank-selector__btn"
                          aria-pressed={idx === ilhaAtiva}
                          data-active={idx === ilhaAtiva ? 'true' : 'false'}
                          onClick={() => setIlhaAtiva(idx)}>
                    {i.contig}
                  </button>
                ))}
              </div>
            ) : null}
            <GenomeTrack length={ilha.end} features={ilha.genes} kinds={TIPOS_ILHA} />
          </div>
        ) : <p className="empty">Nenhuma ilha de defesa nesta seleção.</p>}
      </section>
    </div>
  );
}
