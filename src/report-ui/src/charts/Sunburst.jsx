// Sunburst taxonomico. Nenhum treemap -- decisao explicita do usuario
// registrada no spec. O rank e estado global (useReport): o componente
// nunca fica preso a um nivel fixo, e clicar num arco desce o rank global
// e fixa o filtro de taxon (a navegacao E o grafico, nao um enfeite).
import { hierarchy, partition, arc as d3arc } from 'd3';
import { Chart } from '../viz/Chart.jsx';
import { useTooltip } from '../viz/Tooltip.jsx';
import { PAL, PAL_MUTED } from '../viz/palette.js';
import { useReport } from '../state/store.jsx';
import { RANKS } from '../viz/RankSelector.jsx';

// Linhagem incompleta e normal em virus: a maioria dos contigs virais nao
// chega a genero. Um rank mais fundo que o dado nao pode inventar nivel nem
// gerar arco de nome vazio -- mas quando ALGUMA linha do grupo tem dado
// naquele rank e outra nao, a que nao tem vira "Unclassified" de verdade,
// nao some.
const NAO_CLASSIFICADO = 'Unclassified';

function temValor(v) {
  return v !== undefined && v !== null && v !== '';
}

// Um rank sem NENHUM dado na linhagem toda (ex.: Phylum/Class vazios porque
// a fonte so classificou ate familia) nao pode derrubar a arvore inteira --
// pula esse rank e tenta o proximo, em vez de abortar em depth=0. So aborta
// (devolve []) quando NENHUM rank restante tem dado -- essa lista esgotada e
// o unico caso de "sem taxonomia".
function buildTree(rows, ranks, depth = 0) {
  if (depth >= ranks.length) return [];
  const rankKey = ranks[depth];
  const algumTemDado = rows.some((r) => temValor(r[rankKey]));
  if (!algumTemDado) return buildTree(rows, ranks, depth + 1);
  const grupos = new Map();
  rows.forEach((r) => {
    const nome = temValor(r[rankKey]) ? r[rankKey] : NAO_CLASSIFICADO;
    if (!grupos.has(nome)) grupos.set(nome, []);
    grupos.get(nome).push(r);
  });
  return Array.from(grupos.entries()).map(([name, grupo]) => ({
    name,
    rank: rankKey,
    count: grupo.reduce((soma, r) => soma + (r.count || 0), 0),
    children: buildTree(grupo, ranks, depth + 1),
  }));
}

export function Sunburst({ rows = [], onDrill }) {
  const { rank, setRank, setTaxonFilter } = useReport();
  const { show, hide, node } = useTooltip();
  const ateIndex = Math.max(RANKS.indexOf(rank), 0);
  const ranksAtivos = RANKS.slice(0, ateIndex + 1);
  const arvore = buildTree(rows, ranksAtivos);
  // Duas causas bem diferentes de "vazio": nenhuma linha chegou (fonte
  // ausente/filtrada) vs. linhas existem mas nenhuma tem taxonomia atribuida
  // em rank algum ate o selecionado -- essa segunda precisa dizer isso, ou o
  // usuario le como "a pipeline nao rodou" em vez de "nao classificou".
  const semLinhas = rows.length === 0;
  const semTaxonomiaAtribuida = !semLinhas && arvore.length === 0;
  const vazio = semLinhas || semTaxonomiaAtribuida;

  function handleClick(d) {
    const idxAtual = RANKS.indexOf(d.data.rank);
    const proximo = RANKS[idxAtual + 1] ?? d.data.rank;
    setTaxonFilter({ rank: d.data.rank, name: d.data.name });
    setRank(proximo);
    onDrill?.(d.data.rank, d.data.name);
  }

  return (
    <div data-testid="sunburst" className="chart-wrap">
      <Chart
        height={360}
        empty={vazio}
        emptyLabel={semTaxonomiaAtribuida ? 'Nenhuma taxonomia atribuída nesta rodada' : undefined}
      >
        {({ width, height }) => {
          const raio = Math.min(width, height) / 2;
          const root = hierarchy({ name: 'root', children: arvore }, (d) => d.children)
            .sum((d) => (d.children && d.children.length ? 0 : d.count))
            .sort((a, b) => b.value - a.value);
          partition().size([2 * Math.PI, raio])(root);
          const arcGen = d3arc()
            .startAngle((d) => d.x0)
            .endAngle((d) => d.x1)
            .innerRadius((d) => d.y0)
            .outerRadius((d) => d.y1);
          const nos = root.descendants().filter((d) => d.depth > 0);
          return (
            <g transform={`translate(${width / 2},${height / 2})`}>
              {nos.map((d) => {
                const cor = d.data.name === NAO_CLASSIFICADO ? PAL_MUTED : PAL[(d.depth - 1) % PAL.length];
                return (
                  <path
                    key={`${d.depth}-${d.data.rank}-${d.data.name}`}
                    data-testid={`arc-${d.data.name}`}
                    d={arcGen(d)}
                    fill={cor}
                    stroke="var(--surface, #fff)"
                    strokeWidth={1}
                    onClick={() => handleClick(d)}
                    onMouseMove={(e) => show(e, `${d.data.rank} ${d.data.name}: ${d.value.toLocaleString('pt-BR')}`)}
                    onMouseLeave={hide}
                    style={{ cursor: 'pointer' }}
                  />
                );
              })}
            </g>
          );
        }}
      </Chart>
      {node}
    </div>
  );
}
