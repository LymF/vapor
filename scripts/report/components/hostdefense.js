/* hostdefense.js — Host Prediction (PHIST), Defense/Anti-Defense systems
   (DefenseFinder, with built-in AntiDefenseFinder), AMR (AMRFinderPlus+RGI
   curated vs. DeepARG exploratory), and the Host<->Virus<->Defense<->AMR
   cross-link matrix. */
(function () {
  'use strict';

  window.renderHostDefense = function () {
    const samples = typeof SAMPLES !== 'undefined' ? SAMPLES : [];
    _renderHostCollapseChart(samples);
    _renderDefenseBar(samples);
    _renderAmrBar(samples);
    _renderDefensePhylum();
    _renderDefenseDensity();
    _renderDefenseMechanism(samples);
    _renderDefenseCorrelation();
    _renderDefenseCooccurrence();
    _renderAmrClasses();
    _renderAmrOrigin();
    _renderCrossAnalysis();
    _renderViralAntidefense(samples);
    _renderDefenseIslands();
    makeSampleDropdown('sample-sel-hostpred', _renderHostPrediction, { allSamples: true });
    makeSampleDropdown('sample-sel-hostdefense-matrix', _renderMatrix, { allSamples: true });
  };

  // ── Host Collapse: viral RPKM by predicted host genus ────────────────────
  function _renderHostCollapseChart(samples) {
    const hc = typeof HOST_COLLAPSE !== 'undefined' ? HOST_COLLAPSE : {};

    // Collect top-N genera across all samples (by max total_rpkm in any sample)
    const genusMax = {};
    samples.forEach(s => {
      (hc[s] || []).forEach(d => {
        genusMax[d.genus] = Math.max(genusMax[d.genus] || 0, d.total_rpkm);
      });
    });
    const topGenera = Object.entries(genusMax)
      .sort((a, b) => b[1] - a[1])
      .slice(0, 12)
      .map(e => e[0]);

    if (!topGenera.length) {
      const el = document.getElementById('host-collapse-chart');
      if (el) el.innerHTML = '<p style="color:var(--text-muted);padding:1rem">No host predictions available.</p>';
      return;
    }

    const series = topGenera.map((genus, i) => ({
      name: genus, type: 'bar', stack: 'hc',
      color: PAL[i % PAL.length],
      data: samples.map(s => {
        const row = (hc[s] || []).find(d => d.genus === genus);
        return row ? +row.total_rpkm.toFixed(2) : 0;
      }),
    }));

    mkChart('host-collapse-chart', {
      title: { text: 'Viral Abundance by Predicted Host Genus (RPKM)' },
      tooltip: { trigger: 'axis', axisPointer: { type: 'shadow' } },
      legend: { data: topGenera, type: 'scroll', bottom: 0 },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: 'RPKM' },
      series,
      grid: { bottom: 90 },
    });
  }

  function _byCount(rows, key, samples) {
    const counts = {};
    samples.forEach(s => { counts[s] = 0; });
    (rows || []).forEach(r => { if (r.sample in counts) counts[r.sample] += 1; });
    return samples.map(s => counts[s]);
  }

  // Pearson correlation coefficient; returns null if either series has no variance.
  function _pearson(xs, ys) {
    const n = xs.length;
    if (n < 3) return null;
    const mx = xs.reduce((a, b) => a + b, 0) / n;
    const my = ys.reduce((a, b) => a + b, 0) / n;
    let sxy = 0, sxx = 0, syy = 0;
    for (let i = 0; i < n; i++) {
      const dx = xs[i] - mx, dy = ys[i] - my;
      sxy += dx * dy; sxx += dx * dx; syy += dy * dy;
    }
    if (sxx === 0 || syy === 0) return null;
    return sxy / Math.sqrt(sxx * syy);
  }

  // Abramowitz-Stegun 7.1.26 erf approximation (max error 1.5e-7) -- used
  // to get a standard normal CDF without a stats library (this report is a
  // single offline HTML file, no bundling/imports).
  function _erf(x) {
    const sign = x < 0 ? -1 : 1; x = Math.abs(x);
    const a1=0.254829592, a2=-0.284496736, a3=1.421413741, a4=-1.453152027, a5=1.061405429, p=0.3275911;
    const t = 1 / (1 + p * x);
    const y = 1 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*Math.exp(-x*x);
    return sign * y;
  }
  function _normalCDF(z) { return 0.5 * (1 + _erf(z / Math.SQRT2)); }

  // Two-sided p-value for a Pearson r via Fisher z-transform + normal
  // approximation -- the standard large-sample test (exact for n -> large;
  // close enough here since we're correlating across dozens-hundreds of
  // bins, not a handful).
  function _pearsonPValue(r, n) {
    if (n < 4 || Math.abs(r) >= 1) return n < 4 ? 1 : 0;
    const z = 0.5 * Math.log((1 + r) / (1 - r));
    const se = 1 / Math.sqrt(n - 3);
    return Math.max(0, Math.min(1, 2 * (1 - _normalCDF(Math.abs(z / se)))));
  }

  // Benjamini-Hochberg FDR correction (step-up procedure).
  function _fdrBH(pvals) {
    const m = pvals.length;
    const order = pvals.map((_, i) => i).sort((a, b) => pvals[a] - pvals[b]);
    const q = new Array(m);
    let prevQ = 1;
    for (let rank = m; rank >= 1; rank--) {
      const i = order[rank - 1];
      prevQ = Math.min(prevQ, pvals[i] * m / rank);
      q[i] = prevQ;
    }
    return q;
  }

  // One row per (sample, Bin) with defense/anti-defense/AMR counts + GTDB-Tk
  // taxonomy — built client-side from BIN_ANNOTATIONS + GTDB_DATA so no new
  // Snakemake/Python plumbing is needed for the cross-analysis charts below.
  function _buildBinMatrix() {
    const bins = typeof BIN_ANNOTATIONS !== 'undefined' ? BIN_ANNOTATIONS : {};
    const gtdb = typeof GTDB_DATA !== 'undefined' ? GTDB_DATA : [];
    const gtdbByKey = new Map();
    gtdb.forEach(r => gtdbByKey.set(`${r.sample}::${r.Bin}`, r));
    // CheckM2's own Genome_Size column (bp) -- keyed the same way
    // merge_prok_taxonomy (data_loaders.py) keys it, .fa-suffix stripped.
    const checkm2 = typeof CHECKM2 !== 'undefined' ? CHECKM2 : {};
    const sizeByKey = new Map();
    Object.entries(checkm2).forEach(([sample, rows]) => (rows || []).forEach(r => {
      const name = (r.Name || r.name || '').replace(/\.fa$/, '');
      if (!name) return;
      const size = parseFloat(r.Genome_Size || r.genome_size || 0);
      if (size > 0) sizeByKey.set(`${sample}::${name}`, size);
    }));
    return Object.entries(bins).map(([key, ann]) => {
      const sep = key.indexOf('::');
      const sample = key.slice(0, sep), Bin = key.slice(sep + 2);
      const tax = gtdbByKey.get(key) || {};
      const n_defense = (ann.Defense_systems || []).length;
      const genomeSizeMb = (sizeByKey.get(key) || 0) / 1e6;
      return {
        sample, Bin,
        Domain: tax.Domain || 'Unclassified',
        Phylum: tax.Phylum || 'Unclassified',
        n_defense,
        n_antidefense:    (ann.Antidefense_systems || []).length,
        n_amr:            (ann.AMR_genes || []).length,
        genome_size_mb:        genomeSizeMb,
        defense_density_per_mb: genomeSizeMb > 0 ? n_defense / genomeSizeMb : null,
      };
    });
  }

  // ── Host Prediction (PHIST) ──────────────────────────────────────────────
  function _renderHostPrediction(sample) {
    const all = typeof PHIST_DATA !== 'undefined' ? PHIST_DATA : [];
    const isAll = sample === '__all__';
    const rows = isAll ? all : all.filter(r => r.sample === sample);

    makeTable('phist-table', rows, [
      { key: 'sample', label: 'Sample' },
      { key: 'Virus',  label: 'Virus' },
      { key: 'Host',   label: 'Predicted Host (MAG)' },
      { key: 'Score',  label: '# Common k-mers' },
      { key: 'P_value',label: 'Adj. p-value' },
    ], { searchId: 'phist-search' });

    _renderHostNetwork(rows, isAll ? 'All samples' : sample);
  }

  function _renderHostNetwork(rows, label) {
    const el = document.getElementById('hostpred-network-svg');
    if (!el) return;
    el.innerHTML = '';
    if (!rows.length) {
      el.innerHTML = '<p style="color:var(--text-muted);padding:1rem">No PHIST host predictions for this selection.</p>';
      return;
    }

    const dark  = document.documentElement.dataset.theme === 'dark';
    const W     = el.clientWidth || 900;
    const H     = 520;

    const nodeMap = new Map();
    const edges = [];
    rows.forEach(r => {
      const vId = `v::${r.Virus}`, hId = `h::${r.Host}`;
      if (!nodeMap.has(vId)) nodeMap.set(vId, { id: vId, label: r.Virus, type: 'virus' });
      if (!nodeMap.has(hId)) nodeMap.set(hId, { id: hId, label: r.Host, type: 'host' });
      edges.push({ source: vId, target: hId, score: +r.Score || 0 });
    });
    const nodes = Array.from(nodeMap.values());
    const typeColor = { virus: '#7c3aed', host: '#0d9488' };

    const svg = d3.select(el).append('svg').attr('width', W).attr('height', H);
    const zoomLayer = svg.append('g').attr('class', 'zoom-layer');

    const linkSel = zoomLayer.append('g').selectAll('line').data(edges).join('line')
      .attr('stroke', dark ? '#334155' : '#cbd5e1')
      .attr('stroke-width', 1)
      .attr('stroke-opacity', 0.6);

    const nodeSel = zoomLayer.append('g').selectAll('circle').data(nodes).join('circle')
      .attr('r', d => d.type === 'host' ? 9 : 6)
      .attr('fill', d => typeColor[d.type])
      .attr('fill-opacity', 0.85)
      .attr('stroke', dark ? '#0f172a' : '#fff')
      .attr('stroke-width', 1.5)
      .style('cursor', 'pointer');

    const labels = zoomLayer.append('g').selectAll('text').data(nodes).join('text')
      .attr('font-size', d => d.type === 'host' ? 11 : 9)
      .attr('font-weight', d => d.type === 'host' ? 600 : 400)
      .attr('fill', dark ? '#e2e8f0' : '#334155')
      .attr('dy', '-.7em')
      .style('pointer-events', 'none')
      .text(d => d.label);

    const tip = d3.select(el).append('div').attr('class', 'd3-tooltip').style('display', 'none');
    nodeSel.on('mouseover', (e, d) => {
      tip.style('display', 'block').style('left', (e.offsetX + 12) + 'px').style('top', (e.offsetY - 12) + 'px')
        .html(`<strong>${d.type === 'host' ? 'Host MAG' : 'Virus'}:</strong> ${d.label}`);
    }).on('mouseout', () => tip.style('display', 'none'));

    const sim = d3.forceSimulation(nodes)
      .force('link',    d3.forceLink(edges).id(d => d.id).distance(45).strength(0.5))
      .force('charge',  d3.forceManyBody().strength(-70))
      .force('center',  d3.forceCenter(W / 2, H / 2))
      .force('collide', d3.forceCollide(14));

    const zoom = d3.zoom().scaleExtent([0.1, 8])
      .on('zoom', e => zoomLayer.attr('transform', e.transform));
    svg.call(zoom);

    function fitToView() {
      const b = zoomLayer.node().getBBox();
      if (!b.width || !b.height) return;
      const scale = Math.min(8, Math.max(0.1, 0.9 / Math.max(b.width / W, b.height / H)));
      const tx = W / 2 - scale * (b.x + b.width / 2);
      const ty = H / 2 - scale * (b.y + b.height / 2);
      svg.transition().duration(500).call(zoom.transform, d3.zoomIdentity.translate(tx, ty).scale(scale));
    }

    sim.on('tick', () => {
      linkSel.attr('x1', d => d.source.x).attr('y1', d => d.source.y)
             .attr('x2', d => d.target.x).attr('y2', d => d.target.y);
      nodeSel.attr('cx', d => d.x).attr('cy', d => d.y);
      labels.attr('x', d => d.x).attr('y', d => d.y);
    });
    sim.on('end', fitToView);

    nodeSel.call(d3.drag()
      .on('start', (e, d) => { if (!e.active) sim.alphaTarget(0.3).restart(); d.fx = d.x; d.fy = d.y; })
      .on('drag',  (e, d) => { d.fx = e.x; d.fy = e.y; })
      .on('end',   (e, d) => { if (!e.active) sim.alphaTarget(0); d.fx = null; d.fy = null; })
    );

    const zoomIn  = document.getElementById('hostpred-zoom-in');
    const zoomOut = document.getElementById('hostpred-zoom-out');
    const zoomFit = document.getElementById('hostpred-zoom-fit');
    if (zoomIn)  zoomIn.onclick  = () => svg.transition().duration(250).call(zoom.scaleBy, 1.4);
    if (zoomOut) zoomOut.onclick = () => svg.transition().duration(250).call(zoom.scaleBy, 1 / 1.4);
    if (zoomFit) zoomFit.onclick = fitToView;
  }

  // ── Defense & Anti-Defense ───────────────────────────────────────────────
  function _renderDefenseBar(samples) {
    const def  = typeof DEFENSE_DATA !== 'undefined' ? DEFENSE_DATA : [];
    const anti = typeof ANTIDEFENSE_DATA !== 'undefined' ? ANTIDEFENSE_DATA : [];

    mkChart('defense-bar-chart', {
      title: { text: 'Defense Systems per Sample (DefenseFinder)' },
      tooltip: { trigger: 'axis' },
      legend: { data: ['DefenseFinder', 'AntiDefenseFinder'] },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: 'Systems' },
      series: [
        { name: 'DefenseFinder',      type: 'bar', color: '#0d9488', data: _byCount(def, 'System', samples) },
        { name: 'AntiDefenseFinder',  type: 'bar', color: '#7c3aed', data: _byCount(anti, 'System', samples) },
      ],
      grid: { bottom: 70 },
    });

    const allRows = [
      ...def.map(r => ({ ...r, Tool: 'DefenseFinder', Kind: 'Defense' })),
      ...anti.map(r => ({ ...r, Tool: 'DefenseFinder', Kind: 'Anti-defense' })),
    ];
    makeTable('defense-table', allRows, [
      { key: 'sample', label: 'Sample' },
      { key: 'Bin',    label: 'Bin / Genome unit' },
      { key: 'Tool',   label: 'Tool' },
      { key: 'Kind',   label: 'Kind' },
      { key: 'System', label: 'System' },
      { key: 'System_id', label: 'System ID' },
    ], { searchId: 'defense-search' });
  }

  // ── AMR consensus ─────────────────────────────────────────────────────────
  function _renderAmrBar(samples) {
    const rows = typeof AMR_CONSENSUS !== 'undefined' ? AMR_CONSENSUS : [];

    // Hits per sample, stacked by number of tools that agreed
    const byNTools = s => ({
      three: rows.filter(r => r.sample === s && +r.n_tools === 3).length,
      two:   rows.filter(r => r.sample === s && +r.n_tools === 2).length,
      one:   rows.filter(r => r.sample === s && +r.n_tools === 1).length,
    });
    mkChart('amr-bar-chart', {
      title: { text: 'AMR Loci per Sample — Consensus' },
      tooltip: { trigger: 'axis' },
      legend: { data: ['3 tools', '2 tools', '1 tool'] },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: 'AMR loci' },
      series: [
        { name: '3 tools', type: 'bar', stack: 'amr', color: '#0d9488',
          data: samples.map(s => byNTools(s).three) },
        { name: '2 tools', type: 'bar', stack: 'amr', color: '#0891b2',
          data: samples.map(s => byNTools(s).two) },
        { name: '1 tool',  type: 'bar', stack: 'amr', color: '#94a3b8',
          data: samples.map(s => byNTools(s).one) },
      ],
      grid: { bottom: 70 },
    });

    makeTable('amr-table', rows, [
      { key: 'sample',         label: 'Sample' },
      { key: 'Bin',            label: 'Bin / Genome unit' },
      { key: 'Gene',           label: 'Gene' },
      { key: 'ARO',            label: 'ARO' },
      { key: 'Class',          label: 'Drug class' },
      { key: 'n_tools',        label: 'N tools' },
      { key: 'consensus_score',label: 'Score' },
      { key: 'tools_detected', label: 'Tools' },
    ], { searchId: 'amr-search' });
  }

  // ── Defense Systems by Phylum (Beavogui et al. 2024) ─────────────────────
  function _renderDefensePhylum() {
    const rows = _buildBinMatrix().filter(r => r.n_defense > 0);
    const byPhylum = new Map();
    rows.forEach(r => {
      const e = byPhylum.get(r.Phylum) || { systems: 0, bins: 0 };
      e.systems += r.n_defense; e.bins += 1;
      byPhylum.set(r.Phylum, e);
    });
    const sorted = Array.from(byPhylum.entries()).sort((a, b) => b[1].systems - a[1].systems).slice(0, 12);

    mkChart('defense-phylum-chart', {
      tooltip: { trigger: 'axis' },
      grid: { bottom: 90, left: 50, right: 20, top: 20 },
      xAxis: { type: 'category', data: sorted.map(d => d[0]), axisLabel: { rotate: 40 } },
      yAxis: { type: 'value', name: 'Defense systems' },
      series: [{
        type: 'bar', color: '#0d9488',
        data: sorted.map(d => d[1].systems),
        tooltip: { formatter: p => `${p.name}<br/>Systems: ${p.value}<br/>Bins: ${byPhylum.get(p.name).bins}` },
      }],
    });
  }

  // ── Defense Density by Phylum (systems per Mb, CheckM2 Genome_Size) ───────
  function _renderDefenseDensity() {
    const rows = _buildBinMatrix().filter(r => r.defense_density_per_mb !== null);
    const byPhylum = new Map();
    rows.forEach(r => {
      const e = byPhylum.get(r.Phylum) || { sum: 0, bins: 0 };
      e.sum += r.defense_density_per_mb; e.bins += 1;
      byPhylum.set(r.Phylum, e);
    });
    const sorted = Array.from(byPhylum.entries())
      .map(([phylum, e]) => [phylum, e.sum / e.bins, e.bins])
      .sort((a, b) => b[1] - a[1]).slice(0, 12);

    if (!sorted.length) {
      mkChart('defense-density-chart', { title: { text: 'No bins with both defense hits and CheckM2 genome size', textStyle: { fontSize: 12 } } });
      return;
    }
    mkChart('defense-density-chart', {
      tooltip: { trigger: 'item', formatter: p => `${p.name}<br/>Mean density: ${p.value.toFixed(2)} systems/Mb<br/>Bins: ${sorted.find(d => d[0] === p.name)[2]}` },
      grid: { bottom: 90, left: 50, right: 20, top: 20 },
      xAxis: { type: 'category', data: sorted.map(d => d[0]), axisLabel: { rotate: 40 } },
      yAxis: { type: 'value', name: 'Systems / Mb' },
      series: [{ type: 'bar', color: '#7c3aed', data: sorted.map(d => +d[1].toFixed(3)) }],
    });
  }

  // ── Defense Systems by Mechanism ──────────────────────────────────────────
  // Deliberately conservative: only systems whose mechanism is explicitly
  // confirmed are classified -- RM/CRISPR-Cas (nucleic acid degradation),
  // CBASS/Retron/Abi*-family/Toxin-Antitoxin (abortive infection -- "Abi" is
  // literally named for the mechanism), Viperin/dCTP deaminase (inhibition
  // of replication). Everything else is "Unknown / unclassified" rather
  // than guessed -- most DefenseFinder system types don't have a published
  // mechanism assignment we've verified, and a substantial share of real
  // defense systems genuinely have no characterized mechanism yet.
  function _mechanismFor(systemName) {
    const s = (systemName || '').toLowerCase();
    if (s === 'rm' || s.startsWith('rm_') || s.startsWith('rm type') || s.includes('crispr')) {
      return 'Nucleic acid degradation';
    }
    if (s === 'cbass' || s === 'retron' || s.startsWith('abi') || s.includes('toxin-antitoxin') || s === 'ta') {
      return 'Abortive infection';
    }
    if (s.includes('viperin') || s.includes('dctp deaminase') || s.includes('dctp_deaminase') || s.includes('dxtpase')) {
      return 'Inhibition of replication';
    }
    return 'Unknown / unclassified';
  }

  function _renderDefenseMechanism(samples) {
    const def = typeof DEFENSE_DATA !== 'undefined' ? DEFENSE_DATA : [];
    const order = ['Nucleic acid degradation', 'Abortive infection', 'Inhibition of replication', 'Unknown / unclassified'];
    const colors = { 'Nucleic acid degradation': '#0d9488', 'Abortive infection': '#d97706',
      'Inhibition of replication': '#2563eb', 'Unknown / unclassified': '#94a3b8' };
    const bySampleMech = new Map();
    def.forEach(r => {
      const mech = _mechanismFor(r.System);
      const key = r.sample;
      if (!bySampleMech.has(key)) bySampleMech.set(key, Object.fromEntries(order.map(o => [o, 0])));
      bySampleMech.get(key)[mech] += 1;
    });

    mkChart('defense-mechanism-chart', {
      tooltip: { trigger: 'axis' },
      legend: { data: order, top: 0 },
      grid: { bottom: 70, left: 50, right: 20, top: 50 },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: 'Defense systems' },
      series: order.map(mech => ({
        name: mech, type: 'bar', stack: 'mech', color: colors[mech],
        data: samples.map(s => (bySampleMech.get(s) || {})[mech] || 0),
      })),
    });
  }

  // ── Defense System Co-occurrence (top systems by frequency) ──────────────
  function _renderDefenseCooccurrence() {
    const bins = typeof BIN_ANNOTATIONS !== 'undefined' ? BIN_ANNOTATIONS : {};
    const freq = new Map();
    const perBinSystems = [];
    Object.values(bins).forEach(ann => {
      const systems = ann.Defense_systems || [];
      if (systems.length < 2) return;
      perBinSystems.push(systems);
      systems.forEach(sys => freq.set(sys, (freq.get(sys) || 0) + 1));
    });
    const top = Array.from(freq.entries()).sort((a, b) => b[1] - a[1]).slice(0, 12).map(d => d[0]);
    const idx = new Map(top.map((s, i) => [s, i]));
    const matrix = top.map(() => top.map(() => 0));
    perBinSystems.forEach(systems => {
      const present = systems.filter(s => idx.has(s));
      for (let i = 0; i < present.length; i++) {
        for (let j = 0; j < present.length; j++) {
          matrix[idx.get(present[i])][idx.get(present[j])] += 1;
        }
      }
    });
    const data = [];
    for (let i = 0; i < top.length; i++) {
      for (let j = 0; j < top.length; j++) {
        if (matrix[i][j] > 0) data.push([j, i, matrix[i][j]]);
      }
    }
    const maxV = Math.max(1, ...data.map(d => d[2]));

    if (!top.length) {
      mkChart('defense-cooccur-chart', { title: { text: 'No bins with 2+ defense systems', textStyle: { fontSize: 12 } } });
      return;
    }
    mkChart('defense-cooccur-chart', {
      tooltip: { trigger: 'item', position: 'top', formatter: p => `${top[p.value[1]]} &harr; ${top[p.value[0]]}<br/>Co-occurring bins: ${p.value[2]}` },
      grid: { left: 110, right: 20, top: 10, bottom: 80 },
      xAxis: { type: 'category', data: top, axisLabel: { rotate: 60, fontSize: 10 }, splitArea: { show: true } },
      yAxis: { type: 'category', data: top, axisLabel: { fontSize: 10 }, splitArea: { show: true } },
      visualMap: { min: 0, max: maxV, calculable: true, orient: 'horizontal', left: 'center', bottom: 0,
        inRange: { color: ['#f1f5f9', '#0d9488'] } },
      series: [{ type: 'heatmap', data, label: { show: false } }],
    });
  }

  // ── Defense System Correlation (gene-count proxy, FDR-corrected) ─────────
  // Closer to Han et al. 2026 Fig 4a than the presence/absence co-occurrence
  // heatmap above: a real Pearson r per system pair across every bin (not
  // just "found in the same bin Y/N"), with a two-sided significance test
  // (Fisher z) and Benjamini-Hochberg FDR correction across all pairs
  // tested -- same statistical approach as the paper.
  // Approximation, stated plainly: the paper correlates Salmon-derived gene
  // abundance (GPM); we don't run per-gene read quantification, so this
  // uses DefenseFinder's own per-system gene count (already collected) as
  // the abundance proxy instead. Real signal (multi-copy/expanded systems
  // do get a higher value), just not sequencing-depth-normalized.
  function _renderDefenseCorrelation() {
    const def = typeof DEFENSE_DATA !== 'undefined' ? DEFENSE_DATA : [];
    const byBinSystem = new Map();
    def.forEach(r => {
      const key = `${r.sample}::${r.Bin}`;
      if (!byBinSystem.has(key)) byBinSystem.set(key, new Map());
      const genes = parseInt(r.Genes, 10) || 1;
      const m = byBinSystem.get(key);
      m.set(r.System, (m.get(r.System) || 0) + genes);
    });
    const bins = [...byBinSystem.keys()];
    const freq = new Map();
    byBinSystem.forEach(m => m.forEach((_, sys) => freq.set(sys, (freq.get(sys) || 0) + 1)));
    const top = [...freq.entries()].sort((a, b) => b[1] - a[1]).slice(0, 15).map(d => d[0]);

    const noteEl = document.getElementById('defense-correlation-note');
    if (top.length < 3 || bins.length < 4) {
      mkChart('defense-correlation-chart', { title: { text: 'Not enough systems/bins for a correlation matrix yet', textStyle: { fontSize: 12 } } });
      if (noteEl) noteEl.textContent = '';
      return;
    }

    const vectors = top.map(sys => bins.map(b => byBinSystem.get(b).get(sys) || 0));
    const pairs = [];
    for (let i = 0; i < top.length; i++) {
      for (let j = 0; j < top.length; j++) {
        if (i === j) continue;
        const r = _pearson(vectors[i], vectors[j]);
        pairs.push({ i, j, r: r === null ? 0 : r });
      }
    }
    const pvals = pairs.map(p => _pearsonPValue(p.r, bins.length));
    const qvals = _fdrBH(pvals);
    pairs.forEach((p, k) => { p.p = pvals[k]; p.q = qvals[k]; });

    const sigPairs = pairs.filter(p => p.i < p.j && p.q < 0.05).length;
    if (noteEl) noteEl.textContent =
      `${top.length} systems × ${bins.length} bins — ${sigPairs} pair(s) significant after FDR correction (q < 0.05).`;

    const data = pairs.map(p => [p.j, p.i, +p.r.toFixed(3), +p.q.toFixed(4)]);
    mkChart('defense-correlation-chart', {
      tooltip: { trigger: 'item', formatter: pt => {
        const [x, y, r, q] = pt.data;
        return `${top[y]} &harr; ${top[x]}<br/>r = ${r}<br/>FDR q = ${q}${q < 0.05 ? ' (significant)' : ''}`;
      }},
      grid: { left: 110, right: 20, top: 10, bottom: 80 },
      xAxis: { type: 'category', data: top, axisLabel: { rotate: 60, fontSize: 10 }, splitArea: { show: true } },
      yAxis: { type: 'category', data: top, axisLabel: { fontSize: 10 }, splitArea: { show: true } },
      visualMap: { min: -1, max: 1, calculable: true, orient: 'horizontal', left: 'center', bottom: 0,
        inRange: { color: ['#b91c1c', '#f1f5f9', '#0d9488'] } },
      series: [{ type: 'heatmap', data, label: { show: false } }],
    });
  }

  // ── AMR Drug Classes (consensus) ──────────────────────────────────────────
  function _renderAmrClasses() {
    const rows = typeof AMR_CONSENSUS !== 'undefined' ? AMR_CONSENSUS : [];
    const byClass = new Map();
    rows.forEach(r => {
      const cls = (r.Class || 'Unclassified').trim() || 'Unclassified';
      const e = byClass.get(cls) || { three: 0, two: 0, one: 0 };
      const n = +r.n_tools;
      if (n === 3) e.three += 1; else if (n === 2) e.two += 1; else e.one += 1;
      byClass.set(cls, e);
    });
    const sorted = Array.from(byClass.entries())
      .sort((a, b) => (b[1].three + b[1].two + b[1].one) - (a[1].three + a[1].two + a[1].one))
      .slice(0, 12);

    mkChart('amr-class-chart', {
      tooltip: { trigger: 'axis' },
      legend: { data: ['3 tools', '2 tools', '1 tool'], top: 0 },
      grid: { bottom: 90, left: 50, right: 20, top: 55 },
      xAxis: { type: 'category', data: sorted.map(d => d[0]), axisLabel: { rotate: 40 } },
      yAxis: { type: 'value', name: 'AMR loci' },
      series: [
        { name: '3 tools', type: 'bar', stack: 'amr', color: '#0d9488', data: sorted.map(d => d[1].three) },
        { name: '2 tools', type: 'bar', stack: 'amr', color: '#0891b2', data: sorted.map(d => d[1].two) },
        { name: '1 tool',  type: 'bar', stack: 'amr', color: '#94a3b8', data: sorted.map(d => d[1].one) },
      ],
    });
  }

  // ── AMR by Contig Origin (viral genome vs. prokaryotic MAG vs. unbinned) ──
  function _renderAmrOrigin() {
    const rows   = typeof AMR_CONSENSUS !== 'undefined' ? AMR_CONSENSUS : [];
    const checkv = typeof CHECKV !== 'undefined' ? CHECKV : {};
    const gtdb   = typeof GTDB_DATA !== 'undefined' ? GTDB_DATA : [];

    const viralContigs = new Set();
    Object.entries(checkv).forEach(([s, rows]) => (rows || []).forEach(r => {
      const id = r.contig_id || r.contig;
      if (id) viralContigs.add(`${s}::${id}`);
    }));
    const prokBins = new Set(gtdb.map(r => `${r.sample}::${r.Bin}`));

    function origin(r) {
      const key = `${r.sample}::${r.Bin}`;
      if (viralContigs.has(key)) return 'Viral genome';
      if (prokBins.has(key)) return 'Prokaryotic MAG';
      if (r.Bin === 'contigs_pseudogenome') return 'Unbinned / mixed contigs';
      return 'Other / unbinned';
    }

    const byOrigin = new Map();
    const add = (rows_) => rows_.forEach(r => {
      const o = origin(r);
      const e = byOrigin.get(o) || { count: 0 };
      e.count += 1;
      byOrigin.set(o, e);
    });
    const cats = ['Viral genome', 'Prokaryotic MAG', 'Unbinned / mixed contigs', 'Other / unbinned']
      .filter(c => byOrigin.has(c));

    mkChart('amr-origin-chart', {
      tooltip: { trigger: 'axis' },
      grid: { bottom: 70, left: 50, right: 20, top: 20 },
      xAxis: { type: 'category', data: cats, axisLabel: { rotate: 20 } },
      yAxis: { type: 'value', name: 'AMR loci' },
      series: [
        { type: 'bar', color: '#0891b2', data: cats.map(c => byOrigin.get(c).count) },
      ],
    });
  }

  // ── Cross-Analysis: Defense load vs. AMR load, and vs. PHIST host load ───
  function _renderCrossAnalysis() {
    const matrix = _buildBinMatrix();
    const domainColor = { Bacteria: '#0d9488', Archaea: '#7c3aed', Unclassified: '#94a3b8' };

    // Defense vs. AMR (consensus) per bin
    const withAmr = matrix.filter(r => r.n_defense > 0 || r.n_amr > 0);
    const r1 = _pearson(withAmr.map(r => r.n_defense), withAmr.map(r => r.n_amr));
    _renderCorrBanner('cross-corr-defense-amr', r1,
      'defense systems', 'AMR loci', withAmr.length,
      'A negative trend would suggest bins investing more in anti-phage defense tend to carry fewer AMR genes.');

    mkChart('defense-amr-scatter', {
      tooltip: { trigger: 'item', formatter: p => `${p.data[3]} (${p.data[2]})<br/>Defense systems: ${p.data[0]}<br/>AMR loci: ${p.data[1]}` },
      legend: { data: Object.keys(domainColor) },
      xAxis: { type: 'value', name: 'Defense systems / bin', minInterval: 1 },
      yAxis: { type: 'value', name: 'AMR loci / bin', minInterval: 1 },
      series: Object.keys(domainColor).map(dom => ({
        name: dom, type: 'scatter', symbolSize: 9, color: domainColor[dom],
        data: withAmr.filter(r => r.Domain === dom).map(r => [r.n_defense, r.n_amr, r.Bin, r.sample]),
      })),
    });

    // Defense vs. PHIST host load (distinct viruses predicted per host bin)
    const phist = typeof PHIST_DATA !== 'undefined' ? PHIST_DATA : [];
    const hostViruses = new Map();
    phist.forEach(p => {
      const key = `${p.sample}::${p.Host}`;
      if (!hostViruses.has(key)) hostViruses.set(key, new Set());
      if (p.Virus) hostViruses.get(key).add(p.Virus);
    });
    const withHost = matrix
      .map(r => ({ ...r, n_viruses: (hostViruses.get(`${r.sample}::${r.Bin}`) || new Set()).size }))
      .filter(r => r.n_defense > 0 || r.n_viruses > 0);
    const r2 = _pearson(withHost.map(r => r.n_defense), withHost.map(r => r.n_viruses));
    _renderCorrBanner('cross-corr-defense-host', r2,
      'defense systems', 'predicted phage hosts (PHIST)', withHost.length,
      'A negative trend would suggest bins with richer defense repertoires are predicted as hosts for fewer distinct phages.');

    mkChart('defense-hostload-scatter', {
      tooltip: { trigger: 'item', formatter: p => `${p.data[3]} (${p.data[2]})<br/>Defense systems: ${p.data[0]}<br/>Predicted phages: ${p.data[1]}` },
      legend: { data: Object.keys(domainColor) },
      xAxis: { type: 'value', name: 'Defense systems / bin', minInterval: 1 },
      yAxis: { type: 'value', name: 'Distinct phages predicted (PHIST)', minInterval: 1 },
      series: Object.keys(domainColor).map(dom => ({
        name: dom, type: 'scatter', symbolSize: 9, color: domainColor[dom],
        data: withHost.filter(r => r.Domain === dom).map(r => [r.n_defense, r.n_viruses, r.Bin, r.sample]),
      })),
    });
  }

  function _renderCorrBanner(elId, r, xLabel, yLabel, n, note) {
    const el = document.getElementById(elId);
    if (!el) return;
    if (r === null || n < 3) {
      el.innerHTML = `<span class="corr-text">Not enough bins with both ${xLabel} and ${yLabel} to compute a correlation (n=${n}).</span>`;
      return;
    }
    const dir = r < -0.1 ? 'inverse' : r > 0.1 ? 'positive' : 'negligible';
    el.classList.toggle('warn', Math.abs(r) < 0.1);
    el.innerHTML = `
      <span class="corr-stat">r = ${r.toFixed(2)}</span>
      <span class="corr-text"><strong>${dir.charAt(0).toUpperCase() + dir.slice(1)} relationship</strong> between ${xLabel} and ${yLabel} across ${n} bins.
      <small>${note}</small></span>`;
  }

  // ── Viral anti-defense (Han et al. 2026) — DefenseFinder vs. dbAPIS ───────
  // Same tool/models as the host-side DefenseFinder run, just pointed at
  // viral proteins, plus dbAPIS as a complementary sequence-similarity
  // detector. Kept in separate columns/series, never merged (same
  // never-merge-tiers rule as AMR curated/exploratory).
  function _renderViralAntidefense(samples) {
    const df     = typeof ANTIDEFENSE_VIRAL_DF !== 'undefined' ? ANTIDEFENSE_VIRAL_DF : [];
    const dbapis = typeof ANTIDEFENSE_VIRAL_DBAPIS !== 'undefined' ? ANTIDEFENSE_VIRAL_DBAPIS : [];

    mkChart('viral-antidefense-chart', {
      tooltip: { trigger: 'axis' },
      legend: { data: ['DefenseFinder (viral)', 'dbAPIS'] },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: 'Anti-defense hits' },
      series: [
        { name: 'DefenseFinder (viral)', type: 'bar', color: '#7c3aed', data: _byCount(df, 'System', samples) },
        { name: 'dbAPIS',                type: 'bar', color: '#db2777', data: _byCount(dbapis, 'Virus', samples) },
      ],
      grid: { bottom: 70 },
    });

    const allRows = [
      ...df.map(r => ({ ...r, Hit: r.System, Detail: r.System_id })),
      // Gene/Defense_system_inhibited come from dbAPIS's own
      // seed_and_familyrep_all_infor.tsv (family -> gene name + readable
      // inhibited-system label, e.g. 'APIS331' -> restriction-modification
      // system) -- falls back to the bare family/gene ID from the hit
      // itself if that mapping file isn't downloaded yet.
      ...dbapis.map(r => ({ ...r, Hit: r.Gene,
        Detail: [r.Defense_system_inhibited, `pident=${r.Pident} e=${r.Evalue}`].filter(Boolean).join(' — ') })),
    ];
    makeTable('viral-antidefense-table', allRows, [
      { key: 'sample', label: 'Sample' },
      { key: 'Virus',  label: 'Viral contig' },
      { key: 'Source', label: 'Tool' },
      { key: 'Hit',    label: 'Hit / System' },
      { key: 'Detail', label: 'Detail' },
    ], { searchId: 'viral-antidefense-search' });
  }

  // ── Defense Islands (≥5 genes, ≥3 systems, ≤10-gene window) ───────────────
  // Definition from Beavogui et al. 2024 / Han et al. 2026. Display style
  // borrows the collapsible-card + linear gene-arrow pattern from
  // ViralQuest's sequence viewer (vq-seq-card) — appropriate here because,
  // like a viral genome card, an island is a small genomic *region* best
  // browsed one at a time rather than always-expanded.
  const _ISLAND_COLORS = ['#0d9488', '#7c3aed', '#d97706', '#0891b2', '#db2777', '#16a34a', '#ca8a04', '#2563eb'];

  function _renderDefenseIslands() {
    const islands = typeof DEFENSE_ISLANDS !== 'undefined' ? DEFENSE_ISLANDS : [];
    const cont = document.getElementById('defense-islands-list');
    const countEl = document.getElementById('defense-islands-count');
    if (countEl) countEl.textContent = islands.length;
    if (!cont) return;
    cont.innerHTML = '';

    if (!islands.length) {
      cont.innerHTML = '<p style="color:var(--text-muted);padding:1rem">No defense islands found '
        + '(≥5 genes from ≥3 systems within a 10-gene window).</p>';
      return;
    }

    // Stable color per system name across every island card.
    const sysColor = new Map();
    function colorFor(sys) {
      if (!sys) return 'var(--border-dark)';
      if (!sysColor.has(sys)) sysColor.set(sys, _ISLAND_COLORS[sysColor.size % _ISLAND_COLORS.length]);
      return sysColor.get(sys);
    }

    islands
      .slice().sort((a, b) => b.n_systems - a.n_systems || b.n_genes - a.n_genes)
      .forEach((isl, i) => {
        const card = document.createElement('div');
        card.className = 'region-card';

        const head = document.createElement('div');
        head.className = 'region-card__head';
        head.innerHTML = `
          <span class="region-card__id">${isl.Bin}</span>
          <span class="region-card__meta">
            <span>${isl.Contig}</span>
            <span>${isl.n_genes} genes</span>
            <span>${isl.n_systems} systems: ${isl.Systems.join(', ')}</span>
            <span>${isl.sample}</span>
          </span>
          <svg class="region-card__chevron" width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><polyline points="9 18 15 12 9 6"></polyline></svg>`;
        const body = document.createElement('div');
        body.className = 'region-card__body';
        card.appendChild(head);
        card.appendChild(body);
        cont.appendChild(card);

        head.addEventListener('click', () => {
          const wasOpen = card.classList.contains('open');
          card.classList.toggle('open');
          if (!wasOpen && !body.dataset.rendered) {
            _renderIslandGeneMap(body, isl, colorFor);
            body.dataset.rendered = '1';
          }
        });

        if (i === 0) { card.classList.add('open'); _renderIslandGeneMap(body, isl, colorFor); body.dataset.rendered = '1'; }
      });
  }

  // Genome-browser track for one island. Prodigal writes each gene's locus into
  // the .faa header (`>id # start # end # strand # ...`), so the loader carries
  // real bp coordinates and strand — drawn here as an IGV-style track: a bp
  // ruler, strand-aware arrows whose width is the gene's true length, real
  // intergenic gaps, and a system-span lane underneath. When coordinates are
  // absent (non-Prodigal input) it falls back to evenly spaced gene order.
  function _renderIslandGeneMap(body, isl, colorFor) {
    const wrap = document.createElement('div');
    wrap.className = 'region-genemap';
    body.appendChild(wrap);

    const genes = (isl.window_genes || []).slice();
    if (!genes.length) return;
    const dark = document.documentElement.dataset.theme === 'dark';
    const hasCoords = isl.Start_bp != null && isl.End_bp != null && isl.End_bp > isl.Start_bp;

    // Synthesize coordinates when the input had none, so one code path draws both.
    if (!hasCoords) {
      const SYN = 1000, GAP = 200;
      genes.forEach((g, i) => { g.Start = i * (SYN + GAP) + 1; g.End = g.Start + SYN; g.Strand = g.Strand || 1; });
    }
    const x0 = hasCoords ? isl.Start_bp : genes[0].Start;
    const x1 = hasCoords ? isl.End_bp   : genes[genes.length - 1].End;

    const PAD_L = 14, PAD_R = 14, H = 168;
    const W = Math.max(560, Math.min(1180, (wrap.clientWidth || 900)));
    const innerW = W - PAD_L - PAD_R;
    const span = Math.max(1, x1 - x0);
    const pad = span * 0.02;                       // small breathing room each side
    const x = d3.scaleLinear().domain([x0 - pad, x1 + pad]).range([0, innerW]);

    const ink   = dark ? '#94a3b8' : '#64748b';
    const rule  = dark ? '#1e293b' : '#e2e8f0';
    const other = dark ? '#334155' : '#cbd5e1';

    const svg = d3.select(wrap).append('svg')
      .attr('width', W).attr('height', H)
      .attr('viewBox', `0 0 ${W} ${H}`).style('max-width', '100%');
    const tip = d3.select(wrap).append('div').attr('class', 'd3-tooltip').style('display', 'none');

    const root = svg.append('g').attr('transform', `translate(${PAD_L},0)`);
    const clipId = 'iclip' + Math.random().toString(36).slice(2, 9);
    root.append('clipPath').attr('id', clipId)
      .append('rect').attr('x', -2).attr('y', 0).attr('width', innerW + 4).attr('height', H);
    const view = root.append('g').attr('clip-path', `url(#${clipId})`);

    const AXIS_Y = 30, GENE_Y = 74, GENE_H = 26, SYS_Y = 118;

    // ── bp ruler ────────────────────────────────────────────────────────────
    const axisG = root.append('g').attr('transform', `translate(0,${AXIS_Y})`);
    function fmtBp(v) {
      if (span > 20000) return (v / 1000).toFixed(0) + ' kb';
      if (span > 2000)  return (v / 1000).toFixed(1) + ' kb';
      return Math.round(v).toLocaleString() + ' bp';
    }
    function drawAxis(xs) {
      axisG.selectAll('*').remove();
      axisG.append('line').attr('x1', 0).attr('x2', innerW).attr('stroke', rule).attr('stroke-width', 1);
      if (!hasCoords) {
        axisG.append('text').attr('x', 0).attr('y', -9).attr('font-size', 10).attr('fill', ink)
          .text('gene order (no coordinates in input)');
        return;
      }
      xs.ticks(Math.max(3, Math.round(innerW / 130))).forEach(t => {
        const px = xs(t);
        if (px < -2 || px > innerW + 2) return;
        axisG.append('line').attr('x1', px).attr('x2', px).attr('y1', -5).attr('y2', 0)
          .attr('stroke', ink).attr('stroke-width', 1);
        axisG.append('text').attr('x', px).attr('y', -9).attr('text-anchor', 'middle')
          .attr('font-size', 10).attr('fill', ink)
          .attr('font-family', 'ui-monospace, Menlo, monospace').text(fmtBp(t));
      });
    }

    // ── centre line + genes ─────────────────────────────────────────────────
    view.append('line').attr('class', 'axis-line')
      .attr('x1', 0).attr('x2', innerW).attr('y1', GENE_Y + GENE_H / 2).attr('y2', GENE_Y + GENE_H / 2)
      .attr('stroke', rule).attr('stroke-width', 2);

    const geneG = view.append('g');
    const sysG  = view.append('g');

    function arrowPath(xa, xb, strand, h) {
      const w = Math.max(3, xb - xa);
      const tipw = Math.min(10, w * 0.4);
      const top = -h / 2, bot = h / 2;
      return strand === -1
        ? `M${xa + w},${top} L${xa + tipw},${top} L${xa},0 L${xa + tipw},${bot} L${xa + w},${bot} Z`
        : `M${xa},${top} L${xa + w - tipw},${top} L${xa + w},0 L${xa + w - tipw},${bot} L${xa},${bot} Z`;
    }

    function draw(xs) {
      geneG.selectAll('*').remove();
      sysG.selectAll('*').remove();
      drawAxis(xs);

      genes.forEach(gene => {
        const xa = xs(gene.Start), xb = xs(gene.End);
        if (xb < -20 || xa > innerW + 20) return;
        const isDef = !!gene.System;
        geneG.append('path')
          .attr('transform', `translate(0,${GENE_Y + GENE_H / 2})`)
          .attr('d', arrowPath(xa, xb, gene.Strand, GENE_H))
          .attr('fill', isDef ? colorFor(gene.System) : other)
          .attr('fill-opacity', isDef ? 0.92 : 0.5)
          .attr('stroke', dark ? '#0b1220' : '#ffffff')
          .attr('stroke-width', 2)
          .style('cursor', 'pointer')
          .on('mousemove', (e) => {
            const len = (gene.End - gene.Start + 1);
            tip.style('display', 'block')
              .style('left', (e.offsetX + 14) + 'px').style('top', (e.offsetY - 6) + 'px')
              .html(`<strong>${gene.Protein}</strong>`
                + (isDef ? `<br>System: ${gene.System}` : '<br>No defense system hit')
                + (hasCoords ? `<br>${gene.Start.toLocaleString()}–${gene.End.toLocaleString()} `
                    + `(${gene.Strand === -1 ? '−' : '+'} strand, ${len.toLocaleString()} bp)` : ''));
          })
          .on('mouseout', () => tip.style('display', 'none'));
      });

      // ── system-span lane: one bar covering each system's first→last gene ──
      const spans = new Map();
      genes.forEach(g => {
        if (!g.System) return;
        const key = g.System_id || g.System;
        const cur = spans.get(key);
        if (!cur) spans.set(key, { sys: g.System, s: g.Start, e: g.End });
        else { cur.s = Math.min(cur.s, g.Start); cur.e = Math.max(cur.e, g.End); }
      });
      [...spans.values()].sort((a, b) => a.s - b.s).forEach((sp, i) => {
        const xa = xs(sp.s), xb = xs(sp.e);
        const row = i % 2;                       // stagger to avoid overlap
        const yy = SYS_Y + row * 15;
        sysG.append('rect')
          .attr('x', xa).attr('y', yy).attr('width', Math.max(2, xb - xa)).attr('height', 7)
          .attr('rx', 3).attr('fill', colorFor(sp.sys)).attr('fill-opacity', 0.85);
        if (xb - xa > 44) {
          sysG.append('text').attr('x', xa + 4).attr('y', yy - 3)
            .attr('font-size', 9.5).attr('fill', ink).text(sp.sys);
        }
      });
    }

    draw(x);

    // ── zoom / pan (x only) ─────────────────────────────────────────────────
    if (hasCoords) {
      svg.call(d3.zoom().scaleExtent([1, 40])
        .translateExtent([[0, 0], [innerW, H]]).extent([[0, 0], [innerW, H]])
        .on('zoom', (ev) => draw(ev.transform.rescaleX(x))))
        .style('cursor', 'grab');
      svg.append('text').attr('x', W - PAD_R).attr('y', H - 5).attr('text-anchor', 'end')
        .attr('font-size', 9.5).attr('fill', ink).text('scroll to zoom · drag to pan');
    }

    const legend = document.createElement('div');
    legend.className = 'region-legend';
    isl.Systems.forEach(sys => {
      legend.innerHTML += `<span class="region-legend__item">
        <span class="region-legend__swatch" style="background:${colorFor(sys)}"></span>${sys}</span>`;
    });
    legend.innerHTML += `<span class="region-legend__item">
      <span class="region-legend__swatch" style="background:${dark ? '#334155' : '#cbd5e1'}"></span>Other gene (no defense hit)</span>`;
    body.appendChild(legend);
  }

  // ── Host <-> Virus <-> Defense <-> AMR matrix ────────────────────────────
  // Defense/antidefense/AMR detail lives in BIN_ANNOTATIONS, keyed by
  // 'sample::Host' -- looked up per row here instead of being pre-embedded
  // on every HOST_DEFENSE_LINKS row (a host predicted for hundreds of
  // viruses would otherwise duplicate its full system/gene list hundreds
  // of times in the page's JSON payload).
  function _renderMatrix(sample) {
    const all  = typeof HOST_DEFENSE_LINKS !== 'undefined' ? HOST_DEFENSE_LINKS : [];
    const bins = typeof BIN_ANNOTATIONS !== 'undefined' ? BIN_ANNOTATIONS : {};
    // Viral anti-defense (Han et al. 2026 Fig 6a/6c): grouped by sample::Virus
    // so a virus predicted host to many bins doesn't duplicate its
    // anti-defense list per row -- same dedup rationale as BIN_ANNOTATIONS.
    const viralAnti = new Map();
    const addVirAnti = (rows, label, nameKey) => rows.forEach(r => {
      const key = `${r.sample}::${r.Virus}`;
      if (!viralAnti.has(key)) viralAnti.set(key, { DefenseFinder: new Set(), dbAPIS: new Set() });
      viralAnti.get(key)[label].add(r[nameKey] || '');
    });
    addVirAnti(typeof ANTIDEFENSE_VIRAL_DF !== 'undefined' ? ANTIDEFENSE_VIRAL_DF : [], 'DefenseFinder', 'System');
    addVirAnti(typeof ANTIDEFENSE_VIRAL_DBAPIS !== 'undefined' ? ANTIDEFENSE_VIRAL_DBAPIS : [], 'dbAPIS', 'Gene');

    const isAll = sample === '__all__';
    const rows = (isAll ? all : all.filter(r => r.sample === sample)).map(r => {
      const ann = bins[`${r.sample}::${r.Host}`] || {};
      const va  = viralAnti.get(`${r.sample}::${r.Virus}`) || { DefenseFinder: new Set(), dbAPIS: new Set() };
      return {
        ...r,
        Host_taxonomy_display: r.Host_genus || r.Host_species || r.Host_taxonomy || '—',
        Defense_display:       (ann.Defense_systems || []).join(', ') || '—',
        Antidefense_display:   (ann.Antidefense_systems || []).join(', ') || '—',
        AMR_display:           (ann.AMR_genes || []).join(', ') || '—',
        Viral_antidefense_display: [...va.DefenseFinder, ...va.dbAPIS].filter(Boolean).join(', ') || '—',
        Pair_flag: ((ann.Defense_systems || []).length > 0 && (va.DefenseFinder.size + va.dbAPIS.size) > 0) ? 'Host defends & virus counters' : '',
      };
    });

    makeTable('hostdefense-matrix-table', rows, [
      { key: 'sample',  label: 'Sample' },
      { key: 'Virus',   label: 'Virus' },
      { key: 'Host',    label: 'Host (MAG)' },
      { key: 'Host_taxonomy_display', label: 'Host Taxonomy' },
      { key: 'Defense_display',       label: 'Defense Systems' },
      { key: 'Antidefense_display',   label: 'Anti-Defense Systems (host)' },
      { key: 'Viral_antidefense_display', label: 'Anti-Defense Genes (virus)' },
      { key: 'Pair_flag', label: 'Arms-race pair' },
      { key: 'AMR_display', label: 'AMR genes' },
    ], { searchId: 'hostdefense-matrix-search' });
  }

})();
