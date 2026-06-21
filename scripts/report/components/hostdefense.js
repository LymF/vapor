/* hostdefense.js — Host Prediction (PHIST), Defense/Anti-Defense systems
   (DefenseFinder, with built-in AntiDefenseFinder), AMR (AMRFinderPlus+RGI
   curated vs. DeepARG exploratory), and the Host<->Virus<->Defense<->AMR
   cross-link matrix. */
(function () {
  'use strict';

  window.renderHostDefense = function () {
    const samples = typeof SAMPLES !== 'undefined' ? SAMPLES : [];
    _renderDefenseBar(samples);
    _renderAmrBar(samples);
    _renderDefensePhylum();
    _renderDefenseCooccurrence();
    _renderAmrClasses();
    _renderAmrOrigin();
    _renderCrossAnalysis();
    _renderViralAntidefense(samples);
    _renderDefenseIslands();
    makeSampleDropdown('sample-sel-hostpred', _renderHostPrediction, { allSamples: true });
    makeSampleDropdown('sample-sel-hostdefense-matrix', _renderMatrix, { allSamples: true });
  };

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

  // One row per (sample, Bin) with defense/anti-defense/AMR counts + GTDB-Tk
  // taxonomy — built client-side from BIN_ANNOTATIONS + GTDB_DATA so no new
  // Snakemake/Python plumbing is needed for the cross-analysis charts below.
  function _buildBinMatrix() {
    const bins = typeof BIN_ANNOTATIONS !== 'undefined' ? BIN_ANNOTATIONS : {};
    const gtdb = typeof GTDB_DATA !== 'undefined' ? GTDB_DATA : [];
    const gtdbByKey = new Map();
    gtdb.forEach(r => gtdbByKey.set(`${r.sample}::${r.Bin}`, r));
    return Object.entries(bins).map(([key, ann]) => {
      const sep = key.indexOf('::');
      const sample = key.slice(0, sep), Bin = key.slice(sep + 2);
      const tax = gtdbByKey.get(key) || {};
      return {
        sample, Bin,
        Domain: tax.Domain || 'Unclassified',
        Phylum: tax.Phylum || 'Unclassified',
        n_defense:        (ann.Defense_systems || []).length,
        n_antidefense:    (ann.Antidefense_systems || []).length,
        n_amr_curated:    (ann.AMR_curated || []).length,
        n_amr_exploratory:(ann.AMR_exploratory || []).length,
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

  // ── AMR (curated vs. exploratory) ────────────────────────────────────────
  function _renderAmrBar(samples) {
    const curated     = typeof AMR_DATA !== 'undefined' ? AMR_DATA : [];
    const exploratory  = typeof DEEPARG_DATA !== 'undefined' ? DEEPARG_DATA : [];

    mkChart('amr-bar-chart', {
      title: { text: 'AMR Genes per Sample — Curated vs. Exploratory' },
      tooltip: { trigger: 'axis' },
      legend: { data: ['Curated (AMRFinderPlus + RGI/CARD)', 'Exploratory (DeepARG)'] },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: 'AMR gene hits' },
      series: [
        { name: 'Curated (AMRFinderPlus + RGI/CARD)', type: 'bar', color: '#0891b2',
          data: _byCount(curated, 'Gene', samples) },
        { name: 'Exploratory (DeepARG)', type: 'bar', color: '#ef4444',
          data: _byCount(exploratory, 'Gene', samples) },
      ],
      grid: { bottom: 70 },
    });

    const allRows = [
      ...curated.map(r => ({ ...r, Tier_label: 'Curated' })),
      ...exploratory.map(r => ({ ...r, Tier_label: 'Exploratory' })),
    ];
    makeTable('amr-table', allRows, [
      { key: 'sample', label: 'Sample' },
      { key: 'Bin',    label: 'Bin / Genome unit' },
      { key: 'Source', label: 'Tool' },
      { key: 'Tier_label', label: 'Tier' },
      { key: 'Gene',   label: 'Gene' },
      { key: 'Class',  label: 'Drug class' },
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

  // ── AMR Gene Classes — curated vs. exploratory ────────────────────────────
  function _renderAmrClasses() {
    const curated     = typeof AMR_DATA !== 'undefined' ? AMR_DATA : [];
    const exploratory  = typeof DEEPARG_DATA !== 'undefined' ? DEEPARG_DATA : [];
    const byClass = new Map();
    const add = (rows, tier) => rows.forEach(r => {
      const cls = (r.Class || 'Unclassified').trim() || 'Unclassified';
      const e = byClass.get(cls) || { curated: 0, exploratory: 0 };
      e[tier] += 1;
      byClass.set(cls, e);
    });
    add(curated, 'curated');
    add(exploratory, 'exploratory');
    const sorted = Array.from(byClass.entries())
      .sort((a, b) => (b[1].curated + b[1].exploratory) - (a[1].curated + a[1].exploratory))
      .slice(0, 12);

    mkChart('amr-class-chart', {
      tooltip: { trigger: 'axis' },
      legend: { data: ['Curated', 'Exploratory'], top: 0 },
      grid: { bottom: 90, left: 50, right: 20, top: 55 },
      xAxis: { type: 'category', data: sorted.map(d => d[0]), axisLabel: { rotate: 40 } },
      yAxis: { type: 'value', name: 'AMR gene hits' },
      series: [
        { name: 'Curated',     type: 'bar', stack: 'amr', color: '#0891b2', data: sorted.map(d => d[1].curated) },
        { name: 'Exploratory', type: 'bar', stack: 'amr', color: '#ef4444', data: sorted.map(d => d[1].exploratory) },
      ],
    });
  }

  // ── AMR by Contig Origin (viral genome vs. prokaryotic MAG vs. unbinned) ──
  function _renderAmrOrigin() {
    const curated     = typeof AMR_DATA !== 'undefined' ? AMR_DATA : [];
    const exploratory  = typeof DEEPARG_DATA !== 'undefined' ? DEEPARG_DATA : [];
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
    const add = (rows, tier) => rows.forEach(r => {
      const o = origin(r);
      const e = byOrigin.get(o) || { curated: 0, exploratory: 0 };
      e[tier] += 1;
      byOrigin.set(o, e);
    });
    add(curated, 'curated');
    add(exploratory, 'exploratory');
    const cats = ['Viral genome', 'Prokaryotic MAG', 'Unbinned / mixed contigs', 'Other / unbinned']
      .filter(c => byOrigin.has(c));

    mkChart('amr-origin-chart', {
      tooltip: { trigger: 'axis' },
      legend: { data: ['Curated', 'Exploratory'], top: 0 },
      grid: { bottom: 70, left: 50, right: 20, top: 55 },
      xAxis: { type: 'category', data: cats, axisLabel: { rotate: 20 } },
      yAxis: { type: 'value', name: 'AMR gene hits' },
      series: [
        { name: 'Curated',     type: 'bar', stack: 'amr', color: '#0891b2', data: cats.map(c => byOrigin.get(c).curated) },
        { name: 'Exploratory', type: 'bar', stack: 'amr', color: '#ef4444', data: cats.map(c => byOrigin.get(c).exploratory) },
      ],
    });
  }

  // ── Cross-Analysis: Defense load vs. AMR load, and vs. PHIST host load ───
  function _renderCrossAnalysis() {
    const matrix = _buildBinMatrix();
    const domainColor = { Bacteria: '#0d9488', Archaea: '#7c3aed', Unclassified: '#94a3b8' };

    // Defense vs. AMR (curated) per bin
    const withAmr = matrix.filter(r => r.n_defense > 0 || r.n_amr_curated > 0);
    const r1 = _pearson(withAmr.map(r => r.n_defense), withAmr.map(r => r.n_amr_curated));
    _renderCorrBanner('cross-corr-defense-amr', r1,
      'defense systems', 'curated AMR genes', withAmr.length,
      'Negative r is consistent with Li et al. 2025 (aquifer metagenomes): bins investing more in anti-phage defense tend to carry fewer AMR genes.');

    mkChart('defense-amr-scatter', {
      tooltip: { trigger: 'item', formatter: p => `${p.data[3]} (${p.data[2]})<br/>Defense systems: ${p.data[0]}<br/>AMR genes (curated): ${p.data[1]}` },
      legend: { data: Object.keys(domainColor) },
      xAxis: { type: 'value', name: 'Defense systems / bin', minInterval: 1 },
      yAxis: { type: 'value', name: 'AMR genes (curated) / bin', minInterval: 1 },
      series: Object.keys(domainColor).map(dom => ({
        name: dom, type: 'scatter', symbolSize: 9, color: domainColor[dom],
        data: withAmr.filter(r => r.Domain === dom).map(r => [r.n_defense, r.n_amr_curated, r.Bin, r.sample]),
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
      ...dbapis.map(r => ({ ...r, Hit: r.Hit, Detail: `pident=${r.Pident} e=${r.Evalue}` })),
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

  function _renderIslandGeneMap(body, isl, colorFor) {
    const wrap = document.createElement('div');
    wrap.className = 'region-genemap';
    body.appendChild(wrap);

    const genes = isl.window_genes || [];
    const W_PER_GENE = 46, H = 70, PAD = 10;
    const W = Math.max(200, genes.length * W_PER_GENE + PAD * 2);
    const dark = document.documentElement.dataset.theme === 'dark';

    const svg = d3.select(wrap).append('svg').attr('width', W).attr('height', H);
    const tip = d3.select(wrap).append('div').attr('class', 'd3-tooltip').style('display', 'none');

    const g = svg.append('g').attr('transform', `translate(${PAD},${H / 2})`);
    genes.forEach((gene, idx) => {
      const x = idx * W_PER_GENE;
      const arrow = g.append('path')
        .attr('d', `M${x},-10 L${x + W_PER_GENE - 14},-10 L${x + W_PER_GENE - 4},0 L${x + W_PER_GENE - 14},10 L${x},10 Z`)
        .attr('fill', gene.System ? colorFor(gene.System) : (dark ? '#334155' : '#cbd5e1'))
        .attr('fill-opacity', gene.System ? 0.9 : 0.6)
        .attr('stroke', dark ? '#0f172a' : '#fff')
        .attr('stroke-width', 1)
        .style('cursor', 'pointer');
      arrow.on('mouseover', (e) => {
        tip.style('display', 'block').style('left', (e.offsetX + 12) + 'px').style('top', (e.offsetY) + 'px')
          .html(`<strong>${gene.Protein}</strong>${gene.System ? `<br>System: ${gene.System}` : '<br>No defense system hit'}`);
      }).on('mouseout', () => tip.style('display', 'none'));
    });

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
    addVirAnti(typeof ANTIDEFENSE_VIRAL_DBAPIS !== 'undefined' ? ANTIDEFENSE_VIRAL_DBAPIS : [], 'dbAPIS', 'Hit');

    const isAll = sample === '__all__';
    const rows = (isAll ? all : all.filter(r => r.sample === sample)).map(r => {
      const ann = bins[`${r.sample}::${r.Host}`] || {};
      const va  = viralAnti.get(`${r.sample}::${r.Virus}`) || { DefenseFinder: new Set(), dbAPIS: new Set() };
      return {
        ...r,
        Host_taxonomy_display: r.Host_genus || r.Host_species || r.Host_taxonomy || '—',
        Defense_display:       (ann.Defense_systems || []).join(', ') || '—',
        Antidefense_display:   (ann.Antidefense_systems || []).join(', ') || '—',
        AMR_curated_display:   (ann.AMR_curated || []).join(', ') || '—',
        AMR_exploratory_display: (ann.AMR_exploratory || []).join(', ') || '—',
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
      { key: 'AMR_curated_display',   label: 'AMR (curated)' },
      { key: 'AMR_exploratory_display', label: 'AMR (exploratory)' },
    ], { searchId: 'hostdefense-matrix-search' });
  }

})();
