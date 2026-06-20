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
    makeSampleDropdown('sample-sel-hostpred', _renderHostPrediction, { allSamples: true });
    makeSampleDropdown('sample-sel-hostdefense-matrix', _renderMatrix, { allSamples: true });
  };

  function _byCount(rows, key, samples) {
    const counts = {};
    samples.forEach(s => { counts[s] = 0; });
    (rows || []).forEach(r => { if (r.sample in counts) counts[r.sample] += 1; });
    return samples.map(s => counts[s]);
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

  // ── Host <-> Virus <-> Defense <-> AMR matrix ────────────────────────────
  function _renderMatrix(sample) {
    const all = typeof HOST_DEFENSE_LINKS !== 'undefined' ? HOST_DEFENSE_LINKS : [];
    const isAll = sample === '__all__';
    const rows = (isAll ? all : all.filter(r => r.sample === sample)).map(r => ({
      ...r,
      Host_taxonomy_display: r.Host_genus || r.Host_species || r.Host_taxonomy || '—',
      Defense_display:       (r.Defense_systems || []).join(', ') || '—',
      Antidefense_display:   (r.Antidefense_systems || []).join(', ') || '—',
      AMR_curated_display:   (r.AMR_curated || []).join(', ') || '—',
      AMR_exploratory_display: (r.AMR_exploratory || []).join(', ') || '—',
    }));

    makeTable('hostdefense-matrix-table', rows, [
      { key: 'sample',  label: 'Sample' },
      { key: 'Virus',   label: 'Virus' },
      { key: 'Host',    label: 'Host (MAG)' },
      { key: 'Host_taxonomy_display', label: 'Host Taxonomy' },
      { key: 'Defense_display',       label: 'Defense Systems' },
      { key: 'Antidefense_display',   label: 'Anti-Defense Systems' },
      { key: 'AMR_curated_display',   label: 'AMR (curated)' },
      { key: 'AMR_exploratory_display', label: 'AMR (exploratory)' },
    ], { searchId: 'hostdefense-matrix-search' });
  }

})();
