/* viral.js — Detection, binning, CheckV, taxonomy, vConTACT3 network, lifestyle */
(function () {
  'use strict';

  window.renderViral = function () {
    const samples = typeof SAMPLES !== 'undefined' ? SAMPLES : [];

    _renderDetection(samples);
    _renderBinning(samples);
    _renderLifestyle(samples);
    makeSampleDropdown('sample-sel-viral-tax', _renderTaxonomy);
    makeSampleDropdown('sample-sel-vc3', window.renderVC3Network);
  };

  // ── Detection ─────────────────────────────────────────────────────────────
  function _renderDetection(samples) {
    const vt  = typeof VIRAL_TOOLS !== 'undefined' ? VIRAL_TOOLS : {};
    const sup = typeof SUPPORT     !== 'undefined' ? SUPPORT     : {};

    const tools = ['VirSorter2', 'GeNomad', 'VIBRANT'];
    mkChart('vir-tools-chart', {
      title: { text: 'Viral Contigs per Tool' },
      tooltip: { trigger: 'axis' },
      legend: { data: tools },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: 'Contigs' },
      series: tools.map((t, i) => ({
        name: t, type: 'bar',
        data: samples.map(s => (vt[s] || {})[t] || 0),
        color: PAL[i % PAL.length],
      })),
      grid: { bottom: 70 },
    });

    mkChart('vir-consensus-chart', {
      title: { text: 'Tool Agreement — Consensus Filter' },
      tooltip: { trigger: 'axis' },
      legend: { data: ['1 tool', '2 tools', '≥3 tools'] },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: 'Contigs' },
      series: [
        { name: '1 tool',   type: 'bar', stack: 'c', data: samples.map(s => (sup[s] || {})[1] || 0), color: '#64748b' },
        { name: '2 tools',  type: 'bar', stack: 'c', data: samples.map(s => (sup[s] || {})[2] || 0), color: '#d97706' },
        { name: '≥3 tools', type: 'bar', stack: 'c', data: samples.map(s => (sup[s] || {})[3] || 0), color: '#0d9488' },
      ],
      grid: { bottom: 70 },
    });

    // Tool support heatmap (sample × tool, value = count)
    const heatData = [];
    samples.forEach((s, si) => {
      tools.forEach((t, ti) => {
        heatData.push([ti, si, (vt[s] || {})[t] || 0]);
      });
    });
    mkChart('vir-heatmap-chart', {
      title: { text: 'Tool Detection Matrix' },
      tooltip: {
        formatter: p => `${samples[p.data[1]]} / ${tools[p.data[0]]}: ${p.data[2]} contigs`,
      },
      xAxis: { type: 'category', data: tools },
      yAxis: { type: 'category', data: samples },
      visualMap: {
        min: 0, max: Math.max(1, ...heatData.map(d => d[2])),
        inRange: { color: ['#f0fdfa', '#0d9488', '#0f172a'] },
        show: true, orient: 'horizontal', left: 'right', bottom: 10,
        textStyle: { color: 'inherit' },
      },
      series: [{
        type: 'heatmap',
        data: heatData,
        label: { show: true, formatter: p => p.data[2] },
      }],
      grid: { bottom: 60, top: 40 },
    });
  }

  // ── Binning & Quality ─────────────────────────────────────────────────────
  function _renderBinning(samples) {
    const cv  = typeof CHECKV     !== 'undefined' ? CHECKV     : {};
    const cvr = typeof CHECKV_VRH !== 'undefined' ? CHECKV_VRH : {};
    const vrh = typeof VRHYME     !== 'undefined' ? VRHYME     : {};
    const sup = typeof SUPPORT    !== 'undefined' ? SUPPORT    : {};

    const cvCats = ['Complete', 'High-quality', 'Medium-quality', 'Low-quality', 'Not-determined'];
    const cvCols = ['#16a34a', '#4ade80', '#fbbf24', '#d97706', '#ef4444'];

    // CheckV pie (consensus)
    const pieData = samples.map(s => {
      const counts = {};
      (cv[s] || []).forEach(r => {
        const q = r.checkv_quality || 'Not-determined';
        counts[q] = (counts[q] || 0) + 1;
      });
      return cvCats.map((c, i) => ({ value: counts[c] || 0, name: c, itemStyle: { color: cvCols[i] } }));
    });

    // Multi-pie for consensus
    _multiPie('vir-checkv-pie-chart', 'CheckV — Consensus Contigs Quality', samples, pieData);

    // CheckV pie (vRhyme)
    const pieDataVrh = samples.map(s => {
      const counts = {};
      (cvr[s] || []).forEach(r => {
        const q = r.checkv_quality || 'Not-determined';
        counts[q] = (counts[q] || 0) + 1;
      });
      return cvCats.map((c, i) => ({ value: counts[c] || 0, name: c, itemStyle: { color: cvCols[i] } }));
    });
    _multiPie('vir-checkv-pie-vrh-chart', 'CheckV — vRhyme vMAGs Quality', samples, pieDataVrh);

    // CheckV scatter: length vs completeness
    const scatterSeries = [];
    const qualColorMap = { 'Complete': '#16a34a', 'High-quality': '#4ade80', 'Medium-quality': '#fbbf24', 'Low-quality': '#d97706', 'Not-determined': '#64748b' };
    const seenQual = new Set();
    samples.forEach(s => {
      ['consensus', 'vrhyme'].forEach(src => {
        const rows = src === 'consensus' ? (cv[s] || []) : (cvr[s] || []);
        rows.forEach(r => {
          const len  = +(r.contig_length || 0);
          const comp = +(r.completeness  || 0);
          const qual = r.checkv_quality || 'Not-determined';
          if (!len || !comp) return;
          const key  = `${qual}|${src}`;
          const sym  = src === 'consensus' ? 'circle' : 'diamond';
          if (!seenQual.has(key)) {
            seenQual.add(key);
            scatterSeries.push({ name: `${qual} (${src})`, type: 'scatter', symbolSize: 8,
              symbol: sym, itemStyle: { color: qualColorMap[qual] || '#64748b', opacity: 0.75 },
              data: [] });
          }
          const series = scatterSeries.find(s2 => s2.name === `${qual} (${src})`);
          if (series) series.data.push({ value: [len, comp], name: r.contig_id || r.contig || '' });
        });
      });
    });

    mkChart('vir-checkv-scatter-chart', {
      title: { text: 'CheckV — Length vs Completeness (● consensus  ◆ vRhyme)' },
      tooltip: {
        trigger: 'item',
        formatter: p => `${p.data.name}<br>Length: ${p.data.value[0].toLocaleString()} bp<br>Completeness: ${p.data.value[1].toFixed(1)}%`,
      },
      legend: { type: 'scroll' },
      xAxis: { type: 'log', name: 'Length (bp)', min: 1000, max: 1000000, nameLocation: 'middle', nameGap: 30 },
      yAxis: { type: 'value', name: 'Completeness (%)', min: 0, max: 105 },
      series: scatterSeries,
      grid: { bottom: 60, top: 70 },
    });

    // Add quality zones (markArea + markLine) identical pattern to CheckM2
    const cvChart = window._charts['vir-checkv-scatter-chart'];
    if (cvChart) {
      cvChart.setOption({
        series: [
          ...scatterSeries,
          {
            type: 'scatter', data: [],
            markArea: {
              silent: true,
              data: [
                [{ yAxis: 90,  itemStyle: { color: 'rgba(22,163,74,0.09)' } },  { yAxis: 105 }],
                [{ yAxis: 50,  itemStyle: { color: 'rgba(217,119,6,0.07)' }  },  { yAxis: 90  }],
                [{ yAxis: 0,   itemStyle: { color: 'rgba(239,68,68,0.05)' }  },  { yAxis: 50  }],
              ],
            },
            markLine: {
              silent: true,
              data: [
                { yAxis: 90, lineStyle: { type: 'dashed', color: '#16a34a' }, label: { formatter: '≥90% HQ',  color: '#16a34a', fontSize: 10 } },
                { yAxis: 50, lineStyle: { type: 'dotted', color: '#d97706' }, label: { formatter: '≥50% MQ',  color: '#d97706', fontSize: 10 } },
              ],
            },
          },
        ],
      }, false);
    }

    // vRhyme summary
    mkChart('vir-vrhyme-chart', {
      title: { text: 'vRhyme — vMAG Summary' },
      tooltip: { trigger: 'axis' },
      legend: { data: ['Consensus input', 'vMAGs formed', 'Contigs binned'] },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: 'Count' },
      series: [
        { name: 'Consensus input', type: 'bar', data: samples.map(s => ((sup[s]||{})[3]||0) + ((sup[s]||{})[4]||0)), color: '#0891b2' },
        { name: 'vMAGs formed',    type: 'bar', data: samples.map(s => (vrh[s] || {}).n_bins || 0),         color: '#0d9488' },
        { name: 'Contigs binned',  type: 'bar', data: samples.map(s => (vrh[s] || {}).total_members || 0),  color: '#d97706' },
      ],
      grid: { bottom: 70 },
    });

    // Viral contig length distribution
    const lenData = typeof VIRAL_LENGTHS !== 'undefined' ? VIRAL_LENGTHS : {};
    const binSize = 2000; const maxLen = 60000;
    const binEdges = []; for (let i = 0; i <= maxLen; i += binSize) binEdges.push(i);

    const lenSeries = samples.map((s, i) => {
      const lens = lenData[s] || [];
      const bins = new Array(binEdges.length - 1).fill(0);
      lens.forEach(l => { const b = Math.min(Math.floor(l / binSize), bins.length - 1); bins[b]++; });
      return { name: s, type: 'bar', barMaxWidth: 18, data: bins, color: PAL[i % PAL.length], opacity: 0.75 };
    });

    mkChart('vir-len-chart', {
      title: { text: 'Viral Contig Length Distribution' },
      tooltip: { trigger: 'axis' },
      legend: { data: samples },
      xAxis: { type: 'category', data: binEdges.slice(0, -1).map(v => (v/1000).toFixed(0) + 'k'), axisLabel: { rotate: 45 } },
      yAxis: { type: 'value', name: 'Count' },
      series: lenSeries,
      grid: { bottom: 80 },
    });

    // Viral contig depth distribution
    const depthData = typeof VIRAL_DEPTH !== 'undefined' ? VIRAL_DEPTH : {};
    const dBinSize  = 20; const dMax = 500;
    const dBinEdges = []; for (let i = 0; i <= dMax; i += dBinSize) dBinEdges.push(i);

    const depthSeries = samples.map((s, i) => {
      const vals = depthData[s] || [];
      const bins = new Array(dBinEdges.length - 1).fill(0);
      vals.forEach(v => { const b = Math.min(Math.floor(v / dBinSize), bins.length - 1); bins[b]++; });
      return { name: s, type: 'bar', barMaxWidth: 18, data: bins, color: PAL[i % PAL.length], opacity: 0.75 };
    });

    mkChart('vir-depth-chart', {
      title: { text: 'Viral Contig Coverage Depth Distribution' },
      tooltip: { trigger: 'axis' },
      legend: { data: samples },
      xAxis: { type: 'category', data: dBinEdges.slice(0, -1).map(v => v + '×'), axisLabel: { rotate: 45 } },
      yAxis: { type: 'value', name: 'Count' },
      series: depthSeries,
      grid: { bottom: 80 },
    });
  }

  // ── Multi-pie helper ──────────────────────────────────────────────────────
  // Lays out one donut per sample on a grid (rows x cols) sized so they don't
  // overlap, regardless of how many samples there are.
  function _multiPie(id, title, samples, perSampleData) {
    const n    = Math.max(samples.length, 1);
    const cols = Math.max(1, Math.ceil(Math.sqrt(n)));
    const rows = Math.max(1, Math.ceil(n / cols));
    const el   = document.getElementById(id);
    const W    = (el && el.clientWidth)  || 900;
    const H    = (el && el.clientHeight) || 420;
    const minDim  = Math.min(W, H);
    const outerPx = Math.min(W / cols, H / rows) * 0.40;
    const outerPct = (outerPx / minDim * 100).toFixed(1) + '%';
    const innerPct = (outerPx * 0.5 / minDim * 100).toFixed(1) + '%';
    const dark = document.documentElement.dataset.theme === 'dark';
    const labelColor = dark ? '#94a3b8' : '#64748b';

    const pies = samples.map((s, i) => {
      const col = i % cols, row = Math.floor(i / cols);
      return {
        name: s,
        type: 'pie',
        radius: [innerPct, outerPct],
        center: [`${(col + 0.5) / cols * 100}%`, `${(row + 0.5) / rows * 100}%`],
        data: perSampleData[i],
        label: { show: false },
        tooltip: { formatter: p => `${s}<br>${p.name}: ${p.value} (${p.percent.toFixed(1)}%)` },
      };
    });

    const graphics = samples.map((s, i) => {
      const col = i % cols, row = Math.floor(i / cols);
      const label = s.length > 14 ? s.slice(0, 13) + '…' : s;
      return {
        type: 'text',
        left: `${(col + 0.5) / cols * 100}%`,
        top:  `${Math.min(98, (row + 0.5) / rows * 100 + (outerPx / H * 100) + 2)}%`,
        style: { text: label, fill: labelColor, font: '10px Inter, system-ui, sans-serif', textAlign: 'center' },
      };
    });

    mkChart(id, {
      title: { text: title },
      tooltip: { trigger: 'item' },
      series: pies,
      graphic: graphics,
    });
  }

  // ── Taxonomy ──────────────────────────────────────────────────────────────
  function _renderTaxonomy(sample) {
    const tax = typeof TAX_DATA !== 'undefined' ? TAX_DATA.filter(r => r.sample === sample) : [];

    // Source pie — covers ALL contigs (incl. diamond_custom-only hits and unknown)
    const srcDist = (typeof VIRAL_SOURCE_DIST !== 'undefined' ? VIRAL_SOURCE_DIST : {})[sample] || {};
    const srcPieData = Object.entries(srcDist).map(([name, value]) => ({ name, value }));
    mkChart('vir-tax-source-chart', {
      title: { text: `${sample} — Classification Source` },
      tooltip: { trigger: 'item', formatter: '{b}: {c} ({d}%)' },
      series: [{ type: 'pie', radius: ['35%', '60%'], data: srcPieData, label: { formatter: '{b}\n{d}%' } }],
    });

    // Family bar (top 20) — only contigs with an actual family-level assignment
    const famCount = {};
    tax.forEach(r => {
      const f = r.final_family || r.Family || '';
      if (!f) return;
      famCount[f] = (famCount[f] || 0) + 1;
    });
    const topFams = Object.entries(famCount).sort((a, b) => b[1] - a[1]).slice(0, 20);
    mkChart('vir-tax-family-chart', {
      title: { text: `${sample} — Top Viral Families` },
      tooltip: { trigger: 'axis' },
      xAxis: { type: 'value', name: 'Count' },
      yAxis: { type: 'category', data: topFams.map(x => x[0]).reverse(), axisLabel: { width: 140, overflow: 'truncate' } },
      series: [{ type: 'bar', data: topFams.map(x => x[1]).reverse(), itemStyle: { color: '#0d9488' } }],
      grid: { left: 160, right: 30 },
    });

    // Zoomable icicle (D3) for taxonomy hierarchy
    _renderTaxIcicle(tax, sample);

    // Table
    makeTable('vir-tax-table', tax, [
      { key: 'Genome',       label: 'Contig' },
      { key: 'final_family', label: 'Family' },
      { key: 'final_genus',  label: 'Genus' },
      { key: 'Source',       label: 'Source' },
      { key: 'CheckV_quality', label: 'CheckV', format: qualBadge },
      { key: 'Completeness', label: 'Completeness' },
    ], {
      searchId: 'vir-tax-search',
      format: { CheckV_quality: qualBadge },
    });
  }

  // ── Taxonomy icicle (D3 partition) ────────────────────────────────────────
  function _renderTaxIcicle(tax, sample) {
    const el = document.getElementById('vir-tax-icicle');
    if (!el) return;
    el.innerHTML = '';

    // Build hierarchy: root → family → genus
    const famMap = {};
    tax.forEach(r => {
      const fam = r.final_family || 'Unclassified';
      const gen = r.final_genus  || 'Unclassified';
      if (!famMap[fam]) famMap[fam] = {};
      famMap[fam][gen] = (famMap[fam][gen] || 0) + 1;
    });

    const root = { name: 'Viruses', children: [] };
    Object.entries(famMap).forEach(([fam, gens]) => {
      root.children.push({
        name: fam,
        children: Object.entries(gens).map(([gen, cnt]) => ({ name: gen, value: cnt })),
      });
    });

    if (!root.children.length) {
      el.innerHTML = '<p style="color:var(--text-muted);padding:1rem">No taxonomy data for this sample.</p>';
      return;
    }

    const dark    = document.documentElement.dataset.theme === 'dark';
    const W       = el.clientWidth || 900;
    const H       = 460;
    const color   = d3.scaleOrdinal(PAL);

    const hierarchy = d3.hierarchy(root).sum(d => d.value || 0).sort((a, b) => b.value - a.value);
    d3.partition().size([H, W])(hierarchy);

    const svg = d3.select(el).append('svg').attr('width', W).attr('height', H).style('font', '11px Inter, system-ui, sans-serif');
    const g   = svg.append('g');

    let _current = hierarchy;

    function clicked(event, p) {
      _current = p;
      const t = svg.transition().duration(500);
      g.selectAll('g').transition(t).attr('transform', d => {
        const x = (d.x0 - p.x0) / (p.x1 - p.x0) * H;
        const y = (d.y0 - p.y0) / (p.y1 - p.y0) * W;
        return `translate(${y},${x})`;
      });
      g.selectAll('rect').transition(t)
        .attr('height', d => Math.max(0, (d.x1 - d.x0) / (p.x1 - p.x0) * H - 1))
        .attr('width',  d => Math.max(0, (d.y1 - d.y0) / (p.y1 - p.y0) * W - 1));
      g.selectAll('text').style('display', d => {
        const visH = (d.x1 - d.x0) / (p.x1 - p.x0) * H;
        return visH > 14 ? null : 'none';
      });
    }

    const cell = g.selectAll('g').data(hierarchy.descendants()).join('g')
      .attr('transform', d => `translate(${d.y0},${d.x0})`);

    cell.append('rect')
      .attr('width',  d => Math.max(0, d.y1 - d.y0 - 1))
      .attr('height', d => Math.max(0, d.x1 - d.x0 - 1))
      .attr('fill',   d => { let node = d; while (node.depth > 1) node = node.parent; return color(node.data.name); })
      .attr('fill-opacity', d => 0.6 + d.depth * 0.13)
      .attr('rx', 2)
      .style('cursor', 'pointer')
      .on('click', clicked);

    cell.append('text')
      .attr('x', 4).attr('y', 13)
      .style('font-size', '11px')
      .style('fill', dark ? '#f1f5f9' : '#0f172a')
      .style('pointer-events', 'none')
      .style('display', d => (d.x1 - d.x0) > 14 ? null : 'none')
      .text(d => d.data.name)
      .each(function (d) {
        const w = d.y1 - d.y0 - 8;
        const self = d3.select(this);
        let text = d.data.name;
        self.text(text);
        while (this.getComputedTextLength() > w && text.length > 3) {
          text = text.slice(0, -1);
          self.text(text + '…');
        }
      });

    // Back-click on root
    svg.on('click', e => { if (e.target === svg.node()) clicked(e, hierarchy); });
  }

  // ── vConTACT3 Network (D3 force) ─────────────────────────────────────────
  window._currentVC3Sample = null;

  window.renderVC3Network = function (sample) {
    if (!sample) return;
    window._currentVC3Sample = sample;
    const el = document.getElementById('vc3-network-svg');
    if (!el) return;
    el.innerHTML = '';

    const net  = (typeof VC3_NETWORK !== 'undefined' ? VC3_NETWORK : {})[sample] || { nodes: [], edges: [] };
    const nodes = net.nodes || [];
    const edges = net.edges || [];

    if (!nodes.length) {
      el.innerHTML = '<p style="color:var(--text-muted);padding:1rem">No vConTACT3 network data for this sample.</p>';
      return;
    }

    const dark = document.documentElement.dataset.theme === 'dark';
    const W    = el.clientWidth || 900;
    const H    = 560;
    const color = d3.scaleOrdinal(PAL);

    const svg = d3.select(el).append('svg').attr('width', W).attr('height', H);
    svg.append('defs').append('marker')
      .attr('id', 'arrow').attr('viewBox', '0 -5 10 10').attr('refX', 18)
      .attr('markerWidth', 6).attr('markerHeight', 6).attr('orient', 'auto')
      .append('path').attr('d', 'M0,-5L10,0L0,5').attr('fill', dark ? '#475569' : '#94a3b8');

    const linkSel = svg.append('g').selectAll('line').data(edges).join('line')
      .attr('stroke', dark ? '#334155' : '#cbd5e1')
      .attr('stroke-width', 0.8)
      .attr('stroke-opacity', 0.6);

    const nodeSel = svg.append('g').selectAll('circle').data(nodes).join('circle')
      .attr('r', d => d.is_novel === false ? 5 : 7)
      .attr('fill', d => color(d.cluster || 'Singleton'))
      .attr('fill-opacity', d => d.is_novel === false ? 0.4 : 0.85)
      .attr('stroke', dark ? '#0f172a' : '#fff')
      .attr('stroke-width', 1.5)
      .style('cursor', 'pointer');

    const labels = svg.append('g').selectAll('text').data(nodes.filter(d => d.is_novel !== false)).join('text')
      .attr('font-size', 9)
      .attr('fill', dark ? '#94a3b8' : '#64748b')
      .attr('dy', '-.5em')
      .style('pointer-events', 'none')
      .text(d => (d.id || '').split('_').slice(-2).join('_'));

    // Tooltip
    const tip = d3.select(el).append('div').attr('class', 'd3-tooltip').style('display', 'none');
    nodeSel.on('mouseover', (e, d) => {
      tip.style('display', 'block').style('left', (e.offsetX + 12) + 'px').style('top', (e.offsetY - 12) + 'px')
        .html(`<strong>${d.id || ''}</strong><br>VC: ${d.cluster || 'Singleton'}<br>Family: ${d.family || '—'}`);
    }).on('mouseout', () => tip.style('display', 'none'));

    const sim = d3.forceSimulation(nodes)
      .force('link',    d3.forceLink(edges).id(d => d.id).distance(50).strength(0.3))
      .force('charge',  d3.forceManyBody().strength(-60))
      .force('center',  d3.forceCenter(W / 2, H / 2))
      .force('collide', d3.forceCollide(12));

    sim.on('tick', () => {
      linkSel.attr('x1', d => d.source.x).attr('y1', d => d.source.y)
             .attr('x2', d => d.target.x).attr('y2', d => d.target.y);
      nodeSel.attr('cx', d => d.x).attr('cy', d => d.y);
      labels.attr('x', d => d.x).attr('y', d => d.y);
    });

    // Drag
    nodeSel.call(d3.drag()
      .on('start', (e, d) => { if (!e.active) sim.alphaTarget(0.3).restart(); d.fx = d.x; d.fy = d.y; })
      .on('drag',  (e, d) => { d.fx = e.x; d.fy = e.y; })
      .on('end',   (e, d) => { if (!e.active) sim.alphaTarget(0); d.fx = null; d.fy = null; })
    );
  };

  // ── Lifestyle ─────────────────────────────────────────────────────────────
  function _renderLifestyle(samples) {
    const ls = typeof LIFESTYLE !== 'undefined' ? LIFESTYLE : {};
    const labels = ['Lytic', 'Lysogenic', 'Unknown'];
    const cols   = ['#ef4444', '#7c3aed', '#64748b'];

    const series = labels.map((l, i) => ({
      name: l, type: 'bar', stack: 'ls', color: cols[i],
      data: samples.map(s => (ls[s] || {})[l.toLowerCase()] || 0),
    }));

    mkChart('vir-lifestyle-chart', {
      title: { text: 'Viral Lifestyle Prediction (VIBRANT) — Lytic vs Lysogenic' },
      tooltip: { trigger: 'axis' },
      legend: { data: labels },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: 'Count' },
      series,
      grid: { bottom: 70 },
    });

    // AMG table
    const amg = typeof VIBRANT_AMG !== 'undefined' ? VIBRANT_AMG : [];
    makeTable('amg-table', amg, [
      { key: 'sample',    label: 'Sample' },
      { key: 'Pathway',   label: 'Pathway' },
      { key: 'Metabolism',label: 'Metabolism' },
      { key: 'Total_AMGs',label: 'AMGs' },
      { key: 'KOs',       label: 'KO IDs' },
    ], { searchId: 'amg-search' });
  }

  // Re-render D3 network when its sub-panel becomes visible (so clientWidth is correct)
  document.addEventListener('vapor:subtabshow', function (e) {
    if (e.detail.sub === 'vir-network' && typeof window.renderVC3Network === 'function') {
      setTimeout(() => window.renderVC3Network(window._currentVC3Sample), 50);
    }
  });

})();
