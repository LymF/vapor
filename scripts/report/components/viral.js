/* viral.js — Detection, binning, CheckV, taxonomy, vConTACT3 network, lifestyle */
(function () {
  'use strict';

  window.renderViral = function () {
    const samples = typeof SAMPLES !== 'undefined' ? SAMPLES : [];

    _renderDetection(samples);
    _renderBinning(samples);
    _renderLifestyle(samples);
    makeSampleDropdown('sample-sel-viral-tax', _renderTaxonomy, { allSamples: true });
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

    // CheckV quality donuts — "All samples" (aggregated) or a single sample
    const ALL = 'All samples';
    function _checkvCounts(data, sampleList) {
      const counts = {};
      sampleList.forEach(s => (data[s] || []).forEach(r => {
        const q = r.checkv_quality || 'Not-determined';
        counts[q] = (counts[q] || 0) + 1;
      }));
      return counts;
    }
    function _renderCheckvDonuts(sel) {
      const sampleList = sel === ALL ? samples : [sel];
      const label = sel === ALL ? 'all samples' : sel;

      const counts = _checkvCounts(cv, sampleList);
      mkChart('vir-checkv-pie-chart', {
        title: { text: `CheckV — Consensus Contigs Quality (${label})` },
        tooltip: { trigger: 'item', formatter: p => `${p.name}: ${p.value} (${p.percent.toFixed(1)}%)` },
        legend: { orient: 'vertical', left: 'left', top: 'middle' },
        series: [{
          type: 'pie', radius: ['40%', '70%'],
          data: cvCats.map((c, i) => ({ value: counts[c] || 0, name: c, itemStyle: { color: cvCols[i] } })),
          label: { formatter: '{b}: {c}' },
        }],
      });

      const countsVrh = _checkvCounts(cvr, sampleList);
      mkChart('vir-checkv-pie-vrh-chart', {
        title: { text: `CheckV — vRhyme vMAGs Quality (${label})` },
        tooltip: { trigger: 'item', formatter: p => `${p.name}: ${p.value} (${p.percent.toFixed(1)}%)` },
        legend: { orient: 'vertical', left: 'left', top: 'middle' },
        series: [{
          type: 'pie', radius: ['40%', '70%'],
          data: cvCats.map((c, i) => ({ value: countsVrh[c] || 0, name: c, itemStyle: { color: cvCols[i] } })),
          label: { formatter: '{b}: {c}' },
        }],
      });
    }

    const checkvSel = document.getElementById('sample-sel-vir-checkv');
    if (checkvSel) {
      checkvSel.innerHTML = '';
      [ALL, ...samples].forEach(s => {
        const opt = document.createElement('option');
        opt.value = opt.textContent = s;
        checkvSel.appendChild(opt);
      });
      checkvSel.addEventListener('change', () => _renderCheckvDonuts(checkvSel.value));
    }
    _renderCheckvDonuts(ALL);

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

    // Viral contig length distribution — boxplot per sample
    const lenData = typeof VIRAL_LENGTHS !== 'undefined' ? VIRAL_LENGTHS : {};
    const lenBox = [];
    const lenOutliers = [];
    samples.forEach((s, i) => {
      const { box, outliers: out } = window.boxStats(lenData[s] || []);
      lenBox.push(box);
      out.forEach(v => lenOutliers.push([i, v]));
    });

    mkChart('vir-len-chart', {
      title: { text: 'Viral Contig Length Distribution (bp)' },
      tooltip: { trigger: 'item' },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 45 }, boundaryGap: true },
      yAxis: { type: 'log', name: 'Length (bp)' },
      series: [
        { name: 'Length',  type: 'boxplot', data: lenBox,
          itemStyle: { color: PAL[0], borderColor: PAL[1] } },
        { name: 'Outlier', type: 'scatter', data: lenOutliers,
          symbolSize: 4, itemStyle: { color: PAL[3], opacity: 0.5 } },
      ],
      grid: { bottom: 90 },
    });

    // Viral contig depth distribution — boxplot per sample
    const depthData = typeof VIRAL_DEPTH !== 'undefined' ? VIRAL_DEPTH : {};
    const depthBox = [];
    const depthOutliers = [];
    samples.forEach((s, i) => {
      const { box, outliers: out } = window.boxStats(depthData[s] || []);
      depthBox.push(box);
      out.forEach(v => depthOutliers.push([i, v]));
    });

    mkChart('vir-depth-chart', {
      title: { text: 'Viral Contig Coverage Depth Distribution (×)' },
      tooltip: { trigger: 'item' },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 45 }, boundaryGap: true },
      yAxis: { type: 'log', name: 'Depth (×)' },
      series: [
        { name: 'Depth',   type: 'boxplot', data: depthBox,
          itemStyle: { color: PAL[0], borderColor: PAL[1] } },
        { name: 'Outlier', type: 'scatter', data: depthOutliers,
          symbolSize: 4, itemStyle: { color: PAL[3], opacity: 0.5 } },
      ],
      grid: { bottom: 90 },
    });
  }

  // ── Taxonomy ──────────────────────────────────────────────────────────────
  let _currentTaxSample = null;
  let _currentTax = [];

  function _renderTaxonomy(sample) {
    const allTax = typeof TAX_DATA !== 'undefined' ? TAX_DATA : [];
    const isAll  = sample === '__all__';
    const label  = isAll ? 'All samples' : sample;
    const tax    = isAll ? allTax : allTax.filter(r => r.sample === sample);
    _currentTaxSample = sample;
    _currentTax = tax;

    // Source pie — covers ALL contigs (incl. diamond_custom-only hits and unknown)
    const allSrcDist = typeof VIRAL_SOURCE_DIST !== 'undefined' ? VIRAL_SOURCE_DIST : {};
    let srcDist;
    if (isAll) {
      srcDist = {};
      Object.values(allSrcDist).forEach(d => Object.entries(d).forEach(([k, v]) => { srcDist[k] = (srcDist[k] || 0) + v; }));
    } else {
      srcDist = allSrcDist[sample] || {};
    }
    const srcPieData = Object.entries(srcDist).map(([name, value]) => ({ name, value }));
    mkChart('vir-tax-source-chart', {
      title: { text: `${label} — Classification Source` },
      tooltip: { trigger: 'item', formatter: '{b}: {c} ({d}%)' },
      series: [{ type: 'pie', radius: ['35%', '60%'], data: srcPieData, label: { formatter: '{b}\n{d}%' } }],
    });

    // Family/Genus bar (top 20) — only contigs with an actual assignment at that rank
    function renderRankBar(level) {
      const fields    = level === 'genus' ? ['final_genus', 'Genus'] : ['final_family', 'Family'];
      const rankLabel = level === 'genus' ? 'Genera' : 'Families';
      const count = {};
      tax.forEach(r => {
        const v = r[fields[0]] || r[fields[1]] || '';
        if (!v) return;
        count[v] = (count[v] || 0) + 1;
      });
      const top = Object.entries(count).sort((a, b) => b[1] - a[1]).slice(0, 20);
      mkChart('vir-tax-family-chart', {
        title: { text: `${label} — Top Viral ${rankLabel}` },
        tooltip: { trigger: 'axis' },
        xAxis: { type: 'value', name: 'Count', nameLocation: 'middle', nameGap: 28 },
        yAxis: { type: 'category', data: top.map(x => x[0]).reverse(), axisLabel: { width: 140, overflow: 'truncate' } },
        series: [{ type: 'bar', data: top.map(x => x[1]).reverse(), itemStyle: { color: '#0d9488' } }],
        grid: { left: 160, right: 30, bottom: 50 },
      });
    }

    const famBtn   = document.getElementById('vir-tax-family-btn');
    const genusBtn = document.getElementById('vir-tax-genus-btn');
    const currentLevel = (famBtn && famBtn.classList.contains('active')) || !genusBtn ? 'family'
      : (genusBtn.classList.contains('active') ? 'genus' : 'family');
    if (famBtn && genusBtn) {
      famBtn.onclick = () => { famBtn.classList.add('active'); genusBtn.classList.remove('active'); renderRankBar('family'); };
      genusBtn.onclick = () => { genusBtn.classList.add('active'); famBtn.classList.remove('active'); renderRankBar('genus'); };
      if (!famBtn.classList.contains('active') && !genusBtn.classList.contains('active')) famBtn.classList.add('active');
    }
    renderRankBar(currentLevel);

    // Taxonomy network — disabled for "All samples" (too dense)
    if (isAll) {
      const netEl = document.getElementById('vir-tax-network');
      if (netEl) netEl.innerHTML = '<p style="color:var(--text-muted);padding:1rem">Select a specific sample to view the taxonomy network.</p>';
    } else {
      _renderTaxNetwork(tax, sample);
    }

    // Table — add Sample column when showing all
    const taxCols = isAll
      ? [
          { key: 'sample',       label: 'Sample' },
          { key: 'Genome',       label: 'Contig' },
          { key: 'final_family', label: 'Family' },
          { key: 'final_genus',  label: 'Genus' },
          { key: 'Source',       label: 'Source' },
          { key: 'CheckV_quality', label: 'CheckV' },
          { key: 'Completeness', label: 'Completeness' },
        ]
      : [
          { key: 'Genome',       label: 'Contig' },
          { key: 'final_family', label: 'Family' },
          { key: 'final_genus',  label: 'Genus' },
          { key: 'Source',       label: 'Source' },
          { key: 'CheckV_quality', label: 'CheckV' },
          { key: 'Completeness', label: 'Completeness' },
        ];
    makeTable('vir-tax-table', tax, taxCols, {
      searchId: 'vir-tax-search',
      format: { CheckV_quality: qualBadge },
    });
  }

  // ── Taxonomy network (D3 force) — Order → Family → Genus → Sequence ────────
  function _renderTaxNetwork(tax, sample) {
    const el = document.getElementById('vir-tax-network');
    if (!el) return;
    el.innerHTML = '';

    if (!tax.length) {
      el.innerHTML = '<p style="color:var(--text-muted);padding:1rem">No taxonomy data for this sample.</p>';
      return;
    }

    // Build node/edge sets: each sequence links up through genus → family → order
    const nodesMap = new Map();
    const edgeKeys = new Set();
    const edges = [];

    function addNode(id, type, label) {
      let n = nodesMap.get(id);
      if (!n) { n = { id, type, label, count: 0 }; nodesMap.set(id, n); }
      n.count++;
      return n;
    }
    function addEdge(a, b) {
      const key = a < b ? `${a}|${b}` : `${b}|${a}`;
      if (edgeKeys.has(key)) return;
      edgeKeys.add(key);
      edges.push({ source: a, target: b });
    }

    tax.forEach(r => {
      const order  = r.final_order  || '';
      const family = r.final_family || '';
      const genus  = r.final_genus  || '';
      const seqId  = `seq:${r.Genome}`;
      addNode(seqId, 'sequence', r.Genome);

      let parent = seqId;
      if (genus)  { const id = `genus:${genus}`;   addNode(id, 'genus', genus);   addEdge(parent, id); parent = id; }
      if (family) { const id = `family:${family}`; addNode(id, 'family', family); addEdge(parent, id); parent = id; }
      if (order)  { const id = `order:${order}`;   addNode(id, 'order', order);   addEdge(parent, id); parent = id; }
    });

    const nodes = [...nodesMap.values()];
    if (!nodes.length) {
      el.innerHTML = '<p style="color:var(--text-muted);padding:1rem">No taxonomy data for this sample.</p>';
      return;
    }

    const dark  = document.documentElement.dataset.theme === 'dark';
    const W     = el.clientWidth || 900;
    const H     = 500;

    const typeColor = {
      sequence: dark ? '#475569' : '#cbd5e1',
      genus:    '#d97706',
      family:   '#0891b2',
      order:    '#7c3aed',
    };
    function radius(d) {
      if (d.type === 'sequence') return 3;
      if (d.type === 'genus')    return 5 + Math.sqrt(d.count) * 1.5;
      if (d.type === 'family')   return 7 + Math.sqrt(d.count) * 2;
      return 9 + Math.sqrt(d.count) * 2.5; // order
    }

    const svg = d3.select(el).append('svg').attr('width', W).attr('height', H);
    const zoomLayer = svg.append('g').attr('class', 'zoom-layer');

    const linkSel = zoomLayer.append('g').selectAll('line').data(edges).join('line')
      .attr('stroke', dark ? '#334155' : '#cbd5e1')
      .attr('stroke-width', 0.8)
      .attr('stroke-opacity', 0.6);

    const nodeSel = zoomLayer.append('g').selectAll('circle').data(nodes).join('circle')
      .attr('r', radius)
      .attr('fill', d => typeColor[d.type])
      .attr('fill-opacity', d => d.type === 'sequence' ? 0.6 : 0.9)
      .attr('stroke', dark ? '#0f172a' : '#fff')
      .attr('stroke-width', d => d.type === 'sequence' ? 0.5 : 1.5)
      .style('cursor', 'pointer');

    const labels = zoomLayer.append('g').selectAll('text')
      .data(nodes.filter(d => d.type !== 'sequence')).join('text')
      .attr('font-size', d => d.type === 'order' ? 11 : 10)
      .attr('font-weight', d => d.type === 'order' ? 600 : 400)
      .attr('fill', dark ? '#e2e8f0' : '#334155')
      .attr('dy', '-.7em')
      .style('pointer-events', 'none')
      .text(d => d.label);

    // Invisible larger "hit area" so small sequence nodes are easy to hover
    const hitSel = zoomLayer.append('g').selectAll('circle.hit').data(nodes).join('circle')
      .attr('class', 'hit')
      .attr('r', d => Math.max(radius(d) + 4, 6))
      .attr('fill', 'transparent')
      .style('cursor', 'pointer');

    // Tooltip
    const tip = d3.select(el).append('div').attr('class', 'd3-tooltip').style('display', 'none');
    const rankLabel = { sequence: 'Sequence', genus: 'Genus', family: 'Family', order: 'Order' };
    hitSel.on('mouseover', (e, d) => {
      tip.style('display', 'block').style('left', (e.offsetX + 12) + 'px').style('top', (e.offsetY - 12) + 'px')
        .html(`<strong>${rankLabel[d.type]}:</strong> ${d.label}${d.type !== 'sequence' ? `<br>Sequences: ${d.count}` : ''}`);
    }).on('mouseout', () => tip.style('display', 'none'));

    const sim = d3.forceSimulation(nodes)
      .force('link',    d3.forceLink(edges).id(d => d.id).distance(d =>
        d.source.type === 'sequence' || d.target.type === 'sequence' ? 25 : 60).strength(0.6))
      .force('charge',  d3.forceManyBody().strength(-50))
      .force('center',  d3.forceCenter(W / 2, H / 2))
      .force('collide', d3.forceCollide(d => radius(d) + 2));

    // Zoom / pan
    const zoom = d3.zoom()
      .scaleExtent([0.1, 8])
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
      hitSel.attr('cx', d => d.x).attr('cy', d => d.y);
      labels.attr('x', d => d.x).attr('y', d => d.y);
    });
    sim.on('end', fitToView);

    // Drag
    hitSel.call(d3.drag()
      .on('start', (e, d) => { if (!e.active) sim.alphaTarget(0.3).restart(); d.fx = d.x; d.fy = d.y; })
      .on('drag',  (e, d) => { d.fx = e.x; d.fy = e.y; })
      .on('end',   (e, d) => { if (!e.active) sim.alphaTarget(0); d.fx = null; d.fy = null; })
    );

    // Zoom control buttons
    const zoomIn  = document.getElementById('vir-tax-net-zoom-in');
    const zoomOut = document.getElementById('vir-tax-net-zoom-out');
    const zoomFit = document.getElementById('vir-tax-net-zoom-fit');
    if (zoomIn)  zoomIn.onclick  = () => svg.transition().duration(250).call(zoom.scaleBy, 1.4);
    if (zoomOut) zoomOut.onclick = () => svg.transition().duration(250).call(zoom.scaleBy, 1 / 1.4);
    if (zoomFit) zoomFit.onclick = fitToView;
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

    // All zoomable/pannable content lives inside this group
    const zoomLayer = svg.append('g').attr('class', 'zoom-layer');

    const linkSel = zoomLayer.append('g').selectAll('line').data(edges).join('line')
      .attr('stroke', dark ? '#334155' : '#cbd5e1')
      .attr('stroke-width', 0.8)
      .attr('stroke-opacity', 0.6);

    const nodeSel = zoomLayer.append('g').selectAll('circle').data(nodes).join('circle')
      .attr('r', d => d.is_novel === false ? 5 : 7)
      .attr('fill', d => color(d.cluster || 'Singleton'))
      .attr('fill-opacity', d => d.is_novel === false ? 0.4 : 0.85)
      .attr('stroke', dark ? '#0f172a' : '#fff')
      .attr('stroke-width', 1.5)
      .style('cursor', 'pointer');

    const labels = zoomLayer.append('g').selectAll('text').data(nodes.filter(d => d.is_novel !== false)).join('text')
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

    // Zoom / pan — lets nodes that drift outside the W×H viewport be reached
    const zoom = d3.zoom()
      .scaleExtent([0.1, 8])
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

    // Drag
    nodeSel.call(d3.drag()
      .on('start', (e, d) => { if (!e.active) sim.alphaTarget(0.3).restart(); d.fx = d.x; d.fy = d.y; })
      .on('drag',  (e, d) => { d.fx = e.x; d.fy = e.y; })
      .on('end',   (e, d) => { if (!e.active) sim.alphaTarget(0); d.fx = null; d.fy = null; })
    );

    // Zoom control buttons (rebound each render — old listeners discarded with old svg)
    const zoomIn  = document.getElementById('vc3-zoom-in');
    const zoomOut = document.getElementById('vc3-zoom-out');
    const zoomFit = document.getElementById('vc3-zoom-fit');
    if (zoomIn)  zoomIn.onclick  = () => svg.transition().duration(250).call(zoom.scaleBy, 1.4);
    if (zoomOut) zoomOut.onclick = () => svg.transition().duration(250).call(zoom.scaleBy, 1 / 1.4);
    if (zoomFit) zoomFit.onclick = fitToView;
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
    if (e.detail.sub === 'vir-taxonomy' && _currentTaxSample) {
      setTimeout(() => _renderTaxNetwork(_currentTax, _currentTaxSample), 50);
    }
  });

})();
