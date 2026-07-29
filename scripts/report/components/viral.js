/* viral.js — Detection, binning, CheckV, taxonomy, lifestyle */
(function () {
  'use strict';

  window.renderViral = function () {
    const samples = typeof SAMPLES !== 'undefined' ? SAMPLES : [];

    _renderDetection(samples);
    _renderBinning(samples);
    _renderLifestyle(samples);
    makeSampleDropdown('sample-sel-viral-tax', _renderTaxonomy, { allSamples: true });
    makeSampleDropdown('sample-sel-votu-table', _renderVotuTable);
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
    // Quality tiers are an ORDINAL good->bad ladder (Complete > HQ > MQ > LQ
    // > Not-determined), and each card is really one distribution being
    // read as a composition -- a 100%-stacked horizontal bar (one row) lets
    // the reader compare segment widths directly and shares the same
    // baseline the MIMAG HQ/MQ/LQ bar (overview.js) already uses, instead
    // of asking them to eyeball two separate pie wedges/areas side by side.
    function _checkvBarOption(title, counts) {
      const total = cvCats.reduce((a, c) => a + (counts[c] || 0), 0) || 1;
      return {
        title: { text: title },
        tooltip: { trigger: 'item', formatter: p => `${p.seriesName}: ${p.value} (${(p.value / total * 100).toFixed(1)}%)` },
        legend: { top: 26 },
        grid: { top: 58, bottom: 20, left: 10, right: 10, containLabel: true },
        xAxis: { type: 'value', show: false },
        yAxis: { type: 'category', data: [''], axisLine: { show: false }, axisTick: { show: false } },
        series: cvCats.map((c, i) => ({
          name: c, type: 'bar', stack: 'q', barWidth: 46,
          itemStyle: { color: cvCols[i] },
          label: { show: true, formatter: p => p.value > 0 ? p.value : '', color: '#fff', fontSize: 11 },
          data: [counts[c] || 0],
        })),
      };
    }

    function _renderCheckvDonuts(sel) {
      const sampleList = sel === ALL ? samples : [sel];
      const label = sel === ALL ? 'all samples' : sel;

      mkChart('vir-checkv-pie-chart', _checkvBarOption(
        `CheckV — Consensus Contigs Quality (${label})`, _checkvCounts(cv, sampleList)));
      mkChart('vir-checkv-pie-vrh-chart', _checkvBarOption(
        `CheckV — vRhyme vMAGs Quality (${label})`, _checkvCounts(cvr, sampleList)));
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

    // Viral contig length / depth per sample. Both are heavily right-skewed, so
    // they are drawn on a log axis; distPlot picks strip vs density vs ridgeline
    // from the data (docs/REPORT_VIZ_GUIDE.md §4).
    const lenData = typeof VIRAL_LENGTHS !== 'undefined' ? VIRAL_LENGTHS : {};
    mkChart('vir-len-chart', distPlot({
      groups: samples.map(s => ({ name: s, values: (lenData[s] || []).filter(v => v > 0) })),
      title: 'Viral Contig Length Distribution (bp)',
      xName: 'Length (bp)',
      log: true,
    }));

    const depthData = typeof VIRAL_DEPTH !== 'undefined' ? VIRAL_DEPTH : {};
    mkChart('vir-depth-chart', distPlot({
      groups: samples.map(s => ({ name: s, values: (depthData[s] || []).filter(v => v > 0) })),
      title: 'Viral Contig Coverage Depth Distribution (×)',
      xName: 'Depth (×)',
      log: true,
    }));
  }

  // ── Taxonomy ──────────────────────────────────────────────────────────────
  function _renderTaxonomy(sample) {
    const allTax = typeof TAX_DATA !== 'undefined' ? TAX_DATA : [];
    const isAll  = sample === '__all__';
    const label  = isAll ? 'All samples' : sample;
    const tax    = isAll ? allTax : allTax.filter(r => r.sample === sample);

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

    // Order/Family/Genus bar (top 20) — only contigs with an assignment at that rank.
    // Order is included because GeNomad-only fallback hits (no other tool found a
    // match) typically only resolve to order/class level, never family or genus —
    // without this level they'd silently disappear from every rank chart despite
    // being a large share of "Source".
    const RANK_FIELD = { order: 'final_order', family: 'final_family', genus: 'final_genus' };
    const RANK_LABEL = { order: 'Orders', family: 'Families', genus: 'Genera' };
    function renderRankBar(level) {
      const field = RANK_FIELD[level] || 'final_family';
      const count = {};
      tax.forEach(r => {
        const v = r[field] || '';
        if (!v) return;
        count[v] = (count[v] || 0) + 1;
      });
      const top = Object.entries(count).sort((a, b) => b[1] - a[1]).slice(0, 20);
      mkChart('vir-tax-family-chart', {
        title: { text: `${label} — Top Viral ${RANK_LABEL[level] || 'Families'}` },
        tooltip: { trigger: 'axis' },
        xAxis: { type: 'value', name: 'Count', nameLocation: 'middle', nameGap: 28 },
        yAxis: { type: 'category', data: top.map(x => x[0]).reverse(), axisLabel: { width: 140, overflow: 'truncate' } },
        series: [{ type: 'bar', data: top.map(x => x[1]).reverse(), itemStyle: { color: '#0d9488' } }],
        grid: { left: 160, right: 30, bottom: 50 },
      });
    }

    const rankBtns = {
      order:  document.getElementById('vir-tax-order-btn'),
      family: document.getElementById('vir-tax-family-btn'),
      genus:  document.getElementById('vir-tax-genus-btn'),
    };
    const activeLevel = Object.entries(rankBtns).find(([, btn]) => btn?.classList.contains('active'))?.[0] || 'family';
    Object.entries(rankBtns).forEach(([level, btn]) => {
      if (!btn) return;
      btn.onclick = () => {
        Object.values(rankBtns).forEach(b => b && b.classList.remove('active'));
        btn.classList.add('active');
        renderRankBar(level);
      };
    });
    if (rankBtns.family && !Object.values(rankBtns).some(b => b?.classList.contains('active'))) {
      rankBtns.family.classList.add('active');
    }
    renderRankBar(activeLevel);

    // Table — hidden when "All samples" (too many rows)
    const tableCard = document.querySelector('#vir-tax-table')?.closest('.chart-card');
    if (tableCard) tableCard.style.display = isAll ? 'none' : '';
    if (!isAll) {
      makeTable('vir-tax-table', tax, [
        { key: 'Genome',         label: 'Contig' },
        { key: 'final_order',    label: 'Order' },
        { key: 'final_family',   label: 'Family' },
        { key: 'final_genus',    label: 'Genus' },
        { key: 'Source',         label: 'Source' },
        { key: 'CheckV_quality', label: 'CheckV' },
        { key: 'Completeness',   label: 'Completeness' },
      ], {
        searchId: 'vir-tax-search',
        format: { CheckV_quality: qualBadge },
      });
    }
  }

  // ── Lifestyle ─────────────────────────────────────────────────────────────
  function _renderLifestyle(samples) {
    const ls = typeof LIFESTYLE !== 'undefined' ? LIFESTYLE : {};
    const labels = ['Lytic', 'Lysogenic', 'Unknown'];
    // Lytic/lysogenic is a categorical biological classification, not a
    // good/bad status axis -- red (PAL[7]) means "low quality/problem"
    // everywhere else in this report (CheckV, CheckM2, MIMAG), so using it
    // here would mislead a reader skimming tabs into reading "Lytic" as the
    // concerning outcome. PAL[3]/PAL[2] are plain categorical identity;
    // PAL_MUTED for the true "don't know" bucket.
    const cols = [PAL[3], PAL[2], PAL_MUTED];

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

  // ── vOTU Table ───────────────────────────────────────────────────────────
  function _renderVotuTable(sample) {
    const all  = typeof VOTU_DATA !== 'undefined' ? VOTU_DATA : {};
    const rows = (all[sample] || []).filter(r => String(r.is_rep).toLowerCase() === 'true');
    makeTable('votu-table', rows, [
      { key: 'representative',      label: 'Representative' },
      { key: 'cluster_size',        label: 'Cluster' },
      { key: 'rep_length_bp',       label: 'Length (bp)' },
      { key: 'checkv_quality',      label: 'CheckV' },
      { key: 'checkv_completeness', label: 'Completeness' },
      { key: 'lifestyle',           label: 'Lifestyle' },
      { key: 'n_AMGs',              label: 'AMGs' },
      { key: 'breadth',             label: 'Breadth' },
      { key: 'rpmpm',               label: 'RPMPM' },
    ], {
      searchId: 'votu-table-search',
      format: { checkv_quality: qualBadge },
    });
  }

})();
