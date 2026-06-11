/* prokaryotic.js — Binning, CheckM2, GTDB-Tk taxonomy */
(function () {
  'use strict';

  window.renderProkaryotic = function () {
    const samples = typeof SAMPLES !== 'undefined' ? SAMPLES : [];
    _renderBinning(samples);
    _renderQuality(samples);
    makeSampleDropdown('sample-sel-prok-tax', _renderTaxonomy);
  };

  // ── Binning summary ───────────────────────────────────────────────────────
  function _renderBinning(samples) {
    const binner = typeof BINNER !== 'undefined' ? BINNER : {};
    const tools  = ['MetaBAT2', 'VAMB', 'SemiBin2', 'Binette (final)'];
    const colors  = ['#0891b2', '#7c3aed', '#d97706', '#0d9488'];

    mkChart('prok-binner-chart', {
      title: { text: 'Total Bins per Tool' },
      tooltip: { trigger: 'axis' },
      legend: { data: tools },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: 'Bins' },
      series: tools.map((t, i) => ({
        name: t, type: 'bar', color: colors[i],
        data: samples.map(s => (binner[s] || {})[t]?.total || 0),
      })),
      grid: { bottom: 70 },
    });

    mkChart('prok-domain-chart', {
      title: { text: 'Binette — Domain Composition (GTDB-Tk)' },
      tooltip: { trigger: 'axis' },
      legend: { data: ['Bacteria', 'Archaea', 'Unknown'] },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: 'Bins' },
      series: [
        { name: 'Bacteria', type: 'bar', stack: 'd', color: '#0d9488', data: samples.map(s => (binner[s] || {})['Binette (final)']?.bacteria || 0) },
        { name: 'Archaea',  type: 'bar', stack: 'd', color: '#d97706', data: samples.map(s => (binner[s] || {})['Binette (final)']?.archaea  || 0) },
        { name: 'Unknown',  type: 'bar', stack: 'd', color: '#64748b', data: samples.map(s => (binner[s] || {})['Binette (final)']?.unknown  || 0) },
      ],
      grid: { bottom: 70 },
    });

    // Genome size distribution — boxplot per sample
    const cm = typeof CHECKM2 !== 'undefined' ? CHECKM2 : {};
    const sizeBox = [];
    const sizeOutliers = [];
    samples.forEach((s, i) => {
      const sizes = (cm[s] || []).map(r => +(r.Genome_Size || r.genome_size || 0) / 1e6).filter(v => v > 0);
      const { box, outliers: out } = window.boxStats(sizes);
      sizeBox.push(box);
      out.forEach(v => sizeOutliers.push([i, v]));
    });

    mkChart('prok-size-chart', {
      title: { text: 'MAG Genome Size Distribution (Mb)' },
      tooltip: { trigger: 'item' },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 45 }, boundaryGap: true },
      yAxis: { type: 'value', name: 'Genome size (Mb)' },
      series: [
        { name: 'Genome size', type: 'boxplot', data: sizeBox,
          itemStyle: { color: PAL[0], borderColor: PAL[1] } },
        { name: 'Outlier', type: 'scatter', data: sizeOutliers,
          symbolSize: 4, itemStyle: { color: PAL[3], opacity: 0.5 } },
      ],
      grid: { bottom: 90 },
    });
  }

  // ── Quality (CheckM2) ─────────────────────────────────────────────────────
  function _renderQuality(samples) {
    const cm = typeof CHECKM2 !== 'undefined' ? CHECKM2 : {};

    // Scatter: Contamination vs Completeness + MIMAG zones
    const domColors = { Bacteria: '#0d9488', Archaea: '#d97706', Unknown: '#64748b' };
    const merged = typeof MERGED_PROK !== 'undefined' ? MERGED_PROK : [];
    const domLookup = {};
    merged.forEach(r => { domLookup[`${r.sample}::${r.Bin}`] = r.Domain || 'Unknown'; });
    const seriesMap = {};
    samples.forEach(s => {
      (cm[s] || []).forEach(r => {
        const comp = +(r.Completeness || 0);
        const cont = +(r.Contamination || 0);
        const name = (r.Name || r.name || '').replace(/\.fa$/, '');
        const dom  = domLookup[`${s}::${name}`] || 'Unknown';
        if (!seriesMap[dom]) {
          seriesMap[dom] = { name: dom, type: 'scatter', symbolSize: 9, color: domColors[dom],
            itemStyle: { opacity: 0.8, borderColor: 'white', borderWidth: 1 }, data: [] };
        }
        seriesMap[dom].data.push({ value: [cont, comp], name: `${name} (${s})` });
      });
    });

    mkChart('prok-checkm2-chart', {
      title: { text: 'MAG Quality — Completeness vs Contamination (CheckM2)' },
      tooltip: { formatter: p => `${p.data.name}<br>Comp: ${p.data.value[1].toFixed(1)}%<br>Cont: ${p.data.value[0].toFixed(1)}%` },
      legend: { data: Object.keys(seriesMap) },
      xAxis: { type: 'value', name: 'Contamination (%)', min: -0.5, max: 15, nameLocation: 'middle', nameGap: 30 },
      yAxis: { type: 'value', name: 'Completeness (%)',  min: -2,   max: 105, nameLocation: 'middle', nameGap: 35 },
      series: Object.values(seriesMap),
      // Colored background zones via markArea on the first series
    });

    // Add MIMAG zone shapes via graphic (ECharts markArea on a placeholder series)
    const chart = window._charts['prok-checkm2-chart'];
    if (chart) {
      chart.setOption({
        series: [
          ...Object.values(seriesMap),
          {
            type: 'scatter', data: [], markArea: {
              silent: true, data: [
                [{ coord: [-0.5, 90], itemStyle: { color: 'rgba(22,163,74,0.08)' } }, { coord: [5, 106] }],
                [{ coord: [-0.5, 50], itemStyle: { color: 'rgba(217,119,6,0.07)' } }, { coord: [10, 90] }],
                [{ coord: [-0.5, -3], itemStyle: { color: 'rgba(239,68,68,0.05)' } }, { coord: [15, 50] }],
              ],
            },
            markLine: {
              silent: true,
              data: [
                { yAxis: 90, lineStyle: { type: 'dashed', color: '#16a34a' }, label: { formatter: '≥90% HQ' } },
                { yAxis: 50, lineStyle: { type: 'dotted', color: '#d97706' }, label: { formatter: '≥50% MQ' } },
                { xAxis: 5,  lineStyle: { type: 'dashed', color: '#ef4444' }, label: { formatter: '≤5% cont' } },
                { xAxis: 10, lineStyle: { type: 'dotted', color: '#d97706' }, label: { formatter: '≤10% MQ' } },
              ],
            },
          },
        ],
      }, false);
    }

    // Completeness distribution — boxplot per sample
    const compBox = [];
    const compOutliers = [];
    samples.forEach((s, i) => {
      const vals = (cm[s] || []).map(r => +(r.Completeness || 0));
      const { box, outliers: out } = window.boxStats(vals);
      compBox.push(box);
      out.forEach(v => compOutliers.push([i, v]));
    });
    mkChart('prok-cm2-comp-chart', {
      title: { text: 'Completeness Distribution (%)' },
      tooltip: { trigger: 'item' },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 45 }, boundaryGap: true },
      yAxis: { type: 'value', name: 'Completeness (%)', min: 0, max: 100 },
      series: [
        { name: 'Completeness', type: 'boxplot', data: compBox,
          itemStyle: { color: PAL[0], borderColor: PAL[1] } },
        { name: 'Outlier', type: 'scatter', data: compOutliers,
          symbolSize: 4, itemStyle: { color: PAL[3], opacity: 0.5 } },
      ],
      grid: { bottom: 90 },
    });

    // Contamination distribution — boxplot per sample
    const contBox = [];
    const contOutliers = [];
    samples.forEach((s, i) => {
      const vals = (cm[s] || []).map(r => +(r.Contamination || 0));
      const { box, outliers: out } = window.boxStats(vals);
      contBox.push(box);
      out.forEach(v => contOutliers.push([i, v]));
    });
    mkChart('prok-cm2-cont-chart', {
      title: { text: 'Contamination Distribution (%)' },
      tooltip: { trigger: 'item' },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 45 }, boundaryGap: true },
      yAxis: { type: 'value', name: 'Contamination (%)', min: 0 },
      series: [
        { name: 'Contamination', type: 'boxplot', data: contBox,
          itemStyle: { color: PAL[0], borderColor: PAL[1] } },
        { name: 'Outlier', type: 'scatter', data: contOutliers,
          symbolSize: 4, itemStyle: { color: PAL[3], opacity: 0.5 } },
      ],
      grid: { bottom: 90 },
    });

    // MIMAG table
    const mimag = typeof MIMAG !== 'undefined' ? MIMAG : {};
    const mimagRows = samples.map(s => ({ sample: s, ...mimag[s] }));
    makeTable('mimag-table', mimagRows, [
      { key: 'sample', label: 'Sample' },
      { key: 'HQ',     label: 'HQ (≥90% comp, ≤5% cont)' },
      { key: 'MQ',     label: 'MQ (≥50% comp, ≤10% cont)' },
      { key: 'LQ',     label: 'LQ (other)' },
      { key: 'total',  label: 'Total' },
    ]);
  }

  // ── Taxonomy ──────────────────────────────────────────────────────────────
  function _renderTaxonomy(sample) {
    const mp = typeof MERGED_PROK !== 'undefined' ? MERGED_PROK.filter(r => r.sample === sample) : [];

    // Domain bar
    const domCount = {};
    mp.forEach(r => { const d = r.Domain || 'Unknown'; domCount[d] = (domCount[d] || 0) + 1; });
    const domEntries = Object.entries(domCount);

    mkChart('prok-domain-bar-chart', {
      title: { text: `${sample} — Domain Distribution` },
      tooltip: { trigger: 'axis' },
      xAxis: { type: 'category', data: domEntries.map(x => x[0]) },
      yAxis: { type: 'value', name: 'MAGs' },
      series: [{
        type: 'bar',
        data: domEntries.map(x => x[1]),
        itemStyle: { color: '#0d9488' },
      }],
    });

    // Top phyla bar
    const phylaCount = {};
    mp.forEach(r => { const p = r.Phylum || 'Unknown'; phylaCount[p] = (phylaCount[p] || 0) + 1; });
    const topPhyla = Object.entries(phylaCount).sort((a, b) => b[1] - a[1]).slice(0, 20);

    mkChart('prok-phyla-chart', {
      title: { text: `${sample} — Top Phyla` },
      tooltip: { trigger: 'axis' },
      xAxis: { type: 'value', name: 'MAGs' },
      yAxis: { type: 'category', data: topPhyla.map(x => x[0]).reverse(), axisLabel: { width: 160, overflow: 'truncate' } },
      series: [{ type: 'bar', data: topPhyla.map(x => x[1]).reverse(), itemStyle: { color: '#d97706' } }],
      grid: { left: 180, right: 30 },
    });

    // Table
    makeTable('prok-tax-table', mp, [
      { key: 'Bin',        label: 'Bin' },
      { key: 'Domain',     label: 'Domain' },
      { key: 'Phylum',     label: 'Phylum' },
      { key: 'Genus',      label: 'Genus' },
      { key: 'Source_tax', label: 'Source' },
      { key: 'Completeness',  label: 'Comp %' },
      { key: 'Contamination', label: 'Cont %' },
    ], { searchId: 'prok-tax-search' });
  }

})();
