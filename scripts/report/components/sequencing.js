/* sequencing.js — FastQC, trimming, mapping, assembly charts */
(function () {
  'use strict';

  window.renderSequencing = function () {
    const samples = typeof SAMPLES !== 'undefined' ? SAMPLES : [];
    const fastp   = typeof FASTP   !== 'undefined' ? FASTP   : {};
    const ov      = typeof OVERVIEW !== 'undefined' ? OVERVIEW : {};
    const quast   = typeof QUAST   !== 'undefined' ? QUAST   : {};
    const mapping = typeof MAPPING !== 'undefined' ? MAPPING : {};

    // ── Raw QC ────────────────────────────────────────────────────────────────
    // Grouped bar: total reads per stage/read label
    const stageSet = new Set();
    samples.forEach(s => (fastp[s]?.reads || []).forEach(r => stageSet.add(r.r_label)));
    const stageLabels = [...stageSet];

    const readSeries = stageLabels.map((lbl, i) => ({
      name: lbl,
      type: 'bar',
      data: samples.map(s => {
        const r = (fastp[s]?.reads || []).find(x => x.r_label === lbl);
        return r ? r.total_sequences : 0;
      }),
    }));

    mkChart('seq-reads-chart', samplesBar({
      samples, title: 'Total Reads — Raw vs Trimmed', valueName: 'Reads',
      series: readSeries.map(s => ({ name: s.name, data: s.data, color: s.color || (s.itemStyle || {}).color })),
    }));

    // Mean quality bar. fastp reports one aggregate Q30 per sample, not a
    // per-read distribution, so a bar is the honest form here (there is no
    // distribution to draw) -- see REPORT_VIZ_GUIDE anti-patterns.
    mkChart('seq-qual-chart', samplesBar({
      samples, title: 'Mean Quality (Q30 estimate)', valueName: 'Mean Q',
      valueMin: 0, valueMax: 40,
      series: ['raw', 'trimmed'].map((stage, i) => ({
        name: stage === 'raw' ? 'Raw R1' : 'Trimmed R1',
        data: samples.map(s => {
          const r = (fastp[s]?.reads || []).find(x => x.stage === stage && x.read === 'R1');
          return r ? +r.mean_quality.toFixed(1) : 0;
        }),
      })),
    }));

    // GC content. Same caveat as quality: fastp gives an aggregate GC per
    // sample, so this cannot become a density without new pipeline output.
    mkChart('seq-gc-chart', samplesBar({
      samples, title: 'GC Content (%)', valueName: 'GC %',
      valueMin: 0, valueMax: 100,
      series: ['raw', 'trimmed'].map(stage => ({
        name: stage === 'raw' ? 'Raw R1' : 'Trimmed R1',
        data: samples.map(s => {
          const r = (fastp[s]?.reads || []).find(x => x.stage === stage && x.read === 'R1');
          return r ? +r.gc_percent.toFixed(1) : 0;
        }),
      })),
    }));

    // Retention rate bar
    mkChart('seq-trim-chart', samplesBar({
      samples, title: 'Read Retention After Trimming (%)', valueName: '%',
      valueMin: 0, valueMax: 105,
      refLines: [{ value: 80, label: '80%' }],
      series: [{
        name: 'Retained', color: '#0d9488',
        data: samples.map(s => {
          const t = fastp[s]?.trim || {};
          return +(((t.reads_written || 0) / (t.reads_in || 1)) * 100).toFixed(1);
        }),
      }],
    }));

    // ── Trimming tab ──────────────────────────────────────────────────────────
    mkChart('seq-adapter-chart', samplesBar({
      samples, title: 'Adapter Trimming Rate (%)', valueName: '%', valueMin: 0,
      series: [
        { name: 'R1 adapter %', data: samples.map(s => +((fastp[s]?.trim?.adapter_r1_pct) || 0).toFixed(1)) },
        { name: 'R2 adapter %', data: samples.map(s => +((fastp[s]?.trim?.adapter_r2_pct) || 0).toFixed(1)) },
      ],
    }));

    mkChart('seq-bp-chart', samplesBar({
      samples, title: 'Base-pair Removed (%)', valueName: '%', valueMin: 0,
      series: [{ name: 'bp removed', color: '#d97706',
        data: samples.map(s => +((fastp[s]?.trim?.bp_removed_pct) || 0).toFixed(1)) }],
    }));

    // ── Mapping tab ───────────────────────────────────────────────────────────
    const mapRates = samples.map(s => +(mapping[s] || 0));
    mkChart('seq-mapping-chart', samplesBar({
      samples, title: 'Read Mapping Rate to Assembly (%)', valueName: '%',
      valueMin: 0, valueMax: 105,
      refLines: [{ value: 70, label: '70%' }],
      series: [{ name: 'Mapped', data: mapRates, color: '#0891b2' }],
    }));

    // Contig length distribution per sample (log axis — heavily right-skewed).
    const lenData = typeof VIRAL_LENGTHS !== 'undefined' ? VIRAL_LENGTHS : {};
    mkChart('seq-depth-chart', distPlot({
      groups: samples.map(s => ({ name: s, values: (lenData[s] || []).filter(v => v > 0) })),
      title: 'Viral Contig Length Distribution (bp)',
      xName: 'Length (bp)',
      log: true,
    }));

    // ── Assembly tab ──────────────────────────────────────────────────────────
    // Since item (d) of docs/ROADMAP_SIMPLIFICACAO.md there is a single
    // contig set per sample (the assembly itself -- no merge, no dedup),
    // so QUAST evaluates only one stage. `asmStages` stays an array (rather
    // than a bare string) so the chart-building code below keeps working
    // unchanged if a future stage is ever reintroduced.
    const asmStages = ['assembly'];
    const metrics   = ['N50', 'Total length', '# contigs'];

    const asmRamp = ordinalRamp(asmStages.length);
    const progSeries = asmStages.map((stage, i) => ({
      name: stage,
      type: 'bar',
      data: samples.map(s => +(quast[s]?.[stage]?.['N50'] || 0)),
      color: asmRamp[i],
    }));

    mkChart('seq-asm-prog-chart', samplesBar({
      samples, title: 'Assembly Quality — N50 (bp)', valueName: 'N50 (bp)',
      series: progSeries.map(s => ({ name: s.name, data: s.data, color: s.color })),
    }));

    // Final assembly quality: N50, total length and contig count are three
    // measures on wildly different scales. Plotted as bars on one value axis
    // (as before) the largest silently flattens the other two -- the dual-axis
    // problem in disguise. As a heatmap each metric is normalised to its own
    // max, so samples are comparable within a metric and the grid scales to any
    // sample count (REPORT_VIZ_GUIDE §3).
    const fmtMetric = (m, v) =>
      m === '# contigs' ? Math.round(v).toLocaleString()
                        : (v >= 1e6 ? (v / 1e6).toFixed(2) + ' Mb'
                        : v >= 1e3 ? (v / 1e3).toFixed(1) + ' kb' : Math.round(v));
    const rawVals = metrics.map(m => samples.map(s => +(quast[s]?.assembly?.[m] || 0)));
    const maxVals = rawVals.map(row => Math.max(...row, 1));
    const heatData = [];
    metrics.forEach((m, mi) => samples.forEach((s, si) => {
      heatData.push({ value: [si, mi, +(rawVals[mi][si] / maxVals[mi] * 100).toFixed(1)],
                      raw: rawVals[mi][si] });
    }));
    mkChart('seq-asm-final-chart', {
      __height: Math.max(240, 120 + metrics.length * 34),
      title: { text: 'Final Assembly Quality — % of best sample per metric' },
      tooltip: {
        trigger: 'item',
        formatter: p => `${samples[p.value[0]]}<br>${metrics[p.value[1]]}: `
                      + `${fmtMetric(metrics[p.value[1]], p.data.raw)}<br>${p.value[2]}% of best`,
      },
      legend: { show: false },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'category', data: metrics },
      visualMap: {
        min: 0, max: 100, calculable: false, show: true,
        orient: 'horizontal', left: 'center', bottom: 4, itemHeight: 90,
        inRange: { color: ['#f0fdfa', '#5eead4', '#0d9488', '#134e4a'] },
        textStyle: { fontSize: 10 },
      },
      series: [{ type: 'heatmap', data: heatData,
                 itemStyle: { borderColor: 'transparent', borderWidth: 2 } }],
      grid: { top: 44, bottom: 78, left: 12, right: 20, containLabel: true },
    });

    // Total assembly length, distribution across samples.
    mkChart('seq-asm-len-chart', distPlot({
      groups: asmStages.map(st => ({
        name: st,
        values: samples.map(s => +(quast[s]?.[st]?.['Total length'] || 0)).filter(v => v > 0),
      })),
      title: 'Total Assembly Length (bp)',
      xName: 'Total length (bp)',
      colors: window.ordinalRamp(asmStages.length),
    }));
  };

})();
