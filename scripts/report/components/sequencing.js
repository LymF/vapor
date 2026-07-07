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

    mkChart('seq-reads-chart', {
      title: { text: 'Total Reads — Raw vs Trimmed' },
      tooltip: { trigger: 'axis' },
      legend: { data: stageLabels },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: 'Reads' },
      series: readSeries,
      grid: { bottom: 70 },
    });

    // Mean quality bar
    mkChart('seq-qual-chart', {
      title: { text: 'Mean Quality (Q30 estimate)' },
      tooltip: { trigger: 'axis', formatter: '{b}: {c:.1f}' },
      legend: { data: ['Raw R1', 'Trimmed R1'] },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: 'Mean Q', min: 0, max: 40 },
      series: [
        {
          name: 'Raw R1',
          type: 'bar',
          data: samples.map(s => {
            const r = (fastp[s]?.reads || []).find(x => x.stage === 'raw' && x.read === 'R1');
            return r ? +r.mean_quality.toFixed(1) : 0;
          }),
        },
        {
          name: 'Trimmed R1',
          type: 'bar',
          data: samples.map(s => {
            const r = (fastp[s]?.reads || []).find(x => x.stage === 'trimmed' && x.read === 'R1');
            return r ? +r.mean_quality.toFixed(1) : 0;
          }),
        },
      ],
      grid: { bottom: 70 },
    });

    // GC content bar
    mkChart('seq-gc-chart', {
      title: { text: 'GC Content (%)' },
      tooltip: { trigger: 'axis' },
      legend: { data: ['Raw R1', 'Trimmed R1'] },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: 'GC %', min: 0, max: 100 },
      series: [
        {
          name: 'Raw R1',
          type: 'bar',
          data: samples.map(s => {
            const r = (fastp[s]?.reads || []).find(x => x.stage === 'raw' && x.read === 'R1');
            return r ? +r.gc_percent.toFixed(1) : 0;
          }),
        },
        {
          name: 'Trimmed R1',
          type: 'bar',
          data: samples.map(s => {
            const r = (fastp[s]?.reads || []).find(x => x.stage === 'trimmed' && x.read === 'R1');
            return r ? +r.gc_percent.toFixed(1) : 0;
          }),
        },
      ],
      grid: { bottom: 70 },
    });

    // Retention rate bar
    mkChart('seq-trim-chart', {
      title: { text: 'Read Retention After Trimming (%)' },
      tooltip: { trigger: 'axis', formatter: '{b}: {c:.1f}%' },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: '%', min: 0, max: 105 },
      series: [{
        type: 'bar',
        data: samples.map(s => {
          const t = fastp[s]?.trim || {};
          const rin  = t.reads_in  || 1;
          const rout = t.reads_written || 0;
          return +((rout / rin) * 100).toFixed(1);
        }),
        itemStyle: { color: '#0d9488' },
        markLine: { data: [{ yAxis: 80, name: '80%', lineStyle: { type: 'dashed', color: '#d97706' } }] },
      }],
      grid: { bottom: 70 },
    });

    // ── Trimming tab ──────────────────────────────────────────────────────────
    mkChart('seq-adapter-chart', {
      title: { text: 'Adapter Trimming Rate (%)' },
      tooltip: { trigger: 'axis' },
      legend: { data: ['R1 adapter %', 'R2 adapter %'] },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: '%', min: 0 },
      series: [
        { name: 'R1 adapter %', type: 'bar', data: samples.map(s => +((fastp[s]?.trim?.adapter_r1_pct) || 0).toFixed(1)) },
        { name: 'R2 adapter %', type: 'bar', data: samples.map(s => +((fastp[s]?.trim?.adapter_r2_pct) || 0).toFixed(1)) },
      ],
      grid: { bottom: 70 },
    });

    mkChart('seq-bp-chart', {
      title: { text: 'Base-pair Removed (%)' },
      tooltip: { trigger: 'axis', formatter: '{b}: {c:.1f}%' },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: '%', min: 0 },
      series: [{
        type: 'bar',
        data: samples.map(s => +((fastp[s]?.trim?.bp_removed_pct) || 0).toFixed(1)),
        itemStyle: { color: '#d97706' },
      }],
      grid: { bottom: 70 },
    });

    // ── Mapping tab ───────────────────────────────────────────────────────────
    const mapRates = samples.map(s => +(mapping[s] || 0));
    mkChart('seq-mapping-chart', {
      title: { text: 'Read Mapping Rate to Assembly (%)' },
      tooltip: { trigger: 'axis', formatter: '{b}: {c:.1f}%' },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: '%', min: 0, max: 105 },
      series: [{
        type: 'bar',
        data: mapRates,
        itemStyle: { color: '#0891b2' },
        markLine: { data: [{ yAxis: 70, name: '70%', lineStyle: { type: 'dashed', color: '#d97706' } }] },
      }],
      grid: { bottom: 70 },
    });

    // Contig length distribution (boxplot per sample — scales to any sample count)
    const lenData = typeof VIRAL_LENGTHS !== 'undefined' ? VIRAL_LENGTHS : {};

    const boxData  = [];
    const outliers = [];
    samples.forEach((s, i) => {
      const { box, outliers: out } = window.boxStats(lenData[s] || []);
      boxData.push(box);
      out.forEach(v => outliers.push([i, v]));
    });

    mkChart('seq-depth-chart', {
      title: { text: 'Viral Contig Length Distribution (bp)' },
      tooltip: { trigger: 'item' },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 45 }, boundaryGap: true },
      yAxis: { type: 'log', name: 'Length (bp)' },
      series: [
        { name: 'Length',  type: 'boxplot', data: boxData,
          itemStyle: { color: PAL[0], borderColor: PAL[1] } },
        { name: 'Outlier', type: 'scatter', data: outliers,
          symbolSize: 4, itemStyle: { color: PAL[3], opacity: 0.5 } },
      ],
      grid: { bottom: 90 },
    });

    // ── Assembly tab ──────────────────────────────────────────────────────────
    const asmStages = ['MEGAHIT', 'metaSPAdes', 'metaviralSPAdes', 'merged_filtered', 'deduplicated'];
    const metrics   = ['N50', 'Total length', '# contigs'];

    // Assembly progression per sample, grouped by stage. Stages are ORDINAL
    // (each is a later step of the same pipeline, not an unrelated
    // category) -- a single-hue light->dark ramp reads as progression;
    // 5 arbitrary categorical hues would not.
    const asmRamp = ordinalRamp(asmStages.length);
    const progSeries = asmStages.map((stage, i) => ({
      name: stage,
      type: 'bar',
      data: samples.map(s => +(quast[s]?.[stage]?.['N50'] || 0)),
      color: asmRamp[i],
    }));

    mkChart('seq-asm-prog-chart', {
      title: { text: 'Assembly Progression — N50 (bp)' },
      tooltip: { trigger: 'axis' },
      legend: { data: asmStages },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: 'N50 (bp)' },
      series: progSeries,
      grid: { bottom: 70 },
    });

    // Final assembly (deduplicated) — 3 metrics grouped
    const finalSeries = metrics.map((m, i) => ({
      name: m,
      type: 'bar',
      data: samples.map(s => +(quast[s]?.deduplicated?.[m] || 0)),
      color: PAL[i % PAL.length],
    }));

    mkChart('seq-asm-final-chart', {
      title: { text: 'Final Assembly Quality (deduplicated)' },
      tooltip: { trigger: 'axis' },
      legend: { data: metrics },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: 'Value' },
      series: finalSeries,
      grid: { bottom: 70 },
    });

    // Total assembly length per stage — boxplot across samples (one box per stage)
    const asmLenBox = [];
    const asmLenOutliers = [];
    asmStages.forEach((st, i) => {
      const vals = samples.map(s => +(quast[s]?.[st]?.['Total length'] || 0)).filter(v => v > 0);
      const { box, outliers: out } = window.boxStats(vals);
      asmLenBox.push(box);
      out.forEach(v => asmLenOutliers.push([i, v]));
    });

    mkChart('seq-asm-len-chart', {
      title: { text: 'Total Assembly Length per Stage (bp)' },
      tooltip: { trigger: 'item' },
      xAxis: { type: 'category', data: asmStages, axisLabel: { rotate: 20 }, boundaryGap: true },
      yAxis: { type: 'value', name: 'Total length (bp)' },
      series: [
        { name: 'Total length', type: 'boxplot', data: asmLenBox,
          itemStyle: { color: PAL[0], borderColor: PAL[1] } },
        { name: 'Outlier', type: 'scatter', data: asmLenOutliers,
          symbolSize: 4, itemStyle: { color: PAL[3], opacity: 0.5 } },
      ],
      grid: { bottom: 70 },
    });
  };

})();
