/* overview.js — General Statistics KPI cards + cross-sample summary charts */
(function () {
  'use strict';

  function fmt(n) {
    if (n === undefined || n === null) return 'N/A';
    if (typeof n === 'string' && isNaN(Number(n.replace(/[,%]/g,'')))) return n;
    const num = typeof n === 'number' ? n : Number(String(n).replace(/[,%]/g,''));
    if (isNaN(num)) return String(n);
    if (num >= 1e9) return (num / 1e9).toFixed(2) + 'G';
    if (num >= 1e6) return (num / 1e6).toFixed(2) + 'M';
    if (num >= 1e3) return (num / 1e3).toFixed(1) + 'k';
    return num.toLocaleString();
  }

  // ── KPI summary (cross-sample totals) ────────────────────────────────────
  function buildGeneralStats() {
    const samples = typeof SAMPLES !== 'undefined' ? SAMPLES : [];
    const ov      = typeof OVERVIEW !== 'undefined' ? OVERVIEW : {};

    // Aggregate
    let totalReads = 0, totalViral = 0, totalMAGs = 0;
    let hqViral = 0, totalDefense = 0, totalAmrHq = 0, totalAmgs = 0, totalVmags = 0;
    samples.forEach(s => {
      const d = ov[s] || {};
      totalReads   += Number(d.total_raw_reads)    || 0;
      totalViral   += Number(d.viral_consensus)    || 0;
      totalMAGs    += Number(d.total_bins)         || 0;
      hqViral      += Number(d.complete_viral)     || 0;
      totalDefense += Number(d.total_defense)      || 0;
      totalAmrHq   += Number(d.total_amr_hq)       || 0;
      totalAmgs    += Number(d.total_amgs)         || 0;
      totalVmags   += Number(d.vmags)              || 0;
    });

    const mimag = typeof MIMAG !== 'undefined' ? MIMAG : {};
    let hqMqMags = 0;
    samples.forEach(s => {
      const m = mimag[s] || {};
      hqMqMags += (Number(m.HQ) || 0) + (Number(m.MQ) || 0);
    });

    const kpis = [
      { val: samples.length,    label: 'Samples' },
      { val: fmt(totalReads),   label: 'Total reads',       sub: 'raw' },
      { val: fmt(totalViral),   label: 'Viral vOTUs',       sub: 'consensus' },
      { val: fmt(totalVmags),   label: 'vMAGs',             sub: 'vRhyme bins' },
      { val: fmt(totalMAGs),    label: 'MAGs',              sub: 'Binette final' },
      { val: fmt(hqViral),      label: 'HQ viral',          sub: 'CheckV Complete/HQ' },
      { val: fmt(hqMqMags),     label: 'HQ+MQ MAGs',        sub: 'CheckM2 MIMAG' },
      { val: fmt(totalDefense), label: 'Defense systems',   sub: 'DefenseFinder' },
      { val: fmt(totalAmrHq),   label: 'AMR genes (3-tool)',sub: 'consensus' },
      { val: fmt(totalAmgs),    label: 'AMGs',              sub: 'VIBRANT' },
    ];

    const grid = document.getElementById('kpi-grid');
    if (!grid) return;
    grid.innerHTML = kpis.map(k =>
      `<div class="kpi-card">
         <div class="kpi-val">${k.val}</div>
         <div class="kpi-label">${k.label}</div>
         ${k.sub ? `<div class="kpi-sub">${k.sub}</div>` : ''}
       </div>`
    ).join('');

    // ── Cross-sample charts ────────────────────────────────────────────────
    const readVals    = samples.map(s => Number((ov[s] || {}).total_raw_reads)    || 0);
    const trimVals    = samples.map(s => Number((ov[s] || {}).total_trimmed_reads) || 0);
    const viralVals   = samples.map(s => Number((ov[s] || {}).viral_consensus)     || 0);
    const hqVirVals   = samples.map(s => Number((ov[s] || {}).complete_viral)      || 0);
    const magVals     = samples.map(s => Number((ov[s] || {}).total_bins)          || 0);

    const magHQ  = samples.map(s => (mimag[s] || {}).HQ || 0);
    const magMQ  = samples.map(s => (mimag[s] || {}).MQ || 0);
    const magLQ  = samples.map(s => (mimag[s] || {}).LQ || 0);

    mkChart('ov-reads-chart', {
      title: { text: 'Reads per Sample' },
      tooltip: { trigger: 'axis' },
      legend: { data: ['Raw', 'Trimmed'] },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: 'Reads' },
      series: [
        { name: 'Raw',     type: 'bar', data: readVals,  color: '#0891b2' },
        { name: 'Trimmed', type: 'bar', data: trimVals,  color: '#0d9488' },
      ],
      grid: { bottom: 60 },
    });

    mkChart('ov-votus-chart', {
      title: { text: 'Viral OTUs per Sample' },
      tooltip: { trigger: 'axis' },
      legend: { data: ['Total vOTUs', 'High-quality'] },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: 'vOTUs' },
      series: [
        { name: 'Total vOTUs',  type: 'bar', data: viralVals, color: '#0d9488' },
        { name: 'High-quality', type: 'bar', data: hqVirVals, color: '#16a34a' },
      ],
      grid: { bottom: 60 },
    });

    mkChart('ov-mags-chart', {
      title: { text: 'MAGs per Sample (MIMAG tiers)' },
      tooltip: { trigger: 'axis' },
      legend: { data: ['HQ', 'MQ', 'LQ'] },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: 'MAGs' },
      series: [
        { name: 'HQ', type: 'bar', stack: 'mag', data: magHQ, color: '#16a34a', itemStyle: { borderRadius: [0,0,0,0] } },
        { name: 'MQ', type: 'bar', stack: 'mag', data: magMQ, color: '#d97706' },
        { name: 'LQ', type: 'bar', stack: 'mag', data: magLQ, color: '#ef4444', itemStyle: { borderRadius: [4,4,0,0] } },
      ],
      grid: { bottom: 60 },
    });

    // Funnel (average across all samples)
    const n = samples.length || 1;
    const avgRaw   = readVals.reduce((a,b) => a+b, 0) / n;
    const avgTrim  = trimVals.reduce((a,b) => a+b, 0) / n;
    const avgContig= samples.reduce((a,s) => a + (Number((ov[s]||{}).n_contigs) || 0), 0) / n;
    const avgViral = viralVals.reduce((a,b) => a+b, 0) / n;
    const avgMAG   = magVals.reduce((a,b) => a+b, 0) / n;

    // Pipeline stages are ORDINAL (each stage is a narrower subset of the
    // one before it) -- a single-hue light->dark ramp reads as "funneling
    // down" at a glance; 5 unrelated categorical hues would read as 5
    // disconnected categories instead of one narrowing sequence.
    const funnelRamp = ordinalRamp(5);
    const funnelStages = [
      { name: 'Raw reads',   value: Math.round(avgRaw) },
      { name: 'Trimmed',     value: Math.round(avgTrim) },
      { name: 'Contigs',     value: Math.round(avgContig) },
      { name: 'Viral vOTUs', value: Math.round(avgViral) },
      { name: 'MAGs',        value: Math.round(avgMAG) },
    ];
    mkChart('ov-funnel-chart', {
      title: { text: 'Read Pipeline Funnel (average across samples)' },
      tooltip: { trigger: 'item', formatter: '{b}: {c}' },
      series: [{
        type: 'funnel',
        left: '5%', width: '90%',
        sort: 'none',
        gap: 2,
        label: { show: true, position: 'inside', fontSize: 12, formatter: '{b}\n{c}' },
        data: funnelStages.map((d, i) => ({
          ...d,
          itemStyle: { color: funnelRamp[i] },
          // First two steps of the ramp are light enough that white label
          // text loses contrast -- switch to dark ink for those, white for
          // the rest (same relief the ordinal-ramp rule requires).
          label: { color: i < 2 ? '#0f172a' : '#fff' },
        })),
      }],
    });
  }

  // ── Per-sample drill-down ─────────────────────────────────────────────────
  function buildPerSample(sample) {
    const ov   = (typeof OVERVIEW !== 'undefined' ? OVERVIEW : {})[sample] || {};
    const fastp = (typeof FASTP   !== 'undefined' ? FASTP   : {})[sample] || {};
    const binner= (typeof BINNER  !== 'undefined' ? BINNER  : {})[sample] || {};

    const kpis = [
      // ── Sequencing & Assembly ──────────────────────────────────────────────
      { val: fmt(ov.total_raw_reads),  label: 'Raw reads' },
      { val: ov.mean_qual || 'N/A',    label: 'Mean quality (Q)' },
      { val: ov.gc_pct    || 'N/A',    label: 'GC content' },
      { val: ov.mapping_rate || 'N/A', label: 'Mapping rate' },
      { val: fmt(ov.n_contigs),        label: 'Contigs' },
      { val: fmt(ov.n50),              label: 'N50 (bp)' },
      // ── Viral ─────────────────────────────────────────────────────────────
      { val: fmt(ov.viral_consensus),  label: 'Viral vOTUs' },
      { val: fmt(ov.complete_viral),   label: 'HQ/Complete viral' },
      { val: ov.pct_novel != null
              ? `${ov.pct_novel}%`
              : 'N/A',              label: '% novel vOTUs' },
      { val: fmt(ov.vmags),            label: 'vMAGs (vRhyme)' },
      { val: ov.lytic_count != null
              ? `${ov.lytic_count}L / ${ov.lysogenic_count}T`
              : 'N/A',              label: 'Lytic / Temperate' },
      { val: ov.lytic_ratio != null
              ? `${(ov.lytic_ratio * 100).toFixed(0)}% lytic`
              : 'N/A',              label: 'Lifestyle ratio' },
      { val: fmt(ov.total_amgs),       label: 'AMGs (VIBRANT)' },
      // ── Taxonomy & Hosts ──────────────────────────────────────────────────
      { val: fmt(ov.taxonomy_classified), label: 'Viral taxa classified' },
      { val: fmt(ov.host_pred_total),  label: 'Host predictions (PHIST)' },
      // ── Prokaryotes ───────────────────────────────────────────────────────
      { val: fmt(ov.total_bins),       label: 'MAGs (Binette)' },
      { val: fmt(ov.hq_bins),          label: 'HQ MAGs (≥90/≤5)' },
      { val: fmt(ov.bacteria_bins),    label: 'Bacterial MAGs' },
      { val: fmt(ov.archaea_bins),     label: 'Archaeal MAGs' },
      { val: fmt(ov.gtdb_classified),  label: 'GTDB-classified bins' },
      // ── Defense & AMR ─────────────────────────────────────────────────────
      { val: fmt(ov.total_defense),    label: 'Defense systems' },
      { val: fmt(ov.total_amr_hq),     label: 'AMR genes (3-tool)' },
    ];

    const cards = document.getElementById('ov-sample-cards');
    if (cards) {
      cards.innerHTML = kpis.map(k =>
        `<div class="kpi-card"><div class="kpi-val">${k.val}</div><div class="kpi-label">${k.label}</div></div>`
      ).join('');
    }

    // Reads breakdown (fastp stages)
    const reads = (fastp.reads || []);
    const stageLabels = [...new Set(reads.map(r => r.r_label))];
    const readCounts  = stageLabels.map(l => {
      const r = reads.find(x => x.r_label === l);
      return r ? r.total_sequences : 0;
    });

    mkChart('ov-sample-reads', {
      title: { text: `${sample} — Read Counts by Stage` },
      tooltip: { trigger: 'axis' },
      xAxis: { type: 'category', data: stageLabels, axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: 'Reads' },
      series: [{ type: 'bar', data: readCounts }],
      grid: { bottom: 70 },
    });

    // Viral summary
    const vTools = (typeof VIRAL_TOOLS !== 'undefined' ? VIRAL_TOOLS : {})[sample] || {};
    const toolNames = Object.keys(vTools);
    mkChart('ov-sample-viral', {
      title: { text: `${sample} — Viral Detection by Tool` },
      tooltip: { trigger: 'axis' },
      xAxis: { type: 'category', data: toolNames },
      yAxis: { type: 'value', name: 'Contigs' },
      series: [{ type: 'bar', data: toolNames.map(t => vTools[t] || 0) }],
    });
  }

  window.renderOverview = function () {
    buildGeneralStats();
    makeSampleDropdown('sample-sel-overview', buildPerSample);
  };

})();
