/* coassembly.js — Co-assembly per-group tab.
   Mirrors the per-sample Viral + Prokaryotic tabs for co-assembly GROUPS,
   reusing the shared helpers (mkChart, makeTable, boxStats, qualBadge, PAL).
   Reads COASSEMBLY_DATA (compact counts) + COAS_RICH (per-group checkv/checkm2/
   taxonomy/gtdb). Report-only parity: no vOTU lifestyle/abundance (there is no
   group make_votu_table). A group selector drives the detail charts; boxplots,
   MIMAG and overview counts always compare across all groups. */
(function () {
  'use strict';

  const ALL = '__all__';
  const cvCats = ['Complete', 'High-quality', 'Medium-quality', 'Low-quality', 'Not-determined'];
  const cvCols = ['#16a34a', '#4ade80', '#fbbf24', '#d97706', '#ef4444'];
  const qualColorMap = {
    'Complete': '#16a34a', 'High-quality': '#4ade80', 'Medium-quality': '#fbbf24',
    'Low-quality': '#d97706', 'Not-determined': '#64748b',
  };

  window.renderCoassembly = function () {
    const CA = typeof COASSEMBLY_DATA !== 'undefined' ? COASSEMBLY_DATA : null;
    const R  = typeof COAS_RICH !== 'undefined' ? COAS_RICH : { units: [] };
    const empty   = document.getElementById('coassembly-empty');
    const content = document.getElementById('coas-content');

    if (!CA || !CA.has_data) {
      if (empty)   empty.style.display = '';
      if (content) content.style.display = 'none';
      return;
    }
    if (empty)   empty.style.display = 'none';
    if (content) content.style.display = '';

    const units = (R.units && R.units.length) ? R.units : CA.groups.map(g => g.group);

    _renderOverview(CA);
    _renderStatic(R, units);          // cross-group boxplots + MIMAG (selector-independent)

    _makeGroupDropdown('coas-group-sel', units, unit => {
      _renderViral(R, unit);
      _renderProk(R, unit);
    });
  };

  // ── Group dropdown (mirrors makeSampleDropdown, with "All groups") ──────────
  function _makeGroupDropdown(selId, units, onChange) {
    const sel = document.getElementById(selId);
    if (!sel) return;
    sel.innerHTML = '';
    const optAll = document.createElement('option');
    optAll.value = ALL; optAll.textContent = 'All groups';
    sel.appendChild(optAll);
    units.forEach(u => {
      const o = document.createElement('option');
      o.value = o.textContent = u;
      sel.appendChild(o);
    });
    sel.addEventListener('change', () => onChange(sel.value));
    onChange(sel.value);
  }

  // ── Overview (compact counts + tables from COASSEMBLY_DATA) ─────────────────
  function _renderOverview(CA) {
    const groups = CA.groups || [];

    mkChart('coas-votu-count-chart', {
      title: { text: 'vOTUs & vMAGs per Group' },
      tooltip: { trigger: 'axis' },
      legend: { data: ['vOTUs', 'vMAGs'] },
      xAxis: { type: 'category', data: groups.map(g => g.group), axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: 'Count' },
      series: [
        { name: 'vOTUs', type: 'bar', color: '#0d9488', data: groups.map(g => g.n_votus || 0) },
        { name: 'vMAGs', type: 'bar', color: '#0891b2', data: groups.map(g => g.n_vmags || 0) },
      ],
      grid: { bottom: 70 },
    });

    mkChart('coas-mag-count-chart', {
      title: { text: 'MAGs per Group (CheckM2)' },
      tooltip: { trigger: 'axis' },
      xAxis: { type: 'category', data: groups.map(g => g.group), axisLabel: { rotate: 30 } },
      yAxis: { type: 'value', name: 'MAGs' },
      series: [{ type: 'bar', color: '#7c3aed', data: groups.map(g => (g.mags || []).length) }],
      grid: { bottom: 70 },
    });

    const magRows = [];
    groups.forEach(g => (g.mags || []).forEach(m => magRows.push({
      group: g.group, bin: m.bin,
      completeness: (+m.completeness).toFixed(1),
      contamination: (+m.contamination).toFixed(1),
      classification: m.classification || '',
    })));
    makeTable('coassembly-table', magRows, [
      { key: 'group', label: 'Group' },
      { key: 'bin', label: 'MAG' },
      { key: 'completeness', label: 'Completeness %' },
      { key: 'contamination', label: 'Contamination %' },
      { key: 'classification', label: 'GTDB' },
    ]);

    const votuRows = groups.filter(g => g.n_votus || g.n_vmags).map(g => ({
      group: g.group, n_votus: g.n_votus || 0, n_vmags: g.n_vmags || 0,
      families: (g.votu_families || []).map(f => `${f.family} (${f.count})`).join(', '),
    }));
    makeTable('coassembly-votus-table', votuRows, [
      { key: 'group', label: 'Group' },
      { key: 'n_votus', label: 'vOTUs' },
      { key: 'n_vmags', label: 'vMAGs' },
      { key: 'families', label: 'Top families' },
    ]);
  }

  // ── Static cross-group charts (boxplots + MIMAG table) ──────────────────────
  function _renderStatic(R, units) {
    const cm   = R.checkm2 || {};
    const vlen = R.vlen || {};

    // Distributions per group. distPlot picks strip vs density vs ridgeline from
    // the data and draws the MIMAG cutoffs (docs/REPORT_VIZ_GUIDE.md §4) — group
    // MAG counts are often small, where a KDE would fabricate a smooth curve.
    const byUnit = (fn) => units.map(u => ({ name: u, values: fn(u) }));
    mkChart('coas-vlen-chart', distPlot({
      groups: byUnit(u => (vlen[u] || []).filter(v => v > 0)),
      title: 'vOTU Length Distribution (bp)', xName: 'Length (bp)', log: true,
    }));
    mkChart('coas-size-chart', distPlot({
      groups: byUnit(u => (cm[u] || []).map(r => +(r.Genome_Size || r.genome_size || 0) / 1e6).filter(v => v > 0)),
      title: 'MAG Genome Size Distribution (Mb)', xName: 'Genome size (Mb)', xMin: 0,
    }));
    mkChart('coas-cm2-comp-chart', distPlot({
      groups: byUnit(u => (cm[u] || []).map(r => +(r.Completeness || 0))),
      title: 'MAG Completeness Distribution (%)', xName: 'Completeness (%)',
      xMin: 0, xMax: 100,
      cutoffs: [{ value: 50, label: 'MQ 50%' }, { value: 90, label: 'HQ 90%' }],
    }));
    mkChart('coas-cm2-cont-chart', distPlot({
      groups: byUnit(u => (cm[u] || []).map(r => +(r.Contamination || 0))),
      title: 'MAG Contamination Distribution (%)', xName: 'Contamination (%)', xMin: 0,
      cutoffs: [{ value: 5, label: 'HQ ≤5%' }, { value: 10, label: 'MQ ≤10%' }],
    }));

    const mimag = R.mimag || {};
    const mimagRows = units.map(u => ({ group: u, ...(mimag[u] || { HQ: 0, MQ: 0, LQ: 0, total: 0 }) }));
    makeTable('coas-mimag-table', mimagRows, [
      { key: 'group', label: 'Group' },
      { key: 'HQ', label: 'HQ (≥90% comp, ≤5% cont)' },
      { key: 'MQ', label: 'MQ (≥50% comp, ≤10% cont)' },
      { key: 'LQ', label: 'LQ (other)' },
      { key: 'total', label: 'Total' },
    ]);
  }

  // ── Viral (selection-driven) ────────────────────────────────────────────────
  function _renderViral(R, unit) {
    const isAll = unit === ALL;
    const label = isAll ? 'all groups' : unit;
    const units = (R.units || []);
    const pick  = isAll ? units : [unit];

    const cvRows  = _concat(R.checkv, pick);
    const vrhRows = _concat(R.checkv_vrh, pick);

    mkChart('coas-checkv-bar-chart',
      _checkvBarOption(`CheckV — vOTU / Consensus Quality (${label})`, _tierCounts(cvRows)));
    mkChart('coas-checkv-vrh-chart',
      _checkvBarOption(`CheckV — vRhyme vMAGs Quality (${label})`, _tierCounts(vrhRows)));

    // Scatter: length vs completeness (● consensus  ◆ vRhyme) + MQ/HQ zones
    const series = [];
    function _addPts(rows, srcName, sym) {
      const seen = {};
      rows.forEach(r => {
        const len = +(r.contig_length || 0), comp = +(r.completeness || 0);
        const q = r.checkv_quality || 'Not-determined';
        if (!len || !comp) return;
        const key = `${q} (${srcName})`;
        if (!seen[key]) {
          seen[key] = { name: key, type: 'scatter', symbol: sym, symbolSize: 8,
            itemStyle: { color: qualColorMap[q] || '#64748b', opacity: 0.75 }, data: [] };
          series.push(seen[key]);
        }
        seen[key].data.push({ value: [len, comp], name: r.contig_id || r.contig || '' });
      });
    }
    _addPts(cvRows, 'consensus', 'circle');
    _addPts(vrhRows, 'vRhyme', 'diamond');
    series.push({
      type: 'scatter', data: [],
      markArea: { silent: true, data: [
        [{ yAxis: 90, itemStyle: { color: 'rgba(22,163,74,0.09)' } }, { yAxis: 105 }],
        [{ yAxis: 50, itemStyle: { color: 'rgba(217,119,6,0.07)' } }, { yAxis: 90 }],
        [{ yAxis: 0,  itemStyle: { color: 'rgba(239,68,68,0.05)' } }, { yAxis: 50 }],
      ] },
      markLine: { silent: true, data: [
        { yAxis: 90, lineStyle: { type: 'dashed', color: '#16a34a' }, label: { formatter: '≥90% HQ', color: '#16a34a', fontSize: 10 } },
        { yAxis: 50, lineStyle: { type: 'dotted', color: '#d97706' }, label: { formatter: '≥50% MQ', color: '#d97706', fontSize: 10 } },
      ] },
    });
    mkChart('coas-checkv-scatter-chart', {
      title: { text: `CheckV — Length vs Completeness (● consensus  ◆ vRhyme) — ${label}` },
      tooltip: { trigger: 'item', formatter: p => `${p.data.name}<br>Length: ${(p.data.value[0] || 0).toLocaleString()} bp<br>Completeness: ${(p.data.value[1] || 0).toFixed(1)}%` },
      legend: { type: 'scroll' },
      xAxis: { type: 'log', name: 'Length (bp)', min: 1000, max: 1000000, nameLocation: 'middle', nameGap: 30 },
      yAxis: { type: 'value', name: 'Completeness (%)', min: 0, max: 105 },
      series, grid: { bottom: 60, top: 70 },
    });

    // Taxonomy: source pie + rank bar + master table
    const tax = isAll ? (R.tax || []) : (R.tax || []).filter(r => r.sample === unit);

    const srcDist = {};
    const allSrc = R.source_dist || {};
    (isAll ? units : [unit]).forEach(u => {
      Object.entries(allSrc[u] || {}).forEach(([k, v]) => { srcDist[k] = (srcDist[k] || 0) + v; });
    });
    const srcRows = Object.entries(srcDist).sort((a, b) => a[1] - b[1]);
    const srcTotal = srcRows.reduce((a, [, v]) => a + v, 0) || 1;
    mkChart('coas-tax-source-chart', {
      title: { text: `${label} — Classification Source` },
      tooltip: { trigger: 'item',
                 formatter: p => `${p.name}: ${p.value} (${(p.value / srcTotal * 100).toFixed(1)}%)` },
      legend: { show: false },
      xAxis: { type: 'value', name: 'Contigs', nameLocation: 'middle', nameGap: 28 },
      yAxis: { type: 'category', data: srcRows.map(r => r[0]) },
      series: [{ type: 'bar', data: srcRows.map(r => r[1]), barMaxWidth: 26,
                 itemStyle: { color: PAL[0], borderRadius: [0, 3, 3, 0] },
                 label: { show: true, position: 'right', fontSize: 10,
                          formatter: p => `${(p.value / srcTotal * 100).toFixed(0)}%` } }],
      grid: { left: 12, right: 52, bottom: 44, containLabel: true },
    });

    _renderRankBar(tax, label);
    _renderAccumulation(unit, isAll, units);

    const tableCard = document.querySelector('#coas-tax-table')?.closest('.chart-card');
    if (tableCard) tableCard.style.display = isAll ? 'none' : '';
    if (!isAll) {
      makeTable('coas-tax-table', tax, [
        { key: 'Genome', label: 'Contig' },
        { key: 'final_order', label: 'Order' },
        { key: 'final_family', label: 'Family' },
        { key: 'final_genus', label: 'Genus' },
        { key: 'Source', label: 'Source' },
        { key: 'CheckV_quality', label: 'CheckV' },
        { key: 'Completeness', label: 'Completeness' },
      ], { searchId: 'coas-tax-search', format: { CheckV_quality: qualBadge } });
    }
  }

  // ── vOTU accumulation (collector curve) ─────────────────────────────────────
  // Only meaningful on the co-assembly track: a group's vOTUs are clustered once
  // over the co-assembled contigs, so vOTU identity is shared across the group's
  // samples and can be accumulated. Sample order is arbitrary, so the loader
  // averages over random permutations; the band is the 10-90 percentile.
  // A curve still climbing at the last sample means the group is NOT saturated —
  // more samples would keep yielding unseen vOTUs.
  function _renderAccumulation(unit, isAll, units) {
    const acc = typeof VOTU_ACCUM !== 'undefined' ? VOTU_ACCUM : {};
    const keys = (isAll ? units : [unit]).filter(u => acc[u]);
    const el = document.getElementById('coas-accum-chart');
    if (!el) return;
    if (!keys.length) {
      mkChart('coas-accum-chart', {
        title: { text: 'vOTU Accumulation' },
        graphic: { type: 'text', left: 'center', top: 'middle',
                   style: { text: 'Needs the group abundance matrix (short-read co-binning)',
                            fill: PAL_MUTED, fontSize: 12 } },
      });
      return;
    }

    const maxN = Math.max(...keys.map(k => acc[k].n_samples));
    const xs = Array.from({ length: maxN }, (_, i) => i + 1);
    const series = [];

    // Percentile band, single group only — overlapping bands are unreadable.
    if (keys.length === 1) {
      const c = acc[keys[0]];
      series.push(
        { name: 'lo', type: 'line', stack: 'band', symbol: 'none', silent: true,
          lineStyle: { opacity: 0 }, areaStyle: { opacity: 0 }, data: c.lo },
        { name: 'band', type: 'line', stack: 'band', symbol: 'none', silent: true,
          lineStyle: { opacity: 0 },
          areaStyle: { color: PAL[0], opacity: 0.16 },
          data: c.hi.map((h, i) => h - c.lo[i]) },
      );
    }
    keys.slice(0, window.VIZ.maxSeries).forEach((k, i) => {
      series.push({
        name: k, type: 'line', symbol: 'circle', symbolSize: 6, smooth: false,
        lineStyle: { width: 2, color: PAL[i] }, itemStyle: { color: PAL[i] },
        data: acc[k].mean,
      });
    });

    const one = keys.length === 1 ? acc[keys[0]] : null;
    mkChart('coas-accum-chart', {
      title: { text: one
        ? `vOTU Accumulation — ${keys[0]} (${one.total.toLocaleString()} vOTUs, `
          + `${one.n_samples} samples, ≥${one.min_depth}× depth)`
        : 'vOTU Accumulation per Group' },
      tooltip: { trigger: 'axis',
        formatter: ps => {
          const p = ps.filter(x => x.seriesName !== 'lo' && x.seriesName !== 'band');
          if (!p.length) return '';
          return `${p[0].axisValue} sample(s)<br>`
               + p.map(x => `${x.marker}${x.seriesName}: ${x.value} vOTUs`).join('<br>');
        } },
      legend: keys.length > 1 ? { data: keys.slice(0, window.VIZ.maxSeries) } : { show: false },
      xAxis: { type: 'category', data: xs, name: 'Samples pooled',
               nameLocation: 'middle', nameGap: 28, boundaryGap: false },
      yAxis: { type: 'value', name: 'Distinct vOTUs' },
      series,
      grid: { top: keys.length > 1 ? 58 : 44, bottom: 52, left: 12, right: 24, containLabel: true },
    });
  }

  function _renderRankBar(tax, label) {
    const RANK_FIELD = { order: 'final_order', family: 'final_family', genus: 'final_genus' };
    const RANK_LABEL = { order: 'Orders', family: 'Families', genus: 'Genera' };
    function draw(level) {
      const field = RANK_FIELD[level] || 'final_family';
      const count = {};
      tax.forEach(r => { const v = r[field] || ''; if (v) count[v] = (count[v] || 0) + 1; });
      const top = Object.entries(count).sort((a, b) => b[1] - a[1]).slice(0, 20);
      mkChart('coas-tax-rank-chart', {
        title: { text: `${label} — Top Viral ${RANK_LABEL[level] || 'Families'}` },
        tooltip: { trigger: 'axis' },
        xAxis: { type: 'value', name: 'Count', nameLocation: 'middle', nameGap: 28 },
        yAxis: { type: 'category', data: top.map(x => x[0]).reverse(), axisLabel: { width: 140, overflow: 'truncate' } },
        series: [{ type: 'bar', data: top.map(x => x[1]).reverse(), itemStyle: { color: '#0d9488' } }],
        grid: { left: 160, right: 30, bottom: 50 },
      });
    }
    const btns = {
      order: document.getElementById('coas-tax-order-btn'),
      family: document.getElementById('coas-tax-family-btn'),
      genus: document.getElementById('coas-tax-genus-btn'),
    };
    const active = Object.entries(btns).find(([, b]) => b?.classList.contains('active'))?.[0] || 'family';
    Object.entries(btns).forEach(([level, b]) => {
      if (!b) return;
      b.onclick = () => {
        Object.values(btns).forEach(x => x && x.classList.remove('active'));
        b.classList.add('active');
        draw(level);
      };
    });
    if (btns.family && !Object.values(btns).some(b => b?.classList.contains('active'))) btns.family.classList.add('active');
    draw(active);
  }

  // ── Prokaryotic (selection-driven) ──────────────────────────────────────────
  function _renderProk(R, unit) {
    const isAll = unit === ALL;
    const label = isAll ? 'all groups' : unit;
    const units = (R.units || []);
    const pick  = isAll ? units : [unit];

    // CheckM2 scatter: completeness vs contamination + MIMAG zones (colour by domain)
    const domColors = { Bacteria: '#0d9488', Archaea: '#d97706', Unknown: '#64748b' };
    const merged = (R.merged_prok || []);
    const domLookup = {};
    merged.forEach(r => { domLookup[`${r.sample}::${r.Bin}`] = r.Domain || 'Unknown'; });

    const seriesMap = {};
    pick.forEach(u => {
      ((R.checkm2 || {})[u] || []).forEach(r => {
        const comp = +(r.Completeness || 0), cont = +(r.Contamination || 0);
        const name = (r.Name || r.name || '').replace(/\.(fa|fna)$/, '');
        const dom = domLookup[`${u}::${name}`] || 'Unknown';
        if (!seriesMap[dom]) {
          seriesMap[dom] = { name: dom, type: 'scatter', symbolSize: 9, color: domColors[dom],
            itemStyle: { opacity: 0.8, borderColor: 'white', borderWidth: 1 }, data: [] };
        }
        seriesMap[dom].data.push({ value: [cont, comp], name: `${name} (${u})` });
      });
    });
    const zoneSeries = {
      type: 'scatter', data: [],
      markArea: { silent: true, data: [
        [{ coord: [-0.5, 90], itemStyle: { color: 'rgba(22,163,74,0.08)' } }, { coord: [5, 106] }],
        [{ coord: [-0.5, 50], itemStyle: { color: 'rgba(217,119,6,0.07)' } }, { coord: [10, 90] }],
        [{ coord: [-0.5, -3], itemStyle: { color: 'rgba(239,68,68,0.05)' } }, { coord: [15, 50] }],
      ] },
      markLine: { silent: true, data: [
        { yAxis: 90, lineStyle: { type: 'dashed', color: '#16a34a' }, label: { formatter: '≥90% HQ' } },
        { yAxis: 50, lineStyle: { type: 'dotted', color: '#d97706' }, label: { formatter: '≥50% MQ' } },
        { xAxis: 5,  lineStyle: { type: 'dashed', color: '#ef4444' }, label: { formatter: '≤5% cont' } },
        { xAxis: 10, lineStyle: { type: 'dotted', color: '#d97706' }, label: { formatter: '≤10% MQ' } },
      ] },
    };
    mkChart('coas-checkm2-chart', {
      title: { text: `MAG Quality — Completeness vs Contamination (${label})` },
      tooltip: { trigger: 'item', formatter: p => `${p.data.name}<br>Comp: ${(p.data.value[1] || 0).toFixed(1)}%<br>Cont: ${(p.data.value[0] || 0).toFixed(1)}%` },
      legend: { data: Object.keys(seriesMap) },
      xAxis: { type: 'value', name: 'Contamination (%)', min: -0.5, max: 15, nameLocation: 'middle', nameGap: 30 },
      yAxis: { type: 'value', name: 'Completeness (%)', min: -2, max: 105, nameLocation: 'middle', nameGap: 35 },
      series: [...Object.values(seriesMap), zoneSeries],
    });

    // GTDB taxonomy: domain bar + top phyla + master table
    const mp = isAll ? merged : merged.filter(r => r.sample === unit);

    const domCount = {};
    mp.forEach(r => { const d = r.Domain || 'Unknown'; domCount[d] = (domCount[d] || 0) + 1; });
    const domEntries = Object.entries(domCount);
    mkChart('coas-prok-domain-chart', {
      title: { text: `${label} — Domain Distribution` },
      tooltip: { trigger: 'axis' },
      xAxis: { type: 'category', data: domEntries.map(x => x[0]) },
      yAxis: { type: 'value', name: 'MAGs' },
      series: [{ type: 'bar', data: domEntries.map(x => x[1]), itemStyle: { color: '#0d9488' } }],
    });

    const phylaCount = {};
    mp.forEach(r => { const p = r.Phylum || 'Unknown'; phylaCount[p] = (phylaCount[p] || 0) + 1; });
    const topPhyla = Object.entries(phylaCount).sort((a, b) => b[1] - a[1]).slice(0, 20);
    mkChart('coas-prok-phyla-chart', {
      title: { text: `${label} — Top Phyla` },
      tooltip: { trigger: 'axis' },
      xAxis: { type: 'value', name: 'MAGs', nameLocation: 'middle', nameGap: 28 },
      yAxis: { type: 'category', data: topPhyla.map(x => x[0]).reverse(), axisLabel: { width: 160, overflow: 'truncate' } },
      series: [{ type: 'bar', data: topPhyla.map(x => x[1]).reverse(), itemStyle: { color: '#d97706' } }],
      grid: { left: 180, right: 30, bottom: 50 },
    });

    const tableCard = document.querySelector('#coas-prok-tax-table')?.closest('.chart-card');
    if (tableCard) tableCard.style.display = isAll ? 'none' : '';
    if (!isAll) {
      makeTable('coas-prok-tax-table', mp, [
        { key: 'Bin', label: 'Bin' },
        { key: 'Domain', label: 'Domain' },
        { key: 'Phylum', label: 'Phylum' },
        { key: 'Genus', label: 'Genus' },
        { key: 'Source_tax', label: 'Source' },
        { key: 'Completeness', label: 'Comp %' },
        { key: 'Contamination', label: 'Cont %' },
      ], { searchId: 'coas-prok-tax-search' });
    }
  }

  // ── helpers ─────────────────────────────────────────────────────────────────
  function _concat(dict, units) {
    const out = [];
    (units || []).forEach(u => ((dict || {})[u] || []).forEach(r => out.push(r)));
    return out;
  }
  function _tierCounts(rows) {
    const c = {};
    rows.forEach(r => { const q = r.checkv_quality || 'Not-determined'; c[q] = (c[q] || 0) + 1; });
    return c;
  }
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

})();
