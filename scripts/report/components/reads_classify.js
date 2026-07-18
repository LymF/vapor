/* reads_classify.js — Reads Survey tab (Sylph reads-only classification) */
(function () {
  'use strict';

  const RC = typeof READS_CLASSIFY !== 'undefined' ? READS_CLASSIFY : null;

  window.renderReadsClassify = function () {
    if (!RC || !RC.has_data) {
      const el = document.getElementById('rc-no-data');
      const ct = document.getElementById('rc-viral-content');
      if (el) el.style.display = '';
      if (ct) ct.style.display = 'none';
      return;
    }
    _renderViral();
    _renderProk();
    _renderHost();
    _wireControls();
  };

  /* ── helpers ── */
  function _samples() { return (RC && RC.samples) || []; }

  function _rankField(clade, rank) {
    const map = { family: 'f__', order: 'o__', genus: 'g__', phylum: 'p__', class: 'c__' };
    const prefix = map[rank] || 'f__';
    const parts = (clade || '').split(';');
    for (const p of parts) if (p.startsWith(prefix)) return p.slice(3).trim();
    return '';
  }

  function _aggregateByRank(rows, rank, topN) {
    const samples = _samples();
    const agg = {};
    for (const r of rows) {
      const label = _rankField(r.clade, rank) || 'Unknown';
      if (!agg[label]) { agg[label] = { label }; for (const s of samples) agg[label][s] = 0; }
      for (const s of samples) agg[label][s] += (r[s] || 0);
    }
    const sorted = Object.values(agg).sort((a, b) => {
      const sumA = samples.reduce((t, s) => t + (a[s] || 0), 0);
      const sumB = samples.reduce((t, s) => t + (b[s] || 0), 0);
      return sumB - sumA;
    });
    return sorted.slice(0, topN || 15);
  }

  function _stackedBar(chartId, groups, samples, title) {
    const el = document.getElementById(chartId);
    if (!el) return;
    const chart = echarts.init(el, _theme());
    const COLORS = ['#5b8ff9','#5ad8a6','#5d7092','#f6bd16','#e8684a','#6dc8ec',
                    '#9867bc','#8d684b','#f2b4b8','#87e8de','#ffd591','#b7eb8f'];
    const series = groups.map((g, i) => ({
      name: g.label,
      type: 'bar',
      stack: 'total',
      data: samples.map(s => +(g[s] || 0).toFixed(5)),
      itemStyle: { color: COLORS[i % COLORS.length] },
    }));
    chart.setOption({
      title: { text: title, textStyle: { fontSize: 13 } },
      tooltip: { trigger: 'axis', axisPointer: { type: 'shadow' },
        formatter: params => params.map(p => `${p.marker}${p.seriesName}: ${(p.value*100).toFixed(3)}%`).join('<br>') },
      legend: { type: 'scroll', bottom: 0, textStyle: { fontSize: 11 } },
      grid: { top: 40, bottom: 60, left: 80, right: 20 },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30, fontSize: 10 } },
      yAxis: { type: 'value', name: 'Relative abundance', axisLabel: { formatter: v => `${(v*100).toFixed(1)}%` } },
      series,
    });
  }

  function _barChart(chartId, labels, values, title, color) {
    const el = document.getElementById(chartId);
    if (!el) return;
    const chart = echarts.init(el, _theme());
    chart.setOption({
      title: { text: title, textStyle: { fontSize: 13 } },
      tooltip: { trigger: 'axis' },
      grid: { left: 120, right: 20, top: 40, bottom: 40 },
      xAxis: { type: 'value', axisLabel: { formatter: v => `${(v*100).toFixed(2)}%` } },
      yAxis: { type: 'category', data: labels, axisLabel: { fontSize: 10 } },
      series: [{ type: 'bar', data: values, itemStyle: { color: color || '#5b8ff9' } }],
    });
  }

  function _makeTable(containerId, rows, cols, searchId) {
    const el = document.getElementById(containerId);
    if (!el || !rows.length) { if (el) el.innerHTML = '<p class="muted">No data</p>'; return; }
    const thead = `<thead><tr>${cols.map(c => `<th>${c}</th>`).join('')}</tr></thead>`;
    const renderRows = data => data.map(r => `<tr>${cols.map(c => `<td>${r[c] ?? ''}</td>`).join('')}</tr>`).join('');
    el.innerHTML = `<table class="vapor-table">${thead}<tbody id="${containerId}-body">${renderRows(rows)}</tbody></table>`;
    if (searchId) {
      const inp = document.getElementById(searchId);
      if (inp) inp.addEventListener('input', () => {
        const q = inp.value.toLowerCase();
        const filtered = rows.filter(r => cols.some(c => String(r[c] || '').toLowerCase().includes(q)));
        const body = document.getElementById(`${containerId}-body`);
        if (body) body.innerHTML = renderRows(filtered);
      });
    }
  }

  function _theme() { return document.documentElement.dataset.theme === 'dark' ? 'dark' : 'light'; }

  /* ── Viral ── */
  function _renderViral() {
    const viral = (RC.viral || []);
    const samples = _samples();
    const rank = (document.getElementById('rc-viral-rank-sel') || {}).value || 'family';
    const minPct = parseFloat((document.getElementById('rc-viral-min-abund') || {}).value || 0) / 100;

    const filtered = minPct > 0 ? viral.filter(r => samples.some(s => (r[s] || 0) >= minPct)) : viral;
    const groups = _aggregateByRank(filtered, rank, 15);
    _stackedBar('rc-viral-family-chart', groups, samples, `Top Viral ${rank.charAt(0).toUpperCase()+rank.slice(1)}s — Relative Abundance`);

    // Richness (detected taxa per sample)
    const richness = samples.map(s => viral.filter(r => (r[s] || 0) > 0).length);
    _barChart('rc-viral-richness-chart', samples, richness, 'Detected Viral Taxa per Sample', '#5ad8a6');

    // Domain breakdown (always Viruses here, but sub-realm if available)
    const realms = _aggregateByRank(viral, 'order', 8);
    _stackedBar('rc-viral-domain-chart', realms, samples, 'Viral Orders');

    // Table: top taxa
    const tableRows = filtered
      .map(r => ({
        Family: _rankField(r.clade, 'family') || '—',
        Order:  _rankField(r.clade, 'order')  || '—',
        Genus:  _rankField(r.clade, 'genus')  || '—',
        Species: _rankField(r.clade, 'species') || '—',
        Clade: r.clade,
        ...Object.fromEntries(samples.map(s => [s, `${((r[s]||0)*100).toFixed(4)}%`])),
      }))
      .sort((a, b) => {
        const sa = samples.reduce((t,s) => t + parseFloat(a[s]||0), 0);
        const sb = samples.reduce((t,s) => t + parseFloat(b[s]||0), 0);
        return sb - sa;
      })
      .slice(0, 200);
    const tableCols = ['Family', 'Order', 'Genus', 'Species', ...samples];
    _makeTable('rc-viral-table', tableRows, tableCols, 'rc-viral-search');
  }

  /* ── Prokaryotic ── */
  function _renderProk() {
    const prok = (RC.prok || []).concat(RC.archaea || []);
    const samples = _samples();
    const rank = (document.getElementById('rc-prok-rank-sel') || {}).value || 'phylum';

    const groups = _aggregateByRank(prok, rank, 15);
    _stackedBar('rc-prok-phylum-chart', groups, samples, `Top Prokaryotic ${rank.charAt(0).toUpperCase()+rank.slice(1)}s`);

    // Domain (Bacteria vs Archaea)
    const bact = _aggregateByRank(RC.prok || [], 'phylum', 8);
    const arch = _aggregateByRank(RC.archaea || [], 'phylum', 8);
    const domGroups = [
      ...bact.map(g => ({ ...g, label: `Bacteria — ${g.label}` })),
      ...arch.map(g => ({ ...g, label: `Archaea — ${g.label}` })),
    ];
    _stackedBar('rc-prok-domain-chart', domGroups, samples, 'Bacteria vs. Archaea');

    // Richness
    const richness = samples.map(s => prok.filter(r => (r[s]||0) > 0).length);
    _barChart('rc-prok-richness-chart', samples, richness, 'Prokaryotic Taxa Detected per Sample', '#5d7092');

    // Table
    const tableRows = prok
      .map(r => ({
        Phylum: _rankField(r.clade, 'phylum') || '—',
        Class:  _rankField(r.clade, 'class')  || '—',
        Family: _rankField(r.clade, 'family') || '—',
        Genus:  _rankField(r.clade, 'genus')  || '—',
        Domain: _rankField(r.clade, 'domain') || '—',
        ...Object.fromEntries(samples.map(s => [s, `${((r[s]||0)*100).toFixed(4)}%`])),
      }))
      .sort((a, b) => {
        const sa = samples.reduce((t,s) => t + parseFloat(a[s]||0), 0);
        const sb = samples.reduce((t,s) => t + parseFloat(b[s]||0), 0);
        return sb - sa;
      })
      .slice(0, 200);
    const tableCols = ['Domain', 'Phylum', 'Class', 'Family', 'Genus', ...samples];
    _makeTable('rc-prok-table', tableRows, tableCols, 'rc-prok-search');
  }

  /* ── Host ── */
  function _renderHost() {
    const host = RC.host || [];
    const samples = _samples();
    if (!host.length) {
      const el = document.getElementById('rc-host-chart');
      if (el) el.innerHTML = '<p class="muted" style="padding:20px">Host annotation not available — requires IMGVR or UHGV pre-built database with sylph-tax -a flag.</p>';
      return;
    }

    const top = host.slice(0, 20);
    const labels = top.map(r => r.host_genus || 'Unknown');
    const values = top.map(r => samples.reduce((t, s) => t + (r[s] || 0), 0) / samples.length);
    _barChart('rc-host-chart', labels.reverse(), values.reverse(),
      'Mean Viral Relative Abundance by Predicted Host Genus', '#e8684a');

    const tableRows = host.map(r => ({
      'Host Genus': r.host_genus || 'Unknown',
      'Viral Taxa': r.n_viral_taxa || '',
      ...Object.fromEntries(samples.map(s => [s, `${((r[s]||0)*100).toFixed(4)}%`])),
    }));
    _makeTable('rc-host-table', tableRows, ['Host Genus', 'Viral Taxa', ...samples]);
  }

  /* ── Controls ── */
  function _wireControls() {
    const rankSel = document.getElementById('rc-viral-rank-sel');
    const minInp  = document.getElementById('rc-viral-min-abund');
    if (rankSel) rankSel.addEventListener('change', _renderViral);
    if (minInp)  minInp.addEventListener('change', _renderViral);
    const prokRank = document.getElementById('rc-prok-rank-sel');
    if (prokRank) prokRank.addEventListener('change', _renderProk);
  }

})();
