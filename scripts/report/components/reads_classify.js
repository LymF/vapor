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
    // sylphmpa uses | separator; viral taxonomy uses r__ (realm), k__ (kingdom)
    const map = {
      family: 'f__', order: 'o__', genus: 'g__', phylum: 'p__', class: 'c__',
      species: 's__', realm: 'r__', kingdom: 'k__', domain: 'd__',
    };
    const prefix = map[rank] || 'f__';
    const parts = (clade || '').split('|');
    for (const p of parts) {
      if (p.startsWith(prefix)) {
        const name = p.slice(3).trim();
        return name || 'Unclassified';
      }
    }
    // If the rank prefix doesn't exist but there are UNKNOWN entries at that depth,
    // fall back to 'Unclassified' so viral rows with partial taxonomy still group.
    if (parts.some(p => p === 'UNKNOWN')) return 'Unclassified';
    return '';
  }

  // Standard rank prefixes in order shallowest→deepest (mirrors Python _STD_RANKS)
  const STD_RANK_PREFIXES = ['r__', 'k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__'];

  // Map rank name → prefix
  const RANK_PREFIX = {
    realm: 'r__', kingdom: 'k__', phylum: 'p__', class: 'c__',
    order: 'o__', family: 'f__', genus: 'g__', species: 's__',
  };

  /** Derive effective rank prefix from clade string (mirrors Python _last_std_rank_seg). */
  function _effRank(clade) {
    let last = '';
    for (const seg of (clade || '').split('|')) {
      if (STD_RANK_PREFIXES.some(p => seg.startsWith(p))) last = seg.slice(0, 3);
    }
    return last;
  }

  function _aggregateByRank(rows, rank, topN) {
    // Use only rows at the requested rank level (by recomputing effective rank from
    // the clade string — more robust than relying on the serialized _eff_rank field).
    // Rows at a shallower or deeper rank are skipped to avoid double-counting.
    const prefix = RANK_PREFIX[rank] || 'f__';
    const samples = _samples();
    const agg = {};
    for (const r of rows) {
      if (_effRank(r.clade) !== prefix) continue;
      const label = _rankField(r.clade, rank) || 'Unclassified';
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
    if (!groups.length) {
      el.innerHTML = '<p class="muted" style="padding:20px;text-align:center">No data at this rank level</p>';
      return;
    }
    // Fold to the validated 7 hues + "Other". `groups` arrives sorted by total
    // abundance and could be 15 long -- stacking those cycled PAL, so taxon 1
    // and taxon 9 were painted identically while the legend called them
    // different. Routed through samplesBar so it also stays readable past
    // VIZ.manySamples (REPORT_VIZ_GUIDE §4).
    const keep = groups.slice(0, window.VIZ.maxSeries - 1);
    const rest = groups.slice(window.VIZ.maxSeries - 1);
    const series = keep.map((g, i) => ({
      name: g.label, color: window.PAL[i],
      data: samples.map(s => +(g[s] || 0).toFixed(5)),
    }));
    if (rest.length) {
      const otherData = samples.map(s => +rest.reduce((a, g) => a + (g[s] || 0), 0).toFixed(5));
      if (otherData.some(v => v > 0)) {
        series.push({ name: `Other (${rest.length})`, color: window.PAL_MUTED, data: otherData });
      }
    }
    window.mkChart(chartId, window.samplesBar({
      samples, series, title, stack: true,
      valueName: 'Relative abundance (%)',
    }));
  }

  function _barChart(chartId, labels, values, title, color) {
    const el = document.getElementById(chartId);
    if (!el) return;
    window.mkChart(chartId, {
      title: { text: title, textStyle: { fontSize: 13 } },
      tooltip: { trigger: 'axis' },
      grid: { left: 120 },
      xAxis: { type: 'value', axisLabel: { formatter: v => `${(+v).toFixed(2)}` } },
      yAxis: { type: 'category', data: labels, axisLabel: { fontSize: 10 } },
      series: [{ type: 'bar', data: values, itemStyle: { color: color || window.PAL[0] } }],
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

  /* ── Viral ── */
  function _renderViral() {
    const viral = (RC.viral || []);
    const samples = _samples();
    const rank = (document.getElementById('rc-viral-rank-sel') || {}).value || 'class';
    // values are in % scale (0-100) — threshold is applied directly in %
    const minPct = parseFloat((document.getElementById('rc-viral-min-abund') || {}).value || 0);

    const filtered = minPct > 0 ? viral.filter(r => samples.some(s => (r[s] || 0) >= minPct)) : viral;
    const groups = _aggregateByRank(filtered, rank, 15);
    _stackedBar('rc-viral-family-chart', groups, samples, `Top Viral ${rank.charAt(0).toUpperCase()+rank.slice(1)}s — Relative Abundance`);

    // Richness (detected taxa per sample)
    const richness = samples.map(s => viral.filter(r => (r[s] || 0) > 0).length);
    _barChart('rc-viral-richness-chart', samples, richness, 'Detected Viral Taxa per Sample', window.PAL[4]);

    // Class breakdown — most environmental viromes are Caudoviricetes-dominated
    const classes = _aggregateByRank(viral, 'class', 8);
    _stackedBar('rc-viral-domain-chart', classes, samples, 'Viral Classes');

    // Table: top taxa
    // Table: all effective-rank viral rows sorted by total abundance
    const tableRows = (RC.viral || [])
      .map(r => ({
        Class:  _rankField(r.clade, 'class')  || '—',
        Order:  _rankField(r.clade, 'order')  || '—',
        Family: _rankField(r.clade, 'family') || '—',
        Genus:  _rankField(r.clade, 'genus')  || '—',
        Rank: r._eff_rank || '—',
        ...Object.fromEntries(samples.map(s => [s, `${(+(r[s]||0)).toFixed(4)}%`])),
      }))
      .sort((a, b) => {
        const sa = samples.reduce((t,s) => t + parseFloat(a[s]||0), 0);
        const sb = samples.reduce((t,s) => t + parseFloat(b[s]||0), 0);
        return sb - sa;
      })
      .slice(0, 200);
    const tableCols = ['Rank', 'Class', 'Order', 'Family', 'Genus', ...samples];
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
    _barChart('rc-prok-richness-chart', samples, richness, 'Prokaryotic Taxa Detected per Sample', window.PAL_MUTED);

    // Table
    const tableRows = prok
      .map(r => ({
        Phylum: _rankField(r.clade, 'phylum') || '—',
        Class:  _rankField(r.clade, 'class')  || '—',
        Family: _rankField(r.clade, 'family') || '—',
        Genus:  _rankField(r.clade, 'genus')  || '—',
        Domain: _rankField(r.clade, 'domain') || '—',
        ...Object.fromEntries(samples.map(s => [s, `${(+(r[s]||0)).toFixed(4)}%`])),
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
      'Mean Viral Relative Abundance by Predicted Host Genus', window.PAL[1]);

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
