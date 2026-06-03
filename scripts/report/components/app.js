/* app.js — Tab switching, theme toggle, sub-tab switching, sample selectors */
(function () {
  'use strict';

  // ── ECharts palette (matches design system) ───────────────────────────────
  window.PAL = [
    '#0d9488','#d97706','#7c3aed','#0891b2','#16a34a',
    '#f59e0b','#9333ea','#ef4444','#0e7490','#b45309',
    '#64748b','#2563eb','#db2777','#059669','#ca8a04',
  ];

  // ── ECharts base theme ────────────────────────────────────────────────────
  function echartsTheme() {
    const dark = document.documentElement.dataset.theme === 'dark';
    return {
      backgroundColor: 'transparent',
      textStyle:  { color: dark ? '#f1f5f9' : '#0f172a', fontFamily: 'Inter, system-ui, sans-serif', fontSize: 12 },
      title:      { top: 4, textStyle: { fontSize: 13, fontWeight: 600, color: dark ? '#f1f5f9' : '#0f172a' } },
      legend:     { top: 26, textStyle: { color: dark ? '#94a3b8' : '#64748b', fontSize: 11 } },
      tooltip:    {
        backgroundColor: dark ? '#1e293b' : '#fff',
        borderColor:     dark ? '#334155' : '#e2e8f0',
        textStyle:       { color: dark ? '#f1f5f9' : '#0f172a', fontSize: 12 },
      },
      xAxis:      { axisLine: { lineStyle: { color: dark ? '#334155' : '#e2e8f0' } },
                    splitLine: { lineStyle: { color: dark ? '#1e293b' : '#f1f5f9' } },
                    axisLabel: {
                      color: dark ? '#94a3b8' : '#64748b',
                      // Truncate long category labels (e.g. sample names); tooltip still shows full name
                      formatter: val => String(val).length > 18 ? String(val).slice(0, 17) + '…' : val,
                    } },
      yAxis:      { axisLine: { lineStyle: { color: dark ? '#334155' : '#e2e8f0' } },
                    splitLine: { lineStyle: { color: dark ? '#1e293b' : '#f1f5f9' } },
                    axisLabel: { color: dark ? '#94a3b8' : '#64748b' },
                    nameLocation: 'middle', nameGap: 42,
                    nameTextStyle: { color: dark ? '#94a3b8' : '#64748b', fontSize: 11 } },
      color: window.PAL,
    };
  }
  window.echartsTheme = echartsTheme;

  // Deep-merge axis: preserves theme axisLabel defaults (formatter, color) while
  // allowing per-chart overrides (rotate, width, etc.) to be layered on top.
  function _mergeAxis(themeAxis, optAxis) {
    const merged = { ...themeAxis, ...optAxis };
    if (themeAxis.axisLabel || optAxis.axisLabel) {
      merged.axisLabel = { ...(themeAxis.axisLabel || {}), ...(optAxis.axisLabel || {}) };
    }
    if (themeAxis.nameTextStyle || optAxis.nameTextStyle) {
      merged.nameTextStyle = { ...(themeAxis.nameTextStyle || {}), ...(optAxis.nameTextStyle || {}) };
    }
    return merged;
  }

  // ── Create/update an ECharts instance ────────────────────────────────────
  window._charts = {};
  window.mkChart = function (id, option) {
    const el = document.getElementById(id);
    if (!el) return null;
    let chart = window._charts[id];
    if (!chart) {
      chart = echarts.init(el, null, { renderer: 'canvas' });
      window._charts[id] = chart;
      const ro = new ResizeObserver(() => chart.resize());
      ro.observe(el);
    }
    const theme = echartsTheme();
    const merged = {
      ...theme,
      tooltip: { ...theme.tooltip, trigger: 'axis', ...(option.tooltip || {}) },
      legend:  { ...theme.legend,  ...(option.legend  || {}) },
      xAxis:   option.xAxis ? (Array.isArray(option.xAxis)
                  ? option.xAxis.map(a => _mergeAxis(theme.xAxis, a))
                  : _mergeAxis(theme.xAxis, option.xAxis)) : undefined,
      yAxis:   option.yAxis ? (Array.isArray(option.yAxis)
                  ? option.yAxis.map(a => _mergeAxis(theme.yAxis, a))
                  : _mergeAxis(theme.yAxis, option.yAxis)) : undefined,
      color:   theme.color,
    };
    // merge remaining keys (series, grid, etc.) — grid gets a default top to clear title+legend
    const hasLegend = !!(option.legend && (option.legend.data || option.legend.show !== false));
    const defaultGridTop = hasLegend ? 58 : 36;
    Object.keys(option).forEach(k => {
      if (!['tooltip','legend','xAxis','yAxis','color'].includes(k)) {
        merged[k] = option[k];
      }
    });
    if (!merged.grid) {
      merged.grid = { top: defaultGridTop, bottom: 40, left: 50, right: 20, containLabel: true };
    } else if (!Array.isArray(merged.grid) && merged.grid.top === undefined) {
      merged.grid = { top: defaultGridTop, bottom: 40, left: 50, right: 20, containLabel: true, ...merged.grid };
    }
    chart.setOption(merged, true);
    return chart;
  };

  // Re-render all ECharts on theme change
  function refreshCharts() {
    Object.values(window._charts).forEach(c => {
      if (c && !c.isDisposed()) c.setOption(echartsTheme(), false);
    });
  }

  // ── Theme toggle ──────────────────────────────────────────────────────────
  const html  = document.documentElement;
  const btn   = document.getElementById('theme-toggle');
  function applyTheme(dark) {
    html.dataset.theme = dark ? 'dark' : 'light';
    btn.textContent    = dark ? '☾' : '☀';
    localStorage.setItem('vapor-theme', dark ? 'dark' : 'light');
    refreshCharts();
    if (typeof renderPCoA === 'function') renderPCoA();
    if (typeof renderVC3Network === 'function') renderVC3Network(_currentVC3Sample);
  }
  btn.addEventListener('click', () => applyTheme(html.dataset.theme !== 'dark'));
  // Restore saved preference
  applyTheme(localStorage.getItem('vapor-theme') === 'dark');

  // ── Main tab switching ────────────────────────────────────────────────────
  const tabs   = document.querySelectorAll('.nav-tab');
  const panels = document.querySelectorAll('.tab-panel');

  function showTab(name) {
    tabs.forEach(t => {
      const sel = t.dataset.tab === name;
      t.setAttribute('aria-selected', String(sel));
      t.tabIndex = sel ? 0 : -1;
    });
    panels.forEach(p => {
      const active = p.id === `tab-${name}`;
      p.classList.toggle('active', active);
      p.hidden = !active;
    });
    // Lazy-render on first show
    const evt = new CustomEvent('vapor:tabshow', { detail: { tab: name } });
    document.dispatchEvent(evt);
  }

  tabs.forEach(t => {
    t.addEventListener('click', () => showTab(t.dataset.tab));
    t.addEventListener('keydown', e => {
      if (e.key === 'ArrowRight') { t.nextElementSibling?.click(); }
      if (e.key === 'ArrowLeft')  { t.previousElementSibling?.click(); }
    });
  });

  // ── Sub-tab switching ─────────────────────────────────────────────────────
  document.querySelectorAll('.sub-nav').forEach(nav => {
    const stabs   = nav.querySelectorAll('.sub-tab');
    const section = nav.closest('.tab-panel');
    stabs.forEach(st => {
      st.addEventListener('click', () => {
        stabs.forEach(x => x.classList.remove('active'));
        st.classList.add('active');
        const target = st.dataset.sub;
        section.querySelectorAll('.sub-panel').forEach(p => {
          const show = p.id === `sub-${target}`;
          p.style.display = show ? '' : 'none';
          p.classList.toggle('active', show);
        });
        // Trigger resize for ECharts inside newly visible panel
        setTimeout(() => Object.values(window._charts).forEach(c => c && !c.isDisposed() && c.resize()), 80);
        // Notify components so they can re-render D3 visuals that need visible dimensions
        document.dispatchEvent(new CustomEvent('vapor:subtabshow', { detail: { sub: target } }));
      });
    });
  });

  // ── Generic table builder ─────────────────────────────────────────────────
  window.makeTable = function (containerId, rows, cols, opts = {}) {
    const el = document.getElementById(containerId);
    if (!el || !rows.length) { if (el) el.innerHTML = '<p style="color:var(--text-muted);font-size:.8rem;padding:.5rem">No data</p>'; return; }
    const thead = cols.map(c => `<th>${c.label || c.key}</th>`).join('');
    const tbody = rows.map(r =>
      `<tr>${cols.map(c => `<td>${opts.format?.[c.key]?.(r[c.key], r) ?? (r[c.key] ?? '')}</td>`).join('')}</tr>`
    ).join('');
    el.innerHTML = `<table class="vapor-table"><thead><tr>${thead}</tr></thead><tbody>${tbody}</tbody></table>`;
    // Search
    if (opts.searchId) {
      document.getElementById(opts.searchId)?.addEventListener('input', function () {
        const q = this.value.toLowerCase();
        el.querySelectorAll('tbody tr').forEach(tr => {
          tr.style.display = tr.textContent.toLowerCase().includes(q) ? '' : 'none';
        });
      });
    }
  };

  // ── Quality badge helper ──────────────────────────────────────────────────
  window.qualBadge = function (q) {
    const map = {
      'Complete':        'badge-green',
      'High-quality':    'badge-teal',
      'Medium-quality':  'badge-amber',
      'Low-quality':     'badge-red',
      'Not-determined':  'badge-gray',
    };
    return `<span class="badge ${map[q] || 'badge-gray'}">${q || 'n/a'}</span>`;
  };

  // ── Sample dropdown builder ───────────────────────────────────────────────
  window.makeSampleDropdown = function (selId, onChange) {
    const sel = document.getElementById(selId);
    if (!sel) return;
    sel.innerHTML = '';
    (typeof SAMPLES !== 'undefined' ? SAMPLES : []).forEach(s => {
      const opt = document.createElement('option');
      opt.value = opt.textContent = s;
      sel.appendChild(opt);
    });
    sel.addEventListener('change', () => onChange(sel.value));
    if (sel.options.length) onChange(sel.value);
  };

  // ── Kick off all renders ──────────────────────────────────────────────────
  document.addEventListener('DOMContentLoaded', () => {
    showTab('overview');
    if (typeof renderOverview    === 'function') renderOverview();
    if (typeof renderSequencing  === 'function') renderSequencing();
    if (typeof renderViral       === 'function') renderViral();
    if (typeof renderProkaryotic === 'function') renderProkaryotic();
    if (typeof renderDiversity   === 'function') renderDiversity();
    if (typeof renderAnnotation  === 'function') renderAnnotation();
    if (typeof renderAbout       === 'function') renderAbout();
  });

})();
