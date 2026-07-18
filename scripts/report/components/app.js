/* app.js — Tab switching, theme toggle, sub-tab switching, sample selectors */
(function () {
  'use strict';

  // ── ECharts palette (matches design system) ───────────────────────────────
  // 8 hues, fixed order -- validated with the dataviz-skill validator
  // (scripts/validate_palette.js, light mode): lightness band, chroma floor,
  // CVD adjacent-pair separation and surface contrast all PASS for exactly
  // these 8, in this order. Never index past PAL.length -- a 9th category
  // must fold into "Other" (see window.foldOther) instead of generating or
  // cycling a hue; a validated 9-15 extension was tried and FAILED (two
  // slots below the chroma floor, one adjacent pair at CVD ΔE 2.9).
  window.PAL = [
    '#0d9488','#d97706','#7c3aed','#0891b2',
    '#16a34a','#f59e0b','#9333ea','#ef4444',
  ];
  // Neutral gray for "Other"/"Unknown" buckets -- not a categorical identity
  // slot (never used for a real series), just an escape hatch for the tail.
  window.PAL_MUTED = '#64748b';

  // Fold a {name: count} map down to the top `max` entries + one "Other"
  // bucket summing the rest, sorted descending. Keeps stacked/legend series
  // within the validated 8-hue budget regardless of how many distinct
  // categories the underlying data has (e.g. COG/PHROGS functional letters).
  window.foldOther = function (counts, max = 7) {
    const sorted = Object.entries(counts).sort((a, b) => b[1] - a[1]);
    const top = sorted.slice(0, max);
    const rest = sorted.slice(max).reduce((a, [, v]) => a + v, 0);
    if (rest > 0) top.push(['Other', rest]);
    return top;
  };

  // Single-hue, monotone-lightness ramp for ORDINAL data (pipeline stages,
  // size tiers) -- never use categorical PAL slots for a sequence where
  // order carries meaning. Both ramps validated with --ordinal (lightness
  // monotone, adjacent |dL| >= 0.06, light end clears contrast floor for its
  // surface); dark-mode steps are shifted lighter as a set, not reversed,
  // since the floor is checked against the dark card surface (#1e293b).
  const _RAMP_LIGHT = ['#14b8a6','#0d9488','#0f766e','#134e4a','#042f2e'];
  const _RAMP_DARK  = ['#5eead4','#2dd4bf','#14b8a6','#0d9488','#0f766e'];
  window.ordinalRamp = function (n) {
    const dark = document.documentElement.dataset.theme === 'dark';
    const ramp = dark ? _RAMP_DARK : _RAMP_LIGHT;
    if (n === ramp.length) return ramp;
    // Even subsample for n < 5; ramps this codebase uses today never exceed 5.
    return Array.from({ length: n }, (_, i) => ramp[Math.round(i * (ramp.length - 1) / Math.max(1, n - 1))]);
  };

  // ── Boxplot stats helper (median/IQR/whiskers + 1.5*IQR outliers) ──────────
  // Used by per-sample distribution charts (contig length/depth, genome size,
  // completeness/contamination) to avoid one legend entry per sample.
  window.boxStats = function (values) {
    const sorted = (values || []).slice().sort((a, b) => a - b);
    if (!sorted.length) return { box: [0, 0, 0, 0, 0], outliers: [] };
    const q = p => {
      const pos  = (sorted.length - 1) * p;
      const base = Math.floor(pos);
      const rest = pos - base;
      return sorted[base + 1] !== undefined
        ? sorted[base] + rest * (sorted[base + 1] - sorted[base])
        : sorted[base];
    };
    const q1 = q(0.25), med = q(0.5), q3 = q(0.75);
    const iqr = q3 - q1;
    const lo = q1 - 1.5 * iqr, hi = q3 + 1.5 * iqr;
    const within   = sorted.filter(v => v >= lo && v <= hi);
    const outliers = sorted.filter(v => v < lo || v > hi);
    return {
      box: [within[0] ?? sorted[0], q1, med, q3, within[within.length - 1] ?? sorted[sorted.length - 1]],
      outliers,
    };
  };

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
  }
  btn.addEventListener('click', () => applyTheme(html.dataset.theme !== 'dark'));
  // Restore saved preference
  applyTheme(localStorage.getItem('vapor-theme') === 'dark');

  // ── Back to top ────────────────────────────────────────────────────────────
  const backToTop = document.getElementById('back-to-top');
  if (backToTop) {
    window.addEventListener('scroll', () => {
      backToTop.classList.toggle('visible', window.scrollY > 400);
    });
    backToTop.addEventListener('click', () => window.scrollTo({ top: 0, behavior: 'smooth' }));
  }

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
  // Renders in pages instead of dumping every row into the DOM at once --
  // cross-linked tables (host<->virus<->defense<->AMR) can reach thousands
  // of rows across a multi-sample run and freeze the tab on first paint
  // otherwise. Search re-filters the full row set (not just the rendered
  // page) and restarts pagination from page 1.
  window.makeTable = function (containerId, rows, cols, opts = {}) {
    const el = document.getElementById(containerId);
    if (!el) return;
    const pageSize = opts.pageSize || 100;
    const rowText = r => cols.map(c => String(r[c.key] ?? '')).join(' ').toLowerCase();
    const state = { all: rows || [], filtered: rows || [], shown: 0 };

    function renderRows(start, end) {
      return state.filtered.slice(start, end).map(r =>
        `<tr>${cols.map(c => `<td>${opts.format?.[c.key]?.(r[c.key], r) ?? (r[c.key] ?? '')}</td>`).join('')}</tr>`
      ).join('');
    }

    function renderMoreButton() {
      el.querySelector(':scope > .vapor-table-more')?.remove();
      if (state.shown >= state.filtered.length) return;
      const btn = document.createElement('button');
      btn.type = 'button';
      btn.className = 'vapor-table-more';
      btn.textContent = `Load more (${state.shown} of ${state.filtered.length})`;
      btn.addEventListener('click', () => {
        const start = state.shown;
        state.shown = Math.min(state.filtered.length, state.shown + pageSize);
        el.querySelector('tbody')?.insertAdjacentHTML('beforeend', renderRows(start, state.shown));
        renderMoreButton();
      });
      el.appendChild(btn);
    }

    function renderFromScratch() {
      if (!state.filtered.length) {
        el.innerHTML = '<p style="color:var(--text-muted);font-size:.8rem;padding:.5rem">No data</p>';
        return;
      }
      state.shown = Math.min(state.filtered.length, pageSize);
      const thead = cols.map(c => `<th>${c.label || c.key}</th>`).join('');
      el.innerHTML = `<table class="vapor-table"><thead><tr>${thead}</tr></thead>`
        + `<tbody>${renderRows(0, state.shown)}</tbody></table>`;
      renderMoreButton();
    }

    renderFromScratch();

    if (opts.searchId) {
      const input = document.getElementById(opts.searchId);
      if (input) {
        input.oninput = function () {
          const q = this.value.toLowerCase();
          state.filtered = q ? state.all.filter(r => rowText(r).includes(q)) : state.all;
          renderFromScratch();
        };
      }
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
  window.makeSampleDropdown = function (selId, onChange, opts) {
    const sel = document.getElementById(selId);
    if (!sel) return;
    sel.innerHTML = '';
    if (opts && opts.allSamples) {
      const opt = document.createElement('option');
      opt.value = '__all__';
      opt.textContent = 'All samples';
      sel.appendChild(opt);
    }
    (typeof SAMPLES !== 'undefined' ? SAMPLES : []).forEach(s => {
      const opt = document.createElement('option');
      opt.value = opt.textContent = s;
      sel.appendChild(opt);
    });
    sel.addEventListener('change', () => onChange(sel.value));
    if (sel.options.length) onChange(sel.value);
  };

  // ── Kick off all renders ──────────────────────────────────────────────────
  // Each render runs in its own try/catch so a bug in one tab (e.g. bad data
  // for a single chart) cannot prevent later tabs from rendering at all.
  document.addEventListener('DOMContentLoaded', () => {
    showTab('overview');
    [renderOverview, renderSequencing, renderViral, renderProkaryotic,
     renderHostDefense, renderDiversity, renderAnnotation, renderReadsClassify, renderAbout].forEach(fn => {
      if (typeof fn !== 'function') return;
      try {
        fn();
      } catch (e) {
        console.error(`[VAPOR] render error in ${fn.name}:`, e);
      }
    });
    // Export toolbars (PNG/SVG/PDF) on every chart/table card
    if (window.VaporExport) {
      window.VaporExport.injectAll();
      // Re-scan after lazy tab/sub-tab renders (D3 networks, etc.)
      document.addEventListener('vapor:tabshow',    () => setTimeout(() => window.VaporExport.injectAll(), 100));
      document.addEventListener('vapor:subtabshow', () => setTimeout(() => window.VaporExport.injectAll(), 100));
    }
  });

})();
