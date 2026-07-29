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
  // Validated with the palette checker in BOTH modes (see docs/REPORT_VIZ_GUIDE.md).
  // Slot 6 was '#f59e0b': it FAILED the dark-mode lightness band (L 0.769) and
  // sat at 2.09:1 contrast in light. '#db2777' passes light and dark; its
  // remaining CVD warning (dE 6.1 vs '#16a34a', deutan) is in the 6-8 floor
  // band, which is legal only because every categorical chart here also carries
  // a legend/direct labels and 2px gaps between fills.
  window.PAL = [
    '#0d9488','#d97706','#7c3aed','#0891b2',
    '#16a34a','#db2777','#9333ea','#ef4444',
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

  // ══ Data-shape helpers ═══════════════════════════════════════════════════
  // The numeric triggers from docs/REPORT_VIZ_GUIDE.md §4, in one place. Charts
  // adapt to the SHAPE of the data (how many samples, how many points) so the
  // report works on a 3-sample run and a 300-sample run without per-dataset
  // special cases. Change a threshold here and every chart follows.
  window.VIZ = {
    manySamples:  12,   // above this a per-sample category axis stops being readable
    densityMinN:  20,   // below this a KDE invents structure -> show the points
    manyGroups:    8,   // above this, distributions go to a ridgeline
    denseScatter: 500,  // above this a scatter overplots -> hexbin / 2-D density
    maxSeries:     8,   // validated categorical budget (see foldOther)
  };

  // Gaussian KDE on a fixed grid, Silverman bandwidth. Returns [[x, density], …]
  // normalised so the peak is 1 — ridgeline rows are shape comparisons, not
  // absolute-density comparisons, and per-row normalisation is what makes rows
  // with very different n readable side by side.
  function _kde(values, gridN, lo, hi) {
    const v = values.filter(x => Number.isFinite(x));
    const n = v.length;
    if (!n) return [];
    const mean = v.reduce((a, b) => a + b, 0) / n;
    const sd = Math.sqrt(v.reduce((a, b) => a + (b - mean) ** 2, 0) / Math.max(1, n - 1)) || 1e-9;
    const sorted = v.slice().sort((a, b) => a - b);
    const q = p => sorted[Math.min(n - 1, Math.max(0, Math.round((n - 1) * p)))];
    const iqr = q(0.75) - q(0.25);
    // Silverman: 0.9 * min(sd, IQR/1.34) * n^(-1/5)
    let bw = 0.9 * Math.min(sd, (iqr || sd) / 1.34) * Math.pow(n, -0.2);
    if (!(bw > 0)) bw = (hi - lo) / 50 || 1e-6;
    const out = [];
    let peak = 0;
    for (let i = 0; i < gridN; i++) {
      const x = lo + (hi - lo) * i / (gridN - 1);
      let d = 0;
      for (let j = 0; j < n; j++) {
        const z = (x - v[j]) / bw;
        d += Math.exp(-0.5 * z * z);
      }
      d /= (n * bw * Math.sqrt(2 * Math.PI));
      if (d > peak) peak = d;
      out.push([x, d]);
    }
    if (peak > 0) out.forEach(p => { p[1] /= peak; });
    return out;
  }

  // Round an interval outward to human-readable bounds (1/2/5 x 10^k). Explicit
  // axis min/max are unavoidable for the ridgeline (a `custom` series carries
  // index data, so ECharts cannot infer the value extent from it) and raw
  // float extents would print as "44.228926583" on the axis.
  function _niceBounds(lo, hi) {
    const span = hi - lo;
    if (!(span > 0)) return [lo - 1, hi + 1];
    const raw = span / 4;
    const mag = Math.pow(10, Math.floor(Math.log10(raw)));
    const norm = raw / mag;
    const step = mag * (norm >= 5 ? 5 : norm >= 2 ? 2 : 1);
    return [Math.floor(lo / step) * step, Math.ceil(hi / step) * step];
  }

  function _median(a) {
    const s = a.slice().sort((x, y) => x - y);
    if (!s.length) return 0;
    const m = Math.floor(s.length / 2);
    return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
  }

  /**
   * Distribution of a continuous variable, one row per unit (sample/group).
   * Picks the form from the data, per REPORT_VIZ_GUIDE §4:
   *   median group n < VIZ.densityMinN -> strip plot (real points, no KDE)
   *   otherwise                        -> kernel density
   *   groups > VIZ.manyGroups          -> ridgeline layout
   * `cutoffs` draws threshold lines (MIMAG 50/90, CheckV tiers …) — the whole
   * point of showing the distribution is seeing how close things sit to them.
   */
  window.distPlot = function (opts) {
    const groups  = (opts.groups || []).filter(g => (g.values || []).length);
    const log     = !!opts.log;
    const cutoffs = opts.cutoffs || [];
    const dark    = document.documentElement.dataset.theme === 'dark';
    const axisInk = dark ? '#94a3b8' : '#64748b';
    if (!groups.length) {
      return { title: { text: opts.title || '' },
               graphic: { type: 'text', left: 'center', top: 'middle',
                          style: { text: 'No data', fill: axisInk, fontSize: 12 } } };
    }

    const tx = v => log ? Math.log10(Math.max(v, 1e-9)) : v;
    const inv = v => log ? Math.pow(10, v) : v;
    const all = [];
    groups.forEach(g => g.values.forEach(v => { if (Number.isFinite(v)) all.push(tx(v)); }));
    let lo = Math.min(...all), hi = Math.max(...all);
    if (!(hi > lo)) { hi = lo + 1; lo = lo - 1; }
    const padv = (hi - lo) * 0.06;
    // On a log axis the working values are exponents, so whole numbers are the
    // readable bounds (10^4, 10^5); _niceBounds would land on 10^6.5 = 3162278.
    [lo, hi] = log
      ? [Math.floor(lo - padv), Math.ceil(hi + padv)]
      : _niceBounds(lo - padv, hi + padv);
    // Bounded measures (percentages) should not get an axis running to -20/120
    // just because the rounding padded outward — callers pass their real domain.
    if (opts.xMin !== undefined) lo = Math.max(lo, tx(opts.xMin));
    if (opts.xMax !== undefined) hi = Math.min(hi, tx(opts.xMax));

    // Row colour. Ordinal groups (pipeline stages, tiers) pass their own ramp.
    // Past the categorical budget every row takes ONE hue instead of cycling
    // PAL: identity here is carried by the row label and position, so a cycled
    // hue would claim a distinction it cannot deliver (slot 1 and slot 9 are the
    // same colour) -- see REPORT_VIZ_GUIDE §5.
    const _single = groups.length > window.VIZ.maxSeries;
    const _pal = i => (opts.colors && opts.colors.length)
      ? opts.colors[i % opts.colors.length]
      : (_single ? window.PAL[0] : window.PAL[i]);
    const useDots = _median(groups.map(g => g.values.length)) < window.VIZ.densityMinN;
    const rows    = groups.length;
    const rowH    = rows > window.VIZ.manyGroups ? 26 : 46;
    const height  = Math.max(220, 74 + rows * rowH);
    const names   = groups.map(g => g.name);
    const overlap = rows > window.VIZ.manyGroups ? 1.5 : 0.82;
    // Row labels come from a value axis, so `interval: 1` only lands on integers
    // when min is itself an integer — otherwise ticks fall on 0.4/1.4/… and the
    // formatter looks up names[0.4] and renders nothing.
    const minY = Math.floor(useDots ? -0.8 : -(overlap + 0.35));

    const cutLines = cutoffs.map(c => ({
      xAxis: inv(tx(c.value)),
      lineStyle: { color: c.color || (dark ? '#f87171' : '#ef4444'), type: 'dashed', width: 1.2 },
      label: { formatter: c.label || String(c.value), color: axisInk, fontSize: 10, position: 'insideEndTop' },
    }));

    const base = {
      __height: height,
      title: { text: opts.title || '' },
      legend: { show: false },
      grid: { top: 44, bottom: 46, left: 12, right: 22, containLabel: true },
      xAxis: {
        type: log ? 'log' : 'value',
        name: opts.xName || '', nameLocation: 'middle', nameGap: 28,
        min: inv(lo), max: inv(hi),
        splitLine: { show: true },
      },
      yAxis: {
        type: 'value', min: minY, max: rows - 1 + 0.45, inverse: true,
        interval: 1,
        axisLabel: { formatter: v => names[v] !== undefined ? names[v] : '', color: axisInk, fontSize: 11 },
        splitLine: { show: false }, axisLine: { show: false }, axisTick: { show: false },
      },
    };

    if (useDots) {
      // Strip plot: every observation drawn, deterministic jitter within the row.
      const pts = [];
      groups.forEach((g, i) => g.values.forEach((v, k) => {
        if (!Number.isFinite(v)) return;
        const jitter = ((k * 2654435761 % 1000) / 1000 - 0.5) * 0.36;
        pts.push([inv(tx(v)), i + jitter, g.name, i]);
      }));
      base.tooltip = {
        trigger: 'item',
        formatter: p => `${p.data[2]}<br>${opts.xName || 'value'}: ${(+p.data[0]).toLocaleString()}`,
      };
      base.series = [{
        type: 'scatter', symbolSize: 7, data: pts,
        itemStyle: { color: p => _pal(p.data[3]), opacity: 0.72,
                     borderColor: dark ? '#0b1220' : '#fff', borderWidth: 1 },
        markLine: cutLines.length ? { silent: true, symbol: 'none', data: cutLines } : undefined,
      }];
      base.__mode = 'dots';
      return base;
    }

    // Density / ridgeline: one filled polygon per row, drawn by a single custom
    // series (one series per group would blow past the legend/series budget).
    const curves = groups.map(g => _kde(g.values.map(tx), 72, lo, hi));
    base.tooltip = {
      trigger: 'item',
      formatter: p => {
        const g = groups[p.dataIndex];
        const s = g.values.slice().sort((a, b) => a - b);
        return `<strong>${g.name}</strong><br>n = ${s.length}`
             + `<br>median: ${(+_median(s)).toLocaleString()}`
             + `<br>range: ${(+s[0]).toLocaleString()} – ${(+s[s.length - 1]).toLocaleString()}`;
      },
    };
    base.series = [{
      type: 'custom',
      data: groups.map((_, i) => i),
      renderItem: (params, api) => {
        const i = params.dataIndex;
        const curve = curves[i];
        if (!curve.length) return null;
        const top = curve.map(p => api.coord([inv(p[0]), i - p[1] * overlap]));
        const bot = curve.map(p => api.coord([inv(p[0]), i])).reverse();
        const col = _pal(i);
        return {
          type: 'polygon',
          shape: { points: top.concat(bot) },
          style: { fill: col, opacity: 0.72, stroke: col, lineWidth: 1.4 },
        };
      },
      markLine: cutLines.length ? { silent: true, symbol: 'none', data: cutLines } : undefined,
    }];
    base.__mode = 'density';
    return base;
  };

  /**
   * One value per sample, optionally several series. Vertical bars while the
   * sample count is small; past VIZ.manySamples the axis flips horizontal
   * (names read straight, height grows) and rows sort by total, which keeps
   * every sample visible instead of degrading into an unreadable comb.
   */
  window.samplesBar = function (opts) {
    const samples = opts.samples || [];
    const series  = opts.series  || [];
    const many    = samples.length > window.VIZ.manySamples;
    const stack   = opts.stack ? 'total' : undefined;

    let order = samples.map((_, i) => i);
    if (many && opts.sort !== false) {
      const totals = samples.map((_, i) => series.reduce((a, s) => a + (+s.data[i] || 0), 0));
      order = order.sort((a, b) => totals[a] - totals[b]);   // ascending: biggest on top
    }
    const cats = order.map(i => samples[i]);
    const mk   = s => ({
      name: s.name, type: 'bar', stack,
      itemStyle: { color: s.color, borderRadius: stack ? 0 : [0, 3, 3, 0] },
      barMaxWidth: many ? 18 : 42,
      data: order.map(i => s.data[i]),
    });

    // Reference lines (QC floors like "80% retained") sit on the VALUE axis,
    // which swaps sides when the chart flips horizontal — callers give a value,
    // not an axis, so the line follows the orientation automatically.
    const refs = (opts.refLines || []).map(r => ({
      [many ? 'xAxis' : 'yAxis']: r.value,
      lineStyle: { type: 'dashed', color: r.color || '#d97706' },
      label: { formatter: r.label || String(r.value), fontSize: 10 },
    }));
    const withRefs = arr => (refs.length
      ? arr.map((s, i) => i === 0 ? { ...s, markLine: { silent: true, symbol: 'none', data: refs } } : s)
      : arr);

    const common = {
      title: { text: opts.title || '' },
      legend: series.length > 1 ? { data: series.map(s => s.name) } : { show: false },
      tooltip: { trigger: 'axis', axisPointer: { type: 'shadow' } },
      series: withRefs(series.map(mk)),
    };

    if (!many) {
      return {
        ...common,
        xAxis: { type: 'category', data: cats, axisLabel: { rotate: 30 } },
        yAxis: { type: 'value', name: opts.valueName || '', min: opts.valueMin, max: opts.valueMax },
        grid: { bottom: 70 },
        series: common.series.map(s => ({ ...s, itemStyle: { ...s.itemStyle, borderRadius: stack ? 0 : [3, 3, 0, 0] } })),
      };
    }
    return {
      ...common,
      __height: Math.max(260, 70 + cats.length * 22),
      xAxis: { type: 'value', name: opts.valueName || '', nameLocation: 'middle', nameGap: 28,
               min: opts.valueMin, max: opts.valueMax },
      yAxis: { type: 'category', data: cats, axisLabel: { fontSize: 11 } },
      grid: { top: series.length > 1 ? 58 : 36, bottom: 46, left: 12, right: 26, containLabel: true },
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
    // Forms whose height depends on the data (horizontal per-sample bars,
    // ridgelines) ask for it via __height — the CSS box height is a default for
    // fixed-height charts, not a cap on a chart with 200 rows.
    if (option.__height) {
      const h = Math.round(option.__height);
      if (el.clientHeight !== h) {
        el.style.height = h + 'px';
        if (window._charts[id]) window._charts[id].resize();
      }
    }
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

  // Hide nav tabs whose track did not run (TRACKS injected from the report data).
  (function hideDisabledTrackTabs() {
    if (typeof TRACKS === 'undefined' || !TRACKS) return;
    const TAB_TRACK = {
      viral: 'viral',
      prokaryotic: 'prok',
      hostdefense: 'integration',
      'reads-classify': 'reads',
      coassembly: 'coassembly',
    };
    document.querySelectorAll('.nav-tab').forEach(btn => {
      const key = TAB_TRACK[btn.dataset.tab];
      if (key && TRACKS[key] === false) btn.style.display = 'none';
    });
  })();

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
     renderHostDefense, renderDiversity, renderAnnotation, renderReadsClassify, renderCoassembly, renderAbout].forEach(fn => {
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
