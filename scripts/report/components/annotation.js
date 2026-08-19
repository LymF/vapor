/* annotation.js — Functional annotation (COG/PHROGS) e tabelas de fago */
(function () {
  'use strict';

  window.renderAnnotation = function () {
    const samples = typeof SAMPLES !== 'undefined' ? SAMPLES : [];
    _renderFunctional(samples);
  };

  // ── COG + PHROGS stacked bars ─────────────────────────────────────────────
  function _renderFunctional(samples) {
    const eggnog = typeof EGGNOG_DATA !== 'undefined' ? EGGNOG_DATA : {};
    const phrogs = typeof PHROGS_DATA !== 'undefined' ? PHROGS_DATA : {};

    // COG categories — EggNOG can report ~20 functional letters, well past
    // the validated 8-hue budget. Fold to the top 7 by total count across
    // samples + one "Other" bucket instead of cycling PAL past slot 8
    // (a validated 15-slot extension was tried and failed the CVD/chroma
    // checks — see window.PAL in app.js).
    const cogTotals = {};
    samples.forEach(s => Object.entries(eggnog[s] || {}).forEach(([c, v]) => { cogTotals[c] = (cogTotals[c] || 0) + v; }));
    const cogTop = foldOther(cogTotals, 7).map(d => d[0]);

    const cogSeries = cogTop.map((cat, i) => ({
      name:  cat,
      type:  'bar',
      stack: 'cog',
      color: cat === 'Other' ? PAL_MUTED : PAL[i % PAL.length],
      data:  samples.map(s => cat === 'Other'
        ? Object.entries(eggnog[s] || {}).reduce((a, [c, v]) => a + (cogTop.includes(c) ? 0 : v), 0)
        : (eggnog[s] || {})[cat] || 0),
    }));

    mkChart('ann-cog-chart', samplesBar({
      samples, title: 'COG Functional Categories (EggNOG-mapper)',
      valueName: 'Gene count', stack: true,
      series: cogSeries.map(s => ({ name: s.name, data: s.data,
                                    color: s.color || (s.itemStyle || {}).color })),
    }));

    // PHROGS categories — same fold-to-top-7+Other treatment as COG above.
    // Pharokka moved to the global vOTU catalog on 2026-08-18 ("(h)"), so
    // phrogs[s] is now the SAME catalog-wide counts object for every
    // sample s (see load_phrogs, scripts/report/data_loaders.py) -- the
    // per-sample bars below will all read identical, which is expected,
    // not a bug. Title says so explicitly instead of implying a
    // per-sample recount.
    const phrogsTotals = {};
    samples.forEach(s => Object.entries(phrogs[s] || {}).forEach(([c, v]) => { phrogsTotals[c] = (phrogsTotals[c] || 0) + v; }));
    const phrogsTop = foldOther(phrogsTotals, 7).map(d => d[0]);

    const phrogsSeries = phrogsTop.map((cat, i) => ({
      name:  cat,
      type:  'bar',
      stack: 'ph',
      color: cat === 'Other' ? PAL_MUTED : PAL[i % PAL.length],
      data:  samples.map(s => cat === 'Other'
        ? Object.entries(phrogs[s] || {}).reduce((a, [c, v]) => a + (phrogsTop.includes(c) ? 0 : v), 0)
        : (phrogs[s] || {})[cat] || 0),
    }));

    mkChart('ann-phrogs-chart', samplesBar({
      samples, title: 'PHROGS Functional Categories (Pharokka, vOTU catalog — global, same across samples)',
      valueName: 'Gene count', stack: true,
      series: phrogsSeries.map(s => ({ name: s.name, data: s.data,
                                       color: s.color || (s.itemStyle || {}).color })),
    }));
  }


})();
