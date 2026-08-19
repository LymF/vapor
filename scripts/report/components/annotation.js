/* annotation.js — Functional annotation (COG/PHROGS) e tabelas de fago */
(function () {
  'use strict';

  window.renderAnnotation = function () {
    const samples = typeof SAMPLES !== 'undefined' ? SAMPLES : [];
    _renderFunctional(samples);
    _renderCazy(samples);
    _renderKeggModules();
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


  // ── CAZy ──────────────────────────────────────────────────────────────────
  // A coluna CAZy sempre esteve no eggnog_annotations.tsv (coluna 19) e
  // nenhum consumidor a lia ate 2026-08-19. Nao ha ferramenta nova aqui: e
  // o mesmo arquivo do painel de COG acima, uma coluna adiante.
  function _renderCazy(samples) {
    const cazy = typeof CAZY_DATA !== 'undefined' ? CAZY_DATA : {};
    if (!Object.keys(cazy).length) return;

    // As seis classes cabem no orcamento de 8 matizes sem dobrar em "Other",
    // ao contrario de COG (20 letras) e PHROGS.
    const classTotals = {};
    samples.forEach(s => Object.entries((cazy[s] || {}).by_class || {})
      .forEach(([c, v]) => { classTotals[c] = (classTotals[c] || 0) + v; }));
    const classes = Object.keys(classTotals).sort((a, b) => classTotals[b] - classTotals[a]);

    mkChart('ann-cazy-chart', samplesBar({
      samples, title: 'CAZyme Classes (EggNOG CAZy column)',
      valueName: 'Gene count', stack: true,
      series: classes.map((cls, i) => ({
        name: cls,
        color: cls === 'Other' ? PAL_MUTED : PAL[i % PAL.length],
        data: samples.map(s => ((cazy[s] || {}).by_class || {})[cls] || 0),
      })),
    }));

    // A familia (GH13, GT51) e o que carrega o significado enzimatico; a
    // classe so agrupa. Somada entre amostras, em barra horizontal.
    const famTotals = {};
    samples.forEach(s => ((cazy[s] || {}).top_families || [])
      .forEach(([fam, n]) => { famTotals[fam] = (famTotals[fam] || 0) + n; }));
    const fams = Object.entries(famTotals).sort((a, b) => b[1] - a[1]).slice(0, 15).reverse();

    mkChart('ann-cazy-family-chart', Object.assign(echartsTheme(), {
      title: { text: 'Top CAZyme Families (all samples)', left: 'center' },
      tooltip: { trigger: 'axis', axisPointer: { type: 'shadow' } },
      grid: { left: 90, right: 24, top: 48, bottom: 32 },
      xAxis: { type: 'value', name: 'Gene count' },
      yAxis: { type: 'category', data: fams.map(d => d[0]) },
      series: [{ type: 'bar', data: fams.map(d => d[1]), itemStyle: { color: PAL[0] } }],
    }));
  }

  // ── Completude de modulo KEGG ─────────────────────────────────────────────
  function _renderKeggModules() {
    const mods = typeof KEGG_MODULES !== 'undefined' ? KEGG_MODULES : {};
    if (!Object.keys(mods).length) return;

    makeSampleDropdown('sample-sel-kegg-mod', function (sample) {
      const rows = (mods[sample] || []).map(r => Object.assign({}, r, {
        complete: `${r.n_complete} / ${r.n_mags}`,
      }));
      makeTable('kegg-modules-table', rows, [
        { key: 'module',     label: 'Module' },
        { key: 'name',       label: 'Pathway' },
        { key: 'complete',   label: 'Complete (\u2265 80%)' },
        { key: 'mean',       label: 'Mean %' },
        { key: 'max',        label: 'Max %' },
        { key: 'missing_ko', label: 'Missing KO' },
      ], { pageSize: 30 });
    });
  }


})();
