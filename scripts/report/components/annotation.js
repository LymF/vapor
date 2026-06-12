/* annotation.js — Functional annotation (COG/PHROGS), phage tables, genome maps */
(function () {
  'use strict';

  window.renderAnnotation = function () {
    const samples = typeof SAMPLES !== 'undefined' ? SAMPLES : [];
    _renderFunctional(samples);
    makeSampleDropdown('sample-sel-ann',  _renderPhageAnnotation);
    makeSampleDropdown('sample-sel-maps', _renderMaps);
    _updateMapsSampleLabels();
    // Map mode selector
    const modeEl = document.getElementById('genome-map-mode');
    if (modeEl) modeEl.addEventListener('change', () => {
      _updateMapsSampleLabels();
      const sel = document.getElementById('sample-sel-maps');
      if (sel) _renderMaps(sel.value);
    });
  };

  // ── Annotate the sample dropdown with per-mode genome counts ─────────────
  function _updateMapsSampleLabels() {
    const sel    = document.getElementById('sample-sel-maps');
    const modeEl = document.getElementById('genome-map-mode');
    if (!sel || !modeEl) return;
    const mode = modeEl.value;
    const maps = typeof GENOME_MAPS !== 'undefined' ? GENOME_MAPS : {};
    [...sel.options].forEach(opt => {
      const count = ((maps[opt.value] || {})[mode] || []).length;
      opt.textContent = `${opt.value} (${count})`;
    });
  }

  // ── COG + PHROGS stacked bars ─────────────────────────────────────────────
  function _renderFunctional(samples) {
    const eggnog = typeof EGGNOG_DATA !== 'undefined' ? EGGNOG_DATA : {};
    const phrogs = typeof PHROGS_DATA !== 'undefined' ? PHROGS_DATA : {};

    // COG categories — collect all categories across samples
    const cogCatSet = new Set();
    samples.forEach(s => Object.keys(eggnog[s] || {}).forEach(c => cogCatSet.add(c)));
    const cogCats = [...cogCatSet].sort();

    const cogSeries = cogCats.map((cat, i) => ({
      name:  cat,
      type:  'bar',
      stack: 'cog',
      color: PAL[i % PAL.length],
      data:  samples.map(s => (eggnog[s] || {})[cat] || 0),
    }));

    mkChart('ann-cog-chart', {
      title:   { text: 'COG Functional Categories (EggNOG-mapper)' },
      tooltip: { trigger: 'axis', axisPointer: { type: 'shadow' } },
      legend:  { type: 'scroll', data: cogCats, top: 'bottom' },
      xAxis:   { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis:   { type: 'value', name: 'Gene count' },
      series:  cogSeries,
      grid:    { bottom: 90 },
    });

    // PHROGS categories
    const phrogsCatSet = new Set();
    samples.forEach(s => Object.keys(phrogs[s] || {}).forEach(c => phrogsCatSet.add(c)));
    const phrogsCats = [...phrogsCatSet].sort();

    const phrogsSeries = phrogsCats.map((cat, i) => ({
      name:  cat,
      type:  'bar',
      stack: 'ph',
      color: PAL[i % PAL.length],
      data:  samples.map(s => (phrogs[s] || {})[cat] || 0),
    }));

    mkChart('ann-phrogs-chart', {
      title:   { text: 'PHROGS Functional Categories (Pharokka)' },
      tooltip: { trigger: 'axis', axisPointer: { type: 'shadow' } },
      legend:  { type: 'scroll', data: phrogsCats, top: 'bottom' },
      xAxis:   { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis:   { type: 'value', name: 'Gene count' },
      series:  phrogsSeries,
      grid:    { bottom: 90 },
    });
  }

  // ── Phage annotation table (VIBRANT scaffolds) ────────────────────────────
  function _renderPhageAnnotation(sample) {
    const vibrant = (typeof VIBRANT_DATA !== 'undefined' ? VIBRANT_DATA : [])
      .filter(r => r.sample === sample);
    makeTable('vibrant-table', vibrant, [
      { key: 'Scaffold',    label: 'Scaffold' },
      { key: 'Total_genes', label: 'Total genes' },
      { key: 'VOG_score',   label: 'VOG score' },
      { key: 'Pfam_score',  label: 'Pfam score' },
      { key: 'KEGG_score',  label: 'KEGG score' },
    ]);
  }

  // ── Genome maps ───────────────────────────────────────────────────────────
  function _renderMaps(sample) {
    const modeEl = document.getElementById('genome-map-mode');
    const mode   = modeEl ? modeEl.value : 'phage';
    const maps   = (typeof GENOME_MAPS !== 'undefined' ? GENOME_MAPS : {})[sample] || {};
    const items  = (maps[mode] || []);
    const cont   = document.getElementById('genome-maps-container');
    if (!cont) return;

    if (!items.length) {
      cont.innerHTML = '<p style="color:var(--text-muted);font-size:.85rem;padding:.5rem">No genome maps for this selection.</p>';
      return;
    }

    cont.innerHTML = items.map(m =>
      `<div class="genome-map-item">
         <h4>${m.id}</h4>
         ${m.svg}
       </div>`
    ).join('');

    // Fix SVG sizing for dark/light themes
    cont.querySelectorAll('svg').forEach(svg => {
      svg.style.maxWidth = '100%';
      svg.style.height   = 'auto';
    });

    // Export toolbars (PNG/SVG/PDF) for each genome map
    if (window.VaporExport) {
      cont.querySelectorAll('.genome-map-item').forEach(item => window.VaporExport.attachToSVGHost(item));
    }
  }

})();
