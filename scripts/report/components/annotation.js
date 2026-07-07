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

    mkChart('ann-cog-chart', {
      title:   { text: 'COG Functional Categories (EggNOG-mapper)' },
      tooltip: { trigger: 'axis', axisPointer: { type: 'shadow' } },
      legend:  { type: 'scroll', data: cogTop, top: 'bottom' },
      xAxis:   { type: 'category', data: samples, axisLabel: { rotate: 30 } },
      yAxis:   { type: 'value', name: 'Gene count' },
      series:  cogSeries,
      grid:    { bottom: 90 },
    });

    // PHROGS categories — same fold-to-top-7+Other treatment as COG above.
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

    mkChart('ann-phrogs-chart', {
      title:   { text: 'PHROGS Functional Categories (Pharokka)' },
      tooltip: { trigger: 'axis', axisPointer: { type: 'shadow' } },
      legend:  { type: 'scroll', data: phrogsTop, top: 'bottom' },
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
    const mode   = modeEl ? modeEl.value : 'virus';
    const maps   = (typeof GENOME_MAPS !== 'undefined' ? GENOME_MAPS : {})[sample] || {};
    const items  = (maps[mode] || []);
    const cont   = document.getElementById('genome-maps-container');
    if (!cont) return;

    cont.innerHTML = '';

    if (!items.length) {
      cont.innerHTML = '<p style="color:var(--text-muted);font-size:.85rem;padding:.5rem">No genome maps for this selection.</p>';
      return;
    }

    items.forEach(m => {
      const item = document.createElement('div');
      item.className = 'genome-map-item';

      const h4 = document.createElement('h4');
      h4.textContent = m.id;
      // Phage/Virus badge: the "Virus" view merges both categories into one
      // list -- the badge shows which one PHROGS hallmark-gene evidence
      // assigned, instead of forcing a separate mode toggle per category.
      if (m.category) {
        const badge = document.createElement('span');
        badge.className = 'badge ' + (m.category === 'Phage' ? 'badge-teal' : 'badge-amber');
        badge.style.marginLeft = '.5rem';
        badge.textContent = m.category;
        h4.appendChild(badge);
      }
      item.appendChild(h4);

      // Copy FASTA button
      if (m.seq) {
        const copyBtn = document.createElement('button');
        copyBtn.type = 'button';
        copyBtn.className = 'btn-sm';
        copyBtn.style.cssText = 'margin-bottom:.5rem;display:inline-flex;align-items:center;gap:.3rem';
        copyBtn.textContent = 'Copy FASTA';
        copyBtn.title = 'Copy genome sequence (FASTA) to clipboard';
        copyBtn.addEventListener('click', () => {
          const wrapped = (m.seq.match(/.{1,60}/g) || []).join('\n');
          const fasta = `>${m.id}\n${wrapped}`;
          (navigator.clipboard
            ? navigator.clipboard.writeText(fasta)
            : Promise.reject()
          ).catch(() => {
            const ta = document.createElement('textarea');
            ta.value = fasta;
            document.body.appendChild(ta);
            ta.select();
            document.execCommand('copy');
            document.body.removeChild(ta);
          }).finally(() => {
            copyBtn.textContent = 'Copied!';
            setTimeout(() => { copyBtn.textContent = 'Copy FASTA'; }, 1500);
          });
        });
        item.appendChild(copyBtn);
      }

      // SVG content
      const wrap = document.createElement('div');
      wrap.innerHTML = m.svg;
      const svg = wrap.querySelector('svg');
      if (svg) { svg.style.maxWidth = '100%'; svg.style.height = 'auto'; }
      item.appendChild(svg || wrap);

      cont.appendChild(item);
      if (window.VaporExport) window.VaporExport.attachToSVGHost(item);
    });
  }

})();
