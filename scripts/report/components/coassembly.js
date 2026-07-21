/* coassembly.js — Co-assembly MAGs tab */
(function () {
  'use strict';
  const CA = typeof COASSEMBLY_DATA !== 'undefined' ? COASSEMBLY_DATA : null;

  function esc(s) {
    return String(s == null ? '' : s)
      .replace(/&/g, '&amp;')
      .replace(/</g, '&lt;')
      .replace(/>/g, '&gt;');
  }

  window.renderCoassembly = function () {
    const empty = document.getElementById('coassembly-empty');
    const tbl = document.getElementById('coassembly-table');
    if (!CA || !CA.has_data) { if (empty) empty.style.display = ''; if (tbl) tbl.innerHTML = ''; return; }
    if (empty) empty.style.display = 'none';
    const rows = [];
    for (const g of CA.groups) for (const m of g.mags) {
      rows.push(`<tr><td>${esc(g.group)}</td><td>${esc(m.bin)}</td>` +
        `<td>${(+m.completeness).toFixed(1)}%</td><td>${(+m.contamination).toFixed(1)}%</td>` +
        `<td>${esc(m.classification) || '—'}</td></tr>`);
    }
    tbl.innerHTML = `<table class="vapor-table"><thead><tr>` +
      `<th>Group</th><th>MAG</th><th>Completeness</th><th>Contamination</th><th>GTDB</th>` +
      `</tr></thead><tbody>${rows.join('')}</tbody></table>`;

    const vEmpty = document.getElementById('coassembly-votus-empty');
    const vTbl = document.getElementById('coassembly-votus-table');
    const groupsWithVotus = CA.groups.filter(g => g.n_votus);
    if (!groupsWithVotus.length) {
      if (vEmpty) vEmpty.style.display = '';
      if (vTbl) vTbl.innerHTML = '';
      return;
    }
    if (vEmpty) vEmpty.style.display = 'none';
    const vRows = groupsWithVotus.map(g => {
      const fams = (g.votu_families || [])
        .map(f => `${esc(f.family)} (${f.count})`)
        .join(', ') || '—';
      return `<tr><td>${esc(g.group)}</td><td>${g.n_votus}</td><td>${fams}</td></tr>`;
    });
    vTbl.innerHTML = `<table class="vapor-table"><thead><tr>` +
      `<th>Group</th><th>vOTUs</th><th>Top families</th>` +
      `</tr></thead><tbody>${vRows.join('')}</tbody></table>`;
  };
})();
