/* coassembly.js — Co-assembly MAGs tab */
(function () {
  'use strict';
  const CA = typeof COASSEMBLY_DATA !== 'undefined' ? COASSEMBLY_DATA : null;

  window.renderCoassembly = function () {
    const empty = document.getElementById('coassembly-empty');
    const tbl = document.getElementById('coassembly-table');
    if (!CA || !CA.has_data) { if (empty) empty.style.display = ''; if (tbl) tbl.innerHTML = ''; return; }
    if (empty) empty.style.display = 'none';
    const rows = [];
    for (const g of CA.groups) for (const m of g.mags) {
      rows.push(`<tr><td>${g.group}</td><td>${m.bin}</td>` +
        `<td>${(+m.completeness).toFixed(1)}%</td><td>${(+m.contamination).toFixed(1)}%</td>` +
        `<td>${m.classification || '—'}</td></tr>`);
    }
    tbl.innerHTML = `<table class="vapor-table"><thead><tr>` +
      `<th>Group</th><th>MAG</th><th>Completeness</th><th>Contamination</th><th>GTDB</th>` +
      `</tr></thead><tbody>${rows.join('')}</tbody></table>`;
  };
})();
