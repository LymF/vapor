/* export.js — PNG / SVG / PDF export for charts, tables and genome maps */
(function () {
  'use strict';

  // ── Generic download helpers ───────────────────────────────────────────────
  function downloadDataURL(dataURL, filename) {
    const a = document.createElement('a');
    a.href = dataURL;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
  }

  function downloadBlob(blob, filename) {
    const url = URL.createObjectURL(blob);
    downloadDataURL(url, filename);
    setTimeout(() => URL.revokeObjectURL(url), 2000);
  }

  function escapeXML(s) {
    return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
  }

  function slugify(s) {
    return String(s).toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/(^-|-$)/g, '') || 'export';
  }

  // ── SVG → string / canvas ────────────────────────────────────────────────
  function svgToString(svgEl) {
    const clone = svgEl.cloneNode(true);
    clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    const rect = svgEl.getBoundingClientRect();
    if (!clone.getAttribute('width'))  clone.setAttribute('width',  rect.width  || 800);
    if (!clone.getAttribute('height')) clone.setAttribute('height', rect.height || 500);
    return new XMLSerializer().serializeToString(clone);
  }

  function svgToCanvas(svgEl, scale, bg) {
    return new Promise((resolve, reject) => {
      const rect = svgEl.getBoundingClientRect();
      const w = parseFloat(svgEl.getAttribute('width'))  || rect.width  || 800;
      const h = parseFloat(svgEl.getAttribute('height')) || rect.height || 500;
      const svgStr = svgToString(svgEl);
      const blob = new Blob([svgStr], { type: 'image/svg+xml;charset=utf-8' });
      const url  = URL.createObjectURL(blob);
      const img  = new Image();
      img.onload = () => {
        const canvas = document.createElement('canvas');
        canvas.width  = Math.max(1, Math.round(w * scale));
        canvas.height = Math.max(1, Math.round(h * scale));
        const ctx = canvas.getContext('2d');
        ctx.fillStyle = bg;
        ctx.fillRect(0, 0, canvas.width, canvas.height);
        ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
        URL.revokeObjectURL(url);
        resolve(canvas);
      };
      img.onerror = (e) => { URL.revokeObjectURL(url); reject(e); };
      img.src = url;
    });
  }

  // ── ECharts → SVG string (re-renders the current option with an SVG renderer) ──
  function echartsToSVGString(chart) {
    const opt = chart.getOption();
    const host = document.createElement('div');
    host.style.position = 'absolute';
    host.style.left   = '-99999px';
    host.style.top    = '0';
    host.style.width  = chart.getWidth()  + 'px';
    host.style.height = chart.getHeight() + 'px';
    document.body.appendChild(host);
    const tmp = echarts.init(host, null, { renderer: 'svg' });
    tmp.setOption(opt);
    const svgEl = host.querySelector('svg');
    const svgStr = svgEl ? svgToString(svgEl) : '';
    tmp.dispose();
    document.body.removeChild(host);
    return svgStr;
  }

  // ── Minimal single-page PDF wrapping one JPEG image (no external deps) ──────
  function pdfFromJPEGDataURL(dataURL, widthPx, heightPx) {
    const base64 = dataURL.split(',')[1];
    const bin = atob(base64);
    const jpegBytes = new Uint8Array(bin.length);
    for (let i = 0; i < bin.length; i++) jpegBytes[i] = bin.charCodeAt(i);

    // px → pt at 96dpi
    const ptW = (widthPx  * 72 / 96).toFixed(2);
    const ptH = (heightPx * 72 / 96).toFixed(2);

    const enc = s => new TextEncoder().encode(s);
    const parts = [];
    const objOffsets = [0];
    let pos = 0;
    function add(bytes) { parts.push(bytes); pos += bytes.length; }
    function beginObj(n) { objOffsets[n] = pos; }

    add(enc('%PDF-1.4\n'));

    beginObj(1);
    add(enc('1 0 obj\n<< /Type /Catalog /Pages 2 0 R >>\nendobj\n'));

    beginObj(2);
    add(enc('2 0 obj\n<< /Type /Pages /Kids [3 0 R] /Count 1 >>\nendobj\n'));

    beginObj(3);
    add(enc(`3 0 obj\n<< /Type /Page /Parent 2 0 R /MediaBox [0 0 ${ptW} ${ptH}] /Resources << /XObject << /Im0 4 0 R >> >> /Contents 5 0 R >>\nendobj\n`));

    beginObj(4);
    add(enc(`4 0 obj\n<< /Type /XObject /Subtype /Image /Width ${widthPx} /Height ${heightPx} /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /DCTDecode /Length ${jpegBytes.length} >>\nstream\n`));
    add(jpegBytes);
    add(enc('\nendstream\nendobj\n'));

    beginObj(5);
    const content = `q ${ptW} 0 0 ${ptH} 0 0 cm /Im0 Do Q`;
    const contentBytes = enc(content);
    add(enc(`5 0 obj\n<< /Length ${contentBytes.length} >>\nstream\n`));
    add(contentBytes);
    add(enc('\nendstream\nendobj\n'));

    const xrefStart = pos;
    let xref = 'xref\n0 6\n0000000000 65535 f \n';
    for (let i = 1; i <= 5; i++) xref += String(objOffsets[i]).padStart(10, '0') + ' 00000 n \n';
    add(enc(xref));
    add(enc(`trailer\n<< /Size 6 /Root 1 0 R >>\nstartxref\n${xrefStart}\n%%EOF`));

    return new Blob(parts, { type: 'application/pdf' });
  }

  // ── Table → TSV ──────────────────────────────────────────────────────────
  function tableToTSV(tableEl) {
    return [...tableEl.querySelectorAll('tr')]
      .map(tr => [...tr.children].map(td => td.textContent.replace(/[\t\n\r]/g, ' ').trim()).join('\t'))
      .join('\n');
  }

  // ── Table → layout / canvas / svg ────────────────────────────────────────
  function tableLayout(tableEl) {
    const rows = [...tableEl.querySelectorAll('tr')];
    const colCount = rows.reduce((m, r) => Math.max(m, r.children.length), 0);
    const colWidths = new Array(colCount).fill(0);
    const rowHeights = [];
    const cells = [];
    rows.forEach(tr => {
      const rowCells = [];
      [...tr.children].forEach((td, ci) => {
        const rect = td.getBoundingClientRect();
        colWidths[ci] = Math.max(colWidths[ci], rect.width);
        rowCells.push({ text: td.textContent.trim(), isHeader: td.tagName === 'TH' });
      });
      rowHeights.push(tr.getBoundingClientRect().height || 28);
      cells.push(rowCells);
    });
    return { colWidths, rowHeights, cells };
  }

  function tableToCanvas(tableEl, scale) {
    const { colWidths, rowHeights, cells } = tableLayout(tableEl);
    const dark = document.documentElement.dataset.theme === 'dark';
    const totalW = Math.max(1, Math.ceil(colWidths.reduce((a, b) => a + b, 0)));
    const totalH = Math.max(1, Math.ceil(rowHeights.reduce((a, b) => a + b, 0)));
    const canvas = document.createElement('canvas');
    canvas.width  = Math.ceil(totalW * scale);
    canvas.height = Math.ceil(totalH * scale);
    const ctx = canvas.getContext('2d');
    ctx.scale(scale, scale);
    ctx.fillStyle = dark ? '#1e293b' : '#ffffff';
    ctx.fillRect(0, 0, totalW, totalH);

    let y = 0;
    cells.forEach((row, ri) => {
      let x = 0;
      const h = rowHeights[ri];
      row.forEach((cell, ci) => {
        const w = colWidths[ci] || 0;
        if (cell.isHeader) {
          ctx.fillStyle = dark ? '#0f172a' : '#f1f5f9';
          ctx.fillRect(x, y, w, h);
        }
        ctx.strokeStyle = dark ? '#334155' : '#e2e8f0';
        ctx.lineWidth = 1;
        ctx.strokeRect(x, y, w, h);
        ctx.save();
        ctx.beginPath();
        ctx.rect(x, y, w, h);
        ctx.clip();
        ctx.fillStyle = dark ? '#f1f5f9' : '#0f172a';
        ctx.font = cell.isHeader ? '600 11px Inter, system-ui, sans-serif' : '12px Inter, system-ui, sans-serif';
        ctx.textBaseline = 'middle';
        ctx.fillText(cell.text, x + 8, y + h / 2, Math.max(0, w - 16));
        ctx.restore();
        x += w;
      });
      y += h;
    });
    return canvas;
  }

  function tableToSVG(tableEl) {
    const { colWidths, rowHeights, cells } = tableLayout(tableEl);
    const dark = document.documentElement.dataset.theme === 'dark';
    const totalW = Math.max(1, Math.ceil(colWidths.reduce((a, b) => a + b, 0)));
    const totalH = Math.max(1, Math.ceil(rowHeights.reduce((a, b) => a + b, 0)));
    let svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${totalW}" height="${totalH}" `
            + `viewBox="0 0 ${totalW} ${totalH}" font-family="Inter, system-ui, sans-serif" font-size="12">`;
    svg += `<rect width="${totalW}" height="${totalH}" fill="${dark ? '#1e293b' : '#ffffff'}"/>`;
    let y = 0;
    cells.forEach((row, ri) => {
      let x = 0;
      const h = rowHeights[ri];
      row.forEach((cell, ci) => {
        const w = colWidths[ci] || 0;
        if (cell.isHeader) {
          svg += `<rect x="${x}" y="${y}" width="${w}" height="${h}" fill="${dark ? '#0f172a' : '#f1f5f9'}"/>`;
        }
        svg += `<rect x="${x}" y="${y}" width="${w}" height="${h}" fill="none" stroke="${dark ? '#334155' : '#e2e8f0'}"/>`;
        svg += `<text x="${x + 8}" y="${y + h / 2}" dominant-baseline="middle" fill="${dark ? '#f1f5f9' : '#0f172a'}"`
             + `${cell.isHeader ? ' font-weight="600" font-size="11"' : ''}>${escapeXML(cell.text)}</text>`;
        x += w;
      });
      y += h;
    });
    svg += '</svg>';
    return svg;
  }

  // ── Resolve the exportable thing inside a host element ───────────────────
  function resolveTarget(hostEl) {
    const echartsBox = hostEl.querySelector('[id].echarts-box, .echarts-box[id]');
    if (echartsBox && window._charts && window._charts[echartsBox.id] && !window._charts[echartsBox.id].isDisposed()) {
      return { type: 'echarts', id: echartsBox.id };
    }
    const svg = hostEl.querySelector('svg');
    if (svg) return { type: 'svg', el: svg };
    const table = hostEl.querySelector('table');
    if (table) return { type: 'table', el: table };
    return null;
  }

  function defaultFilename(hostEl, target) {
    const titleEl = hostEl.querySelector('h3.section-title, h4');
    let base;
    if (titleEl) base = titleEl.textContent.trim();
    else if (target.type === 'echarts') base = target.id;
    else if (target.el && target.el.id) base = target.el.id;
    else base = 'vapor-export';
    return slugify(base);
  }

  // ── Perform the actual export ────────────────────────────────────────────
  async function handleExport(target, fmt, name) {
    const dark = document.documentElement.dataset.theme === 'dark';
    const bg = dark ? '#1e293b' : '#ffffff';
    try {
      if (fmt === 'TSV') {
        if (target.type === 'table') {
          // Leading BOM so Excel detects UTF-8 instead of guessing the
          // system codepage -- without it, any non-ASCII char (e.g. the
          // '—' placeholder used for missing values) comes out as
          // mojibake ("â€”") when opened in Excel on Windows.
          downloadBlob(new Blob(['﻿' + tableToTSV(target.el)], { type: 'text/tab-separated-values' }), `${name}.tsv`);
        }
        return;
      }
      if (target.type === 'echarts') {
        const chart = window._charts[target.id];
        if (fmt === 'PNG') {
          downloadDataURL(chart.getDataURL({ type: 'png', pixelRatio: 2, backgroundColor: bg }), `${name}.png`);
        } else if (fmt === 'SVG') {
          downloadBlob(new Blob([echartsToSVGString(chart)], { type: 'image/svg+xml' }), `${name}.svg`);
        } else {
          const pixelRatio = 2;
          const url = chart.getDataURL({ type: 'jpeg', pixelRatio, backgroundColor: bg });
          downloadBlob(pdfFromJPEGDataURL(url, chart.getWidth() * pixelRatio, chart.getHeight() * pixelRatio), `${name}.pdf`);
        }
      } else if (target.type === 'svg') {
        if (fmt === 'SVG') {
          downloadBlob(new Blob([svgToString(target.el)], { type: 'image/svg+xml' }), `${name}.svg`);
        } else if (fmt === 'PNG') {
          const canvas = await svgToCanvas(target.el, 2, bg);
          downloadDataURL(canvas.toDataURL('image/png'), `${name}.png`);
        } else {
          const canvas = await svgToCanvas(target.el, 2, bg);
          downloadBlob(pdfFromJPEGDataURL(canvas.toDataURL('image/jpeg', 0.92), canvas.width, canvas.height), `${name}.pdf`);
        }
      } else if (target.type === 'table') {
        if (fmt === 'SVG') {
          downloadBlob(new Blob([tableToSVG(target.el)], { type: 'image/svg+xml' }), `${name}.svg`);
        } else if (fmt === 'PNG') {
          downloadDataURL(tableToCanvas(target.el, 2).toDataURL('image/png'), `${name}.png`);
        } else {
          const canvas = tableToCanvas(target.el, 2);
          downloadBlob(pdfFromJPEGDataURL(canvas.toDataURL('image/jpeg', 0.92), canvas.width, canvas.height), `${name}.pdf`);
        }
      }
    } catch (err) {
      console.error('[VAPOR export]', err);
      alert('Export failed: ' + err.message);
    }
  }

  // ── Toolbar injection ─────────────────────────────────────────────────────
  function attachToolbar(hostEl, resolver, filenameFn, formats) {
    if (hostEl.querySelector(':scope > .export-toolbar')) return;
    formats = formats || ['PNG', 'SVG', 'PDF'];
    const bar = document.createElement('div');
    bar.className = 'export-toolbar';
    formats.forEach(fmt => {
      const btn = document.createElement('button');
      btn.type = 'button';
      btn.className = 'export-btn';
      btn.textContent = fmt;
      btn.title = `Download as ${fmt}`;
      btn.addEventListener('click', e => {
        e.stopPropagation();
        const target = resolver();
        if (!target) return;
        const name = filenameFn ? filenameFn(target) : defaultFilename(hostEl, target);
        handleExport(target, fmt, name);
      });
      bar.appendChild(btn);
    });
    hostEl.appendChild(bar);
  }

  // Attach export toolbars to every chart/table card on the page (idempotent).
  function injectAll() {
    document.querySelectorAll('.chart-card').forEach(card => {
      const formats = card.querySelector('.table-wrap') ? ['TSV'] : ['PNG', 'SVG', 'PDF'];
      attachToolbar(card, () => resolveTarget(card), null, formats);
    });
  }

  // Used by annotation.js for dynamically-generated genome map items.
  function attachToSVGHost(hostEl, filenameFn) {
    attachToolbar(hostEl, () => resolveTarget(hostEl), filenameFn);
  }

  window.VaporExport = { injectAll, attachToSVGHost };
})();
