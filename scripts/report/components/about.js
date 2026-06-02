/* about.js — Static About tab: tool info, citation, GitHub, pipeline config */
(function () {
  'use strict';

  window.renderAbout = function () {
    const el = document.getElementById('about-content');
    if (!el) return;

    const cfg     = typeof CFG_PARAMS    !== 'undefined' ? CFG_PARAMS    : {};
    const vers    = typeof TOOL_VERSIONS !== 'undefined' ? TOOL_VERSIONS : {};
    const bench   = typeof BENCH_DATA    !== 'undefined' ? BENCH_DATA    : [];
    const samples = typeof SAMPLES       !== 'undefined' ? SAMPLES       : [];

    const cfgRows = Object.entries(cfg).map(([k, v]) =>
      `<tr><td>${k}</td><td>${v}</td></tr>`).join('');

    // Total runtime from benchmarks
    let totalSec = 0;
    bench.forEach(r => { totalSec += +(r.s || r.cpu_time || 0); });
    const runtimeStr = totalSec > 0
      ? `${(totalSec / 3600).toFixed(1)} CPU-hours`
      : 'N/A';

    // Tool version rows (2-column grid)
    const toolRows = Object.entries(vers)
      .map(([t, v]) => `<div class="tool-item"><span>${t}</span><span class="tv">${v}</span></div>`)
      .join('');

    el.innerHTML = `
      <div class="about-wrap">
        <h1>VAPOR</h1>
        <p class="about-subtitle">Virome and metagenome Analysis Pipeline for Omics Research</p>

        <div class="about-section">
          <h2>About</h2>
          <p>
            VAPOR is a comprehensive Snakemake-based pipeline for analysing viral and prokaryotic
            communities from environmental metagenomic samples. It supports both short-read
            (Illumina PE) and long-read (Nanopore / PacBio HiFi) sequencing data.
          </p>
          <p style="margin-top:.6rem">
            The pipeline integrates quality control, assembly, deduplication, viral detection
            (VirSorter2, GeNomad, VIBRANT), binning, taxonomy assignment (Diamond, vConTACT3,
            GTDB-Tk), host prediction (PHIST), and functional annotation (Pharokka, EggNOG-mapper).
          </p>
        </div>

        <div class="about-section">
          <h2>Citation</h2>
          <p>
            <strong>Paper in preparation.</strong><br>
            If you use VAPOR in your research, please check the GitHub repository for the most
            up-to-date citation information.
          </p>
          <p style="margin-top:.6rem">
            <a class="about-link" href="https://github.com/lucasmelo/vapor" target="_blank" rel="noopener">
              github.com/lucasmelo/vapor
            </a>
          </p>
        </div>

        <div class="about-section">
          <h2>This Run</h2>
          <table class="config-table">
            <tr><td>Samples</td><td>${samples.join(', ') || 'N/A'}</td></tr>
            <tr><td>Total CPU time</td><td>${runtimeStr}</td></tr>
            ${cfgRows}
          </table>
        </div>

        <div class="about-section">
          <h2>Tool Versions</h2>
          <div class="tool-grid">${toolRows || '<p style="color:var(--text-muted);font-size:.8rem">Version detection requires conda environments.</p>'}</div>
        </div>

        <div class="about-section">
          <h2>Benchmark Summary</h2>
          ${bench.length
            ? _benchTable(bench)
            : '<p style="color:var(--text-muted);font-size:.85rem">No benchmark data available.</p>'}
        </div>
      </div>`;
  };

  function _benchTable(bench) {
    const cols = Object.keys(bench[0] || {}).slice(0, 8);
    const thead = cols.map(c => `<th>${c}</th>`).join('');
    const tbody = bench.map(r =>
      `<tr>${cols.map(c => `<td>${r[c] ?? ''}</td>`).join('')}</tr>`
    ).join('');
    return `<div class="table-wrap"><table class="vapor-table"><thead><tr>${thead}</tr></thead><tbody>${tbody}</tbody></table></div>`;
  }

})();
