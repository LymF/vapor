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
        <p class="about-subtitle">Viral And Prokaryotic mOdular pipelineR</p>

        <div class="about-section">
          <h2>About</h2>
          <p>
            VAPOR is a modular Snakemake pipeline for metagenomics and viromics, designed to
            characterise viral and prokaryotic communities from environmental samples.
            It supports <strong>short reads</strong> (Illumina PE and SE) and
            <strong>long reads</strong> (Nanopore ONT / PacBio HiFi).
          </p>
          <p style="margin-top:.6rem">
            The pipeline covers the full analysis workflow: quality control, assembly,
            deduplication, viral detection (VirSorter2, GeNomad), binning and quality
            assessment (CheckV, vRhyme, CheckM2, GUNC), taxonomy (MMseqs2, GTDB-Tk),
            host prediction (PHIST), functional annotation (Pharokka, Phold, Bakta, EggNOG-mapper),
            abundance quantification (CoverM), and alpha/beta diversity analysis.
          </p>
          <p style="margin-top:.6rem">
            VAPOR auto-detects the available runtime (<strong>Apptainer → Singularity → Conda</strong>)
            and derives container bind-mounts automatically from <code>config.yaml</code>.
          </p>
          <p style="margin-top:.6rem">
            Lytic/lysogenic lifestyle is predicted once, globally, by <strong>BACPHLIP</strong>
            over the vOTU catalog's representative sequences. Auxiliary metabolic gene
            candidates come from <strong>EggNOG-mapper</strong> run on the same representatives
            and are reported as <em>putative</em> AMGs — annotation-based AMG calls are
            known to be prone to mis-annotation and require genomic-context inspection the
            pipeline does not perform.
          </p>
        </div>

        <div class="about-section">
          <h2>Citation</h2>
          <p>
            <strong>Paper in preparation.</strong><br>
            If you use VAPOR in your research, please check the GitHub repository for the most
            up-to-date citation information.
          </p>
          <p style="margin-top:.8rem">
            <a class="about-link" href="https://github.com/LymF/vapor" target="_blank" rel="noopener">
              &#128279; github.com/LymF/vapor
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
