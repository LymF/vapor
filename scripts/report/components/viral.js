/* viral.js — Detection, binning, CheckV, taxonomy, lifestyle */
(function () {
  'use strict';

  // ── Global vOTU catalog ──────────────────────────────────────────────────
  // Per-sample counts and the total come from the SAME catalog, so they are
  // on one scale. Summing per-sample counts is not the total and never was.
  function buildVotuCatalog() {
    const cat = typeof VOTU_CATALOG !== 'undefined' ? VOTU_CATALOG : null;
    const box = document.getElementById('votu-catalog-summary');
    if (!box || !cat || !cat.n_votus) { if (box) box.innerHTML = ''; return; }
    box.innerHTML =
      `<div class="chart-card"><h3 style="margin-top:0">Global vOTU Catalog</h3>` +
      `<p><strong>${cat.n_votus.toLocaleString()}</strong> vOTUs from ` +
      `<strong>${cat.n_pool.toLocaleString()}</strong> viral contigs ` +
      `(${cat.reduction_pct}% redundancy removed). ICTV clustering: ` +
      `95% ANI + 85% AF, in a single pass over the complete set.</p></div>`;
  }

  function buildVotuPresence() {
    const pres = typeof VOTU_PRESENCE !== 'undefined' ? VOTU_PRESENCE : null;
    if (!pres) return;
    const samples = typeof SAMPLES !== 'undefined' ? SAMPLES : [];
    const assembled = samples.map(s => (pres.per_sample[s] || {}).assembled || 0);
    const recruited = samples.map(s => (pres.per_sample[s] || {}).recruited || 0);

    mkChart('votu-presence-chart', {
      title: { text: 'vOTUs per Sample — Assembly vs Read Recruitment' },
      tooltip: { trigger: 'axis' },
      legend: { data: ['Assembled', 'Recruited'], top: 28 },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 45 } },
      yAxis: { type: 'value', name: 'vOTUs' },
      grid: { bottom: 110, top: 70 },
      series: [
        { name: 'Assembled',  type: 'bar', data: assembled },
        { name: 'Recruited', type: 'bar', data: recruited },
      ],
    });

    // "Which sample(s) each virus is present in" -- capped for page weight.
    const box = document.getElementById('votu-presence-table');
    if (!box) return;
    const shown = pres.votus.slice(0, 500);
    const head = '<th>vOTU</th>' + samples.map(s => `<th>${s}</th>`).join('');
    const rows = shown.map(v => {
      const cells = samples.map(s => {
        const st = v.samples[s];
        const mark = st === 'both' ? '●' : st === 'assembled' ? '◐'
                   : st === 'recruited' ? '○' : '';
        return `<td title="${st}">${mark}</td>`;
      }).join('');
      return `<tr><td>${v.votu_id}</td>${cells}</tr>`;
    }).join('');
    box.innerHTML =
      `<div class="chart-card"><h3 style="margin-top:0">Presence by Sample</h3>` +
      `<p>● assembled and recruited · ◐ assembled only · ○ recruited only. ` +
      `Showing ${shown.length} of ${pres.votus.length} vOTUs.</p>` +
      `<div class="table-wrap"><table class="vapor-table"><thead><tr>${head}</tr></thead>` +
      `<tbody>${rows}</tbody></table></div></div>`;
  }

  window.renderViral = function () {
    const samples = typeof SAMPLES !== 'undefined' ? SAMPLES : [];

    buildVotuCatalog();
    buildVotuPresence();
    _renderDetection(samples);
    _renderBinning(samples);
    makeSampleDropdown('sample-sel-viral-tax', _renderTaxonomy, { allSamples: true });
    _renderLifestyle();
    makeSampleDropdown('sample-sel-votu-table', _renderVotuTable);
  };

  // ── Detection ─────────────────────────────────────────────────────────────
  function _renderDetection(samples) {
    const vt  = typeof VIRAL_TOOLS !== 'undefined' ? VIRAL_TOOLS : {};
    const sup = typeof SUPPORT     !== 'undefined' ? SUPPORT     : {};

    const tools = ['VirSorter2', 'GeNomad'];
    mkChart('vir-tools-chart', samplesBar({
      samples, title: 'Viral Contigs per Tool', valueName: 'Contigs',
      series: tools.map((t, i) => ({
        name: t, color: PAL[i], data: samples.map(s => (vt[s] || {})[t] || 0),
      })),
    }));

    // Support level (1 / 2 / >=3 tools) is ORDINAL -- more agreement is
    // stronger evidence -- so it takes a graded ramp, not identity hues.
    mkChart('vir-consensus-chart', samplesBar({
      samples, title: 'Tool Agreement — Consensus Filter', valueName: 'Contigs', stack: true,
      series: [
        { name: '1 tool',   color: '#94a3b8', data: samples.map(s => (sup[s] || {})[1] || 0) },
        { name: '2 tools',  color: '#5eead4', data: samples.map(s => (sup[s] || {})[2] || 0) },
        { name: '≥3 tools', color: '#0d9488', data: samples.map(s => (sup[s] || {})[3] || 0) },
      ],
    }));

    // Tool agreement as an UpSet: which COMBINATION of detectors supports each
    // contig. The old sample x tool heatmap showed per-tool totals, which cannot
    // answer "which contigs does each detector find that the other misses" -- the question a
    // consensus filter exists to settle. Aggregated across samples; the
    // per-sample stacked bar above keeps the sample-level view.
    const combosAll = {};
    const sc = typeof SUPPORT_COMBOS !== 'undefined' ? SUPPORT_COMBOS : {};
    samples.forEach(s => Object.entries(sc[s] || {}).forEach(([k, v]) => {
      combosAll[k] = (combosAll[k] || 0) + v;
    }));
    mkChart('vir-heatmap-chart', upsetPlot({
      sets: tools, combos: combosAll,
      title: 'Detector Agreement — Contigs per Tool Combination (all samples)',
      valueName: 'Contigs', unit: 'contigs',
    }));
  }

  // ── Binning & Quality ─────────────────────────────────────────────────────
  function _renderBinning(samples) {
    const cv  = typeof CHECKV     !== 'undefined' ? CHECKV     : {};
    const cvr = typeof CHECKV_VRH !== 'undefined' ? CHECKV_VRH : {};
    const vrh = typeof VRHYME     !== 'undefined' ? VRHYME     : {};
    const sup = typeof SUPPORT    !== 'undefined' ? SUPPORT    : {};

    const cvCats = ['Complete', 'High-quality', 'Medium-quality', 'Low-quality', 'Not-determined'];
    const cvCols = ['#16a34a', '#4ade80', '#fbbf24', '#d97706', '#ef4444'];

    // CheckV quality donuts — "All samples" (aggregated) or a single sample
    const ALL = 'All samples';
    function _checkvCounts(data, sampleList) {
      const counts = {};
      sampleList.forEach(s => (data[s] || []).forEach(r => {
        const q = r.checkv_quality || 'Not-determined';
        counts[q] = (counts[q] || 0) + 1;
      }));
      return counts;
    }
    // Quality tiers are an ORDINAL good->bad ladder (Complete > HQ > MQ > LQ
    // > Not-determined), and each card is really one distribution being
    // read as a composition -- a 100%-stacked horizontal bar (one row) lets
    // the reader compare segment widths directly and shares the same
    // baseline the MIMAG HQ/MQ/LQ bar (overview.js) already uses, instead
    // of asking them to eyeball two separate pie wedges/areas side by side.
    function _checkvBarOption(title, counts) {
      const total = cvCats.reduce((a, c) => a + (counts[c] || 0), 0) || 1;
      return {
        title: { text: title },
        tooltip: { trigger: 'item', formatter: p => `${p.seriesName}: ${p.value} (${(p.value / total * 100).toFixed(1)}%)` },
        legend: { top: 26 },
        grid: { top: 58, bottom: 20, left: 10, right: 10, containLabel: true },
        xAxis: { type: 'value', show: false },
        yAxis: { type: 'category', data: [''], axisLine: { show: false }, axisTick: { show: false } },
        series: cvCats.map((c, i) => ({
          name: c, type: 'bar', stack: 'q', barWidth: 46,
          itemStyle: { color: cvCols[i] },
          label: { show: true, formatter: p => p.value > 0 ? p.value : '', color: '#fff', fontSize: 11 },
          data: [counts[c] || 0],
        })),
      };
    }

    function _renderCheckvDonuts(sel) {
      const sampleList = sel === ALL ? samples : [sel];
      const label = sel === ALL ? 'all samples' : sel;

      mkChart('vir-checkv-pie-chart', _checkvBarOption(
        `CheckV — Consensus Contigs Quality (${label})`, _checkvCounts(cv, sampleList)));
      mkChart('vir-checkv-pie-vrh-chart', _checkvBarOption(
        `CheckV — vRhyme vMAGs Quality (${label})`, _checkvCounts(cvr, sampleList)));
    }

    const checkvSel = document.getElementById('sample-sel-vir-checkv');
    if (checkvSel) {
      checkvSel.innerHTML = '';
      [ALL, ...samples].forEach(s => {
        const opt = document.createElement('option');
        opt.value = opt.textContent = s;
        checkvSel.appendChild(opt);
      });
      checkvSel.addEventListener('change', () => _renderCheckvDonuts(checkvSel.value));
    }
    _renderCheckvDonuts(ALL);

    // CheckV scatter: length vs completeness
    const scatterSeries = [];
    const qualColorMap = { 'Complete': '#16a34a', 'High-quality': '#4ade80', 'Medium-quality': '#fbbf24', 'Low-quality': '#d97706', 'Not-determined': '#64748b' };
    const seenQual = new Set();
    samples.forEach(s => {
      ['consensus', 'vrhyme'].forEach(src => {
        const rows = src === 'consensus' ? (cv[s] || []) : (cvr[s] || []);
        rows.forEach(r => {
          const len  = +(r.contig_length || 0);
          const comp = +(r.completeness  || 0);
          const qual = r.checkv_quality || 'Not-determined';
          if (!len || !comp) return;
          const key  = `${qual}|${src}`;
          const sym  = src === 'consensus' ? 'circle' : 'diamond';
          if (!seenQual.has(key)) {
            seenQual.add(key);
            scatterSeries.push({ name: `${qual} (${src})`, type: 'scatter', symbolSize: 8,
              symbol: sym, itemStyle: { color: qualColorMap[qual] || '#64748b', opacity: 0.75 },
              data: [] });
          }
          const series = scatterSeries.find(s2 => s2.name === `${qual} (${src})`);
          if (series) series.data.push({ value: [len, comp], name: r.contig_id || r.contig || '' });
        });
      });
    });

    // Past VIZ.denseScatter the per-contig marks overlap into a blob, so the
    // same data is shown as a hex-binned density instead — the quality zones and
    // cutoff lines stay either way, since they are what the chart is for.
    const allPts = [];
    scatterSeries.forEach(s => s.data.forEach(d => allPts.push([Math.log10(d.value[0]), d.value[1]])));
    const hb = window.hexbin(allPts);
    const scatterOpt = hb ? {
      title: { text: `CheckV — Length vs Completeness (${hb.total.toLocaleString()} contigs, hex-binned)` },
      tooltip: { trigger: 'item',
                 formatter: p => `${Math.round(Math.pow(10, p.value[0])).toLocaleString()} bp`
                               + `<br>${p.value[1].toFixed(0)}% complete<br><strong>${p.value[2]}</strong> contigs` },
      legend: { show: false },
      xAxis: { type: 'value', name: 'Length (log10 bp)', nameLocation: 'middle', nameGap: 30,
               min: v => Math.floor(v.min * 2) / 2, max: v => Math.ceil(v.max * 2) / 2 },
      yAxis: { type: 'value', name: 'Completeness (%)', min: 0, max: 105 },
      visualMap: { min: 1, max: hb.max, calculable: false, orient: 'horizontal',
                   left: 'right', bottom: 6, itemHeight: 80, text: ['dense', 'sparse'],
                   inRange: { color: ['#ccfbf1', '#5eead4', '#0d9488', '#134e4a'] },
                   textStyle: { fontSize: 10 } },
      series: [{ type: 'scatter', data: hb.cells, symbol: 'circle', symbolSize: 9,
                 itemStyle: { opacity: 0.9 },
                 // The MQ/HQ completeness cutoffs are the reason this chart
                 // exists, so they are drawn in both the point and binned forms.
                 markLine: { silent: true, symbol: 'none', data: [
                   { yAxis: 90, lineStyle: { type: 'dashed', color: '#16a34a' },
                     label: { formatter: '≥90% HQ', color: '#16a34a', fontSize: 10 } },
                   { yAxis: 50, lineStyle: { type: 'dotted', color: '#d97706' },
                     label: { formatter: '≥50% MQ', color: '#d97706', fontSize: 10 } },
                 ] } }],
      grid: { bottom: 70, top: 50 },
    } : {
      title: { text: 'CheckV — Length vs Completeness (● consensus  ◆ vRhyme)' },
      tooltip: {
        trigger: 'item',
        formatter: p => `${p.data.name}<br>Length: ${p.data.value[0].toLocaleString()} bp<br>Completeness: ${p.data.value[1].toFixed(1)}%`,
      },
      legend: { type: 'scroll' },
      xAxis: { type: 'log', name: 'Length (bp)', min: 1000, max: 1000000, nameLocation: 'middle', nameGap: 30 },
      yAxis: { type: 'value', name: 'Completeness (%)', min: 0, max: 105 },
      series: scatterSeries,
      grid: { bottom: 60, top: 70 },
    };
    mkChart('vir-checkv-scatter-chart', scatterOpt);

    // Add quality zones (markArea + markLine) identical pattern to CheckM2
    const cvChart = window._charts['vir-checkv-scatter-chart'];
    if (cvChart && !hb) {
      cvChart.setOption({
        series: [
          ...scatterSeries,
          {
            type: 'scatter', data: [],
            markArea: {
              silent: true,
              data: [
                [{ yAxis: 90,  itemStyle: { color: 'rgba(22,163,74,0.09)' } },  { yAxis: 105 }],
                [{ yAxis: 50,  itemStyle: { color: 'rgba(217,119,6,0.07)' }  },  { yAxis: 90  }],
                [{ yAxis: 0,   itemStyle: { color: 'rgba(239,68,68,0.05)' }  },  { yAxis: 50  }],
              ],
            },
            markLine: {
              silent: true,
              data: [
                { yAxis: 90, lineStyle: { type: 'dashed', color: '#16a34a' }, label: { formatter: '≥90% HQ',  color: '#16a34a', fontSize: 10 } },
                { yAxis: 50, lineStyle: { type: 'dotted', color: '#d97706' }, label: { formatter: '≥50% MQ',  color: '#d97706', fontSize: 10 } },
              ],
            },
          },
        ],
      }, false);
    }

    // vRhyme summary
    mkChart('vir-vrhyme-chart', samplesBar({
      samples, title: 'vRhyme — vMAG Summary', valueName: 'Count',
      series: [
        { name: 'Consensus input', color: '#0891b2', data: samples.map(s => ((sup[s]||{})[3]||0) + ((sup[s]||{})[4]||0)) },
        { name: 'vMAGs formed',    color: '#0d9488', data: samples.map(s => (vrh[s] || {}).n_bins || 0) },
        { name: 'Contigs binned',  color: '#d97706', data: samples.map(s => (vrh[s] || {}).total_members || 0) },
      ],
    }));

    // Viral contig length / depth per sample. Both are heavily right-skewed, so
    // they are drawn on a log axis; distPlot picks strip vs density vs ridgeline
    // from the data (docs/REPORT_VIZ_GUIDE.md §4).
    const lenData = typeof VIRAL_LENGTHS !== 'undefined' ? VIRAL_LENGTHS : {};
    mkChart('vir-len-chart', distPlot({
      groups: samples.map(s => ({ name: s, values: (lenData[s] || []).filter(v => v > 0) })),
      title: 'Viral Contig Length Distribution (bp)',
      xName: 'Length (bp)',
      log: true,
    }));

    const depthData = typeof VIRAL_DEPTH !== 'undefined' ? VIRAL_DEPTH : {};
    mkChart('vir-depth-chart', distPlot({
      groups: samples.map(s => ({ name: s, values: (depthData[s] || []).filter(v => v > 0) })),
      title: 'Viral Contig Coverage Depth Distribution (×)',
      xName: 'Depth (×)',
      log: true,
    }));
  }

  // ── Taxonomy ──────────────────────────────────────────────────────────────
  function _renderTaxonomy(sample) {
    const allTax = typeof TAX_DATA !== 'undefined' ? TAX_DATA : [];
    const isAll  = sample === '__all__';
    const label  = isAll ? 'All samples' : sample;
    const tax    = isAll ? allTax : allTax.filter(r => r.sample === sample);

    // Source pie — covers ALL contigs (incl. diamond_custom-only hits and unknown)
    const allSrcDist = typeof VIRAL_SOURCE_DIST !== 'undefined' ? VIRAL_SOURCE_DIST : {};
    let srcDist;
    if (isAll) {
      srcDist = {};
      Object.values(allSrcDist).forEach(d => Object.entries(d).forEach(([k, v]) => { srcDist[k] = (srcDist[k] || 0) + v; }));
    } else {
      srcDist = allSrcDist[sample] || {};
    }
    // Sorted horizontal bar, not a pie: sources routinely exceed the ~4 slices
    // an angle comparison can carry, and bar length is read directly.
    const srcRows = Object.entries(srcDist).sort((a, b) => a[1] - b[1]);
    const srcTotal = srcRows.reduce((a, [, v]) => a + v, 0) || 1;
    mkChart('vir-tax-source-chart', {
      title: { text: `${label} — Classification Source` },
      tooltip: { trigger: 'item',
                 formatter: p => `${p.name}: ${p.value} (${(p.value / srcTotal * 100).toFixed(1)}%)` },
      legend: { show: false },
      xAxis: { type: 'value', name: 'Contigs', nameLocation: 'middle', nameGap: 28 },
      yAxis: { type: 'category', data: srcRows.map(r => r[0]) },
      series: [{ type: 'bar', data: srcRows.map(r => r[1]), barMaxWidth: 26,
                 itemStyle: { color: PAL[0], borderRadius: [0, 3, 3, 0] },
                 label: { show: true, position: 'right', fontSize: 10,
                          formatter: p => `${(p.value / srcTotal * 100).toFixed(0)}%` } }],
      grid: { left: 12, right: 52, bottom: 44, containLabel: true },
    });

    // Order/Family/Genus bar (top 20) — only contigs with an assignment at that rank.
    // Order is included because GeNomad-only fallback hits (no other tool found a
    // match) typically only resolve to order/class level, never family or genus —
    // without this level they'd silently disappear from every rank chart despite
    // being a large share of "Source".
    const RANK_FIELD = { order: 'final_order', family: 'final_family', genus: 'final_genus' };
    const RANK_LABEL = { order: 'Orders', family: 'Families', genus: 'Genera' };
    function renderRankBar(level) {
      const field = RANK_FIELD[level] || 'final_family';
      const count = {};
      tax.forEach(r => {
        const v = r[field] || '';
        if (!v) return;
        count[v] = (count[v] || 0) + 1;
      });
      const top = Object.entries(count).sort((a, b) => b[1] - a[1]).slice(0, 20);
      mkChart('vir-tax-family-chart', {
        title: { text: `${label} — Top Viral ${RANK_LABEL[level] || 'Families'}` },
        tooltip: { trigger: 'axis' },
        xAxis: { type: 'value', name: 'Count', nameLocation: 'middle', nameGap: 28 },
        yAxis: { type: 'category', data: top.map(x => x[0]).reverse(), axisLabel: { width: 140, overflow: 'truncate' } },
        series: [{ type: 'bar', data: top.map(x => x[1]).reverse(), itemStyle: { color: '#0d9488' } }],
        grid: { left: 160, right: 30, bottom: 50 },
      });
    }

    const rankBtns = {
      order:  document.getElementById('vir-tax-order-btn'),
      family: document.getElementById('vir-tax-family-btn'),
      genus:  document.getElementById('vir-tax-genus-btn'),
    };
    const activeLevel = Object.entries(rankBtns).find(([, btn]) => btn?.classList.contains('active'))?.[0] || 'family';
    Object.entries(rankBtns).forEach(([level, btn]) => {
      if (!btn) return;
      btn.onclick = () => {
        Object.values(rankBtns).forEach(b => b && b.classList.remove('active'));
        btn.classList.add('active');
        renderRankBar(level);
      };
    });
    if (rankBtns.family && !Object.values(rankBtns).some(b => b?.classList.contains('active'))) {
      rankBtns.family.classList.add('active');
    }
    renderRankBar(activeLevel);

    // Table — hidden when "All samples" (too many rows)
    const tableCard = document.querySelector('#vir-tax-table')?.closest('.chart-card');
    if (tableCard) tableCard.style.display = isAll ? 'none' : '';
    if (!isAll) {
      makeTable('vir-tax-table', tax, [
        { key: 'Genome',         label: 'Contig' },
        { key: 'final_order',    label: 'Order' },
        { key: 'final_family',   label: 'Family' },
        { key: 'final_genus',    label: 'Genus' },
        { key: 'Source',         label: 'Source' },
        { key: 'CheckV_quality', label: 'CheckV' },
        { key: 'Completeness',   label: 'Completeness' },
      ], {
        searchId: 'vir-tax-search',
        format: { CheckV_quality: qualBadge },
      });
    }
  }

  // ── Lifestyle (BACPHLIP, global) + Putative AMGs (eggNOG, global) ──────────
  function _renderLifestyle() {
    const life = typeof VOTU_LIFESTYLE !== 'undefined' ? VOTU_LIFESTYLE : { rows: [], counts: {}, quality_tier_used: '' };
    const dark = document.documentElement.dataset.theme === 'dark';
    const ink  = dark ? '#94a3b8' : '#64748b';

    if (!life.rows || !life.rows.length) {
      mkChart('vir-lifestyle-chart', {
        title: { text: 'vOTU Lifestyle — BACPHLIP' },
        graphic: { type: 'text', left: 'center', top: 'middle',
                   style: { text: 'No data', fill: ink, fontSize: 12 } },
      });
    } else {
      const lytic = life.counts.lytic || 0;
      const lysogenic = life.counts.lysogenic || 0;
      const total = lytic + lysogenic || 1;
      mkChart('vir-lifestyle-chart', {
        title: { text: `vOTU Lifestyle — BACPHLIP (${life.rows.length.toLocaleString()} representatives)` },
        tooltip: { trigger: 'item', formatter: p => `${p.seriesName}: ${p.value} (${(p.value / total * 100).toFixed(1)}%)` },
        legend: { top: 26 },
        grid: { top: 58, bottom: 20, left: 10, right: 10, containLabel: true },
        xAxis: { type: 'value', show: false },
        yAxis: { type: 'category', data: [''], axisLine: { show: false }, axisTick: { show: false } },
        series: [
          { name: 'Lytic', type: 'bar', stack: 'q', barWidth: 46, itemStyle: { color: '#0d9488' },
            label: { show: true, formatter: p => p.value > 0 ? p.value : '', color: '#fff', fontSize: 11 },
            data: [lytic] },
          { name: 'Lysogenic', type: 'bar', stack: 'q', barWidth: 46, itemStyle: { color: '#d97706' },
            label: { show: true, formatter: p => p.value > 0 ? p.value : '', color: '#fff', fontSize: 11 },
            data: [lysogenic] },
        ],
      });
    }

    // Which CheckV quality tier BACPHLIP was restricted to -- lysogeny
    // domains can be missing from an incomplete genome, so the caveat
    // travels with the chart rather than being buried in About.
    const chartCard = document.getElementById('vir-lifestyle-chart')?.closest('.chart-card');
    let noteEl = chartCard?.querySelector('.lifestyle-tier-note');
    if (chartCard) {
      if (life.quality_tier_used) {
        if (!noteEl) {
          noteEl = document.createElement('p');
          noteEl.className = 'lifestyle-tier-note';
          noteEl.style.cssText = 'color:var(--text-muted);font-size:.8rem;margin-top:.5rem';
          chartCard.appendChild(noteEl);
        }
        noteEl.textContent = `Lifestyle calculado apenas em vOTUs ${life.quality_tier_used} — o BACPHLIP depende de domínios de lisogenia que podem faltar em genoma incompleto.`;
      } else if (noteEl) {
        noteEl.remove();
      }
    }

    const amg = typeof PUTATIVE_AMGS !== 'undefined' ? PUTATIVE_AMGS : { rows: [] };
    makeTable('amg-table', amg.rows, [
      { key: 'votu_id',       label: 'vOTU' },
      { key: 'protein_id',    label: 'Protein' },
      { key: 'KEGG_ko',       label: 'KEGG KO' },
      { key: 'KEGG_Pathway',  label: 'KEGG Pathway' },
      { key: 'COG_category',  label: 'COG' },
      { key: 'Description',   label: 'Description' },
    ], {
      searchId: 'amg-search',
    });
  }

  // ── vOTU Table ───────────────────────────────────────────────────────────
  function _renderVotuTable(sample) {
    const all  = typeof VOTU_DATA !== 'undefined' ? VOTU_DATA : {};
    const rows = (all[sample] || []).filter(r => String(r.is_rep).toLowerCase() === 'true');
    makeTable('votu-table', rows, [
      { key: 'representative',      label: 'Representative' },
      { key: 'cluster_size',        label: 'Cluster' },
      { key: 'rep_length_bp',       label: 'Length (bp)' },
      { key: 'checkv_quality',      label: 'CheckV' },
      { key: 'checkv_completeness', label: 'Completeness' },
      { key: 'lifestyle',           label: 'Lifestyle' },
      { key: 'n_putative_AMGs',     label: 'Putative AMGs' },
      { key: 'breadth',             label: 'Breadth' },
      { key: 'rpmpm',               label: 'RPMPM' },
    ], {
      searchId: 'votu-table-search',
      format: { checkv_quality: qualBadge },
    });
  }

})();
