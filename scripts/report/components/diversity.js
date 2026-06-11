/* diversity.js — Alpha diversity ECharts bars + Beta PCoA D3 scatter */
(function () {
  'use strict';

  window.renderDiversity = function () {
    _renderAlpha();
    _setupPCoA();
  };

  // ── Alpha diversity ───────────────────────────────────────────────────────
  function _renderAlpha() {
    const rows = typeof ALPHA_DATA !== 'undefined' ? ALPHA_DATA : [];

    // Group by domain and index
    const byDomain = {};
    rows.forEach(r => {
      const dom = r.domain || 'combined';
      const idx = r.index  || 'shannon';
      const s   = r.sample || '';
      const v   = +r.value || 0;
      if (!byDomain[dom]) byDomain[dom] = {};
      if (!byDomain[dom][idx]) byDomain[dom][idx] = {};
      byDomain[dom][idx][s] = v;
    });

    const indices = ['shannon', 'simpson', 'observed', 'chao1'];
    const samples = typeof SAMPLES !== 'undefined' ? SAMPLES : [];

    ['viral', 'prok', 'combined'].forEach(dom => {
      const chartId = `div-alpha-${dom}-chart`;
      const label   = dom === 'prok' ? 'Prokaryotic' : dom.charAt(0).toUpperCase() + dom.slice(1);
      const domData = byDomain[dom] || {};

      const series = indices.filter(idx => domData[idx]).map((idx, i) => ({
        name:   idx.charAt(0).toUpperCase() + idx.slice(1),
        type:   'bar',
        color:  PAL[i % PAL.length],
        data:   samples.map(s => +(domData[idx]?.[s] || 0).toFixed(3)),
      }));

      if (!series.length) {
        const el = document.getElementById(chartId);
        if (el) el.innerHTML = '<p style="color:var(--text-muted);font-size:.8rem;padding:1rem">No alpha diversity data for this domain.</p>';
        return;
      }

      mkChart(chartId, {
        title: { text: `Alpha Diversity — ${label}` },
        tooltip: { trigger: 'axis' },
        legend: { data: series.map(s => s.name) },
        xAxis: { type: 'category', data: samples, axisLabel: { rotate: 30 } },
        yAxis: { type: 'value', name: 'Index value' },
        series,
        grid: { bottom: 70 },
      });
    });
  }

  // ── Beta diversity / PCoA (D3) ────────────────────────────────────────────
  function _setupPCoA() {
    const sel = document.getElementById('pcoa-dataset-sel');
    if (!sel) return;
    sel.addEventListener('change', _renderPCoA);
    _renderPCoA();
  }

  window.renderPCoA = _renderPCoA;

  function _renderPCoA() {
    const sel     = document.getElementById('pcoa-dataset-sel');
    const dataset = sel ? sel.value : 'viral';
    const dataMap = {
      viral:    typeof PCOA_VIRAL    !== 'undefined' ? PCOA_VIRAL    : [],
      prok:     typeof PCOA_PROK     !== 'undefined' ? PCOA_PROK     : [],
      combined: typeof PCOA_COMBINED !== 'undefined' ? PCOA_COMBINED : [],
    };
    const rows = dataMap[dataset] || [];
    const el   = document.getElementById('pcoa-svg');
    if (!el) return;
    el.innerHTML = '';

    if (!rows.length) {
      el.innerHTML = '<p style="color:var(--text-muted);padding:1rem">No PCoA data available.</p>';
      return;
    }

    const dark   = document.documentElement.dataset.theme === 'dark';
    const W      = el.clientWidth  || 900;
    const H      = 500;
    const margin = { top: 40, right: 40, bottom: 70, left: 70 };
    const iW     = W - margin.left - margin.right;
    const iH     = H - margin.top  - margin.bottom;

    const pc1Vals  = rows.map(r => +r.pc1);
    const pc2Vals  = rows.map(r => +r.pc2);
    const var1     = rows[0]?.var_pc1 ? (+rows[0].var_pc1 * 100).toFixed(1) : '?';
    const var2     = rows[0]?.var_pc2 ? (+rows[0].var_pc2 * 100).toFixed(1) : '?';

    const xScale = d3.scaleLinear()
      .domain([d3.min(pc1Vals) * 1.2, d3.max(pc1Vals) * 1.2]).range([0, iW]);
    const yScale = d3.scaleLinear()
      .domain([d3.min(pc2Vals) * 1.2, d3.max(pc2Vals) * 1.2]).range([iH, 0]);

    const color  = d3.scaleOrdinal(PAL);

    const svg = d3.select(el).append('svg').attr('width', W).attr('height', H);
    const g   = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`);

    // Grid lines
    g.append('g').attr('class', 'grid').call(
      d3.axisLeft(yScale).tickSize(-iW).tickFormat('')
    ).selectAll('line').attr('stroke', dark ? '#1e293b' : '#f1f5f9');

    // Axes
    g.append('g').attr('transform', `translate(0,${iH})`).call(d3.axisBottom(xScale))
      .selectAll('text').attr('fill', dark ? '#94a3b8' : '#64748b');

    g.append('g').call(d3.axisLeft(yScale))
      .selectAll('text').attr('fill', dark ? '#94a3b8' : '#64748b');

    // Axis labels
    g.append('text').attr('x', iW / 2).attr('y', iH + 50).attr('text-anchor', 'middle')
      .attr('fill', dark ? '#94a3b8' : '#64748b').attr('font-size', 12)
      .text(`PC1 (${var1}% variance explained)`);
    g.append('text').attr('transform', 'rotate(-90)').attr('x', -iH / 2).attr('y', -55)
      .attr('text-anchor', 'middle').attr('fill', dark ? '#94a3b8' : '#64748b').attr('font-size', 12)
      .text(`PC2 (${var2}% variance explained)`);

    // Zero lines
    g.append('line').attr('x1', xScale(0)).attr('x2', xScale(0)).attr('y1', 0).attr('y2', iH)
      .attr('stroke', dark ? '#334155' : '#cbd5e1').attr('stroke-dasharray', '4 3');
    g.append('line').attr('x1', 0).attr('x2', iW).attr('y1', yScale(0)).attr('y2', yScale(0))
      .attr('stroke', dark ? '#334155' : '#cbd5e1').attr('stroke-dasharray', '4 3');

    // Points
    const tip = d3.select(el).append('div').attr('class', 'd3-tooltip').style('display', 'none');

    g.selectAll('circle').data(rows).join('circle')
      .attr('cx', d => xScale(+d.pc1))
      .attr('cy', d => yScale(+d.pc2))
      .attr('r',  8)
      .attr('fill', d => color(d.sample))
      .attr('fill-opacity', 0.85)
      .attr('stroke', dark ? '#0f172a' : '#fff')
      .attr('stroke-width', 1.5)
      .on('mouseover', (e, d) => {
        tip.style('display', 'block')
           .style('left', (e.offsetX + 12) + 'px')
           .style('top',  (e.offsetY - 12) + 'px')
           .html(`<strong>${d.sample}</strong><br>PC1: ${(+d.pc1).toFixed(4)}<br>PC2: ${(+d.pc2).toFixed(4)}`);
      })
      .on('mouseout', () => tip.style('display', 'none'));

    // Labels next to points (each point is already labeled — no separate legend needed,
    // since a per-sample legend with ~32 entries would overflow the SVG and be clipped)
    g.selectAll('.pcoa-label').data(rows).join('text')
      .attr('class', 'pcoa-label')
      .attr('x', d => xScale(+d.pc1) + 10)
      .attr('y', d => yScale(+d.pc2) + 4)
      .attr('font-size', 10)
      .attr('fill', dark ? '#94a3b8' : '#64748b')
      .text(d => d.sample);
  }

})();
