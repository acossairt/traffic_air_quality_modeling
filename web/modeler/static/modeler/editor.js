// Live editor for the traffic model. Mirrors julia_model/traffic_configurable.jl
// so the per-edge demand-curve preview shows the same γ(t) the solver will compute.

(function () {
  // ---------------------------------------------------------------------------
  // Model math (mirrors the helpers in traffic_configurable.jl)
  // ---------------------------------------------------------------------------
  // The diurnal clock in the Julia model is a Hopf-like oscillator with
  // u(0)=1, v(0)=0, period P. Its closed-form solution is
  //   u(t) = cos(2π t / P), v(t) = sin(2π t / P).
  // Therefore v_shifted(s, t) = v*cos(2π s/P) - u*sin(2π s/P) = sin(2π (t-s)/P).

  function vShifted(s, t, period) {
    return Math.sin((2 * Math.PI * (t - s)) / period);
  }
  function fSig(x, r, x0, m) {
    return m / (1 + Math.exp(-r * (x - x0)));
  }
  function gWeib(x, a, b, c) {
    if (b === 0) return 0;
    return a * Math.exp(-Math.LN2 * Math.pow(x / b, c));
  }
  // The form exposes "peak_time" (24h clock) instead of the raw shift.
  // shift = peak_time - 6 keeps the underlying Julia math identical.
  const PEAK_TIME_OFFSET = 6;
  function demandAtTime(t, params, fromPop, period) {
    const sat = 1 - gWeib(fromPop, 1, params.demhf, 4);
    const shift = (params.peak_time ?? 0) - PEAK_TIME_OFFSET;
    const diurnal = fSig(vShifted(shift, t, period), params.dsharp, params.dur, 1);
    return sat * diurnal + params.bgd;
  }

  // ---------------------------------------------------------------------------
  // Demand curve SVG renderer (one per edge card)
  // ---------------------------------------------------------------------------
  function renderDemandCurve(svg, params, fromPop, period) {
    const W = 320, H = 160, ML = 36, MR = 10, MT = 10, MB = 24;
    const N = 200;
    const tMax = 2 * period;
    const dt = tMax / N;
    const yMax = 1.0;
    const pts = [];
    for (let i = 0; i <= N; i++) {
      const t = i * dt;
      const y = demandAtTime(t, params, fromPop, period);
      pts.push([t, y]);
    }
    const xMap = (t) => ML + ((W - ML - MR) * t) / tMax;
    const yMap = (y) => H - MB - ((H - MT - MB) * Math.min(Math.max(y, 0), yMax)) / yMax;
    const path = pts.map(([t, y], i) =>
      (i === 0 ? "M" : "L") + xMap(t).toFixed(1) + "," + yMap(y).toFixed(1)
    ).join(" ");

    const gridX = [];
    for (let h = 0; h <= tMax; h += period / 4) {
      gridX.push(`<line x1="${xMap(h)}" y1="${MT}" x2="${xMap(h)}" y2="${H - MB}" stroke="#eee"/>`);
    }
    const labels = [0, period, tMax].map(h =>
      `<text x="${xMap(h)}" y="${H - MB + 14}" font-size="10" text-anchor="middle">${h}h</text>`
    ).join("");

    svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
    svg.innerHTML = `
      <rect x="${ML}" y="${MT}" width="${W - ML - MR}" height="${H - MT - MB}" fill="#fafafa" stroke="#ccc"/>
      ${gridX.join("")}
      <path d="${path}" fill="none" stroke="#2563eb" stroke-width="2"/>
      <text x="${ML - 4}" y="${MT + 4}" font-size="10" text-anchor="end" alignment-baseline="hanging">1</text>
      <text x="${ML - 4}" y="${H - MB}" font-size="10" text-anchor="end">0</text>
      <text x="${ML - 4}" y="${(MT + (H - MB)) / 2}" font-size="10" text-anchor="end" transform="rotate(-90 ${ML - 22} ${(MT + (H - MB)) / 2})">γ</text>
      ${labels}
    `;
  }

  // ---------------------------------------------------------------------------
  // Edge card wiring
  // ---------------------------------------------------------------------------
  function readEdgeParams(card) {
    const params = {};
    card.querySelectorAll('[data-edge-param]').forEach(el => {
      params[el.dataset.edgeParam] = parseFloat(el.value || "0") || 0;
    });
    return params;
  }
  // Map from patch id (string) -> initial_pop. Built once at boot from #network-data.
  let PATCH_POPS = {};

  function readFromPop(card) {
    const sel = card.querySelector('select[data-edge-from]');
    if (!sel || !sel.value) return 0;
    return parseFloat(PATCH_POPS[sel.value] ?? 0) || 0;
  }
  function getPeriod() {
    const el = document.querySelector('input[name="period"]');
    return el ? parseFloat(el.value) || 24 : 24;
  }

  // ---------------------------------------------------------------------------
  // Corridor curves: normalized velocity & ramp-throttle vs. fractional density
  // ---------------------------------------------------------------------------
  // Render a single g(x, a, b, c) = a * exp(-ln2 * (x/b)^c) curve over
  // fractional density x = nc/(ψL) ∈ [0, 1.2]. y-axis is in real units 0..yMax.
  function renderGCurve(svg, { a, b, c, color }) {
    const W = 320, H = 160, ML = 44, MR = 10, MT = 10, MB = 24;
    const N = 200;
    const xMaxFrac = 1.2;
    const dx = xMaxFrac / N;
    const yMax = a > 0 ? a : 1;

    const pts = [];
    for (let i = 0; i <= N; i++) {
      const x = i * dx;
      const y = a * Math.exp(-Math.LN2 * Math.pow(x / b, c));
      pts.push([x, y]);
    }

    const xMap = (x) => ML + ((W - ML - MR) * x) / xMaxFrac;
    const yMap = (y) => H - MB - ((H - MT - MB) * y) / (yMax * 1.05);
    const toPath = (p) => p
      .filter(([_, y]) => Number.isFinite(y))
      .map(([x, y], i) => (i === 0 ? "M" : "L") + xMap(x).toFixed(1) + "," + yMap(y).toFixed(1))
      .join(" ");

    const gridX = [0.3, 0.6, 0.9, 1.2].map(x =>
      `<line x1="${xMap(x)}" y1="${MT}" x2="${xMap(x)}" y2="${H - MB}" stroke="#eee"/>`
    ).join("");
    const xLabels = [0, 0.3, 0.6, 0.9, 1.2].map(x =>
      `<text x="${xMap(x)}" y="${H - MB + 14}" font-size="10" text-anchor="middle">${x}</text>`
    ).join("");

    // y-axis ticks at 0, yMax/2, yMax — formatted compactly.
    const fmt = (v) => v >= 100 ? Math.round(v).toString() : v.toFixed(v < 1 ? 2 : 1);
    const yTicks = [0, yMax / 2, yMax].map(y =>
      `<text x="${ML - 4}" y="${yMap(y) + 3}" font-size="10" text-anchor="end">${fmt(y)}</text>`
    ).join("");

    svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
    svg.innerHTML = `
      <rect x="${ML}" y="${MT}" width="${W - ML - MR}" height="${H - MT - MB}" fill="#fafafa" stroke="#ccc"/>
      ${gridX}
      <line x1="${xMap(1)}" y1="${MT}" x2="${xMap(1)}" y2="${H - MB}" stroke="#bbb" stroke-dasharray="3,3"/>
      <path d="${toPath(pts)}" fill="none" stroke="${color}" stroke-width="2"/>
      ${yTicks}
      <text x="${(ML + W - MR) / 2}" y="${H - 2}" font-size="10" text-anchor="middle" fill="#555">nc / (ψL)</text>
      ${xLabels}
    `;
  }

  function initEdgeCard(card) {
    const demandSvg = card.querySelector('svg.demand-preview');
    const velSvg    = card.querySelector('svg.velocity-preview');
    const rampSvg   = card.querySelector('svg.ramp-preview');

    function updatePreview() {
      const params = readEdgeParams(card);
      if (demandSvg) renderDemandCurve(demandSvg, params, readFromPop(card), getPeriod());
      if (velSvg) renderGCurve(velSvg, {
        a: params.vff, b: params.nc_half_ff, c: params.vsharp, color: '#2563eb',
      });
      if (rampSvg) renderGCurve(rampSvg, {
        a: params.onff, b: params.on_half, c: params.onsharp, color: '#dc2626',
      });
    }

    // Sync each slider with its sibling number input, both ways.
    card.querySelectorAll('.slider-field').forEach(field => {
      const slider = field.querySelector('input[type="range"]');
      const number = field.querySelector('input[type="number"]');
      if (!slider || !number) return;
      slider.addEventListener('input', () => { number.value = slider.value; updatePreview(); });
      number.addEventListener('input', () => { slider.value = number.value; updatePreview(); });
    });

    // From-patch select changes the saturation factor, so re-render.
    const fromSel = card.querySelector('select[data-edge-from]');
    if (fromSel) fromSel.addEventListener('change', updatePreview);

    // Collapse / expand toggle.
    const toggle = card.querySelector('.edge-card-toggle');
    if (toggle) {
      toggle.addEventListener('click', () => card.classList.toggle('collapsed'));
    }

    updatePreview();
  }

  // ---------------------------------------------------------------------------
  // Network preview (patches in a circle, edges as curved arrows)
  // ---------------------------------------------------------------------------
  function renderNetwork(svg, patches, edges) {
    const W = svg.clientWidth || 600;
    const H = svg.clientHeight || 360;
    svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
    const cx = W / 2, cy = H / 2;
    const radius = Math.min(W, H) * 0.36;
    const n = patches.length;
    if (n === 0) {
      svg.innerHTML = '<text x="50%" y="50%" text-anchor="middle" fill="#888">No patches yet</text>';
      return;
    }

    const pos = {};
    patches.forEach((p, i) => {
      const a = (-Math.PI / 2) + (2 * Math.PI * i) / n;
      pos[p.id] = { x: cx + radius * Math.cos(a), y: cy + radius * Math.sin(a) };
    });

    const edgeSvg = edges.map(e => {
      const a = pos[e.from], b = pos[e.to];
      if (!a || !b) return "";
      // curve out to one side so reverse-direction edges don't overlap
      const dx = b.x - a.x, dy = b.y - a.y;
      const len = Math.hypot(dx, dy) || 1;
      const nx = -dy / len, ny = dx / len;
      const offset = 18;
      const mx = (a.x + b.x) / 2 + nx * offset;
      const my = (a.y + b.y) / 2 + ny * offset;
      // shorten endpoints so the arrow head doesn't sit inside the node circle
      const r = 28;
      const t1 = r / len, t2 = 1 - r / len;
      const sx = a.x + dx * t1, sy = a.y + dy * t1;
      const ex = a.x + dx * t2, ey = a.y + dy * t2;
      return `<path d="M${sx.toFixed(1)},${sy.toFixed(1)} Q${mx.toFixed(1)},${my.toFixed(1)} ${ex.toFixed(1)},${ey.toFixed(1)}"
                    fill="none" stroke="#666" stroke-width="1.5" marker-end="url(#arrow)"/>`;
    }).join("");

    const nodeSvg = patches.map(p => {
      const c = pos[p.id];
      return `
        <g>
          <circle cx="${c.x.toFixed(1)}" cy="${c.y.toFixed(1)}" r="26" fill="#e0ecff" stroke="#2563eb" stroke-width="1.5"/>
          <text x="${c.x.toFixed(1)}" y="${(c.y - 2).toFixed(1)}" text-anchor="middle" font-size="11" font-weight="600">${escapeHtml(p.label)}</text>
          <text x="${c.x.toFixed(1)}" y="${(c.y + 11).toFixed(1)}" text-anchor="middle" font-size="9" fill="#555">${formatPop(p.initial_pop)}</text>
        </g>`;
    }).join("");

    svg.innerHTML = `
      <defs>
        <marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto">
          <path d="M0,0 L10,5 L0,10 z" fill="#666"/>
        </marker>
      </defs>
      ${edgeSvg}
      ${nodeSvg}
    `;
  }

  function escapeHtml(s) {
    return String(s).replace(/[&<>"']/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
  }
  function formatPop(n) {
    n = Number(n) || 0;
    return n.toLocaleString(undefined, { maximumFractionDigits: 0 });
  }

  // ---------------------------------------------------------------------------
  // Add-form buttons (Django formset row cloning)
  // ---------------------------------------------------------------------------
  function initFormsetButtons() {
    document.querySelectorAll('[data-add-form]').forEach(btn => {
      btn.addEventListener('click', () => addForm(btn.dataset.addForm));
    });
  }
  function addForm(prefix) {
    const totalInput = document.querySelector(`#id_${prefix}-TOTAL_FORMS`);
    if (!totalInput) return;
    const total = parseInt(totalInput.value, 10);
    const tmplEl = document.querySelector(`#${prefix}-empty-form`);
    if (!tmplEl) return;
    const html = tmplEl.innerHTML.replaceAll('__prefix__', String(total));
    const wrapper = document.createElement('div');
    wrapper.innerHTML = html.trim();
    const newNode = wrapper.firstElementChild;
    document.querySelector(`#${prefix}-list`).appendChild(newNode);
    totalInput.value = total + 1;
    if (newNode.classList.contains('edge-card')) initEdgeCard(newNode);
  }

  // ---------------------------------------------------------------------------
  // Boot
  // ---------------------------------------------------------------------------
  document.addEventListener('DOMContentLoaded', () => {
    const netData = document.querySelector('#network-data');
    let networkData = { patches: [], edges: [] };
    if (netData) {
      try { networkData = JSON.parse(netData.textContent); }
      catch (e) { console.error("network data parse failed", e); }
    }
    PATCH_POPS = {};
    (networkData.patches || []).forEach(p => { PATCH_POPS[String(p.id)] = p.initial_pop; });

    document.querySelectorAll('#edges-list .edge-card').forEach(initEdgeCard);
    // Default existing edge cards to collapsed; freshly added cards stay open.
    document.querySelectorAll('#edges-list .edge-card').forEach(c => c.classList.add('collapsed'));
    initFormsetButtons();

    const collapseAll = document.querySelector('#edges-collapse-all');
    const expandAll = document.querySelector('#edges-expand-all');
    if (collapseAll) collapseAll.addEventListener('click', () => {
      document.querySelectorAll('#edges-list .edge-card').forEach(c => c.classList.add('collapsed'));
    });
    if (expandAll) expandAll.addEventListener('click', () => {
      document.querySelectorAll('#edges-list .edge-card').forEach(c => c.classList.remove('collapsed'));
    });

    const netSvg = document.querySelector('#network-preview');
    if (netSvg) renderNetwork(netSvg, networkData.patches || [], networkData.edges || []);
  });
})();
