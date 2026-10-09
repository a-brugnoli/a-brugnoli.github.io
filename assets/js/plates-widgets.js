/* =====================================================================
   Plates course — interactive widgets (vanilla JS + Observable Plot)
   Requires d3.min.js and plot.umd.min.js loaded before this file.
   Usage: <div class="pw" data-widget="buckling"></div>
   Widgets: buckling, navier, platebeam, diskritz, ritz1, rrconv
   ===================================================================== */
(function () {
  "use strict";
  const PI = Math.PI;
  const TEAL = "#0b7a75", AMBER = "#d68a1c", RED = "#c0392b", INK = "#0e3b3f";
  const SERIES = ["#0b7a75", "#d68a1c", "#2b7bb9", "#a63d6b", "#6a8f2b", "#7a5bb8"];

  // ---------- helpers -------------------------------------------------
  const fr = (x, d = 3) => (Number.isFinite(x) ? x.toLocaleString("fr-FR", { maximumFractionDigits: d, minimumFractionDigits: d }) : "—");
  const sci = (x, d = 2) => {
    if (!Number.isFinite(x)) return "—";
    if (x === 0) return "0";
    const e = Math.floor(Math.log10(Math.abs(x)));
    const m = x / Math.pow(10, e);
    return `${m.toLocaleString("fr-FR", { maximumFractionDigits: d, minimumFractionDigits: d })}×10<sup>${e}</sup>`;
  };
  function el(tag, attrs = {}, html = "") {
    const e = document.createElement(tag);
    for (const [k, v] of Object.entries(attrs)) e.setAttribute(k, v);
    if (html) e.innerHTML = html;
    return e;
  }
  function frame(root, badge, title) {
    root.innerHTML = "";
    const head = el("div", { class: "pw-head" }, `<span class="pw-badge">${badge}</span><span class="pw-title">${title}</span>`);
    const controls = el("div", { class: "pw-controls" });
    const body = el("div", { class: "pw-body" });
    const readout = el("div", { class: "pw-readout" });
    root.append(head, controls, body, readout);
    return { controls, body, readout };
  }
  function slider(parent, label, { min, max, step, value, fmt = (v) => v }) {
    const lab = el("label");
    lab.innerHTML = `<span>${label}</span>`;
    const inp = el("input", { type: "range", min, max, step, value });
    const out = el("output");
    out.innerHTML = fmt(+value);
    lab.append(inp, out);
    parent.append(lab);
    inp.addEventListener("input", () => (out.innerHTML = fmt(+inp.value)));
    return inp;
  }
  function select(parent, label, options, value) {
    const lab = el("label");
    lab.innerHTML = `<span>${label}</span>`;
    const s = el("select");
    for (const [v, t] of options) {
      const o = el("option", { value: v }, t);
      if (v === value) o.selected = true;
      s.append(o);
    }
    lab.append(s);
    parent.append(lab);
    return s;
  }
  function checkbox(parent, label, checked) {
    const lab = el("label");
    const c = el("input", { type: "checkbox" });
    c.checked = !!checked;
    lab.append(c);
    lab.insertAdjacentHTML("beforeend", `<span>${label}</span>`);
    parent.append(lab);
    return c;
  }
  function kv(k, v) { return `<span class="kv">${k} = <b>${v}</b></span>`; }
  function panel(body) { const p = el("div", { class: "pw-panel" }); body.append(p); return p; }
  function put(p, node) { p.innerHTML = ""; p.append(node); }

  // Gauss–Legendre nodes/weights on [0,1]
  function gauss01(n) {
    const x = new Float64Array(n), w = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      let z = Math.cos(PI * (i + 0.75) / (n + 0.5)), pp = 0;
      for (let it = 0; it < 100; it++) {
        let p1 = 1, p2 = 0;
        for (let j = 1; j <= n; j++) { const p3 = p2; p2 = p1; p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j; }
        pp = n * (z * p1 - p2) / (z * z - 1);
        const z1 = z; z = z1 - p1 / pp;
        if (Math.abs(z - z1) < 1e-15) break;
      }
      x[i] = 0.5 * (1 - z); w[i] = 1 / ((1 - z * z) * pp * pp);
    }
    return { x, w };
  }
  // Dense linear solve (Gaussian elimination with partial pivoting)
  function solve(A, b) {
    const n = b.length, M = A.map((r, i) => [...r, b[i]]);
    for (let k = 0; k < n; k++) {
      let p = k; for (let i = k + 1; i < n; i++) if (Math.abs(M[i][k]) > Math.abs(M[p][k])) p = i;
      [M[k], M[p]] = [M[p], M[k]];
      for (let i = k + 1; i < n; i++) { const f = M[i][k] / M[k][k]; for (let j = k; j <= n; j++) M[i][j] -= f * M[k][j]; }
    }
    const x = new Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) { let s = M[i][n]; for (let j = i + 1; j < n; j++) s -= M[i][j] * x[j]; x[i] = s / M[i][i]; }
    return x;
  }
  // Symmetric generalized eigenproblem K v = λ B v (B SPD) — smallest pair
  function geneigMin(K, B) {
    const n = K.length;
    const L = Array.from({ length: n }, () => new Float64Array(n));
    for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) {
      let s = B[i][j]; for (let k = 0; k < j; k++) s -= L[i][k] * L[j][k];
      if (i === j) { if (s <= 0) return null; L[i][i] = Math.sqrt(s); } else L[i][j] = s / L[j][j];
    }
    // C = L^-1 K L^-T
    const Y = Array.from({ length: n }, () => new Float64Array(n)); // Y = L^-1 K
    for (let c = 0; c < n; c++) for (let i = 0; i < n; i++) { let s = K[i][c]; for (let k = 0; k < i; k++) s -= L[i][k] * Y[k][c]; Y[i][c] = s / L[i][i]; }
    const C = Array.from({ length: n }, () => new Float64Array(n)); // C = Y L^-T  -> C^T = L^-1 Y^T
    for (let r = 0; r < n; r++) for (let i = 0; i < n; i++) { let s = Y[r][i]; for (let k = 0; k < i; k++) s -= L[i][k] * C[r][k]; C[r][i] = s / L[i][i]; }
    for (let i = 0; i < n; i++) for (let j = 0; j < i; j++) { const m = 0.5 * (C[i][j] + C[j][i]); C[i][j] = C[j][i] = m; }
    // Cyclic Jacobi
    const V = Array.from({ length: n }, (_, i) => { const r = new Float64Array(n); r[i] = 1; return r; });
    for (let sweep = 0; sweep < 60; sweep++) {
      let off = 0; for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) off += C[i][j] * C[i][j];
      if (off < 1e-22) break;
      for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) {
        if (Math.abs(C[p][q]) < 1e-300) continue;
        const th = (C[q][q] - C[p][p]) / (2 * C[p][q]);
        const t = Math.sign(th || 1) / (Math.abs(th) + Math.sqrt(th * th + 1));
        const c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < n; k++) { const a = C[k][p], b = C[k][q]; C[k][p] = c * a - s * b; C[k][q] = s * a + c * b; }
        for (let k = 0; k < n; k++) { const a = C[p][k], b = C[q][k]; C[p][k] = c * a - s * b; C[q][k] = s * a + c * b; }
        for (let k = 0; k < n; k++) { const a = V[k][p], b = V[k][q]; V[k][p] = c * a - s * b; V[k][q] = s * a + c * b; }
      }
    }
    let im = 0; for (let i = 1; i < n; i++) if (C[i][i] < C[im][im]) im = i;
    const y = V.map((r) => r[im]);
    const v = new Float64Array(n); // solve L^T v = y
    for (let i = n - 1; i >= 0; i--) { let s = y[i]; for (let k = i + 1; k < n; k++) s -= L[k][i] * v[k]; v[i] = s / L[i][i]; }
    return { lambda: C[im][im], v };
  }

  // Colour bar legend (HTML)
  function cbar(lo, hi, scheme) {
    const stops = d3.range(0, 1.0001, 0.1).map((t) => scheme(t)).join(",");
    return `<div class="pw-cbar"><span>${lo}</span><div class="bar" style="background:linear-gradient(90deg,${stops})"></div><span>${hi}</span></div>`;
  }
  const RdBu = (t) => d3.interpolateRdBu(1 - t); // blue = negative, red = positive
  const Defl = (t) => d3.interpolateYlGnBu(0.08 + 0.92 * t);

  // Heat map of f(x,y) on [0,Lx]x[0,Ly] with optional contours
  function fieldPlot(f, Lx, Ly, { width = 360, scheme = "rdbu", domain, contours = 0, marks = [], title } = {}) {
    const height = Math.max(90, Math.round(width * Ly / Lx));
    const opts = {
      width, height, marginLeft: 26, marginBottom: 26, marginTop: title ? 22 : 8, marginRight: 8,
      x: { domain: [0, Lx], label: "x", ticks: 4, labelAnchor: "right" },
      y: { domain: [0, Ly], label: "y", ticks: Math.max(2, Math.round(4 * Ly / Lx)), labelAnchor: "top" },
      color: { type: "linear", scheme, domain, reverse: scheme === "rdbu" ? false : undefined },
      title,
      marks: [
        Plot.raster({ fill: f, x1: 0, x2: Lx, y1: 0, y2: Ly, width: 90, height: Math.max(20, Math.round(90 * Ly / Lx)), interpolate: "nearest", imageRendering: "auto" }),
        ...(contours ? [Plot.contour({ value: f, x1: 0, x2: Lx, y1: 0, y2: Ly, width: 70, height: Math.max(20, Math.round(70 * Ly / Lx)), stroke: "white", strokeOpacity: 0.55, strokeWidth: 0.8, thresholds: contours })] : []),
        Plot.frame({ stroke: "currentColor", strokeOpacity: 0.6 }),
        ...marks,
      ],
    };
    return Plot.plot(opts);
  }

  // =====================================================================
  // 1. Buckling coefficient of an SSSS plate under uniaxial compression
  // =====================================================================
  function buckling(root) {
    const { controls, body, readout } = frame(root, "Simulation", "Flambement d'une plaque SSSS en compression uniaxiale");
    const sA = slider(controls, "α = a/b", { min: 0.3, max: 5, step: 0.01, value: +(root.dataset.alpha || 1.0), fmt: (v) => fr(v, 2) });
    const cN = checkbox(controls, "afficher n = 2", false);
    const pL = panel(body), pR = panel(body);
    const kf = (m, n, a) => Math.pow(m / a + n * n * a / m, 2);
    const grid = d3.range(0.3, 5.0001, 0.01);
    function draw() {
      const a = +sA.value;
      let best = { m: 1, n: 1, k: Infinity };
      for (let m = 1; m <= 8; m++) for (let n = 1; n <= 2; n++) { const k = kf(m, n, a); if (k < best.k) best = { m, n, k }; }
      const lines = [];
      for (let m = 1; m <= 5; m++) lines.push(Plot.line(grid.map((x) => [x, kf(m, 1, x)]), { clip: true, stroke: SERIES[m - 1], strokeWidth: m === best.m ? 2.6 : 1.4, strokeOpacity: m === best.m ? 1 : 0.6 }));
      if (cN.checked) lines.push(Plot.line(grid.map((x) => [x, kf(1, 2, x)]), { clip: true, stroke: "currentColor", strokeDasharray: "4,3", strokeOpacity: 0.6 }));
      const env = grid.map((x) => [x, d3.min(d3.range(1, 9), (m) => kf(m, 1, x))]);
      put(pL, Plot.plot({
        width: 440, height: 300, marginLeft: 40, marginBottom: 36,
        x: { domain: [0.3, 5], label: "α = a/b →" }, y: { domain: [0, 12], label: "↑ k", grid: true },
        marks: [
          Plot.areaY(env, { clip: true, x: (d) => d[0], y1: 0, y2: (d) => d[1], fill: TEAL, fillOpacity: 0.07 }),
          ...lines,
          Plot.ruleY([4], { stroke: RED, strokeDasharray: "2,3" }),
          Plot.ruleX([a], { stroke: "currentColor", strokeOpacity: 0.35 }),
          Plot.dot([[a, best.k]], { r: 5, fill: AMBER, stroke: "white" }),
          Plot.text([1, 2, 3, 4, 5].map((m) => [m, kf(m, 1, m) + 0.6, `m=${m}`]), { x: (d) => d[0], y: (d) => d[1], text: (d) => d[2], fill: (d, i) => SERIES[i], fontSize: 10 }),
        ],
      }));
      const b = 1, L = a * b;
      put(pR, fieldPlot((x, y) => Math.sin(best.m * PI * x / L) * Math.sin(best.n * PI * y / b), L, b,
        { width: 360, scheme: "rdbu", domain: [-1, 1], contours: 7, title: `Mode critique (m, n) = (${best.m}, ${best.n})`,
          marks: [Plot.arrow([[-0.12 * L, b / 2, 0.0, b / 2], [1.12 * L, b / 2, L, b / 2]], { x1: (d) => d[0], y1: (d) => d[1], x2: (d) => d[2], y2: (d) => d[3], stroke: RED, strokeWidth: 2, headLength: 7 })] }));
      readout.innerHTML = kv("mode critique", `(${best.m}, ${best.n})`) + kv("k<sub>min</sub>", fr(best.k, 3)) +
        kv("N<sub>cr</sub>", `${fr(best.k, 2)} π²D/b²`) + kv("transition suivante α", fr(Math.sqrt(best.m * (best.m + 1)), 3));
    }
    sA.addEventListener("input", draw); cN.addEventListener("change", draw); draw();
  }

  // =====================================================================
  // 2. Navier series for an SSSS plate (uniform or point load)
  // =====================================================================
  function navier(root) {
    const { controls, body, readout } = frame(root, "Simulation", "Solution de Navier — plaque rectangulaire appuyée");
    const sLoad = select(controls, "charge", [["u", "uniforme p₀"], ["p", "ponctuelle F₀"]], root.dataset.load || "u");
    const sN = slider(controls, "termes K", { min: 1, max: 15, step: 1, value: 3, fmt: (v) => `${v}` });
    const sB = slider(controls, "b/a", { min: 0.5, max: 3, step: 0.05, value: 1, fmt: (v) => fr(v, 2) });
    const sX = slider(controls, "x*/a", { min: 0.05, max: 0.95, step: 0.01, value: 0.5, fmt: (v) => fr(v, 2) });
    const sY = slider(controls, "y*/b", { min: 0.05, max: 0.95, step: 0.01, value: 0.5, fmt: (v) => fr(v, 2) });
    const pL = panel(body), pR = panel(body);
    const note = el("div", { class: "pw-note" });
    root.append(note);
    // w normalised: uniform -> p0 a^4 / D ; point -> F0 a^2 / D (positive = downwards)
    function coeffs(load, K, a, b, xs, ys) {
      const c = [];
      for (let m = 1; m <= 2 * K - 1; m++) for (let n = 1; n <= 2 * K - 1; n++) {
        const den = Math.pow((m / a) ** 2 + (n / b) ** 2, 2);
        if (load === "u") { if (m % 2 && n % 2) c.push([m, n, 16 / (Math.pow(PI, 6) * m * n * den)]); }
        else c.push([m, n, 4 * Math.sin(m * PI * xs / a) * Math.sin(n * PI * ys / b) / (a * b * Math.pow(PI, 4) * den)]);
      }
      return c;
    }
    const evalW = (c, a, b, x, y) => { let s = 0; for (const [m, n, w] of c) s += w * Math.sin(m * PI * x / a) * Math.sin(n * PI * y / b); return s; };
    function draw() {
      const load = sLoad.value, K = +sN.value, a = 1, b = +sB.value;
      const xs = +sX.value * a, ys = +sY.value * b;
      sX.parentElement.style.display = sY.parentElement.style.display = load === "p" ? "" : "none";
      const c = coeffs(load, K, a, b, xs, ys);
      const px = load === "u" ? a / 2 : xs, py = load === "u" ? b / 2 : ys;
      const wref = evalW(coeffs(load, load === "u" ? 40 : 60, a, b, xs, ys), a, b, px, py);
      const wK = evalW(c, a, b, px, py);
      // field (sampled on a grid, cached sin tables)
      const NX = 64, NY = Math.max(16, Math.round(64 * b / a));
      const G = new Float64Array(NX * NY);
      let gmax = 1e-300;
      const sx = c.map(([m]) => Float64Array.from({ length: NX }, (_, i) => Math.sin(m * PI * ((i + 0.5) / NX))));
      const sy = c.map(([, n]) => Float64Array.from({ length: NY }, (_, j) => Math.sin(n * PI * ((j + 0.5) / NY))));
      for (let k = 0; k < c.length; k++) { const w = c[k][2]; for (let j = 0; j < NY; j++) { const t = w * sy[k][j]; for (let i = 0; i < NX; i++) G[j * NX + i] += t * sx[k][i]; } }
      for (const v of G) gmax = Math.max(gmax, v);
      const f = (x, y) => { // bilinear interpolation on the cell-centred grid
        const gx = Math.max(0, Math.min(NX - 1.001, x / a * NX - 0.5)), gy = Math.max(0, Math.min(NY - 1.001, y / b * NY - 0.5));
        const i = Math.floor(gx), j = Math.floor(gy), tx = gx - i, ty = gy - j;
        const v = (1 - tx) * (1 - ty) * G[j * NX + i] + tx * (1 - ty) * G[j * NX + i + 1] + (1 - tx) * ty * G[(j + 1) * NX + i] + tx * ty * G[(j + 1) * NX + i + 1];
        return v / gmax;
      };
      put(pL, fieldPlot(f, a, b, {
        width: b > 1.6 ? 240 : 330, scheme: "ylgnbu", domain: [0, 1], contours: 8,
        title: load === "u" ? "Flèche w / w_max" : "Flèche w / w_max — cliquez pour déplacer F₀",
        marks: load === "p" ? [Plot.dot([[xs, ys]], { r: 5, fill: RED, stroke: "white" })] : [Plot.dot([[a / 2, b / 2]], { r: 3, fill: RED })],
      }));
      const svg = pL.querySelector("svg");
      if (load === "p" && svg) {
        svg.style.cursor = "crosshair";
        svg.addEventListener("click", (ev) => {
          const r = svg.getBoundingClientRect(), vb = svg.viewBox.baseVal;
          const X = (ev.clientX - r.left) * vb.width / r.width, Y = (ev.clientY - r.top) * vb.height / r.height;
          const W = vb.width - 26 - 8, H = vb.height - 26 - 22;
          const fx = (X - 26) / W, fy = 1 - (Y - 22) / H;
          if (fx > 0 && fx < 1 && fy > 0 && fy < 1) { sX.value = fx.toFixed(2); sY.value = fy.toFixed(2); sX.dispatchEvent(new Event("input")); sY.dispatchEvent(new Event("input")); }
        });
      }
      // convergence curve
      const conv = d3.range(1, 16).map((k) => {
        const wk = evalW(coeffs(load, k, a, b, xs, ys), a, b, px, py);
        return { k, err: Math.max(Math.abs(wk - wref) / Math.abs(wref), 1e-12) };
      });
      put(pR, Plot.plot({
        width: 330, height: 250, marginLeft: 50, marginBottom: 36, title: load === "u" ? "Erreur relative sur w(a/2, b/2)" : "Erreur relative sur w(x*, y*)",
        x: { label: "K (m, n ≤ 2K−1) →", domain: [1, 15] }, y: { type: "log", label: "↑ erreur", grid: true, tickFormat: (d) => { const e = Math.round(Math.log10(d)); return Math.abs(Math.log10(d) - e) < 1e-6 ? `1e${e}` : ""; } },
        marks: [Plot.line(conv, { x: "k", y: "err", stroke: TEAL }), Plot.dot(conv, { x: "k", y: "err", r: 2.5, fill: TEAL }), Plot.dot(conv.filter((d) => d.k === K), { x: "k", y: "err", r: 6, fill: AMBER, stroke: "white" })],
      }));
      const unit = load === "u" ? "p₀a⁴/D" : "F₀a²/D";
      readout.innerHTML = kv("termes non nuls", c.filter((d) => d[2] !== 0).length) + kv("w", `${fr(wK, 6)} ${unit}`) + kv("référence", `${fr(wref, 6)} ${unit}`) +
        kv("erreur", `${fr(100 * Math.abs(wK - wref) / Math.abs(wref), 3)} %`);
      note.innerHTML = load === "u"
        ? "Seuls les termes m, n impairs sont non nuls : la charge uniforme est symétrique. La série converge en 1/(mn)⁵ : quelques termes suffisent pour la flèche."
        : "Sous une charge ponctuelle, les moments sont singuliers ; la flèche converge mais plus lentement (≈ 1/(m²+n²)²).";
    }
    for (const s of [sLoad, sN, sB, sX, sY]) s.addEventListener("input", draw);
    draw();
  }

  // =====================================================================
  // 3. Plate vs beam model (PC1, exercise 1)
  // =====================================================================
  function platebeam(root) {
    const { controls, body, readout } = frame(root, "Visualisation", "Déformée : modèle plaque vs modèle poutre");
    body.classList.add("single");
    const sK = slider(controls, "amplification K/E", { min: 0, max: 0.12, step: 0.002, value: 0.08, fmt: (v) => fr(v, 3) });
    const sNu = slider(controls, "ν", { min: 0, max: 0.5, step: 0.01, value: 0.33, fmt: (v) => fr(v, 2) });
    const p = panel(body);
    const A = 4, B = 1.5;
    function draw() {
      const k = +sK.value, nu = +sNu.value;
      const plate = (x, y) => [x + k * (x * y + B * x), y - k * (x * x / 2 + nu * y * y / 2 + nu * B * y)];
      const beam = (x, y) => [x + k * (x * y + B * x), y - k * x * x / 2];
      const outline = (F) => [
        ...d3.range(0, A + 1e-9, A / 60).map((x) => F(x, -B)), ...d3.range(-B, B + 1e-9, 2 * B / 30).map((y) => F(A, y)),
        ...d3.range(A, -1e-9, -A / 60).map((x) => F(x, B)), ...d3.range(B, -B - 1e-9, -2 * B / 30).map((y) => F(0, y)),
      ];
      const gridLines = [];
      for (const y of d3.range(-B, B + 1e-9, B / 3)) gridLines.push(d3.range(0, A + 1e-9, A / 40).map((x) => plate(x, y)));
      for (const x of d3.range(0, A + 1e-9, A / 6)) gridLines.push(d3.range(-B, B + 1e-9, B / 15).map((y) => plate(x, y)));
      put(p, Plot.plot({
        width: 680, height: 330, marginLeft: 36, marginBottom: 30,
        x: { domain: [-0.4, 6.0], label: "x" }, y: { domain: [-2.6, 1.9], label: "y" }, aspectRatio: 1,
        marks: [
          Plot.line([[0, -B], [A, -B], [A, B], [0, B], [0, -B]], { stroke: "currentColor", strokeOpacity: 0.35, strokeDasharray: "4,3" }),
          ...gridLines.map((l) => Plot.line(l, { stroke: TEAL, strokeOpacity: 0.25 })),
          Plot.line(outline(plate), { stroke: TEAL, strokeWidth: 2.2 }),
          Plot.line(outline(beam), { stroke: AMBER, strokeWidth: 2.2, strokeDasharray: "6,4" }),
          Plot.dot([[0, 0]], { r: 4, fill: RED }),
          Plot.text([[A * 1.02, B + 0.35, "— plaque"], [A * 1.02, B + 0.05, "- - poutre"]], { x: (d) => d[0], y: (d) => d[1], text: (d) => d[2], fill: (d, i) => (i ? AMBER : TEAL), textAnchor: "start", fontSize: 12, fontWeight: 600 }),
        ],
      }));
      const dv = k * nu * (B * B / 2 + B * B); // |v_plate - v_beam| max at y=-b? compute numerically
      let maxd = 0; for (const y of d3.range(-B, B + 1e-9, 0.01)) maxd = Math.max(maxd, Math.abs(k * (nu * y * y / 2 + nu * B * y)));
      readout.innerHTML = kv("θ<sub>z</sub>(a)", `−${fr(k * A, 3)} rad (identique)`) + kv("écart max sur v (effet Poisson)", fr(maxd, 4)) + kv("v(a,0) flèche", fr(-k * A * A / 2, 3));
      void dv;
    }
    sK.addEventListener("input", draw); sNu.addEventListener("input", draw); draw();
  }

  // =====================================================================
  // 4. Clamped disk under a central point load — Ritz with N functions
  // =====================================================================
  function diskritz(root) {
    const { controls, body, readout } = frame(root, "Simulation", "Disque encastré, force ponctuelle — méthode de Ritz");
    const sN = slider(controls, "fonctions N", { min: 1, max: 6, step: 1, value: 1, fmt: (v) => `${v}` });
    const pL = panel(body), pR = panel(body);
    const note = el("div", { class: "pw-note" }, "Base φ<sub>k</sub>(r) = (1 − r²/a²)<sup>k+1</sup>, k = 1…N. Unités : a = 1, D = 1, F = 1.");
    root.append(note);
    const q = gauss01(48);
    const lap = (p, r) => -4 * p * Math.pow(1 - r * r, p - 1) + 4 * p * (p - 1) * r * r * Math.pow(1 - r * r, Math.max(p - 2, 0)) * (p >= 2 ? 1 : 0);
    const exact = (r) => (r <= 0 ? 1 / (16 * PI) : (1 - r * r + 2 * r * r * Math.log(r)) / (16 * PI));
    function ritz(N) {
      const K = Array.from({ length: N }, () => new Array(N).fill(0)), f = new Array(N).fill(1);
      for (let i = 0; i < N; i++) for (let j = 0; j < N; j++) {
        let s = 0; for (let g = 0; g < q.x.length; g++) { const r = q.x[g]; s += q.w[g] * lap(i + 2, r) * lap(j + 2, r) * r; }
        K[i][j] = 2 * PI * s;
      }
      return solve(K, f);
    }
    const wN = (c, r) => c.reduce((s, ci, i) => s + ci * Math.pow(1 - r * r, i + 2), 0);
    const ratios = d3.range(1, 9).map((n) => ({ n, ratio: wN(ritz(n), 0) / exact(0) }));
    function draw() {
      const N = +sN.value, c = ritz(N);
      const rr = d3.range(-1, 1.0001, 0.01);
      put(pL, Plot.plot({
        width: 380, height: 260, marginLeft: 46, marginBottom: 34, title: "Profil de flèche w(r) (vers le bas)",
        x: { label: "r / a →" }, y: { label: "↑ w · D/(F a²)", reverse: true, grid: true },
        marks: [
          Plot.line(rr.map((r) => [r, exact(Math.abs(r))]), { stroke: "currentColor", strokeWidth: 2 }),
          Plot.line(rr.map((r) => [r, wN(c, Math.abs(r))]), { stroke: TEAL, strokeWidth: 2.4, strokeDasharray: "6,3" }),
          Plot.text([[0.35, exact(0) * 0.12, "exacte"], [0.35, exact(0) * 0.3, `Ritz N = ${N}`]], { x: (d) => d[0], y: (d) => d[1], text: (d) => d[2], fill: (d, i) => (i ? TEAL : "currentColor"), textAnchor: "start", fontWeight: 600 }),
          Plot.ruleY([0]),
        ],
      }));
      put(pR, Plot.plot({
        width: 300, height: 260, marginLeft: 46, marginBottom: 34, title: "w_Ritz(0) / w_exacte(0)",
        x: { label: "N →", ticks: 8 }, y: { domain: [0.7, 1.0], grid: true, label: null },
        marks: [Plot.ruleY([1], { stroke: RED, strokeDasharray: "3,3" }), Plot.line(ratios, { x: "n", y: "ratio", stroke: TEAL }),
          Plot.dot(ratios, { x: "n", y: "ratio", r: 3, fill: TEAL }), Plot.dot(ratios.filter((d) => d.n === N), { x: "n", y: "ratio", r: 6, fill: AMBER, stroke: "white" })],
      }));
      readout.innerHTML = kv("w(0) Ritz", `${fr(wN(c, 0) * 64 * PI, 4)} · Fa²/(64πD)`) + kv("w(0) exacte", "4 · Fa²/(64πD)") +
        kv("rapport", `${fr(100 * wN(c, 0) / exact(0), 2)} %`) + kv("c", "[" + c.map((v) => fr(v * 64 * PI, 3)).join(" ; ") + "] · F a²/(64πD)");
    }
    sN.addEventListener("input", draw); draw();
  }

  // =====================================================================
  // 5. One-term Ritz with separated envelopes (any C/S/F combination)
  // =====================================================================
  // polynomial helpers with integer coefficients (ascending powers)
  const pmul = (a, b) => { const r = new Array(a.length + b.length - 1).fill(0); a.forEach((x, i) => b.forEach((y, j) => (r[i + j] += x * y))); return r; };
  const pder = (a) => a.slice(1).map((x, i) => x * (i + 1));
  const gcd = (a, b) => { a = Math.abs(a); b = Math.abs(b); while (b) [a, b] = [b, a % b]; return a || 1; };
  function pint01(a) { // exact rational ∫_0^1
    let num = 0, den = 1;
    a.forEach((c, i) => { const d = i + 1; num = num * d + c * den; den *= d; const g = gcd(num, den); num /= g; den /= g; });
    return [num, den];
  }
  const fracStr = ([n, d]) => (n === 0 ? "0" : d === 1 ? `${n}` : `${n < 0 ? "−" : ""}${Math.abs(n)}/${d}`);
  const R = { C: 2, S: 1, F: 0 };
  function envelope(b0, b1) { // ξ^r0 (1-ξ)^r1
    let p = [1];
    for (let i = 0; i < R[b0]; i++) p = pmul(p, [0, 1]);
    for (let i = 0; i < R[b1]; i++) p = pmul(p, [1, -1]);
    return p;
  }
  function integrals1D(b0, b1, trig) {
    if (trig) return { g: (t) => Math.sin(PI * t), S: [2 / PI, "2/π"], C: [0.5, "1/2"], B: [PI * PI / 2, "π²/2"], Dd: [-PI * PI / 2, "−π²/2"], A: [Math.pow(PI, 4) / 2, "π⁴/2"], name: "sin πξ" };
    const g = envelope(b0, b1), g1 = pder(g), g2 = pder(g1);
    const I = (p) => { const f = pint01(p); return [f[0] / f[1], fracStr(f)]; };
    const term = (c, i) => (c === 0 ? "" : `${c < 0 ? " − " : " + "}${Math.abs(c) === 1 && i > 0 ? "" : Math.abs(c)}${i ? "ξ" + (i > 1 ? "<sup>" + i + "</sup>" : "") : ""}`);
    let name = g.map(term).join("").replace(/^ \+ /, "").replace(/^ − /, "−");
    return { g: (t) => g.reduce((s, c, i) => s + c * Math.pow(t, i), 0), S: I(g), C: I(pmul(g, g)), B: I(pmul(g1, g1)), Dd: I(pmul(g, g2)), A: I(pmul(g2, g2)), name: name || "1" };
  }
  function ritz1(root) {
    const { controls, body, readout } = frame(root, "Calculateur", "Ritz à un terme : w = c · g(x/a) · g(y/b)");
    const opts = [["C", "C (encastré)"], ["S", "S (appui simple)"], ["F", "F (libre)"]];
    const s0 = select(controls, "x = 0", opts, root.dataset.bc ? root.dataset.bc[0] : "S");
    const s1 = select(controls, "x = a", opts, root.dataset.bc ? root.dataset.bc[1] : "S");
    const s2 = select(controls, "y = 0", opts, root.dataset.bc ? root.dataset.bc[2] : "S");
    const s3 = select(controls, "y = b", opts, root.dataset.bc ? root.dataset.bc[3] : "S");
    const cT = checkbox(controls, "sinus si S–S", false);
    const sA = slider(controls, "a/b", { min: 0.5, max: 3, step: 0.05, value: 1, fmt: (v) => fr(v, 2) });
    const sNu = slider(controls, "ν", { min: 0, max: 0.5, step: 0.01, value: 0.3, fmt: (v) => fr(v, 2) });
    const pL = panel(body), pR = panel(body);
    const REF = { SSSS: 0.00406, CCCC: 0.00126 }; // square plate, uniform load (Timoshenko)
    function draw() {
      const bx = [s0.value, s1.value], by = [s2.value, s3.value];
      const a = +sA.value, b = 1, nu = +sNu.value;
      const X = integrals1D(bx[0], bx[1], cT.checked && bx.join("") === "SS");
      const Y = integrals1D(by[0], by[1], cT.checked && by.join("") === "SS");
      const Ax = X.A[0] / a ** 3, Bx = X.B[0] / a, Cx = X.C[0] * a, Dx = X.Dd[0] / a, Sx = X.S[0] * a;
      const Ay = Y.A[0] / b ** 3, By = Y.B[0] / b, Cy = Y.C[0] * b, Dy = Y.Dd[0] / b, Sy = Y.S[0] * b;
      const K = Ax * Cy + Cx * Ay + 2 * nu * Dx * Dy + 2 * (1 - nu) * Bx * By;
      const F = Sx * Sy;
      const singular = K < 1e-12;
      const c = singular ? NaN : F / K;
      // max deflection over the plate
      let wmax = 0, at = [0, 0];
      for (let i = 0; i <= 40; i++) for (let j = 0; j <= 40; j++) { const v = c * X.g(i / 40) * Y.g(j / 40); if (Math.abs(v) > Math.abs(wmax)) { wmax = v; at = [i / 40, j / 40]; } }
      const ncr = singular || Bx * Cy === 0 ? NaN : K / (Bx * Cy);
      const key = bx.join("") + by.join("");
      const exactK = key === "SSSS" ? d3.min(d3.range(1, 9), (m) => Math.pow(m / a + a / m, 2)) * PI * PI : NaN;
      // edge marks
      const style = (t) => (t === "C" ? { strokeWidth: 5, stroke: "currentColor" } : t === "S" ? { strokeWidth: 2.5, stroke: "currentColor", strokeDasharray: "7,4" } : { strokeWidth: 1.2, stroke: "currentColor", strokeDasharray: "1,4" });
      const edges = [
        Plot.line([[0, 0], [0, b]], style(bx[0])), Plot.line([[a, 0], [a, b]], style(bx[1])),
        Plot.line([[0, 0], [a, 0]], style(by[0])), Plot.line([[0, b], [a, b]], style(by[1])),
      ];
      const gmax = singular ? 1 : Math.abs(wmax / c) || 1;
      put(pL, fieldPlot((x, y) => (singular ? 0 : X.g(x / a) * Y.g(y / b) / gmax), a, b, { width: a > 1.7 ? 360 : 300, scheme: "ylgnbu", domain: [0, 1], contours: 8, marks: edges, title: "Fonction d'essai g(x/a)·g(y/b)" }));
      const row = (n, I) => `<tr><td>${n}</td><td>${I.name}</td><td>${I.S[1]}</td><td>${I.C[1]}</td><td>${I.B[1]}</td><td>${I.Dd[1]}</td><td>${I.A[1]}</td></tr>`;
      pR.innerHTML = `<table class="table table-sm"><thead><tr><th></th><th>g(ξ)</th><th>∫g</th><th>∫g²</th><th>∫g′²</th><th>∫gg″</th><th>∫g″²</th></tr></thead><tbody>${row("x", X)}${row("y", Y)}</tbody></table>
        <div class="pw-note">Trait épais : encastré · tirets : appui simple · pointillés : libre.<br>K = D[A<sub>x</sub>C<sub>y</sub> + C<sub>x</sub>A<sub>y</sub> + 2νD<sub>x</sub>D<sub>y</sub> + 2(1−ν)B<sub>x</sub>B<sub>y</sub>], F = p₀S<sub>x</sub>S<sub>y</sub>.</div>`;
      if (singular) { readout.innerHTML = `<span class="kv" style="color:${RED}">K = 0 : mode rigide (aucune condition essentielle suffisante) — choisir d'autres bords.</span>`; return; }
      let ref = "";
      if (Math.abs(a - 1) < 1e-9 && REF[key] && !(cT.checked && key !== "SSSS")) ref = kv("référence", `${fr(REF[key], 5)} p₀a⁴/D`) + kv("écart", `${fr(100 * (Math.abs(wmax) - REF[key]) / REF[key], 1)} %`);
      readout.innerHTML = kv("K/D", fr(K, 4)) + kv("F/p₀", fr(F, 5)) + kv("w<sub>max</sub>", `${fr(Math.abs(wmax), 5)} p₀b⁴/D`) + ref +
        (Number.isFinite(ncr) ? kv("Rayleigh N<sub>cr</sub> (compression x)", `${fr(ncr, 2)} D/b²`) : "") +
        (Number.isFinite(exactK) ? kv("exact SSSS", `${fr(exactK, 2)} D/b²`) : "");
    }
    for (const s of [s0, s1, s2, s3, cT, sA, sNu]) s.addEventListener("input", draw);
    draw();
  }

  // =====================================================================
  // 6. Rayleigh–Ritz convergence for SSSS buckling (polynomial basis)
  // =====================================================================
  function rrconv(root) {
    const { controls, body, readout } = frame(root, "Simulation", "Rayleigh–Ritz : flambement SSSS, base ξ(1−ξ)ξᵐ · η(1−η)ηⁿ");
    const sM = slider(controls, "M (base M×M)", { min: 1, max: 6, step: 1, value: 2, fmt: (v) => `${v} → ${v * v} ddl` });
    const sA = slider(controls, "α = a/b", { min: 0.5, max: 3, step: 0.05, value: +(root.dataset.alpha || 1), fmt: (v) => fr(v, 2) });
    const sNu = slider(controls, "ν", { min: 0, max: 0.5, step: 0.05, value: 0.3, fmt: (v) => fr(v, 2) });
    const pL = panel(body), pR = panel(body);
    const q = gauss01(24);
    function mats1D(M, L) { // X_m(x) = ξ(1-ξ)ξ^m, ξ = x/L
      const vals = [], d1 = [], d2 = [];
      for (let m = 0; m < M; m++) {
        const p = pmul([0, 1, -1], [...new Array(m).fill(0), 1]), p1 = pder(p), p2 = pder(p1);
        const ev = (pp, t) => pp.reduce((s, c, i) => s + c * Math.pow(t, i), 0);
        vals.push(Array.from(q.x, (t) => ev(p, t))); d1.push(Array.from(q.x, (t) => ev(p1, t) / L)); d2.push(Array.from(q.x, (t) => ev(p2, t) / (L * L)));
      }
      const mk = (f, g) => Array.from({ length: M }, (_, i) => Array.from({ length: M }, (_, j) => { let s = 0; for (let k = 0; k < q.x.length; k++) s += q.w[k] * f[i][k] * g[j][k]; return s * L; }));
      return { A0: mk(vals, vals), A1: mk(d1, d1), A2: mk(d2, d2), P: mk(d2, vals), vals };
    }
    function solveRR(M, a, nu) {
      const X = mats1D(M, a), Y = mats1D(M, 1);
      const idx = []; for (let m = 0; m < M; m++) for (let n = 0; n < M; n++) idx.push([m, n]);
      const N = idx.length;
      const K = Array.from({ length: N }, () => new Array(N).fill(0)), B = Array.from({ length: N }, () => new Array(N).fill(0));
      idx.forEach(([m, n], I) => idx.forEach(([k, l], J) => {
        K[I][J] = X.A2[m][k] * Y.A0[n][l] + X.A0[m][k] * Y.A2[n][l] + nu * (X.P[m][k] * Y.P[l][n] + X.P[k][m] * Y.P[n][l]) + 2 * (1 - nu) * X.A1[m][k] * Y.A1[n][l];
        B[I][J] = X.A1[m][k] * Y.A0[n][l];
      }));
      const r = geneigMin(K, B);
      return { ...r, idx, N };
    }
    const kExact = (a) => d3.min(d3.range(1, 9), (m) => Math.pow(m / a + a / m, 2));
    function draw() {
      const M = +sM.value, a = +sA.value, nu = +sNu.value;
      const ke = kExact(a);
      const conv = d3.range(1, 7).map((m) => { const r = solveRR(m, a, nu); return { m, dofs: m * m, k: r ? r.lambda / (PI * PI) : NaN }; });
      const cur = solveRR(M, a, nu);
      put(pL, Plot.plot({
        width: 340, height: 260, marginLeft: 48, marginBottom: 36, title: "k Ritz / k exact",
        x: { label: "ddl (M²) →", type: "log", ticks: [1, 4, 9, 16, 25, 36], tickFormat: (d) => `${d}` }, y: { grid: true, label: null, domain: [0.98, Math.max(1.15, d3.max(conv, (d) => d.k / ke) * 1.01)] },
        marks: [Plot.ruleY([1], { stroke: RED, strokeDasharray: "3,3" }), Plot.line(conv, { x: "dofs", y: (d) => d.k / ke, stroke: TEAL }),
          Plot.dot(conv, { x: "dofs", y: (d) => d.k / ke, r: 3, fill: TEAL }), Plot.dot(conv.filter((d) => d.m === M), { x: "dofs", y: (d) => d.k / ke, r: 6, fill: AMBER, stroke: "white" })],
      }));
      // mode shape
      const v = cur.v, idx = cur.idx;
      const ev = (m, t) => t * (1 - t) * Math.pow(t, m);
      let mx = 0; const samp = (x, y) => { let s = 0; idx.forEach(([m, n], I) => (s += v[I] * ev(m, x / a) * ev(n, y))); return s; };
      for (let i = 0; i <= 30; i++) for (let j = 0; j <= 30; j++) mx = Math.max(mx, Math.abs(samp(a * i / 30, j / 30)));
      put(pR, fieldPlot((x, y) => samp(x, y) / mx, a, 1, { width: a > 1.7 ? 360 : 280, scheme: "rdbu", domain: [-1, 1], contours: 7, title: "Mode approché" }));
      readout.innerHTML = kv("k Ritz", fr(cur.lambda / (PI * PI), 4)) + kv("k exact", fr(ke, 4)) + kv("écart", `+${fr(100 * (cur.lambda / (PI * PI) / ke - 1), 3)} %`) +
        kv("N<sub>cr</sub>", `${fr(cur.lambda, 2)} D/b²`);
    }
    for (const s of [sM, sA, sNu]) s.addEventListener("input", draw);
    draw();
  }

  // =====================================================================
  // 7. Energy landscape (Trefftz criterion)
  // =====================================================================
  function stability(root) {
    const { controls, body, readout } = frame(root, "Visualisation", "Critère énergétique : Π(ηw) = (U_b − λ W_m) η²");
    body.classList.add("single");
    const sL = slider(controls, "λ / λ<sub>cr</sub>", { min: 0, max: 1.6, step: 0.01, value: 0.6, fmt: (v) => fr(v, 2) });
    const p = panel(body);
    function draw() {
      const l = +sL.value, eta = d3.range(-1, 1.0001, 0.02);
      const col = l < 0.995 ? TEAL : l > 1.005 ? RED : "currentColor";
      put(p, Plot.plot({
        width: 560, height: 230, marginLeft: 44, x: { label: "amplitude η →" }, y: { domain: [-0.65, 1.05], label: "↑ Π", grid: true },
        marks: [Plot.ruleY([0]), Plot.line(eta.map((e) => [e, (1 - l) * e * e]), { stroke: col, strokeWidth: 3 }), Plot.dot([[0, 0]], { r: 7, fill: col, stroke: "white" })],
      }));
      readout.innerHTML = kv("δ²Π[0]", `${fr(2 * (1 - l), 2)} · U_b`) + kv("état plan", l < 0.995 ? "stable" : l > 1.005 ? "instable" : "critique (neutre)");
    }
    sL.addEventListener("input", draw); draw();
  }

  const REGISTRY = { buckling, navier, platebeam, diskritz, ritz1, rrconv, stability };

  function mountAll(scope) {
    (scope || document).querySelectorAll("[data-widget]").forEach((node) => {
      if (node.dataset.mounted) return;
      const fn = REGISTRY[node.dataset.widget];
      if (!fn) return;
      try { fn(node); node.dataset.mounted = "1"; } catch (e) { node.innerHTML = `<div class="pw-note">Erreur du widget : ${e.message}</div>`; console.error(e); }
    });
  }
  // "show / hide all solutions" toolbar on exercise pages
  function solToolbar() {
    document.querySelectorAll("[data-sol]").forEach((btn) => {
      btn.addEventListener("click", () => {
        const open = btn.dataset.sol === "open";
        document.querySelectorAll(".callout-tip .callout-collapse").forEach((el) => {
          if (window.bootstrap && bootstrap.Collapse) {
            const c = bootstrap.Collapse.getOrCreateInstance(el, { toggle: false });
            open ? c.show() : c.hide();
          } else el.classList.toggle("show", open);
          const head = el.parentElement && el.parentElement.querySelector(".callout-header");
          if (head) head.classList.toggle("collapsed", !open);
        });
      });
    });
  }
  window.PlatesWidgets = { mountAll, REGISTRY };
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", solToolbar); else solToolbar();
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", () => mountAll());
  else mountAll();
})();
