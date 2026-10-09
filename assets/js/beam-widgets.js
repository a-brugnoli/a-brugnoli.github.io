/* =====================================================================
   Structural dynamics — interactive widgets for the beam FE course
   Requires d3.min.js, plot.umd.min.js, beam-core.js (loaded before).
   Usage: <div class="pw" data-widget="hermite"></div>
   Widgets: kinematics, hermite, elemK, assembly, cantilever, modes, modalsum, frf
   ===================================================================== */
(function () {
  "use strict";
  const C = window.BeamCore;
  const PI = Math.PI;
  const TEAL = "#0b7a75", AMBER = "#d68a1c", RED = "#c0392b", INK = "#0e3b3f", BLUE = "#2b7bb9", PLUM = "#a63d6b", GREEN = "#6a8f2b";
  const SERIES = [TEAL, AMBER, BLUE, PLUM, GREEN, "#7a5bb8"];
  const EL_COL = [RED, BLUE, GREEN, AMBER, "#7a5bb8"];
  const SOFT = "var(--ab-soft, #eef6f5)";

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
    inp.fmtOut = (v) => { inp.value = v; out.innerHTML = fmt(+inp.value); };
    return inp;
  }
  function select(parent, label, options, value) {
    const lab = el("label");
    lab.innerHTML = `<span>${label}</span>`;
    const s = el("select");
    for (const [v, t] of options) {
      const o = el("option", { value: v }, t);
      if (String(v) === String(value)) o.selected = true;
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
  function button(parent, label) {
    const b = el("button", { type: "button", class: "pw-btn" }, label);
    parent.append(b);
    return b;
  }
  const kv = (k, v) => `<span class="kv">${k} = <b>${v}</b></span>`;
  function panel(body, cls = "") { const p = el("div", { class: "pw-panel " + cls }); body.append(p); return p; }
  function put(p, node) { p.innerHTML = ""; p.append(node); }
  const CUR = "currentColor";

  // animation loop that stops itself when the widget leaves the DOM / is hidden
  function animate(root, fn) {
    let raf = 0, t0 = performance.now();
    const loop = (now) => {
      if (!root.isConnected) return;
      if (root.offsetParent !== null) fn((now - t0) / 1000);
      raf = requestAnimationFrame(loop);
    };
    raf = requestAnimationFrame(loop);
    return () => cancelAnimationFrame(raf);
  }

  // cached modal computations (they are slow-ish for large Ne)
  const modalCache = new Map();
  function modalOf(Ne, bc) {
    const k = Ne + "|" + bc;
    if (!modalCache.has(k)) {
      const m = C.modal(Ne, bc, { L: 1, EI: 2e5, rhoS: 0.78 });
      // normalise shapes so that max |v| = 1 and is positive
      m.shapes = m.shapes.map((s) => {
        let mx = 0, sg = 1;
        for (let i = 0; i < s.length; i += 2) if (Math.abs(s[i]) > mx) { mx = Math.abs(s[i]); sg = Math.sign(s[i]); }
        // tip rotation can dominate for SS: also look at the interpolated curve
        const curve = C.sample(s, Ne, 1, 1, 6);
        let cm = 0, cs = 1;
        curve.v.forEach((v) => { if (Math.abs(v) > cm) { cm = Math.abs(v); cs = Math.sign(v); } });
        const out = new Float64Array(s.length);
        for (let i = 0; i < s.length; i++) out[i] = (cs / cm) * s[i];
        return out;
      });
      modalCache.set(k, m);
    }
    return modalCache.get(k);
  }

  // =====================================================================
  // 1. Kinematics (Euler–Bernoulli: sections stay perpendicular to the axis)
  // =====================================================================
  function kinematics(root) {
    const { controls, body, readout } = frame(root, "Simulation", "Cinématique d'Euler–Bernoulli : les sections restent perpendiculaires à l'axe");
    const sD = slider(controls, "flèche en bout δ/L", { min: 0, max: 0.45, step: 0.005, value: 0.25, fmt: (v) => fr(v, 2) });
    const sX = slider(controls, "section x₀/L", { min: 0, max: 1, step: 0.02, value: 0.6, fmt: (v) => fr(v, 2) });
    const mode = select(controls, "champ", [["exact", "rotation exacte"], ["lin", "petits angles"], ["both", "les deux"]], "exact");
    const cb = checkbox(controls, "colorer σₓₓ", true);
    const p = panel(body);
    body.classList.add("single");
    const W = 700, H = 330, s = 560, x0 = 60, y0 = 58, h = 0.13, Nx = 28, Ny = 6;

    function draw() {
      const d = +sD.value, xs = +sX.value, m = mode.value;
      const v = (x) => -d * (3 * x * x - x ** 3) / 2;
      const ph = (x) => -d * (6 * x - 3 * x * x) / 2;
      const kap = (x) => -3 * d * (1 - x);
      const P = (X, Y) => [x0 + s * X, y0 - s * Y];
      const exact = (x, y) => [x - y * Math.sin(ph(x)), v(x) + y * Math.cos(ph(x))];
      const lin = (x, y) => [x - y * ph(x), v(x) + y];
      const svg = d3.create("svg").attr("viewBox", `0 0 ${W} ${H}`).style("width", "100%").style("font", "13px Inter, system-ui, sans-serif");
      // wall
      svg.append("line").attr("x1", x0).attr("x2", x0).attr("y1", y0 - 60).attr("y2", y0 + 60).attr("stroke", CUR).attr("stroke-width", 4);
      for (let k = -60; k < 60; k += 10) svg.append("line").attr("x1", x0 - 12).attr("x2", x0).attr("y1", y0 + k + 10).attr("y2", y0 + k).attr("stroke", CUR).attr("opacity", 0.6);
      // undeformed ghost
      svg.append("rect").attr("x", x0).attr("y", y0 - s * h / 2).attr("width", s).attr("height", s * h).attr("fill", "none").attr("stroke", CUR).attr("stroke-dasharray", "4 4").attr("opacity", 0.35);
      const field = m === "lin" ? lin : exact;
      // cells
      const cells = svg.append("g");
      for (let i = 0; i < Nx; i++) for (let j = 0; j < Ny; j++) {
        const xa = i / Nx, xb = (i + 1) / Nx, ya = (j / Ny - 0.5) * h, yb = ((j + 1) / Ny - 0.5) * h;
        const pts = [[xa, ya], [xb, ya], [xb, yb], [xa, yb]].map(([x, y]) => P(...field(x, y)));
        const xm = (xa + xb) / 2, ym = (ya + yb) / 2;
        const sig = (ym * (1 - xm)) / (h / 2); // normalised σxx ∈ [-1,1]
        const fill = cb.checked ? d3.interpolateRdBu(0.5 - 0.46 * sig) : "rgba(11,122,117,.18)";
        cells.append("polygon").attr("points", pts.map((q) => q.join(",")).join(" ")).attr("fill", fill).attr("stroke", "rgba(0,0,0,.18)").attr("stroke-width", 0.6);
      }
      // linear boundary overlay
      if (m === "both") {
        const top = d3.range(0, Nx + 1).map((i) => P(...lin(i / Nx, h / 2)).join(",")).join(" ");
        const bot = d3.range(0, Nx + 1).map((i) => P(...lin(i / Nx, -h / 2)).join(",")).join(" ");
        [top, bot].forEach((pts) => svg.append("polyline").attr("points", pts).attr("fill", "none").attr("stroke", AMBER).attr("stroke-width", 2.2).attr("stroke-dasharray", "6 4"));
        const sl = [P(...lin(xs, -h / 2)), P(...lin(xs, h / 2))];
        svg.append("line").attr("x1", sl[0][0]).attr("y1", sl[0][1]).attr("x2", sl[1][0]).attr("y2", sl[1][1]).attr("stroke", AMBER).attr("stroke-width", 2.2);
      }
      // neutral axis
      const axis = d3.range(0, 101).map((i) => P(i / 100, v(i / 100)).join(",")).join(" ");
      svg.append("polyline").attr("points", axis).attr("fill", "none").attr("stroke", INK).attr("stroke-width", 2).attr("stroke-dasharray", "7 4");
      // highlighted section
      const f2 = m === "lin" ? lin : exact;
      const a = P(...f2(xs, -h / 2)), b = P(...f2(xs, h / 2)), c0 = P(...f2(xs, 0));
      svg.append("line").attr("x1", a[0]).attr("y1", a[1]).attr("x2", b[0]).attr("y2", b[1]).attr("stroke", RED).attr("stroke-width", 3.2);
      // right-angle marker between tangent and section (screen coordinates)
      const ang = ph(xs); // tangent angle (y up) ; in screen y is flipped
      const t = [Math.cos(ang), -Math.sin(ang)];           // tangent, screen coords
      const u = [t[0] * 13, t[1] * 13];
      const up = [-Math.sin(ang) * 13, -Math.cos(ang) * 13]; // section direction (towards +y), screen coords
      if (m !== "lin") {
        svg.append("polyline").attr("fill", "none").attr("stroke", RED).attr("stroke-width", 1.6)
          .attr("points", [[c0[0] + u[0], c0[1] + u[1]], [c0[0] + u[0] + up[0], c0[1] + u[1] + up[1]], [c0[0] + up[0], c0[1] + up[1]]].map((q) => q.join(",")).join(" "));
      }
      const tl = 70;
      svg.append("line").attr("x1", c0[0] - t[0] * tl).attr("y1", c0[1] - t[1] * tl).attr("x2", c0[0] + t[0] * tl).attr("y2", c0[1] + t[1] * tl).attr("stroke", TEAL).attr("stroke-width", 1.8).attr("opacity", 0.9);
      svg.append("text").attr("x", c0[0] + t[0] * tl + 6).attr("y", c0[1] + t[1] * tl + 4).attr("fill", TEAL).text("tangente");
      svg.append("text").attr("x", b[0] + 8).attr("y", b[1] - 6).attr("fill", RED).text("section");
      // axes
      svg.append("line").attr("x1", x0 + 6).attr("y1", H - 34).attr("x2", x0 + 56).attr("y2", H - 34).attr("stroke", CUR).attr("marker-end", "url(#a)");
      svg.append("text").attr("x", x0 + 62).attr("y", H - 30).attr("fill", CUR).text("x");
      svg.append("line").attr("x1", x0 + 6).attr("y1", H - 34).attr("x2", x0 + 6).attr("y2", H - 82).attr("stroke", CUR);
      svg.append("text").attr("x", x0 + 2).attr("y", H - 88).attr("fill", CUR).text("y");
      // colour legend
      const gid = "g" + Math.random().toString(36).slice(2, 7);
      const lg = svg.append("defs").append("linearGradient").attr("id", gid);
      d3.range(0, 1.001, 0.1).forEach((tt) => lg.append("stop").attr("offset", tt * 100 + "%").attr("stop-color", d3.interpolateRdBu(0.5 - 0.46 * (1 - 2 * tt))));
      if (cb.checked) {
        svg.append("rect").attr("x", W - 220).attr("y", H - 40).attr("width", 180).attr("height", 10).attr("rx", 5).attr("fill", `url(#${gid})`);
        svg.append("text").attr("x", W - 220).attr("y", H - 46).attr("fill", CUR).attr("font-size", 12).text("σₓₓ : traction");
        svg.append("text").attr("x", W - 40).attr("y", H - 46).attr("text-anchor", "end").attr("fill", CUR).attr("font-size", 12).text("compression");
      }
      put(p, svg.node());
      const yt = h / 2, E1 = exact(1, yt), L1 = lin(1, yt);
      readout.innerHTML = kv("φ(x₀)", fr(ph(xs) * 180 / PI, 1) + "°") + kv("v(x₀)/L", fr(v(xs), 3)) + kv("κ = v″", fr(kap(xs), 3) + " /L") +
        kv("ε<sub>xx</sub>(x₀, h/2)", fr(-yt * kap(xs), 4)) + kv("écart exact / linéaire (coin du bout)", fr(Math.hypot(E1[0] - L1[0], E1[1] - L1[1]), 4) + " L");
    }
    [sD, sX, mode, cb].forEach((c) => c.addEventListener("input", draw));
    mode.addEventListener("change", draw); cb.addEventListener("change", draw);
    draw();
  }

  // =====================================================================
  // 2. Hermite cubic shape functions
  // =====================================================================
  function hermite(root) {
    const { controls, body, readout } = frame(root, "Simulation", "Élément de poutre : interpolation d'Hermite  v(x)=N(x)·q<sub>e</sub>");
    const preset = select(controls, "préréglage", [
      ["free", "libre"], ["1", "N₁ : v₁ = 1"], ["2", "N₂ : φ₁ = 1"], ["3", "N₃ : v₂ = 1"], ["4", "N₄ : φ₂ = 1"], ["rt", "translation rigide"], ["rr", "rotation rigide"], ["bend", "flexion pure (courbure const.)"],
    ], "free");
    const s1 = slider(controls, "v₁", { min: -1, max: 1, step: 0.05, value: 0, fmt: (v) => fr(v, 2) });
    const s2 = slider(controls, "φ₁", { min: -2, max: 2, step: 0.05, value: 1, fmt: (v) => fr(v, 2) });
    const s3 = slider(controls, "v₂", { min: -1, max: 1, step: 0.05, value: 0.5, fmt: (v) => fr(v, 2) });
    const s4 = slider(controls, "φ₂", { min: -2, max: 2, step: 0.05, value: -0.5, fmt: (v) => fr(v, 2) });
    const pA = panel(body), pB = panel(body);
    const sl = [s1, s2, s3, s4];
    const PRE = { 1: [1, 0, 0, 0], 2: [0, 1, 0, 0], 3: [0, 0, 1, 0], 4: [0, 0, 0, 1], rt: [0.5, 0, 0.5, 0], rr: [0, 1, 1, 1], bend: [0, 0, 0.5, 1] };
    preset.addEventListener("change", () => { if (PRE[preset.value]) PRE[preset.value].forEach((v, i) => sl[i].fmtOut(v)); draw(); });
    const xs = d3.range(0, 1.0001, 0.01);
    const names = ["N₁", "N₂", "N₃", "N₄"];
    function draw() {
      const q = sl.map((s) => +s.value);
      const data = [];
      xs.forEach((x) => { const N = C.hermite(x, 1); N.forEach((y, k) => data.push({ x, y, k: names[k] })); });
      put(pA, Plot.plot({
        width: 430, height: 270, marginLeft: 42, grid: true,
        x: { label: "ξ = x / Lₑ →" }, y: { label: "↑ N_i", domain: [-0.25, 1.1] },
        color: { domain: names, range: SERIES.slice(0, 4), legend: true },
        marks: [Plot.ruleY([0]), Plot.line(data, { x: "x", y: "y", stroke: "k", strokeWidth: 2.4 })],
      }));
      const v = xs.map((x) => { const N = C.hermite(x, 1); return N[0] * q[0] + N[1] * q[1] + N[2] * q[2] + N[3] * q[3]; });
      const dx = 0.14;
      const tang = [[0, q[0], q[1]], [1, q[2], q[3]]].map(([x, y, m]) => ({ x1: x - dx, y1: y - m * dx, x2: x + dx, y2: y + m * dx }));
      put(pB, Plot.plot({
        width: 430, height: 270, marginLeft: 42, grid: true,
        x: { label: "ξ →", domain: [-0.2, 1.2] }, y: { label: "↑ v(ξ)", domain: [-2.2, 2.2] },
        marks: [Plot.ruleY([0]), Plot.line(xs.map((x, i) => [x, v[i]]), { stroke: TEAL, strokeWidth: 3.2 }),
          Plot.link(tang, { x1: "x1", y1: "y1", x2: "x2", y2: "y2", stroke: AMBER, strokeWidth: 2.2, strokeDasharray: "5 3" }),
          Plot.dot([[0, q[0]], [1, q[2]]], { r: 6, fill: INK, stroke: "white" }),
          Plot.text([[0, q[0]], [1, q[2]]], { text: ["nœud 1", "nœud 2"], dy: 18, fill: CUR })],
      }));
      const a0 = q[0], a1 = q[1], a2 = -3 * q[0] - 2 * q[1] + 3 * q[2] - q[3], a3 = 2 * q[0] + q[1] - 2 * q[2] + q[3];
      readout.innerHTML = kv("a₀,a₁", `${fr(a0, 2)} ; ${fr(a1, 2)}`) + kv("a₂", fr(a2, 2)) + kv("a₃", fr(a3, 2)) +
        kv("v″(0) ∝ M₁", fr(2 * a2, 2)) + kv("v″(1) ∝ M₂", fr(2 * a2 + 6 * a3, 2)) + kv("v‴ ∝ T (constant)", fr(6 * a3, 2));
    }
    sl.forEach((s) => s.addEventListener("input", () => { preset.value = "free"; draw(); }));
    draw();
  }

  // =====================================================================
  // 3. Element matrices: column j of K_e = nodal forces for a unit dof
  // =====================================================================
  function elemK(root) {
    const { controls, body, readout } = frame(root, "Simulation", "Matrice de raideur élémentaire : que représente une colonne de K<sub>e</sub> ?");
    const sel = select(controls, "déformée imposée", [
      ["0", "v₁ = 1"], ["1", "φ₁ = 1"], ["2", "v₂ = 1"], ["3", "φ₂ = 1"], ["rt", "translation rigide (1,0,1,0)"], ["rr", "rotation rigide (0,1,Lₑ,1)"], ["mix", "(v₁,φ₁,v₂,φ₂) = (0, 0.5, 0.3, −0.4)"],
    ], "0");
    const sL = slider(controls, "Lₑ", { min: 0.5, max: 2, step: 0.1, value: 1, fmt: (v) => fr(v, 1) });
    const pA = panel(body), pB = panel(body);
    const lab = ["v₁", "φ₁", "v₂", "φ₂"];
    function draw() {
      const Le = +sL.value, EI = 1;
      const ke = C.Ke(EI, Le), me = C.Me(1, Le);
      let q = [0, 0, 0, 0], hi = -1;
      if (/^\d$/.test(sel.value)) { q[+sel.value] = 1; hi = +sel.value; }
      else if (sel.value === "rt") q = [1, 0, 1, 0];
      else if (sel.value === "rr") q = [0, 1, Le, 1];
      else q = [0, 0.5, 0.3, -0.4];
      const f = ke.map((r) => r.reduce((s, kij, j) => s + kij * q[j], 0));
      // --- left: element drawing
      const W = 520, H = 330, x0 = 70, X1 = 450, ym = 105, yb = 250;
      const px = (xi) => x0 + xi * (X1 - x0);
      const svg = d3.create("svg").attr("viewBox", `0 0 ${W} ${H}`).style("width", "100%").style("font", "13px Inter, system-ui, sans-serif");
      const xiv = d3.range(0, 1.001, 0.02);
      const vv = xiv.map((xi) => { const N = C.hermite(xi, Le); return N[0] * q[0] + N[1] * q[1] + N[2] * q[2] + N[3] * q[3]; });
      const ys = 62 / Math.max(0.35, d3.max(vv, Math.abs));
      svg.append("line").attr("x1", px(0)).attr("x2", px(1)).attr("y1", ym).attr("y2", ym).attr("stroke", CUR).attr("stroke-dasharray", "4 4").attr("opacity", 0.4);
      svg.append("polyline").attr("points", xiv.map((xi, k) => `${px(xi)},${ym - ys * vv[k]}`).join(" ")).attr("fill", "none").attr("stroke", TEAL).attr("stroke-width", 5).attr("stroke-linecap", "round");
      [0, 1].forEach((k) => {
        svg.append("circle").attr("cx", px(k)).attr("cy", ym - ys * q[2 * k]).attr("r", 6).attr("fill", INK);
        svg.append("line").attr("x1", px(k)).attr("x2", px(k)).attr("y1", ym - ys * q[2 * k] + 8).attr("y2", yb - 70).attr("stroke", CUR).attr("opacity", 0.18).attr("stroke-dasharray", "2 4");
        svg.append("text").attr("x", px(k)).attr("y", ym - ys * q[2 * k] + 26).attr("text-anchor", "middle").attr("fill", CUR).attr("opacity", 0.8).text("nœud " + (k + 1));
      });
      svg.append("text").attr("x", 10).attr("y", 20).attr("fill", CUR).attr("font-weight", 600).text("q = (" + q.map((v) => fr(v, 2)).join(" ; ") + ")  →  f = K_e q");
      svg.append("text").attr("x", 10).attr("y", 190).attr("fill", CUR).attr("opacity", 0.7).attr("font-size", 12).text("efforts nodaux équivalents (échelle commune) :");
      // nodal forces (arrows) and moments (arcs) on a baseline
      const fm = Math.max(1e-9, ...f.map(Math.abs)), tol = 1e-9 * Math.max(1, fm);
      svg.append("line").attr("x1", x0 - 30).attr("x2", X1 + 30).attr("y1", yb).attr("y2", yb).attr("stroke", CUR).attr("opacity", 0.15);
      const arrow = (x, val, color) => {
        if (Math.abs(val) < tol) { svg.append("circle").attr("cx", x).attr("cy", yb).attr("r", 3).attr("fill", color); return; }
        const len = 10 + 56 * Math.abs(val) / fm, dir = val > 0 ? -1 : 1; // positive force points up
        svg.append("line").attr("x1", x).attr("y1", yb).attr("x2", x).attr("y2", yb + dir * (len - 4)).attr("stroke", color).attr("stroke-width", 3.4);
        svg.append("polygon").attr("points", `${x - 6},${yb + dir * (len - 11)} ${x + 6},${yb + dir * (len - 11)} ${x},${yb + dir * len}`).attr("fill", color);
      };
      const arc = (x, val, color) => {
        if (Math.abs(val) < tol) return;
        const r = 14 + 14 * Math.abs(val) / fm, ccw = val > 0; // positive moment = counter-clockwise
        const ts = ccw ? -0.75 * PI : 0.75 * PI, te = -ts;
        const P = (a) => [x + r * Math.cos(a), yb - r * Math.sin(a)];
        const S = P(ts), E = P(te);
        svg.append("path").attr("d", `M${S[0]},${S[1]} A${r},${r} 0 1 ${ccw ? 0 : 1} ${E[0]},${E[1]}`).attr("fill", "none").attr("stroke", color).attr("stroke-width", 3);
        const sg = ccw ? 1 : -1, tg = [-Math.sin(te) * sg, -Math.cos(te) * sg], pp = [-tg[1], tg[0]];
        svg.append("polygon").attr("fill", color).attr("points", [[E[0] + 9 * tg[0], E[1] + 9 * tg[1]], [E[0] + 6 * pp[0], E[1] + 6 * pp[1]], [E[0] - 6 * pp[0], E[1] - 6 * pp[1]]].map((q2) => q2.join(",")).join(" "));
      };
      arrow(px(0), f[0], RED); arc(px(0) + 62, f[1], BLUE);
      arrow(px(1), f[2], RED); arc(px(1) - 62, f[3], BLUE);
      const lbl = (x, a, b) => {
        svg.append("text").attr("x", x).attr("y", H - 22).attr("text-anchor", "middle").attr("fill", RED).attr("font-weight", 600).text(a);
        svg.append("text").attr("x", x).attr("y", H - 6).attr("text-anchor", "middle").attr("fill", BLUE).attr("font-weight", 600).text(b);
      };
      lbl(px(0) + 20, `F₁ = ${fr(f[0], 2)}`, `M₁ = ${fr(f[1], 2)}`); lbl(px(1) - 20, `F₂ = ${fr(f[2], 2)}`, `M₂ = ${fr(f[3], 2)}`);
      put(pA, svg.node());
      // --- right: matrices
      const tbl = (A, title, hiCol, fmtc) => {
        let h = `<div style="font-weight:600;margin:.1em 0">${title}</div><table style="border-collapse:collapse;font-variant-numeric:tabular-nums;margin-bottom:.4em"><thead><tr><th></th>${lab.map((l) => `<th style="padding:.1em .45em;color:${AMBER}">${l}</th>`).join("")}</tr></thead><tbody>`;
        A.forEach((r, i) => { h += `<tr><th style="padding:.1em .45em;color:${AMBER}">${["F₁", "M₁", "F₂", "M₂"][i]}</th>` + r.map((v, j) => `<td style="text-align:right;padding:.1em .5em;${j === hiCol ? "background:rgba(11,122,117,.18);font-weight:700;" : ""}">${fmtc(v)}</td>`).join("") + "</tr>"; });
        return h + "</tbody></table>";
      };
      pB.innerHTML = tbl(ke, "K<sub>e</sub>  (EI = 1)", hi, (v) => fr(v, 2)) + tbl(me.map((r) => r.map((v) => v * 420 / Le)), "M<sub>e</sub> · 420/(ρS Lₑ)", hi, (v) => fr(v, 2));
      const nrm = Math.hypot(...f);
      readout.innerHTML = kv("‖K<sub>e</sub> q‖", fr(nrm, 3)) + (nrm < 1e-9 ? kv("conclusion", "mouvement rigide : aucun effort ⇒ K<sub>e</sub> singulière") : "") +
        kv("ΣF", fr(f[0] + f[2], 3)) + kv("équilibre des moments", fr(f[1] + f[3] + f[2] * Le, 3));
    }
    sel.addEventListener("change", draw); sL.addEventListener("input", draw);
    draw();
  }

  // =====================================================================
  // 4. Assembly of K and f
  // =====================================================================
  function assembly(root) {
    const { controls, body, readout } = frame(root, "Simulation", "Assemblage de la matrice de raideur et du vecteur force");
    const sN = slider(controls, "nombre d'éléments", { min: 2, max: 4, step: 1, value: 3, fmt: (v) => v });
    const sK = slider(controls, "éléments assemblés", { min: 0, max: 3, step: 1, value: 3, fmt: (v) => v });
    const p = panel(body);
    body.classList.add("single");
    function draw() {
      const Ne = +sN.value; sK.max = Ne; if (+sK.value > Ne) sK.fmtOut(Ne);
      const k = +sK.value, n = 2 * (Ne + 1), cs = Ne >= 4 ? 40 : 46;
      const Ke1 = [[12, 6, -12, 6], [6, 4, -6, 2], [-12, -6, 12, -6], [6, 2, -6, 4]];
      const fe1 = [6, 1, 6, -1];
      const mx = 56, my = 44, W = mx + n * cs + 120 + 70, H = my + n * cs + 28;
      const K = Array.from({ length: n }, () => new Array(n).fill(0)), f = new Array(n).fill(0), cnt = Array.from({ length: n }, () => new Array(n).fill(0));
      const svg = d3.create("svg").attr("viewBox", `0 0 ${W} ${H}`).style("width", "100%").style("font", "13px Inter, system-ui, sans-serif");
      for (let e = 0; e < k; e++) for (let a = 0; a < 4; a++) { f[2 * e + a] += fe1[a]; for (let b = 0; b < 4; b++) { K[2 * e + a][2 * e + b] += Ke1[a][b]; cnt[2 * e + a][2 * e + b]++; } }
      // element backgrounds
      for (let e = 0; e < k; e++) {
        svg.append("rect").attr("x", mx + 2 * e * cs).attr("y", my + 2 * e * cs).attr("width", 4 * cs).attr("height", 4 * cs).attr("fill", EL_COL[e]).attr("opacity", e === k - 1 ? 0.3 : 0.18);
      }
      // grid
      for (let i = 0; i <= n; i++) {
        svg.append("line").attr("x1", mx).attr("x2", mx + n * cs).attr("y1", my + i * cs).attr("y2", my + i * cs).attr("stroke", CUR).attr("opacity", 0.14);
        svg.append("line").attr("y1", my).attr("y2", my + n * cs).attr("x1", mx + i * cs).attr("x2", mx + i * cs).attr("stroke", CUR).attr("opacity", 0.14);
      }
      svg.append("rect").attr("x", mx).attr("y", my).attr("width", n * cs).attr("height", n * cs).attr("fill", "none").attr("stroke", CUR).attr("stroke-width", 1.4);
      for (let e = 0; e < k; e++) svg.append("rect").attr("x", mx + 2 * e * cs).attr("y", my + 2 * e * cs).attr("width", 4 * cs).attr("height", 4 * cs).attr("fill", "none").attr("stroke", EL_COL[e]).attr("stroke-width", 2.4);
      // labels
      for (let i = 0; i < n; i++) {
        const t = (i % 2 ? "φ" : "v") + "₁₂₃₄₅"[Math.floor(i / 2)];
        svg.append("text").attr("x", mx + i * cs + cs / 2).attr("y", my - 10).attr("text-anchor", "middle").attr("fill", AMBER).text(t);
        svg.append("text").attr("x", mx - 10).attr("y", my + i * cs + cs / 2 + 5).attr("text-anchor", "end").attr("fill", AMBER).text(t);
      }
      // numbers
      for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) if (cnt[i][j]) {
        svg.append("text").attr("x", mx + j * cs + cs / 2).attr("y", my + i * cs + cs / 2 + 5).attr("text-anchor", "middle").attr("fill", CUR)
          .attr("font-weight", cnt[i][j] > 1 ? 800 : 500).attr("font-size", cnt[i][j] > 1 ? 15 : 13).text(K[i][j]);
        if (cnt[i][j] > 1) svg.append("rect").attr("x", mx + j * cs + 3).attr("y", my + i * cs + 3).attr("width", cs - 6).attr("height", cs - 6).attr("fill", "none").attr("stroke", RED).attr("stroke-width", 1.6).attr("rx", 5);
      }
      // f vector
      const fx = mx + n * cs + 36;
      svg.append("text").attr("x", fx + 22).attr("y", my - 10).attr("text-anchor", "middle").attr("fill", CUR).attr("font-weight", 600).text("f");
      for (let i = 0; i < n; i++) {
        svg.append("rect").attr("x", fx).attr("y", my + i * cs).attr("width", 44).attr("height", cs).attr("fill", "none").attr("stroke", CUR).attr("opacity", 0.25);
        const e = Math.floor(i / 2), contrib = [];
        for (let ee = Math.max(0, e - 1); ee <= Math.min(k - 1, e); ee++) contrib.push(ee);
        if (f[i] !== 0 || contrib.length) svg.append("text").attr("x", fx + 22).attr("y", my + i * cs + cs / 2 + 5).attr("text-anchor", "middle").attr("fill", CUR).attr("font-weight", 600).text(f[i]);
      }
      svg.append("text").attr("x", fx).attr("y", my + n * cs + 20).attr("fill", CUR).attr("opacity", 0.7).attr("font-size", 12).text("f en unités −pLₑ/12");
      svg.append("text").attr("x", mx).attr("y", my + n * cs + 20).attr("fill", CUR).attr("opacity", 0.7).attr("font-size", 12).text("K en unités EI/Lₑ³ (Lₑ = 1)");
      put(p, svg.node());
      readout.innerHTML = kv("ddl globaux", n) + kv("élément ajouté", k || "—") + kv("recouvrements", `${Math.max(0, k - 1)} nœud${k > 2 ? "s" : ""} partagé${k > 2 ? "s" : ""} (cases cerclées : somme)`) + kv("bande de K", "6 (largeur limitée)");
    }
    sN.addEventListener("input", draw); sK.addEventListener("input", draw);
    draw();
  }

  // =====================================================================
  // 5. Static cantilever FE solver + convergence
  // =====================================================================
  function cantilever(root) {
    const show = (root.dataset.show || "defl,moment,conv").split(",");
    const { controls, body, readout } = frame(root, "Simulation", "Poutre encastrée : éléments finis de Hermite vs solution exacte");
    const sN = slider(controls, "éléments N", { min: 1, max: 32, step: 1, value: +(root.dataset.n || 2), fmt: (v) => v });
    const sP = slider(controls, "charge p [N/m]", { min: 0, max: 200, step: 5, value: 100, fmt: (v) => v });
    const sF = slider(controls, "force F [N]", { min: 0, max: 40, step: 1, value: 20, fmt: (v) => v });
    const cN = checkbox(controls, "nœuds", true);
    const panels = {};
    show.forEach((k) => (panels[k] = panel(body)));
    const L = 1, EI = 2e5;
    const cache = {};
    function conv(p, F) {
      const key = p + "|" + F;
      if (cache[key]) return cache[key];
      const out = [];
      for (let Ne = 1; Ne <= 32; Ne++) out.push([Ne, errOf(Ne, p, F).rel]);
      return (cache[key] = out);
    }
    function errOf(Ne, p, F) {
      const kF = Math.max(1, Math.round(Ne / 2)), a = kF * L / Ne;
      const r = C.staticCantilever(Ne, { L, EI, p, F, kF });
      const s = C.sample(r.q, Ne, L, EI, 10);
      let e = 0, mx = 1e-30;
      s.x.forEach((x, i) => { const ex = C.exactCantilever(x, { L, EI, p, F, a }).v; e = Math.max(e, Math.abs(s.v[i] - ex)); mx = Math.max(mx, Math.abs(ex)); });
      return { rel: e / mx, r, s, a, kF };
    }
    function draw() {
      const Ne = +sN.value, p = +sP.value, F = +sF.value;
      const { rel, r, s, a, kF } = errOf(Ne, p, F);
      const ex = d3.range(0, 1.0001, 0.01).map((x) => ({ x, ...C.exactCantilever(x, { L, EI, p, F, a }) }));
      const mm = 1e3; // display in mm and N·m
      if (panels.defl) put(panels.defl, Plot.plot({
        width: 440, height: 250, marginLeft: 52, grid: true, x: { label: "x / L →" }, y: { label: "↑ flèche v [mm]" },
        marks: [Plot.ruleY([0]),
          Plot.line(ex, { x: "x", y: (d) => d.v * mm, stroke: CUR, strokeDasharray: "6 4", strokeWidth: 2, opacity: 0.8 }),
          Plot.line(s.x.map((x, i) => [x, s.v[i] * mm]), { stroke: TEAL, strokeWidth: 3.2 }),
          cN.checked ? Plot.dot(d3.range(0, Ne + 1).map((k) => [k * L / Ne, r.q[2 * k] * mm]), { r: 5, fill: AMBER, stroke: "white" }) : null,
          F > 0 ? Plot.dot([[a, r.q[2 * kF] * mm]], { r: 8, symbol: "diamond", fill: RED, stroke: "white" }) : null],
      }));
      if (panels.moment) {
        const data = s.x.map((x, i) => ({ x, M: s.M[i], e: s.elem[i] }));
        put(panels.moment, Plot.plot({
          width: 440, height: 250, marginLeft: 52, grid: true, x: { label: "x / L →" }, y: { label: "↑ moment M = EI v″ [N·m]" },
          marks: [Plot.ruleY([0]),
            Plot.line(ex, { x: "x", y: "M", stroke: CUR, strokeDasharray: "6 4", strokeWidth: 2, opacity: 0.8 }),
            Plot.line(data, { x: "x", y: "M", z: "e", stroke: AMBER, strokeWidth: 3 })],
        }));
      }
      if (panels.conv) {
        const cv = conv(p, F).map(([n, e]) => ({ n, e: Math.max(e, 1e-14) }));
        const ref = [1, 32].map((n) => ({ n, e: Math.max(cv[0].e, 1e-14) / n ** 4 }));
        put(panels.conv, Plot.plot({
          width: 440, height: 250, marginLeft: 56, grid: true,
          x: { type: "log", label: "nombre d'éléments N →", ticks: [1, 2, 4, 8, 16, 32], tickFormat: (d) => d },
          y: { type: "log", label: "↑ erreur relative max sur v", domain: [1e-14, 1] },
          marks: [Plot.line(ref, { x: "n", y: "e", stroke: CUR, strokeDasharray: "5 4", opacity: 0.6 }),
            Plot.text([ref[1]], { x: "n", y: "e", text: ["pente −4"], dy: -9, dx: -22, fill: CUR, opacity: 0.7 }),
            Plot.line(cv, { x: "n", y: "e", stroke: TEAL, strokeWidth: 2.4 }),
            Plot.dot(cv, { x: "n", y: "e", r: 3.5, fill: TEAL }),
            Plot.dot(cv.filter((d) => d.n === Ne), { x: "n", y: "e", r: 8, fill: AMBER, stroke: "white" })],
        }));
      }
      const exTip = C.exactCantilever(L, { L, EI, p, F, a }).v;
      const R = p * L + F, M0 = p * L * L / 2 + F * a;
      readout.innerHTML = kv("v(L) EF", fr(r.q[2 * Ne] * mm, 4) + " mm") + kv("v(L) exact", fr(exTip * mm, 4) + " mm") +
        kv("erreur max", rel < 1e-12 ? "≈ 0 (précision machine)" : sci(rel, 1)) +
        kv("réaction F₁", `${fr(r.reaction[0], 2)} N (exact ${fr(R, 2)})`) + kv("moment M₁", `${fr(r.reaction[1], 2)} N·m (exact ${fr(M0, 2)})`) +
        kv("EI", "2·10⁵ Pa·m⁴");
    }
    [sN, sP, sF].forEach((s) => s.addEventListener("input", draw)); cN.addEventListener("change", draw);
    draw();
  }

  // =====================================================================
  // 6. Modal analysis with animation
  // =====================================================================
  function drawBeamSVG(svg, W, H, s, x0, yc, bc) {
    if (bc === "cantilever" || bc === "cc") {
      svg.append("line").attr("x1", x0).attr("x2", x0).attr("y1", yc - 40).attr("y2", yc + 40).attr("stroke", CUR).attr("stroke-width", 4);
      for (let k = -40; k < 40; k += 10) svg.append("line").attr("x1", x0 - 11).attr("x2", x0).attr("y1", yc + k + 10).attr("y2", yc + k).attr("stroke", CUR).attr("opacity", 0.55);
    }
    if (bc === "ss") svg.append("polygon").attr("points", `${x0},${yc} ${x0 - 12},${yc + 22} ${x0 + 12},${yc + 22}`).attr("fill", "none").attr("stroke", CUR).attr("stroke-width", 2.2);
    if (bc === "cc" || bc === "ss") {
      const xe = x0 + s;
      if (bc === "cc") { svg.append("line").attr("x1", xe).attr("x2", xe).attr("y1", yc - 40).attr("y2", yc + 40).attr("stroke", CUR).attr("stroke-width", 4); for (let k = -40; k < 40; k += 10) svg.append("line").attr("x1", xe + 11).attr("x2", xe).attr("y1", yc + k + 10).attr("y2", yc + k).attr("stroke", CUR).attr("opacity", 0.55); }
      else svg.append("polygon").attr("points", `${xe},${yc} ${xe - 12},${yc + 22} ${xe + 12},${yc + 22}`).attr("fill", "none").attr("stroke", CUR).attr("stroke-width", 2.2);
    }
  }

  function modes(root) {
    const { controls, body, readout } = frame(root, "Simulation", "Analyse modale : modes propres et fréquences (EF vs exact)");
    const bcS = select(controls, "conditions aux limites", [["cantilever", "encastrée–libre"], ["ss", "appuyée–appuyée"], ["cc", "encastrée–encastrée"]], root.dataset.bc || "cantilever");
    const sN = slider(controls, "éléments N", { min: 2, max: 30, step: 1, value: +(root.dataset.n || 8), fmt: (v) => v });
    const sM = slider(controls, "mode i", { min: 1, max: 8, step: 1, value: 1, fmt: (v) => v });
    const sA = slider(controls, "amplitude", { min: 0.03, max: 0.2, step: 0.01, value: 0.14, fmt: (v) => fr(v, 2) });
    const play = checkbox(controls, "animer", true);
    const pA = panel(body), pB = panel(body);
    const W = 520, H = 250, s = 440, x0 = 40, yc = 130;
    const svg = d3.create("svg").attr("viewBox", `0 0 ${W} ${H}`).style("width", "100%").style("font", "13px Inter, system-ui, sans-serif");
    put(pA, svg.node());
    let gDyn = null, cur = null;
    function setup() {
      const bc = bcS.value, Ne = +sN.value;
      const m = modalOf(Ne, bc);
      sM.max = Math.min(8, m.omega.length);
      if (+sM.value > +sM.max) sM.fmtOut(sM.max);
      const i = +sM.value;
      const shape = m.shapes[i - 1];
      const cv = C.sample(shape, Ne, 1, 1, 8);
      cur = { cv, omega: m.omega[i - 1], Ne, bc, i, m };
      svg.selectAll("*").remove();
      drawBeamSVG(svg, W, H, s, x0, yc, bc);
      svg.append("line").attr("x1", x0).attr("x2", x0 + s).attr("y1", yc).attr("y2", yc).attr("stroke", CUR).attr("stroke-dasharray", "3 5").attr("opacity", 0.35);
      gDyn = {
        env1: svg.append("polyline").attr("fill", "none").attr("stroke", AMBER).attr("stroke-dasharray", "4 4").attr("opacity", 0.55),
        env2: svg.append("polyline").attr("fill", "none").attr("stroke", AMBER).attr("stroke-dasharray", "4 4").attr("opacity", 0.55),
        line: svg.append("polyline").attr("fill", "none").attr("stroke", TEAL).attr("stroke-width", 7).attr("stroke-linecap", "round").attr("stroke-linejoin", "round"),
        dots: svg.append("g"),
        txt: svg.append("text").attr("x", W - 6).attr("y", 20).attr("text-anchor", "end").attr("fill", CUR).attr("font-weight", 600),
      };
      gDyn.dots.selectAll("circle").data(d3.range(0, Ne + 1)).join("circle").attr("r", 4.5).attr("fill", AMBER).attr("stroke", "white");
      render(0);
      // ratio plot
      const nm = Math.min(m.omega.length, 24), data = [];
      for (let k = 1; k <= nm; k++) data.push({ k, r: m.omega[k - 1] / C.omegaExact(bc, k) });
      put(pB, Plot.plot({
        width: 430, height: 250, marginLeft: 50, grid: true,
        x: { label: "numéro du mode i →", ticks: 8 }, y: { type: "log", label: "↑ ω_EF / ω_exact", domain: [0.995, Math.max(1.05, d3.max(data, (d) => d.r) * 1.1)], tickFormat: (d) => fr(d, 2) },
        marks: [Plot.ruleY([1], { stroke: CUR, strokeDasharray: "4 4" }), Plot.line(data, { x: "k", y: "r", stroke: TEAL, strokeWidth: 2 }),
          Plot.dot(data, { x: "k", y: "r", r: 4, fill: TEAL }), Plot.dot(data.filter((d) => d.k === i), { x: "k", y: "r", r: 9, fill: AMBER, stroke: "white" })],
      }));
      const fe = cur.omega / (2 * PI), fx = C.omegaExact(bc, i) / (2 * PI);
      readout.innerHTML = kv(`f<sub>${i}</sub> EF`, fr(fe, 2) + " Hz") + kv(`f<sub>${i}</sub> exact`, fr(fx, 2) + " Hz") + kv("écart", fr(100 * (fe / fx - 1), 3) + " %") +
        kv("ddl libres", m.omega.length) + kv("EI, ρS, L", "2·10⁵, 0,78, 1");
    }
    function render(t) {
      if (!cur) return;
      const { cv, i } = cur;
      const A = +sA.value, c = play.checked ? Math.cos(2 * PI * t / 3.2) : 1;
      const mk = (a) => cv.x.map((x, k) => `${x0 + s * x},${yc - s * A * a * cv.v[k]}`).join(" ");
      gDyn.line.attr("points", mk(c));
      gDyn.env1.attr("points", mk(1)); gDyn.env2.attr("points", mk(-1));
      gDyn.dots.selectAll("circle").attr("cx", (_, k) => x0 + s * k / cur.Ne).attr("cy", (_, k) => yc - s * A * c * cv.v[k * 8]);
      gDyn.txt.text(`mode ${i}  ·  ω = ${fr(cur.omega, 1)} rad/s`);
    }
    [bcS, sN, sM].forEach((c) => c.addEventListener("input", setup)); bcS.addEventListener("change", setup);
    sA.addEventListener("input", () => render(performance.now() / 1000));
    animate(root, (t) => { if (play.checked) render(t); });
    setup();
  }

  // =====================================================================
  // 7. Modal superposition — mode invariance
  // =====================================================================
  function modalsum(root) {
    const { controls, body, readout } = frame(root, "Simulation", "Invariance des modes : q⁰ = Σ c<sub>i</sub> φ<sub>i</sub> ⇒ q(t) = Σ c<sub>i</sub> φ<sub>i</sub> cos ω<sub>i</sub>t");
    const c1 = slider(controls, "c₁", { min: -1, max: 1, step: 0.1, value: 1, fmt: (v) => fr(v, 1) });
    const c2 = slider(controls, "c₂", { min: -1, max: 1, step: 0.1, value: 0, fmt: (v) => fr(v, 1) });
    const c3 = slider(controls, "c₃", { min: -1, max: 1, step: 0.1, value: 0, fmt: (v) => fr(v, 1) });
    const sV = slider(controls, "vitesse", { min: 0.1, max: 1.5, step: 0.05, value: 0.4, fmt: (v) => fr(v, 2) });
    const play = checkbox(controls, "animer", true);
    const pA = panel(body), pB = panel(body);
    const Ne = 10, m = modalOf(Ne, "cantilever");
    const cvs = [0, 1, 2].map((i) => C.sample(m.shapes[i], Ne, 1, 1, 6));
    const w = [0, 1, 2].map((i) => m.omega[i] / m.omega[0]);
    const cs = [c1, c2, c3];
    const W = 520, H = 250, S = 440, x0 = 40, yc = 130;
    const svg = d3.create("svg").attr("viewBox", `0 0 ${W} ${H}`).style("width", "100%").style("font", "13px Inter, system-ui, sans-serif");
    drawBeamSVG(svg, W, H, S, x0, yc, "cantilever");
    svg.append("line").attr("x1", x0).attr("x2", x0 + S).attr("y1", yc).attr("y2", yc).attr("stroke", CUR).attr("stroke-dasharray", "3 5").attr("opacity", 0.35);
    const ghosts = d3.range(7).map(() => svg.append("polyline").attr("fill", "none").attr("stroke", TEAL).attr("stroke-linecap", "round"));
    const main = svg.append("polyline").attr("fill", "none").attr("stroke", TEAL).attr("stroke-width", 7).attr("stroke-linecap", "round").attr("stroke-linejoin", "round");
    put(pA, svg.node());
    // modal-coordinate chart with cursor
    const CW = 430, CH = 250, cm = { l: 40, r: 10, t: 12, b: 34 };
    const csvg = d3.create("svg").attr("viewBox", `0 0 ${CW} ${CH}`).style("width", "100%").style("font", "12px Inter, system-ui, sans-serif");
    const xT = d3.scaleLinear([0, 2], [cm.l, CW - cm.r]), yT = d3.scaleLinear([-1.1, 1.1], [CH - cm.b, cm.t]);
    csvg.append("g").attr("transform", `translate(0,${yT(0)})`).call(d3.axisBottom(xT).ticks(4).tickSizeOuter(0)).call((g) => g.selectAll("text").attr("fill", CUR)).call((g) => g.selectAll("path,line").attr("stroke", CUR).attr("opacity", 0.5));
    csvg.append("g").attr("transform", `translate(${cm.l},0)`).call(d3.axisLeft(yT).ticks(5).tickSizeOuter(0)).call((g) => g.selectAll("text").attr("fill", CUR)).call((g) => g.selectAll("path,line").attr("stroke", CUR).attr("opacity", 0.5));
    csvg.append("text").attr("x", CW - cm.r).attr("y", CH - 4).attr("text-anchor", "end").attr("fill", CUR).text("t / T₁ →");
    csvg.append("text").attr("x", cm.l + 6).attr("y", CH - 4).attr("fill", CUR).text("coordonnées modales η_i(t)");
    const paths = [0, 1, 2].map((i) => csvg.append("path").attr("fill", "none").attr("stroke", SERIES[i]).attr("stroke-width", 2));
    const cursor = csvg.append("line").attr("y1", cm.t).attr("y2", CH - cm.b).attr("stroke", RED).attr("stroke-width", 1.6);
    put(pB, csvg.node());
    const T1 = 1;
    function drawStatic() {
      const ts = d3.range(0, 2.0005, 0.002);
      paths.forEach((p, i) => p.attr("d", d3.line().x((t) => xT(t)).y((t) => yT(+cs[i].value * Math.cos(2 * PI * w[i] * t / T1)))(ts)).attr("opacity", Math.abs(+cs[i].value) < 1e-9 ? 0.15 : 1));
      const act = [0, 1, 2].filter((i) => Math.abs(+cs[i].value) > 1e-9).map((i) => i + 1);
      readout.innerHTML = kv("modes actifs", act.length ? act.join(", ") : "aucun") + kv("ω₂/ω₁", fr(w[1], 2)) + kv("ω₃/ω₁", fr(w[2], 2)) +
        kv("mouvement", act.length === 1 ? "synchrone : forme fixe (sous-espace invariant)" : act.length === 0 ? "repos" : "somme de modes : forme variable");
    }
    function frame_(tau) { // tau in units of T1
      const q = (t) => cvs[0].x.map((x, k) => cs.reduce((sum, c, i) => sum + +c.value * cvs[i].v[k] * Math.cos(2 * PI * w[i] * t), 0));
      const A = 0.14 * S;
      const pts = (t) => { const v = q(t); return cvs[0].x.map((x, k) => `${x0 + S * x},${yc - A * v[k]}`).join(" "); };
      ghosts.forEach((g, j) => g.attr("points", pts(tau - (j + 1) * 0.012)).attr("stroke-width", 7).attr("opacity", 0.22 * (1 - j / 8)));
      main.attr("points", pts(tau));
      cursor.attr("x1", xT(tau % 2)).attr("x2", xT(tau % 2));
    }
    [c1, c2, c3].forEach((c) => c.addEventListener("input", () => { drawStatic(); frame_(tauNow); }));
    let tauNow = 0;
    animate(root, (t) => { if (play.checked) { tauNow = t * +sV.value * 0.31; frame_(tauNow); } });
    drawStatic(); frame_(0);
  }

  // =====================================================================
  // 8. Frequency response & modal truncation
  // =====================================================================
  function frf(root) {
    const { controls, body, readout } = frame(root, "Simulation", "Réduction modale : réponse en fréquence en bout de poutre");
    const sM = slider(controls, "modes conservés m", { min: 1, max: 40, step: 1, value: 3, fmt: (v) => v });
    const sLo = slider(controls, "f<sub>min</sub> de la bande", { min: 1, max: 4.2, step: 0.02, value: 1.6, fmt: (v) => fr(Math.pow(10, v), 0) + " Hz" });
    const sHi = slider(controls, "f<sub>max</sub> de la bande", { min: 1, max: 4.2, step: 0.02, value: 2.9, fmt: (v) => fr(Math.pow(10, v), 0) + " Hz" });
    const p = panel(body);
    body.classList.add("single");
    const Ne = 20, L = 1, EI = 2e5, rhoS = 0.78, zeta = 0.005;
    const raw = C.modal(Ne, "cantilever", { L, EI, rhoS }); // un-normalised (M-orthonormal in non-dim units)
    const psi = raw.shapes.map((s) => s[2 * Ne] / Math.sqrt(rhoS * L)); // tip value of the physically mass-normalised shape
    const om = raw.omega, n = om.length;
    const fs = d3.range(1, 4.3, 0.01).map((e) => Math.pow(10, e));
    const H = (f, mm) => { const w = 2 * PI * f; let re = 0, im = 0; for (let i = 0; i < mm; i++) { const a = om[i] * om[i] - w * w, b = 2 * zeta * om[i] * w, d = a * a + b * b, g = psi[i] * psi[i]; re += g * a / d; im -= g * b / d; } return Math.hypot(re, im); };
    const full = fs.map((f) => ({ f, h: H(f, n) }));
    const slideMode = !!root.closest(".reveal"), W = slideMode ? 1060 : 760, Hh = slideMode ? 350 : 330;
    function draw() {
      const mm = +sM.value, lo = Math.pow(10, Math.min(+sLo.value, +sHi.value)), hi = Math.pow(10, Math.max(+sLo.value, +sHi.value));
      const red = fs.map((f) => ({ f, h: H(f, mm) }));
      let err = 0; fs.forEach((f, k) => { if (f >= lo && f <= hi) err = Math.max(err, Math.abs(red[k].h - full[k].h) / full[k].h); });
      const fm = om.slice(0, 8).map((w) => w / (2 * PI));
      put(p, Plot.plot({
        width: W, height: Hh, marginLeft: 62, marginBottom: 40, marginRight: 24, grid: true,
        x: { type: "log", label: "fréquence f [Hz] →", domain: [10, 2e4] }, y: { type: "log", label: "↑ |H(f)| déplacement / force [m/N]", domain: [1e-8, 1e-2] },
        marks: [Plot.rectX([[lo, hi]], { x1: (d) => d[0], x2: (d) => d[1], fill: AMBER, fillOpacity: 0.16 }),
          Plot.ruleX(fm.filter((f) => f < 2e4), { stroke: CUR, strokeOpacity: 0.15 }),
          Plot.line(full, { x: "f", y: "h", stroke: CUR, strokeWidth: 2.4, opacity: 0.85 }),
          Plot.line(red, { x: "f", y: "h", stroke: TEAL, strokeWidth: 3, strokeDasharray: mm >= n ? null : "7 4" }),
          Plot.ruleX(om.slice(0, mm).map((w) => w / (2 * PI)).filter((f) => f < 2e4), { stroke: TEAL, strokeOpacity: 0.55, strokeDasharray: "2 3", y1: 1e-8, y2: 4e-8 }),
          Plot.text([[lo, 6e-3]], { text: ["bande de la sollicitation"], textAnchor: "start", dx: 6, fill: AMBER, fontWeight: 600 })],
      }));
      const inBand = om.map((w) => w / (2 * PI)).filter((f) => f >= lo && f <= hi).length;
      const staticErr = Math.abs(red[0].h - full[0].h) / full[0].h;
      readout.innerHTML = kv("taille du système", `${n} ddl → ${mm} modes`) + kv("modes dans la bande", inBand) + kv("erreur max dans la bande", err < 1e-3 ? "< 0,1 %" : fr(100 * err, 1) + " %") +
        kv("compliance statique", sci(full[0].h, 2) + " m/N (exact L³/3EI = " + sci(L ** 3 / (3 * EI), 2) + ")") + kv("ω₁", fr(om[0], 0) + " rad/s");
    }
    [sM, sLo, sHi].forEach((s) => s.addEventListener("input", draw));
    draw();
  }

  const REGISTRY = { kinematics, hermite, elemK, assembly, cantilever, modes, modalsum, frf };
  function mountAll(scope) {
    (scope || document).querySelectorAll("[data-widget]").forEach((node) => {
      if (node.dataset.mounted) return;
      const fn = REGISTRY[node.dataset.widget];
      if (!fn) return;
      try { fn(node); node.dataset.mounted = "1"; } catch (e) { node.innerHTML = `<div class="pw-note">Erreur du widget : ${e.message}</div>`; console.error(e); }
    });
  }
  window.BeamWidgets = { mountAll, REGISTRY };
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", () => mountAll());
  else mountAll();
})();
