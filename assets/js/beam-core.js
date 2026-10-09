/* =====================================================================
   Euler–Bernoulli beam finite elements — numerical core (no DOM).
   Hermite cubic elements, assembly, static solve, generalized eigenproblem.
   Works in the browser (window.BeamCore) and in Node (module.exports).
   Degrees of freedom per node: (v, φ).  Global vector q = (v1,φ1,v2,φ2,…).
   ===================================================================== */
(function (root) {
  "use strict";

  // ---------- small dense linear algebra --------------------------------
  const zeros = (n, m = n) => Array.from({ length: n }, () => new Float64Array(m));

  function solve(A, b) {
    const n = b.length, M = A.map((r, i) => [...r, b[i]]);
    for (let k = 0; k < n; k++) {
      let p = k;
      for (let i = k + 1; i < n; i++) if (Math.abs(M[i][k]) > Math.abs(M[p][k])) p = i;
      [M[k], M[p]] = [M[p], M[k]];
      for (let i = k + 1; i < n; i++) {
        const f = M[i][k] / M[k][k];
        for (let j = k; j <= n; j++) M[i][j] -= f * M[k][j];
      }
    }
    const x = new Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) {
      let s = M[i][n];
      for (let j = i + 1; j < n; j++) s -= M[i][j] * x[j];
      x[i] = s / M[i][i];
    }
    return x;
  }

  // All eigenpairs of K v = λ M v (K symmetric, M SPD). Cholesky + cyclic Jacobi.
  // Returns { lambda: sorted ascending, vecs: array of Float64Array, M-orthonormal }.
  function eigGen(K, M) {
    const n = K.length;
    const L = zeros(n);
    for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) {
      let s = M[i][j];
      for (let k = 0; k < j; k++) s -= L[i][k] * L[j][k];
      if (i === j) { if (s <= 0) throw new Error("M non définie positive"); L[i][i] = Math.sqrt(s); }
      else L[i][j] = s / L[j][j];
    }
    const Y = zeros(n); // Y = L^-1 K
    for (let c = 0; c < n; c++) for (let i = 0; i < n; i++) {
      let s = K[i][c];
      for (let k = 0; k < i; k++) s -= L[i][k] * Y[k][c];
      Y[i][c] = s / L[i][i];
    }
    const C = zeros(n); // C = Y L^-T
    for (let r = 0; r < n; r++) for (let i = 0; i < n; i++) {
      let s = Y[r][i];
      for (let k = 0; k < i; k++) s -= L[i][k] * C[r][k];
      C[r][i] = s / L[i][i];
    }
    for (let i = 0; i < n; i++) for (let j = 0; j < i; j++) { const m = 0.5 * (C[i][j] + C[j][i]); C[i][j] = C[j][i] = m; }
    const V = zeros(n);
    for (let i = 0; i < n; i++) V[i][i] = 1;
    let tot = 0;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) tot += C[i][j] * C[i][j];
    for (let sweep = 0; sweep < 100; sweep++) {
      let off = 0;
      for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) off += C[i][j] * C[i][j];
      if (off < 1e-26 * tot) break;
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
    const idx = Array.from({ length: n }, (_, i) => i).sort((a, b) => C[a][a] - C[b][b]);
    const lambda = idx.map((i) => C[i][i]);
    const vecs = idx.map((col) => {
      const y = new Float64Array(n);
      for (let r = 0; r < n; r++) y[r] = V[r][col];
      const v = new Float64Array(n); // solve L^T v = y
      for (let i = n - 1; i >= 0; i--) {
        let s = y[i];
        for (let k = i + 1; k < n; k++) s -= L[k][i] * v[k];
        v[i] = s / L[i][i];
      }
      return v;
    });
    return { lambda, vecs };
  }

  // ---------- Hermite element ---------------------------------------------
  // N(ξ), N'(x), N''(x) for an element of length Le, ξ = x/Le in [0,1]
  function hermite(xi, Le) {
    const x2 = xi * xi, x3 = x2 * xi;
    return [1 - 3 * x2 + 2 * x3, Le * (xi - 2 * x2 + x3), 3 * x2 - 2 * x3, Le * (-x2 + x3)];
  }
  function hermiteD1(xi, Le) { // dN/dx
    const x2 = xi * xi;
    return [(-6 * xi + 6 * x2) / Le, 1 - 4 * xi + 3 * x2, (6 * xi - 6 * x2) / Le, -2 * xi + 3 * x2];
  }
  function hermiteD2(xi, Le) { // d²N/dx² = B
    return [(-6 + 12 * xi) / (Le * Le), (-4 + 6 * xi) / Le, (6 - 12 * xi) / (Le * Le), (-2 + 6 * xi) / Le];
  }

  function Ke(EI, Le) {
    const c = EI / Le ** 3, l = Le;
    return [
      [12 * c, 6 * l * c, -12 * c, 6 * l * c],
      [6 * l * c, 4 * l * l * c, -6 * l * c, 2 * l * l * c],
      [-12 * c, -6 * l * c, 12 * c, -6 * l * c],
      [6 * l * c, 2 * l * l * c, -6 * l * c, 4 * l * l * c],
    ];
  }
  function Me(rhoS, Le) {
    const c = rhoS * Le / 420, l = Le;
    return [
      [156 * c, 22 * l * c, 54 * c, -13 * l * c],
      [22 * l * c, 4 * l * l * c, 13 * l * c, -3 * l * l * c],
      [54 * c, 13 * l * c, 156 * c, -22 * l * c],
      [-13 * l * c, -3 * l * l * c, -22 * l * c, 4 * l * l * c],
    ];
  }
  // consistent load vector of a uniform load p (positive = downward): −p·Le/2·(1, Le/6, 1, −Le/6)
  function fe(p, Le) { return [-p * Le / 2, -p * Le * Le / 12, -p * Le / 2, p * Le * Le / 12]; }

  // ---------- assembly ------------------------------------------------------
  // Uniform mesh of Ne elements on [0, L]
  function assemble(Ne, L, EI, rhoS, p = 0) {
    const n = 2 * (Ne + 1), Le = L / Ne;
    const K = zeros(n), M = zeros(n), f = new Float64Array(n);
    const ke = Ke(EI, Le), me = Me(rhoS, Le), fl = fe(p, Le);
    for (let e = 0; e < Ne; e++) {
      const dof = [2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3];
      for (let a = 0; a < 4; a++) {
        f[dof[a]] += fl[a];
        for (let b = 0; b < 4; b++) { K[dof[a]][dof[b]] += ke[a][b]; M[dof[a]][dof[b]] += me[a][b]; }
      }
    }
    return { K, M, f, Le, n };
  }

  // Boundary conditions → list of blocked dofs
  function blocked(Ne, bc) {
    const last = 2 * Ne;
    if (bc === "cantilever") return [0, 1];
    if (bc === "ss") return [0, last];
    if (bc === "cc") return [0, 1, last, last + 1];
    if (bc === "free") return [];
    throw new Error("bc inconnue " + bc);
  }
  const reduce = (A, free) => free.map((i) => free.map((j) => A[i][j]));

  // ---------- static cantilever ---------------------------------------------
  // p: uniform load (downward positive); F: point force (downward) at node index kF (0..Ne) or null
  function staticCantilever(Ne, { L = 1, EI = 2e5, p = 0, F = 0, kF = null } = {}) {
    const { K, f, n } = assemble(Ne, L, EI, 1, p);
    if (F && kF != null) f[2 * kF] += -F;
    const blk = blocked(Ne, "cantilever");
    const free = Array.from({ length: n }, (_, i) => i).filter((i) => !blk.includes(i));
    const ql = solve(reduce(K, free), free.map((i) => f[i]));
    const q = new Float64Array(n);
    free.forEach((i, k) => (q[i] = ql[k]));
    // reactions r = K_bl q_l − f_b
    const r = blk.map((i) => { let s = -f[i]; for (let j = 0; j < n; j++) s += K[i][j] * q[j]; return s; });
    return { q, reaction: r, Le: L / Ne };
  }

  // Exact Euler–Bernoulli cantilever solution (displacement, moment M = EI v'')
  function exactCantilever(x, { L = 1, EI = 2e5, p = 0, F = 0, a = L / 2 } = {}) {
    let v = -p * x * x * (6 * L * L - 4 * L * x + x * x) / (24 * EI);
    let M = -p * (L - x) ** 2 / 2;
    if (F) {
      if (x <= a) { v += -F * x * x * (3 * a - x) / (6 * EI); M += -F * (a - x); }
      else { v += -F * a * a * (3 * x - a) / (6 * EI); }
    }
    return { v, M };
  }

  // Sample a FE solution (displacement and bending moment) with m points / element
  function sample(q, Ne, L, EI, m = 12) {
    const Le = L / Ne, x = [], v = [], M = [], Msecond = [];
    for (let e = 0; e < Ne; e++) {
      const qe = [q[2 * e], q[2 * e + 1], q[2 * e + 2], q[2 * e + 3]];
      for (let k = 0; k <= m; k++) {
        const xi = k / m;
        const N = hermite(xi, Le), B = hermiteD2(xi, Le);
        x.push(e * Le + xi * Le);
        v.push(N[0] * qe[0] + N[1] * qe[1] + N[2] * qe[2] + N[3] * qe[3]);
        M.push(EI * (B[0] * qe[0] + B[1] * qe[1] + B[2] * qe[2] + B[3] * qe[3]));
        Msecond.push(e);
      }
    }
    return { x, v, M, elem: Msecond };
  }

  // ---------- modal analysis -------------------------------------------------
  // Returns frequencies ω_i [rad/s] (ascending) and mass-normalised shapes on the full dof vector
  function modal(Ne, bc, { L = 1, EI = 2e5, rhoS = 0.78 } = {}) {
    // work in non-dimensional form (EI=ρS=L=1) for conditioning, rescale ω afterwards
    const { K, M, n } = assemble(Ne, 1, 1, 1);
    const blk = blocked(Ne, bc);
    const free = Array.from({ length: n }, (_, i) => i).filter((i) => !blk.includes(i));
    const { lambda, vecs } = eigGen(reduce(K, free), reduce(M, free));
    const scale = Math.sqrt(EI / (rhoS * L ** 4));
    const omega = lambda.map((l) => Math.sqrt(Math.max(l, 0)) * scale);
    const shapes = vecs.map((v) => { const full = new Float64Array(n); free.forEach((i, k) => (full[i] = v[k])); return full; });
    return { omega, shapes, free, Le: L / Ne, n };
  }

  // Analytical natural frequencies  ω_i = (β_i L)² √(EI/ρS L⁴)
  function betaL(bc, i) { // i = 1,2,…
    const cant = [1.8751040687, 4.6940911330, 7.8547574382, 10.9955407349, 14.1371683911];
    const cc = [4.7300407449, 7.8532046241, 10.9956078380, 14.1371654913, 17.2787596574];
    if (bc === "ss") return i * Math.PI;
    if (bc === "cantilever") return i <= 5 ? cant[i - 1] : (2 * i - 1) * Math.PI / 2;
    if (bc === "cc") return i <= 5 ? cc[i - 1] : (2 * i + 1) * Math.PI / 2;
    return NaN;
  }
  function omegaExact(bc, i, { L = 1, EI = 2e5, rhoS = 0.78 } = {}) {
    return betaL(bc, i) ** 2 * Math.sqrt(EI / (rhoS * L ** 4));
  }

  // Interpolated shape (Hermite) of a global dof vector, sampled with m points per element
  function shapeCurve(vec, Ne, m = 8) {
    return sample(vec, Ne, 1, 1, m);
  }

  root.BeamCore = {
    solve, eigGen, hermite, hermiteD1, hermiteD2, Ke, Me, fe, assemble, blocked, reduce,
    staticCantilever, exactCantilever, sample, modal, betaL, omegaExact, shapeCurve,
  };
  if (typeof module !== "undefined" && module.exports) module.exports = root.BeamCore;
})(typeof window !== "undefined" ? window : globalThis);
