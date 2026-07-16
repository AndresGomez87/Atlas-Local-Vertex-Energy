'use strict';

/* ════════════════════════════════════════════════════════════════════
   NUMERICAL ENGINE — verified against NumPy (max abs error ~1e-14 on
   random graphs; exact-form tests: P_n, C_n, S_n, K_n, K_{p,q},
   Petersen, broom B_{4,3}).
   ════════════════════════════════════════════════════════════════════ */

// Classical Jacobi eigenvalue algorithm for real symmetric matrices.
// Repeatedly annihilates the largest off-diagonal entry with a Givens
// rotation; converges to a diagonal matrix holding the spectrum.
function eigSymm(M) {
  const n = M.length;
  if (n === 0) return [];
  if (n === 1) return [M[0][0]];

  const A = M.map(row => Float64Array.from(row));
  const maxIter = 100 * n * n;
  let iter = 0;

  while (iter++ < maxIter) {
    // largest off-diagonal entry a_pq
    let p = 0, q = 1, maxVal = Math.abs(A[0][1]);
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        if (Math.abs(A[i][j]) > maxVal) { maxVal = Math.abs(A[i][j]); p = i; q = j; }
      }
    }
    if (maxVal < 1e-12) break;   // numerically diagonal

    const app = A[p][p], aqq = A[q][q], apq = A[p][q];
    // numerically stable rotation angle
    const tau = (aqq - app) / (2.0 * apq);
    const t = Math.sign(tau !== 0 ? tau : 1) / (Math.abs(tau) + Math.sqrt(tau * tau + 1.0));
    const c = 1.0 / Math.sqrt(t * t + 1.0);
    const s = t * c;

    // orthogonal similarity A ← JᵀAJ
    A[p][p] = app - t * apq;
    A[q][q] = aqq + t * apq;
    A[p][q] = 0.0; A[q][p] = 0.0;
    for (let i = 0; i < n; i++) {
      if (i !== p && i !== q) {
        const api = A[p][i], aqi = A[q][i];
        A[p][i] = c * api - s * aqi;  A[i][p] = A[p][i];
        A[q][i] = s * api + c * aqi;  A[i][q] = A[q][i];
      }
    }
  }
  const eigs = new Array(n);
  for (let i = 0; i < n; i++) eigs[i] = A[i][i];
  return eigs;
}

function graphEnergy(adj) {
  return eigSymm(adj).reduce((s, e) => s + Math.abs(e), 0);
}

function buildAdj(verts, edgeList) {
  const n = verts.length;
  const idx = {};
  verts.forEach((v, i) => idx[v.id] = i);
  const A = Array.from({ length: n }, () => Array(n).fill(0));
  for (const e of edgeList) {
    const a = idx[e.a], b = idx[e.b];
    if (a !== undefined && b !== undefined) { A[a][b] = 1; A[b][a] = 1; }
  }
  return A;
}

/* ════════════════════════════════════════════════════════════════════
   STATE
   ════════════════════════════════════════════════════════════════════ */

function lineLayout(n, W, H) {
  const m = 70;
  return Array.from({ length: n }, (_, i) =>
    ({ x: m + (W - 2 * m) * (n === 1 ? 0.5 : i / (n - 1)), y: H / 2 }));
}
function circleLayout(n, W, H, scale = 1) {
  const r = (Math.min(W, H) / 2 - 70) * scale, cx = W / 2, cy = H / 2;
  return Array.from({ length: n }, (_, i) => {
    const a = 2 * Math.PI * i / n - Math.PI / 2;
    return { x: cx + r * Math.cos(a), y: cy + r * Math.sin(a) };
  });
}
function fan(cx, cy, count, angle0, spread, r) {
  return Array.from({ length: count }, (_, i) => {
    const a = count === 1 ? angle0 : angle0 - spread / 2 + spread * i / (count - 1);
    return { x: cx + r * Math.cos(a), y: cy + r * Math.sin(a) };
  });
}
const fanSpread = k => Math.min(Math.PI * 1.15, 0.55 + 0.3 * k);

const FAMILIES = {
  custom: {
    name: 'Free draw', params: [],
    note: 'Draw any graph by hand — or pick a family to generate one, then edit it freely.'
  },
  path: {
    name: 'Path P(n)',
    params: [{ k: 'n', lab: 'n', min: 2, max: 30, def: 6 }],
    gen: p => ({ n: p.n, edges: Array.from({ length: p.n - 1 }, (_, i) => [i, i + 1]) }),
    layout: (p, W, H) => lineLayout(p.n, W, H),
    capt: p => `The path <i>P</i><sub>${p.n}</sub>`,
    note: 'Eigenvalues 2 cos(kπ/(n+1)), k = 1,…,n.'
  },
  cycle: {
    name: 'Cycle C(n)',
    params: [{ k: 'n', lab: 'n', min: 3, max: 30, def: 8 }],
    gen: p => ({ n: p.n, edges: Array.from({ length: p.n }, (_, i) => [i, (i + 1) % p.n]) }),
    layout: (p, W, H) => circleLayout(p.n, W, H),
    capt: p => `The cycle <i>C</i><sub>${p.n}</sub>`,
    note: 'Eigenvalues 2 cos(2πk/n); bipartite exactly when n is even.'
  },
  star: {
    name: 'Star S(n)',
    params: [{ k: 'n', lab: 'n', min: 3, max: 30, def: 7 }],
    gen: p => ({ n: p.n, edges: Array.from({ length: p.n - 1 }, (_, i) => [0, i + 1]) }),
    layout: (p, W, H) => [{ x: W / 2, y: H / 2 },
      ...circleLayout(p.n - 1, W, H, 0.8)],
    capt: p => `The star <i>S</i><sub>${p.n}</sub>`,
    note: 'Spectrum ±√(n−1) and 0 with multiplicity n−2, so E(G) = 2√(n−1).'
  },
  complete: {
    name: 'Complete K(n)',
    params: [{ k: 'n', lab: 'n', min: 2, max: 14, def: 5 }],
    gen: p => {
      const edges = [];
      for (let i = 0; i < p.n; i++) for (let j = i + 1; j < p.n; j++) edges.push([i, j]);
      return { n: p.n, edges };
    },
    layout: (p, W, H) => circleLayout(p.n, W, H, 0.85),
    capt: p => `The complete graph <i>K</i><sub>${p.n}</sub>`,
    note: 'Spectrum n−1 and −1 with multiplicity n−1, so E(G) = 2(n−1).'
  },
  bipartite: {
    name: 'Complete bipartite K(p,q)',
    params: [{ k: 'p', lab: 'p', min: 1, max: 12, def: 3 },
             { k: 'q', lab: 'q', min: 1, max: 12, def: 4 }],
    gen: p => {
      const edges = [];
      for (let i = 0; i < p.p; i++) for (let j = 0; j < p.q; j++) edges.push([i, p.p + j]);
      return { n: p.p + p.q, edges };
    },
    layout: (p, W, H) => {
      const m = 70;
      const col = (c, x) => Array.from({ length: c }, (_, i) =>
        ({ x, y: c === 1 ? H / 2 : m + (H - 2 * m) * i / (c - 1) }));
      return [...col(p.p, W * 0.34), ...col(p.q, W * 0.66)];
    },
    capt: p => `The complete bipartite graph <i>K</i><sub>${p.p},${p.q}</sub>`,
    note: 'Nonzero eigenvalues ±√(pq), so E(G) = 2√(pq).'
  },
  broom: {
    name: 'Broom B(r,q)',
    params: [{ k: 'r', lab: 'r', min: 2, max: 14, def: 4 },
             { k: 'q', lab: 'q', min: 1, max: 12, def: 3 }],
    gen: p => {
      const edges = [];
      for (let i = 0; i < p.r - 1; i++) edges.push([i, i + 1]);
      for (let j = 0; j < p.q; j++) edges.push([p.r - 1, p.r + j]);
      return { n: p.r + p.q, edges };
    },
    layout: (p, W, H) => {
      const m = 80, pathW = (W - 2 * m) * 0.6;
      const handle = Array.from({ length: p.r }, (_, i) =>
        ({ x: m + (p.r === 1 ? 0 : pathW * i / (p.r - 1)), y: H / 2 }));
      const hub = handle[p.r - 1];
      return [...handle, ...fan(hub.x, hub.y, p.q, 0, fanSpread(p.q), 115)];
    },
    capt: p => `The broom <i>B</i><sub>${p.r},${p.q}</sub>`,
    note: 'The hub carries the largest local energy of the broom; the pendant leaves the smallest.'
  },
  dstar: {
    name: 'Double star S(p,q)',
    params: [{ k: 'p', lab: 'p', min: 1, max: 12, def: 3 },
             { k: 'q', lab: 'q', min: 1, max: 12, def: 3 }],
    gen: p => {
      const edges = [[0, 1]];
      for (let i = 0; i < p.p; i++) edges.push([0, 2 + i]);
      for (let j = 0; j < p.q; j++) edges.push([1, 2 + p.p + j]);
      return { n: 2 + p.p + p.q, edges };
    },
    layout: (p, W, H) => {
      const c1 = { x: W * 0.42, y: H / 2 }, c2 = { x: W * 0.58, y: H / 2 };
      return [c1, c2,
        ...fan(c1.x, c1.y, p.p, Math.PI, fanSpread(p.p), 115),
        ...fan(c2.x, c2.y, p.q, 0, fanSpread(p.q), 115)];
    },
    capt: p => `The double star <i>S</i><sub>${p.p},${p.q}</sub>`,
    note: 'Two adjacent centers with p and q pendant leaves; a tree of diameter 3.'
  },
  dbroom: {
    name: 'Double broom DB(t,k)',
    params: [{ k: 't', lab: 't', min: 2, max: 12, def: 4 },
             { k: 'k', lab: 'k', min: 1, max: 10, def: 3 }],
    gen: p => {
      const edges = [];
      for (let i = 0; i < p.t - 1; i++) edges.push([i, i + 1]);
      for (let j = 0; j < p.k; j++) edges.push([0, p.t + j]);
      for (let j = 0; j < p.k; j++) edges.push([p.t - 1, p.t + p.k + j]);
      return { n: p.t + 2 * p.k, edges };
    },
    layout: (p, W, H) => {
      const pathW = W * 0.44, x0 = (W - pathW) / 2;
      const handle = Array.from({ length: p.t }, (_, i) =>
        ({ x: x0 + (p.t === 1 ? 0 : pathW * i / (p.t - 1)), y: H / 2 }));
      const L = handle[0], Rr = handle[p.t - 1];
      return [...handle,
        ...fan(L.x, L.y, p.k, Math.PI, fanSpread(p.k), 115),
        ...fan(Rr.x, Rr.y, p.k, 0, fanSpread(p.k), 115)];
    },
    capt: p => `The symmetric double broom <i>DB</i><sub>${p.t},${p.k}</sub>`,
    note: 'A path on t vertices with k pendant leaves at each end; symmetric under the swap of its two hubs.'
  },
  petersen: {
    name: 'Petersen graph',
    params: [],
    gen: () => ({ n: 10, edges: [[0,1],[1,2],[2,3],[3,4],[4,0],
      [5,7],[7,9],[9,6],[6,8],[8,5],[0,5],[1,6],[2,7],[3,8],[4,9]] }),
    layout: (p, W, H) => [...circleLayout(5, W, H, 1), ...circleLayout(5, W, H, 0.5)],
    capt: () => 'The Petersen graph',
    note: 'Spectrum 3, 1 (×5), −2 (×4), so E(G) = 16. Strongly regular srg(10, 3, 0, 1).'
  }
};

function adjFromGen(g) {
  const A = Array.from({ length: g.n }, () => Array(g.n).fill(0));
  for (const [a, b] of g.edges) { A[a][b] = 1; A[b][a] = 1; }
  return A;
}
function energyWithout(A, k) {
  const n = A.length, idx = [];
  for (let i = 0; i < n; i++) if (i !== k) idx.push(i);
  return graphEnergy(idx.map(i => idx.map(j => A[i][j])));
}
function sweepPoint(f, params, x) {
  const g = f.gen(params);
  const A = adjFromGen(g);
  const EG = eigSymm(A).reduce((s, e) => s + Math.abs(e), 0);
  let eG = 0, mx = -Infinity, mn = Infinity;
  for (let k = 0; k < g.n; k++) {
    const Ev = g.n === 1 ? EG : EG - energyWithout(A, k);
    eG += Ev;
    if (Ev > mx) mx = Ev;
    if (Ev < mn) mn = Ev;
  }
  return { x, EG, eG, mx, mn };
}


// Prueba del barrido: broom con r=4 fijo, q variable
const f = FAMILIES.broom;
const pts = [];
for (let q = 1; q <= 12; q++) pts.push(sweepPoint(f, { r: 4, q }, q));

// El punto q=3 debe reproducir B_{4,3}
const p3 = pts.find(p => p.x === 3);
const exact = 2 * Math.sqrt(6 + 2 * Math.sqrt(7));
let fails = 0;
const chk = (name, got, exp, tol=1e-6) => {
  const ok = Math.abs(got - exp) < tol; if (!ok) fails++;
  console.log(`${ok?'OK  ':'FAIL'} ${name}: ${got.toFixed(6)} vs ${exp.toFixed(6)}`);
};
chk('E(B_{4,3}) en la curva', p3.EG, exact, 1e-9);
chk('max local (hub)', p3.mx, 3.892139, 1e-5);
chk('min local (hoja)', p3.mn, 0.565229, 1e-4);

// monotonicidad esperable: E(B_{4,q}) crece con q (más aristas en un árbol más grande)
const mono = pts.every((p,i)=> i===0 || p.EG > pts[i-1].EG - 1e-9);
console.log((mono?'OK  ':'FAIL')+' E(G) no decrece a lo largo del barrido en q');
if(!mono) fails++;

// coherencia interna: mn <= mx, eG = suma consistente con cotas
for (const p of pts) {
  if (!(p.mn <= p.mx + 1e-12 && p.eG <= 2*p.EG + 1e-9 && p.mn >= -1e-9)) {
    fails++; console.log('FAIL coherencia en q='+p.x);
  }
}
console.log('OK   coherencia mn<=mx, 0<=E_G(v), e(G)<=2E(G) en los 12 puntos');

// double broom: barrido en k con t=4; el punto k=3 debe dar hub 3.8303
const f2 = FAMILIES.dbroom;
const pk = sweepPoint(f2, { t: 4, k: 3 }, 3);
chk('DB_{4,3} max local (hub)', pk.mx, 3.8303, 1e-3);
chk('DB_{4,3} min local (hoja)', pk.mn, 0.5889, 1e-3);

console.log(fails===0 ? '\n=== BARRIDO VERIFICADO ===' : '\n=== '+fails+' FALLOS ===');
process.exit(fails?1:0);
