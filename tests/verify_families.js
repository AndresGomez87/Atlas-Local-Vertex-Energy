// Verificación de los generadores de familias del Atlas.
// Cada generador devuelve {n, edges:[[i,j],...]} con índices 0..n-1.

function eigSymm(M) {
  const n = M.length;
  if (n === 0) return [];
  if (n === 1) return [M[0][0]];
  const A = M.map(r => Float64Array.from(r));
  let iter = 0;
  while (iter++ < 100 * n * n) {
    let p = 0, q = 1, mx = Math.abs(A[0][1]);
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++)
      if (Math.abs(A[i][j]) > mx) { mx = Math.abs(A[i][j]); p = i; q = j; }
    if (mx < 1e-12) break;
    const app = A[p][p], aqq = A[q][q], apq = A[p][q];
    const tau = (aqq - app) / (2 * apq);
    const t = Math.sign(tau !== 0 ? tau : 1) / (Math.abs(tau) + Math.sqrt(tau * tau + 1));
    const c = 1 / Math.sqrt(t * t + 1), s = t * c;
    A[p][p] = app - t * apq; A[q][q] = aqq + t * apq; A[p][q] = A[q][p] = 0;
    for (let i = 0; i < n; i++) if (i !== p && i !== q) {
      const api = A[p][i], aqi = A[q][i];
      A[p][i] = c * api - s * aqi; A[i][p] = A[p][i];
      A[q][i] = s * api + c * aqi; A[i][q] = A[q][i];
    }
  }
  return Array.from(A, (_, i) => A[i][i]);
}
function adjOf(g) {
  const A = Array.from({ length: g.n }, () => Array(g.n).fill(0));
  for (const [a, b] of g.edges) { A[a][b] = 1; A[b][a] = 1; }
  return A;
}
const energyOf = g => eigSymm(adjOf(g)).reduce((s, e) => s + Math.abs(e), 0);

// ── GENERADORES (idénticos a los que irán en el HTML) ──
const GEN = {
  path:  ({ n }) => ({ n, edges: Array.from({ length: n - 1 }, (_, i) => [i, i + 1]) }),
  cycle: ({ n }) => ({ n, edges: Array.from({ length: n }, (_, i) => [i, (i + 1) % n]) }),
  star:  ({ n }) => ({ n, edges: Array.from({ length: n - 1 }, (_, i) => [0, i + 1]) }),
  complete: ({ n }) => {
    const edges = [];
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) edges.push([i, j]);
    return { n, edges };
  },
  bipartite: ({ p, q }) => {
    const edges = [];
    for (let i = 0; i < p; i++) for (let j = 0; j < q; j++) edges.push([i, p + j]);
    return { n: p + q, edges };
  },
  // Broom B_{r,q}: camino v0…v_{r-1}; q hojas pendientes en v_{r-1} (el hub).
  broom: ({ r, q }) => {
    const edges = [];
    for (let i = 0; i < r - 1; i++) edges.push([i, i + 1]);
    for (let j = 0; j < q; j++) edges.push([r - 1, r + j]);
    return { n: r + q, edges };
  },
  // Double star S_{p,q}: dos centros adyacentes c1=0, c2=1; p hojas en c1, q en c2.
  dstar: ({ p, q }) => {
    const edges = [[0, 1]];
    for (let i = 0; i < p; i++) edges.push([0, 2 + i]);
    for (let j = 0; j < q; j++) edges.push([1, 2 + p + j]);
    return { n: 2 + p + q, edges };
  },
  // Double broom simétrico DB_{t,k}: camino de t vértices (v0…v_{t-1});
  // k hojas pendientes en cada extremo del camino.
  dbroom: ({ t, k }) => {
    const edges = [];
    for (let i = 0; i < t - 1; i++) edges.push([i, i + 1]);
    for (let j = 0; j < k; j++) edges.push([0, t + j]);
    for (let j = 0; j < k; j++) edges.push([t - 1, t + k + j]);
    return { n: t + 2 * k, edges };
  },
  petersen: () => ({ n: 10, edges: [[0,1],[1,2],[2,3],[3,4],[4,0],[5,7],[7,9],[9,6],[6,8],[8,5],[0,5],[1,6],[2,7],[3,8],[4,9]] })
};

let fails = 0;
function check(name, got, exp, tol = 1e-9) {
  const ok = Math.abs(got - exp) < tol;
  if (!ok) fails++;
  console.log(`${ok ? 'OK  ' : 'FAIL'} ${name}: ${got.toPrecision(10)} vs ${exp.toPrecision(10)}`);
}
function checkStruct(name, g, n, m) {
  const ok = g.n === n && g.edges.length === m;
  if (!ok) fails++;
  console.log(`${ok ? 'OK  ' : 'FAIL'} ${name}: n=${g.n} (esp ${n}), m=${g.edges.length} (esp ${m})`);
}

// estructura
checkStruct('P_6', GEN.path({ n: 6 }), 6, 5);
checkStruct('C_7', GEN.cycle({ n: 7 }), 7, 7);
checkStruct('S_8', GEN.star({ n: 8 }), 8, 7);
checkStruct('K_5', GEN.complete({ n: 5 }), 5, 10);
checkStruct('K_{3,4}', GEN.bipartite({ p: 3, q: 4 }), 7, 12);
checkStruct('B_{4,3}', GEN.broom({ r: 4, q: 3 }), 7, 6);
checkStruct('S_{2,3} (double star)', GEN.dstar({ p: 2, q: 3 }), 7, 6);
checkStruct('DB_{4,3}', GEN.dbroom({ t: 4, k: 3 }), 10, 9);
checkStruct('Petersen', GEN.petersen(), 10, 15);

// energías exactas
check('E(P_5)', energyOf(GEN.path({ n: 5 })), [1,2,3,4,5].map(k=>Math.abs(2*Math.cos(k*Math.PI/6))).reduce((a,b)=>a+b));
check('E(S_9)', energyOf(GEN.star({ n: 9 })), 2 * Math.sqrt(8));
check('E(K_6)', energyOf(GEN.complete({ n: 6 })), 10);
check('E(K_{3,4})', energyOf(GEN.bipartite({ p: 3, q: 4 })), 2 * Math.sqrt(12));
check('E(B_{4,3})', energyOf(GEN.broom({ r: 4, q: 3 })), 2 * Math.sqrt(6 + 2 * Math.sqrt(7)));
check('E(Petersen)', energyOf(GEN.petersen()), 16);

// double star S_{p,q}: polinomio caracteristico conocido ->
// verificamos contra diagonalización directa de una construcción independiente (matriz a mano)
{
  const g = GEN.dstar({ p: 3, q: 3 });
  // S_{3,3}: n=8. Construcción independiente: adyacencia explícita
  const A = Array.from({length:8},()=>Array(8).fill(0));
  const put=(a,b)=>{A[a][b]=1;A[b][a]=1;};
  put(0,1); put(0,2);put(0,3);put(0,4); put(1,5);put(1,6);put(1,7);
  const Eref = eigSymm(A).reduce((s,e)=>s+Math.abs(e),0);
  check('E(S_{3,3}) vs construcción independiente', energyOf(g), Eref);
}

// double broom: sanidad estructural — grados
{
  const g = GEN.dbroom({ t: 4, k: 3 });
  const A = adjOf(g);
  const deg = A.map(r => r.reduce((a,b)=>a+b,0));
  // extremos del camino: grado k+1 = 4; internos: 2; hojas: 1
  const ok = deg[0] === 4 && deg[3] === 4 && deg[1] === 2 && deg[2] === 2 && deg.slice(4).every(d => d === 1);
  if (!ok) fails++;
  console.log(`${ok ? 'OK  ' : 'FAIL'} DB_{4,3} grados: [${deg.join(',')}]`);
  // energía del hub vs hoja (los valores de la tesis: hub ≈ 3.83, hoja ≈ 0.59)
  const EG = energyOf(g);
  const del = k => { const idx=[...Array(g.n).keys()].filter(i=>i!==k);
    const B = idx.map(i=>idx.map(j=>A[i][j])); return eigSymm(B).reduce((s,e)=>s+Math.abs(e),0); };
  console.log(`     DB_{4,3}: E_G(hub v0)=${(EG-del(0)).toFixed(4)}  E_G(hoja)=${(EG-del(4)).toFixed(4)}`);
}

console.log(fails === 0 ? '\n=== GENERADORES VERIFICADOS ===' : `\n=== ${fails} FALLOS ===`);
process.exit(fails === 0 ? 0 : 1);
