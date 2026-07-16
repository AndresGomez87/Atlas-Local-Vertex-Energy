// Verificación del motor de eigenvalores (Jacobi clásico) extraído del HTML.
// Compara contra valores exactos de forma cerrada y contra referencias NumPy.

function eigSymm(M) {
  const n = M.length;
  if (n === 0) return [];
  if (n === 1) return [M[0][0]];
  const A = M.map(row => Float64Array.from(row));
  const maxIter = 100 * n * n;
  let iter = 0;
  while (iter++ < maxIter) {
    let p = 0, q = 1, maxVal = Math.abs(A[0][1]);
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        if (Math.abs(A[i][j]) > maxVal) { maxVal = Math.abs(A[i][j]); p = i; q = j; }
      }
    }
    if (maxVal < 1e-12) break;
    const app = A[p][p], aqq = A[q][q], apq = A[p][q];
    const tau = (aqq - app) / (2.0 * apq);
    const t = Math.sign(tau !== 0 ? tau : 1) / (Math.abs(tau) + Math.sqrt(tau * tau + 1.0));
    const c = 1.0 / Math.sqrt(t * t + 1.0);
    const s = t * c;
    A[p][p] = app - t * apq;
    A[q][q] = aqq + t * apq;
    A[p][q] = 0.0; A[q][p] = 0.0;
    for (let i = 0; i < n; i++) {
      if (i !== p && i !== q) {
        const api = A[p][i], aqi = A[q][i];
        A[p][i] = c * api - s * aqi; A[i][p] = A[p][i];
        A[q][i] = s * api + c * aqi; A[i][q] = A[q][i];
      }
    }
  }
  return Array.from(A, (_, i) => A[i][i]);
}

const energy = A => eigSymm(A).reduce((s, e) => s + Math.abs(e), 0);

// ---- constructores de familias ----
const zeros = n => Array.from({length: n}, () => Array(n).fill(0));
function fromEdges(n, E) { const A = zeros(n); for (const [a,b] of E) { A[a][b]=1; A[b][a]=1; } return A; }
function path(n){ const E=[]; for(let i=0;i<n-1;i++) E.push([i,i+1]); return fromEdges(n,E); }
function cycle(n){ const E=[]; for(let i=0;i<n;i++) E.push([i,(i+1)%n]); return fromEdges(n,E); }
function star(n){ const E=[]; for(let i=1;i<n;i++) E.push([0,i]); return fromEdges(n,E); }
function complete(n){ const E=[]; for(let i=0;i<n;i++) for(let j=i+1;j<n;j++) E.push([i,j]); return fromEdges(n,E); }
function completeBipartite(p,q){ const E=[]; for(let i=0;i<p;i++) for(let j=0;j<q;j++) E.push([i,p+j]); return fromEdges(p+q,E); }
// Broom B_{r,q}: camino de r vértices, con q hojas pegadas al último vértice del camino.
// B_{4,3}: camino w3-w2-w1-v (4 vértices), v con 3 hojas.  n = 7.
function broom(r,q){ const E=[]; for(let i=0;i<r-1;i++) E.push([i,i+1]); for(let j=0;j<q;j++) E.push([r-1, r+j]); return fromEdges(r+q,E); }
function petersen(){
  const E=[[0,1],[1,2],[2,3],[3,4],[4,0],[5,7],[7,9],[9,6],[6,8],[8,5],[0,5],[1,6],[2,7],[3,8],[4,9]];
  return fromEdges(10,E);
}

// Energía local: E_G(v) = E(G) - E(G - v)
function deleteVertex(A, k){
  const n = A.length, B = zeros(n-1);
  for(let i=0,bi=0;i<n;i++){ if(i===k) continue; for(let j=0,bj=0;j<n;j++){ if(j===k) continue; B[bi][bj]=A[i][j]; bj++; } bi++; }
  return B;
}
const localEnergies = A => A.map((_, k) => energy(A) - energy(deleteVertex(A, k)));

// ---- pruebas contra formas cerradas ----
let fails = 0;
function check(name, got, expected, tol = 1e-9) {
  const err = Math.abs(got - expected);
  const ok = err < tol;
  if (!ok) fails++;
  console.log(`${ok ? 'OK  ' : 'FAIL'} ${name}: got=${got.toPrecision(12)} exp=${expected.toPrecision(12)} err=${err.toExponential(2)}`);
}

// Espectros exactos
// P_n: eigenvalores 2cos(kπ/(n+1)), E(P_n) = suma |.|
for (const n of [2,3,5,10,20]) {
  const exact = Array.from({length:n},(_,k)=>2*Math.cos((k+1)*Math.PI/(n+1))).reduce((s,x)=>s+Math.abs(x),0);
  check(`E(P_${n})`, energy(path(n)), exact);
}
// C_n: 2cos(2πk/n)
for (const n of [3,4,6,11]) {
  const exact = Array.from({length:n},(_,k)=>2*Math.cos(2*Math.PI*k/n)).reduce((s,x)=>s+Math.abs(x),0);
  check(`E(C_${n})`, energy(cycle(n)), exact, 1e-8);
}
// S_n (estrella con n vértices): E = 2√(n−1)
for (const n of [4,7,12]) check(`E(S_${n})`, energy(star(n)), 2*Math.sqrt(n-1));
// K_n: E = 2(n−1)
for (const n of [3,5,8]) check(`E(K_${n})`, energy(complete(n)), 2*(n-1));
// K_{p,q}: E = 2√(pq)
for (const [p,q] of [[2,3],[3,3],[2,6],[4,5]]) check(`E(K_{${p},${q}})`, energy(completeBipartite(p,q)), 2*Math.sqrt(p*q));
// Petersen: espectro 3, 1^5, (−2)^4 → E = 3+5+8 = 16
check('E(Petersen)', energy(petersen()), 16);
// B_{4,3}: E = 2√(6+2√7)  (verificado en la tesis)
check('E(B_{4,3})', energy(broom(4,3)), 2*Math.sqrt(6+2*Math.sqrt(7)));

// Energías locales de B_{4,3} (valores verificados previamente con NumPy en la tesis):
// hoja u ≈ 0.565, hub v ≈ 3.892, extremo del camino w ≈ 0.902
{
  const le = localEnergies(broom(4,3));
  // vértices: 0..3 camino (3 = hub v), 4..6 hojas
  const EB = 2*Math.sqrt(6+2*Math.sqrt(7));
  check('E_v(B43) hub', le[3], EB - 2*Math.sqrt(2));           // B−v = P3 ∪ 3K1, E = 2√2
  console.log(`     hub=${le[3].toFixed(4)} hoja=${le[4].toFixed(4)} extremo=${le[0].toFixed(4)}`);
}

// Propiedades estructurales sobre grafos aleatorios:
// (1) E_G(v) ≥ 0 (submatriz principal: norma traza no aumenta)
// (2) E_G(v) ≤ 2√(d_v)  (Espinal–Rada, Teorema 2)
// (3) e(G) ≤ 2 E(G)
function randomGraph(n, p, seed) {
  let s = seed;
  const rnd = () => (s = (s*1103515245+12345) & 0x7fffffff) / 0x7fffffff;
  const E=[]; for(let i=0;i<n;i++) for(let j=i+1;j<n;j++) if(rnd()<p) E.push([i,j]);
  return fromEdges(n,E);
}
let propFails = 0;
for (let trial = 0; trial < 60; trial++) {
  const n = 3 + (trial % 10);
  const A = randomGraph(n, 0.4, 1000+trial);
  const EG = energy(A);
  const le = localEnergies(A);
  const eG = le.reduce((s,x)=>s+x,0);
  le.forEach((E_v, v) => {
    const d = A[v].reduce((s,x)=>s+x,0);
    if (E_v < -1e-8) { propFails++; console.log(`FAIL E_G(v)>=0 en trial ${trial}, v=${v}: ${E_v}`); }
    if (E_v > 2*Math.sqrt(d) + 1e-8) { propFails++; console.log(`FAIL cota 2√d en trial ${trial}, v=${v}`); }
  });
  if (eG > 2*EG + 1e-8) { propFails++; console.log(`FAIL e(G)<=2E(G) en trial ${trial}`); }
}
console.log(propFails === 0 ? 'OK   60 grafos aleatorios: 0 ≤ E_G(v) ≤ 2√d_v y e(G) ≤ 2E(G) en todos' : `FAIL propiedades: ${propFails}`);

// Exportar casos aleatorios para comparación cruzada con NumPy
const cases = [];
for (let t = 0; t < 25; t++) {
  const n = 3 + (t % 12);
  const A = randomGraph(n, 0.45, 777+t);
  cases.push({ A, E: energy(A), le: localEnergies(A) });
}
require('fs').writeFileSync(require('path').join(__dirname, 'js_results.json'), JSON.stringify(cases));
console.log(`\nExportados 25 casos aleatorios a js_results.json para cruce con NumPy.`);
console.log(fails === 0 && propFails === 0 ? '\n=== MOTOR VERIFICADO: todas las pruebas de forma cerrada pasan ===' : `\n=== ${fails+propFails} FALLOS ===`);
