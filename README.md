# Atlas of Local Vertex Energy

[![tests](https://github.com/AndresGomez87/atlas-local-vertex-energy/actions/workflows/tests.yml/badge.svg)](https://github.com/AndresGomez87/atlas-local-vertex-energy/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-0b4a2f.svg)](LICENSE)

**[Versión en español →](README.es.md)**

An interactive atlas for exploring the **local energy of vertices** in graphs — a recent invariant from spectral graph theory. Draw any graph or generate one from a parametric family, and watch the energy distribute itself across the vertices as a live heatmap, with the spectrum plotted on the real line and energy curves traced across the whole family.

**Live demo:** https://AndresGomez87.github.io/atlas-local-vertex-energy/

![The broom B(6,5) with its local energy heatmap](assets/hero.png)

*The broom B(6,5). The hub concentrates the largest local energy (deep green) while the pendant leaves carry the least. Note the alternating intensity along the handle.*

## The mathematics

The **energy** of a graph *G* is the sum of the absolute values of the eigenvalues of its adjacency matrix, E(G) = Σᵢ |λᵢ|, a classical invariant introduced by Gutman with roots in Hückel molecular orbital theory. The **local energy of a vertex**, introduced by Espinal and Rada (2024), measures how much of that energy the vertex is responsible for:

> E_G(v) = E(G) − E(G − v)

Since the adjacency matrix of G − v is a principal submatrix of that of G, interlacing gives 0 ≤ E_G(v), and Espinal–Rada prove the sharp upper bound E_G(v) ≤ 2√(deg v). The **total local energy** e(G) = Σᵥ E_G(v) satisfies e(G) ≤ 2·E(G). The Atlas verifies and displays all three quantities live.

## Features

- **Free drawing** — click to place vertices, join them with edges, drag, delete. Everything recomputes as you draw (up to n = 30; on-demand beyond that).
- **Parametric families** — paths, cycles, stars, complete and complete bipartite graphs, brooms B(r,q), double stars S(p,q), symmetric double brooms DB(t,k), and the Petersen graph. One slider per parameter; the graph regenerates and lays itself out automatically. Editing a generated graph by hand detaches it into free-draw mode.
- **Local energy heatmap** — every vertex is coloured by its E_G(v), with a colour bar. Structure becomes visible: hubs glow, leaves fade.
- **Spectrum on the real line** — eigenvalues plotted as points on ℝ, stacked by multiplicity, with hollow points at zero. The symmetry of bipartite spectra is visible at a glance, and the figure caption notes it.
- **Family energy curves (Fig. 2)** — E(G), e(G), maxᵥ E_G(v) and minᵥ E_G(v) traced across the family as a parameter sweeps its range, computed progressively and cached. A dotted rule marks where the currently displayed graph sits inside its own family.
- **Shareable URLs** — `index.html#f=broom&r=4&q=3` reconstructs that exact graph. **PNG export** of the plate with one click.
- **Zero dependencies** — a single HTML file. The whole numerical stack is implemented from scratch.

## Parametric families

| Family | Parameters | Spectral fact shown in-app |
|---|---|---|
| Path P(n) | n | eigenvalues 2cos(kπ/(n+1)) |
| Cycle C(n) | n | eigenvalues 2cos(2πk/n); bipartite iff n even |
| Star S(n) | n | spectrum ±√(n−1), 0ⁿ⁻²; E = 2√(n−1) |
| Complete K(n) | n | spectrum n−1, (−1)ⁿ⁻¹; E = 2(n−1) |
| Complete bipartite K(p,q) | p, q | nonzero eigenvalues ±√(pq); E = 2√(pq) |
| Broom B(r,q) | r, q | hub maximises local energy within the graph |
| Double star S(p,q) | p, q | tree of diameter 3 |
| Symmetric double broom DB(t,k) | t, k | symmetric under the swap of its two hubs |
| Petersen graph | — | spectrum 3, 1⁵, (−2)⁴; E = 16 |

## Verified numerics

The eigenvalue engine is the **classical Jacobi rotation algorithm** for real symmetric matrices, implemented from scratch in ~40 lines of JavaScript — provably convergent and highly accurate for the matrix sizes the Atlas handles. It is tested three ways, and the full suite runs in CI on every push:

1. **Closed forms** (`tests/verify_engine.js`, `tests/verify_families.js`) — energies of P_n, C_n, S_n, K_n, K_{p,q}, the Petersen graph (E = 16) and the broom B(4,3) (E = 2√(6+2√7)) are checked against their exact values, along with the structural bounds 0 ≤ E_G(v) ≤ 2√(deg v) and e(G) ≤ 2E(G) on random graphs.
2. **Cross-check against NumPy** (`tests/cross_check_numpy.py`) — the engine's output on 25 random graphs is compared with `numpy.linalg.eigvalsh`; the maximum discrepancy observed is ~1e-14.
3. **Curve sweeps** (`tests/verify_curves.js`) — points of the Fig. 2 family curves are validated against known values and internal-consistency checks.

```bash
npm test              # JavaScript suites (no dependencies)
npm run test:numpy    # cross-check against NumPy (requires python3 + numpy)
```

## Run locally & deploy

It is a single file — open `index.html` in a browser and that's it. To publish your own copy on GitHub Pages:

```bash
git init && git add . && git commit -m "Atlas of local vertex energy"
git branch -M main
git remote add origin https://github.com/AndresGomez87/atlas-local-vertex-energy.git
git push -u origin main
```

Then on GitHub: **Settings → Pages → Source: Deploy from a branch → main / (root)**. The site appears at `https://AndresGomez87.github.io/atlas-local-vertex-energy/` within a minute or two.

## References

- C. Espinal, J. Rada, *Graph energy change on vertex deletion* (2024) — definition of local vertex energy and the bound E_G(v) ≤ 2√(deg v).
- I. Gutman, *The energy of a graph* (1978) — the classical graph energy.

## Author

Developed by **Andres Gomez**, mathematics graduate of Universidad de Antioquia (Medellín, Colombia), as a companion to an undergraduate thesis on local vertex energy in structured families of trees.

Licensed under the [MIT License](LICENSE).
