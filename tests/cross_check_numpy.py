"""Cross-check: the Atlas JavaScript eigenvalue engine vs NumPy.

Run `node tests/verify_engine.js` first (it exports tests/js_results.json
with 25 random graphs and the energies computed by the JS engine), then:

    python3 tests/cross_check_numpy.py

Passes if the maximum absolute discrepancy is below 1e-9 (in practice
it is ~1e-14).
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "js_results.json")

if not os.path.exists(RESULTS):
    sys.exit("tests/js_results.json not found - run `node tests/verify_engine.js` first.")

with open(RESULTS) as f:
    cases = json.load(f)


def energy(A: np.ndarray) -> float:
    return float(np.sum(np.abs(np.linalg.eigvalsh(A))))


def local_energies(A: np.ndarray) -> list[float]:
    n = A.shape[0]
    EG = energy(A)
    out = []
    for k in range(n):
        idx = [i for i in range(n) if i != k]
        out.append(EG - energy(A[np.ix_(idx, idx)]))
    return out


max_err_E = 0.0
max_err_le = 0.0
for c in cases:
    A = np.array(c["A"], dtype=float)
    max_err_E = max(max_err_E, abs(energy(A) - c["E"]))
    for a, b in zip(local_energies(A), c["le"]):
        max_err_le = max(max_err_le, abs(a - b))

print(f"{len(cases)} random graphs, JS engine vs numpy.linalg.eigvalsh:")
print(f"  max |Delta E(G)|   = {max_err_E:.2e}")
print(f"  max |Delta E_G(v)| = {max_err_le:.2e}")

if max(max_err_E, max_err_le) < 1e-9:
    print("CROSS-CHECK OK")
else:
    sys.exit("CROSS-CHECK FAILED")
