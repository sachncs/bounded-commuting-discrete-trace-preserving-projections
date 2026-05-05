# Research Paper Review: BCDTPP Implementation vs. Ern–Guzmán–Potu (2026)

**Paper:** "Bounded, Commuting, Discrete-trace Preserving Projections" (arXiv:2604.28103v1)  
**Authors:** Alexandre Ern, Johnny Guzmán, Pratyush Potu  
**Repository:** `/Users/sachin/Research/math/bcdtpp`  
**Review Date:** 2026-05-05

> **Note on evidence limits:** The full paper text was not extractable from the PDF (compressed binary streams only). Therefore, findings requiring specific theorem or algorithm verification are labeled **suspected** or **needs paper verification**. All other findings are based on direct code evidence.

---

## 1. Fidelity Verdict

**Partially faithful.**

The implementation faithfully reproduces the *lowest-order* boundary-aware projection operators for the 3D de Rham complex (`l = 0,1,2,3`, `p = 0`). The Alfeld and Worsey-Farin splits, local surface stiffness solves, and trace-preserving boundary DoFs are all present and well-structured. However, the code does not implement the full decomposition `Pi^l = Pi_partial^l + Pi_ring^l (I - Pi_partial^l)` advertised in the README — `Pi_ring^l` is entirely missing. Higher-order projections (`p > 0`) are only supported for `l = 0` (scalar bubble) and `l = 3` (monomial L2); `l = 1` and `l = 2` higher-order are explicitly unimplemented. No numerical experiment pipeline, convergence verification, or stability-norm computation from the paper is present.

---

## 2. Executive Summary

The repository is a clean, well-tested JavaScript implementation (~2,700 LOC) of the *novel* contribution from Ern–Guzmán–Potu (2026): the boundary-correction operators `Pi_partial^l` using Alfeld/Worsey-Farin subdivided patches. Topology, Whitney forms, quadrature, and local solvers are all implemented natively. Recent refactoring improved encapsulation and safety.

**Risks:**
- The full projection formula from the paper is not implemented (`Pi_ring^l` missing).
- Higher-order vector-valued projections (`H(curl)` and `H(div)` for `p > 0`) are absent.
- No experimental validation or convergence study is included.
- Several algorithmic shortcuts (simplex row replacement for constrained solve, midpoint-split fallback in AABB tree) deviate from best-practice numerics.

**Bottom line:** The code is a solid *reference implementation* of the paper's boundary-correction machinery, but it is not a complete reproduction of the full theoretical framework or experimental section.

---

## 3. Paper-to-Code Mapping Table

| Paper Component (inferred from code/comments) | Implementation Location | Status |
|---|---|---|
| Tetrahedral mesh topology & boundary extraction | `src/lib/mesh.js` | **Fully implemented** |
| Alfeld split of boundary faces (Sec 6.1.3) | `src/lib/mesh_refinement.js:computeAlfeldSplit` | **Fully implemented** |
| Worsey-Farin split of tetrahedra (Sec 6.1.4) | `src/lib/mesh_refinement.js:computeWorseyFarinSplit` | **Fully implemented** |
| Lowest-order vertex weights (Sec 6.3.1) | `src/lib/boundary_weight_computer.js:#computeVertexWeights` | **Fully implemented** |
| `H^1` projection (`l = 0`, `p = 0`) | `src/lib/bcdtpp.js:projectH1` | **Fully implemented** |
| `H(curl)` projection (`l = 1`, `p = 0`) | `src/lib/bcdtpp.js:projectHcurl` | **Fully implemented** |
| `H(div)` projection (`l = 2`, `p = 0`) | `src/lib/bcdtpp.js:projectHdiv` | **Fully implemented** |
| `L^2` projection (`l = 3`, `p = 0`) | `src/lib/bcdtpp.js:projectL2` | **Fully implemented** |
| Higher-order scalar bubble basis (Sec 7) | `src/lib/higher_order_projection.js` | **Partial** (only `l=0`, `p>=4`) |
| Higher-order L2 monomial basis | `src/lib/higher_order_projection.js:solveL2Projection` | **Partial** (only `l=3`) |
| `Pi_ring^l` operator | — | **Missing** |
| Discrete extension / stable lifting operators | — | **Missing** |
| Numerical experiments / convergence rates | — | **Missing** |
| Graph-norm stability verification | — | **Missing** |
| Commuting-diagram property verification | `tests/BCDTPP.test.js` | **Partial** (smoke tests only) |

---

## 4. Critical Gaps and Deviations

### G1: `Pi_ring^l` is not implemented
- **Category:** missing
- **Location:** README defines `Pi^l = Pi_partial^l + Pi_ring^l (I - Pi_partial^l)`; no code file implements `Pi_ring^l`.
- **Impact:** The code computes what appears to be the *full* projection directly (boundary corrections + standard interior interpolation), but the README claims only `Pi_partial^l` is implemented. The relationship between the code and the paper's operator decomposition is **unclear**.
- **Fix:** Either implement `Pi_ring^l` as a separate operator and compose it in `projectHp`, or update the README to clarify that `projectH1`/`projectHcurl`/etc. already represent the complete operator and `Pi_ring^l` is not needed for the implemented scope.

### G2: Higher-order `H(curl)` and `H(div)` missing
- **Category:** missing
- **Location:** `src/lib/bcdtpp.js:projectHp` lines 296-299.
- **Impact:** For `l = 1` or `l = 2` with `p > 0`, the code throws. The paper likely covers full polynomial de Rham sequences.
- **Fix:** Implement Nédélec and Raviart-Thomas higher-order basis extensions.

### G3: `p = 1,2,3` bubble space is empty
- **Category:** partial
- **Location:** `src/lib/higher_order_projection.js:#getBubbleExponents`.
- **Impact:** For `l = 0` and `p < 4`, `projectHp` returns the lowest-order base without any enrichment. The bubble space `b * P^{p-4}` is only non-empty for `p >= 4`. This means degrees 1, 2, 3 are not actually higher-order.
- **Fix:** Add hierarchical enrichment (e.g., edge/face/bubble DoFs) for `p = 1,2,3`.

### G4: No numerical experiment pipeline
- **Category:** missing
- **Location:** Entire repo.
- **Impact:** The paper's results cannot be reproduced. No convergence plots, no stability norm measurements, no discrete extension tests.
- **Fix:** Add a batch evaluation harness with manufactured solutions and mesh refinement loops.

### G5: Local solver uses textbook row-replacement
- **Category:** divergent
- **Location:** `src/lib/local_solver.js:solveWithConstraint` lines 75-96.
- **Impact:** The constraint `sum(x) = 0` is enforced by overwriting the last row of `K` with ones. The JSDoc itself admits this "destroys symmetry and can degrade conditioning." This is a known numerics shortcut.
- **Fix:** Use a Lagrange-multiplier formulation or projected-gradient solve for larger patches.

### G6: `HigherOrderProjection.evaluateBubble` ignores `tIdx`
- **Category:** suspected
- **Location:** `src/lib/higher_order_projection.js:evaluateBubble` line 170.
- **Impact:** The method signature accepts `tIdx` but only uses it to get barycentric coordinates via `this.whitney.getBarycentric`. The bubble basis is purely barycentric, so this is mathematically correct, but the parameter is misleading.
- **Fix:** Remove unused `tIdx` from signature or document why it is retained.

---

## 5. Algorithm Correctness Risks

| Risk | Location | Severity | Evidence |
|---|---|---|---|
| **Constrained solve fallback** | `local_solver.js:92-97` | Medium | Falls back to uniform zero on failure. Mathematically satisfies the constraint, but may silently degrade accuracy without user awareness. |
| **AABB tree balance** | `point_locator.js:42-53` | Low | Fixed in recent refactoring: now uses median-split. Previously midpoint-split could degrade to O(N). |
| **Face permutation parity** | `mesh.js:303-309` | Low | Fixed: now uses inversion count instead of product-of-differences. |
| **Type inconsistency in projections** | `bcdtpp.js:137-238` | Medium | Fixed: scalar/vector paths are now determined at evaluation points, not just at the barycenter. |
| **0^0 in polynomial basis** | `higher_order_projection.js:54-60` | Low | Fixed: explicit `#pow` helper documents the convention. |
| **luSolve scaling** | `math_utils.js:109-163` | Low | Fixed: row equilibration + relative tolerance added. |
| **Missing tIdx bounds checks** | `bcdtpp.js:111-270` | Low | Fixed: `#validateTetIdx` guards all projection entry points. |

---

## 6. Reproducibility Assessment

**What works:**
- Unit tests verify exactness for constants, linear polynomials, and basic commuting properties on a single tetrahedron.
- Boundary-weight computation runs without errors on the test mesh.
- Point location and safe expression evaluation are tested.

**What is missing:**
- No mesh refinement studies (uniform refinement, adaptive refinement).
- No convergence rate measurement (L2, H1, H(curl), H(div) norms).
- No comparison against standard interpolation operators.
- No stability constant estimation (graph norm bounds).
- No discrete extension experiment (the paper's "optimality result").
- No batch processing harness or CLI.

**What blocks replication:**
- The absence of `Pi_ring^l` and higher-order vector projections means anyone trying to reproduce the full theoretical construction will find the implementation incomplete.
- No datasets or mesh generators are provided; users must supply their own tetrahedral meshes in a specific JSON format.

---

## 7. Dead Code and Unused Concepts

| Item | Status | Notes |
|---|---|---|
| `Pi_ring^l` | **Confirmed missing concept** | Referenced in README, absent in code. |
| `safe_eval.js` recursion depth limits | **Recently added** | Previously unlimited; now bounded. |
| `yargs` dependency | **Removed** | Was unused. |
| `getBubbleBasis` in `whitney.js` | **Removed** | Was dead code. |
| `integrateTriangle` | **Restored after accidental deletion** | Used by tests. |
| `Mesh.computeAlfeldSplit` / `computeWorseyFarinSplit` | **Moved to `MeshRefinement`** | Old methods removed from `Mesh`. |

---

## 8. Architecture Assessment

**Strengths:**
- **Good separation of concerns:** `Mesh` (data), `MeshRefinement` (mutation), `BoundaryWeightComputer` (local solves), `Whitney` (basis), `Bcdtpp` (projections), `PointLocator` (geometry query).
- **Private fields:** Recent refactoring made `edgeKeyToIndex`, `tetEdgeSigns`, `tetFaceSigns`, and boundary-weight Maps private.
- **Injectable warnings:** `onWarning` option allows library consumers to intercept warnings.
- **Pure JavaScript:** No native dependencies; easy to audit.

**Weaknesses:**
- **Tight coupling:** `Bcdtpp` directly accesses `this.mesh.vertices`, `this.mesh.faces`, etc. Documented in JSDoc, but still means the class cannot be used with a different mesh implementation without extensive rewriting.
- **Mesh mutability:** `MeshRefinement` mutates `mesh.vertices` and `mesh.vertexCount` in-place. This makes the mesh object stateful and harder to reason about after refinement.
- **UI coupling:** `ui.js` hardcodes DOM IDs. Acceptable for a demo, but documented as non-reusable.
- **No plugin architecture:** Higher-order support is hardcoded; adding new polynomial spaces requires editing `HigherOrderProjection`.

---

## 9. Improvements (Prioritized)

### P0 — Critical for correctness
1. **Clarify or implement `Pi_ring^l`.** The README claims the library implements `Pi_partial^l` only, but the user-facing methods appear to return full projections. Reconcile the README with the code, or implement the missing operator.
2. **Higher-order `H(curl)` and `H(div)`.** Implement Nédélec and Raviart-Thomas `p > 0` basis functions and projection solves.

### P1 — Important for fidelity
3. **Convergence experiment harness.** Add a script that runs projections on a sequence of refined meshes and computes L2/H1/H(curl)/H(div) errors against manufactured solutions.
4. **Replace row-replacement constraint solver.** Use a symmetric constrained solver (e.g., augmented Lagrangian or null-space method) in `LocalSolver.solveWithConstraint`.
5. **Bubble space for `p = 1,2,3`.** The current code returns the base projection for `p < 4`. Add hierarchical DoFs (edge bubbles, face bubbles) for low polynomial degrees.

### P2 — Quality of life
6. **CLI / batch mode.** Add a Node.js CLI that accepts mesh files and expression strings, outputs projections to JSON.
7. **Structured error types.** Replace generic `Error` throws with specific error classes (`MeshValidationError`, `ProjectionError`, `SingularMatrixError`) for better observability.
8. **Metrics hooks.** Add optional timing callbacks for `computeBoundaryWeights`, `solveBubbleProjection`, and `findTetrahedron`.

---

## 10. Feature Additions

| Feature | Research Alignment | Practical Value |
|---|---|---|
| **Convergence plotting** | High (paper has convergence results) | Essential for validation |
| **Adaptive mesh refinement loop** | Medium | Research extension |
| **JSON schema for mesh I/O** | Low | Production robustness |
| **WebWorker batch evaluator** | Low | UI responsiveness |
| **Export to VTK/Paraview** | Medium | Visualization of projected fields |
| **Stability constant estimator** | High (paper's graph-norm stability) | Theoretical verification |
| **Discrete extension operator** | High (paper's second application) | Practical for domain decomposition |

---

## 11. Final Recommendation

**Minimum changes required to achieve fidelity:**
1. **Document the scope gap.** Update the README to explicitly state that `Pi_ring^l` is not implemented and that the code provides the boundary-correction operators plus standard interior interpolation.
2. **Implement higher-order vector projections** (`l = 1, 2`, `p > 0`). Without these, the de Rham complex is incomplete for polynomial degrees above 0.
3. **Add a convergence experiment.** Even a single script with one manufactured solution on a refined mesh would validate the implementation against the paper's claims.

**Steps to reach production readiness:**
1. Add structured logging and error types.
2. Implement a CLI for batch processing.
3. Add VTK export for field visualization.
4. Profile and optimize `Bcdtpp.#computeBoundaryIntegralH1` for large meshes (the current O(F·A) filtering is now improved but still not cached across vertices).
5. Add input mesh generators (unit cube, sphere) so users can run experiments without external meshing tools.

**Overall:** The repository is a **solid, well-tested reference implementation** of the paper's *boundary-correction machinery*. It is suitable for educational use, algorithm verification on simple meshes, and as a foundation for a fuller implementation. It is **not yet a complete reproduction** of the paper's theoretical framework or experimental results.
