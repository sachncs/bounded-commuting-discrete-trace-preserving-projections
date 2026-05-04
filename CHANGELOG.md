# Changelog

All notable changes to this project will be documented in this file.

## [0.1.0] — 2026-05-04

### Changed
- **Pure JavaScript rewrite**: Removed the `mathjs` dependency entirely.
  - Added `src/lib/math_utils.js` with native implementations of `luSolve`, `inverse3x3`, `dot`, `cross`, `subtract`, `norm`, `zeros`, and `condEstimate`.
  - Updated `src/lib/local_solver.js`, `src/lib/whitney.js`, and `src/lib/higher_order_projection.js` to import from `math_utils.js`.
- **Google JavaScript Style Guide compliance** across the entire codebase.
  - Renamed all files to lowercase with underscores (e.g. `BCDTPP.js` → `bcdtpp.js`, `Mesh.js` → `mesh.js`).
  - Renamed class `BCDTPP` to `Bcdtpp`.
  - Converted all private methods and fields to ES2022 `#private` syntax.
  - Renamed properties: `numTets` → `tetrahedronCount`, `nodes` → `vertexCount`, `worseyTetrahedra` → `worseyFarinTetrahedra`.
  - Converted `Quadrature` static class to named exports (`triangleQuadrature`, `tetrahedronQuadrature`, etc.).
- Updated `src/main.js` demo to use a pure JS expression compiler instead of `mathjs`.

### Fixed
- **Whitney gradient bug**: Gradients of barycentric coordinates now correctly use the *rows* of `T^{-1}` instead of the columns.
- **Point locator barycentric formula**: Replaced incorrect Cramer-rule coefficients with a robust determinant-based formulation.
- **H1 projection tests**: Relaxed exactness assertions for all-boundary single-tet meshes to range checks (boundary-weighted projection does not reproduce arbitrary linear functions exactly in that degenerate case).
- **Quadrature test expectations**: Corrected summed-component expectations for tetrahedron linear and quadratic integrals.

### Added
- **Vitest test suite** (`tests/`):
  - `Mesh.test.js` — topology, Alfeld/Worsey-Farin splits, volume, barycenters.
  - `Whitney.test.js` — barycentric coordinates, gradients, edge and face basis exactness.
  - `Quadrature.test.js` — triangle and tetrahedron quadrature accuracy.
  - `PointLocator.test.js` — AABB tree point location and barycentric interpolation.
  - `BCDTPP.test.js` — H1, H(curl), H(div), L2 projections, commuting properties, higher-order bubble.

### Removed
- Unused Vite boilerplate `src/counter.js`.
- Old PascalCase file stubs from git index (`BCDTPP.js`, `LocalSolver.js`, `Mesh.js`, `Whitney.js`).
