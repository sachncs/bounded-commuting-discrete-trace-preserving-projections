# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

### Added

- Add Tech Stack section to README
- Add Security section to README
- Add Code of Conduct section to README
- Add centered header with aligned badges to README
- Add license copyright to README

### Changed

- Restructure README header with centered layout and aligned badges

## [0.0.7](https://github.com/sachncs/bounded-commuting-discrete-trace-preserving-projections/commit/e9b643c) — 2026-06-20

- Add comprehensive test suite: interior projector, convergence harness, local solver, mesh generator, point locator, Whitney basis, quadrature, coverage edge cases
- Add convergence test framework (h-refinement and p-refinement on all four form degrees)
- Add mesh generator for cube subdivision
- Add PointLocator AABB tree for O(log N) point-in-tet queries
- Add LocalSolver with surface stiffness assembly and constrained solve
- Add BoundaryWeightComputer for trace-preserving projection weights
- Implement Worsey-Farin split (§6.1.4) with idempotency
- Implement Alfeld split (§6.1.3) with idempotency
- Fix UMD global name to Bcdtpp in rollup config
- Fix linting issues across all source files

## [0.0.6](https://github.com/sachncs/bounded-commuting-discrete-trace-preserving-projections/commit/403d048) — 2026-05-12

- Add Mesh class with topology (faces, edges, boundary), geometry (volume, normals), and adjacency
- Add Whitney class for barycentric coordinates and Whitney finite-element basis functions
- Add math utilities: LU solver, 3x3 inverse, numerical gradient, vector operations
- Add higher-order projection framework (§7) with bubble corrections
- Add L2 projector (l=3) for piecewise constants
- Add MeshRefinement class for Alfeld/Worsey-Farin splits
- Fix orientation sign computation for edges and faces
- Fix tetDeterminant (scalar triple product) and add tetVolume
- Refactor projectors into separate files under src/lib/projectors/

## [0.0.5](https://github.com/sachncs/bounded-commuting-discrete-trace-preserving-projections/commit/8daf318) — 2026-05-12

- Add H(curl) projector (l=1) for Nédélec first-kind (Whitney 1-form) space
- Add H(div) projector (l=2) for Raviart-Thomas (Whitney 2-form) space
- Implement interior projector Π_ring^l and discrete extension operator E^l
- Implement boundary correction part Π_partial^l
- Implement global projector decomposition: Π^l = Π_partial^l + Π_ring^l (I − Π_partial^l)
- Add computeBoundaryWeights for trace-preserving projections

## [0.0.4](https://github.com/sachncs/bounded-commuting-discrete-trace-preserving-projections/commit/399c333) — 2026-05-05

- Add H1 projector (l=0) for continuous piecewise-linear space
- Implement boundary-aware projection with vertex weights
- Add projectAtPoint convenience method with AABB point locator
- Add extractBoundaryDofs for boundary degree-of-freedom extraction

## [0.0.3](https://github.com/sachncs/bounded-commuting-discrete-trace-preserving-projections/commit/593e4d2) — 2026-05-05

- Add Bcdtpp main class with project method dispatching to form degree l
- Add basic error classes: MeshValidationError, ProjectionError, SingularMatrixError
- Implement first version of H1 projection
- Add rollup UMD bundle configuration

## [0.0.2](https://github.com/sachncs/bounded-commuting-discrete-trace-preserving-projections/commit/cc2313d) — 2026-05-04

- Add Mesh class with tetrahedral mesh construction and validation
- Add Whitney class with barycentric coordinates and edge/face basis functions
- Add math utilities: dot, cross, norm, subtract, matrix operations
- Add Gaussian quadrature for triangles, tetrahedra, and lines
- Add Mocha/Chai test infrastructure

## [0.0.1](https://github.com/sachncs/bounded-commuting-discrete-trace-preserving-projections/commit/197b86b) — 2026-05-03

- Initial project setup with npm, StandardJS, Mocha/Chai, c8 coverage
- Add README, CONTRIBUTING, LICENSE, and project configuration
- Add GitHub Actions CI workflow
- Add documentation: architecture, math background, setup, testing, publishing

---

## Dependency Updates

### 2026-07-07

- Merge pull request #1 — Bump actions/checkout from 4 to 7
- Merge pull request #2 — Bump actions/setup-node from 4 to 6
- Merge pull request #3 — Bump dev-dependencies group (4 updates)
- Merge pull request #4 — Bump chai from 5.3.3 to 6.2.2
- Merge pull request #5 — Bump c8 from 10.1.3 to 11.0.0
- Merge pull request #6 — Bump @babel/preset-env from 7.29.5 to 8.0.2
- Merge pull request #7 — Bump @babel/core from 7.29.0 to 8.0.1
