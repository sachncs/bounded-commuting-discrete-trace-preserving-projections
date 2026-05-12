# Bcdtpp.js

[![CI](https://github.com/sachin/bcdtpp/actions/workflows/ci.yml/badge.svg)](https://github.com/sachin/bcdtpp/actions/workflows/ci.yml)
[![JavaScript Style Guide](https://img.shields.io/badge/code_style-standard-brightgreen.svg)](https://standardjs.com)

JavaScript implementation of **Bounded, Commuting, Discrete-trace Preserving Projections** for the 3D de Rham complex on simplicial meshes.

Based on the paper: [*Ern, Guzmán, Potu (2026) arXiv:2604.28103v1*](https://arxiv.org/abs/2604.28103).

> **Disclaimer:** I am not an author of the paper above. This repository is an independent JavaScript implementation of the algorithm described in that work.

## Features

- **Pure JavaScript**: No external math dependencies; all linear algebra implemented natively.
- **Simplicial Mesh Support**: 3D tetrahedral mesh processing with full topology (faces, edges, incidence).
- **Geometric Splits**:
  - **Alfeld Split**: Boundary face subdivision for local solvers.
  - **Worsey-Farin Split**: Bulk tetrahedron subdivision for weight extensions.
- **Full de Rham Complex (Lowest-Order)**:
  - `l = 0` (`H^1`): Vertex-based projections using surface Poisson solvers.
  - `l = 1` (`H(curl)`): Edge-based projections using edge Whitney forms.
  - `l = 2` (`H(div)`): Face-based projections using face Whitney forms.
  - `l = 3` (`L^2`): Cell-based projections via cell averages.
- **Higher-Order Projections**: Scalar bubble basis support for `p >= 4` on `H^1`; L2 monomial basis for `p >= 1` on `l = 3`. Vector-valued higher-order (`l = 1, 2`, `p > 0`) is not yet implemented.
- **Commuting Properties**: Designed to commute with exterior derivatives while preserving discrete traces.
- **Fast Point Location**: AABB tree for `O(log N)` point-in-tetrahedron queries.

## Installation

```bash
npm install bcdtpp
```

### From source

```bash
git clone https://github.com/sachin/bcdtpp.git
cd bcdtpp
npm install
```

## Usage

### 1. Define a Mesh

```javascript
import { Mesh } from 'bcdtpp'

const mesh = new Mesh(vertices, tetrahedra)
```

### 2. Initialize Projections

```javascript
import { Whitney, Bcdtpp } from 'bcdtpp'

const whitney = new Whitney(mesh)
const bcdtpp = new Bcdtpp(mesh, whitney, { quadratureOrder: 3 })

// Compute boundary weights (Alfeld splits + local solves)
bcdtpp.computeBoundaryWeights()
bcdtpp.buildPointLocator()
```

### 3. Project a Function

```javascript
const u = (p) => Math.sin(p[0])
const val = bcdtpp.projectH1(u, [0.5, 0.5, 0.5], 0)
```

### 4. Higher-Order Projection

```javascript
const val = bcdtpp.projectHp(u, point, tetIdx, /* l */ 0, /* p */ 2)
```

### 5. Point Location

```javascript
const result = bcdtpp.projectAtPoint(u, [0.1, 0.2, 0.3], /* l */ 0, /* p */ 0)
// result = { value, tIdx, bary }
```

## Development

This project uses a zero-config toolchain. Every task is a single `npm run` command away.

| Script | Purpose |
|--------|---------|
| `npm run lint` | Check StandardJS style |
| `npm run lint:fix` | Auto-fix StandardJS violations |
| `npm test` | Run the full Mocha test suite |
| `npm run test:watch` | Run tests in watch mode |
| `npm run test:coverage` | Run tests with c8 coverage |
| `npm run build` | Build ESM, CJS, and UMD bundles |
| `npm run docs` | Regenerate API docs from JSDoc |

## Project Structure

```
src/
  lib/
    bcdtpp.js                 — Main projection class (Bcdtpp)
    mesh.js                   — Tetrahedral mesh topology & splits
    whitney.js                — Whitney forms and barycentric utilities
    quadrature.js             — Gaussian quadrature for triangles & tetrahedra
    math_utils.js             — Pure JS linear algebra (LU, inverse, vector ops)
    local_solver.js           — Boundary-patch stiffness assembly & solves
    point_locator.js          — AABB tree point-in-tet locator
    higher_order_projection.js — Bubble basis for p >= 1
    boundary_weight_computer.js — Boundary weight computation
    mesh_refinement.js        — Alfeld and Worsey-Farin mesh refinement
    errors.js                 — Custom error classes
    projectors/
      h1_projector.js         — Pi^0 (vertex-based)
      hcurl_projector.js      — Pi^1 (edge-based)
      hdiv_projector.js       — Pi^2 (face-based)
      l2_projector.js         — Pi^3 (cell-based)

tests/
  *.test.js                  — Mocha + Chai test suites

docs/
  api.md                     — Generated API reference
  architecture.md            — Module map and data flow
  math.md                    — Mathematical background
  exceptions.md              — Error taxonomy and recovery
  setup.md                   — Installation and usage
  development.md             — Development workflow
  testing.md                 — Testing philosophy and guide
  publishing.md              — Versioning and release process
```

## Mathematical Implementation Details

This library implements the construction of boundary correction operators `Pi_partial^l`.
The final projection is defined as:

```
Pi^l = Pi_partial^l + Pi_ring^l (I - Pi_partial^l)
```

The novel contribution of the paper is the construction of `Pi_partial^l` using local problems on subdivided patches to ensure stability and trace-preservation.

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct and the process for submitting pull requests.

## License

[MIT](LICENSE)
