<p align="center">
  <h1 align="center">bcdtpp</h1>
  <p align="center">Bounded, Commuting, Discrete-trace Preserving Projections for the 3D de Rham complex on simplicial meshes.</p>
  <p align="center">
    <a href="#installation"><img src="https://img.shields.io/npm/v/bcdtpp.svg" alt="npm version"></a>
    <a href="LICENSE"><img src="https://img.shields.io/badge/License-MIT-green" alt="License"></a>
    <a href="https://github.com/sachncs/bounded-commuting-discrete-trace-preserving-projections/actions"><img src="https://img.shields.io/github/actions/workflow/status/sachncs/bounded-commuting-discrete-trace-preserving-projections/ci.yml?branch=master" alt="CI"></a>
    <a href="https://www.npmjs.com/package/bcdtpp"><img src="https://img.shields.io/npm/dm/bcdtpp" alt="Downloads"></a>
    <a href="https://standardjs.com"><img src="https://img.shields.io/badge/code_style-standard-brightgreen" alt="StandardJS"></a>
  </p>
</p>

Based on the paper: [*Ern, Guzmán, Potu (2026) arXiv:2604.28103v1*](https://arxiv.org/abs/2604.28103).

> **Disclaimer:** I am not an author of the paper above. This repository is an independent JavaScript implementation of the algorithm described in that work.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
- [API Overview](#api-overview)
- [Error Handling](#error-handling)
- [TypeScript](#typescript)
- [Browser Usage](#browser-usage)
- [Performance](#performance)
- [Development](#development)
- [Project Structure](#project-structure)
- [Mathematical Background](#mathematical-background)
- [What's Next](#whats-next)
- [Changelog](#changelog)
- [Contributing](#contributing)
- [License](#license)

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
git clone https://github.com/sachncs/bounded-commuting-discrete-trace-preserving-projections.git
cd bounded-commuting-discrete-trace-preserving-projections
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

## API Overview

The package exports the following named entities:

| Export | Subpath | Description |
|--------|---------|-------------|
| `Bcdtpp` | `bcdtpp` | Main projection orchestrator |
| `Mesh` | `bcdtpp/mesh` | Tetrahedral mesh topology and geometry |
| `Whitney` | `bcdtpp/whitney` | Barycentric coordinates and Whitney basis |
| `PointLocator` | `bcdtpp/point-locator` | AABB tree point-in-tet locator |
| `quadrature` | `bcdtpp/quadrature` | Gaussian quadrature utilities |
| `mathUtils` | `bcdtpp/math-utils` | Linear algebra primitives |

All subpaths are exposed via the `exports` field in `package.json` and include hand-written `.d.ts` files.

## Error Handling

The library throws specific error subclasses so callers can distinguish mesh problems from projection problems:

- `MeshValidationError` — Invalid mesh data (degenerate tets, bad indices, non-finite vertices).
- `ProjectionError` — Invalid arguments to projection methods or unsupported configurations.
- `SingularMatrixError` — Numerical linear algebra failure (singular Jacobian, zero pivot).

Recoverable boundary-weight failures emit warnings rather than throwing, so a single bad element does not halt the entire mesh projection. See [docs/exceptions.md](docs/exceptions.md) for the full taxonomy.

## TypeScript

Hand-written `.d.ts` declaration files are co-located with every source module in `src/lib/`. If your bundler does not resolve them automatically, reference the package in your `tsconfig.json`:

```json
{
  "compilerOptions": {
    "types": ["bcdtpp"]
  }
}
```

## Browser Usage

The UMD bundle is available via CDN:

```html
<script src="https://cdn.jsdelivr.net/npm/bcdtpp/dist/bcdtpp.umd.js"></script>
<script>
  const { Mesh, Whitney, Bcdtpp } = window.Bcdtpp
</script>
```

## Performance

- **Mesh construction**: `O(V + T)` where `V` is vertices and `T` is tetrahedra.
- **Point location**: `O(log T)` per query after AABB tree construction.
- **Boundary weight computation**: `O(V * k)` where `k` is the typical boundary-patch size; this is the dominant upfront cost.
- **Per-element projection**: `O(1)` for lowest-order; `O(p^3)` for higher-order scalar enrichment.

## Tech Stack

| Category | Technology |
|----------|------------|
| Language | JavaScript (ES2022) |
| Build | [Rollup](https://rollupjs.org/) (ESM + CJS + UMD) |
| Lint | [StandardJS](https://standardjs.com/) |
| Testing | [Mocha](https://mochajs.org/) + [Chai](https://www.chaijs.com/) |
| Coverage | [c8](https://github.com/bcoe/c8) |
| CI | [GitHub Actions](https://github.com/features/actions) |

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

## Mathematical Background

This library implements the construction of boundary correction operators `Pi_partial^l`.
The final projection is defined as:

```
Pi^l = Pi_partial^l + Pi_ring^l (I - Pi_partial^l)
```

The novel contribution of the paper is the construction of `Pi_partial^l` using local problems on subdivided patches to ensure stability and trace-preservation.

See [docs/math.md](docs/math.md) for the full mathematical exposition (de Rham complex, Whitney forms, commuting diagrams, and mesh refinement theory).

## What's Next

The following items are planned but not yet implemented. Contributions are welcome.

### High Priority

- **Vector-valued higher-order projections** (`l = 1, 2`, `p > 0`): Nédélec and Raviart-Thomas enrichment for `H(curl)` and `H(div)`.
- **Adaptive mesh refinement support**: Integrate `MeshRefinement` APIs into `Bcdtpp` so boundary weights can be recomputed incrementally as the mesh refines.

### Medium Priority

- **Web Worker parallelization**: Offload per-tet projection and boundary weight solves to workers for large meshes.
- **WASM acceleration**: Port the dense linear algebra routines (`luSolve`, `inverse3x3`) to WebAssembly for a 2-5x speedup.

### Low Priority / Research

- **Anisotropic mesh support**: Generalize the point locator and quadrature to handle highly stretched tets without loss of precision.
- **Time-dependent projections**: Cache-friendly APIs for projecting fields that evolve between time steps.
- **Arbitrary polynomial degree `p`**: Unify the scalar bubble and L2 monomial bases into a single hierarchical basis generator.

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for a detailed history of releases and breaking changes.

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct and the process for submitting pull requests.

## Code of Conduct

This project follows the [Contributor Covenant v2.1](CODE_OF_CONDUCT.md).
By participating you agree to abide by its terms.

## Security

Report vulnerabilities to **sachncs@gmail.com** — see [SECURITY.md](SECURITY.md).

## License

[MIT](LICENSE) © 2026 Sachin
