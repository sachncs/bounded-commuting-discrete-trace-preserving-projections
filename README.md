# Bcdtpp.js

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
- **Higher-Order Projections**: Scalar bubble basis support for `p >= 1` on `H^1`.
- **Commuting Properties**: Designed to commute with exterior derivatives while preserving discrete traces.
- **Fast Point Location**: AABB tree for `O(log N)` point-in-tetrahedron queries.

## Installation

```bash
npm install
```

## Usage

### 1. Define a Mesh

```javascript
import {Mesh} from './src/lib/mesh.js';

const mesh = new Mesh(vertices, tetrahedra);
```

### 2. Initialize Projections

```javascript
import {Whitney} from './src/lib/whitney.js';
import {Bcdtpp} from './src/lib/bcdtpp.js';

const whitney = new Whitney(mesh);
const bcdtpp = new Bcdtpp(mesh, whitney, {quadratureOrder: 3});

// Compute boundary weights (Alfeld splits + local solves)
await bcdtpp.computeBoundaryWeights();
bcdtpp.buildPointLocator();
```

### 3. Project a Function

```javascript
const u = (p) => Math.sin(p[0]);
const val = bcdtpp.projectH1(u, [0.5, 0.5, 0.5], 0);
```

### 4. Higher-Order Projection

```javascript
const val = bcdtpp.projectHp(u, point, tetIdx, /* l */ 0, /* p */ 2);
```

### 5. Point Location

```javascript
const result = bcdtpp.projectAtPoint(u, [0.1, 0.2, 0.3], /* l */ 0, /* p */ 0);
// result = {value, tIdx, bary}
```

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
  main.js                     — Interactive demo (Vite + Three.js)

tests/
  *.test.js                  — Vitest suite (38 tests)
```

## Mathematical Implementation Details

This library implements the construction of boundary correction operators `Pi_partial^l`.
The final projection is defined as:

```
Pi^l = Pi_partial^l + Pi_ring^l (I - Pi_partial^l)
```

The novel contribution of the paper is the construction of `Pi_partial^l` using local problems on subdivided patches to ensure stability and trace-preservation.

## Style & Conventions

This codebase follows the [Google JavaScript Style Guide](https://google.github.io/styleguide/jsguide.html):

- `camelCase` for variables, functions, and methods
- `PascalCase` for classes
- `CONSTANT_CASE` for constants
- ES2022 private class fields and methods (`#field`, `#method()`)
- Lowercase filenames with underscores (e.g., `math_utils.js`)

## License

MIT
