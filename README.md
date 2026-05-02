# BCDTPP.js

JavaScript implementation of **Bounded, Commuting, Discrete-trace Preserving Projections** for the 3D de Rham complex on simplicial meshes.

Based on the paper: *Ern, Guzmán, Potu (2026) arXiv:2604.28103v1*.

## Features

- **Simplicial Mesh Support**: 3D tetrahedral mesh processing with full topology (faces, edges, incidence).
- **Geometric Splits**: 
  - **Alfeld Split**: Boundary faces subdivision for local solvers.
  - **Worsey-Farin Split**: Bulk tetrahedron subdivision for weight extensions.
- **Full de Rham Complex (Lowest-Order)**:
  - $l=0$ ($H^1$): Vertex-based projections using surface Poisson solvers.
  - $l=1$ ($H(\text{curl})$): Edge-based projections using edge Whitney forms.
  - $l=2$ ($H(\text{div})$): Face-based projections using face Whitney forms.
  - $l=3$ ($L^2$): Cell-based projections via cell averages.
- **Commuting Properties**: Designed to commute with exterior derivatives while preserving discrete traces.
- **High-Fidelity Solvers**: Uses `mathjs` for sparse matrix assembly and linear solves on boundary patches.

## Installation

```bash
npm install
```

## Usage

### 1. Define a Mesh
```javascript
const mesh = new Mesh(vertices, tetrahedra);
```

### 2. Initialize Projections
```javascript
const whitney = new Whitney(mesh);
const bcdtpp = new BCDTPP(mesh, whitney);

// Compute boundary weights (Alfeld splits + local solves)
await bcdtpp.computeAllWeights();
```

### 3. Project a function
```javascript
const u = (p) => Math.sin(p[0]);
const val = bcdtpp.projectH1(u, [0.5, 0.5, 0.5], 0);
```

## Mathematical Implementation Details

This library implements the construction of boundary correction operators $\Pi_\partial^l$.
The final projection is defined as:
$$\Pi^l = \Pi_\partial^l + \mathring{\Pi}^l(I - \Pi_\partial^l)$$

The novel contribution of this paper is the construction of $\Pi_\partial^l$ using local problems on subdivided patches to ensure stability and trace-preservation.

## License
MIT
