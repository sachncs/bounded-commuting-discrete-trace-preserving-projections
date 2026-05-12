# Architecture

## Module Map

```
src/lib/
  bcdtpp.js                   — Main API: Bcdtpp class
  mesh.js                     — Tetrahedral mesh topology, geometry, adjacency
  whitney.js                  — Barycentric coords, Whitney edge/face basis
  quadrature.js               — Gaussian quadrature on triangles, tetrahedra, lines
  math_utils.js               — Linear algebra: LU, inverse, vector ops
  local_solver.js             — Boundary-patch stiffness + constrained solve
  point_locator.js            — AABB tree for O(log N) point location
  higher_order_projection.js  — Bubble basis assembly, L2 monomial projection
  boundary_weight_computer.js — Alfeld/Worsey-Farin split + weight solves
  mesh_refinement.js          — Alfeld and Worsey-Farin mesh splits
  errors.js                   — Custom error classes
  projectors/
    h1_projector.js           — Pi^0 (vertex-based)
    hcurl_projector.js        — Pi^1 (edge-based)
    hdiv_projector.js         — Pi^2 (face-based)
    l2_projector.js           — Pi^3 (cell-based)
```

## Data Flow

```
User Input
    |
    v
Mesh (vertices, tets)
    |
    +---> Whitney (barycentric coords, basis)
    |         |
    |         v
    |     Bcdtpp (API facade)
    |         |
    |         +---> BoundaryWeightComputer
    |         |         +---> MeshRefinement (Alfeld / Worsey-Farin)
    |         |         +---> LocalSolver
    |         |
    |         +---> PointLocator (AABB tree)
    |         |
    |         +---> HigherOrderProjection (bubble / L2 enrichment)
    |         |
    |         +---> Projectors (H1, Hcurl, Hdiv, L2)
    |                   +---> math_utils (LU, vector ops)
    |
    +---> Quadrature (integrals)
```

## Projector Hierarchy

Each projector implements the same interface:

```
project(u, point, tIdx, boundaryFaceSet) -> value
```

- `H1Projector`: Vertex DoFs + boundary-weighted vertex values.
- `HcurlProjector`: Edge DoFs + line-integral constraints on boundary edges.
- `HdivProjector`: Face DoFs + normal-flux constraints on boundary faces.
- `L2Projector`: Cell average via quadrature.

## Key Design Decisions

1. **Pure ES modules**: No CommonJS in source; bundlers handle multi-format output.
2. **Zero external runtime dependencies**: All linear algebra is native JavaScript.
3. **Immutable mesh inputs**: `Mesh` validates and freezes topology at construction.
4. **Lazy caching**: `Whitney` caches per-tet barycentric gradients; `Bcdtpp` caches boundary weights on demand.
5. **Warning instead of throwing for local failures**: `BoundaryWeightComputer` and `HigherOrderProjection` warn on singular matrices so that a single bad element does not crash the entire mesh projection.
