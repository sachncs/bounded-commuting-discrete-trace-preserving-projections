# BCDTPP API Reference

Bounded, Commuting, Discrete-trace Preserving Projections on simplicial meshes.

## Table of Contents

- [Bcdtpp](#bcdtpp) — Main entry point
- [Mesh](#mesh) — Tetrahedral mesh data structure
- [Whitney](#whitney) — Barycentric coordinates and Whitney basis functions
- [BoundaryWeightComputer](#boundaryweightcomputer) — Boundary patch weights
- [HigherOrderProjection](#higherorderprojection) — Higher-order projection framework
- [LocalSolver](#localsolver) — Surface stiffness assembly and constrained solve
- [MeshRefinement](#meshrefinement) — Alfeld and Worsey-Farin splits
- [PointLocator](#pointlocator) — AABB tree point-in-tet queries
- [H1Projector](#h1projector) — H¹ (l=0) vertex-based projector
- [HcurlProjector](#hcurlprojector) — H(curl) (l=1) edge-based projector
- [HdivProjector](#hdivprojector) — H(div) (l=2) face-based projector
- [L2Projector](#l2projector) — L² (l=3) cell-based projector
- [Math Utilities](#math-utilities) — Pure JS linear algebra and vector functions
- [Error Classes](#error-classes)

---

## Bcdtpp

Main entry point implementing the global de Rham projection operator Πˡ.

```js
import { Bcdtpp } from 'bcdtpp/src/lib/bcdtpp.js'
```

### Constructor

```js
new Bcdtpp(mesh, whitney, options?)
```

| Param | Type | Default | Description |
|-------|------|---------|-------------|
| `mesh` | `Mesh` | — | Tetrahedral mesh. |
| `whitney` | `Whitney` | — | Whitney basis instance. |
| `options.quadratureOrder` | `number` | `3` | Gaussian quadrature order for integrations. |
| `options.onWarning` | `(ctx) => void` | `console.warn` | Callback for local projection failures. |

Validates that the mesh contains no degenerate (zero-volume) tetrahedra.

### Methods

#### `buildPointLocator()`

Builds an AABB tree for O(log N) point-in-tetrahedron queries. Called automatically by `projectAtPoint` if not already built.

#### `computeBoundaryWeights()`

Computes boundary patch weights (§6.3.1). Triggers a Worsey-Farin split if the mesh lacks face barycenter vertices. Must be called before `projectH1`, `projectHcurl`, or `projectHdiv`.

#### `project(u, point, tIdx, l, p?)`

Global projector Πˡ implementing the decomposition formula:

```
Πˡ = Π_partialˡ + Π_ringˡ (I − Π_partialˡ)
```

| Param | Type | Description |
|-------|------|-------------|
| `u` | `(pt) => number \| Array<number>` | Input function (scalar for l=0,3; vector for l=1,2). |
| `point` | `[x, y, z]` | Evaluation point inside the tetrahedron. |
| `tIdx` | `number` | Tetrahedron index. |
| `l` | `number` | Form degree: 0 (H¹), 1 (H(curl)), 2 (H(div)), 3 (L²). |
| `p` | `number` | Polynomial degree (default 0). |

Returns `number` (l=0,3) or `Array<number>` (l=1,2).

#### `projectAtPoint(u, point, l?, p?)`

Convenience wrapper that finds the containing tetrahedron via the AABB tree, then projects. Auto-builds the point locator on first call.

Returns `{ value, tIdx, bary }` where `value` is the projected value and `tIdx`/`bary` identify the containing tetrahedron.

#### `projectHp(u, point, tIdx, l, p)`

Higher-order projection for form degree l and polynomial degree p ≥ 1. Delegates to the specialized `HigherOrderProjection` instance.

#### `projectH1(u, point, tIdx)` / `projectHcurl(u, point, tIdx)` / `projectHdiv(u, point, tIdx)` / `projectL2(u, tIdx)`

Lowest-order (p=0) projections for each form degree. `projectH1`, `projectHcurl`, and `projectHdiv` require `computeBoundaryWeights()` to have been called.

#### `projectRing(u, point, tIdx, l)` / `extendBoundary(boundaryData, point, tIdx, l)` / `projectPartial(u, point, tIdx, l)`

Decomposition components: interior projector, discrete extension operator, and boundary correction part.

#### `extractBoundaryDofs(u, l)`

Extracts boundary degrees of freedom for form degree l.

Returns `Map<number, number>` mapping boundary entity indices to coefficient values.

#### `quadratureOrder` (getter)

Returns the quadrature order used for integrations.

---

## Mesh

Tetrahedral mesh data structure with topological connectivity, geometry, and boundary extraction.

```js
import { Mesh } from 'bcdtpp/src/lib/mesh.js'
```

### Constructor

```js
new Mesh(vertices, tetrahedra)
```

| Param | Type | Description |
|-------|------|-------------|
| `vertices` | `Array<[x, y, z]>` | Vertex coordinates. |
| `tetrahedra` | `Array<[i0, i1, i2, i3]>` | Tetrahedra as vertex index arrays. |

Validates orientation (positive signed volume), bounds, and uniqueness. Throws `MeshValidationError` on invalid input.

### Topology Fields

| Field | Type | Description |
|-------|------|-------------|
| `vertices` | `Array<Array<number>>` | Vertex coordinates. |
| `tetrahedra` | `Array<Array<number>>` | Tetrahedra as vertex index arrays. |
| `originalVertexCount` | `number` | Vertex count before refinement. |
| `vertexCount` | `number` | Current vertex count. |
| `tetrahedronCount` | `number` | Number of tetrahedra. |
| `faces` | `Array<Array<number>>` | All faces as vertex index arrays. |
| `edges` | `Array<Array<number>>` | All edges as vertex index pairs. |
| `boundaryFaces` | `Array<number>` | Indices of boundary faces. |
| `boundaryEdges` | `Array<number>` | Indices of boundary edges. |
| `boundaryNodes` | `Set<number>` | Boundary vertex indices. |
| `faceToTets` | `Array<Array<number>>` | Face → incident tetrahedra adjacency. |
| `edgeToFaces` | `Array<Array<number>>` | Edge → incident faces adjacency. |
| `vertexToTets` | `Array<Array<number>>` | Vertex → incident tetrahedra. |
| `vertexToBoundaryFaces` | `Array<Array<number>>` | Vertex → incident boundary faces. |
| `vertexToEdges` | `Array<Array<number>>` | Vertex → incident edges. |

### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `getFaceArea(fIdx)` | `number` | Geometric area of a face. |
| `getFaceOutwardNormal(fIdx)` | `[x, y, z]` | Outward-pointing unit normal for a boundary face. |
| `getEdgeIndex(edgeKey)` | `number \| undefined` | Global edge index for an edge key. |
| `getFaceIndex(faceKey)` | `number \| undefined` | Global face index for a face key. |
| `getTetEdgeSign(tIdx, e)` | `number` | Edge orientation sign (±1) for local edge e (0–5). |
| `getTetFaceSign(tIdx, f)` | `number` | Face orientation sign (±1) for local face f (0–3). |
| `getFaceBarycenter(fIdx)` | `[x, y, z]` | Cached centroid of a triangle face. |
| `getTetrahedronBarycenter(tIdx)` | `[x, y, z]` | Cached centroid of a tetrahedron. |
| `getStar(vIdx)` | `Array<number>` | Tetrahedra incident on a vertex. |
| `getBoundaryStar(vIdx)` | `Array<number>` | Boundary faces incident on a vertex. |
| `getEdgeStar(vIdx)` | `Array<number>` | Edges incident on a vertex. |
| `getVolume(tIdx)` | `number` | Cached volume of a tetrahedron. |
| `getTetrahedronFaces(tIdx)` | `Array<number>` | Four global face indices (opposite-vertex convention). |
| `getVertices()` | `Array<Array<number>>` | All vertex coordinates. |
| `getTetrahedra()` | `Array<Array<number>>` | All tetrahedra. |
| `getFaces()` | `Array<Array<number>>` | All faces. |
| `getEdges()` | `Array<Array<number>>` | All edges. |
| `getBoundaryNodes()` | `Set<number>` | Boundary node indices. |
| `getBoundaryFaces()` | `Array<number>` | Boundary face indices. |
| `getBoundaryEdges()` | `Array<number>` | Boundary edge indices. |
| `getOriginalVertexCount()` | `number` | Original vertex count before refinement. |

### Static Methods

#### `Mesh.computeEdgeKey(a, b, vertexCount)`

Returns a canonical integer key for an edge. Safe for `vertexCount < 2²⁶`.

---

## Whitney

Barycentric coordinate computation and Whitney finite-element basis functions.

```js
import { Whitney } from 'bcdtpp/src/lib/whitney.js'
```

### Constructor

```js
new Whitney(mesh)
```

Precomputes per-tet edge matrices, inverses, determinants, and barycentric gradients.

### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `getBarycentric(tIdx, point)` | `[λ0, λ1, λ2, λ3]` | Barycentric coordinates of a point w.r.t. a tetrahedron. |
| `getGradBarycentric(tIdx)` | `Array<[x,y,z]>` | Gradients of the four barycentric coordinate functions. |
| `getEdgeBasis(tIdx, bary)` | `Array<[x,y,z]>` | Six Nédélec (Whitney 1-form) edge basis vectors in order `[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]`. |
| `getFaceBasis(tIdx, bary)` | `Array<[x,y,z]>` | Four Raviart-Thomas (Whitney 2-form) face basis vectors indexed by opposite vertex 0–3. |

---

## BoundaryWeightComputer

Computes boundary patch weights used by trace-preserving projection operators.

```js
import { BoundaryWeightComputer } from 'bcdtpp/src/lib/boundary_weight_computer.js'
```

### Constructor

```js
new BoundaryWeightComputer(mesh, meshRefinement, onWarning?)
```

| Param | Type | Default | Description |
|-------|------|---------|-------------|
| `mesh` | `Mesh` | — | Tetrahedral mesh. |
| `meshRefinement` | `MeshRefinement` | — | Mesh refinement instance. |
| `onWarning` | `(ctx) => void` | `console.warn` | Callback for local computation failures. |

### Methods

#### `compute()`

Returns:
```js
{
  vertexBoundaryData: Map<number, { nodeMap, invNodeMap, psi }>,
  edgeBoundaryData:  Map<number, { v0, v1, tangent, length }>,
  faceBoundaryData:  Map<number, { normal, area }>
}
```

---

## HigherOrderProjection

Higher-order projection framework (§7). Builds projections for p ≥ 1 from lowest-order via bubble corrections on Alfeld-split patches.

```js
import { HigherOrderProjection } from 'bcdtpp/src/lib/higher_order_projection.js'
```

### Constructor

```js
new HigherOrderProjection(mesh, whitney, quadratureOrder?, onWarning?)
```

| Param | Type | Default | Description |
|-------|------|---------|-------------|
| `mesh` | `Mesh` | — | Tetrahedral mesh. |
| `whitney` | `Whitney` | — | Whitney basis instance. |
| `quadratureOrder` | `number` | `3` | Quadrature order. |
| `onWarning` | `(ctx) => void` | `console.warn` | Callback for projection failures. |

### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `evaluateBubbleBasis(bary, p)` | `Array<number>` | Evaluates scalar bubble basis at barycentric coordinates for degree p. |
| `assembleBubbleMass(tIdx, p)` | `Array<Array<number>>` | Assembles local mass matrix for the bubble space. |
| `assembleBubbleRhs(tIdx, p, residualFn)` | `Array<number>` | Assembles RHS for bubble projection. |
| `solveBubbleProjection(tIdx, p, residualFn)` | `Array<number> \| null` | Solves for bubble coefficients. Returns null if p < 4 or solve fails. |
| `evaluateBubble(tIdx, p, coeffs, point)` | `number` | Evaluates bubble correction at a point. |
| `evaluateBernsteinBasis(bary, p)` | `Array<number>` | Evaluates P^p Bernstein basis at barycentric coordinates. |
| `solveL2Projection(tIdx, p, u)` | `Array<number>` | Solves L2 projection of u onto P^p(T). Returns Bernstein coefficients. |
| `evaluateL2Projection(coeffs, bary, p)` | `number` | Evaluates L2 projection at a point. |

---

## LocalSolver

Static utility for surface stiffness assembly and constrained linear systems.

```js
import { LocalSolver } from 'bcdtpp/src/lib/local_solver.js'
```

### Static Methods

#### `LocalSolver.assembleSurfaceStiffness(vertices, triangles)`

Assembles the surface stiffness matrix for −Δ_Γ.

| Param | Type | Description |
|-------|------|-------------|
| `vertices` | `Array<Array<number>>` | Vertex coordinates. |
| `triangles` | `Array<Array<number>>` | Triangle connectivity. |

Returns `Array<Array<number>>` — the N×N stiffness matrix.

#### `LocalSolver.solveWithConstraint(K, b, onWarning?)`

Solves Kx = b with a mean-zero constraint Σx = 0. Uses Gaussian elimination with partial pivoting and detects ill-conditioning.

| Param | Type | Default | Description |
|-------|------|---------|-------------|
| `K` | `Array<Array<number>>` | — | Stiffness matrix. |
| `b` | `Array<number>` | — | Right-hand side. |
| `onWarning` | `(ctx) => void` | `console.warn` | Callback for singularity warnings. |

Returns `Array<number>` — the solution vector.

---

## MeshRefinement

Mesh refinement operators implementing Alfeld face splitting (§6.1.3) and Worsey-Farin tetrahedron splitting (§6.1.4).

```js
import { MeshRefinement } from 'bcdtpp/src/lib/mesh_refinement.js'
```

### Constructor

```js
new MeshRefinement(mesh)
```

### Methods

#### `computeAlfeldSplit()`

Alfeld split: inserts a face barycenter vertex for each boundary face and subdivides. Mutates `mesh.vertices`. Idempotent.

#### `computeWorseyFarinSplit()`

Worsey-Farin split: splits each tetrahedron using face barycenters. Calls Alfeld first if needed. Mutates `mesh.vertices`. Idempotent.

### Fields

| Field | Description |
|-------|-------------|
| `faceBarycenters` | Computed face barycenter coordinates. |
| `tetBarycenters` | Computed tet barycenter coordinates. |
| `alfeldTriangles` | Triangles from Alfeld split. |
| `faceToAlfeld` | Map from face index to Alfeld split data. |
| `worseyFarinTetrahedra` | Tets from Worsey-Farin split. |

---

## PointLocator

AABB tree for O(log N) point-in-tetrahedron queries. Builds a balanced tree via median-split along the longest axis.

```js
import { PointLocator } from 'bcdtpp/src/lib/point_locator.js'
```

### Constructor

```js
new PointLocator(mesh, maxLeafSize?)
```

| Param | Type | Default | Description |
|-------|------|---------|-------------|
| `mesh` | `Mesh` | — | Tetrahedral mesh. |
| `maxLeafSize` | `number` | `8` | Maximum tets per leaf node. |

### Methods

#### `findTetrahedron(point)`

Returns `{ tIdx, bary }` containing the tetrahedron index and barycentric coordinates, or `null` if not found.

---

## H1Projector

Lowest-order H¹ (l=0) vertex-based projector. Projects scalar functions onto continuous piecewise-linear (P1 Lagrange) space.

```js
import { H1Projector } from 'bcdtpp/src/lib/projectors/h1_projector.js'
```

### Constructor

```js
new H1Projector(mesh, whitney, meshRefinement)
```

### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `project(u, point, tIdx, vertexBoundaryData)` | `number` | H¹ projection of scalar u at a point. |
| `computeBoundaryIntegralH1(vIdx, u, vertexBoundaryData)` | `number` | Computes boundary integral DoF for a boundary vertex. |
| `projectRing(u, point, tIdx)` | `number` | Π_ring^0: interior projector with zero boundary trace. |
| `extendBoundary(boundaryData, point, tIdx)` | `number` | E^0: discrete extension of vertex boundary data. |

---

## HcurlProjector

Lowest-order H(curl) (l=1) edge-based projector. Projects vector functions onto the Nédélec first-kind (Whitney 1-form) space.

```js
import { HcurlProjector } from 'bcdtpp/src/lib/projectors/hcurl_projector.js'
```

### Constructor

```js
new HcurlProjector(mesh, whitney, quadratureOrder)
```

### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `project(u, point, tIdx, boundaryEdgeSet)` | `Array<number>` | H(curl) projection at a point. |
| `computeEdgeDof(u, eIdx)` | `number` | Exact edge DoF: ∫_e u·t ds. For scalar u reduces to u(v1)−u(v0). |
| `computeInteriorEdgeCoeff(u, tIdx, i, j, isScalar)` | `number` | Interior edge coefficient in local orientation. |
| `projectRing(u, point, tIdx, boundaryEdgeSet)` | `Array<number>` | Π_ring^1: interior projector with zero boundary trace. |
| `extendBoundary(boundaryData, point, tIdx, boundaryEdgeSet)` | `Array<number>` | E^1: discrete extension of edge boundary data. |

---

## HdivProjector

Lowest-order H(div) (l=2) face-based projector. Projects vector functions onto the Raviart-Thomas (Whitney 2-form) space.

```js
import { HdivProjector } from 'bcdtpp/src/lib/projectors/hdiv_projector.js'
```

### Constructor

```js
new HdivProjector(mesh, whitney, quadratureOrder)
```

### Methods

| Method | Returns | Description |
|--------|---------|-------------|
| `project(u, point, tIdx, boundaryFaceSet)` | `Array<number>` | H(div) projection at a point. |
| `computeFaceDof(u, fIdx)` | `number` | Exact face DoF: ∫_f u·n dA. For scalar u integrates grad(u)·n. |
| `computeInteriorFaceCoeff(u, tIdx, f, isScalar)` | `number` | Interior face coefficient in local outward normal orientation. |
| `projectRing(u, point, tIdx, boundaryFaceSet)` | `Array<number>` | Π_ring^2: interior projector with zero boundary trace. |
| `extendBoundary(boundaryData, point, tIdx, boundaryFaceSet)` | `Array<number>` | E^2: discrete extension of face boundary data. |

---

## L2Projector

Lowest-order L² (l=3) cell-based projector. Projects scalar functions onto piecewise constants (P0).

```js
import { L2Projector } from 'bcdtpp/src/lib/projectors/l2_projector.js'
```

### Constructor

```js
new L2Projector(mesh, whitney, quadratureOrder)
```

### Methods

#### `project(u, tIdx)`

Returns the L² projection: volume-weighted average of u over the tetrahedron.

---

## Math Utilities

Pure JavaScript linear algebra and vector utilities. No external dependencies.

```js
import { dot, cross, norm, ... } from 'bcdtpp/src/lib/math_utils.js'
```

### Vector Operations

| Function | Signature | Returns | Description |
|----------|-----------|---------|-------------|
| `dot(a, b)` | `(number[], number[]) → number` | `number` | Dot product of two 3D vectors. |
| `cross(a, b)` | `(number[], number[]) → number[]` | `number[]` | Cross product of two 3D vectors. |
| `subtract(a, b)` | `(number[], number[]) → number[]` | `number[]` | Vector subtraction a − b. |
| `subtractInto(a, b, out)` | `(number[], number[], number[]) → number[]` | `number[]` | In-place subtraction: out = a − b. |
| `crossInto(a, b, out)` | `(number[], number[], number[]) → number[]` | `number[]` | In-place cross product: out = a × b. |
| `norm(v)` | `(number[]) → number` | `number` | Euclidean norm of a vector. |

### Geometry

| Function | Signature | Returns | Description |
|----------|-----------|---------|-------------|
| `tetDeterminant(v0, v1, v2, v3)` | `(number[], number[], number[], number[]) → number` | `number` | Scalar triple product; \|det\|/6 = volume. |
| `tetVolume(v0, v1, v2, v3)` | `(number[], number[], number[], number[]) → number` | `number` | Geometric volume of a tetrahedron. |
| `triangleArea(p1, p2, p3)` | `(number[], number[], number[]) → number` | `number` | Area of a triangle. |
| `factorial(n)` | `(number) → number` | `number` | Computes n!. Throws if n > 170. |

### Matrix Operations

| Function | Signature | Returns | Description |
|----------|-----------|---------|-------------|
| `zeros(rows, cols)` | `(number, number) → number[][]` | `number[][]` | Zero-initialized dense matrix. |
| `infinityNorm(a)` | `(number[][]) → number` | `number` | Infinity norm (max absolute row sum). |
| `luSolve(a, b)` | `(number[][], number[]) → number[]` | `number[]` | Solves Ax = b via Gaussian elimination with partial pivoting and row equilibration. |
| `inverse3x3(m)` | `(number[][]) → number[][]` | `number[][]` | Inverse of a 3×3 matrix. |
| `solve3x3(a, b)` | `(number[][], number[]) → number[]` | `number[]` | Solves a 3×3 system via Cramer's rule. |

### Numerical Methods

| Function | Signature | Returns | Description |
|----------|-----------|---------|-------------|
| `numericalGradient(u, pt, h?)` | `((number[]→number), number[], number?) → number[]` | `number[]` | Numerical gradient via central differences. Default h = 1e-6. |

### Constants

| Name | Type | Value | Description |
|------|------|-------|-------------|
| `MAX_SAFE_FACTORIAL` | `number` | `170` | Maximum n for which n! fits in a JavaScript Number. |

---

## Error Classes

All error classes extend `Error` and are exported from `bcdtpp/src/lib/errors.js`.

### MeshValidationError

Thrown when mesh input data fails validation (e.g., degenerate tetrahedra, out-of-bounds indices, non-unique vertices).

```js
new MeshValidationError(message, index?)
```

| Param | Type | Description |
|-------|------|-------------|
| `message` | `string` | Description of the validation failure. |
| `index` | `number` | The offending element index, if applicable. |

### ProjectionError

Thrown when a projection cannot be computed (e.g., singular matrix, invalid arguments).

```js
new ProjectionError(message)
```

### SingularMatrixError

Thrown when a linear system is singular or numerically ill-conditioned.

```js
new SingularMatrixError(message)
```
