# Error Taxonomy

This library uses a small hierarchy of custom error classes to provide actionable diagnostics.

## MeshValidationError

Thrown by `Mesh` when input data is geometrically or topologically invalid.

| Condition | Message pattern |
|-----------|-----------------|
| Empty tetrahedra array | `Mesh must contain at least one tetrahedron` |
| Non-finite vertex coordinate | `Vertex ... contains non-finite value` |
| Duplicate vertex indices in a tet | `Tetrahedron ... contains duplicate vertex indices` |
| Out-of-bounds vertex index | `Tetrahedron ... references out-of-bounds vertex index` |
| Non-integer tetrahedron index | `Tetrahedron ... contains non-integer vertex index` |
| Negative orientation (det <= 0) | `Tetrahedron ... is degenerate or negatively oriented` |

**Recovery**: Fix the input mesh data and re-instantiate `Mesh`.

## ProjectionError

Thrown by `Bcdtpp` and projector classes when arguments are invalid or operations cannot proceed.

| Condition | Message pattern |
|-----------|-----------------|
| Invalid `tIdx` type | `tIdx must be a non-negative integer` |
| Out-of-range `tIdx` | `tIdx ... is out of range` |
| Invalid point type | `point must be an array of 3 finite numbers` |
| Point outside mesh | `Point ... not found in mesh` |
| Unimplemented higher-order vector projection | `Higher-order vector projection (l=1 or l=2, p>0) is not yet implemented` |

**Recovery**: Validate inputs before calling projection methods; ensure the point locator is built for global queries.

## SingularMatrixError

Thrown by linear algebra routines in `math_utils.js` when a matrix is singular or numerically rank-deficient.

| Condition | Message pattern |
|-----------|-----------------|
| Zero pivot during LU | `Singular matrix encountered` |
| Zero determinant in `inverse3x3` | `Singular matrix encountered` |
| Zero determinant in `solve3x3` | `Singular matrix encountered` |

**Recovery**: Check mesh validity (degenerate elements cause singular Jacobians) or verify the input matrix.

## BoundaryWeightComputer Warnings

`BoundaryWeightComputer` logs warnings via an injected `warn` function rather than throwing, because per-vertex failures should not halt the entire computation.

| Code | Message | Meaning |
|------|---------|---------|
| `BWC_ZERO_STAR_AREA` | `Vertex ... has zero star area; skipping weight computation` | The vertex patch has no geometric area, likely due to a degenerate mesh. |
| `BWC_VERTEX_FAILURE` | `Failed to compute weights for vertex ...` | An unexpected error occurred during the local solve for this vertex. |

**Recovery**: Inspect the mesh near the reported vertex for degenerate or inverted elements.

## HigherOrderProjection Warnings

`HigherOrderProjection` also warns rather than throwing for singular mass matrices.

| Code | Message | Meaning |
|------|---------|---------|
| `HOP_SINGULAR_MASS` | `Singular mass matrix for bubble projection ...` | The bubble mass matrix is singular; bubble coefficients are skipped. |
| `HOP_SINGULAR_L2` | `Singular L2 mass matrix ...` | The L2 mass matrix is singular (often zero volume); coefficients are empty. |

**Recovery**: Ensure the tetrahedron has strictly positive volume.
