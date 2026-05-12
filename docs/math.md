# Mathematical Background

This document provides a concise mathematical introduction to the operators implemented in this library.

## The 3D de Rham Complex

On a contractible domain, the de Rham complex is the exact sequence:

```
H^1  --grad-->  H(curl)  --curl-->  H(div)  --div-->  L^2  -->  0
```

For a bounded Lipschitz domain, the complex is:

```
H^1_0  --grad-->  H_0(curl)  --curl-->  H_0(div)  --div-->  L^2_0  -->  0
```

## Discrete Spaces

On a simplicial mesh, the lowest-order finite-element spaces that form a discrete de Rham complex are:

- **P^1_0 (H^1)**: continuous piecewise-linear Lagrange elements (vertex DoFs).
- **N^0_1 (H(curl))**: Nédélec edge elements (edge DoFs).
- **RT^0_1 (H(div))**: Raviart-Thomas face elements (face DoFs).
- **P^0_3 (L^2)**: piecewise constants (cell DoFs).

## Projection Operators

The library implements bounded, commuting, discrete-trace preserving projections `Pi^l` for `l = 0,1,2,3`.

### Decomposition

Each projection is decomposed as:

```
Pi^l = Pi_partial^l + Pi_ring^l (I - Pi_partial^l)
```

Where:

- **Pi_partial^l**: Prescribes boundary data exactly and extends it into the interior.
- **Pi_ring^l**: The interior projector with vanishing trace on the boundary.

### Trace Preservation

The operators satisfy:

```
tr^l(Pi^l v) = tr^l(v)    on boundary faces/edges/vertices
```

Where `tr^l` is the canonical trace operator for the `l`-form space.

### Commuting Diagram

The projections commute with the exterior derivative:

```
d^l Pi^l = Pi^{l+1} d^l
```

In concrete terms:
- `grad Pi^0 = Pi^1 grad`
- `curl Pi^1 = Pi^2 curl`
- `div Pi^2 = Pi^3 div`

## Whitney Forms

Barycentric coordinates `lambda_i` form the basis for `P^1_0`. Their gradients are constant per tetrahedron.

The Whitney edge basis for `N^0_1` is:

```
phi_{ij} = lambda_i grad(lambda_j) - lambda_j grad(lambda_i)
```

The Whitney face basis for `RT^0_1` is:

```
psi_{ijk} = 2 (lambda_i grad(lambda_j) x grad(lambda_k) + cyclic)
```

## Mesh Refinement

Boundary weights require solving local problems on subdivided patches:

- **Alfeld split**: A tetrahedron is split into 4 sub-tetrahedra by connecting the barycenter to each face.
- **Worsey-Farin split**: A tetrahedron is split into 12 sub-tetrahedra by connecting face barycenters to the tet barycenter.

These splits are used to construct stable local solvers for the boundary correction operators.

## Higher-Order Extensions

For polynomial degree `p > 0`:

- **H^1 (l=0)**: Enriched with bubble functions that vanish on the element boundary.
- **L^2 (l=3)**: Enriched with monomials of total degree up to `p`.

Vector-valued higher-order (`l = 1, 2`, `p > 0`) is not yet implemented.
