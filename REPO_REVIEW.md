# Repository Review: BCDTPP Implementation

This document lists the technical debt, missing features, and numerical issues identified in the current implementation of the BCDTPP projection framework.

## 1. Missing Mathematical Features

### Higher-Order Projections ($p \ge 1$)
- **Status**: Not implemented.
- **Details**: The current implementation only handles the lowest-order case ($p=0$). Section 7 of the paper defines a recursive construction for $p \ge 1$ using bubble functions and higher-degree local problems. This is a significant missing piece for a "complete" implementation of the paper.

### Full $H(\text{curl})$ and $H(\text{div})$ Integration
- **Status**: Simplified / Placeholder.
- **Details**: While the methods `projectHcurl` and `projectHdiv` exist, the computation of the boundary coefficients $\alpha_e$ and $\alpha_f$ currently uses simplified integration. A faithful implementation requires solving the specific surface-curl and surface-divergence problems on the Alfeld split patches as defined in §6.3.2 and §6.3.3.

### Interior Projections ($\mathring{\Pi}^l$)
- **Status**: Proxy implementation.
- **Details**: The code uses a standard nodal interpolation as a proxy for the interior projection $\mathring{\Pi}^l$. To be fully faithful, a commuting zero-trace projection (e.g., from [11]) should be integrated.

## 2. Numerical & Algorithmic Issues

### Integration Accuracy
- **Status**: Low order.
- **Details**: Integrals are currently computed using midpoint quadrature. For de Rham complexes, exactness properties (like $\operatorname{grad} P^0 = P^1 \operatorname{grad}$) often rely on higher-order Gaussian quadrature (e.g., Keast or Dunavant rules for tetrahedra).

### Point Location Performance
- **Status**: Naive ($O(1)$/Fixed).
- **Details**: Evaluation currently assumes tetrahedron index 0 or requires the caller to know the tetrahedron index. A production library requires an **AABB Tree** or **Spatial Hash** to perform point-in-tetrahedron searches in $O(\log N)$ time.

### Computational Performance
- **Status**: Single-threaded.
- **Details**: The Alfeld/Worsey-Farin splitting and the local matrix solves are computationally expensive. For meshes with $>1000$ tetrahedra, these operations should be moved to a **WebWorker** to prevent UI freezing.

## 3. Engineering & Code Quality

### Lack of Automated Tests
- **Status**: Critical.
- **Details**: There is no `vitest` or `jest` suite to verify the commuting properties ($\text{d} \Pi^l = \Pi^{l+1} \text{d}$). These tests are standard for validating de Rham complex implementations.

### Dependency Management
- **Status**: Minimal.
- **Details**: The repo relies heavily on `mathjs`. While powerful, it may be slow for large sparse systems. A specialized WASM-based solver (like `Eigen`) would be more appropriate for large-scale FEM.

### Missing Error Handling
- **Status**: Moderate.
- **Details**: Degenerate tetrahedra (zero volume) will cause `mathjs` matrix inversions to fail without graceful recovery.

---
*Review generated on: 2026-05-03*
