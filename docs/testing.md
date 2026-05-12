# Testing Guide

## Philosophy

This project maintains high test coverage (currently >98% statements, >93% branches). Tests are grouped by module and behavior, not by implementation detail.

## Test Organization

| File | Scope |
|------|-------|
| `bcdtpp.test.js` | Core projections (H1, Hcurl, Hdiv, L2, point location) |
| `coverage.test.js` | Multi-tet meshes, mesh validation, bubble projection, boundary weight faults |
| `convergence.test.js` | Convergence rates on refined cube meshes |
| `interior_extension.test.js` | Ring projector, extension operator, decomposition properties |
| `local_solver.test.js` | Surface stiffness assembly, constrained solves |
| `math_utils.test.js` | Vector operations, LU solver, matrix inversion |
| `mesh.test.js` | Mesh topology, volume, boundary stars |
| `point_locator.test.js` | AABB tree point-in-tet queries |
| `quadrature.test.js` | Gaussian quadrature exactness |
| `whitney.test.js` | Barycentric coordinates, Whitney basis functions |

## Running Tests

```bash
npm test              # Single run
npm run test:watch    # Watch mode
npm run test:coverage # With c8 coverage report
```

## Writing New Tests

Use Mocha's `describe`/`it` and Chai's `expect`:

```javascript
import { expect } from 'chai'
import { Mesh } from '../src/lib/mesh.js'

describe('Mesh', () => {
  it('computes positive volume', () => {
    const mesh = new Mesh(...)
    expect(mesh.getVolume(0)).to.be.above(0)
  })
})
```

## Coverage Thresholds

While not enforced by CI, the project aims for:
- Statements: >= 95%
- Branches: >= 90%
- Functions: >= 95%
- Lines: >= 95%

If you add new code, add tests that exercise both success and error paths.
