/**
 * Local patch stiffness assembly and constrained solves for boundary
 * weight computation.
 */

import { dot, cross, subtract, norm, zeros, luSolve, infinityNorm } from './math_utils.js'
import { SingularMatrixError } from './errors.js'

export class LocalSolver {
  /**
   * Assembles the surface stiffness matrix for -Delta_Gamma.
   * @param {!Array<!Array<number>>} vertices
   * @param {!Array<!Array<number>>} triangles
   * @return {!Array<!Array<number>>}
   */
  static assembleSurfaceStiffness (vertices, triangles) {
    const n = vertices.length
    const K = zeros(n, n)

    triangles.forEach((tri) => {
      const v = tri.map((i) => vertices[i])
      const ke = LocalSolver.#triangleStiffness(v)

      for (let i = 0; i < 3; i++) {
        for (let j = 0; j < 3; j++) {
          K[tri[i]][tri[j]] += ke[i][j]
        }
      }
    })

    return K
  }

  /**
   * Local stiffness matrix for -Delta_Gamma on a single triangle.
   * @param {!Array<!Array<number>>} v - Triangle vertices.
   * @return {!Array<!Array<number>>}
   * @private
   */
  static #triangleStiffness (v) {
    const v1 = subtract(v[1], v[0])
    const v2 = subtract(v[2], v[0])

    const c = cross(v1, v2)
    const area = 0.5 * norm(c)

    if (area < 1e-12) {
      throw new SingularMatrixError(`Degenerate triangle in stiffness assembly: area=${area}`)
    }

    const G = [subtract(v[2], v[1]), subtract(v[0], v[2]), subtract(v[1], v[0])]

    const ke = zeros(3, 3)
    for (let i = 0; i < 3; i++) {
      for (let j = 0; j < 3; j++) {
        ke[i][j] = dot(G[i], G[j]) / (4 * area)
      }
    }
    return ke
  }

  /**
   * Solves K x = b with a mean-zero constraint sum(x) = 0.
   *
   * This implementation enforces the constraint by replacing the last row of K
   * with ones and setting the last entry of b to zero. This is a simple
   * textbook approach that works well for small patch sizes (valence < 20).
   * It destroys symmetry and can degrade conditioning for larger systems; for
   * production-scale patches a Lagrange-multiplier or projected-gradient method
   * is preferable.
   *
   * @param {!Array<!Array<number>>} K
   * @param {!Array<number>} b
   * @return {!Array<number>}
   */
  static solveWithConstraint (K, b, onWarning = console.warn) {
    const n = K.length
    if (n === 0) {
      return []
    }
    const KDense = K.map((row) => [...row])
    for (let j = 0; j < n; j++) {
      KDense[n - 1][j] = 1.0
    }
    const bMod = [...b]
    bMod[n - 1] = 0

    const normEst = infinityNorm(KDense)
    if (normEst > 1e12) {
      onWarning({
        code: 'LOCAL_SOLVER_ILL_CONDITIONED',
        severity: 'warn',
        message:
          `LocalSolver: matrix is ill-conditioned (norm=${normEst}). ` +
          'Results may be inaccurate.'
      })
    }

    try {
      return luSolve(KDense, bMod)
    } catch (e) {
      throw new SingularMatrixError(
        `LocalSolver: constrained solve failed (${e.message}). ` +
          'Patch matrix may be singular or ill-conditioned.'
      )
    }
  }
}
