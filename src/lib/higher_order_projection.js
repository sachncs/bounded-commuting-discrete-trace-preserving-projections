/**
 * Higher-Order Projection Framework (Section 7).
 *
 * For p >= 1, the projection Pi^l_p is built recursively from the lowest-order
 * projection Pi^l_0 by adding bubble corrections on Alfeld-split patches.
 */

import { zeros, luSolve, factorial } from './math_utils.js'
import {
  barycentricToCartesian,
  tetrahedronQuadrature,
  compositeTetrahedronQuadrature
} from './quadrature.js'

/**
 * Higher-order projection framework implementing Section 7 of the paper.
 *
 * For polynomial degree p >= 4 on H^1 (l=0), the projection Pi^l_p is built
 * from the lowest-order projection by adding bubble corrections on
 * Alfeld-split patches.  For L2 (l=3) with p >= 1, a Bernstein-basis
 * L2 projection is used instead.
 *
 * For p = 1, 2, 3 on H^1 the bubble space is empty (the product
 * b = lambda_0 * lambda_1 * lambda_2 * lambda_3 already has degree 4),
 * so projectHp returns the lowest-order projection without enrichment.
 */
export class HigherOrderProjection {
  /**
   * @param {!Mesh} mesh
   * @param {!Whitney} whitney
   * @param {number=} quadratureOrder
   * @param {function=} onWarning - Callback invoked with a warning context
   *   object when a local solve fails (e.g. singular mass matrix).
   */
  constructor (mesh, whitney, quadratureOrder = 3, onWarning = console.warn) {
    this.mesh = mesh
    this.whitney = whitney
    this.quadratureOrder = quadratureOrder
    this.onWarning =
      typeof onWarning === 'function'
        ? onWarning
        : (ctx) => console.warn(ctx.message ?? ctx)
    this.bubbleExponentCache = new Map()
  }

  /**
   * Generates the scalar bubble basis for polynomial degree p on a tetrahedron.
   * The bubble space consists of polynomials vanishing on the boundary.
   * For degree p, the bubble space is b * P^{p-4} where b = prod(lambda_i).
   *
   * **Note:** The bubble space is empty for p = 1, 2, 3 because the product
   * b = lambda0 * lambda1 * lambda2 * lambda3 already has degree 4. Therefore
   * projectHp with l=0 and p=1,2,3 returns the lowest-order projection without
   * enrichment. True higher-order enrichment begins at p >= 4.
   *
   * @param {number} p
   * @return {!Array<!Array<number>>} List of exponent tuples [a,b,c,d]
   *   such that basis function = lambda0^a * lambda1^b * lambda2^c * lambda3^d * b.
   */
  #getBubbleExponents (p) {
    if (this.bubbleExponentCache.has(p)) {
      return this.bubbleExponentCache.get(p)
    }
    const exponents = []
    if (p >= 4) {
      const m = p - 4
      for (let a = 0; a <= m; a++) {
        for (let b = 0; b <= m - a; b++) {
          for (let c = 0; c <= m - a - b; c++) {
            const d = m - a - b - c
            // Exponents for lambda_i include the +1 from b.
            exponents.push([a + 1, b + 1, c + 1, d + 1])
          }
        }
      }
    }
    this.bubbleExponentCache.set(p, exponents)
    return exponents
  }

  /**
   * Polynomial power with explicit 0^0 = 1 convention.
   * @param {number} base
   * @param {number} exp
   * @return {number}
   * @private
   */
  static #pow (base, exp) {
    if (base === 0 && exp === 0) return 1
    return Math.pow(base, exp)
  }

  /**
   * Evaluates the scalar bubble basis at a point given in barycentric coords.
   * @param {!Array<number>} bary
   * @param {number} p
   * @return {!Array<number>}
   */
  evaluateBubbleBasis (bary, p) {
    const exponents = this.#getBubbleExponents(p)
    if (exponents.length === 0) {
      return []
    }
    return exponents.map(([a, b, c, d]) =>
      HigherOrderProjection.#pow(bary[0], a) *
      HigherOrderProjection.#pow(bary[1], b) *
      HigherOrderProjection.#pow(bary[2], c) *
      HigherOrderProjection.#pow(bary[3], d)
    )
  }

  /**
   * Assembles the local mass matrix for the bubble space.
   * @param {number} tIdx
   * @param {number} p
   * @return {!Array<!Array<number>>}
   */
  assembleBubbleMass (tIdx, p) {
    const exponents = this.#getBubbleExponents(p)
    if (exponents.length === 0) {
      return []
    }

    const { bary, weights } = tetrahedronQuadrature(this.quadratureOrder)
    const vol = this.mesh.getVolume(tIdx)

    // Precompute basis values at quadrature points.
    const basisAtPoints = bary.map((lam) => this.evaluateBubbleBasis(lam, p))
    const n = basisAtPoints[0].length
    const M = zeros(n, n)

    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        let sum = 0
        for (let q = 0; q < bary.length; q++) {
          sum += weights[q] * basisAtPoints[q][i] * basisAtPoints[q][j]
        }
        M[i][j] = sum * vol
      }
    }
    return M
  }

  /**
   * Assembles the RHS for bubble projection.
   * @param {number} tIdx
   * @param {number} p
   * @param {function(!Array<number>): number} residualFn
   * @return {!Array<number>}
   */
  assembleBubbleRhs (tIdx, p, residualFn) {
    const exponents = this.#getBubbleExponents(p)
    if (exponents.length === 0) {
      return []
    }

    const { bary, weights } = tetrahedronQuadrature(this.quadratureOrder)
    const vol = this.mesh.getVolume(tIdx)

    const n = exponents.length
    const f = new Array(n).fill(0)

    for (let q = 0; q < bary.length; q++) {
      const basis = this.evaluateBubbleBasis(bary[q], p)
      const verts = this.mesh.tetrahedra[tIdx].map((i) => this.mesh.vertices[i])
      const pt = barycentricToCartesian(verts, bary[q])
      const r = residualFn(pt)
      for (let i = 0; i < n; i++) {
        f[i] += weights[q] * r * basis[i]
      }
    }
    return f.map((x) => x * vol)
  }

  /**
   * Solves for bubble coefficients.
   * @param {number} tIdx
   * @param {number} p
   * @param {function(!Array<number>): number} residualFn
   * @return {?Array<number>}
   */
  solveBubbleProjection (tIdx, p, residualFn) {
    if (p < 4) {
      return null
    }
    try {
      const M = this.assembleBubbleMass(tIdx, p)
      const f = this.assembleBubbleRhs(tIdx, p, residualFn)
      return luSolve(M, f)
    } catch (e) {
      this.onWarning({
        code: 'HOP_BUBBLE_SOLVE_FAILED',
        severity: 'warn',
        message:
          `HigherOrderProjection: bubble solve failed for tet ${tIdx}, p=${p}: ${e.message}`
      })
      return null
    }
  }

  /**
   * Evaluates the bubble correction at a point.
   * @param {number} tIdx
   * @param {number} p
   * @param {!Array<number>} coeffs
   * @param {!Array<number>} point
   * @return {number}
   */
  evaluateBubble (tIdx, p, coeffs, point) {
    if (!coeffs || coeffs.length === 0) {
      return 0
    }
    const bary = this.whitney.getBarycentric(tIdx, point)
    const basis = this.evaluateBubbleBasis(bary, p)
    return basis.reduce((sum, b, i) => sum + (coeffs[i] || 0) * b, 0)
  }

  /**
   * Generates Bernstein exponents for P^p on a tetrahedron.
   * Basis: {B_{i,j,k,l} = (p!/(i!j!k!l!)) lambda_0^i lambda_1^j lambda_2^k lambda_3^l
   * | i+j+k+l = p}.
   * This basis is well-conditioned and forms a partition of unity.
   * @param {number} p
   * @return {!Array<!Array<number>>}
   * @private
   */
  #getBernsteinExponents (p) {
    const exponents = []
    for (let i = 0; i <= p; i++) {
      for (let j = 0; j <= p - i; j++) {
        for (let k = 0; k <= p - i - j; k++) {
          const l = p - i - j - k
          exponents.push([i, j, k, l])
        }
      }
    }
    return exponents
  }

  /**
   * Computes the multinomial coefficient p!/(i!j!k!l!).
   * @param {number} p
   * @param {!Array<number>} exps
   * @return {number}
   * @private
   */
  static #multinomial (p, exps) {
    let result = 1
    let remaining = p
    for (const e of exps) {
      for (let i = 1; i <= e; i++) {
        result *= remaining
        result /= i
        remaining--
      }
    }
    return result
  }

  /**
   * Evaluates the P^p Bernstein basis at barycentric coordinates.
   * @param {!Array<number>} bary
   * @param {number} p
   * @return {!Array<number>}
   */
  evaluateBernsteinBasis (bary, p) {
    const exponents = this.#getBernsteinExponents(p)
    return exponents.map(([i, j, k, l]) => {
      const coeff = HigherOrderProjection.#multinomial(p, [i, j, k, l])
      return (
        coeff *
        HigherOrderProjection.#pow(bary[0], i) *
        HigherOrderProjection.#pow(bary[1], j) *
        HigherOrderProjection.#pow(bary[2], k) *
        HigherOrderProjection.#pow(bary[3], l)
      )
    })
  }

  /**
   * Computes the analytical mass matrix for the Bernstein basis on a tetrahedron.
   * M[(i,j,k,l), (i',j',k',l')] = vol * (p!)^2 * (i+i')! (j+j')! (k+k')! (l+l')!
   *   / (i! j! k! l! i'! j'! k'! l'! * (2p+3)!)
   * @param {number} p
   * @param {number} vol
   * @param {!Array<!Array<number>>} exponents
   * @return {!Array<!Array<number>>}
   * @private
   */
  #assembleBernsteinMass (p, vol, exponents) {
    const n = exponents.length
    const pFact = factorial(p)
    const denom = factorial(2 * p + 3)
    const M = zeros(n, n)
    for (let i = 0; i < n; i++) {
      const [a, b, c, d] = exponents[i]
      const factI = [
        factorial(a),
        factorial(b),
        factorial(c),
        factorial(d)
      ]
      for (let j = 0; j < n; j++) {
        const [a2, b2, c2, d2] = exponents[j]
        const num =
          factorial(a + a2) *
          factorial(b + b2) *
          factorial(c + c2) *
          factorial(d + d2)
        const den =
          factI[0] *
          factI[1] *
          factI[2] *
          factI[3] *
          factorial(a2) *
          factorial(b2) *
          factorial(c2) *
          factorial(d2)
        M[i][j] = (6 * vol * pFact * pFact * num) / (den * denom)
      }
    }
    return M
  }

  /**
   * Solves the L2 projection of u onto P^p(T) using the Bernstein basis.
   * The mass matrix is computed analytically; only the RHS uses quadrature.
   * @param {number} tIdx
   * @param {number} p
   * @param {function(!Array<number>): number} u
   * @return {!Array<number>} Coefficients of the Bernstein basis.
   */
  solveL2Projection (tIdx, p, u) {
    const exponents = this.#getBernsteinExponents(p)
    const n = exponents.length
    if (n === 0) {
      return []
    }

    const vol = this.mesh.getVolume(tIdx)
    const verts = this.mesh.tetrahedra[tIdx].map((i) => this.mesh.vertices[i])

    // Use composite quadrature for RHS to ensure enough sample points.
    const { bary, weights } =
      p >= 2
        ? compositeTetrahedronQuadrature(this.quadratureOrder)
        : tetrahedronQuadrature(this.quadratureOrder)
    const basisAtPoints = bary.map((lam) => this.evaluateBernsteinBasis(lam, p))
    const f = new Array(n).fill(0)

    for (let q = 0; q < bary.length; q++) {
      const pt = barycentricToCartesian(verts, bary[q])
      const uVal = u(pt)
      const bq = basisAtPoints[q]
      for (let i = 0; i < n; i++) {
        f[i] += weights[q] * uVal * bq[i]
      }
    }
    for (let i = 0; i < n; i++) {
      f[i] *= vol
    }

    const M = this.#assembleBernsteinMass(p, vol, exponents)
    try {
      return luSolve(M, f)
    } catch (e) {
      this.onWarning({
        code: 'HOP_L2_SOLVE_FAILED',
        severity: 'warn',
        message:
          `HigherOrderProjection: L2 solve failed for tet ${tIdx}, p=${p}: ${e.message}`
      })
      return []
    }
  }

  /**
   * Evaluates the L2 projection at a point.
   * @param {!Array<number>} coeffs
   * @param {!Array<number>} bary
   * @param {number} p
   * @return {number}
   */
  evaluateL2Projection (coeffs, bary, p) {
    if (!coeffs || coeffs.length === 0) {
      return 0
    }
    const basis = this.evaluateBernsteinBasis(bary, p)
    return basis.reduce((sum, b, i) => sum + (coeffs[i] || 0) * b, 0)
  }
}
