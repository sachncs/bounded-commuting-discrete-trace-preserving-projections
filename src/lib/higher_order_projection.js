/**
 * Higher-Order Projection Framework (Section 7).
 *
 * For p >= 1, the projection Pi^l_p is built recursively from the lowest-order
 * projection Pi^l_0 by adding bubble corrections on Alfeld-split patches.
 */

import {zeros, luSolve} from './math_utils.js';
import {
  barycentricToCartesian,
  tetrahedronQuadrature,
} from './quadrature.js';

export class HigherOrderProjection {
  /**
   * @param {!Mesh} mesh
   * @param {!Whitney} whitney
   * @param {number=} quadratureOrder
   */
  constructor(mesh, whitney, quadratureOrder = 3, onWarning = console.warn) {
    this.mesh = mesh;
    this.whitney = whitney;
    this.quadratureOrder = quadratureOrder;
    this.onWarning = onWarning;
    this.bubbleExponentCache = new Map();
  }

  /**
   * Generates the scalar bubble basis for polynomial degree p on a tetrahedron.
   * The bubble space consists of polynomials vanishing on the boundary.
   * For degree p, the bubble space is b * P^{p-4} where b = prod(lambda_i).
   * @param {number} p
   * @return {!Array<!Array<number>>} List of exponent tuples [a,b,c,d]
   *   such that basis function = lambda0^a * lambda1^b * lambda2^c * lambda3^d * b.
   */
  #getBubbleExponents(p) {
    if (this.bubbleExponentCache.has(p)) {
      return this.bubbleExponentCache.get(p);
    }
    const exponents = [];
    if (p >= 4) {
      const m = p - 4;
      for (let a = 0; a <= m; a++) {
        for (let b = 0; b <= m - a; b++) {
          for (let c = 0; c <= m - a - b; c++) {
            const d = m - a - b - c;
            // Exponents for lambda_i include the +1 from b.
            exponents.push([a + 1, b + 1, c + 1, d + 1]);
          }
        }
      }
    }
    this.bubbleExponentCache.set(p, exponents);
    return exponents;
  }

  /**
   * Polynomial power with explicit 0^0 = 1 convention.
   * @param {number} base
   * @param {number} exp
   * @return {number}
   * @private
   */
  static #pow(base, exp) {
    if (base === 0 && exp === 0) return 1;
    return Math.pow(base, exp);
  }

  /**
   * Evaluates the scalar bubble basis at a point given in barycentric coords.
   * @param {!Array<number>} bary
   * @param {number} p
   * @return {!Array<number>}
   */
  evaluateBubbleBasis(bary, p) {
    const exponents = this.#getBubbleExponents(p);
    if (exponents.length === 0) {
      return [];
    }
    return exponents.map(([a, b, c, d]) =>
      HigherOrderProjection.#pow(bary[0], a) *
      HigherOrderProjection.#pow(bary[1], b) *
      HigherOrderProjection.#pow(bary[2], c) *
      HigherOrderProjection.#pow(bary[3], d),
    );
  }

  /**
   * Assembles the local mass matrix for the bubble space.
   * @param {number} tIdx
   * @param {number} p
   * @return {!Array<!Array<number>>}
   */
  assembleBubbleMass(tIdx, p) {
    const exponents = this.#getBubbleExponents(p);
    if (exponents.length === 0) {
      return [];
    }

    const {bary, weights} = tetrahedronQuadrature(this.quadratureOrder);
    const vol = this.mesh.getVolume(tIdx);

    // Precompute basis values at quadrature points.
    const basisAtPoints = bary.map((lam) => this.evaluateBubbleBasis(lam, p));
    const n = basisAtPoints[0].length;
    const M = zeros(n, n);

    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        let sum = 0;
        for (let q = 0; q < bary.length; q++) {
          sum += weights[q] * basisAtPoints[q][i] * basisAtPoints[q][j];
        }
        M[i][j] = sum * vol;
      }
    }
    return M;
  }

  /**
   * Assembles the RHS for bubble projection.
   * @param {number} tIdx
   * @param {number} p
   * @param {function(!Array<number>): number} residualFn
   * @return {!Array<number>}
   */
  assembleBubbleRhs(tIdx, p, residualFn) {
    const exponents = this.#getBubbleExponents(p);
    if (exponents.length === 0) {
      return [];
    }

    const {bary, weights} = tetrahedronQuadrature(this.quadratureOrder);
    const vol = this.mesh.getVolume(tIdx);

    const n = exponents.length;
    const f = new Array(n).fill(0);

    for (let q = 0; q < bary.length; q++) {
      const basis = this.evaluateBubbleBasis(bary[q], p);
      const verts = this.mesh.tetrahedra[tIdx].map((i) => this.mesh.vertices[i]);
      const pt = barycentricToCartesian(verts, bary[q]);
      const r = residualFn(pt);
      for (let i = 0; i < n; i++) {
        f[i] += weights[q] * r * basis[i];
      }
    }
    return f.map((x) => x * vol);
  }

  /**
   * Solves for bubble coefficients.
   * @param {number} tIdx
   * @param {number} p
   * @param {function(!Array<number>): number} residualFn
   * @return {?Array<number>}
   */
  solveBubbleProjection(tIdx, p, residualFn) {
    if (p < 4) {
      return null;
    }
    try {
      const M = this.assembleBubbleMass(tIdx, p);
      const f = this.assembleBubbleRhs(tIdx, p, residualFn);
      return luSolve(M, f);
    } catch (e) {
      this.onWarning(
        `HigherOrderProjection: bubble solve failed for tet ${tIdx}, p=${p}: ${e.message}`,
      );
      return null;
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
  evaluateBubble(tIdx, p, coeffs, point) {
    if (!coeffs || coeffs.length === 0) {
      return 0;
    }
    const bary = this.whitney.getBarycentric(tIdx, point);
    const basis = this.evaluateBubbleBasis(bary, p);
    return basis.reduce((sum, b, i) => sum + (coeffs[i] || 0) * b, 0);
  }

  /**
   * Generates monomial exponents for P^p on a tetrahedron.
   * Basis: {lambda_0^a lambda_1^b lambda_2^c | a+b+c <= p}.
   * Lambda_3 is omitted because lambda_3 = 1 - lambda_0 - lambda_1 - lambda_2,
   * so the three independent barycentric coordinates fully describe P^p on T.
   * @param {number} p
   * @return {!Array<!Array<number>>}
   * @private
   */
  #getMonomialExponents(p) {
    const exponents = [];
    for (let a = 0; a <= p; a++) {
      for (let b = 0; b <= p - a; b++) {
        for (let c = 0; c <= p - a - b; c++) {
          exponents.push([a, b, c]);
        }
      }
    }
    return exponents;
  }

  /**
   * Evaluates the P^p monomial basis at barycentric coordinates.
   * @param {!Array<number>} bary
   * @param {number} p
   * @return {!Array<number>}
   */
  evaluateMonomialBasis(bary, p) {
    const exponents = this.#getMonomialExponents(p);
    const l0 = bary[0];
    const l1 = bary[1];
    const l2 = bary[2];
    return exponents.map(([a, b, c]) =>
      HigherOrderProjection.#pow(l0, a) *
      HigherOrderProjection.#pow(l1, b) *
      HigherOrderProjection.#pow(l2, c),
    );
  }

  /**
   * Solves the L2 projection of u onto P^p(T).
   * @param {number} tIdx
   * @param {number} p
   * @param {function(!Array<number>): number} u
   * @return {!Array<number>} Coefficients of the monomial basis.
   */
  solveL2Projection(tIdx, p, u) {
    const exponents = this.#getMonomialExponents(p);
    const n = exponents.length;
    if (n === 0) {
      return [];
    }

    const {bary, weights} = tetrahedronQuadrature(this.quadratureOrder);
    const vol = this.mesh.getVolume(tIdx);
    const verts = this.mesh.tetrahedra[tIdx].map((i) => this.mesh.vertices[i]);

    const basisAtPoints = bary.map((lam) => this.evaluateMonomialBasis(lam, p));
    const M = zeros(n, n);
    const f = new Array(n).fill(0);

    for (let q = 0; q < bary.length; q++) {
      const pt = barycentricToCartesian(verts, bary[q]);
      const uVal = u(pt);
      const bq = basisAtPoints[q];
      for (let i = 0; i < n; i++) {
        f[i] += weights[q] * uVal * bq[i];
        for (let j = 0; j < n; j++) {
          M[i][j] += weights[q] * bq[i] * bq[j];
        }
      }
    }

    for (let i = 0; i < n; i++) {
      f[i] *= vol;
      for (let j = 0; j < n; j++) {
        M[i][j] *= vol;
      }
    }

    return luSolve(M, f);
  }

  /**
   * Evaluates the L2 projection at a point.
   * @param {!Array<number>} coeffs
   * @param {!Array<number>} bary
   * @param {number} p
   * @return {number}
   */
  evaluateL2Projection(coeffs, bary, p) {
    if (!coeffs || coeffs.length === 0) {
      return 0;
    }
    const basis = this.evaluateMonomialBasis(bary, p);
    return basis.reduce((sum, b, i) => sum + (coeffs[i] || 0) * b, 0);
  }
}
