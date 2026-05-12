/**
 * Whitney forms and barycentric coordinate utilities on a tetrahedral mesh.
 */

import { cross, inverse3x3, tetDeterminant } from './math_utils.js'
import { ProjectionError } from './errors.js'

export class Whitney {
  /** @type {!Array<!Array<!Array<number>>|undefined>} */
  #gradLCache
  /** @type {!Array<number>} */
  #tetDetCache
  /** @type {!Array<!Array<!Array<number>>>} */
  #tetMatrixCache
  /** @type {!Array<!Array<!Array<number>>>} */
  #tetTinvCache

  /**
   * @param {!Mesh} mesh
   */
  constructor (mesh) {
    this.mesh = mesh
    this.#gradLCache = new Array(mesh.tetrahedronCount)
    this.#tetDetCache = new Array(mesh.tetrahedronCount)
    this.#tetMatrixCache = new Array(mesh.tetrahedronCount)
    this.#tetTinvCache = new Array(mesh.tetrahedronCount)

    for (let tIdx = 0; tIdx < mesh.tetrahedronCount; tIdx++) {
      const tet = mesh.tetrahedra[tIdx]
      const v = tet.map((i) => mesh.vertices[i])
      const T = [
        [v[0][0] - v[3][0], v[1][0] - v[3][0], v[2][0] - v[3][0]],
        [v[0][1] - v[3][1], v[1][1] - v[3][1], v[2][1] - v[3][1]],
        [v[0][2] - v[3][2], v[1][2] - v[3][2], v[2][2] - v[3][2]]
      ]

      const det = tetDeterminant(v[3], v[0], v[1], v[2])

      this.#tetDetCache[tIdx] = det
      this.#tetMatrixCache[tIdx] = T

      if (Math.abs(det) >= 1e-12) {
        const Tinv = inverse3x3(T)
        const g0 = [Tinv[0][0], Tinv[0][1], Tinv[0][2]]
        const g1 = [Tinv[1][0], Tinv[1][1], Tinv[1][2]]
        const g2 = [Tinv[2][0], Tinv[2][1], Tinv[2][2]]
        const g3 = [
          -g0[0] - g1[0] - g2[0],
          -g0[1] - g1[1] - g2[1],
          -g0[2] - g1[2] - g2[2]
        ]
        this.#gradLCache[tIdx] = [g0, g1, g2, g3]
        this.#tetTinvCache[tIdx] = Tinv
      }
    }
  }

  /**
   * Computes barycentric coordinates of a point with respect to a tetrahedron.
   * @param {number} tIdx
   * @param {!Array<number>} point
   * @return {!Array<number>}
   */
  getBarycentric (tIdx, point) {
    if (!Array.isArray(point) || point.length !== 3 ||
        !point.every((n) => typeof n === 'number' && Number.isFinite(n))) {
      throw new ProjectionError(
        `point must be an array of 3 finite numbers, got ${JSON.stringify(point)}`
      )
    }
    const v = this.mesh.tetrahedra[tIdx].map((i) => this.mesh.vertices[i])
    const T = this.#tetMatrixCache[tIdx] ?? [
      [v[0][0] - v[3][0], v[1][0] - v[3][0], v[2][0] - v[3][0]],
      [v[0][1] - v[3][1], v[1][1] - v[3][1], v[2][1] - v[3][1]],
      [v[0][2] - v[3][2], v[1][2] - v[3][2], v[2][2] - v[3][2]]
    ]

    const det = this.#tetDetCache[tIdx] ?? tetDeterminant(v[3], v[0], v[1], v[2])

    if (Math.abs(det) < 1e-12) {
      throw new ProjectionError('Degenerate tetrahedron')
    }

    const b = [point[0] - v[3][0], point[1] - v[3][1], point[2] - v[3][2]]
    const Tinv = this.#tetTinvCache[tIdx] ?? inverse3x3(T)
    const L = [
      Tinv[0][0] * b[0] + Tinv[0][1] * b[1] + Tinv[0][2] * b[2],
      Tinv[1][0] * b[0] + Tinv[1][1] * b[1] + Tinv[1][2] * b[2],
      Tinv[2][0] * b[0] + Tinv[2][1] * b[1] + Tinv[2][2] * b[2]
    ]
    return [L[0], L[1], L[2], 1 - L[0] - L[1] - L[2]]
  }

  /**
   * Computes gradients of the barycentric coordinate functions.
   * @param {number} tIdx
   * @return {!Array<!Array<number>>}
   */
  getGradBarycentric (tIdx) {
    const gradL = this.#gradLCache[tIdx]
    if (!gradL) {
      throw new ProjectionError(`Degenerate tetrahedron ${tIdx}`)
    }
    return gradL
  }

  /**
   * Edge-based Whitney 1-form (Nedelec) basis functions.
   *
   * For edge e = (i, j):
   *   phi_e = lambda_i * grad(lambda_j) - lambda_j * grad(lambda_i)
   *
   * @param {number} tIdx
   * @param {!Array<number>} bary
   * @return {!Array<!Array<number>>} Six edge basis vectors in order
   *   [0,1], [0,2], [0,3], [1,2], [1,3], [2,3].
   */
  getEdgeBasis (tIdx, bary) {
    if (!Array.isArray(bary) || bary.length !== 4 ||
        !bary.every((n) => typeof n === 'number' && Number.isFinite(n))) {
      throw new ProjectionError(
        `bary must be an array of 4 finite numbers, got ${JSON.stringify(bary)}`
      )
    }
    const gradL = this.getGradBarycentric(tIdx)
    const edges = [
      [0, 1],
      [0, 2],
      [0, 3],
      [1, 2],
      [1, 3],
      [2, 3]
    ]
    return edges.map(([i, j]) => {
      return [
        bary[i] * gradL[j][0] - bary[j] * gradL[i][0],
        bary[i] * gradL[j][1] - bary[j] * gradL[i][1],
        bary[i] * gradL[j][2] - bary[j] * gradL[i][2]
      ]
    })
  }

  /**
   * Face-based Raviart-Thomas (Whitney 2-form) basis functions.
   *
   * @param {number} tIdx
   * @param {!Array<number>} bary
   * @return {!Array<!Array<number>>} Four face basis vectors indexed by
   *   opposite vertex: 0, 1, 2, 3.
   */
  getFaceBasis (tIdx, bary) {
    if (!Array.isArray(bary) || bary.length !== 4 ||
        !bary.every((n) => typeof n === 'number' && Number.isFinite(n))) {
      throw new ProjectionError(
        `bary must be an array of 4 finite numbers, got ${JSON.stringify(bary)}`
      )
    }
    const gradL = this.getGradBarycentric(tIdx)
    // Vertices ordered so that the basis function points outward from the tet.
    const faces = [
      [1, 2, 3],
      [0, 3, 2],
      [0, 1, 3],
      [0, 2, 1]
    ]
    return faces.map(([j, k, l]) => {
      const c1 = cross(gradL[k], gradL[l])
      const c2 = cross(gradL[l], gradL[j])
      const c3 = cross(gradL[j], gradL[k])
      return [
        2 * (bary[j] * c1[0] + bary[k] * c2[0] + bary[l] * c3[0]),
        2 * (bary[j] * c1[1] + bary[k] * c2[1] + bary[l] * c3[1]),
        2 * (bary[j] * c1[2] + bary[k] * c2[2] + bary[l] * c3[2])
      ]
    })
  }
}
