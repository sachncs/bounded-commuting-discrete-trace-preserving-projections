/**
 * Lowest-order H1 (l=0) projector.
 *
 * Vertex-based scalar projection with boundary-aware trace DoFs.
 */

import { triangleQuadrature, barycentricToCartesian } from '../quadrature.js'
import { triangleArea } from '../math_utils.js'

export class H1Projector {
  /**
   * @param {!Mesh} mesh
   * @param {!Whitney} whitney
   * @param {!MeshRefinement} meshRefinement
   */
  constructor (mesh, whitney, meshRefinement) {
    this.mesh = mesh
    this.whitney = whitney
    this.meshRefinement = meshRefinement
  }

  /**
   * @param {function(!Array<number>): number} u
   * @param {!Array<number>} point
   * @param {number} tIdx
   * @param {!Map<number, {nodeMap: !Array<number>, psi: !Array<number>}>} vertexBoundaryData
   * @return {number}
   */
  project (u, point, tIdx, vertexBoundaryData) {
    const tet = this.mesh.getTetrahedra()[tIdx]
    const bary = this.whitney.getBarycentric(tIdx, point)

    let result = 0
    for (let i = 0; i < 4; i++) {
      const vIdx = tet[i]
      if (this.mesh.getBoundaryNodes().has(vIdx)) {
        const alpha = this.computeBoundaryIntegralH1(vIdx, u, vertexBoundaryData)
        result += alpha * bary[i]
      } else {
        result += u(this.mesh.getVertices()[vIdx]) * bary[i]
      }
    }
    return result
  }

  /**
   * @param {number} vIdx
   * @param {function(!Array<number>): number} u
   * @param {!Map<number, {nodeMap: !Array<number>, psi: !Array<number>}>} vertexBoundaryData
   * @return {number}
   */
  computeBoundaryIntegralH1 (vIdx, u, vertexBoundaryData) {
    const data = vertexBoundaryData.get(vIdx)
    if (!data) {
      return u(this.mesh.getVertices()[vIdx])
    }

    const { invNodeMap, psi } = data

    let integral = 0
    let zetaIntegral = 0
    const starFaces = this.mesh.getBoundaryStar(vIdx)

    starFaces.forEach((fIdx) => {
      const at = this.meshRefinement.faceToAlfeld.get(fIdx)
      if (!at) return
      at.triangles.forEach((tri) => {
        const verts = tri.map((i) => this.mesh.getVertices()[i])
        const area = triangleArea(verts[0], verts[1], verts[2])
        const localPsi = tri.map((v) => psi[invNodeMap.get(v)])
        const { bary, weights } = triangleQuadrature(2)
        let triIntegral = 0
        let triZetaIntegral = 0
        for (let q = 0; q < bary.length; q++) {
          const pt = barycentricToCartesian(verts, bary[q])
          const zeta =
            localPsi[0] * bary[q][0] +
            localPsi[1] * bary[q][1] +
            localPsi[2] * bary[q][2]
          triIntegral += weights[q] * u(pt) * zeta
          triZetaIntegral += weights[q] * zeta
        }
        integral += triIntegral * area
        zetaIntegral += triZetaIntegral * area
      })
    })

    if (Math.abs(zetaIntegral) < 1e-12) {
      return u(this.mesh.getVertices()[vIdx])
    }
    return integral / zetaIntegral
  }

  /**
   * Pi_ring^0: interior H1 projector with zero boundary trace.
   * @param {function(!Array<number>): number} u
   * @param {!Array<number>} point
   * @param {number} tIdx
   * @return {number}
   */
  projectRing (u, point, tIdx) {
    const tet = this.mesh.getTetrahedra()[tIdx]
    const bary = this.whitney.getBarycentric(tIdx, point)
    let result = 0
    for (let i = 0; i < 4; i++) {
      const vIdx = tet[i]
      if (!this.mesh.getBoundaryNodes().has(vIdx)) {
        result += u(this.mesh.getVertices()[vIdx]) * bary[i]
      }
    }
    return result
  }

  /**
   * E^0: discrete extension of vertex boundary data.
   * @param {!Map<number, number>} boundaryData
   * @param {!Array<number>} point
   * @param {number} tIdx
   * @return {number}
   */
  extendBoundary (boundaryData, point, tIdx) {
    const tet = this.mesh.getTetrahedra()[tIdx]
    const bary = this.whitney.getBarycentric(tIdx, point)
    let result = 0
    for (let i = 0; i < 4; i++) {
      const vIdx = tet[i]
      if (this.mesh.getBoundaryNodes().has(vIdx) && boundaryData.has(vIdx)) {
        result += boundaryData.get(vIdx) * bary[i]
      }
    }
    return result
  }
}
