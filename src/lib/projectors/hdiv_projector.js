/**
 * Lowest-order H(div) (l=2) projector.
 *
 * Face-based vector projection with boundary-aware normal trace DoFs.
 */

import { dot, cross, subtract, norm, numericalGradient, triangleArea } from '../math_utils.js'
import { triangleQuadrature, barycentricToCartesian } from '../quadrature.js'
import { ProjectionError } from '../errors.js'

/**
 * Lowest-order H(div) (l=2) face-based projector implementing Pi^2.
 *
 * Projects vector functions onto the Raviart-Thomas (Whitney 2-form) space.
 * Boundary faces use exact normal-flux degrees of freedom (∫_f u·n dA);
 * interior faces use barycenter evaluation of the normal component.
 */
export class HdivProjector {
  /**
   * @param {!Mesh} mesh
   * @param {!Whitney} whitney
   * @param {number} quadratureOrder
   */
  constructor (mesh, whitney, quadratureOrder) {
    this.mesh = mesh
    this.whitney = whitney
    this.quadratureOrder = quadratureOrder
  }

  /**
   * @param {function(!Array<number>): (number|!Array<number>)} u
   * @param {!Array<number>} point
   * @param {number} tIdx
   * @param {!Set<number>} boundaryFaceSet
   * @return {!Array<number>}
   */
  project (u, point, tIdx, boundaryFaceSet) {
    const bary = this.whitney.getBarycentric(tIdx, point)
    const faceBasis = this.whitney.getFaceBasis(tIdx, bary)

    const tFaces = this.mesh.getTetrahedronFaces(tIdx)
    const result = [0, 0, 0]

    const samplePt = this.mesh.getTetrahedronBarycenter(tIdx)
    const isScalar = typeof u(samplePt) === 'number'

    for (let f = 0; f < 4; f++) {
      const fIdx = tFaces[f]

      let coefficient
      if (boundaryFaceSet.has(fIdx)) {
        coefficient = this.computeFaceDof(u, fIdx)
      } else {
        coefficient = this.computeInteriorFaceCoeff(u, tIdx, f, isScalar)
      }

      result[0] += coefficient * faceBasis[f][0]
      result[1] += coefficient * faceBasis[f][1]
      result[2] += coefficient * faceBasis[f][2]
    }

    return result
  }

  /**
   * Computes the exact face DoF for H(div): ∫_f u·n dA.
   * For scalar u, this integrates grad(u)·n over the face.
   * @param {function(!Array<number>): (number|!Array<number>)} u
   * @param {number} fIdx
   * @return {number}
   */
  computeFaceDof (u, fIdx) {
    const f = this.mesh.getFaces()[fIdx]
    const verts = f.map((i) => this.mesh.getVertices()[i])
    const normal = this.mesh.getFaceOutwardNormal(fIdx)
    const area = triangleArea(verts[0], verts[1], verts[2])
    const { bary, weights } = triangleQuadrature(this.quadratureOrder)

    const isScalar = typeof u(verts[0]) === 'number'
    let integral = 0
    for (let q = 0; q < bary.length; q++) {
      const pt = barycentricToCartesian(verts, bary[q])
      if (isScalar) {
        const grad = numericalGradient(u, pt)
        integral += weights[q] * dot(grad, normal)
      } else {
        integral += weights[q] * dot(u(pt), normal)
      }
    }
    return integral * area
  }

  /**
   * Computes the interior face coefficient for H(div) in local outward normal.
   * @param {function(!Array<number>): (number|!Array<number>)} u
   * @param {number} tIdx
   * @param {number} f
   * @param {boolean} isScalar
   * @return {number}
   */
  computeInteriorFaceCoeff (u, tIdx, f, isScalar) {
    const tet = this.mesh.getTetrahedra()[tIdx]
    const fIdx = this.mesh.getTetrahedronFaces(tIdx)[f]
    const faceBary = this.mesh.getFaceBarycenter(fIdx)
    const localFaces = [
      [tet[1], tet[2], tet[3]],
      [tet[0], tet[2], tet[3]],
      [tet[0], tet[1], tet[3]],
      [tet[0], tet[1], tet[2]]
    ]
    const lf = localFaces[f]
    const verts = lf.map((v) => this.mesh.getVertices()[v])
    const e1 = subtract(verts[1], verts[0])
    const e2 = subtract(verts[2], verts[0])
    const nRaw = cross(e1, e2)
    const area = 0.5 * norm(nRaw)
    if (area < 1e-12) {
      throw new ProjectionError(
        `Degenerate face in tetrahedron ${tIdx}: face area=${area}`
      )
    }
    const normal = nRaw.map((c) => c / (2 * area))
    const oppV = tet.find((v) => !lf.includes(v))
    const oppP = this.mesh.getVertices()[oppV]
    const centroid = [
      (verts[0][0] + verts[1][0] + verts[2][0]) / 3,
      (verts[0][1] + verts[1][1] + verts[2][1]) / 3,
      (verts[0][2] + verts[1][2] + verts[2][2]) / 3
    ]
    const toOpp = subtract(oppP, centroid)
    if (dot(normal, toOpp) > 0) {
      normal[0] *= -1
      normal[1] *= -1
      normal[2] *= -1
    }
    if (isScalar) {
      return dot(numericalGradient(u, faceBary), normal) * area
    }
    return dot(u(faceBary), normal) * area
  }

  /**
   * Pi_ring^2: interior H(div) projector with zero boundary trace.
   * @param {function(!Array<number>): (number|!Array<number>)} u
   * @param {!Array<number>} point
   * @param {number} tIdx
   * @param {!Set<number>} boundaryFaceSet
   * @return {!Array<number>}
   */
  projectRing (u, point, tIdx, boundaryFaceSet) {
    const bary = this.whitney.getBarycentric(tIdx, point)
    const faceBasis = this.whitney.getFaceBasis(tIdx, bary)
    const tFaces = this.mesh.getTetrahedronFaces(tIdx)
    const result = [0, 0, 0]
    const samplePt = this.mesh.getTetrahedronBarycenter(tIdx)
    const isScalar = typeof u(samplePt) === 'number'

    for (let f = 0; f < 4; f++) {
      const fIdx = tFaces[f]
      if (boundaryFaceSet.has(fIdx)) {
        continue
      }
      const coefficient = this.computeInteriorFaceCoeff(u, tIdx, f, isScalar)
      result[0] += coefficient * faceBasis[f][0]
      result[1] += coefficient * faceBasis[f][1]
      result[2] += coefficient * faceBasis[f][2]
    }
    return result
  }

  /**
   * E^2: discrete extension of face boundary data.
   * @param {!Map<number, number>} boundaryData
   * @param {!Array<number>} point
   * @param {number} tIdx
   * @param {!Set<number>} boundaryFaceSet
   * @return {!Array<number>}
   */
  extendBoundary (boundaryData, point, tIdx, boundaryFaceSet) {
    const bary = this.whitney.getBarycentric(tIdx, point)
    const faceBasis = this.whitney.getFaceBasis(tIdx, bary)
    const tFaces = this.mesh.getTetrahedronFaces(tIdx)
    const result = [0, 0, 0]

    for (let f = 0; f < 4; f++) {
      const fIdx = tFaces[f]
      if (!boundaryFaceSet.has(fIdx) || !boundaryData.has(fIdx)) {
        continue
      }
      const coefficient = boundaryData.get(fIdx)
      result[0] += coefficient * faceBasis[f][0]
      result[1] += coefficient * faceBasis[f][1]
      result[2] += coefficient * faceBasis[f][2]
    }
    return result
  }
}
