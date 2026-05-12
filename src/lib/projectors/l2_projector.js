/**
 * Lowest-order L2 (l=3) projector.
 *
 * Cell-based scalar projection (average over tetrahedron).
 */

import { integrateTetrahedron } from '../quadrature.js'
import { ProjectionError } from '../errors.js'

export class L2Projector {
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
   * @param {function(!Array<number>): number} u
   * @param {number} tIdx
   * @return {number}
   */
  project (u, tIdx) {
    const tet = this.mesh.getTetrahedra()[tIdx]
    const verts = tet.map((i) => this.mesh.getVertices()[i])
    const vol = this.mesh.getVolume(tIdx)
    if (vol < 1e-12) {
      throw new ProjectionError(
        `Cannot perform L2 projection on degenerate tetrahedron ${tIdx}`
      )
    }

    const integrand = (pt) => u(pt)
    const integral = integrateTetrahedron(verts, integrand, this.quadratureOrder)
    return integral / vol
  }
}
