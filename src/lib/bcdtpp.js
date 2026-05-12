/**
 * BCDTPP: Bounded, Commuting, Discrete-trace Preserving Projections.
 *
 * Implements the de Rham projection operators Pi^l for l = 0,1,2,3 on
 * tetrahedral meshes with boundary-aware trace preservation.
 */

import {
  numericalGradient
} from './math_utils.js'
import { ProjectionError } from './errors.js'
import { BoundaryWeightComputer } from './boundary_weight_computer.js'
import { PointLocator } from './point_locator.js'
import { HigherOrderProjection } from './higher_order_projection.js'
import { MeshRefinement } from './mesh_refinement.js'
import { H1Projector } from './projectors/h1_projector.js'
import { HcurlProjector } from './projectors/hcurl_projector.js'
import { HdivProjector } from './projectors/hdiv_projector.js'
import { L2Projector } from './projectors/l2_projector.js'

/**
 * BCDTPP: Bounded, Commuting, Discrete-trace Preserving Projections.
 *
 * Implements the de Rham projection operators Pi^l for l = 0,1,2,3 on
 * tetrahedral meshes with boundary-aware trace preservation.
 *
 * **Coupling note:** This class accesses mesh data through the {@link Mesh}
 * public API (getters for vertices, faces, edges, boundary flags, orientation
 * signs, etc.).  Swapping in a different mesh implementation requires only that
 * the new class implements the same getter interface.
 */
export class Bcdtpp {
  /** @type {!Mesh} */
  #mesh
  /** @type {!Whitney} */
  #whitney
  /** @type {number} */
  #quadratureOrder
  /** @type {(msg: string) => void} */
  #onWarning
  /** @type {PointLocator|null} */
  #pointLocator
  /** @type {!MeshRefinement} */
  #meshRefinement
  /** @type {!BoundaryWeightComputer} */
  #boundaryWeightComputer
  /** @type {!HigherOrderProjection} */
  #higherOrderProjector
  /** @type {!H1Projector} */
  #h1Projector
  /** @type {!HcurlProjector} */
  #hcurlProjector
  /** @type {!HdivProjector} */
  #hdivProjector
  /** @type {!L2Projector} */
  #l2Projector
  /** @type {!Map<number, {nodeMap: !Array<number>, psi: !Array<number>}>} */
  #vertexBoundaryData
  /** @type {!Set<number>} */
  #boundaryEdgeSet
  /** @type {!Set<number>} */
  #boundaryFaceSet

  /**
   * @param {!Mesh} mesh
   * @param {!Whitney} whitney
   * @param {!Object=} options
   * @param {number=} options.quadratureOrder - Quadrature order for integration (default 3).
   */
  constructor (mesh, whitney, options = {}) {
    this.#mesh = mesh
    this.#whitney = whitney
    this.#quadratureOrder = options.quadratureOrder || 3
    this.#onWarning =
      options.onWarning || ((ctx) => console.warn(ctx.message ?? ctx))
    this.#vertexBoundaryData = new Map()
    this.#pointLocator = null

    this.#meshRefinement = new MeshRefinement(this.#mesh)
    this.#boundaryWeightComputer = new BoundaryWeightComputer(
      this.#mesh,
      this.#meshRefinement,
      this.#onWarning
    )

    this.#boundaryEdgeSet = new Set(this.#mesh.getBoundaryEdges())
    this.#boundaryFaceSet = new Set(this.#mesh.getBoundaryFaces())

    this.#higherOrderProjector = new HigherOrderProjection(
      this.#mesh,
      this.#whitney,
      this.#quadratureOrder,
      this.#onWarning
    )

    this.#h1Projector = new H1Projector(this.#mesh, this.#whitney, this.#meshRefinement)
    this.#hcurlProjector = new HcurlProjector(this.#mesh, this.#whitney, this.#quadratureOrder)
    this.#hdivProjector = new HdivProjector(this.#mesh, this.#whitney, this.#quadratureOrder)
    this.#l2Projector = new L2Projector(this.#mesh, this.#whitney, this.#quadratureOrder)

    this.#validateMesh()
  }

  /** @private */
  #validateMesh () {
    let degenerateCount = 0
    for (let tIdx = 0; tIdx < this.#mesh.tetrahedronCount; tIdx++) {
      const vol = this.#mesh.getVolume(tIdx)
      if (vol < 1e-12) {
        degenerateCount++
      }
    }
    if (degenerateCount > 0) {
      this.#onWarning({
        code: 'BCDTPP_DEGENERATE_MESH',
        severity: 'warn',
        message:
          `Bcdtpp: mesh contains ${degenerateCount} degenerate or ` +
          'near-degenerate tetrahedra. Projections may fail.'
      })
    }
  }

  /** Quadrature order used for integrations. */
  get quadratureOrder () {
    return this.#quadratureOrder
  }

  buildPointLocator () {
    this.#pointLocator = new PointLocator(this.#mesh)
  }

  /** @private */
  #validateTetIdx (tIdx) {
    if (typeof tIdx !== 'number' || !Number.isInteger(tIdx)) {
      throw new ProjectionError(`tIdx must be an integer, got ${tIdx}`)
    }
    if (tIdx < 0 || tIdx >= this.#mesh.tetrahedronCount) {
      throw new ProjectionError(
        `tIdx=${tIdx} out of range [0, ${this.#mesh.tetrahedronCount - 1}]`
      )
    }
  }

  /** @private */
  static #validatePoint (point) {
    if (!Array.isArray(point) || point.length !== 3 ||
        !point.every((n) => typeof n === 'number' && Number.isFinite(n))) {
      throw new ProjectionError(
        `point must be an array of 3 finite numbers, got ${JSON.stringify(point)}`
      )
    }
  }

  /** Section 6.3.1: Construction of lowest-order vertex weights. */
  computeBoundaryWeights () {
    if (this.#meshRefinement.alfeldTriangles.length === 0) {
      this.#meshRefinement.computeWorseyFarinSplit()
    }
    const weights = this.#boundaryWeightComputer.compute()
    this.#vertexBoundaryData = weights.vertexBoundaryData
  }

  /**
   * H1 projection (l=0) of a scalar field.
   * @param {function(!Array<number>): number} u
   * @param {!Array<number>} point
   * @param {number} tIdx
   * @return {number}
   */
  projectH1 (u, point, tIdx) {
    Bcdtpp.#validatePoint(point)
    this.#validateTetIdx(tIdx)
    return this.#h1Projector.project(u, point, tIdx, this.#vertexBoundaryData)
  }

  /**
   * H(curl) projection (l=1).
   * Boundary edges use exact trace DoFs ∫_e u·t ds.
   * Interior edges use midpoint evaluation.
   * @param {function(!Array<number>): (number|!Array<number>)} u
   * @param {!Array<number>} point
   * @param {number} tIdx
   * @return {!Array<number>}
   */
  projectHcurl (u, point, tIdx) {
    Bcdtpp.#validatePoint(point)
    this.#validateTetIdx(tIdx)
    return this.#hcurlProjector.project(u, point, tIdx, this.#boundaryEdgeSet)
  }

  /**
   * H(div) projection (l=2).
   * Boundary faces use exact trace DoFs ∫_f u·n dA.
   * Interior faces use barycenter evaluation.
   * @param {function(!Array<number>): (number|!Array<number>)} u
   * @param {!Array<number>} point
   * @param {number} tIdx
   * @return {!Array<number>}
   */
  projectHdiv (u, point, tIdx) {
    Bcdtpp.#validatePoint(point)
    this.#validateTetIdx(tIdx)
    return this.#hdivProjector.project(u, point, tIdx, this.#boundaryFaceSet)
  }

  /**
   * L2 projection (l=3).
   * @param {function(!Array<number>): number} u
   * @param {number} tIdx
   * @return {number}
   */
  projectL2 (u, tIdx) {
    this.#validateTetIdx(tIdx)
    return this.#l2Projector.project(u, tIdx)
  }

  /**
   * Higher-order projection.
   * @param {function(!Array<number>): (number|!Array<number>)} u
   * @param {!Array<number>} point
   * @param {number} tIdx
   * @param {number} l - Form degree (0=H1 scalar, 1=Hcurl vector, 2=Hdiv vector, 3=L2 scalar).
   * @param {number} p
   * @return {(number|!Array<number>)}
   */
  projectHp (u, point, tIdx, l, p) {
    Bcdtpp.#validatePoint(point)
    this.#validateTetIdx(tIdx)
    if (p === 0) {
      const dispatch = {
        0: () => this.projectH1(u, point, tIdx),
        1: () => this.projectHcurl(u, point, tIdx),
        2: () => this.projectHdiv(u, point, tIdx),
        3: () => this.projectL2(u, tIdx)
      }
      if (!Object.hasOwn(dispatch, l)) {
        throw new ProjectionError(`Invalid form degree l=${l}`)
      }
      return dispatch[l]()
    }

    if (l === 0) {
      if (p === 1) {
        return this.projectH1(u, point, tIdx)
      }
      const bary = this.#whitney.getBarycentric(tIdx, point)
      const coeffs = this.#higherOrderProjector.solveL2Projection(
        tIdx,
        p,
        u
      )
      return this.#higherOrderProjector.evaluateL2Projection(
        coeffs,
        bary,
        p
      )
    }

    if (l === 3) {
      const bary = this.#whitney.getBarycentric(tIdx, point)
      const coeffs = this.#higherOrderProjector.solveL2Projection(
        tIdx,
        p,
        u
      )
      return this.#higherOrderProjector.evaluateL2Projection(
        coeffs,
        bary,
        p
      )
    }

    throw new ProjectionError(
      `Higher-order projections for l=${l}, p=${p} not yet implemented. ` +
        'Only l=0 and l=3 support p>0.'
    )
  }

  /**
   * Global projector Pi^l.
   * @param {function(!Array<number>): (number|!Array<number>)} u
   * @param {!Array<number>} point
   * @param {number} tIdx
   * @param {number} l
   * @param {number=} p
   * @return {(number|!Array<number>)}
   */
  project (u, point, tIdx, l, p = 0) {
    Bcdtpp.#validatePoint(point)
    this.#validateTetIdx(tIdx)
    if (p !== 0) {
      return this.projectHp(u, point, tIdx, l, p)
    }
    if (l === 3) {
      return this.projectL2(u, tIdx)
    }
    const boundaryData = this.extractBoundaryDofs(u, l)
    const partial = this.extendBoundary(boundaryData, point, tIdx, l)
    const samplePt = this.#mesh.getTetrahedronBarycenter(tIdx)
    const isScalar = typeof u(samplePt) === 'number'

    const w = (l === 0 || !isScalar) ? u : numericalGradient.bind(this, u)

    const v = (pt) => {
      const valW = w(pt)
      const valP = this.extendBoundary(boundaryData, pt, tIdx, l)
      if (typeof valW === 'number') {
        return valW - valP
      }
      return [
        valW[0] - valP[0],
        valW[1] - valP[1],
        valW[2] - valP[2]
      ]
    }
    const ring = this.projectRing(v, point, tIdx, l)
    if (typeof ring === 'number') {
      return ring + partial
    }
    return [
      ring[0] + partial[0],
      ring[1] + partial[1],
      ring[2] + partial[2]
    ]
  }

  /**
   * Extracts the boundary degrees of freedom for a given function.
   * @param {function(!Array<number>): (number|!Array<number>)} u
   * @param {number} l
   * @return {!Map<number, number>}
   */
  extractBoundaryDofs (u, l) {
    const result = new Map()
    if (l === 0) {
      for (const vIdx of this.#mesh.getBoundaryNodes()) {
        result.set(vIdx, this.#h1Projector.computeBoundaryIntegralH1(vIdx, u, this.#vertexBoundaryData))
      }
    } else if (l === 1) {
      for (const eIdx of this.#mesh.getBoundaryEdges()) {
        result.set(eIdx, this.#hcurlProjector.computeEdgeDof(u, eIdx))
      }
    } else if (l === 2) {
      for (const fIdx of this.#mesh.getBoundaryFaces()) {
        result.set(fIdx, this.#hdivProjector.computeFaceDof(u, fIdx))
      }
    }
    return result
  }

  /**
   * Pi_ring^l: interior projector with zero boundary trace.
   * @param {function(!Array<number>): (number|!Array<number>)} u
   * @param {!Array<number>} point
   * @param {number} tIdx
   * @param {number} l
   * @return {(number|!Array<number>)}
   */
  projectRing (u, point, tIdx, l) {
    Bcdtpp.#validatePoint(point)
    this.#validateTetIdx(tIdx)
    const dispatch = {
      0: () => this.#h1Projector.projectRing(u, point, tIdx),
      1: () => this.#hcurlProjector.projectRing(u, point, tIdx, this.#boundaryEdgeSet),
      2: () => this.#hdivProjector.projectRing(u, point, tIdx, this.#boundaryFaceSet),
      3: () => this.#l2Projector.project(u, tIdx)
    }
    if (!Object.hasOwn(dispatch, l)) {
      throw new ProjectionError(`Invalid form degree l=${l}`)
    }
    return dispatch[l]()
  }

  /**
   * E^l: discrete extension operator.
   * @param {!Map<number, number>} boundaryData
   * @param {!Array<number>} point
   * @param {number} tIdx
   * @param {number} l
   * @return {(number|!Array<number>)}
   */
  extendBoundary (boundaryData, point, tIdx, l) {
    Bcdtpp.#validatePoint(point)
    this.#validateTetIdx(tIdx)
    const dispatch = {
      0: () => this.#h1Projector.extendBoundary(boundaryData, point, tIdx),
      1: () => this.#hcurlProjector.extendBoundary(boundaryData, point, tIdx, this.#boundaryEdgeSet),
      2: () => this.#hdivProjector.extendBoundary(boundaryData, point, tIdx, this.#boundaryFaceSet)
    }
    if (!Object.hasOwn(dispatch, l)) {
      throw new ProjectionError(
        `Discrete extension not defined for form degree l=${l}`
      )
    }
    return dispatch[l]()
  }

  /**
   * Pi_partial^l: boundary correction part of the projection.
   * @param {function(!Array<number>): (number|!Array<number>)} u
   * @param {!Array<number>} point
   * @param {number} tIdx
   * @param {number} l
   * @return {(number|!Array<number>)}
   */
  projectPartial (u, point, tIdx, l) {
    Bcdtpp.#validatePoint(point)
    this.#validateTetIdx(tIdx)
    const boundaryData = this.extractBoundaryDofs(u, l)
    return this.extendBoundary(boundaryData, point, tIdx, l)
  }

  /**
   * Finds the tetrahedron containing the point and projects.
   * @param {function(!Array<number>): (number|!Array<number>)} u
   * @param {!Array<number>} point
   * @param {number=} l
   * @param {number=} p
   * @return {{value: (number|!Array<number>), tIdx: number, bary: !Array<number>}}
   */
  projectAtPoint (u, point, l = 0, p = 0) {
    Bcdtpp.#validatePoint(point)
    if (!this.#pointLocator) {
      this.buildPointLocator()
    }
    const found = this.#pointLocator.findTetrahedron(point)
    if (!found) {
      throw new ProjectionError(
        `Point [${point.join(', ')}] not found in any tetrahedron`
      )
    }
    const { tIdx, bary } = found
    return {
      value: this.projectHp(u, point, tIdx, l, p),
      tIdx,
      bary
    }
  }
}
