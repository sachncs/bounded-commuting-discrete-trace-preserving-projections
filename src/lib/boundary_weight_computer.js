/**
 * Boundary weight computation for BCDTPP projections.
 *
 * Computes vertex patch weights, edge tangents/lengths, and face normals/areas
 * used by trace-preserving boundary DoFs.
 */

import { subtract, norm, triangleArea } from './math_utils.js'
import { LocalSolver } from './local_solver.js'

export class BoundaryWeightComputer {
  /**
   * @param {!Mesh} mesh
   * @param {!MeshRefinement} meshRefinement
   */
  constructor (mesh, meshRefinement, onWarning = console.warn) {
    this.mesh = mesh
    this.meshRefinement = meshRefinement
    this.onWarning =
      typeof onWarning === 'function'
        ? onWarning
        : (ctx) => console.warn(ctx.message ?? ctx)
  }

  /**
   * Computes all boundary weights.
   * @return {{
   *   vertexBoundaryData: !Map<number, {nodeMap: !Array<number>, invNodeMap: !Map<number, number>, psi: !Array<number>}>,
   *   edgeBoundaryData: !Map<number, {v0: number, v1: number, tangent: !Array<number>, length: number}>,
   *   faceBoundaryData: !Map<number, {normal: !Array<number>, area: number}>,
   * }}
   */
  compute () {
    const vertexBoundaryData = this.#computeVertexWeights()
    const edgeBoundaryData = this.#computeEdgeData()
    const faceBoundaryData = this.#computeFaceData()
    return { vertexBoundaryData, edgeBoundaryData, faceBoundaryData }
  }

  /** @private */
  #computeVertexWeights () {
    const zeta0Vertex = new Map()
    for (const vIdx of this.mesh.boundaryNodes) {
      try {
        const starFaces = this.mesh.getBoundaryStar(vIdx)
        const alfeldTris = this.meshRefinement.alfeldTriangles.filter((at) =>
          starFaces.includes(at.parentFaceIdx)
        )

        const triangles = alfeldTris.flatMap((at) => at.triangles)
        const starNodes = new Set(triangles.flat())
        const nodeMap = Array.from(starNodes)
        const invNodeMap = new Map(nodeMap.map((id, i) => [id, i]))

        const localTris = triangles.map((t) => t.map((v) => invNodeMap.get(v)))
        const localVerts = nodeMap.map((v) => this.mesh.vertices[v])

        const K = LocalSolver.assembleSurfaceStiffness(localVerts, localTris)
        const b = new Array(nodeMap.length).fill(0)

        const starArea = starFaces.reduce(
          (acc, fIdx) => acc + this.mesh.getFaceArea(fIdx),
          0
        )
        if (starArea < 1e-12) {
          this.onWarning({
            code: 'BWC_ZERO_STAR_AREA',
            severity: 'warn',
            message:
              `BoundaryWeightComputer: vertex ${vIdx} has zero or ` +
              `near-zero star area (${starArea}). Skipping weight computation.`
          })
          continue
        }
        const eta = 1.0 / starArea

        localTris.forEach((tri) => {
          const area = triangleArea(
            localVerts[tri[0]],
            localVerts[tri[1]],
            localVerts[tri[2]]
          )
          tri.forEach((nodeIdx) => {
            b[nodeIdx] += eta * (area / 3.0)
          })
        })

        const psi = LocalSolver.solveWithConstraint(K, b, this.onWarning)
        zeta0Vertex.set(vIdx, { nodeMap, invNodeMap, psi })
      } catch (e) {
        this.onWarning({
          code: 'BWC_VERTEX_FAILURE',
          severity: 'warn',
          message:
            'BoundaryWeightComputer: failed to compute weights for ' +
            `vertex ${vIdx}: ${e.message}`
        })
      }
    }
    return zeta0Vertex
  }

  /** @private */
  #computeEdgeData () {
    const zeta1Edge = new Map()
    for (const eIdx of this.mesh.boundaryEdges) {
      const e = this.mesh.edges[eIdx]
      const edgeVec = subtract(this.mesh.vertices[e[1]], this.mesh.vertices[e[0]])
      const edgeLen = norm(edgeVec)
      if (edgeLen < 1e-12) {
        continue
      }
      zeta1Edge.set(eIdx, {
        v0: e[0],
        v1: e[1],
        tangent: edgeVec.map((x) => x / edgeLen),
        length: edgeLen
      })
    }
    return zeta1Edge
  }

  /** @private */
  #computeFaceData () {
    const zeta2Face = new Map()
    for (const fIdx of this.mesh.boundaryFaces) {
      const normal = this.mesh.getFaceOutwardNormal(fIdx)
      const area = this.mesh.getFaceArea(fIdx)
      zeta2Face.set(fIdx, { normal, area })
    }
    return zeta2Face
  }
}
