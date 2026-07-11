/**
 * Tetrahedral mesh data structure with topological connectivity, geometry,
 * and boundary extraction.
 *
 * This module provides the core Mesh class used throughout the BCDTPP library.
 * The Mesh is constructed once and treated as immutable for the duration of a
 * projection computation.  Mesh refinement operators (Alfeld split,
 * Worsey-Farin split) live in a separate {@link MeshRefinement} class.
 */

import { triangleArea, tetDeterminant, tetVolume } from './math_utils.js'
import { MeshValidationError } from './errors.js'

/**
 * Sorts three numbers in ascending order.
 * @param {number} a
 * @param {number} b
 * @param {number} c
 * @return {!Array<number>}
 */
function sort3 (a, b, c) {
  if (a > b) [a, b] = [b, a]
  if (b > c) [b, c] = [c, b]
  if (a > b) [a, b] = [b, a]
  return [a, b, c]
}

/**
 * Tetrahedral mesh data structure with topological connectivity and boundary
 * extraction.  This class is intentionally a *pure data structure*; mesh
 * refinement operators (Alfeld split, Worsey-Farin split) live in
 * {@link MeshRefinement}.
 */

export class Mesh {
  /** @type {!Map<number, number>} */
  #edgeKeyToIndex
  /** @type {!Map<number, number>} */
  #faceKeyToIndex
  /** @type {!Array<!Array<number>>} */
  #tetEdgeSigns
  /** @type {!Array<!Array<number>>} */
  #tetFaceSigns
  /** @type {!Array<number>} */
  #volumeCache
  /** @type {!Array<number>} */
  #faceAreaCache
  /** @type {!Array<!Array<number>>} */
  #faceNormalCache
  /** @type {!Array<!Array<number>>} */
  #faceBarycenterCache
  /** @type {!Array<!Array<number>>} */
  #tetBarycenterCache

  /**
   * @param {!Array<!Array<number>>} vertices - Vertex coordinates.
   * @param {!Array<!Array<number>>} tetrahedra - Tetrahedron vertex indices.
   */
  constructor (vertices, tetrahedra) {
    Mesh.#validateInput(vertices, tetrahedra)
    this.vertices = vertices
    this.tetrahedra = tetrahedra
    this.originalVertexCount = vertices.length
    this.vertexCount = vertices.length
    this.tetrahedronCount = tetrahedra.length

    this.faces = []
    this.edges = []
    this.boundaryFaces = []
    this.boundaryEdges = []
    this.boundaryNodes = new Set()

    this.faceToTets = []
    this.edgeToFaces = []

    this.vertexToTets = []
    this.vertexToBoundaryFaces = []
    this.vertexToEdges = []

    this.#edgeKeyToIndex = new Map()
    this.#faceKeyToIndex = new Map()
    this.#tetEdgeSigns = []
    this.#tetFaceSigns = []

    this.#volumeCache = []
    this.#faceAreaCache = []
    this.#faceNormalCache = []
    this.#faceBarycenterCache = []
    this.#tetBarycenterCache = []

    this.#buildTopology()
    this.#computeOrientationSigns()
    this.#precomputeGeometry()
  }

  /**
   * Validates mesh input before construction.
   * @param {!Array<!Array<number>>} vertices
   * @param {!Array<!Array<number>>} tetrahedra
   * @private
   */
  static #validateInput (vertices, tetrahedra) {
    if (!Array.isArray(vertices) || vertices.length === 0) {
      throw new MeshValidationError('vertices must be a non-empty array')
    }
    if (!Array.isArray(tetrahedra) || tetrahedra.length === 0) {
      throw new MeshValidationError('tetrahedra must be a non-empty array')
    }
    for (let i = 0; i < vertices.length; i++) {
      const v = vertices[i]
      if (!Array.isArray(v) || v.length !== 3 || !v.every((n) => typeof n === 'number' && Number.isFinite(n))) {
        throw new MeshValidationError(`Vertex ${i} must be an array of 3 finite numbers`)
      }
    }
    for (let i = 0; i < tetrahedra.length; i++) {
      const t = tetrahedra[i]
      if (!Array.isArray(t) || t.length !== 4 || !t.every((n) => typeof n === 'number' && Number.isInteger(n))) {
        throw new MeshValidationError(`Tetrahedron ${i} must be an array of 4 integer indices`)
      }
      for (let j = 0; j < 4; j++) {
        if (t[j] < 0 || t[j] >= vertices.length) {
          throw new MeshValidationError(
            `Tetrahedron ${i} contains out-of-bounds vertex index ${t[j]} (valid range: 0-${vertices.length - 1})`
          )
        }
      }
      if (new Set(t).size !== 4) {
        throw new MeshValidationError(`Tetrahedron ${i} contains duplicate vertex indices`)
      }

      // Orientation check: signed volume must be strictly positive.
      const v0 = vertices[t[0]]
      const v1 = vertices[t[1]]
      const v2 = vertices[t[2]]
      const v3 = vertices[t[3]]
      const signedVol = tetDeterminant(v0, v1, v2, v3)
      if (signedVol <= 1e-12) {
        throw new MeshValidationError(
          `Tetrahedron ${i} has non-positive or degenerate volume (signedVol=${signedVol}). ` +
            'Reorder vertices to satisfy the right-hand rule, or check for collinear/coplanar vertices.'
        )
      }
    }
  }

  /** @private */
  #precomputeGeometry () {
    for (let tIdx = 0; tIdx < this.tetrahedronCount; tIdx++) {
      this.#volumeCache[tIdx] = this.#computeVolume(tIdx)
      this.#tetBarycenterCache[tIdx] = this.#computeTetrahedronBarycenter(tIdx)
    }
    for (let fIdx = 0; fIdx < this.faces.length; fIdx++) {
      this.#faceAreaCache[fIdx] = this.#computeFaceArea(fIdx)
      this.#faceNormalCache[fIdx] = this.#computeFaceOutwardNormal(fIdx)
      this.#faceBarycenterCache[fIdx] = this.#computeFaceBarycenter(fIdx)
    }
  }

  /** @private */
  #computeFaceArea (fIdx) {
    const f = this.faces[fIdx]
    return triangleArea(
      this.vertices[f[0]],
      this.vertices[f[1]],
      this.vertices[f[2]]
    )
  }

  /** @private */
  #computeFaceBarycenter (fIdx) {
    const f = this.faces[fIdx]
    const v = f.map((i) => this.vertices[i])
    return [
      (v[0][0] + v[1][0] + v[2][0]) / 3,
      (v[0][1] + v[1][1] + v[2][1]) / 3,
      (v[0][2] + v[1][2] + v[2][2]) / 3
    ]
  }

  /** @private */
  #computeTetrahedronBarycenter (tIdx) {
    const t = this.tetrahedra[tIdx]
    const v = t.map((i) => this.vertices[i])
    return [
      (v[0][0] + v[1][0] + v[2][0] + v[3][0]) / 4,
      (v[0][1] + v[1][1] + v[2][1] + v[3][1]) / 4,
      (v[0][2] + v[1][2] + v[2][2] + v[3][2]) / 4
    ]
  }

  /** @private */
  #computeVolume (tIdx) {
    const tet = this.tetrahedra[tIdx]
    const v = tet.map((i) => this.vertices[i])
    return tetVolume(v[0], v[1], v[2], v[3])
  }

  /** @private */
  #computeFaceOutwardNormal (fIdx) {
    const f = this.faces[fIdx]
    const v0 = this.vertices[f[0]]
    const v1 = this.vertices[f[1]]
    const v2 = this.vertices[f[2]]
    const e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]]
    const e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]]
    const n = [
      e1[1] * e2[2] - e1[2] * e2[1],
      e1[2] * e2[0] - e1[0] * e2[2],
      e1[0] * e2[1] - e1[1] * e2[0]
    ]

    const tets = this.faceToTets[fIdx]
    if (tets.length === 1) {
      const tet = this.tetrahedra[tets[0]]
      const oppV = tet.find((v) => !f.includes(v))
      const oppP = this.vertices[oppV]
      const centroid = [
        (v0[0] + v1[0] + v2[0]) / 3,
        (v0[1] + v1[1] + v2[1]) / 3,
        (v0[2] + v1[2] + v2[2]) / 3
      ]
      const toOpp = [oppP[0] - centroid[0], oppP[1] - centroid[1], oppP[2] - centroid[2]]
      const dotProduct = n[0] * toOpp[0] + n[1] * toOpp[1] + n[2] * toOpp[2]
      if (dotProduct > 0) {
        n[0] = -n[0]
        n[1] = -n[1]
        n[2] = -n[2]
      }
    }

    const len = Math.sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2])
    if (len < 1e-12) {
      return [0, 0, 1]
    }
    return [n[0] / len, n[1] / len, n[2] / len]
  }

  /**
   * Returns the geometric area of a face.
   * @param {number} fIdx
   * @return {number}
   */
  getFaceArea (fIdx) {
    return this.#faceAreaCache[fIdx]
  }

  /**
   * Looks up the global edge index for an edge key.
   * @param {number} edgeKey
   * @return {number|undefined}
   */
  getEdgeIndex (edgeKey) {
    return this.#edgeKeyToIndex.get(edgeKey)
  }

  /**
   * Looks up the global face index for a face key.
   * @param {number} faceKey
   * @return {number|undefined}
   */
  getFaceIndex (faceKey) {
    return this.#faceKeyToIndex.get(faceKey)
  }

  /**
   * Returns the edge orientation sign for a local edge within a tetrahedron.
   * @param {number} tIdx
   * @param {number} e - Local edge index (0-5).
   * @return {number}
   */
  getTetEdgeSign (tIdx, e) {
    return this.#tetEdgeSigns[tIdx][e]
  }

  /**
   * Returns all vertex coordinates.
   * @return {!Array<!Array<number>>}
   */
  getVertices () {
    return this.vertices
  }

  /**
   * Returns all tetrahedra as vertex index arrays.
   * @return {!Array<!Array<number>>}
   */
  getTetrahedra () {
    return this.tetrahedra
  }

  /**
   * Returns all faces as vertex index arrays.
   * @return {!Array<!Array<number>>}
   */
  getFaces () {
    return this.faces
  }

  /**
   * Returns all edges as vertex index pairs.
   * @return {!Array<!Array<number>>}
   */
  getEdges () {
    return this.edges
  }

  /**
   * Returns the set of boundary node indices.
   * @return {!Set<number>}
   */
  getBoundaryNodes () {
    return this.boundaryNodes
  }

  /**
   * Returns the array of boundary face indices.
   * @return {!Array<number>}
   */
  getBoundaryFaces () {
    return this.boundaryFaces
  }

  /**
   * Returns the array of boundary edge indices.
   * @return {!Array<number>}
   */
  getBoundaryEdges () {
    return this.boundaryEdges
  }

  /**
   * Returns the original vertex count before refinement.
   * @return {number}
   */
  getOriginalVertexCount () {
    return this.originalVertexCount
  }

  /**
   * Returns the face orientation sign for a local face within a tetrahedron.
   * @param {number} tIdx
   * @param {number} f - Local face index (0-3).
   * @return {number}
   */
  getTetFaceSign (tIdx, f) {
    return this.#tetFaceSigns[tIdx][f]
  }

  /**
   * Computes a canonical integer key for an edge given two vertex indices.
   *
   * Safe while vertexCount < 2^26 (~67 million); beyond that the product
   * a * vertexCount may exceed Number.MAX_SAFE_INTEGER.
   *
   * @param {number} a
   * @param {number} b
   * @param {number} vertexCount
   * @return {number}
   */
  static computeEdgeKey (a, b, vertexCount) {
    return a < b ? a * vertexCount + b : b * vertexCount + a
  }

  /** @private */
  #buildTopology () {
    const faceMap = new Map()
    const edgeMap = new Map()
    const vc = this.originalVertexCount

    for (let tIdx = 0; tIdx < this.tetrahedronCount; tIdx++) {
      const tet = this.tetrahedra[tIdx]
      const localFaces = [
        [tet[1], tet[2], tet[3]],
        [tet[0], tet[2], tet[3]],
        [tet[0], tet[1], tet[3]],
        [tet[0], tet[1], tet[2]]
      ]

      for (const f of localFaces) {
        const key = this.#faceKey(f)
        if (!faceMap.has(key)) {
          faceMap.set(key, { verts: f, tets: [tIdx] })
        } else {
          faceMap.get(key).tets.push(tIdx)
        }
      }
    }

    let fIdx = 0
    for (const [key, data] of faceMap.entries()) {
      this.faces.push(data.verts)
      this.faceToTets.push(data.tets)
      this.#faceKeyToIndex.set(key, fIdx)

      const isBoundary = data.tets.length === 1
      if (isBoundary) {
        this.boundaryFaces.push(fIdx)
        data.verts.forEach((v) => this.boundaryNodes.add(v))
      }

      const fEdges = [
        [data.verts[0], data.verts[1]],
        [data.verts[1], data.verts[2]],
        [data.verts[2], data.verts[0]]
      ]

      for (const e of fEdges) {
        const eKey = Mesh.computeEdgeKey(e[0], e[1], vc)
        if (!edgeMap.has(eKey)) {
          edgeMap.set(eKey, { verts: e, isBoundary, faces: [fIdx] })
        } else {
          const eData = edgeMap.get(eKey)
          eData.faces.push(fIdx)
          if (isBoundary) {
            eData.isBoundary = true
          }
        }
      }
      fIdx++
    }

    let eIdx = 0
    for (const [key, data] of edgeMap.entries()) {
      this.edges.push(data.verts)
      this.edgeToFaces.push(data.faces)
      this.#edgeKeyToIndex.set(key, eIdx)
      if (data.isBoundary) {
        this.boundaryEdges.push(eIdx)
      }
      eIdx++
    }

    // Build adjacency lists.
    this.vertexToTets = Array.from({ length: this.vertexCount }, () => [])
    for (let tIdx = 0; tIdx < this.tetrahedronCount; tIdx++) {
      for (const v of this.tetrahedra[tIdx]) {
        this.vertexToTets[v].push(tIdx)
      }
    }

    this.vertexToBoundaryFaces = Array.from({ length: this.vertexCount }, () => [])
    for (const fIdx of this.boundaryFaces) {
      for (const v of this.faces[fIdx]) {
        this.vertexToBoundaryFaces[v].push(fIdx)
      }
    }

    this.vertexToEdges = Array.from({ length: this.vertexCount }, () => [])
    for (let eIdx = 0; eIdx < this.edges.length; eIdx++) {
      for (const v of this.edges[eIdx]) {
        this.vertexToEdges[v].push(eIdx)
      }
    }
  }

  /**
   * Computes orientation signs for edges and faces within each tetrahedron.
   * @private
   */
  #computeOrientationSigns () {
    this.#tetEdgeSigns = new Array(this.tetrahedronCount)
    this.#tetFaceSigns = new Array(this.tetrahedronCount)

    for (let tIdx = 0; tIdx < this.tetrahedronCount; tIdx++) {
      const tet = this.tetrahedra[tIdx]

      const edges = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]]
      this.#tetEdgeSigns[tIdx] = edges.map(([i, j]) => {
        const gEdge = [tet[i], tet[j]]
        const key = Mesh.computeEdgeKey(gEdge[0], gEdge[1], this.originalVertexCount)
        const gIdx = this.#edgeKeyToIndex.get(key)
        const stored = this.edges[gIdx]
        return stored[0] === gEdge[0] && stored[1] === gEdge[1] ? 1 : -1
      })

      const faces = [
        [tet[1], tet[2], tet[3]],
        [tet[0], tet[2], tet[3]],
        [tet[0], tet[1], tet[3]],
        [tet[0], tet[1], tet[2]]
      ]
      this.#tetFaceSigns[tIdx] = faces.map((f) => {
        const key = this.#faceKey(f)
        const gIdx = this.#faceKeyToIndex.get(key)
        const stored = this.faces[gIdx]
        return this.#facePermutationSign(stored, f)
      })
    }
  }

  /**
   * Returns +1 if local is an even permutation of stored, -1 if odd.
   * @param {!Array<number>} stored
   * @param {!Array<number>} local
   * @return {number}
   * @private
   */
  #facePermutationSign (stored, local) {
    const idx = local.map((v) => stored.indexOf(v))
    // Count inversions in the permutation (idx[0], idx[1], idx[2]) of (0,1,2).
    // Even parity → +1, odd parity → -1.
    let inversions = 0
    for (let i = 0; i < idx.length; i++) {
      for (let j = i + 1; j < idx.length; j++) {
        if (idx[i] > idx[j]) inversions++
      }
    }
    return inversions % 2 === 0 ? 1 : -1
  }

  /**
   * Returns the cached barycenter (centroid) of a triangle face.
   * @param {number} fIdx
   * @return {!Array<number>} [x, y, z]
   */
  getFaceBarycenter (fIdx) {
    return this.#faceBarycenterCache[fIdx]
  }

  /**
   * Returns the cached barycenter (centroid) of a tetrahedron.
   * @param {number} tIdx
   * @return {!Array<number>} [x, y, z]
   */
  getTetrahedronBarycenter (tIdx) {
    return this.#tetBarycenterCache[tIdx]
  }

  /**
   * Returns the set of tetrahedra incident on a vertex (the vertex star).
   * @param {number} vIdx
   * @return {!Array<number>} Array of tetrahedron indices.
   */
  getStar (vIdx) {
    return this.vertexToTets[vIdx] || []
  }

  /**
   * Returns the boundary faces incident on a vertex (the boundary star).
   * @param {number} vIdx
   * @return {!Array<number>} Array of boundary face indices.
   */
  getBoundaryStar (vIdx) {
    return this.vertexToBoundaryFaces[vIdx] || []
  }

  /**
   * Returns the edges incident on a vertex.
   * @param {number} vIdx
   * @return {!Array<number>} Array of edge indices.
   */
  getEdgeStar (vIdx) {
    return this.vertexToEdges[vIdx] || []
  }

  /**
   * Returns the cached volume of a tetrahedron.
   * @param {number} tIdx
   * @return {number} Volume (absolute value of signed determinant / 6).
   */
  getVolume (tIdx) {
    return this.#volumeCache[tIdx]
  }

  /**
   * Returns the four global face indices of a tetrahedron in local order
   * (opposite-vertex convention: face 0 is opposite vertex 0, etc.).
   * @param {number} tIdx
   * @return {!Array<number>} Array of four face indices.
   */
  getTetrahedronFaces (tIdx) {
    const tet = this.tetrahedra[tIdx]
    const local = [
      [tet[1], tet[2], tet[3]],
      [tet[0], tet[2], tet[3]],
      [tet[0], tet[1], tet[3]],
      [tet[0], tet[1], tet[2]]
    ]
    return local.map((lf) => {
      const key = this.#faceKey(lf)
      return this.#faceKeyToIndex.get(key)
    })
  }

  /**
   * Integer face key using sorted vertex ids.
   *
   * Safe while originalVertexCount < 2^17 (~131 thousand); beyond that the
   * quadratic term s0 * vc^2 may exceed Number.MAX_SAFE_INTEGER.
   *
   * @param {!Array<number>} f
   * @return {number}
   * @private
   */
  #faceKey (f) {
    const [s0, s1, s2] = sort3(f[0], f[1], f[2])
    const vc = this.originalVertexCount
    return ((s0 * vc) + s1) * vc + s2
  }

  /**
   * Computes the outward-pointing unit normal for a boundary face.
   * For interior faces the orientation is arbitrary.
   * @param {number} fIdx
   * @return {!Array<number>}
   */
  getFaceOutwardNormal (fIdx) {
    return this.#faceNormalCache[fIdx]
  }
}
