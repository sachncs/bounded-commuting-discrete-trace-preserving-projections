/**
 * Sorts three numbers in ascending order.
 * @param {number} a
 * @param {number} b
 * @param {number} c
 * @return {!Array<number>}
 */
function sort3(a, b, c) {
  if (a > b) [a, b] = [b, a];
  if (b > c) [b, c] = [c, b];
  if (a > b) [a, b] = [b, a];
  return [a, b, c];
}

/**
 * Tetrahedral mesh with topological connectivity, boundary extraction,
 * Alfeld face splitting, and Worsey-Farin tetrahedron splitting.
 */

export class Mesh {
  /**
   * @param {!Array<!Array<number>>} vertices - Vertex coordinates.
   * @param {!Array<!Array<number>>} tetrahedra - Tetrahedron vertex indices.
   */
  constructor(vertices, tetrahedra) {
    this.vertices = vertices;
    this.tetrahedra = tetrahedra;
    this.originalVertexCount = vertices.length;
    this.vertexCount = vertices.length;
    this.tetrahedronCount = tetrahedra.length;

    this.faces = [];
    this.edges = [];
    this.boundaryFaces = [];
    this.boundaryEdges = [];
    this.boundaryNodes = new Set();

    this.faceToTets = [];
    this.edgeToFaces = [];
    /** @type {!Map<number, number>} Maps edge integer key to global edge index. */
    this.edgeKeyToIndex = new Map();

    this.vertexToTets = [];
    this.vertexToBoundaryFaces = [];
    this.vertexToEdges = [];

    /**
     * Orientation signs for edges within each tetrahedron.
     * tetEdgeSigns[tIdx][e] = +1 if local edge order matches global storage, -1 otherwise.
     * @type {!Array<!Array<number>>}
     */
    this.tetEdgeSigns = [];
    /**
     * Orientation signs for faces within each tetrahedron.
     * tetFaceSigns[tIdx][f] = +1 if local face is an even permutation of storage, -1 otherwise.
     * @type {!Array<!Array<number>>}
     */
    this.tetFaceSigns = [];

    this.#buildTopology();
    this.#computeOrientationSigns();

    this.faceBarycenters = [];
    this.tetBarycenters = [];
    this.alfeldTriangles = [];
    this.worseyFarinTetrahedra = [];
  }

  /** @private */
  #buildTopology() {
    const faceMap = new Map();
    const edgeMap = new Map();
    const vc = this.originalVertexCount;

    const edgeKey = (e) => {
      const a = e[0];
      const b = e[1];
      return a < b ? a * vc + b : b * vc + a;
    };

    for (let tIdx = 0; tIdx < this.tetrahedronCount; tIdx++) {
      const tet = this.tetrahedra[tIdx];
      const localFaces = [
        [tet[1], tet[2], tet[3]],
        [tet[0], tet[2], tet[3]],
        [tet[0], tet[1], tet[3]],
        [tet[0], tet[1], tet[2]],
      ];

      for (const f of localFaces) {
        const key = this.#faceKey(f);
        if (!faceMap.has(key)) {
          faceMap.set(key, {verts: f, tets: [tIdx]});
        } else {
          faceMap.get(key).tets.push(tIdx);
        }
      }
    }

    this.faceKeyToIndex = new Map();

    let fIdx = 0;
    for (const [key, data] of faceMap.entries()) {
      this.faces.push(data.verts);
      this.faceToTets.push(data.tets);
      this.faceKeyToIndex.set(key, fIdx);

      const isBoundary = data.tets.length === 1;
      if (isBoundary) {
        this.boundaryFaces.push(fIdx);
        data.verts.forEach((v) => this.boundaryNodes.add(v));
      }

      const fEdges = [
        [data.verts[0], data.verts[1]],
        [data.verts[1], data.verts[2]],
        [data.verts[2], data.verts[0]],
      ];

      for (const e of fEdges) {
        const eKey = edgeKey(e);
        if (!edgeMap.has(eKey)) {
          edgeMap.set(eKey, {verts: e, isBoundary, faces: [fIdx]});
        } else {
          const eData = edgeMap.get(eKey);
          eData.faces.push(fIdx);
          if (isBoundary) {
            eData.isBoundary = true;
          }
        }
      }
      fIdx++;
    }

    let eIdx = 0;
    for (const [key, data] of edgeMap.entries()) {
      this.edges.push(data.verts);
      this.edgeToFaces.push(data.faces);
      this.edgeKeyToIndex.set(key, eIdx);
      if (data.isBoundary) {
        this.boundaryEdges.push(eIdx);
      }
      eIdx++;
    }

    // Build adjacency lists.
    this.vertexToTets = Array.from({length: this.vertexCount}, () => []);
    for (let tIdx = 0; tIdx < this.tetrahedronCount; tIdx++) {
      for (const v of this.tetrahedra[tIdx]) {
        this.vertexToTets[v].push(tIdx);
      }
    }

    this.vertexToBoundaryFaces = Array.from({length: this.vertexCount}, () => []);
    for (const fIdx of this.boundaryFaces) {
      for (const v of this.faces[fIdx]) {
        this.vertexToBoundaryFaces[v].push(fIdx);
      }
    }

    this.vertexToEdges = Array.from({length: this.vertexCount}, () => []);
    for (let eIdx = 0; eIdx < this.edges.length; eIdx++) {
      for (const v of this.edges[eIdx]) {
        this.vertexToEdges[v].push(eIdx);
      }
    }
  }

  /**
   * Computes orientation signs for edges and faces within each tetrahedron.
   * tetEdgeSigns[tIdx][e] = +1 if local edge order matches global storage.
   * tetFaceSigns[tIdx][f] = +1 if local face is an even permutation of storage.
   * @private
   */
  #computeOrientationSigns() {
    const edgeKey = (e) => {
      const a = e[0];
      const b = e[1];
      return a < b ? a * this.originalVertexCount + b : b * this.originalVertexCount + a;
    };

    this.tetEdgeSigns = new Array(this.tetrahedronCount);
    this.tetFaceSigns = new Array(this.tetrahedronCount);

    for (let tIdx = 0; tIdx < this.tetrahedronCount; tIdx++) {
      const tet = this.tetrahedra[tIdx];

      const edges = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];
      this.tetEdgeSigns[tIdx] = edges.map(([i, j]) => {
        const gEdge = [tet[i], tet[j]];
        const key = edgeKey(gEdge);
        const gIdx = this.edgeKeyToIndex.get(key);
        const stored = this.edges[gIdx];
        return stored[0] === gEdge[0] && stored[1] === gEdge[1] ? 1 : -1;
      });

      const faces = [
        [tet[1], tet[2], tet[3]],
        [tet[0], tet[2], tet[3]],
        [tet[0], tet[1], tet[3]],
        [tet[0], tet[1], tet[2]],
      ];
      this.tetFaceSigns[tIdx] = faces.map((f) => {
        const key = this.#faceKey(f);
        const gIdx = this.faceKeyToIndex.get(key);
        const stored = this.faces[gIdx];
        return this.#facePermutationSign(stored, f);
      });
    }
  }

  /**
   * Returns +1 if local is an even permutation of stored, -1 if odd.
   * @param {!Array<number>} stored
   * @param {!Array<number>} local
   * @return {number}
   * @private
   */
  #facePermutationSign(stored, local) {
    const idx = local.map((v) => stored.indexOf(v));
    // Parity of permutation (idx[0], idx[1], idx[2]) of (0,1,2).
    // The sign of the permutation is the sign of the product of differences.
    const inv = (idx[1] - idx[0]) * (idx[2] - idx[0]) * (idx[2] - idx[1]);
    return inv > 0 ? 1 : -1;
  }

  /** Section 6.1.3: Alfeld split of boundary faces. */
  computeAlfeldSplit() {
    this.faceBarycenters = new Array(this.faces.length).fill(-1);
    this.alfeldTriangles = [];

    this.boundaryFaces.forEach((fIdx) => {
      const f = this.faces[fIdx];
      const bary = this.getFaceBarycenter(fIdx);
      const vBaryIdx = this.vertices.length;
      this.vertices.push(bary);
      this.faceBarycenters[fIdx] = vBaryIdx;

      const subTriangles = [
        [f[0], f[1], vBaryIdx],
        [f[1], f[2], vBaryIdx],
        [f[2], f[0], vBaryIdx],
      ];
      this.alfeldTriangles.push({parentFaceIdx: fIdx, triangles: subTriangles});
    });
    this.vertexCount = this.vertices.length;
  }

  /** Section 6.1.4: Worsey-Farin split of tetrahedra. */
  computeWorseyFarinSplit() {
    if (this.faceBarycenters.length === 0) {
      this.computeAlfeldSplit();
    }

    this.tetBarycenters = new Array(this.tetrahedronCount).fill(-1);
    this.worseyFarinTetrahedra = [];

    for (let tIdx = 0; tIdx < this.tetrahedronCount; tIdx++) {
      const tet = this.tetrahedra[tIdx];
      const bary = this.getTetrahedronBarycenter(tIdx);
      const vTetBaryIdx = this.vertices.length;
      this.vertices.push(bary);
      this.tetBarycenters[tIdx] = vTetBaryIdx;

      const tetSubTets = [];
      const tFaces = this.getTetrahedronFaces(tIdx);

      tFaces.forEach((fIdx) => {
        const f = this.faces[fIdx];
        let fvBaryIdx = this.faceBarycenters[fIdx];
        if (fvBaryIdx === -1) {
          const fbary = this.getFaceBarycenter(fIdx);
          fvBaryIdx = this.vertices.length;
          this.vertices.push(fbary);
          this.faceBarycenters[fIdx] = fvBaryIdx;
        }

        const subTris = [
          [f[0], f[1], fvBaryIdx],
          [f[1], f[2], fvBaryIdx],
          [f[2], f[0], fvBaryIdx],
        ];

        subTris.forEach((tri) => {
          tetSubTets.push([...tri, vTetBaryIdx]);
        });
      });
      this.worseyFarinTetrahedra.push({
        parentTetIdx: tIdx,
        tetrahedra: tetSubTets,
      });
    }
    this.vertexCount = this.vertices.length;
  }

  /**
   * @param {number} fIdx
   * @return {!Array<number>}
   */
  getFaceBarycenter(fIdx) {
    const f = this.faces[fIdx];
    const v = f.map((i) => this.vertices[i]);
    return [
      (v[0][0] + v[1][0] + v[2][0]) / 3,
      (v[0][1] + v[1][1] + v[2][1]) / 3,
      (v[0][2] + v[1][2] + v[2][2]) / 3,
    ];
  }

  /**
   * @param {number} tIdx
   * @return {!Array<number>}
   */
  getTetrahedronBarycenter(tIdx) {
    const t = this.tetrahedra[tIdx];
    const v = t.map((i) => this.vertices[i]);
    return [
      (v[0][0] + v[1][0] + v[2][0] + v[3][0]) / 4,
      (v[0][1] + v[1][1] + v[2][1] + v[3][1]) / 4,
      (v[0][2] + v[1][2] + v[2][2] + v[3][2]) / 4,
    ];
  }

  /**
   * @param {number} vIdx
   * @return {!Array<number>}
   */
  getStar(vIdx) {
    return this.vertexToTets[vIdx] || [];
  }

  /**
   * @param {number} vIdx
   * @return {!Array<number>}
   */
  getBoundaryStar(vIdx) {
    return this.vertexToBoundaryFaces[vIdx] || [];
  }

  /**
   * @param {number} vIdx
   * @return {!Array<number>}
   */
  getEdgeStar(vIdx) {
    return this.vertexToEdges[vIdx] || [];
  }

  /**
   * @param {number} tIdx
   * @return {number}
   */
  getVolume(tIdx) {
    const tet = this.tetrahedra[tIdx];
    const v = tet.map((i) => this.vertices[i]);
    const v0 = [v[1][0] - v[0][0], v[1][1] - v[0][1], v[1][2] - v[0][2]];
    const v1 = [v[2][0] - v[0][0], v[2][1] - v[0][1], v[2][2] - v[0][2]];
    const v2 = [v[3][0] - v[0][0], v[3][1] - v[0][1], v[3][2] - v[0][2]];
    return (
      Math.abs(
        v0[0] * (v1[1] * v2[2] - v1[2] * v2[1]) -
          v0[1] * (v1[0] * v2[2] - v1[2] * v2[0]) +
          v0[2] * (v1[0] * v2[1] - v1[1] * v2[0]),
      ) / 6.0
    );
  }

  /**
   * @param {number} tIdx
   * @return {!Array<number>}
   */
  getTetrahedronFaces(tIdx) {
    const tet = this.tetrahedra[tIdx];
    const local = [
      [tet[1], tet[2], tet[3]],
      [tet[0], tet[2], tet[3]],
      [tet[0], tet[1], tet[3]],
      [tet[0], tet[1], tet[2]],
    ];
    return local.map((lf) => {
      const key = this.#faceKey(lf);
      return this.faceKeyToIndex.get(key);
    });
  }

  /**
   * Integer face key using sorted vertex ids.
   * @param {!Array<number>} f
   * @return {number}
   * @private
   */
  #faceKey(f) {
    const [s0, s1, s2] = sort3(f[0], f[1], f[2]);
    const vc = this.originalVertexCount;
    return ((s0 * vc) + s1) * vc + s2;
  }

  /**
   * Computes the outward-pointing unit normal for a boundary face.
   * For interior faces the orientation is arbitrary.
   * @param {number} fIdx
   * @return {!Array<number>}
   */
  getFaceOutwardNormal(fIdx) {
    const f = this.faces[fIdx];
    const v0 = this.vertices[f[0]];
    const v1 = this.vertices[f[1]];
    const v2 = this.vertices[f[2]];
    const e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
    const e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
    const n = [
      e1[1] * e2[2] - e1[2] * e2[1],
      e1[2] * e2[0] - e1[0] * e2[2],
      e1[0] * e2[1] - e1[1] * e2[0],
    ];

    const tets = this.faceToTets[fIdx];
    if (tets.length === 1) {
      const tet = this.tetrahedra[tets[0]];
      const oppV = tet.find((v) => !f.includes(v));
      const oppP = this.vertices[oppV];
      const centroid = [
        (v0[0] + v1[0] + v2[0]) / 3,
        (v0[1] + v1[1] + v2[1]) / 3,
        (v0[2] + v1[2] + v2[2]) / 3,
      ];
      const toOpp = [oppP[0] - centroid[0], oppP[1] - centroid[1], oppP[2] - centroid[2]];
      const dotProduct = n[0] * toOpp[0] + n[1] * toOpp[1] + n[2] * toOpp[2];
      if (dotProduct > 0) {
        n[0] = -n[0];
        n[1] = -n[1];
        n[2] = -n[2];
      }
    }

    const len = Math.sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    if (len < 1e-12) {
      return [0, 0, 1];
    }
    return [n[0] / len, n[1] / len, n[2] / len];
  }
}
