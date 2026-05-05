/**
 * BCDTPP: Bounded, Commuting, Discrete-trace Preserving Projections.
 *
 * Implements the de Rham projection operators Pi^l for l = 0,1,2,3 on
 * tetrahedral meshes with boundary-aware trace preservation.
 */

import {
  subtract, norm, dot, triangleArea,
} from './math_utils.js';
import {BoundaryWeightComputer} from './boundary_weight_computer.js';
import {
  barycentricToCartesian,
  integrateTetrahedron,
  triangleQuadrature,
  lineQuadrature,
} from './quadrature.js';
import {PointLocator} from './point_locator.js';
import {HigherOrderProjection} from './higher_order_projection.js';
import {MeshRefinement} from './mesh_refinement.js';

/**
 * BCDTPP: Bounded, Commuting, Discrete-trace Preserving Projections.
 *
 * Implements the de Rham projection operators Pi^l for l = 0,1,2,3 on
 * tetrahedral meshes with boundary-aware trace preservation.
 *
 * **Coupling note:** This class is tightly coupled to the {@link Mesh} data
 * structure.  It accesses mesh internals (vertices, faces, edges, boundary
 * flags, orientation signs, etc.) directly.  Swapping in a different mesh
 * implementation requires updating every property access inside this class.
 */
export class Bcdtpp {
  /** @type {!Map<number, {nodeMap: !Array<number>, psi: !Array<number>}>} */
  #zeta0Vertex;
  /** @type {!Map<number, {v0: number, v1: number, tangent: !Array<number>, length: number}>} */
  #zeta1Edge;
  /** @type {!Map<number, {normal: !Array<number>, area: number}>} */
  #zeta2Face;

  /**
   * @param {!Mesh} mesh
   * @param {!Whitney} whitney
   * @param {!Object=} options
   * @param {number=} options.quadratureOrder - Quadrature order for integration (default 3).
   */
  constructor(mesh, whitney, options = {}) {
    this.mesh = mesh;
    this.whitney = whitney;
    this.quadratureOrder = options.quadratureOrder || 3;
    this.onWarning = options.onWarning || ((msg) => console.warn(msg));
    this.#zeta0Vertex = new Map();
    this.#zeta1Edge = new Map();
    this.#zeta2Face = new Map();
    this.pointLocator = null;

    this.meshRefinement = new MeshRefinement(this.mesh);
    this.meshRefinement.computeWorseyFarinSplit();
    this.boundaryWeightComputer = new BoundaryWeightComputer(
      this.mesh,
      this.meshRefinement,
      this.onWarning,
    );

    this.higherOrderProjector = new HigherOrderProjection(
      this.mesh,
      this.whitney,
      this.quadratureOrder,
      this.onWarning,
    );

    this.#validateMesh();
  }

  /** @private */
  #validateMesh() {
    let degenerateCount = 0;
    for (let tIdx = 0; tIdx < this.mesh.tetrahedronCount; tIdx++) {
      const vol = this.mesh.getVolume(tIdx);
      if (vol < 1e-12) {
        degenerateCount++;
      }
    }
    if (degenerateCount > 0) {
      this.onWarning(
        `Mesh contains ${degenerateCount} degenerate or near-degenerate ` +
          `tetrahedra. Projections may fail.`,
      );
    }
  }

  buildPointLocator() {
    this.pointLocator = new PointLocator(this.mesh);
  }

  /** @private */
  #validateTetIdx(tIdx) {
    if (typeof tIdx !== 'number' || !Number.isInteger(tIdx)) {
      throw new Error(`tIdx must be an integer, got ${tIdx}`);
    }
    if (tIdx < 0 || tIdx >= this.mesh.tetrahedronCount) {
      throw new Error(
        `tIdx=${tIdx} out of range [0, ${this.mesh.tetrahedronCount - 1}]`,
      );
    }
  }

  /** Section 6.3.1: Construction of lowest-order vertex weights. */
  computeBoundaryWeights() {
    const weights = this.boundaryWeightComputer.compute();
    this.#zeta0Vertex = weights.zeta0Vertex;
    this.#zeta1Edge = weights.zeta1Edge;
    this.#zeta2Face = weights.zeta2Face;
  }

  /**
   * H1 projection (l=0) of a scalar field.
   * @param {function(!Array<number>): number} u
   * @param {!Array<number>} point
   * @param {number} tIdx
   * @return {number}
   */
  projectH1(u, point, tIdx) {
    this.#validateTetIdx(tIdx);
    const tet = this.mesh.tetrahedra[tIdx];
    const bary = this.whitney.getBarycentric(tIdx, point);

    let result = 0;
    for (let i = 0; i < 4; i++) {
      const vIdx = tet[i];
      if (this.mesh.boundaryNodes.has(vIdx)) {
        const alpha = this.#computeBoundaryIntegralH1(vIdx, u);
        result += alpha * bary[i];
      } else {
        result += u(this.mesh.vertices[vIdx]) * bary[i];
      }
    }
    return result;
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
  projectHcurl(u, point, tIdx) {
    this.#validateTetIdx(tIdx);
    const tet = this.mesh.tetrahedra[tIdx];
    const bary = this.whitney.getBarycentric(tIdx, point);
    const edgeBasis = this.whitney.getEdgeBasis(tIdx, bary);

    const localEdges = [
      [0, 1],
      [0, 2],
      [0, 3],
      [1, 2],
      [1, 3],
      [2, 3],
    ];
    let result = [0, 0, 0];
    const vc = this.mesh.originalVertexCount;
    const edgeKey = (e) => {
      const a = e[0];
      const b = e[1];
      return a < b ? a * vc + b : b * vc + a;
    };

    for (let e = 0; e < 6; e++) {
      const [i, j] = localEdges[e];
      const globalEdge = [tet[i], tet[j]];
      const eKey = edgeKey(globalEdge);
      const eIdx = this.mesh.getEdgeIndex(eKey);
      const sigma = this.mesh.getTetEdgeSign(tIdx, e);

      let coefficient;
      if (this.mesh.boundaryEdges.includes(eIdx)) {
        coefficient = this.#computeEdgeDof(u, eIdx);
      } else {
        const mid = this.#midpoint(tet[i], tet[j]);
        const edgeVec = subtract(this.mesh.vertices[tet[j]], this.mesh.vertices[tet[i]]);
        const val = u(mid);
        if (typeof val === 'number') {
          const grad = this.#numericalGradient(u, mid);
          coefficient = dot(grad, edgeVec);
        } else {
          coefficient = dot(val, edgeVec);
        }
      }

      coefficient *= sigma;
      result[0] += coefficient * edgeBasis[e][0];
      result[1] += coefficient * edgeBasis[e][1];
      result[2] += coefficient * edgeBasis[e][2];
    }

    return result;
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
  projectHdiv(u, point, tIdx) {
    this.#validateTetIdx(tIdx);
    const tet = this.mesh.tetrahedra[tIdx];
    const bary = this.whitney.getBarycentric(tIdx, point);
    const faceBasis = this.whitney.getFaceBasis(tIdx, bary);

    const tFaces = this.mesh.getTetrahedronFaces(tIdx);
    let result = [0, 0, 0];

    for (let f = 0; f < 4; f++) {
      const fIdx = tFaces[f];
      const sigma = this.mesh.getTetFaceSign(tIdx, f);

      let coefficient;
      if (this.mesh.boundaryFaces.includes(fIdx)) {
        coefficient = this.#computeFaceDof(u, fIdx);
      } else {
        const faceBary = this.mesh.getFaceBarycenter(fIdx);
        const normal = this.mesh.getFaceOutwardNormal(fIdx);
        const area = this.mesh.getFaceArea(fIdx);
        const val = u(faceBary);
        if (typeof val === 'number') {
          const grad = this.#numericalGradient(u, faceBary);
          coefficient = dot(grad, normal) * area;
        } else {
          coefficient = dot(val, normal) * area;
        }
      }

      coefficient *= sigma;
      result[0] += coefficient * faceBasis[f][0];
      result[1] += coefficient * faceBasis[f][1];
      result[2] += coefficient * faceBasis[f][2];
    }

    return result;
  }

  /**
   * L2 projection (l=3).
   * @param {function(!Array<number>): number} u
   * @param {number} tIdx
   * @return {number}
   */
  projectL2(u, tIdx) {
    this.#validateTetIdx(tIdx);
    const tet = this.mesh.tetrahedra[tIdx];
    const verts = tet.map((i) => this.mesh.vertices[i]);
    const vol = this.mesh.getVolume(tIdx);
    if (vol < 1e-12) {
      throw new Error(
        `Cannot perform L2 projection on degenerate tetrahedron ${tIdx}`,
      );
    }

    const integrand = (pt) => u(pt);
    const integral = integrateTetrahedron(verts, integrand, this.quadratureOrder);
    return integral / vol;
  }

  /**
   * Higher-order projection.
   * @param {function(!Array<number>): (number|!Array<number>)} u
   * @param {!Array<number>} point
   * @param {number} tIdx
   * @param {number} l - Form degree (0=H1 scalar, 1=Hcurl vector, 2=Hdiv vector, 3=L2 scalar).
   * @param {number} p
   * @return {(number|!Array<number>)} Returns a number for l=0 or l=3, and a 3-vector for l=1 or l=2.
   */
  projectHp(u, point, tIdx, l, p) {
    if (p === 0) {
      const dispatch = {
        0: () => this.projectH1(u, point, tIdx),
        1: () => this.projectHcurl(u, point, tIdx),
        2: () => this.projectHdiv(u, point, tIdx),
        3: () => this.projectL2(u, tIdx),
      };
      if (!dispatch.hasOwnProperty(l)) {
        throw new Error(`Invalid form degree l=${l}`);
      }
      return dispatch[l]();
    }

    if (l === 0) {
      const base = this.projectH1(u, point, tIdx);
      if (p < 4) {
        return base;
      }
      const residualFn = (pt) => u(pt) - base;
      const coeffs = this.higherOrderProjector.solveBubbleProjection(
        tIdx,
        p,
        residualFn,
      );
      if (coeffs) {
        const bubble = this.higherOrderProjector.evaluateBubble(
          tIdx,
          p,
          coeffs,
          point,
        );
        return base + bubble;
      }
      return base;
    }

    if (l === 3) {
      const bary = this.whitney.getBarycentric(tIdx, point);
      const coeffs = this.higherOrderProjector.solveL2Projection(
        tIdx,
        p,
        u,
      );
      return this.higherOrderProjector.evaluateL2Projection(
        coeffs,
        bary,
        p,
      );
    }

    throw new Error(
      `Higher-order projections for l=${l}, p=${p} not yet implemented. ` +
        `Only l=0 and l=3 support p>0.`,
    );
  }

  /**
   * Finds the tetrahedron containing the point and projects.
   * @param {function(!Array<number>): (number|!Array<number>)} u
   * @param {!Array<number>} point
   * @param {number=} l
   * @param {number=} p
   * @return {{value: (number|!Array<number>), tIdx: number, bary: !Array<number>}}
   */
  projectAtPoint(u, point, l = 0, p = 0) {
    if (!this.pointLocator) {
      this.buildPointLocator();
    }
    const found = this.pointLocator.findTetrahedron(point);
    if (!found) {
      throw new Error(
        `Point [${point.join(', ')}] not found in any tetrahedron`,
      );
    }
    const {tIdx, bary} = found;
    const physPoint = barycentricToCartesian(
      this.mesh.tetrahedra[tIdx].map((i) => this.mesh.vertices[i]),
      bary,
    );
    return {
      value: this.projectHp(u, physPoint, tIdx, l, p),
      tIdx,
      bary,
    };
  }

  /** @private */
  #computeBoundaryIntegralH1(vIdx, u) {
    const data = this.#zeta0Vertex.get(vIdx);
    if (!data) {
      return u(this.mesh.vertices[vIdx]);
    }

    const {nodeMap, psi} = data;
    const invNodeMap = new Map(nodeMap.map((id, i) => [id, i]));

    let integral = 0;
    let zetaIntegral = 0;
    const starFaces = this.mesh.getBoundaryStar(vIdx);

    starFaces.forEach((fIdx) => {
      const at = this.meshRefinement.faceToAlfeld.get(fIdx);
      if (!at) return;
      at.triangles.forEach((tri) => {
        const verts = tri.map((i) => this.mesh.vertices[i]);
        const area = triangleArea(verts[0], verts[1], verts[2]);
        const localPsi = tri.map((v) => psi[invNodeMap.get(v)]);
        const {bary, weights} = triangleQuadrature(2);
        let triIntegral = 0;
        let triZetaIntegral = 0;
        for (let q = 0; q < bary.length; q++) {
          const pt = barycentricToCartesian(verts, bary[q]);
          const zeta =
            localPsi[0] * bary[q][0] +
            localPsi[1] * bary[q][1] +
            localPsi[2] * bary[q][2];
          triIntegral += weights[q] * u(pt) * zeta;
          triZetaIntegral += weights[q] * zeta;
        }
        integral += triIntegral * area;
        zetaIntegral += triZetaIntegral * area;
      });
    });

    if (Math.abs(zetaIntegral) < 1e-12) {
      return u(this.mesh.vertices[vIdx]);
    }
    return integral / zetaIntegral;
  }

  /**
   * Computes the exact edge DoF for H(curl): ∫_e u·t ds.
   * For scalar u, this reduces to u(v1) - u(v0) (the fundamental theorem).
   * @param {function(!Array<number>): (number|!Array<number>)} u
   * @param {number} eIdx
   * @return {number}
   * @private
   */
  #computeEdgeDof(u, eIdx) {
    const e = this.mesh.edges[eIdx];
    const v0 = this.mesh.vertices[e[0]];
    const v1 = this.mesh.vertices[e[1]];
    const edgeVec = subtract(v1, v0);

    if (typeof u(v0) === 'number') {
      return u(v1) - u(v0);
    }

    const {points, weights} = lineQuadrature(this.quadratureOrder);
    let integral = 0;
    for (let q = 0; q < points.length; q++) {
      const s = points[q];
      const pt = [
        v0[0] + s * edgeVec[0],
        v0[1] + s * edgeVec[1],
        v0[2] + s * edgeVec[2],
      ];
      integral += weights[q] * dot(u(pt), edgeVec);
    }
    return integral;
  }

  /**
   * Computes the exact face DoF for H(div): ∫_f u·n dA.
   * For scalar u, this integrates grad(u)·n over the face.
   * @param {function(!Array<number>): (number|!Array<number>)} u
   * @param {number} fIdx
   * @return {number}
   * @private
   */
  #computeFaceDof(u, fIdx) {
    const f = this.mesh.faces[fIdx];
    const verts = f.map((i) => this.mesh.vertices[i]);
    const normal = this.mesh.getFaceOutwardNormal(fIdx);
    const area = triangleArea(verts[0], verts[1], verts[2]);
    const {bary, weights} = triangleQuadrature(this.quadratureOrder);

    const isScalar = typeof u(verts[0]) === 'number';
    let integral = 0;
    for (let q = 0; q < bary.length; q++) {
      const pt = barycentricToCartesian(verts, bary[q]);
      if (isScalar) {
        const grad = this.#numericalGradient(u, pt);
        integral += weights[q] * dot(grad, normal);
      } else {
        integral += weights[q] * dot(u(pt), normal);
      }
    }
    return integral * area;
  }

  /** @private */
  #numericalGradient(u, pt, h = 1e-6) {
    return [
      (u([pt[0] + h, pt[1], pt[2]]) - u([pt[0] - h, pt[1], pt[2]])) / (2 * h),
      (u([pt[0], pt[1] + h, pt[2]]) - u([pt[0], pt[1] - h, pt[2]])) / (2 * h),
      (u([pt[0], pt[1], pt[2] + h]) - u([pt[0], pt[1], pt[2] - h])) / (2 * h),
    ];
  }

  /** @private */
  #midpoint(vI, vJ) {
    const a = this.mesh.vertices[vI];
    const b = this.mesh.vertices[vJ];
    return [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2];
  }

}
