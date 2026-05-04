/**
 * BCDTPP: Bounded, Commuting, Discrete-trace Preserving Projections.
 *
 * Implements the de Rham projection operators Pi^l for l = 0,1,2,3 on
 * tetrahedral meshes with boundary-aware trace preservation.
 */

import {
  subtract, cross, norm, dot, luSolve,
} from './math_utils.js';
import {LocalSolver} from './local_solver.js';
import {
  barycentricToCartesian,
  integrateTetrahedron,
  triangleQuadrature,
  lineQuadrature,
} from './quadrature.js';
import {PointLocator} from './point_locator.js';
import {HigherOrderProjection} from './higher_order_projection.js';

export class Bcdtpp {
  /**
   * @param {!Mesh} mesh
   * @param {!Whitney} whitney
   * @param {!Object=} options
   */
  constructor(mesh, whitney, options = {}) {
    this.mesh = mesh;
    this.whitney = whitney;
    this.quadratureOrder = options.quadratureOrder || 3;
    this.zeta0Vertex = new Map();
    this.zeta1Edge = new Map();
    this.zeta2Face = new Map();
    this.pointLocator = null;

    this.mesh.computeWorseyFarinSplit();
    this.higherOrderProjector = new HigherOrderProjection(
      this.mesh,
      this.whitney,
      this.quadratureOrder,
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
      console.warn(
        `Mesh contains ${degenerateCount} degenerate or near-degenerate ` +
          `tetrahedra. Projections may fail.`,
      );
    }
  }

  buildPointLocator() {
    this.pointLocator = new PointLocator(this.mesh);
  }

  /** Section 6.3.1: Construction of lowest-order vertex weights. */
  computeBoundaryWeights() {
    for (const vIdx of this.mesh.boundaryNodes) {
      const starFaces = this.mesh.getBoundaryStar(vIdx);
      const alfeldTris = this.mesh.alfeldTriangles.filter((at) =>
        starFaces.includes(at.parentFaceIdx),
      );

      const triangles = alfeldTris.flatMap((at) => at.triangles);
      const starNodes = new Set(triangles.flat());
      const nodeMap = Array.from(starNodes);
      const invNodeMap = new Map(nodeMap.map((id, i) => [id, i]));

      const localTris = triangles.map((t) => t.map((v) => invNodeMap.get(v)));
      const localVerts = nodeMap.map((v) => this.mesh.vertices[v]);

      const K = LocalSolver.assembleSurfaceStiffness(localVerts, localTris);
      const b = new Array(nodeMap.length).fill(0);

      const starArea = starFaces.reduce(
        (acc, fIdx) => acc + this.#getFaceArea(fIdx),
        0,
      );
      const eta = 1.0 / starArea;

      localTris.forEach((tri) => {
        const area = this.#getTriangleArea(
          localVerts[tri[0]],
          localVerts[tri[1]],
          localVerts[tri[2]],
        );
        tri.forEach((nodeIdx) => {
          b[nodeIdx] += eta * (area / 3.0);
        });
      });

      const psi = LocalSolver.solveWithConstraint(K, b);
      this.zeta0Vertex.set(vIdx, {nodeMap, psi});
    }

    this.#computeEdgeBoundaryWeights();
    this.#computeFaceBoundaryWeights();
  }

  /** @private */
  #computeEdgeBoundaryWeights() {
    for (const eIdx of this.mesh.boundaryEdges) {
      const e = this.mesh.edges[eIdx];
      const edgeVec = subtract(this.mesh.vertices[e[1]], this.mesh.vertices[e[0]]);
      const edgeLen = norm(edgeVec);
      if (edgeLen < 1e-12) {
        continue;
      }
      this.zeta1Edge.set(eIdx, {
        v0: e[0],
        v1: e[1],
        tangent: edgeVec.map((x) => x / edgeLen),
        length: edgeLen,
      });
    }
  }

  /** @private */
  #computeFaceBoundaryWeights() {
    for (const fIdx of this.mesh.boundaryFaces) {
      const normal = this.mesh.getFaceOutwardNormal(fIdx);
      const area = this.#getFaceArea(fIdx);
      this.zeta2Face.set(fIdx, {normal, area});
    }
  }

  /**
   * H1 projection (l=0) of a scalar field.
   * @param {function(!Array<number>): number} u
   * @param {!Array<number>} point
   * @param {number} tIdx
   * @return {number}
   */
  projectH1(u, point, tIdx) {
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
    const tet = this.mesh.tetrahedra[tIdx];
    const bary = this.whitney.getBarycentric(tIdx, point);
    const edgeBasis = this.whitney.getEdgeBasis(tIdx, bary);

    const testVal = u(this.mesh.getTetrahedronBarycenter(tIdx));
    const isScalar = typeof testVal === 'number';

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
      const eIdx = this.mesh.edgeKeyToIndex.get(eKey);
      const sigma = this.mesh.tetEdgeSigns[tIdx][e];

      let coefficient;
      if (this.mesh.boundaryEdges.includes(eIdx)) {
        coefficient = this.#computeEdgeDof(u, eIdx, isScalar);
      } else {
        const mid = this.#midpoint(tet[i], tet[j]);
        const edgeVec = subtract(this.mesh.vertices[tet[j]], this.mesh.vertices[tet[i]]);
        if (isScalar) {
          const grad = this.#numericalGradient(u, mid);
          coefficient = dot(grad, edgeVec);
        } else {
          coefficient = dot(u(mid), edgeVec);
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
    const tet = this.mesh.tetrahedra[tIdx];
    const bary = this.whitney.getBarycentric(tIdx, point);
    const faceBasis = this.whitney.getFaceBasis(tIdx, bary);

    const testVal = u(this.mesh.getTetrahedronBarycenter(tIdx));
    const isScalar = typeof testVal === 'number';

    const tFaces = this.mesh.getTetrahedronFaces(tIdx);
    let result = [0, 0, 0];

    for (let f = 0; f < 4; f++) {
      const fIdx = tFaces[f];
      const sigma = this.mesh.tetFaceSigns[tIdx][f];

      let coefficient;
      if (this.mesh.boundaryFaces.includes(fIdx)) {
        coefficient = this.#computeFaceDof(u, fIdx, isScalar);
      } else {
        const fVerts = this.mesh.faces[fIdx];
        const faceBary = this.#getTriangleBarycenter(fVerts);
        const normal = this.#getFaceNormal(fIdx);
        const area = this.#getFaceArea(fIdx);
        if (isScalar) {
          const grad = this.#numericalGradient(u, faceBary);
          coefficient = dot(grad, normal) * area;
        } else {
          coefficient = dot(u(faceBary), normal) * area;
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
   * @param {number} l
   * @param {number} p
   * @return {(number|!Array<number>)}
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

    const base = this.projectHp(u, point, tIdx, l, 0);
    console.warn(
      `Higher-order projections for l=${l}, p=${p} not yet fully ` +
        `implemented. Returning lowest-order part with exact boundary DoFs.`,
    );
    return base;
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
    const data = this.zeta0Vertex.get(vIdx);
    if (!data) {
      return u(this.mesh.vertices[vIdx]);
    }

    const {nodeMap, psi} = data;
    const invNodeMap = new Map(nodeMap.map((id, i) => [id, i]));

    let integral = 0;
    let zetaIntegral = 0;
    const starFaces = this.mesh.getBoundaryStar(vIdx);
    const alfeldTris = this.mesh.alfeldTriangles.filter((at) =>
      starFaces.includes(at.parentFaceIdx),
    );

    alfeldTris.forEach((at) => {
      at.triangles.forEach((tri) => {
        const verts = tri.map((i) => this.mesh.vertices[i]);
        const area = this.#getTriangleArea(verts[0], verts[1], verts[2]);
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
   * @param {boolean} isScalar
   * @return {number}
   * @private
   */
  #computeEdgeDof(u, eIdx, isScalar) {
    const e = this.mesh.edges[eIdx];
    const v0 = this.mesh.vertices[e[0]];
    const v1 = this.mesh.vertices[e[1]];
    const edgeVec = subtract(v1, v0);

    if (isScalar) {
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
   * @param {boolean} isScalar
   * @return {number}
   * @private
   */
  #computeFaceDof(u, fIdx, isScalar) {
    const f = this.mesh.faces[fIdx];
    const verts = f.map((i) => this.mesh.vertices[i]);
    const normal = this.mesh.getFaceOutwardNormal(fIdx);
    const area = this.#getTriangleArea(verts[0], verts[1], verts[2]);
    const {bary, weights} = triangleQuadrature(this.quadratureOrder);

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

  /** @private */
  #getFaceArea(fIdx) {
    const f = this.mesh.faces[fIdx];
    return this.#getTriangleArea(
      this.mesh.vertices[f[0]],
      this.mesh.vertices[f[1]],
      this.mesh.vertices[f[2]],
    );
  }

  /** @private */
  #getFaceNormal(fIdx) {
    const f = this.mesh.faces[fIdx];
    const v0 = this.mesh.vertices[f[0]];
    const v1 = this.mesh.vertices[f[1]];
    const v2 = this.mesh.vertices[f[2]];
    const e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
    const e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
    const n = [
      e1[1] * e2[2] - e1[2] * e2[1],
      e1[2] * e2[0] - e1[0] * e2[2],
      e1[0] * e2[1] - e1[1] * e2[0],
    ];
    const len = Math.sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    if (len < 1e-12) {
      return [0, 0, 1];
    }
    return n.map((x) => x / len);
  }

  /** @private */
  #getTriangleArea(p1, p2, p3) {
    const v1 = subtract(p2, p1);
    const v2 = subtract(p3, p1);
    return 0.5 * norm(cross(v1, v2));
  }

  /** @private */
  #getTriangleBarycenter(tri) {
    const v = tri.map((i) => this.mesh.vertices[i]);
    return [
      (v[0][0] + v[1][0] + v[2][0]) / 3,
      (v[0][1] + v[1][1] + v[2][1]) / 3,
      (v[0][2] + v[1][2] + v[2][2]) / 3,
    ];
  }

  /** @private */
  #getStarArea(vIdx) {
    return this.mesh
      .getBoundaryStar(vIdx)
      .reduce((acc, f) => acc + this.#getFaceArea(f), 0);
  }

  /**
   * Backward-compatible alias for projectH1.
   * @param {function(!Array<number>): number} u
   * @param {!Array<number>} point
   * @param {number} tIdx
   * @return {number}
   */
  project(u, point, tIdx) {
    return this.projectH1(u, point, tIdx);
  }
}
