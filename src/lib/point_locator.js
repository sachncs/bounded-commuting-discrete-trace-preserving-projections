/**
 * Axis-aligned bounding box (AABB) tree for fast point-in-tetrahedron queries.
 *
 * Recursively partitions the mesh bounding box to achieve O(log N) search.
 */

export class PointLocator {
  /**
   * @param {!Mesh} mesh
   * @param {number=} maxLeafSize - Maximum tets per leaf node (default 8).
   */
  constructor(mesh, maxLeafSize = 8) {
    this.mesh = mesh;
    this.maxLeafSize = maxLeafSize;
    this.root = this.#buildTree(mesh.tetrahedra.map((_, i) => i));
  }

  /**
   * Builds the AABB tree recursively.
   * @param {!Array<number>} tetIndices
   * @return {!Object}
   */
  #buildTree(tetIndices) {
    const aabb = this.#computeAabb(tetIndices);

    if (tetIndices.length <= this.maxLeafSize) {
      return {aabb, tets: tetIndices, left: null, right: null};
    }

    const dims = [
      aabb.max[0] - aabb.min[0],
      aabb.max[1] - aabb.min[1],
      aabb.max[2] - aabb.min[2],
    ];
    let axis = 0;
    if (dims[1] > dims[axis]) {
      axis = 1;
    }
    if (dims[2] > dims[axis]) {
      axis = 2;
    }
    const mid = (aabb.min[axis] + aabb.max[axis]) * 0.5;

    const leftTets = [];
    const rightTets = [];
    for (const tIdx of tetIndices) {
      const center = this.#getTetrahedronCenter(tIdx);
      if (center[axis] < mid) {
        leftTets.push(tIdx);
      } else {
        rightTets.push(tIdx);
      }
    }

    if (leftTets.length === 0 || rightTets.length === 0) {
      const midIdx = Math.floor(tetIndices.length / 2);
      leftTets.length = 0;
      rightTets.length = 0;
      for (let i = 0; i < tetIndices.length; i++) {
        if (i < midIdx) {
          leftTets.push(tetIndices[i]);
        } else {
          rightTets.push(tetIndices[i]);
        }
      }
    }

    return {
      aabb,
      tets: null,
      left: this.#buildTree(leftTets),
      right: this.#buildTree(rightTets),
    };
  }

  /**
   * @param {!Array<number>} tetIndices
   * @return {{min: !Array<number>, max: !Array<number>}}
   */
  #computeAabb(tetIndices) {
    let minX = Infinity;
    let minY = Infinity;
    let minZ = Infinity;
    let maxX = -Infinity;
    let maxY = -Infinity;
    let maxZ = -Infinity;
    for (const tIdx of tetIndices) {
      const tet = this.mesh.tetrahedra[tIdx];
      for (const vIdx of tet) {
        const v = this.mesh.vertices[vIdx];
        if (v[0] < minX) {
          minX = v[0];
        }
        if (v[1] < minY) {
          minY = v[1];
        }
        if (v[2] < minZ) {
          minZ = v[2];
        }
        if (v[0] > maxX) {
          maxX = v[0];
        }
        if (v[1] > maxY) {
          maxY = v[1];
        }
        if (v[2] > maxZ) {
          maxZ = v[2];
        }
      }
    }
    return {min: [minX, minY, minZ], max: [maxX, maxY, maxZ]};
  }

  /**
   * @param {number} tIdx
   * @return {!Array<number>}
   */
  #getTetrahedronCenter(tIdx) {
    const tet = this.mesh.tetrahedra[tIdx];
    let cx = 0;
    let cy = 0;
    let cz = 0;
    for (const vIdx of tet) {
      const v = this.mesh.vertices[vIdx];
      cx += v[0];
      cy += v[1];
      cz += v[2];
    }
    return [cx / 4, cy / 4, cz / 4];
  }

  /**
   * @param {!Array<number>} point
   * @param {{min: !Array<number>, max: !Array<number>}} aabb
   * @return {boolean}
   */
  #pointInAabb(point, aabb) {
    return (
      point[0] >= aabb.min[0] &&
      point[0] <= aabb.max[0] &&
      point[1] >= aabb.min[1] &&
      point[1] <= aabb.max[1] &&
      point[2] >= aabb.min[2] &&
      point[2] <= aabb.max[2]
    );
  }

  /**
   * Computes barycentric coordinates and checks if the point lies inside.
   * @param {number} tIdx
   * @param {!Array<number>} point
   * @return {?Array<number>}
   */
  #barycentricInTetrahedron(tIdx, point) {
    const tet = this.mesh.tetrahedra[tIdx];
    const v = tet.map((i) => this.mesh.vertices[i]);
    const e0 = [v[0][0] - v[3][0], v[0][1] - v[3][1], v[0][2] - v[3][2]];
    const e1 = [v[1][0] - v[3][0], v[1][1] - v[3][1], v[1][2] - v[3][2]];
    const e2 = [v[2][0] - v[3][0], v[2][1] - v[3][1], v[2][2] - v[3][2]];
    const b = [point[0] - v[3][0], point[1] - v[3][1], point[2] - v[3][2]];

    const det3 = (a, b_, c) =>
      a[0] * (b_[1] * c[2] - b_[2] * c[1]) -
      a[1] * (b_[0] * c[2] - b_[2] * c[0]) +
      a[2] * (b_[0] * c[1] - b_[1] * c[0]);

    const det = det3(e0, e1, e2);
    if (Math.abs(det) < 1e-12) {
      return null;
    }

    const lambda0 = det3(b, e1, e2) / det;
    const lambda1 = det3(e0, b, e2) / det;
    const lambda2 = det3(e0, e1, b) / det;
    const lambda3 = 1 - lambda0 - lambda1 - lambda2;
    const tol = -1e-9;
    if (lambda0 >= tol && lambda1 >= tol && lambda2 >= tol && lambda3 >= tol) {
      return [lambda0, lambda1, lambda2, lambda3];
    }
    return null;
  }

  /**
   * Finds the tetrahedron containing the given point.
   * @param {!Array<number>} point - [x, y, z].
   * @return {{tIdx: number, bary: !Array<number>}|null}
   */
  findTetrahedron(point) {
    return this.#findInNode(this.root, point);
  }

  /**
   * @param {!Object} node
   * @param {!Array<number>} point
   * @return {{tIdx: number, bary: !Array<number>}|null}
   */
  #findInNode(node, point) {
    if (!this.#pointInAabb(point, node.aabb)) {
      return null;
    }

    if (node.tets) {
      for (const tIdx of node.tets) {
        const bary = this.#barycentricInTetrahedron(tIdx, point);
        if (bary) {
          return {tIdx, bary};
        }
      }
      return null;
    }

    const leftResult = this.#findInNode(node.left, point);
    if (leftResult) {
      return leftResult;
    }
    return this.#findInNode(node.right, point);
  }
}
