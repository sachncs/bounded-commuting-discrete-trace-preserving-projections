import type {Mesh} from './mesh.js';

/**
 * Whitney forms and barycentric coordinate utilities on a tetrahedral mesh.
 */
export class Whitney {
  constructor(mesh: Mesh);

  /**
   * Computes barycentric coordinates of a point with respect to a tetrahedron.
   * @param tIdx - Tetrahedron index.
   * @param point - Cartesian point [x, y, z].
   * @returns Barycentric coordinates [lambda0, lambda1, lambda2, lambda3].
   */
  getBarycentric(tIdx: number, point: number[]): number[];

  /**
   * Computes gradients of the barycentric coordinate functions.
   * @param tIdx - Tetrahedron index.
   * @returns Four gradient vectors, one per vertex.
   */
  getGradBarycentric(tIdx: number): number[][];

  /**
   * Edge-based Whitney 1-form (Nedelec) basis functions.
   *
   * For edge e = (i, j):
   *   phi_e = lambda_i * grad(lambda_j) - lambda_j * grad(lambda_i)
   *
   * @param tIdx - Tetrahedron index.
   * @param bary - Barycentric coordinates at which to evaluate.
   * @returns Six edge basis vectors in order [0,1], [0,2], [0,3], [1,2], [1,3], [2,3].
   */
  getEdgeBasis(tIdx: number, bary: number[]): number[][];

  /**
   * Face-based Raviart-Thomas (Whitney 2-form) basis functions.
   *
   * @param tIdx - Tetrahedron index.
   * @param bary - Barycentric coordinates at which to evaluate.
   * @returns Four face basis vectors indexed by opposite vertex: 0, 1, 2, 3.
   */
  getFaceBasis(tIdx: number, bary: number[]): number[][];
}
