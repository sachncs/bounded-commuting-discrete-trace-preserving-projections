import type {Mesh} from './mesh.js';
import type {Whitney} from './whitney.js';

/**
 * Higher-Order Projection Framework (Section 7).
 *
 * For p >= 1, the projection Pi^l_p is built recursively from the lowest-order
 * projection Pi^l_0 by adding bubble corrections on Alfeld-split patches.
 */
export class HigherOrderProjection {
  /**
   * @param mesh - The mesh.
   * @param whitney - Whitney form utilities.
   * @param quadratureOrder - Quadrature order (default 3).
   * @param onWarning - Optional warning callback.
   */
  constructor(
    mesh: Mesh,
    whitney: Whitney,
    quadratureOrder?: number,
    onWarning?: (ctx: {code: string; severity: 'warn' | 'error'; message: string}) => void,
  );

  /**
   * Evaluates the scalar bubble basis at a point given in barycentric coords.
   * @param bary - Barycentric coordinates [lambda0, lambda1, lambda2, lambda3].
   * @param p - Polynomial degree.
   * @returns Bubble basis values.
   */
  evaluateBubbleBasis(bary: number[], p: number): number[];

  /**
   * Assembles the local mass matrix for the bubble space.
   * @param tIdx - Tetrahedron index.
   * @param p - Polynomial degree.
   * @returns Mass matrix.
   */
  assembleBubbleMass(tIdx: number, p: number): number[][];

  /**
   * Assembles the RHS for bubble projection.
   * @param tIdx - Tetrahedron index.
   * @param p - Polynomial degree.
   * @param residualFn - Residual function.
   * @returns RHS vector.
   */
  assembleBubbleRhs(
    tIdx: number,
    p: number,
    residualFn: (point: number[]) => number,
  ): number[];

  /**
   * Solves for bubble coefficients.
   * @param tIdx - Tetrahedron index.
   * @param p - Polynomial degree.
   * @param residualFn - Residual function.
   * @returns Bubble coefficients, or null if p < 4 or solve fails.
   */
  solveBubbleProjection(
    tIdx: number,
    p: number,
    residualFn: (point: number[]) => number,
  ): number[] | null;

  /**
   * Evaluates the bubble correction at a point.
   * @param tIdx - Tetrahedron index.
   * @param p - Polynomial degree.
   * @param coeffs - Bubble coefficients.
   * @param point - Cartesian point [x, y, z].
   * @returns Bubble correction value.
   */
  evaluateBubble(
    tIdx: number,
    p: number,
    coeffs: number[],
    point: number[],
  ): number;

  /**
   * Evaluates the P^p Bernstein basis at barycentric coordinates.
   * @param bary - Barycentric coordinates.
   * @param p - Polynomial degree.
   * @returns Bernstein basis values.
   */
  evaluateBernsteinBasis(bary: number[], p: number): number[];

  /**
   * Solves the L2 projection of u onto P^p(T) using the Bernstein basis.
   * @param tIdx - Tetrahedron index.
   * @param p - Polynomial degree.
   * @param u - Scalar function to project.
   * @returns Coefficients of the Bernstein basis.
   */
  solveL2Projection(
    tIdx: number,
    p: number,
    u: (point: number[]) => number,
  ): number[];

  /**
   * Evaluates the L2 projection at a point.
   * @param coeffs - Bernstein coefficients.
   * @param bary - Barycentric coordinates.
   * @param p - Polynomial degree.
   * @returns Projected value.
   */
  evaluateL2Projection(coeffs: number[], bary: number[], p: number): number;
}
