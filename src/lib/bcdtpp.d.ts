import type {Mesh} from './mesh.js';
import type {Whitney} from './whitney.js';

/**
 * BCDTPP: Bounded, Commuting, Discrete-trace Preserving Projections.
 *
 * Implements the de Rham projection operators Pi^l for l = 0,1,2,3 on
 * tetrahedral meshes with boundary-aware trace preservation.
 */
export class Bcdtpp {
  constructor(
    mesh: Mesh,
    whitney: Whitney,
    options?: {
      quadratureOrder?: number;
      onWarning?: (ctx: {code: string; severity: 'warn' | 'error'; message: string}) => void;
    },
  );

  /**
   * Builds the AABB tree for point location.
   * Must be called before projectAtPoint if the tree was not built automatically.
   */
  buildPointLocator(): void;

  /**
   * Section 6.3.1: Construction of lowest-order vertex weights.
   * Computes boundary-aware weights for trace-preserving projections.
   */
  computeBoundaryWeights(): void;

  /**
   * H1 projection (l=0) of a scalar field.
   * @param u - Scalar function.
   * @param point - Cartesian point [x, y, z].
   * @param tIdx - Tetrahedron index.
   * @returns Projected scalar value.
   */
  projectH1(u: (point: number[]) => number, point: number[], tIdx: number): number;

  /**
   * H(curl) projection (l=1).
   * Boundary edges use exact trace DoFs ∫_e u·t ds.
   * Interior edges use midpoint evaluation.
   * @param u - Scalar or vector function.
   * @param point - Cartesian point [x, y, z].
   * @param tIdx - Tetrahedron index.
   * @returns Projected vector [vx, vy, vz].
   */
  projectHcurl(
    u: (point: number[]) => number | number[],
    point: number[],
    tIdx: number,
  ): number[];

  /**
   * H(div) projection (l=2).
   * Boundary faces use exact trace DoFs ∫_f u·n dA.
   * Interior faces use barycenter evaluation.
   * @param u - Scalar or vector function.
   * @param point - Cartesian point [x, y, z].
   * @param tIdx - Tetrahedron index.
   * @returns Projected vector [vx, vy, vz].
   */
  projectHdiv(
    u: (point: number[]) => number | number[],
    point: number[],
    tIdx: number,
  ): number[];

  /**
   * L2 projection (l=3).
   * @param u - Scalar function.
   * @param tIdx - Tetrahedron index.
   * @returns Projected scalar value (average over tet).
   */
  projectL2(u: (point: number[]) => number, tIdx: number): number;

  /**
   * Higher-order projection.
   * @param u - Scalar or vector function.
   * @param point - Cartesian point [x, y, z].
   * @param tIdx - Tetrahedron index.
   * @param l - Form degree (0=H1 scalar, 1=Hcurl vector, 2=Hdiv vector, 3=L2 scalar).
   * @param p - Polynomial degree (0 returns lowest-order projection).
   * @returns Projected value (number for l=0 or l=3, 3-vector for l=1 or l=2).
   */
  projectHp(
    u: (point: number[]) => number | number[],
    point: number[],
    tIdx: number,
    l: number,
    p: number,
  ): number | number[];

  /**
   * Extracts the boundary degrees of freedom for a given function.
   * @param u - Scalar or vector function.
   * @param l - Form degree (0,1,2,3).
   * @returns Map from boundary entity index to DoF value.
   */
  extractBoundaryDofs(
    u: (point: number[]) => number | number[],
    l: number,
  ): Map<number, number>;

  /**
   * Pi_ring^l: interior projector with zero boundary trace.
   * @param u - Scalar or vector function.
   * @param point - Cartesian point [x, y, z].
   * @param tIdx - Tetrahedron index.
   * @param l - Form degree (0,1,2,3).
   * @returns Projected value (number for l=0 or l=3, 3-vector for l=1 or l=2).
   */
  projectRing(
    u: (point: number[]) => number | number[],
    point: number[],
    tIdx: number,
    l: number,
  ): number | number[];

  /**
   * E^l: discrete extension operator.
   * Takes boundary DoFs and produces a Whitney form whose trace on the
   * boundary equals the given data.
   * @param boundaryData - Map from boundary entity index to DoF.
   * @param point - Cartesian point [x, y, z].
   * @param tIdx - Tetrahedron index.
   * @param l - Form degree (0,1,2).
   * @returns Extended value.
   */
  extendBoundary(
    boundaryData: Map<number, number>,
    point: number[],
    tIdx: number,
    l: number,
  ): number | number[];

  /**
   * Pi_partial^l: boundary correction part of the projection.
   * Computes the discrete extension of the boundary trace of u.
   * @param u - Scalar or vector function.
   * @param point - Cartesian point [x, y, z].
   * @param tIdx - Tetrahedron index.
   * @param l - Form degree (0,1,2,3).
   * @returns Boundary correction value.
   */
  projectPartial(
    u: (point: number[]) => number | number[],
    point: number[],
    tIdx: number,
    l: number,
  ): number | number[];

  /**
   * Finds the tetrahedron containing the point and projects.
   * @param u - Scalar or vector function.
   * @param point - Cartesian point [x, y, z].
   * @param l - Form degree (default 0).
   * @param p - Polynomial degree (default 0).
   * @returns Object with projected value, tetrahedron index, and barycentric coordinates.
   */
  projectAtPoint(
    u: (point: number[]) => number | number[],
    point: number[],
    l?: number,
    p?: number,
  ): {
    value: number | number[];
    tIdx: number;
    bary: number[];
  };
}
