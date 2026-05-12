/**
 * Local patch stiffness assembly and constrained solves for boundary
 * weight computation.
 */
export class LocalSolver {
  /**
   * Assembles the surface stiffness matrix for -Delta_Gamma.
   * @param vertices - Patch vertex coordinates.
   * @param triangles - Patch triangles as local vertex indices.
   * @returns Dense stiffness matrix.
   */
  static assembleSurfaceStiffness(
    vertices: number[][],
    triangles: number[][],
  ): number[][];

  /**
   * Solves K x = b with a mean-zero constraint sum(x) = 0.
   * @param K - Stiffness matrix.
   * @param b - Right-hand side vector.
   * @param onWarning - Optional warning callback.
   * @returns Solution vector.
   */
  static solveWithConstraint(
    K: number[][],
    b: number[],
    onWarning?: (ctx: {code: string; severity: 'warn' | 'error'; message: string}) => void,
  ): number[];
}
