/**
 * Returns quadrature points and weights on the reference triangle.
 * Weights are normalized so that sum(w_i) = 1.
 * @param order - Target polynomial exactness (1, 2, or 3).
 */
export function triangleQuadrature(order: number): {
  bary: number[][];
  weights: number[];
};

/**
 * Returns quadrature points and weights on the reference tetrahedron.
 * Weights are normalized so that sum(w_i) = 1.
 * @param order - Target polynomial exactness (1, 2, or 3).
 */
export function tetrahedronQuadrature(order: number): {
  bary: number[][];
  weights: number[];
};

/**
 * Integrates a scalar function over a triangle using quadrature.
 * @param vertices - Triangle vertices.
 * @param fn - Scalar function.
 * @param order - Quadrature order (default 2).
 */
export function integrateTriangle(
  vertices: number[][],
  fn: (point: number[]) => number,
  order?: number,
): number;

/**
 * Maps barycentric coordinates to a Cartesian point.
 * @param vertices - Vertices of the element.
 * @param bary - Barycentric coordinates.
 */
export function barycentricToCartesian(
  vertices: number[][],
  bary: number[],
): number[];

/**
 * Returns quadrature points and weights on the reference interval [0, 1].
 * @param order - Target polynomial exactness (1, 2, or 3).
 */
export function lineQuadrature(order: number): {
  points: number[];
  weights: number[];
};

/**
 * Returns a composite quadrature rule by subdividing the reference tetrahedron
 * into 4 sub-tetrahedra via the centroid and applying the base rule on each.
 * @param order - Base quadrature order.
 */
export function compositeTetrahedronQuadrature(order?: number): {
  bary: number[][];
  weights: number[];
};

/**
 * Integrates a scalar function over a tetrahedron using quadrature.
 * @param vertices - Tetrahedron vertices.
 * @param fn - Scalar function.
 * @param order - Quadrature order (default 2).
 */
export function integrateTetrahedron(
  vertices: number[][],
  fn: (point: number[]) => number,
  order?: number,
): number;
