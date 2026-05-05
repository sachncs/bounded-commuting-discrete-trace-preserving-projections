/**
 * Gaussian quadrature rules for triangles and tetrahedra.
 *
 * Weights are normalized so that sum(w_i) = 1.  The caller multiplies by the
 * element area or volume:
 *   integral ≈ area * Σ w_i * f(x_i)
 */

/**
 * Returns quadrature points and weights on the reference triangle.
 * @param {number} order - Target polynomial exactness (1, 2, or 3).
 * @return {{bary: !Array<!Array<number>>, weights: !Array<number>}}
 */
export function triangleQuadrature(order) {
  switch (order) {
    case 1:
      return {
        bary: [[1 / 3, 1 / 3, 1 / 3]],
        weights: [1.0],
      };
    case 2:
      return {
        bary: [
          [2 / 3, 1 / 6, 1 / 6],
          [1 / 6, 2 / 3, 1 / 6],
          [1 / 6, 1 / 6, 2 / 3],
        ],
        weights: [1 / 3, 1 / 3, 1 / 3],
      };
    case 3:
    default:
      return {
        bary: [
          [1 / 3, 1 / 3, 1 / 3],
          [3 / 5, 1 / 5, 1 / 5],
          [1 / 5, 3 / 5, 1 / 5],
          [1 / 5, 1 / 5, 3 / 5],
        ],
        weights: [-9 / 16, 25 / 48, 25 / 48, 25 / 48],
      };
  }
}

/**
 * Returns quadrature points and weights on the reference tetrahedron.
 * @param {number} order - Target polynomial exactness (1, 2, or 3).
 * @return {{bary: !Array<!Array<number>>, weights: !Array<number>}}
 */
export function tetrahedronQuadrature(order) {
  switch (order) {
    case 1:
      return {
        bary: [[1 / 4, 1 / 4, 1 / 4, 1 / 4]],
        weights: [1.0],
      };
    case 2: {
      const a = 0.5854101966249685;
      const b = 0.1381966011250105;
      return {
        bary: [
          [a, b, b, b],
          [b, a, b, b],
          [b, b, a, b],
          [b, b, b, a],
        ],
        weights: [1 / 4, 1 / 4, 1 / 4, 1 / 4],
      };
    }
    case 3:
    default:
      return {
        bary: [
          [1 / 4, 1 / 4, 1 / 4, 1 / 4],
          [1 / 3, 1 / 3, 1 / 3, 0],
          [1 / 3, 1 / 3, 0, 1 / 3],
          [1 / 3, 0, 1 / 3, 1 / 3],
          [0, 1 / 3, 1 / 3, 1 / 3],
        ],
        weights: [-4 / 5, 9 / 20, 9 / 20, 9 / 20, 9 / 20],
      };
  }
}

/**
 * Integrates a scalar function over a triangle using quadrature.
 * @param {!Array<!Array<number>>} vertices - Triangle vertices.
 * @param {function(!Array<number>): number} fn - Scalar function.
 * @param {number=} order - Quadrature order (default 2).
 * @return {number}
 */
export function integrateTriangle(vertices, fn, order = 2) {
  const {bary, weights} = triangleQuadrature(order);
  const v0 = [
    vertices[1][0] - vertices[0][0],
    vertices[1][1] - vertices[0][1],
    vertices[1][2] - vertices[0][2],
  ];
  const v1 = [
    vertices[2][0] - vertices[0][0],
    vertices[2][1] - vertices[0][1],
    vertices[2][2] - vertices[0][2],
  ];
  const cross = [
    v0[1] * v1[2] - v0[2] * v1[1],
    v0[2] * v1[0] - v0[0] * v1[2],
    v0[0] * v1[1] - v0[1] * v1[0],
  ];
  const area =
    0.5 *
    Math.sqrt(
      cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2],
    );
  let sum = 0;
  for (let i = 0; i < bary.length; i++) {
    const point = barycentricToCartesian(vertices, bary[i]);
    sum += weights[i] * fn(point);
  }
  return area * sum;
}

/**
 * Maps barycentric coordinates to a Cartesian point.
 * @param {!Array<!Array<number>>} vertices - Vertices of the element.
 * @param {!Array<number>} bary - Barycentric coordinates.
 * @return {!Array<number>}
 */
export function barycentricToCartesian(vertices, bary) {
  const dim = vertices[0].length;
  const point = new Array(dim).fill(0);
  for (let i = 0; i < vertices.length; i++) {
    for (let d = 0; d < dim; d++) {
      point[d] += bary[i] * vertices[i][d];
    }
  }
  return point;
}

/**
 * Returns quadrature points and weights on the reference interval [0, 1].
 * @param {number} order - Target polynomial exactness (1, 2, or 3).
 * @return {{points: !Array<number>, weights: !Array<number>}}
 * @example
 * const {points, weights} = lineQuadrature(2);
 * const integral = points.reduce((s, x, i) => s + weights[i] * f(x), 0);
 */
export function lineQuadrature(order) {
  switch (order) {
    case 1:
      return {
        points: [0.5],
        weights: [1.0],
      };
    case 2:
      return {
        points: [
          0.5 * (1.0 - 1.0 / Math.sqrt(3.0)),
          0.5 * (1.0 + 1.0 / Math.sqrt(3.0)),
        ],
        weights: [0.5, 0.5],
      };
    case 3:
    default:
      return {
        points: [
          0.5,
          0.5 * (1.0 - Math.sqrt(0.6)),
          0.5 * (1.0 + Math.sqrt(0.6)),
        ],
        weights: [4.0 / 9.0, 5.0 / 18.0, 5.0 / 18.0],
      };
  }
}

/**
 * Integrates a scalar function over a tetrahedron using quadrature.
 * @param {!Array<!Array<number>>} vertices - Tetrahedron vertices.
 * @param {function(!Array<number>): number} fn - Scalar function.
 * @param {number=} order - Quadrature order (default 2).
 * @return {number}
 */
export function integrateTetrahedron(vertices, fn, order = 2) {
  const {bary, weights} = tetrahedronQuadrature(order);
  const v0 = [
    vertices[1][0] - vertices[0][0],
    vertices[1][1] - vertices[0][1],
    vertices[1][2] - vertices[0][2],
  ];
  const v1 = [
    vertices[2][0] - vertices[0][0],
    vertices[2][1] - vertices[0][1],
    vertices[2][2] - vertices[0][2],
  ];
  const v2 = [
    vertices[3][0] - vertices[0][0],
    vertices[3][1] - vertices[0][1],
    vertices[3][2] - vertices[0][2],
  ];
  const volume =
    Math.abs(
      v0[0] * (v1[1] * v2[2] - v1[2] * v2[1]) -
        v0[1] * (v1[0] * v2[2] - v1[2] * v2[0]) +
        v0[2] * (v1[0] * v2[1] - v1[1] * v2[0]),
    ) / 6.0;
  let sum = 0;
  for (let i = 0; i < bary.length; i++) {
    const point = barycentricToCartesian(vertices, bary[i]);
    sum += weights[i] * fn(point);
  }
  return volume * sum;
}
