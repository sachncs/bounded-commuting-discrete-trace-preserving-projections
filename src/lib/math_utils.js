/**
 * Pure JavaScript linear algebra and vector utilities.
 *
 * This module replaces the mathjs dependency for the core BCDTPP library,
 * ensuring a pure JavaScript implementation with no external runtime math
 * dependencies.
 */

/** @const {number} */
const EPSILON = 1e-12;

/**
 * Computes the dot product of two 3D vectors.
 * @param {!Array<number>} a
 * @param {!Array<number>} b
 * @return {number}
 */
export function dot(a, b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/**
 * Computes the cross product of two 3D vectors.
 * @param {!Array<number>} a
 * @param {!Array<number>} b
 * @return {!Array<number>}
 */
export function cross(a, b) {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ];
}

/**
 * Subtracts vector b from vector a.
 * @param {!Array<number>} a
 * @param {!Array<number>} b
 * @return {!Array<number>}
 */
export function subtract(a, b) {
  return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
}

/**
 * Computes the Euclidean norm of a vector.
 * @param {!Array<number>} v
 * @return {number}
 */
export function norm(v) {
  return Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

/**
 * Computes the area of a triangle given its three vertices.
 * @param {!Array<number>} p1
 * @param {!Array<number>} p2
 * @param {!Array<number>} p3
 * @return {number}
 */
export function triangleArea(p1, p2, p3) {
  const v1 = subtract(p2, p1);
  const v2 = subtract(p3, p1);
  return 0.5 * norm(cross(v1, v2));
}

/**
 * Creates a zero-initialized dense matrix.
 * @param {number} rows
 * @param {number} cols
 * @return {!Array<!Array<number>>}
 */
export function zeros(rows, cols) {
  return Array.from({length: rows}, () => new Array(cols).fill(0));
}

/**
 * Computes the infinity norm (maximum absolute row sum) of a matrix.
 * @param {!Array<!Array<number>>} a
 * @return {number}
 */
export function infinityNorm(a) {
  const n = a.length;
  let maxRowSum = 0;
  for (let i = 0; i < n; i++) {
    let sum = 0;
    for (let j = 0; j < n; j++) {
      sum += Math.abs(a[i][j]);
    }
    if (sum > maxRowSum) {
      maxRowSum = sum;
    }
  }
  return maxRowSum;
}

/**
 * Solves Ax = b using Gaussian elimination with partial pivoting and row
 * equilibration (pivot scaling).
 *
 * The matrix is scaled so that each row has max-norm 1 before elimination,
 * reducing the risk of false singularity claims on poorly scaled systems.
 *
 * @param {!Array<!Array<number>>} a - Dense square matrix (not modified).
 * @param {!Array<number>} b
 * @return {!Array<number>}
 */
export function luSolve(a, b) {
  const n = a.length;

  // Row equilibration: scale each row so its max element is 1.
  const scale = new Array(n);
  for (let i = 0; i < n; i++) {
    let rowMax = 0;
    for (let j = 0; j < n; j++) {
      const val = Math.abs(a[i][j]);
      if (val > rowMax) rowMax = val;
    }
    scale[i] = rowMax === 0 ? 1 : rowMax;
  }

  const aug = a.map((row, i) => [...row.map((v) => v / scale[i]), b[i] / scale[i]]);

  // Tolerance relative to matrix size and typical entry magnitude.
  const normEst = infinityNorm(a);
  const tol = EPSILON * Math.max(1, normEst) * n;

  for (let col = 0; col < n; col++) {
    let maxRow = col;
    let maxVal = Math.abs(aug[col][col]);
    for (let row = col + 1; row < n; row++) {
      const val = Math.abs(aug[row][col]);
      if (val > maxVal) {
        maxVal = val;
        maxRow = row;
      }
    }
    if (maxVal < tol) {
      throw new Error('Singular or near-singular matrix');
    }
    if (maxRow !== col) {
      [aug[col], aug[maxRow]] = [aug[maxRow], aug[col]];
    }

    for (let row = col + 1; row < n; row++) {
      const factor = aug[row][col] / aug[col][col];
      for (let j = col; j <= n; j++) {
        aug[row][j] -= factor * aug[col][j];
      }
    }
  }

  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let sum = aug[i][n];
    for (let j = i + 1; j < n; j++) {
      sum -= aug[i][j] * x[j];
    }
    x[i] = sum / aug[i][i];
  }
  return x;
}

/**
 * Computes the inverse of a 3x3 matrix.
 * @param {!Array<!Array<number>>} m
 * @return {!Array<!Array<number>>}
 */
export function inverse3x3(m) {
  const det =
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
    m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
    m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
  if (Math.abs(det) < EPSILON) {
    throw new Error('Singular 3x3 matrix');
  }
  const invDet = 1 / det;
  return [
    [
      (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * invDet,
      (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invDet,
      (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invDet,
    ],
    [
      (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invDet,
      (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invDet,
      (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * invDet,
    ],
    [
      (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * invDet,
      (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * invDet,
      (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * invDet,
    ],
  ];
}

/**
 * Solves a 3x3 linear system Ax = b using Cramer's rule.
 * @param {!Array<!Array<number>>} a - 3x3 matrix.
 * @param {!Array<number>} b - 3-vector.
 * @return {!Array<number>}
 */
export function solve3x3(a, b) {
  const det =
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
    a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
    a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
  if (Math.abs(det) < EPSILON) {
    throw new Error('Singular 3x3 matrix');
  }

  const x0 =
    (b[0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) -
      a[0][1] * (b[1] * a[2][2] - a[1][2] * b[2]) +
      a[0][2] * (b[1] * a[2][1] - a[1][1] * b[2])) /
    det;
  const x1 =
    (a[0][0] * (b[1] * a[2][2] - a[1][2] * b[2]) -
      b[0] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) +
      a[0][2] * (a[1][0] * b[2] - b[1] * a[2][0])) /
    det;
  const x2 =
    (a[0][0] * (a[1][1] * b[2] - b[1] * a[2][1]) -
      a[0][1] * (a[1][0] * b[2] - b[1] * a[2][0]) +
      b[0] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])) /
    det;

  return [x0, x1, x2];
}
