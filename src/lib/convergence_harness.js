/**
 * Convergence experiment harness for BCDTPP projection operators.
 *
 * Computes discrete error norms and convergence rates on a sequence of meshes.
 */

import { integrateTetrahedron } from './quadrature.js'
import { ProjectionError } from './errors.js'
import { Bcdtpp } from './bcdtpp.js'
import { Whitney } from './whitney.js'

/**
 * Computes the L2 error between an exact function and its projection.
 *
   * err_L2^2 = Σ_T ∫_T |u_exact - u_proj|^2 dx
   *
 * @param {!Mesh} mesh
 * @param {!Bcdtpp} bcdtpp
 * @param {function(!Array<number>): number} exactFn
 * @param {function(number, !Array<number>): number} projFn
 *   Function taking (tIdx, point) and returning the projected value at that point.
 * @return {number}
 */
export function computeL2ErrorScalar (mesh, bcdtpp, exactFn, projFn) {
  let errSq = 0
  for (let tIdx = 0; tIdx < mesh.tetrahedronCount; tIdx++) {
    const verts = mesh.tetrahedra[tIdx].map((i) => mesh.vertices[i])
    const integrand = (pt) => {
      const exact = exactFn(pt)
      const proj = projFn(tIdx, pt)
      return (exact - proj) * (exact - proj)
    }
    errSq += integrateTetrahedron(verts, integrand, bcdtpp.quadratureOrder)
  }
  return Math.sqrt(Math.max(0, errSq))
}

/**
 * Computes the L2 error for a vector-valued projection.
 *
   * err_L2^2 = Σ_T ∫_T |v_exact - v_proj|^2 dx
   *
 * @param {!Mesh} mesh
 * @param {!Bcdtpp} bcdtpp
 * @param {function(!Array<number>): !Array<number>} exactFn
 * @param {function(number, !Array<number>): !Array<number>} projFn
 * @return {number}
 */
export function computeL2ErrorVector (mesh, bcdtpp, exactFn, projFn) {
  let errSq = 0
  for (let tIdx = 0; tIdx < mesh.tetrahedronCount; tIdx++) {
    const verts = mesh.tetrahedra[tIdx].map((i) => mesh.vertices[i])
    const integrand = (pt) => {
      const exact = exactFn(pt)
      const proj = projFn(tIdx, pt)
      const dx = exact[0] - proj[0]
      const dy = exact[1] - proj[1]
      const dz = exact[2] - proj[2]
      return dx * dx + dy * dy + dz * dz
    }
    errSq += integrateTetrahedron(verts, integrand, bcdtpp.quadratureOrder)
  }
  return Math.sqrt(Math.max(0, errSq))
}

/**
 * Computes the H1 semi-norm error (L2 error of the gradient) for scalar projections.
 * Uses numerical differentiation of the exact function for comparison.
 *
   * err_H1^2 = Σ_T ∫_T |grad(u_exact) - grad(u_proj)|^2 dx
   *
 * @param {!Mesh} mesh
 * @param {!Bcdtpp} bcdtpp
 * @param {function(!Array<number>): number} exactFn
 * @param {function(number, !Array<number>): number} projFn
 * @return {number}
 */
export function computeH1SemiError (mesh, bcdtpp, exactFn, projFn) {
  const h = 1e-6
  const gradExact = (pt) => [
    (exactFn([pt[0] + h, pt[1], pt[2]]) - exactFn([pt[0] - h, pt[1], pt[2]])) / (2 * h),
    (exactFn([pt[0], pt[1] + h, pt[2]]) - exactFn([pt[0], pt[1] - h, pt[2]])) / (2 * h),
    (exactFn([pt[0], pt[1], pt[2] + h]) - exactFn([pt[0], pt[1], pt[2] - h])) / (2 * h)
  ]

  let errSq = 0
  for (let tIdx = 0; tIdx < mesh.tetrahedronCount; tIdx++) {
    const verts = mesh.tetrahedra[tIdx].map((i) => mesh.vertices[i])
    const gradProj = (pt) => {
      const valXm = projFn(tIdx, [pt[0] - h, pt[1], pt[2]])
      const valXp = projFn(tIdx, [pt[0] + h, pt[1], pt[2]])
      const valYm = projFn(tIdx, [pt[0], pt[1] - h, pt[2]])
      const valYp = projFn(tIdx, [pt[0], pt[1] + h, pt[2]])
      const valZm = projFn(tIdx, [pt[0], pt[1], pt[2] - h])
      const valZp = projFn(tIdx, [pt[0], pt[1], pt[2] + h])
      return [
        (valXp - valXm) / (2 * h),
        (valYp - valYm) / (2 * h),
        (valZp - valZm) / (2 * h)
      ]
    }

    const integrand = (pt) => {
      const gE = gradExact(pt)
      const gP = gradProj(pt)
      const dx = gE[0] - gP[0]
      const dy = gE[1] - gP[1]
      const dz = gE[2] - gP[2]
      return dx * dx + dy * dy + dz * dz
    }
    errSq += integrateTetrahedron(verts, integrand, bcdtpp.quadratureOrder)
  }
  return Math.sqrt(Math.max(0, errSq))
}

/**
 * Estimates the mesh size h as the cube root of the average tetrahedron volume
 * scaled to unit volume, or more simply the maximum edge length.
 * @param {!Mesh} mesh
 * @return {number}
 */
export function estimateMeshSize (mesh) {
  let maxEdgeLen = 0
  for (let eIdx = 0; eIdx < mesh.edges.length; eIdx++) {
    const e = mesh.edges[eIdx]
    const v0 = mesh.vertices[e[0]]
    const v1 = mesh.vertices[e[1]]
    const dx = v1[0] - v0[0]
    const dy = v1[1] - v0[1]
    const dz = v1[2] - v0[2]
    const len = Math.sqrt(dx * dx + dy * dy + dz * dz)
    if (len > maxEdgeLen) {
      maxEdgeLen = len
    }
  }
  return maxEdgeLen
}

/**
 * Computes the observed convergence rate between two successive error measurements.
 *
   * rate = log(err_1 / err_2) / log(h_1 / h_2)
   *
 * @param {number} err1 - Error on finer mesh.
 * @param {number} err2 - Error on coarser mesh.
 * @param {number} h1 - Mesh size of finer mesh.
 * @param {number} h2 - Mesh size of coarser mesh.
 * @return {number}
 */
export function computeRate (err1, err2, h1, h2) {
  if (err1 === 0 || err2 === 0 || h1 === 0 || h2 === 0 || h2 === h1) {
    return 0
  }
  return Math.log(err1 / err2) / Math.log(h1 / h2)
}

/**
 * Runs a convergence study on a sequence of meshes.
 *
 * @param {!Array<!Mesh>} meshes - Sequence of progressively finer meshes.
 * @param {Object} config
 * @param {function(!Array<number>): number} config.exactScalar - Exact scalar function.
 * @param {function(!Array<number>): !Array<number>=} config.exactVector - Exact vector function.
 * @param {number=} config.l - Form degree (default 0).
 * @param {number=} config.p - Polynomial degree (default 0).
 * @param {number=} config.quadratureOrder - Quadrature order (default 3).
 * @return {!Array<{h: number, l2Err: number, h1Err: number|undefined, rateL2: number|undefined, rateH1: number|undefined}>}
 */
export function runConvergenceStudy (meshes, config) {
  const results = []
  const {
    exactScalar,
    exactVector,
    l = 0,
    p = 0,
    quadratureOrder = 3
  } = config

  for (let i = 0; i < meshes.length; i++) {
    const mesh = meshes[i]
    // Re-create Bcdtpp for each mesh since it holds mesh-specific data.
    const w = new Whitney(mesh)
    const b = new Bcdtpp(mesh, w, { quadratureOrder })
    b.computeBoundaryWeights()
    b.buildPointLocator()

    const h = estimateMeshSize(mesh)
    let l2Err
    let h1Err

    if (l === 0) {
      const projFn = (tIdx, pt) => b.projectHp(exactScalar, pt, tIdx, l, p)
      l2Err = computeL2ErrorScalar(mesh, b, exactScalar, projFn)
      h1Err = computeH1SemiError(mesh, b, exactScalar, projFn)
    } else if (l === 3) {
      const projFn = (tIdx, pt) => b.projectHp(exactScalar, pt, tIdx, l, p)
      // For L2, H1 semi-norm is not meaningful; just compute L2 error.
      l2Err = computeL2ErrorScalar(mesh, b, exactScalar, projFn)
    } else if (l === 1 || l === 2) {
      if (!exactVector) {
        throw new ProjectionError('exactVector is required for l=1 or l=2')
      }
      const projFn = (tIdx, pt) => b.projectHp(exactVector, pt, tIdx, l, p)
      l2Err = computeL2ErrorVector(mesh, b, exactVector, projFn)
    }

    const prev = results.length > 0 ? results[results.length - 1] : null
    const rateL2 = prev ? computeRate(l2Err, prev.l2Err, h, prev.h) : undefined
    const rateH1 = prev && h1Err !== undefined ? computeRate(h1Err, prev.h1Err, h, prev.h) : undefined

    results.push({ h, l2Err, h1Err, rateL2, rateH1 })
  }

  return results
}
