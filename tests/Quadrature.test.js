/**
 * Tests for Gaussian quadrature rules on triangles, tetrahedra, and lines.
 * Verifies exactness at each quadrature order for polynomial integrands
 * of the appropriate degree, and barycentric-to-Cartesian coordinate
 * conversion.
 */
import { expect } from 'chai'
import {
  integrateTriangle,
  integrateTetrahedron,
  barycentricToCartesian,
  lineQuadrature
} from '../src/lib/quadrature.js'

// Verifies exactness of Gaussian quadrature at each supported order.
// Tolerance 1e-10 is used for exact analytical values (machine epsilon
// level for well-conditioned quadrature weights).
describe('Quadrature', () => {
  it('triangle order 1 integrates constants exactly', () => {
    const verts = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
    const result = integrateTriangle(verts, () => 1, 1)
    expect(result).to.be.closeTo(0.5, Math.pow(10, -10))
  })

  it('triangle order 2 integrates linear functions exactly', () => {
    const verts = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
    const result = integrateTriangle(verts, (pt) => pt[0] + pt[1], 2)
    expect(result).to.be.closeTo(1 / 3, Math.pow(10, -10))
  })

  it('triangle order 3 integrates quadratics exactly', () => {
    const verts = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
    const result = integrateTriangle(verts, (pt) => pt[0] * pt[0], 3)
    expect(result).to.be.closeTo(1 / 12, Math.pow(10, -10))
  })

  it('tetrahedron order 1 integrates constants exactly', () => {
    const verts = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
    const result = integrateTetrahedron(verts, () => 1, 1)
    expect(result).to.be.closeTo(1 / 6, Math.pow(10, -10))
  })

  it('tetrahedron order 2 integrates linear functions exactly', () => {
    const verts = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
    const result = integrateTetrahedron(
      verts,
      (pt) => pt[0] + pt[1] + pt[2],
      2
    )
    expect(result).to.be.closeTo(1 / 8, Math.pow(10, -10))
  })

  it('tetrahedron order 3 integrates quadratics exactly', () => {
    const verts = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
    const result = integrateTetrahedron(
      verts,
      (pt) => pt[0] * pt[0] + pt[1] * pt[1] + pt[2] * pt[2],
      3
    )
    expect(result).to.be.closeTo(1 / 20, Math.pow(10, -10))
  })

  it('barycentricToCartesian works for triangle', () => {
    const verts = [[0, 0, 0], [1, 0, 0], [0, 1, 0]]
    const pt = barycentricToCartesian(verts, [0.5, 0.5, 0])
    expect(pt[0]).to.be.closeTo(0.5, Math.pow(10, -10))
    expect(pt[1]).to.be.closeTo(0, Math.pow(10, -10))
    expect(pt[2]).to.be.closeTo(0, Math.pow(10, -10))
  })

  it('lineQuadrature order 1 integrates constants exactly', () => {
    const { points, weights } = lineQuadrature(1)
    let sum = 0
    for (let q = 0; q < points.length; q++) {
      sum += weights[q] * 1
    }
    expect(sum).to.be.closeTo(1, Math.pow(10, -10))
  })

  it('lineQuadrature order 2 integrates linear functions exactly', () => {
    const { points, weights } = lineQuadrature(2)
    let sum = 0
    for (let q = 0; q < points.length; q++) {
      sum += weights[q] * points[q]
    }
    expect(sum).to.be.closeTo(0.5, Math.pow(10, -10))
  })

  it('lineQuadrature order 3 integrates quadratics exactly', () => {
    const { points, weights } = lineQuadrature(3)
    let sum = 0
    for (let q = 0; q < points.length; q++) {
      sum += weights[q] * points[q] * points[q]
    }
    expect(sum).to.be.closeTo(1 / 3, Math.pow(10, -10))
  })
})
