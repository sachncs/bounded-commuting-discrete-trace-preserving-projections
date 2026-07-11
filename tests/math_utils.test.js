/**
 * Tests for pure JavaScript linear algebra and vector utilities:
 * dot/cross products, norms, matrix operations, LU solver with partial
 * pivoting, 3x3 inverse and solve, tetrahedron determinant/volume,
 * factorial bounds, and numerical gradient approximation.
 */
import { expect } from 'chai'
import {
  dot,
  cross,
  subtract,
  norm,
  subtractInto,
  crossInto,
  triangleArea,
  zeros,
  infinityNorm,
  luSolve,
  inverse3x3,
  solve3x3,
  tetDeterminant,
  tetVolume,
  factorial,
  numericalGradient
} from '../src/lib/math_utils.js'
import { SingularMatrixError } from '../src/lib/errors.js'

describe('Math Utils', () => {
  it('dot product is correct', () => {
    expect(dot([1, 2, 3], [4, 5, 6])).to.equal(32)
  })

  it('cross product is correct', () => {
    expect(cross([1, 0, 0], [0, 1, 0])).to.deep.equal([0, 0, 1])
  })

  it('subtract works', () => {
    expect(subtract([3, 3, 3], [1, 2, 3])).to.deep.equal([2, 1, 0])
  })

  it('norm works', () => {
    expect(norm([3, 4, 0])).to.equal(5)
  })

  it('subtractInto modifies output array', () => {
    const out = [0, 0, 0]
    subtractInto([3, 3, 3], [1, 2, 3], out)
    expect(out).to.deep.equal([2, 1, 0])
  })

  it('crossInto modifies output array', () => {
    const out = [0, 0, 0]
    crossInto([1, 0, 0], [0, 1, 0], out)
    expect(out).to.deep.equal([0, 0, 1])
  })

  it('triangleArea works', () => {
    const area = triangleArea([0, 0, 0], [1, 0, 0], [0, 1, 0])
    expect(area).to.be.closeTo(0.5, Math.pow(10, -10))
  })

  it('zeros creates correct matrix', () => {
    const z = zeros(2, 3)
    expect(z.length).to.equal(2)
    expect(z[0].length).to.equal(3)
    expect(z[0][0]).to.equal(0)
  })

  it('infinityNorm works', () => {
    const m = [
      [1, 2],
      [3, -4]
    ]
    expect(infinityNorm(m)).to.equal(7)
  })

  it('luSolve solves a simple system', () => {
    const a = [
      [2, 1],
      [1, 3]
    ]
    const b = [5, 8]
    const x = luSolve(a, b)
    expect(x[0]).to.be.closeTo(1.4, Math.pow(10, -6))
    expect(x[1]).to.be.closeTo(2.2, Math.pow(10, -6))
  })

  it('luSolve throws SingularMatrixError for singular matrix', () => {
    const a = [
      [1, 2],
      [2, 4]
    ]
    const b = [3, 6]
    expect(() => luSolve(a, b)).to.throw(SingularMatrixError)
  })

  it('inverse3x3 works', () => {
    const m = [
      [1, 0, 0],
      [0, 2, 0],
      [0, 0, 3]
    ]
    const inv = inverse3x3(m)
    expect(inv[0][0]).to.be.closeTo(1, Math.pow(10, -10))
    expect(inv[1][1]).to.be.closeTo(0.5, Math.pow(10, -10))
    expect(inv[2][2]).to.be.closeTo(1 / 3, Math.pow(10, -10))
  })

  it('inverse3x3 throws SingularMatrixError for singular matrix', () => {
    const m = [
      [1, 2, 3],
      [2, 4, 6],
      [3, 6, 9]
    ]
    expect(() => inverse3x3(m)).to.throw(SingularMatrixError)
  })

  it('solve3x3 works', () => {
    const a = [
      [2, -1, 0],
      [-1, 2, -1],
      [0, -1, 2]
    ]
    const b = [1, 0, 1]
    const x = solve3x3(a, b)
    expect(x[0]).to.be.closeTo(1, Math.pow(10, -6))
    expect(x[1]).to.be.closeTo(1, Math.pow(10, -6))
    expect(x[2]).to.be.closeTo(1, Math.pow(10, -6))
  })

  it('solve3x3 throws SingularMatrixError for singular matrix', () => {
    const a = [
      [1, 2, 3],
      [2, 4, 6],
      [3, 6, 9]
    ]
    const b = [1, 2, 3]
    expect(() => solve3x3(a, b)).to.throw(SingularMatrixError)
  })

  it('tetDeterminant is positive for right-handed tet', () => {
    const v0 = [0, 0, 0]
    const v1 = [1, 0, 0]
    const v2 = [0, 1, 0]
    const v3 = [0, 0, 1]
    expect(tetDeterminant(v0, v1, v2, v3)).to.equal(1)
  })

  it('tetDeterminant is negative for left-handed tet', () => {
    const v0 = [0, 0, 0]
    const v1 = [1, 0, 0]
    const v2 = [0, 1, 0]
    const v3 = [0, 0, 1]
    // Swapping v1 and v2 gives negative orientation.
    expect(tetDeterminant(v0, v2, v1, v3)).to.equal(-1)
  })

  it('tetVolume of unit tet is 1/6', () => {
    const v0 = [0, 0, 0]
    const v1 = [1, 0, 0]
    const v2 = [0, 1, 0]
    const v3 = [0, 0, 1]
    expect(tetVolume(v0, v1, v2, v3)).to.be.closeTo(1 / 6, Math.pow(10, -10))
  })

  it('luSolve handles zero row via equilibration', () => {
    const a = [
      [0, 0, 0],
      [0, 1, 0],
      [0, 0, 1]
    ]
    const b = [0, 2, 3]
    expect(() => luSolve(a, b)).to.throw(SingularMatrixError)
  })

  it('factorial(0) is 1', () => {
    expect(factorial(0)).to.equal(1)
  })

  it('factorial(5) is 120', () => {
    expect(factorial(5)).to.equal(120)
  })

  it('factorial throws for negative input', () => {
    expect(() => factorial(-1)).to.throw(/non-negative/)
  })

  it('factorial throws for non-integer input', () => {
    expect(() => factorial(3.5)).to.throw(/non-negative/)
  })

  it('factorial throws for n > 170', () => {
    expect(() => factorial(171)).to.throw(/overflows/)
  })

  it('factorial(170) is finite', () => {
    expect(Number.isFinite(factorial(170))).to.equal(true)
  })

  it('numericalGradient approximates gradient of x^2 + y^2', () => {
    const u = (pt) => pt[0] * pt[0] + pt[1] * pt[1]
    const grad = numericalGradient(u, [2, 3, 0])
    expect(grad[0]).to.be.closeTo(4, Math.pow(10, -4))
    expect(grad[1]).to.be.closeTo(6, Math.pow(10, -4))
    expect(grad[2]).to.be.closeTo(0, Math.pow(10, -4))
  })
})
