import { expect } from 'chai'
import { Whitney } from '../src/lib/whitney.js'
import { Bcdtpp } from '../src/lib/bcdtpp.js'
import {
  generateUnitCubeMesh,
  generateSingleTetMesh
} from '../src/lib/mesh_generator.js'
import {
  computeL2ErrorScalar,
  computeL2ErrorVector,
  computeH1SemiError,
  estimateMeshSize,
  computeRate,
  runConvergenceStudy
} from '../src/lib/convergence_harness.js'

describe('Mesh Generator', () => {
  it('generates a valid single-cube mesh with 6 tets', () => {
    const mesh = generateUnitCubeMesh(1)
    expect(mesh.tetrahedronCount).to.equal(6)
    expect(mesh.vertexCount).to.equal(8)
    expect(mesh.boundaryFaces.length).to.be.above(0)
  })

  it('generates a valid 2x2x2 cube mesh with 48 tets', () => {
    const mesh = generateUnitCubeMesh(2)
    expect(mesh.tetrahedronCount).to.equal(48)
    expect(mesh.vertexCount).to.equal(27)
  })

  it('generates positively oriented tets', () => {
    // Mesh constructor validates orientation; no throw means success.
    expect(() => generateUnitCubeMesh(3)).to.not.throw()
  })

  it('rejects invalid n values', () => {
    expect(() => generateUnitCubeMesh(0)).to.throw()
    expect(() => generateUnitCubeMesh(-1)).to.throw()
    expect(() => generateUnitCubeMesh(1.5)).to.throw()
  })
})

describe('Convergence Harness Utilities', () => {
  it('estimateMeshSize returns positive for a cube mesh', () => {
    const mesh = generateUnitCubeMesh(2)
    const h = estimateMeshSize(mesh)
    expect(h).to.be.above(0)
    // For n=2, cube side = 0.5, but Freudenthal tets contain body
    // diagonals of length 0.5 * sqrt(3) ≈ 0.866.
    expect(h).to.be.closeTo(Math.sqrt(3) / 2, Math.pow(10, -6))
  })

  it('computeRate calculates expected slope', () => {
    // If error halves when h halves, rate = 1.
    expect(computeRate(0.5, 1.0, 0.5, 1.0)).to.be.closeTo(1, Math.pow(10, -6))
    // If error quarters when h halves, rate = 2.
    expect(computeRate(0.25, 1.0, 0.5, 1.0)).to.be.closeTo(2, Math.pow(10, -6))
  })

  it('computeRate returns 0 when err2 is zero', () => {
    expect(computeRate(0.5, 0, 0.5, 1.0)).to.equal(0)
  })

  it('computeRate returns 0 when h2 equals h1', () => {
    expect(computeRate(0.5, 1.0, 0.5, 0.5)).to.equal(0)
  })

  it('computeRate returns 0 when err1 is zero', () => {
    expect(computeRate(0, 1.0, 0.5, 1.0)).to.equal(0)
  })

  it('computeRate returns 0 when h1 is zero', () => {
    expect(computeRate(0.5, 1.0, 0, 1.0)).to.equal(0)
  })

  it('computeRate returns 0 when h2 is zero', () => {
    expect(computeRate(0.5, 1.0, 0.5, 0)).to.equal(0)
  })

  it('runConvergenceStudy returns results for scalar l=0', () => {
    const meshes = [1, 2].map((n) => generateUnitCubeMesh(n))
    const results = runConvergenceStudy(meshes, {
      exactScalar: (pt) => pt[0] + pt[1] + pt[2],
      l: 0,
      p: 0,
      quadratureOrder: 3
    })
    expect(results.length).to.equal(2)
    expect(results[0]).to.have.property('h')
    expect(results[0]).to.have.property('l2Err')
    expect(results[0]).to.have.property('h1Err')
    expect(results[0]).to.have.property('rateL2')
    expect(results[0]).to.have.property('rateH1')
    expect(results[0].rateL2).to.equal(undefined)
    expect(typeof results[1].rateL2).to.equal('number')
  })

  it('runConvergenceStudy returns results for vector l=1', () => {
    const meshes = [1, 2].map((n) => generateUnitCubeMesh(n))
    const results = runConvergenceStudy(meshes, {
      exactScalar: () => 0,
      exactVector: (pt) => [pt[0], pt[1], pt[2]],
      l: 1,
      p: 0,
      quadratureOrder: 3
    })
    expect(results.length).to.equal(2)
    expect(results[0]).to.have.property('l2Err')
    expect(results[0].h1Err).to.equal(undefined)
  })

  it('runConvergenceStudy throws for vector l without exactVector', () => {
    const meshes = [1, 2].map((n) => generateUnitCubeMesh(n))
    expect(() =>
      runConvergenceStudy(meshes, {
        exactScalar: () => 0,
        l: 1,
        p: 0
      })
    ).to.throw(/exactVector is required/)
  })

  it('runConvergenceStudy works for l=3', () => {
    const meshes = [1, 2].map((n) => generateUnitCubeMesh(n))
    const results = runConvergenceStudy(meshes, {
      exactScalar: (pt) => pt[0] + pt[1] + pt[2],
      l: 3,
      p: 0,
      quadratureOrder: 3
    })
    expect(results.length).to.equal(2)
    expect(results[0]).to.have.property('l2Err')
    expect(results[0].h1Err).to.equal(undefined)
  })
})

describe('p-Convergence on Single Tet', () => {
  const mesh = generateSingleTetMesh()
  const whitney = new Whitney(mesh)
  const bcdtpp = new Bcdtpp(mesh, whitney, { quadratureOrder: 3 })
  bcdtpp.computeBoundaryWeights()

  // Quadratic exact function: u(x,y,z) = x^2 + y^2 + z^2.
  // P^1 cannot represent quadratics exactly; P^2 can.
  const exactQuadratic = (pt) =>
    pt[0] * pt[0] + pt[1] * pt[1] + pt[2] * pt[2]

  it('l=0 p=1 L2 error is non-zero for quadratic', () => {
    const projFn = (tIdx, pt) => bcdtpp.projectHp(exactQuadratic, pt, tIdx, 0, 1)
    const err = computeL2ErrorScalar(mesh, bcdtpp, exactQuadratic, projFn)
    expect(err).to.be.above(1e-3)
  })

  it('l=0 p=2 L2 error is smaller than p=1 for quadratic', () => {
    const projP1 = (tIdx, pt) => bcdtpp.projectHp(exactQuadratic, pt, tIdx, 0, 1)
    const projP2 = (tIdx, pt) => bcdtpp.projectHp(exactQuadratic, pt, tIdx, 0, 2)
    const errP1 = computeL2ErrorScalar(mesh, bcdtpp, exactQuadratic, projP1)
    const errP2 = computeL2ErrorScalar(mesh, bcdtpp, exactQuadratic, projP2)
    expect(errP2).to.be.below(errP1)
  })

  it('l=3 p=2 differs from p=0 for quadratic at barycenter', () => {
    const pt = mesh.getTetrahedronBarycenter(0)
    const p0 = bcdtpp.projectHp(exactQuadratic, pt, 0, 3, 0)
    const p2 = bcdtpp.projectHp(exactQuadratic, pt, 0, 3, 2)
    expect(typeof p2).to.equal('number')
    expect(Number.isFinite(p2)).to.equal(true)
    expect(p2).to.not.equal(p0)
  })
})

describe('h-Convergence on Cube Meshes', () => {
  const exactScalar = (pt) =>
    Math.sin(Math.PI * pt[0]) *
    Math.sin(Math.PI * pt[1]) *
    Math.sin(Math.PI * pt[2])

  const exactVector = (pt) => [
    Math.sin(Math.PI * pt[0]),
    Math.sin(Math.PI * pt[1]),
    Math.sin(Math.PI * pt[2])
  ]

  it('l=0 p=0: H1 semi-norm converges at rate ~1 (asymptotic)', () => {
    // Note: on coarse all-boundary meshes (open cube), boundary integral
    // weights are not yet exact, so pre-asymptotic rates are lower.
    const meshes = [1, 2, 4].map((n) => generateUnitCubeMesh(n))
    const results = []

    for (const mesh of meshes) {
      const whitney = new Whitney(mesh)
      const bcdtpp = new Bcdtpp(mesh, whitney, { quadratureOrder: 3 })
      bcdtpp.computeBoundaryWeights()

      const projFn = (tIdx, pt) => bcdtpp.projectHp(exactScalar, pt, tIdx, 0, 0)
      const l2Err = computeL2ErrorScalar(mesh, bcdtpp, exactScalar, projFn)
      const h1Err = computeH1SemiError(mesh, bcdtpp, exactScalar, projFn)
      const h = estimateMeshSize(mesh)
      results.push({ h, l2Err, h1Err })
    }

    // H1 semi-norm rate approaches 1 as the mesh refines.
    const rateH1One = computeRate(
      results[1].h1Err,
      results[0].h1Err,
      results[1].h,
      results[0].h
    )
    const rateH1Two = computeRate(
      results[2].h1Err,
      results[1].h1Err,
      results[2].h,
      results[1].h
    )

    expect(rateH1One).to.be.above(0.3)
    expect(rateH1Two).to.be.above(0.6)
  })

  it('l=0 p=0: L2 error converges at rate ~1+ (asymptotic)', () => {
    // On all-boundary meshes, boundary-weighted vertex DoFs deviate from
    // exact nodal values on coarse meshes, lowering pre-asymptotic rates.
    const meshes = [1, 2, 4].map((n) => generateUnitCubeMesh(n))
    const results = []

    for (const mesh of meshes) {
      const whitney = new Whitney(mesh)
      const bcdtpp = new Bcdtpp(mesh, whitney, { quadratureOrder: 3 })
      bcdtpp.computeBoundaryWeights()

      const projFn = (tIdx, pt) => bcdtpp.projectHp(exactScalar, pt, tIdx, 0, 0)
      const l2Err = computeL2ErrorScalar(mesh, bcdtpp, exactScalar, projFn)
      const h = estimateMeshSize(mesh)
      results.push({ h, l2Err })
    }

    const rateL2One = computeRate(
      results[1].l2Err,
      results[0].l2Err,
      results[1].h,
      results[0].h
    )
    const rateL2Two = computeRate(
      results[2].l2Err,
      results[1].l2Err,
      results[2].h,
      results[1].h
    )

    expect(rateL2One).to.be.above(0.8)
    expect(rateL2Two).to.be.above(1.2)
  })

  it('l=3 p=0: L2 error converges at rate ~1', () => {
    const meshes = [1, 2, 4].map((n) => generateUnitCubeMesh(n))
    const results = []

    for (const mesh of meshes) {
      const whitney = new Whitney(mesh)
      const bcdtpp = new Bcdtpp(mesh, whitney, { quadratureOrder: 3 })
      bcdtpp.computeBoundaryWeights()

      const projFn = (tIdx, pt) => bcdtpp.projectHp(exactScalar, pt, tIdx, 3, 0)
      const l2Err = computeL2ErrorScalar(mesh, bcdtpp, exactScalar, projFn)
      const h = estimateMeshSize(mesh)
      results.push({ h, l2Err })
    }

    const rateL2One = computeRate(
      results[1].l2Err,
      results[0].l2Err,
      results[1].h,
      results[0].h
    )
    const rateL2Two = computeRate(
      results[2].l2Err,
      results[1].l2Err,
      results[2].h,
      results[1].h
    )

    // Piecewise constant L2 projection has L2 error O(h).
    expect(rateL2One).to.be.above(0.6)
    expect(rateL2Two).to.be.above(0.6)
  })

  it('l=1 p=0: L2 error converges at rate ~1 for vector field', () => {
    const meshes = [1, 2, 4].map((n) => generateUnitCubeMesh(n))
    const results = []

    for (const mesh of meshes) {
      const whitney = new Whitney(mesh)
      const bcdtpp = new Bcdtpp(mesh, whitney, { quadratureOrder: 3 })
      bcdtpp.computeBoundaryWeights()

      const projFn = (tIdx, pt) => bcdtpp.projectHp(exactVector, pt, tIdx, 1, 0)
      const l2Err = computeL2ErrorVector(mesh, bcdtpp, exactVector, projFn)
      const h = estimateMeshSize(mesh)
      results.push({ h, l2Err })
    }

    const rateL2One = computeRate(
      results[1].l2Err,
      results[0].l2Err,
      results[1].h,
      results[0].h
    )
    const rateL2Two = computeRate(
      results[2].l2Err,
      results[1].l2Err,
      results[2].h,
      results[1].h
    )

    // Lowest-order Nedelec has L2 error O(h).
    expect(rateL2One).to.be.above(0.5)
    expect(rateL2Two).to.be.above(0.8)
  })
})
