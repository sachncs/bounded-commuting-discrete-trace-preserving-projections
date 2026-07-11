/**
 * Broad coverage and edge-case tests: multi-tet projection correctness,
 * mesh validation errors, higher-order edge cases, singular-matrix handling,
 * BoundaryWeightComputer fault isolation, and scalar-to-vector projection
 * interoperability.
 */
import { expect } from 'chai'
import sinon from 'sinon'
import { Mesh } from '../src/lib/mesh.js'
import { Whitney } from '../src/lib/whitney.js'
import { Bcdtpp } from '../src/lib/bcdtpp.js'
import { HigherOrderProjection } from '../src/lib/higher_order_projection.js'
import { BoundaryWeightComputer } from '../src/lib/boundary_weight_computer.js'
import { MeshRefinement } from '../src/lib/mesh_refinement.js'
import { MeshValidationError } from '../src/lib/errors.js'
import { factorial } from '../src/lib/math_utils.js'
import { generateUnitCubeMesh } from '../src/lib/mesh_generator.js'

// Multi-tet mesh: exercises projection across shared faces/edges and
// verifies commuting properties hold element-by-element.
describe('Multi-tet Bcdtpp Projections', () => {
  const twoTet = {
    vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]],
    tetrahedra: [[0, 1, 2, 3], [1, 2, 3, 4]]
  }

  let mesh
  let whitney
  let bcdtpp

  before(() => {
    mesh = new Mesh(twoTet.vertices, twoTet.tetrahedra)
    whitney = new Whitney(mesh)
    bcdtpp = new Bcdtpp(mesh, whitney, { quadratureOrder: 3 })
    bcdtpp.computeBoundaryWeights()
    bcdtpp.buildPointLocator()
  })

  it('projectH1 is exact for constant on multi-tet mesh', () => {
    const u = () => 5
    for (let tIdx = 0; tIdx < mesh.tetrahedronCount; tIdx++) {
      const pt = mesh.getTetrahedronBarycenter(tIdx)
      const result = bcdtpp.projectH1(u, pt, tIdx)
      expect(result).to.be.closeTo(5, Math.pow(10, -6))
    }
  })

  it('projectL2 integrates constants exactly on multi-tet mesh', () => {
    const u = () => 7
    for (let tIdx = 0; tIdx < mesh.tetrahedronCount; tIdx++) {
      const result = bcdtpp.projectL2(u, tIdx)
      expect(result).to.be.closeTo(7, Math.pow(10, -6))
    }
  })

  it('projectAtPoint finds correct tet in multi-tet mesh', () => {
    const u = (pt) => pt[0] + pt[1] + pt[2]
    const result = bcdtpp.projectAtPoint(u, [0.1, 0.1, 0.1], 0)
    expect(result.tIdx).to.equal(0)
    expect(Number.isFinite(result.value)).to.equal(true)
  })

  it('projectHdiv returns finite vector on multi-tet mesh', () => {
    const v = (pt) => [2 * pt[0], 2 * pt[1], 2 * pt[2]]
    for (let tIdx = 0; tIdx < mesh.tetrahedronCount; tIdx++) {
      const pt = mesh.getTetrahedronBarycenter(tIdx)
      const projV = bcdtpp.projectHdiv(v, pt, tIdx)
      expect(Array.isArray(projV)).to.equal(true)
      expect(projV.length).to.equal(3)
      projV.forEach((c) => expect(Number.isFinite(c)).to.equal(true))
    }
  })

  // Tolerance 1e-1: numerical divergence via finite differences (h=1e-5)
  // on the projected field introduces O(h) discretization error; 1e-1
  // accommodates the compounded error from projection + finite difference.
  it('commutes div Pi^2 = Pi^3 div on multi-tet mesh', () => {
    const v = (pt) => [2 * pt[0], 2 * pt[1], 2 * pt[2]]
    const divV = () => 6
    const h = 1e-5
    for (let tIdx = 0; tIdx < mesh.tetrahedronCount; tIdx++) {
      const pt = mesh.getTetrahedronBarycenter(tIdx)
      const proj = bcdtpp.projectHdiv(v, pt, tIdx)
      const projVx = bcdtpp.projectHdiv(v, [pt[0] + h, pt[1], pt[2]], tIdx)
      const projVy = bcdtpp.projectHdiv(v, [pt[0], pt[1] + h, pt[2]], tIdx)
      const projVz = bcdtpp.projectHdiv(v, [pt[0], pt[1], pt[2] + h], tIdx)
      const numDiv =
        (projVx[0] - proj[0]) / h +
        (projVy[1] - proj[1]) / h +
        (projVz[2] - proj[2]) / h
      const l2Div = bcdtpp.projectL2(divV, tIdx)
      expect(numDiv).to.be.closeTo(l2Div, Math.pow(10, -1))
    }
  })

  it('Pi^1 reproduces gradient of linear function on multi-tet mesh', () => {
    const u = (pt) => 2 * pt[0] - 3 * pt[1] + 5 * pt[2]
    const gradU = [2, -3, 5]
    for (let tIdx = 0; tIdx < mesh.tetrahedronCount; tIdx++) {
      const pt = mesh.getTetrahedronBarycenter(tIdx)
      const proj = bcdtpp.projectHcurl(u, pt, tIdx)
      expect(proj[0]).to.be.closeTo(gradU[0], Math.pow(10, -5))
      expect(proj[1]).to.be.closeTo(gradU[1], Math.pow(10, -5))
      expect(proj[2]).to.be.closeTo(gradU[2], Math.pow(10, -5))
    }
  })
})

// Verifies that Mesh throws MeshValidationError for various invalid inputs.
describe('Mesh Validation', () => {
  it('throws MeshValidationError for negative orientation', () => {
    expect(() => new Mesh(
      [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
      [[0, 2, 1, 3]] // inverted: 0,2,1 instead of 0,1,2
    )).to.throw(MeshValidationError)
  })

  it('throws MeshValidationError for out-of-bounds index', () => {
    expect(() => new Mesh(
      [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
      [[0, 1, 2, 4]]
    )).to.throw(MeshValidationError)
  })

  it('throws MeshValidationError for duplicate vertices in tet', () => {
    expect(() => new Mesh(
      [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
      [[0, 1, 1, 3]]
    )).to.throw(MeshValidationError)
  })

  it('throws MeshValidationError for non-finite vertex', () => {
    expect(() => new Mesh(
      [[0, 0, NaN], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
      [[0, 1, 2, 3]]
    )).to.throw(MeshValidationError)
  })

  it('throws MeshValidationError for non-integer tetrahedron index', () => {
    expect(() => new Mesh(
      [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
      [[0, 1, 2, 3.5]]
    )).to.throw(MeshValidationError)
  })

  it('throws MeshValidationError for empty tetrahedra', () => {
    expect(() => new Mesh([[0, 0, 0]], [])).to.throw(MeshValidationError)
  })

  it('rejects degenerate (zero-volume) tetrahedra', () => {
    expect(() => new Mesh(
      [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]],
      [[0, 1, 2, 3]]
    )).to.throw(MeshValidationError)
  })
})

// Tests edge cases in the higher-order projection framework: p=1 fallback
// to lowest-order, p=2,3 L2 enrichment for H^1, and unimplemented l=1,2 p>0.
describe('Higher-Order Edge Cases', () => {
  const singleTet = {
    vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
    tetrahedra: [[0, 1, 2, 3]]
  }

  let mesh
  let whitney
  let bcdtpp

  before(() => {
    mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    whitney = new Whitney(mesh)
    bcdtpp = new Bcdtpp(mesh, whitney, { quadratureOrder: 3 })
    bcdtpp.computeBoundaryWeights()
  })

  it('p=1 H1 returns exact base projection', () => {
    const u = (pt) => pt[0] * pt[0] + pt[1] * pt[1]
    const pt = mesh.getTetrahedronBarycenter(0)
    const p0 = bcdtpp.projectHp(u, pt, 0, 0, 0)
    const p1 = bcdtpp.projectHp(u, pt, 0, 0, 1)
    expect(p1).to.equal(p0)
  })

  it('p=2,3 H1 uses L2 enrichment', () => {
    const u = (pt) => pt[0] * pt[0] + pt[1] * pt[1]
    const pt = mesh.getTetrahedronBarycenter(0)
    const p0 = bcdtpp.projectHp(u, pt, 0, 0, 0)
    const p2 = bcdtpp.projectHp(u, pt, 0, 0, 2)
    const p3 = bcdtpp.projectHp(u, pt, 0, 0, 3)
    expect(typeof p2).to.equal('number')
    expect(Number.isFinite(p2)).to.equal(true)
    expect(typeof p3).to.equal('number')
    expect(Number.isFinite(p3)).to.equal(true)
    // Higher-order should differ from lowest-order for non-linear functions.
    expect(p2).to.not.equal(p0)
    expect(p3).to.not.equal(p0)
  })

  it('throws for unimplemented l=1, p>0', () => {
    const u = (pt) => [pt[0], pt[1], pt[2]]
    const pt = mesh.getTetrahedronBarycenter(0)
    expect(() => bcdtpp.projectHp(u, pt, 0, 1, 1)).to.throw(/not yet implemented/)
  })

  it('throws for unimplemented l=2, p>0', () => {
    const u = (pt) => [pt[0], pt[1], pt[2]]
    const pt = mesh.getTetrahedronBarycenter(0)
    expect(() => bcdtpp.projectHp(u, pt, 0, 2, 1)).to.throw(/not yet implemented/)
  })
})

// Verifies ProjectionError is thrown for invalid arguments to projection methods.
describe('Projection Error Handling', () => {
  const singleTet = {
    vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
    tetrahedra: [[0, 1, 2, 3]]
  }

  let mesh
  let whitney
  let bcdtpp

  before(() => {
    mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    whitney = new Whitney(mesh)
    bcdtpp = new Bcdtpp(mesh, whitney, { quadratureOrder: 3 })
    bcdtpp.computeBoundaryWeights()
    bcdtpp.buildPointLocator()
  })

  it('throws for invalid tIdx type', () => {
    expect(() => bcdtpp.projectH1(() => 1, [0, 0, 0], '0')).to.throw()
  })

  it('throws for out-of-range tIdx', () => {
    expect(() => bcdtpp.projectH1(() => 1, [0, 0, 0], 5)).to.throw()
  })

  it('throws for negative tIdx', () => {
    expect(() => bcdtpp.projectH1(() => 1, [0, 0, 0], -1)).to.throw()
  })

  it('throws for invalid point type', () => {
    expect(() => bcdtpp.projectH1(() => 1, 'bad', 0)).to.throw(/point must be/)
  })

  it('throws for NaN in point', () => {
    expect(() => bcdtpp.projectH1(() => 1, [0, NaN, 0], 0)).to.throw(/point must be/)
  })
})

// Tests the HigherOrderProjection bubble solve path, including singular
// matrix fault injection via sinon stubs.
describe('Bubble Projection', () => {
  afterEach(() => sinon.restore())

  const singleTet = {
    vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
    tetrahedra: [[0, 1, 2, 3]]
  }

  let mesh
  let whitney

  before(() => {
    mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    whitney = new Whitney(mesh)
  })

  it('solveBubbleProjection returns null for p < 4', () => {
    const hop = new HigherOrderProjection(mesh, whitney)
    const result = hop.solveBubbleProjection(0, 3, () => 1)
    expect(result).to.equal(null)
  })

  it('factorial overflow throws for n > 170', () => {
    expect(() => factorial(171)).to.throw(/overflows/)
  })

  it('solveBubbleProjection computes coefficients for p = 4', () => {
    const hop = new HigherOrderProjection(mesh, whitney)
    const coeffs = hop.solveBubbleProjection(0, 4, (pt) => pt[0] * pt[0])
    expect(coeffs).to.not.equal(null)
    expect(coeffs.length).to.be.above(0)
  })

  it('evaluateBubble returns 0 when coeffs are empty', () => {
    const hop = new HigherOrderProjection(mesh, whitney)
    const val = hop.evaluateBubble(0, 4, [], [0.1, 0.1, 0.1])
    expect(val).to.equal(0)
  })

  it('evaluateBubble evaluates a bubble correction', () => {
    const hop = new HigherOrderProjection(mesh, whitney)
    const coeffs = hop.solveBubbleProjection(0, 4, (pt) => pt[0] * pt[0])
    const pt = [0.1, 0.1, 0.1]
    const val = hop.evaluateBubble(0, 4, coeffs, pt)
    expect(Number.isFinite(val)).to.equal(true)
  })

  it('evaluateL2Projection returns 0 for empty coeffs', () => {
    const hop = new HigherOrderProjection(mesh, whitney)
    const val = hop.evaluateL2Projection([], [0.25, 0.25, 0.25, 0.25], 1)
    expect(val).to.equal(0)
  })

  it('solveL2Projection returns empty array for negative p', () => {
    const hop = new HigherOrderProjection(mesh, whitney)
    const coeffs = hop.solveL2Projection(0, -1, () => 1)
    expect(coeffs).to.deep.equal([])
  })

  it('solveBubbleProjection warns on singular mass matrix', () => {
    const warnSpy = sinon.spy()
    const mesh = new Mesh(
      [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
      [[0, 1, 2, 3]]
    )
    const whitney = new Whitney(mesh)
    const hop = new HigherOrderProjection(mesh, whitney, 3, warnSpy)
    // Inject a singular mass matrix to force the solver to fail.
    hop.assembleBubbleMass = () => [[0]]
    const coeffs = hop.solveBubbleProjection(0, 4, () => 1)
    expect(coeffs).to.equal(null)
    expect(warnSpy.called).to.equal(true)
  })

  it('solveL2Projection warns on singular mass matrix', () => {
    const warnSpy = sinon.spy()
    const mesh = new Mesh(
      [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
      [[0, 1, 2, 3]]
    )
    const whitney = new Whitney(mesh)
    const hop = new HigherOrderProjection(mesh, whitney, 3, warnSpy)
    // Force zero volume to make the mass matrix singular.
    const origGetVolume = mesh.getVolume.bind(mesh)
    mesh.getVolume = () => 0
    const coeffs = hop.solveL2Projection(0, 2, () => 1)
    expect(coeffs).to.deep.equal([])
    expect(warnSpy.called).to.equal(true)
    mesh.getVolume = origGetVolume
  })
})

// Verifies that BoundaryWeightComputer isolates per-vertex failures:
// zero star area and adjacency exceptions emit warnings instead of throwing.
describe('BoundaryWeightComputer fault isolation', () => {
  afterEach(() => sinon.restore())

  it('warns and skips vertex with zero star area', () => {
    const warnSpy = sinon.spy()
    const mesh = new Mesh(
      [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
      [[0, 1, 2, 3]]
    )
    const refinement = new MeshRefinement(mesh)
    refinement.computeWorseyFarinSplit()
    // Force all face areas to zero so starArea falls below the threshold.
    const origGetFaceArea = mesh.getFaceArea.bind(mesh)
    mesh.getFaceArea = (fIdx) => 0
    const bwc = new BoundaryWeightComputer(mesh, refinement, warnSpy)
    bwc.compute()
    expect(warnSpy.called).to.equal(true)
    expect(warnSpy.getCall(0).args[0].code).to.equal('BWC_ZERO_STAR_AREA')
    expect(warnSpy.getCall(0).args[0].message).to.match(/star area/)
    // Restore for downstream tests.
    mesh.getFaceArea = origGetFaceArea
  })

  it('catches per-vertex solver failures and warns', () => {
    const warnSpy = sinon.spy()
    const mesh = new Mesh(
      [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
      [[0, 1, 2, 3]]
    )
    const refinement = new MeshRefinement(mesh)
    refinement.computeWorseyFarinSplit()
    const bwc = new BoundaryWeightComputer(mesh, refinement, warnSpy)
    // Inject a fault by temporarily corrupting the mesh adjacency.
    const orig = mesh.getBoundaryStar
    mesh.getBoundaryStar = (vIdx) => {
      throw new Error('simulated adjacency failure')
    }
    expect(() => bwc.compute()).to.not.throw()
    expect(warnSpy.called).to.equal(true)
    expect(warnSpy.getCall(0).args[0].code).to.equal('BWC_VERTEX_FAILURE')
    expect(warnSpy.getCall(0).args[0].message).to.match(/failed to compute weights/)
    mesh.getBoundaryStar = orig
  })
})

// Tests that H(curl) and H(div) projectors accept scalar functions
// (which are internally promoted to vector fields via numerical gradient).
describe('Scalar inputs to vector projections', () => {
  const cubeMesh = generateUnitCubeMesh(2)
  let whitney
  let bcdtpp

  before(() => {
    whitney = new Whitney(cubeMesh)
    bcdtpp = new Bcdtpp(cubeMesh, whitney, { quadratureOrder: 3 })
    bcdtpp.computeBoundaryWeights()
  })

  it('projectHcurl accepts scalar function on multi-tet mesh', () => {
    const u = (pt) => pt[0] + pt[1] + pt[2]
    const pt = [0.5, 0.5, 0.5]
    const result = bcdtpp.projectAtPoint(u, pt, 1)
    expect(Array.isArray(result.value)).to.equal(true)
    expect(result.value.length).to.equal(3)
    result.value.forEach((v) => expect(Number.isFinite(v)).to.equal(true))
  })

  it('projectHdiv accepts scalar function on multi-tet mesh', () => {
    const u = (pt) => pt[0] + pt[1] + pt[2]
    const pt = [0.5, 0.5, 0.5]
    const result = bcdtpp.projectAtPoint(u, pt, 2)
    expect(Array.isArray(result.value)).to.equal(true)
    expect(result.value.length).to.equal(3)
    result.value.forEach((v) => expect(Number.isFinite(v)).to.equal(true))
  })
})

// Verifies that Mesh rejects degenerate (coplanar/collinear) tetrahedra.
describe('Whitney degenerate tet handling', () => {
  it('rejects degenerate tet at mesh construction', () => {
    expect(() => new Mesh(
      [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]],
      [[0, 1, 2, 3]]
    )).to.throw(MeshValidationError)
  })
})
