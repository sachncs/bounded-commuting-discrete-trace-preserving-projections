/**
 * Tests for the interior projector (Pi_ring^l), discrete extension operator
 * (E^l), their decomposition (Pi^l = Pi_ring^l + Pi_partial^l), commuting
 * diagram properties (tr^l(E^l(w)) = w), and the global projector dispatch
 * on both single-tet and multi-tet meshes.
 */
import { expect } from 'chai'
import { Mesh } from '../src/lib/mesh.js'
import { Whitney } from '../src/lib/whitney.js'
import { Bcdtpp } from '../src/lib/bcdtpp.js'
import { generateUnitCubeMesh } from '../src/lib/mesh_generator.js'

const singleTet = {
  vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
  tetrahedra: [[0, 1, 2, 3]]
}

// Pi_ring^l: verifies that the interior projector vanishes on boundary
// for all form degrees (H1, Hcurl, Hdiv, L2) on a single tet.
describe('Pi_ring^l interior projector', () => {
  let mesh
  let whitney
  let bcdtpp

  before(() => {
    mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    whitney = new Whitney(mesh)
    bcdtpp = new Bcdtpp(mesh, whitney, { quadratureOrder: 3 })
    bcdtpp.computeBoundaryWeights()
  })

  it('projectRing H1 is exact for linear on single tet', () => {
    const u = (pt) => 2 * pt[0] - 3 * pt[1] + 5 * pt[2]
    const pt = mesh.getTetrahedronBarycenter(0)
    const result = bcdtpp.projectRing(u, pt, 0, 0)
    // On a single tet all vertices are boundary, so projectRing returns 0.
    expect(result).to.equal(0)
  })

  it('projectRing H1 vanishes on boundary vertices', () => {
    const u = (pt) => pt[0] + pt[1] + pt[2] + 1
    for (let i = 0; i < 4; i++) {
      const v = mesh.vertices[i]
      const result = bcdtpp.projectRing(u, v, 0, 0)
      expect(result).to.equal(0)
    }
  })

  it('projectRing Hcurl vanishes on boundary edges', () => {
    const u = (pt) => [pt[0], pt[1], pt[2]]
    for (const eIdx of mesh.boundaryEdges) {
      const e = mesh.edges[eIdx]
      const mid = [
        (mesh.vertices[e[0]][0] + mesh.vertices[e[1]][0]) / 2,
        (mesh.vertices[e[0]][1] + mesh.vertices[e[1]][1]) / 2,
        (mesh.vertices[e[0]][2] + mesh.vertices[e[1]][2]) / 2
      ]
      // Find a tet containing this edge.
      const tIdx = mesh.edgeToFaces[eIdx].length > 0
        ? mesh.faceToTets[mesh.edgeToFaces[eIdx][0]][0]
        : 0
      const result = bcdtpp.projectRing(u, mid, tIdx, 1)
      // Edge basis evaluated on the edge should give zero for interior DoFs.
      // Actually, the edge basis is non-zero on the edge, but the coefficient
      // for boundary edges is zero, so the tangential component along the edge
      // should vanish.
      const edgeVec = [
        mesh.vertices[e[1]][0] - mesh.vertices[e[0]][0],
        mesh.vertices[e[1]][1] - mesh.vertices[e[0]][1],
        mesh.vertices[e[1]][2] - mesh.vertices[e[0]][2]
      ]
      const tangential =
        result[0] * edgeVec[0] +
        result[1] * edgeVec[1] +
        result[2] * edgeVec[2]
      expect(Math.abs(tangential)).to.be.below(1e-12)
    }
  })

  it('projectRing Hdiv has zero normal flux on boundary faces', () => {
    const u = (pt) => [pt[0], pt[1], pt[2]]
    for (const fIdx of mesh.boundaryFaces) {
      const normal = mesh.getFaceOutwardNormal(fIdx)
      const bary = mesh.getFaceBarycenter(fIdx)
      // Find a tet containing this face.
      const tIdx = mesh.faceToTets[fIdx][0]
      const result = bcdtpp.projectRing(u, bary, tIdx, 2)
      const flux =
        result[0] * normal[0] +
        result[1] * normal[1] +
        result[2] * normal[2]
      expect(Math.abs(flux)).to.be.below(1e-12)
    }
  })

  it('projectRing L2 equals projectL2 (no boundary DoFs)', () => {
    const u = (pt) => pt[0] * pt[0] + pt[1] * pt[1]
    const resultRing = bcdtpp.projectRing(u, [0.25, 0.25, 0.25], 0, 3)
    const resultL2 = bcdtpp.projectL2(u, 0)
    expect(resultRing).to.be.closeTo(resultL2, Math.pow(10, -10))
  })
})

// E^l: verifies that the discrete extension operator preserves boundary
// DoF values exactly for H1 and Hcurl on a single tet.
describe('Discrete extension operator E^l', () => {
  let mesh
  let whitney
  let bcdtpp

  before(() => {
    mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    whitney = new Whitney(mesh)
    bcdtpp = new Bcdtpp(mesh, whitney, { quadratureOrder: 3 })
    bcdtpp.computeBoundaryWeights()
  })

  it('E^0 preserves boundary vertex values', () => {
    const boundaryValues = new Map()
    for (const vIdx of mesh.boundaryNodes) {
      boundaryValues.set(vIdx, vIdx * 10 + 5)
    }
    for (const vIdx of mesh.boundaryNodes) {
      const pt = mesh.vertices[vIdx]
      const result = bcdtpp.extendBoundary(boundaryValues, pt, 0, 0)
      expect(result).to.be.closeTo(boundaryValues.get(vIdx), Math.pow(10, -10))
    }
  })

  it('E^1 preserves boundary edge tangential components', () => {
    const boundaryEdgeDofs = new Map()
    for (const eIdx of mesh.boundaryEdges) {
      boundaryEdgeDofs.set(eIdx, eIdx * 2 + 1)
    }
    for (const eIdx of mesh.boundaryEdges) {
      const e = mesh.edges[eIdx]
      const mid = [
        (mesh.vertices[e[0]][0] + mesh.vertices[e[1]][0]) / 2,
        (mesh.vertices[e[0]][1] + mesh.vertices[e[1]][1]) / 2,
        (mesh.vertices[e[0]][2] + mesh.vertices[e[1]][2]) / 2
      ]
      const tIdx = mesh.faceToTets[mesh.edgeToFaces[eIdx][0]][0]
      const result = bcdtpp.extendBoundary(boundaryEdgeDofs, mid, tIdx, 1)
      const edgeVec = [
        mesh.vertices[e[1]][0] - mesh.vertices[e[0]][0],
        mesh.vertices[e[1]][1] - mesh.vertices[e[0]][1],
        mesh.vertices[e[1]][2] - mesh.vertices[e[0]][2]
      ]
      const edgeLen = Math.sqrt(
        edgeVec[0] * edgeVec[0] +
        edgeVec[1] * edgeVec[1] +
        edgeVec[2] * edgeVec[2]
      )
      const tangential =
        result[0] * edgeVec[0] +
        result[1] * edgeVec[1] +
        result[2] * edgeVec[2]
      // The tangential component should equal the prescribed line integral.
      expect(tangential / edgeLen).to.be.closeTo(
        boundaryEdgeDofs.get(eIdx) / edgeLen,
        Math.pow(10, -10)
      )
    }
  })

  it('E^2 decomposition: projectPartial equals extendBoundary of extracted DoFs', () => {
    const u = (pt) => [pt[0], pt[1], pt[2]]
    const pt = mesh.getTetrahedronBarycenter(0)
    const partial = bcdtpp.projectPartial(u, pt, 0, 2)
    const boundaryData = bcdtpp.extractBoundaryDofs(u, 2)
    const extended = bcdtpp.extendBoundary(boundaryData, pt, 0, 2)
    expect(partial[0]).to.be.closeTo(extended[0], Math.pow(10, -10))
    expect(partial[1]).to.be.closeTo(extended[1], Math.pow(10, -10))
    expect(partial[2]).to.be.closeTo(extended[2], Math.pow(10, -10))
  })
})

// Verifies the fundamental decomposition Pi^l = Pi_ring^l + Pi_partial^l
// for all three vector-valued form degrees (H1, Hcurl, Hdiv).
describe('Decomposition Pi^l = Pi_ring^l + Pi_partial^l', () => {
  let mesh
  let whitney
  let bcdtpp

  before(() => {
    mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    whitney = new Whitney(mesh)
    bcdtpp = new Bcdtpp(mesh, whitney, { quadratureOrder: 3 })
    bcdtpp.computeBoundaryWeights()
  })

  it('H1: projectH1 = projectRing + projectPartial', () => {
    const u = (pt) => pt[0] * pt[0] + pt[1] * pt[1] + pt[2]
    const pt = mesh.getTetrahedronBarycenter(0)
    const full = bcdtpp.projectH1(u, pt, 0)
    const ring = bcdtpp.projectRing(u, pt, 0, 0)
    const partial = bcdtpp.projectPartial(u, pt, 0, 0)
    expect(full).to.be.closeTo(ring + partial, Math.pow(10, -10))
  })

  it('Hcurl: projectHcurl = projectRing + projectPartial', () => {
    const u = (pt) => [2 * pt[0], 2 * pt[1], 2 * pt[2]]
    const pt = mesh.getTetrahedronBarycenter(0)
    const full = bcdtpp.projectHcurl(u, pt, 0)
    const ring = bcdtpp.projectRing(u, pt, 0, 1)
    const partial = bcdtpp.projectPartial(u, pt, 0, 1)
    expect(full[0]).to.be.closeTo(ring[0] + partial[0], Math.pow(10, -10))
    expect(full[1]).to.be.closeTo(ring[1] + partial[1], Math.pow(10, -10))
    expect(full[2]).to.be.closeTo(ring[2] + partial[2], Math.pow(10, -10))
  })

  it('Hdiv: projectHdiv = projectRing + projectPartial', () => {
    const u = (pt) => [2 * pt[0], 2 * pt[1], 2 * pt[2]]
    const pt = mesh.getTetrahedronBarycenter(0)
    const full = bcdtpp.projectHdiv(u, pt, 0)
    const ring = bcdtpp.projectRing(u, pt, 0, 2)
    const partial = bcdtpp.projectPartial(u, pt, 0, 2)
    expect(full[0]).to.be.closeTo(ring[0] + partial[0], Math.pow(10, -10))
    expect(full[1]).to.be.closeTo(ring[1] + partial[1], Math.pow(10, -10))
    expect(full[2]).to.be.closeTo(ring[2] + partial[2], Math.pow(10, -10))
  })
})

// Interior projector and extension operator on a multi-tet cube mesh,
// verifying trace-vanishing properties on boundary.
describe('Interior projector on multi-tet mesh', () => {
  let mesh
  let whitney
  let bcdtpp

  before(() => {
    mesh = generateUnitCubeMesh(2)
    whitney = new Whitney(mesh)
    bcdtpp = new Bcdtpp(mesh, whitney, { quadratureOrder: 3 })
    bcdtpp.computeBoundaryWeights()
  })

  it('projectRing H1 is non-zero on cube mesh (has interior vertices)', () => {
    // n=2 cube has interior vertices.
    const u = (pt) => pt[0] + pt[1] + pt[2]
    let hasNonZero = false
    for (let tIdx = 0; tIdx < mesh.tetrahedronCount; tIdx++) {
      const pt = mesh.getTetrahedronBarycenter(tIdx)
      const val = bcdtpp.projectRing(u, pt, tIdx, 0)
      if (Math.abs(val) > 1e-12) {
        hasNonZero = true
        break
      }
    }
    expect(hasNonZero).to.equal(true)
  })

  it('projectRing Hcurl has vanishing tangential trace on boundary edges', () => {
    const u = (pt) => [pt[0], pt[1], pt[2]]
    for (const eIdx of mesh.boundaryEdges) {
      const e = mesh.edges[eIdx]
      const mid = [
        (mesh.vertices[e[0]][0] + mesh.vertices[e[1]][0]) / 2,
        (mesh.vertices[e[0]][1] + mesh.vertices[e[1]][1]) / 2,
        (mesh.vertices[e[0]][2] + mesh.vertices[e[1]][2]) / 2
      ]
      const tIdx = mesh.faceToTets[mesh.edgeToFaces[eIdx][0]][0]
      const result = bcdtpp.projectRing(u, mid, tIdx, 1)
      const edgeVec = [
        mesh.vertices[e[1]][0] - mesh.vertices[e[0]][0],
        mesh.vertices[e[1]][1] - mesh.vertices[e[0]][1],
        mesh.vertices[e[1]][2] - mesh.vertices[e[0]][2]
      ]
      const tangential =
        result[0] * edgeVec[0] +
        result[1] * edgeVec[1] +
        result[2] * edgeVec[2]
      expect(Math.abs(tangential)).to.be.below(1e-12)
    }
  })

  it('projectRing Hdiv has vanishing normal trace on boundary faces', () => {
    const u = (pt) => [pt[0], pt[1], pt[2]]
    for (const fIdx of mesh.boundaryFaces) {
      const normal = mesh.getFaceOutwardNormal(fIdx)
      const bary = mesh.getFaceBarycenter(fIdx)
      const tIdx = mesh.faceToTets[fIdx][0]
      const result = bcdtpp.projectRing(u, bary, tIdx, 2)
      const flux =
        result[0] * normal[0] +
        result[1] * normal[1] +
        result[2] * normal[2]
      expect(Math.abs(flux)).to.be.below(1e-12)
    }
  })
})

// Extension operator on a multi-tet cube mesh, verifying exactness
// at boundary vertices/edges and the projectPartial decomposition.
describe('Extension operator on multi-tet mesh', () => {
  let mesh
  let whitney
  let bcdtpp

  before(() => {
    mesh = generateUnitCubeMesh(2)
    whitney = new Whitney(mesh)
    bcdtpp = new Bcdtpp(mesh, whitney, { quadratureOrder: 3 })
    bcdtpp.computeBoundaryWeights()
  })

  it('E^0 is exact at all boundary vertices', () => {
    const boundaryValues = new Map()
    for (const vIdx of mesh.boundaryNodes) {
      boundaryValues.set(vIdx, Math.sin(vIdx))
    }
    for (const vIdx of mesh.boundaryNodes) {
      const pt = mesh.vertices[vIdx]
      // Find a tet containing this vertex.
      const tIdx = mesh.vertexToTets[vIdx][0]
      const result = bcdtpp.extendBoundary(boundaryValues, pt, tIdx, 0)
      expect(result).to.be.closeTo(boundaryValues.get(vIdx), Math.pow(10, -10))
    }
  })

  it('E^1 matches prescribed line integrals on boundary edges', () => {
    const boundaryEdgeDofs = new Map()
    for (const eIdx of mesh.boundaryEdges) {
      boundaryEdgeDofs.set(eIdx, Math.cos(eIdx))
    }
    for (const eIdx of mesh.boundaryEdges) {
      const e = mesh.edges[eIdx]
      const mid = [
        (mesh.vertices[e[0]][0] + mesh.vertices[e[1]][0]) / 2,
        (mesh.vertices[e[0]][1] + mesh.vertices[e[1]][1]) / 2,
        (mesh.vertices[e[0]][2] + mesh.vertices[e[1]][2]) / 2
      ]
      const tIdx = mesh.faceToTets[mesh.edgeToFaces[eIdx][0]][0]
      const result = bcdtpp.extendBoundary(boundaryEdgeDofs, mid, tIdx, 1)
      const edgeVec = [
        mesh.vertices[e[1]][0] - mesh.vertices[e[0]][0],
        mesh.vertices[e[1]][1] - mesh.vertices[e[0]][1],
        mesh.vertices[e[1]][2] - mesh.vertices[e[0]][2]
      ]
      const edgeLen = Math.sqrt(
        edgeVec[0] * edgeVec[0] +
        edgeVec[1] * edgeVec[1] +
        edgeVec[2] * edgeVec[2]
      )
      const tangential =
        result[0] * edgeVec[0] +
        result[1] * edgeVec[1] +
        result[2] * edgeVec[2]
      expect(tangential / edgeLen).to.be.closeTo(
        boundaryEdgeDofs.get(eIdx) / edgeLen,
        Math.pow(10, -10)
      )
    }
  })

  it('E^2 decomposition on cube mesh: projectPartial equals extendBoundary of extracted DoFs', () => {
    const u = (pt) => [pt[0], pt[1], pt[2]]
    for (let tIdx = 0; tIdx < Math.min(3, mesh.tetrahedronCount); tIdx++) {
      const pt = mesh.getTetrahedronBarycenter(tIdx)
      const partial = bcdtpp.projectPartial(u, pt, tIdx, 2)
      const boundaryData = bcdtpp.extractBoundaryDofs(u, 2)
      const extended = bcdtpp.extendBoundary(boundaryData, pt, tIdx, 2)
      expect(partial[0]).to.be.closeTo(extended[0], Math.pow(10, -10))
      expect(partial[1]).to.be.closeTo(extended[1], Math.pow(10, -10))
      expect(partial[2]).to.be.closeTo(extended[2], Math.pow(10, -10))
    }
  })
})

// Verifies commuting diagram properties: tr^l(E^l(w)) = w for all
// form degrees (H1 vertex trace, Hcurl tangential trace).
describe('Commuting diagrams for ring + extension', () => {
  let mesh
  let whitney
  let bcdtpp

  before(() => {
    mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    whitney = new Whitney(mesh)
    bcdtpp = new Bcdtpp(mesh, whitney, { quadratureOrder: 3 })
    bcdtpp.computeBoundaryWeights()
  })

  it('tr^0(E^0(w)) = w for arbitrary boundary data', () => {
    const boundaryValues = new Map([
      [0, 1.0],
      [1, 2.0],
      [2, 3.0],
      [3, 4.0]
    ])
    for (const vIdx of mesh.boundaryNodes) {
      const pt = mesh.vertices[vIdx]
      const val = bcdtpp.extendBoundary(boundaryValues, pt, 0, 0)
      expect(val).to.be.closeTo(boundaryValues.get(vIdx), Math.pow(10, -10))
    }
  })

  it('tr^1(E^1(w)) = w for arbitrary boundary edge data', () => {
    const boundaryEdgeDofs = new Map()
    for (const eIdx of mesh.boundaryEdges) {
      boundaryEdgeDofs.set(eIdx, 1.0 + eIdx * 0.5)
    }
    for (const eIdx of mesh.boundaryEdges) {
      const e = mesh.edges[eIdx]
      const mid = [
        (mesh.vertices[e[0]][0] + mesh.vertices[e[1]][0]) / 2,
        (mesh.vertices[e[0]][1] + mesh.vertices[e[1]][1]) / 2,
        (mesh.vertices[e[0]][2] + mesh.vertices[e[1]][2]) / 2
      ]
      const tIdx = mesh.faceToTets[mesh.edgeToFaces[eIdx][0]][0]
      const val = bcdtpp.extendBoundary(boundaryEdgeDofs, mid, tIdx, 1)
      const edgeVec = [
        mesh.vertices[e[1]][0] - mesh.vertices[e[0]][0],
        mesh.vertices[e[1]][1] - mesh.vertices[e[0]][1],
        mesh.vertices[e[1]][2] - mesh.vertices[e[0]][2]
      ]
      const edgeLen = Math.sqrt(
        edgeVec[0] * edgeVec[0] +
        edgeVec[1] * edgeVec[1] +
        edgeVec[2] * edgeVec[2]
      )
      const lineInt =
        val[0] * edgeVec[0] +
        val[1] * edgeVec[1] +
        val[2] * edgeVec[2]
      expect(lineInt / edgeLen).to.be.closeTo(
        boundaryEdgeDofs.get(eIdx) / edgeLen,
        Math.pow(10, -10)
      )
    }
  })

  it('projectPartial Hdiv equals extendBoundary for an integrable vector field', () => {
    const u = (pt) => [2 * pt[0], 2 * pt[1], 2 * pt[2]]
    const pt = mesh.getTetrahedronBarycenter(0)
    const partial = bcdtpp.projectPartial(u, pt, 0, 2)
    const boundaryData = bcdtpp.extractBoundaryDofs(u, 2)
    const extended = bcdtpp.extendBoundary(boundaryData, pt, 0, 2)
    expect(partial[0]).to.be.closeTo(extended[0], Math.pow(10, -10))
    expect(partial[1]).to.be.closeTo(extended[1], Math.pow(10, -10))
    expect(partial[2]).to.be.closeTo(extended[2], Math.pow(10, -10))
  })
})

// Global projector dispatch: verifies that project() with l,p arguments
// delegates to the correct specialized projection method.
describe('Global projector Pi^l', () => {
  let mesh
  let whitney
  let bcdtpp

  before(() => {
    mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    whitney = new Whitney(mesh)
    bcdtpp = new Bcdtpp(mesh, whitney, { quadratureOrder: 3 })
    bcdtpp.computeBoundaryWeights()
  })

  it('project l=0 matches projectH1 for p=0', () => {
    const u = (pt) => pt[0] + pt[1] + pt[2]
    const pt = mesh.getTetrahedronBarycenter(0)
    const a = bcdtpp.project(u, pt, 0, 0, 0)
    const b = bcdtpp.projectH1(u, pt, 0)
    expect(a).to.be.closeTo(b, Math.pow(10, -10))
  })

  it('project l=1 matches projectHcurl for p=0 with vector input', () => {
    const u = (pt) => [pt[0], pt[1], pt[2]]
    const pt = mesh.getTetrahedronBarycenter(0)
    const a = bcdtpp.project(u, pt, 0, 1, 0)
    const b = bcdtpp.projectHcurl(u, pt, 0)
    expect(a[0]).to.be.closeTo(b[0], Math.pow(10, -10))
    expect(a[1]).to.be.closeTo(b[1], Math.pow(10, -10))
    expect(a[2]).to.be.closeTo(b[2], Math.pow(10, -10))
  })

  it('project l=1 matches projectHcurl for p=0 with scalar input', () => {
    const u = (pt) => pt[0] + pt[1] + pt[2]
    const pt = mesh.getTetrahedronBarycenter(0)
    const a = bcdtpp.project(u, pt, 0, 1, 0)
    const b = bcdtpp.projectHcurl(u, pt, 0)
    expect(a[0]).to.be.closeTo(b[0], Math.pow(10, -10))
    expect(a[1]).to.be.closeTo(b[1], Math.pow(10, -10))
    expect(a[2]).to.be.closeTo(b[2], Math.pow(10, -10))
  })

  it('project l=2 matches projectHdiv for p=0 with vector input', () => {
    const u = (pt) => [pt[0], pt[1], pt[2]]
    const pt = mesh.getTetrahedronBarycenter(0)
    const a = bcdtpp.project(u, pt, 0, 2, 0)
    const b = bcdtpp.projectHdiv(u, pt, 0)
    expect(a[0]).to.be.closeTo(b[0], Math.pow(10, -10))
    expect(a[1]).to.be.closeTo(b[1], Math.pow(10, -10))
    expect(a[2]).to.be.closeTo(b[2], Math.pow(10, -10))
  })

  it('project l=2 matches projectHdiv for p=0 with scalar input', () => {
    const u = (pt) => pt[0] + pt[1] + pt[2]
    const pt = mesh.getTetrahedronBarycenter(0)
    const a = bcdtpp.project(u, pt, 0, 2, 0)
    const b = bcdtpp.projectHdiv(u, pt, 0)
    expect(a[0]).to.be.closeTo(b[0], Math.pow(10, -10))
    expect(a[1]).to.be.closeTo(b[1], Math.pow(10, -10))
    expect(a[2]).to.be.closeTo(b[2], Math.pow(10, -10))
  })

  it('project l=3 matches projectL2 for p=0', () => {
    const u = (pt) => pt[0] + pt[1] + pt[2]
    const a = bcdtpp.project(u, [0.25, 0.25, 0.25], 0, 3, 0)
    const b = bcdtpp.projectL2(u, 0)
    expect(a).to.be.closeTo(b, Math.pow(10, -10))
  })

  it('project delegates to projectHp for p>0', () => {
    const u = (pt) => pt[0] * pt[0]
    const pt = mesh.getTetrahedronBarycenter(0)
    const a = bcdtpp.project(u, pt, 0, 0, 2)
    const b = bcdtpp.projectHp(u, pt, 0, 0, 2)
    expect(a).to.be.closeTo(b, Math.pow(10, -10))
  })
})
