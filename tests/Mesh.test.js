/**
 * Tests for the Mesh and MeshRefinement classes: topology construction
 * (faces, edges, boundary), geometry (volume, normals), orientation
 * signs, vertex stars, and Alfeld/Worsey-Farin split idempotency.
 */
import { expect } from 'chai'
import { Mesh } from '../src/lib/mesh.js'
import { MeshRefinement } from '../src/lib/mesh_refinement.js'

// Verifies Mesh topology construction: faces, edges, boundary sets,
// volume, orientation signs, vertex stars, and multi-tet connectivity.
describe('Mesh', () => {
  const singleTet = {
    vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
    tetrahedra: [[0, 1, 2, 3]]
  }

  it('builds topology for a single tetrahedron', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    expect(mesh.tetrahedronCount).to.equal(1)
    expect(mesh.faces.length).to.equal(4)
    expect(mesh.edges.length).to.equal(6)
    expect(mesh.boundaryFaces.length).to.equal(4)
    expect(mesh.boundaryEdges.length).to.equal(6)
    expect(mesh.boundaryNodes.size).to.equal(4)
  })

  it('computes correct volume', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    const vol = mesh.getVolume(0)
    expect(vol).to.be.closeTo(1 / 6, Math.pow(10, -10))
  })

  it('rejects degenerate tetrahedra', () => {
    expect(() => new Mesh(
      [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]],
      [[0, 1, 2, 3]]
    )).to.throw(/degenerate/)
  })

  it('getBoundaryStar returns correct faces for a vertex', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    const star = mesh.getBoundaryStar(0)
    expect(star.length).to.equal(3)
  })

  it('getStar returns all tets containing a vertex', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    expect(mesh.getStar(0)).to.deep.equal([0])
  })

  it('handles two-tetrahedron mesh', () => {
    const twoTet = {
      vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]],
      tetrahedra: [[0, 1, 2, 3], [1, 2, 3, 4]]
    }
    const mesh = new Mesh(twoTet.vertices, twoTet.tetrahedra)
    expect(mesh.tetrahedronCount).to.equal(2)
    expect(mesh.faces.length).to.equal(7)
    expect(mesh.boundaryFaces.length).to.equal(6)
  })

  it('getTetFaceSign returns +1 or -1', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    const sign = mesh.getTetFaceSign(0, 0)
    expect(sign === 1 || sign === -1).to.equal(true)
  })

  it('getEdgeStar returns edges for a boundary vertex', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    const star = mesh.getEdgeStar(0)
    expect(star.length).to.be.above(0)
    expect(star.every((e) => Number.isInteger(e))).to.equal(true)
  })
})

// Verifies Alfeld and Worsey-Farin mesh refinement: split counts,
// sub-triangle/tet counts, and idempotency of both splits.
describe('MeshRefinement', () => {
  const singleTet = {
    vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
    tetrahedra: [[0, 1, 2, 3]]
  }

  it('computes Alfeld split', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    const refinement = new MeshRefinement(mesh)
    refinement.computeAlfeldSplit()
    expect(refinement.alfeldTriangles.length).to.equal(4)
    const totalSubTris = refinement.alfeldTriangles.reduce(
      (s, at) => s + at.triangles.length,
      0
    )
    expect(totalSubTris).to.equal(12)
  })

  it('computes Worsey-Farin split', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    const refinement = new MeshRefinement(mesh)
    refinement.computeWorseyFarinSplit()
    expect(refinement.worseyFarinTetrahedra.length).to.equal(1)
    expect(refinement.worseyFarinTetrahedra[0].tetrahedra.length).to.equal(12)
  })

  it('Alfeld split is idempotent', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    const refinement = new MeshRefinement(mesh)
    refinement.computeAlfeldSplit()
    const vCountAfterFirst = mesh.vertexCount
    refinement.computeAlfeldSplit()
    expect(mesh.vertexCount).to.equal(vCountAfterFirst)
  })

  it('Worsey-Farin split is idempotent', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    const refinement = new MeshRefinement(mesh)
    refinement.computeWorseyFarinSplit()
    const vCountAfterFirst = mesh.vertexCount
    refinement.computeWorseyFarinSplit()
    expect(mesh.vertexCount).to.equal(vCountAfterFirst)
  })
})
