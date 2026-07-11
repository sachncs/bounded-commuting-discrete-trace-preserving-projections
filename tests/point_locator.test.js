/**
 * Tests for the PointLocator AABB tree: interior/vertex/boundary queries,
 * multi-tet meshes, tree partitioning on larger meshes, and degenerate
 * tet rejection at construction time.
 */
import { expect } from 'chai'
import { PointLocator } from '../src/lib/point_locator.js'
import { Mesh } from '../src/lib/mesh.js'
import { MeshValidationError } from '../src/lib/errors.js'
import { generateUnitCubeMesh } from '../src/lib/mesh_generator.js'

// Tests the AABB tree point-in-tet locator on single-tet, two-tet,
// and multi-tet (cube) meshes, plus tree partitioning behavior.
describe('PointLocator', () => {
  const singleTet = {
    vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
    tetrahedra: [[0, 1, 2, 3]]
  }

  it('finds point inside tetrahedron', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    const locator = new PointLocator(mesh)
    const result = locator.findTetrahedron([0.1, 0.1, 0.1])
    expect(result).to.not.equal(null)
    expect(result.tIdx).to.equal(0)
  })

  it('returns null for point outside', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    const locator = new PointLocator(mesh)
    const result = locator.findTetrahedron([1, 1, 1])
    expect(result).to.equal(null)
  })

  it('finds point on vertex', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    const locator = new PointLocator(mesh)
    const result = locator.findTetrahedron([0, 0, 0])
    expect(result).to.not.equal(null)
    expect(result.tIdx).to.equal(0)
  })

  it('handles two-tetrahedron mesh', () => {
    const twoTet = {
      vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]],
      tetrahedra: [[0, 1, 2, 3], [1, 2, 3, 4]]
    }
    const mesh = new Mesh(twoTet.vertices, twoTet.tetrahedra)
    const locator = new PointLocator(mesh)
    const r1 = locator.findTetrahedron([0.1, 0.1, 0.1])
    expect(r1).to.not.equal(null)
    const r2 = locator.findTetrahedron([0.8, 0.8, 0.8])
    expect(r2).to.not.equal(null)
  })

  it('partitions a mesh with more than 8 tets', () => {
    const mesh = generateUnitCubeMesh(2)
    const locator = new PointLocator(mesh)
    // Cube [0,1]^3 center should be inside some tet.
    const result = locator.findTetrahedron([0.5, 0.5, 0.5])
    expect(result).to.not.equal(null)
    expect(Number.isInteger(result.tIdx)).to.equal(true)
  })

  it('rejects degenerate tet at mesh construction', () => {
    expect(() => new Mesh(
      [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]],
      [[0, 1, 2, 3]]
    )).to.throw(MeshValidationError)
  })

  it('searches right subtree when point not in left', () => {
    const mesh = generateUnitCubeMesh(2)
    const locator = new PointLocator(mesh)
    // Query a point near a corner to exercise both subtrees.
    const result = locator.findTetrahedron([0.9, 0.9, 0.9])
    expect(result).to.not.equal(null)
  })
})
