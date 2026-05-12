import { expect } from 'chai'
import { Whitney } from '../src/lib/whitney.js'
import { Mesh } from '../src/lib/mesh.js'

describe('Whitney', () => {
  const singleTet = {
    vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
    tetrahedra: [[0, 1, 2, 3]]
  }

  it('computes correct barycentric coordinates', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    const w = new Whitney(mesh)
    const bary = w.getBarycentric(0, [0, 0, 0])
    expect(bary[0]).to.be.closeTo(1, Math.pow(10, -10))
    expect(bary[1]).to.be.closeTo(0, Math.pow(10, -10))
    expect(bary[2]).to.be.closeTo(0, Math.pow(10, -10))
    expect(bary[3]).to.be.closeTo(0, Math.pow(10, -10))
  })

  it('barycentric coordinates sum to 1', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    const w = new Whitney(mesh)
    const bary = w.getBarycentric(0, [0.1, 0.1, 0.1])
    const sum = bary.reduce((a, b) => a + b, 0)
    expect(sum).to.be.closeTo(1, Math.pow(10, -10))
  })

  it('gradients of barycentric coords sum to zero', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    const w = new Whitney(mesh)
    const grads = w.getGradBarycentric(0)
    const sum = [
      grads[0][0] + grads[1][0] + grads[2][0] + grads[3][0],
      grads[0][1] + grads[1][1] + grads[2][1] + grads[3][1],
      grads[0][2] + grads[1][2] + grads[2][2] + grads[3][2]
    ]
    expect(sum[0]).to.be.closeTo(0, Math.pow(10, -10))
    expect(sum[1]).to.be.closeTo(0, Math.pow(10, -10))
    expect(sum[2]).to.be.closeTo(0, Math.pow(10, -10))
  })

  it('edge basis has 6 components', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    const w = new Whitney(mesh)
    const bary = [0.25, 0.25, 0.25, 0.25]
    const edgeBasis = w.getEdgeBasis(0, bary)
    expect(edgeBasis.length).to.equal(6)
    edgeBasis.forEach((phi) => {
      expect(phi.length).to.equal(3)
    })
  })

  it('face basis has 4 components', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra)
    const w = new Whitney(mesh)
    const bary = [0.25, 0.25, 0.25, 0.25]
    const faceBasis = w.getFaceBasis(0, bary)
    expect(faceBasis.length).to.equal(4)
    faceBasis.forEach((psi) => {
      expect(psi.length).to.equal(3)
    })
  })

  it('rejects degenerate tetrahedron at mesh level', () => {
    expect(() => new Mesh(
      [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]],
      [[0, 1, 2, 3]]
    )).to.throw(/degenerate/)
  })
})
