import { expect } from 'chai'
import sinon from 'sinon'
import { LocalSolver } from '../src/lib/local_solver.js'
import { SingularMatrixError } from '../src/lib/errors.js'

describe('LocalSolver', () => {
  afterEach(() => sinon.restore())

  it('assembleSurfaceStiffness returns symmetric matrix', () => {
    const vertices = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
    const triangles = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
    const K = LocalSolver.assembleSurfaceStiffness(vertices, triangles)
    const n = K.length
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        expect(K[i][j]).to.be.closeTo(K[j][i], Math.pow(10, -10))
      }
    }
  })

  it('solveWithConstraint warns for ill-conditioned matrix', () => {
    const warnSpy = sinon.spy()
    // Construct an ill-conditioned matrix by scaling one row/column very large.
    const K = [
      [1e14, 0, 0],
      [0, 1, 0],
      [0, 0, 1]
    ]
    const b = [0, 0, 0]
    LocalSolver.solveWithConstraint(K, b, warnSpy)
    expect(warnSpy.called).to.equal(true)
    expect(warnSpy.getCall(0).args[0].code).to.equal('LOCAL_SOLVER_ILL_CONDITIONED')
    expect(warnSpy.getCall(0).args[0].message).to.match(/ill-conditioned/)
  })

  it('solveWithConstraint throws SingularMatrixError for singular matrix', () => {
    const K = [
      [1, 1, 1],
      [1, 1, 1],
      [1, 1, 1]
    ]
    const b = [0, 0, 0]
    expect(() => LocalSolver.solveWithConstraint(K, b)).to.throw(SingularMatrixError)
  })

  it('assembleSurfaceStiffness throws for degenerate triangle', () => {
    const vertices = [[0, 0, 0], [1, 0, 0], [2, 0, 0]]
    const triangles = [[0, 1, 2]]
    expect(() => LocalSolver.assembleSurfaceStiffness(vertices, triangles)).to.throw(
      SingularMatrixError
    )
  })

  it('solveWithConstraint returns empty array for n=0', () => {
    const result = LocalSolver.solveWithConstraint([], [])
    expect(result).to.deep.equal([])
  })
})
