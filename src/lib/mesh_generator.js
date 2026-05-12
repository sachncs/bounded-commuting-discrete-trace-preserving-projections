/**
 * Simple structured tetrahedral mesh generators for convergence experiments.
 */

import { Mesh } from './mesh.js'
import { MeshValidationError } from './errors.js'

import { tetDeterminant } from './math_utils.js'

/**
 * Generates a uniform tetrahedral mesh of the unit cube [0,1]^3 using the
 * Freudenthal (Kuhn) triangulation: each cube is split into 6 tets along
 * the body diagonal from (0,0,0) to (1,1,1).
 *
 * @param {number} n - Number of cubes per axis (creates n^3 cubes).
 * @return {!Mesh}
 */
export function generateUnitCubeMesh (n) {
  if (!Number.isInteger(n) || n < 1) {
    throw new MeshValidationError('n must be a positive integer')
  }

  const vertices = []
  for (let k = 0; k <= n; k++) {
    for (let j = 0; j <= n; j++) {
      for (let i = 0; i <= n; i++) {
        vertices.push([i / n, j / n, k / n])
      }
    }
  }

  const vertexIndex = (i, j, k) => k * (n + 1) * (n + 1) + j * (n + 1) + i

  // Freudenthal chains from (0,0,0) to (1,1,1): each path flips one bit at a time.
  // Local cube vertices: 0=(0,0,0), 1=(1,0,0), 2=(0,1,0), 3=(1,1,0),
  //                     4=(0,0,1), 5=(1,0,1), 6=(0,1,1), 7=(1,1,1)
  const chains = [
    [0, 1, 3, 7],
    [0, 1, 5, 7],
    [0, 2, 3, 7],
    [0, 2, 6, 7],
    [0, 4, 5, 7],
    [0, 4, 6, 7]
  ]

  const tetrahedra = []
  for (let k = 0; k < n; k++) {
    for (let j = 0; j < n; j++) {
      for (let i = 0; i < n; i++) {
        const cubeVerts = [
          vertexIndex(i, j, k),
          vertexIndex(i + 1, j, k),
          vertexIndex(i, j + 1, k),
          vertexIndex(i + 1, j + 1, k),
          vertexIndex(i, j, k + 1),
          vertexIndex(i + 1, j, k + 1),
          vertexIndex(i, j + 1, k + 1),
          vertexIndex(i + 1, j + 1, k + 1)
        ]

        for (const chain of chains) {
          const tet = chain.map((lv) => cubeVerts[lv])
          const verts = tet.map((vi) => vertices[vi])
          const vol = tetDeterminant(verts[0], verts[1], verts[2], verts[3])
          if (vol < 0) {
            // Reverse last two vertices to fix orientation.
            [tet[2], tet[3]] = [tet[3], tet[2]]
          }
          tetrahedra.push(tet)
        }
      }
    }
  }

  return new Mesh(vertices, tetrahedra)
}

/**
 * Generates a single reference tetrahedron mesh.
 * @return {!Mesh}
 */
export function generateSingleTetMesh () {
  return new Mesh(
    [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
    [[0, 1, 2, 3]]
  )
}
