import type {Mesh} from './mesh.js';

/**
 * Axis-aligned bounding box (AABB) tree for fast point-in-tetrahedron queries.
 *
 * Recursively partitions the mesh bounding box to achieve O(log N) search.
 */
export class PointLocator {
  /**
   * @param mesh - The mesh to index.
   * @param maxLeafSize - Maximum tets per leaf node (default 8).
   */
  constructor(mesh: Mesh, maxLeafSize?: number);

  /**
   * Finds the tetrahedron containing the given point.
   * @param point - Cartesian point [x, y, z].
   * @returns Object with tetrahedron index and barycentric coordinates, or null if not found.
   */
  findTetrahedron(point: number[]): {
    tIdx: number;
    bary: number[];
  } | null;
}
