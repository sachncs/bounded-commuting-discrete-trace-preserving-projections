import type {Mesh} from './mesh.js';

/**
 * Mesh refinement operators: Alfeld face splitting and Worsey-Farin tetrahedron
 * splitting.
 *
 * This class mutates the underlying Mesh by appending barycenter vertices.
 * It stores the refinement data (sub-triangles, sub-tetrahedra) separately so
 * that Mesh remains a pure data structure.
 */
export class MeshRefinement {
  constructor(mesh: Mesh);

  /** Face barycenter vertex indices (mutated during split). */
  faceBarycenters: number[];
  /** Tetrahedron barycenter vertex indices (mutated during split). */
  tetBarycenters: number[];
  /** Alfeld-split triangles for each boundary face. */
  alfeldTriangles: Array<{
    parentFaceIdx: number;
    triangles: number[][];
  }>;
  /** Maps parent face index to its Alfeld split data. */
  faceToAlfeld: Map<number, {parentFaceIdx: number; triangles: number[][]}>;
  /** Worsey-Farin sub-tetrahedra for each parent tet. */
  worseyFarinTetrahedra: Array<{
    parentTetIdx: number;
    tetrahedra: number[][];
  }>;

  /**
   * Section 6.1.3: Alfeld split of boundary faces.
   *
   * Mutates mesh.vertices by appending face barycenters.
   */
  computeAlfeldSplit(): void;

  /**
   * Section 6.1.4: Worsey-Farin split of tetrahedra.
   *
   * Mutates mesh.vertices by appending face and tetrahedron barycenters.
   * Calls computeAlfeldSplit() if not already done.
   */
  computeWorseyFarinSplit(): void;
}
