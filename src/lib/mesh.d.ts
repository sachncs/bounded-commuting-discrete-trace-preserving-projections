/**
 * Tetrahedral mesh data structure with topological connectivity and boundary
 * extraction.
 */
export class Mesh {
  /** Vertex coordinates. */
  vertices: number[][];
  /** Tetrahedron vertex indices. */
  tetrahedra: number[][];
  /** Number of vertices in the original (unrefined) mesh. */
  readonly originalVertexCount: number;
  /** Current total number of vertices (includes refinements). */
  vertexCount: number;
  /** Number of tetrahedra. */
  readonly tetrahedronCount: number;

  /** Global faces as vertex index triples. */
  faces: number[][];
  /** Global edges as vertex index pairs. */
  edges: number[][];
  /** Global indices of boundary faces. */
  boundaryFaces: number[];
  /** Global indices of boundary edges. */
  boundaryEdges: number[];
  /** Global vertex indices on the boundary. */
  boundaryNodes: Set<number>;

  /** Maps each face index to the tetrahedra that contain it. */
  faceToTets: number[][];
  /** Maps each edge index to the faces that contain it. */
  edgeToFaces: number[][];

  /** Maps each vertex to the tetrahedra that contain it. */
  vertexToTets: number[][];
  /** Maps each vertex to the boundary faces that contain it. */
  vertexToBoundaryFaces: number[][];
  /** Maps each vertex to the edges that contain it. */
  vertexToEdges: number[][];

  /**
   * @param vertices - Vertex coordinates.
   * @param tetrahedra - Tetrahedron vertex indices.
   */
  constructor(vertices: number[][], tetrahedra: number[][]);

  /** Returns the geometric area of a face. */
  getFaceArea(fIdx: number): number;

  /** Looks up the global edge index for an edge key. */
  getEdgeIndex(edgeKey: number): number | undefined;

  /** Looks up the global face index for a face key. */
  getFaceIndex(faceKey: number): number | undefined;

  /** Returns the edge orientation sign for a local edge within a tetrahedron. */
  getTetEdgeSign(tIdx: number, e: number): number;

  /** Returns the face orientation sign for a local face within a tetrahedron. */
  getTetFaceSign(tIdx: number, f: number): number;

  /** Returns all vertex coordinates. */
  getVertices(): number[][];

  /** Returns all tetrahedra as vertex index arrays. */
  getTetrahedra(): number[][];

  /** Returns all faces as vertex index triples. */
  getFaces(): number[][];

  /** Returns all edges as vertex index pairs. */
  getEdges(): number[][];

  /** Returns the set of boundary node indices. */
  getBoundaryNodes(): Set<number>;

  /** Returns the array of boundary face indices. */
  getBoundaryFaces(): number[];

  /** Returns the array of boundary edge indices. */
  getBoundaryEdges(): number[];

  /** Returns the original vertex count before refinement. */
  getOriginalVertexCount(): number;

  /** Returns the barycenter of a face. */
  getFaceBarycenter(fIdx: number): number[];

  /** Returns the barycenter of a tetrahedron. */
  getTetrahedronBarycenter(tIdx: number): number[];

  /** Returns the star (incident tets) of a vertex. */
  getStar(vIdx: number): number[];

  /** Returns the boundary faces incident to a vertex. */
  getBoundaryStar(vIdx: number): number[];

  /** Returns the edges incident to a vertex. */
  getEdgeStar(vIdx: number): number[];

  /** Returns the volume of a tetrahedron. */
  getVolume(tIdx: number): number;

  /** Returns the global face indices of a tetrahedron's four faces. */
  getTetrahedronFaces(tIdx: number): number[];

  /**
   * Computes the outward-pointing unit normal for a boundary face.
   * For interior faces the orientation is arbitrary.
   */
  getFaceOutwardNormal(fIdx: number): number[];
}
