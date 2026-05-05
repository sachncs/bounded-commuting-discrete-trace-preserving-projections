/**
 * Mesh refinement operators: Alfeld face splitting and Worsey-Farin tetrahedron
 * splitting.
 *
 * This class mutates the underlying Mesh by appending barycenter vertices.
 * It stores the refinement data (sub-triangles, sub-tetrahedra) separately so
 * that Mesh remains a pure data structure.
 */

export class MeshRefinement {
  /**
   * @param {!Mesh} mesh
   */
  constructor(mesh) {
    this.mesh = mesh;
    this.faceBarycenters = [];
    this.tetBarycenters = [];
    this.alfeldTriangles = [];
    /** @type {!Map<number, {parentFaceIdx: number, triangles: !Array<!Array<number>>}>} */
    this.faceToAlfeld = new Map();
    this.worseyFarinTetrahedra = [];
  }

  /**
   * Section 6.1.3: Alfeld split of boundary faces.
   *
   * Mutates mesh.vertices by appending face barycenters.
   */
  computeAlfeldSplit() {
    this.faceBarycenters = new Array(this.mesh.faces.length).fill(-1);
    this.alfeldTriangles = [];

    this.mesh.boundaryFaces.forEach((fIdx) => {
      const f = this.mesh.faces[fIdx];
      const bary = this.mesh.getFaceBarycenter(fIdx);
      const vBaryIdx = this.mesh.vertices.length;
      this.mesh.vertices.push(bary);
      this.faceBarycenters[fIdx] = vBaryIdx;

      const subTriangles = [
        [f[0], f[1], vBaryIdx],
        [f[1], f[2], vBaryIdx],
        [f[2], f[0], vBaryIdx],
      ];
      const entry = {parentFaceIdx: fIdx, triangles: subTriangles};
      this.alfeldTriangles.push(entry);
      this.faceToAlfeld.set(fIdx, entry);
    });
    this.mesh.vertexCount = this.mesh.vertices.length;
  }

  /**
   * Section 6.1.4: Worsey-Farin split of tetrahedra.
   *
   * Mutates mesh.vertices by appending face and tetrahedron barycenters.
   * Calls computeAlfeldSplit() if not already done.
   */
  computeWorseyFarinSplit() {
    if (this.faceBarycenters.length === 0) {
      this.computeAlfeldSplit();
    }

    this.tetBarycenters = new Array(this.mesh.tetrahedronCount).fill(-1);
    this.worseyFarinTetrahedra = [];

    for (let tIdx = 0; tIdx < this.mesh.tetrahedronCount; tIdx++) {
      const tet = this.mesh.tetrahedra[tIdx];
      const bary = this.mesh.getTetrahedronBarycenter(tIdx);
      const vTetBaryIdx = this.mesh.vertices.length;
      this.mesh.vertices.push(bary);
      this.tetBarycenters[tIdx] = vTetBaryIdx;

      const tetSubTets = [];
      const tFaces = this.mesh.getTetrahedronFaces(tIdx);

      tFaces.forEach((fIdx) => {
        const f = this.mesh.faces[fIdx];
        let fvBaryIdx = this.faceBarycenters[fIdx];
        if (fvBaryIdx === -1) {
          const fbary = this.mesh.getFaceBarycenter(fIdx);
          fvBaryIdx = this.mesh.vertices.length;
          this.mesh.vertices.push(fbary);
          this.faceBarycenters[fIdx] = fvBaryIdx;
        }

        const subTris = [
          [f[0], f[1], fvBaryIdx],
          [f[1], f[2], fvBaryIdx],
          [f[2], f[0], fvBaryIdx],
        ];

        subTris.forEach((tri) => {
          tetSubTets.push([...tri, vTetBaryIdx]);
        });
      });
      this.worseyFarinTetrahedra.push({
        parentTetIdx: tIdx,
        tetrahedra: tetSubTets,
      });
    }
    this.mesh.vertexCount = this.mesh.vertices.length;
  }
}
