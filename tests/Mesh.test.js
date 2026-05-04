import {describe, it, expect} from 'vitest';
import {Mesh} from '../src/lib/mesh.js';

describe('Mesh', () => {
  const singleTet = {
    vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
    tetrahedra: [[0, 1, 2, 3]],
  };

  it('builds topology for a single tetrahedron', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra);
    expect(mesh.tetrahedronCount).toBe(1);
    expect(mesh.faces.length).toBe(4);
    expect(mesh.edges.length).toBe(6);
    expect(mesh.boundaryFaces.length).toBe(4);
    expect(mesh.boundaryEdges.length).toBe(6);
    expect(mesh.boundaryNodes.size).toBe(4);
  });

  it('computes correct volume', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra);
    const vol = mesh.getVolume(0);
    expect(vol).toBeCloseTo(1 / 6, 10);
  });

  it('detects degenerate tetrahedra', () => {
    const badMesh = {
      vertices: [[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]],
      tetrahedra: [[0, 1, 2, 3]],
    };
    const mesh = new Mesh(badMesh.vertices, badMesh.tetrahedra);
    expect(mesh.getVolume(0)).toBeLessThan(1e-12);
  });

  it('computes Alfeld split', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra);
    mesh.computeAlfeldSplit();
    expect(mesh.alfeldTriangles.length).toBe(4);
    const totalSubTris = mesh.alfeldTriangles.reduce(
      (s, at) => s + at.triangles.length,
      0,
    );
    expect(totalSubTris).toBe(12);
  });

  it('computes Worsey-Farin split', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra);
    mesh.computeWorseyFarinSplit();
    expect(mesh.worseyFarinTetrahedra.length).toBe(1);
    expect(mesh.worseyFarinTetrahedra[0].tetrahedra.length).toBe(12);
  });

  it('getBoundaryStar returns correct faces for a vertex', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra);
    const star = mesh.getBoundaryStar(0);
    expect(star.length).toBe(3);
  });

  it('getStar returns all tets containing a vertex', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra);
    expect(mesh.getStar(0)).toEqual([0]);
  });

  it('handles two-tetrahedron mesh', () => {
    const twoTet = {
      vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]],
      tetrahedra: [[0, 1, 2, 3], [1, 2, 3, 4]],
    };
    const mesh = new Mesh(twoTet.vertices, twoTet.tetrahedra);
    expect(mesh.tetrahedronCount).toBe(2);
    expect(mesh.faces.length).toBe(7);
    expect(mesh.boundaryFaces.length).toBe(6);
  });
});
