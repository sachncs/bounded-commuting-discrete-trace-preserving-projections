import {describe, it, expect} from 'vitest';
import {PointLocator} from '../src/lib/point_locator.js';
import {Mesh} from '../src/lib/mesh.js';

describe('PointLocator', () => {
  const singleTet = {
    vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
    tetrahedra: [[0, 1, 2, 3]],
  };

  it('finds point inside tetrahedron', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra);
    const locator = new PointLocator(mesh);
    const result = locator.findTetrahedron([0.1, 0.1, 0.1]);
    expect(result).not.toBeNull();
    expect(result.tIdx).toBe(0);
  });

  it('returns null for point outside', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra);
    const locator = new PointLocator(mesh);
    const result = locator.findTetrahedron([1, 1, 1]);
    expect(result).toBeNull();
  });

  it('finds point on vertex', () => {
    const mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra);
    const locator = new PointLocator(mesh);
    const result = locator.findTetrahedron([0, 0, 0]);
    expect(result).not.toBeNull();
    expect(result.tIdx).toBe(0);
  });

  it('handles two-tetrahedron mesh', () => {
    const twoTet = {
      vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]],
      tetrahedra: [[0, 1, 2, 3], [1, 2, 3, 4]],
    };
    const mesh = new Mesh(twoTet.vertices, twoTet.tetrahedra);
    const locator = new PointLocator(mesh);
    const r1 = locator.findTetrahedron([0.1, 0.1, 0.1]);
    expect(r1).not.toBeNull();
    const r2 = locator.findTetrahedron([0.8, 0.8, 0.8]);
    expect(r2).not.toBeNull();
  });
});
