import {describe, it, expect} from 'vitest';
import {validateMeshData} from '../src/ui.js';

describe('validateMeshData', () => {
  it('accepts valid mesh data', () => {
    const data = {
      vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
      tetrahedra: [[0, 1, 2, 3]],
    };
    expect(() => validateMeshData(data)).not.toThrow();
  });

  it('rejects missing vertices', () => {
    expect(() => validateMeshData({tetrahedra: [[0, 1, 2, 3]]})).toThrow('vertices');
  });

  it('rejects missing tetrahedra', () => {
    expect(() => validateMeshData({vertices: [[0, 0, 0]]})).toThrow('tetrahedra');
  });

  it('rejects bad vertex shape', () => {
    const data = {
      vertices: [[0, 0]],
      tetrahedra: [[0, 1, 2, 3]],
    };
    expect(() => validateMeshData(data)).toThrow('Vertex 0');
  });

  it('rejects out-of-bounds indices', () => {
    const data = {
      vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
      tetrahedra: [[0, 1, 2, 4]],
    };
    expect(() => validateMeshData(data)).toThrow('out-of-bounds');
  });

  it('rejects non-integer tetrahedron indices', () => {
    const data = {
      vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
      tetrahedra: [[0, 1, 2, 3.5]],
    };
    expect(() => validateMeshData(data)).toThrow('integer');
  });
});
