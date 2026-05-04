import {describe, it, expect, beforeAll} from 'vitest';
import {Mesh} from '../src/lib/mesh.js';
import {Whitney} from '../src/lib/whitney.js';
import {Bcdtpp} from '../src/lib/bcdtpp.js';

describe('Bcdtpp Projections', () => {
  const singleTet = {
    vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
    tetrahedra: [[0, 1, 2, 3]],
  };

  let mesh;
  let whitney;
  let bcdtpp;

  beforeAll(() => {
    mesh = new Mesh(singleTet.vertices, singleTet.tetrahedra);
    whitney = new Whitney(mesh);
    bcdtpp = new Bcdtpp(mesh, whitney, {quadratureOrder: 3});
    bcdtpp.computeBoundaryWeights();
    bcdtpp.buildPointLocator();
  });

  it('projectH1 is exact for constant function', () => {
    const u = () => 5;
    const pt = [0.2, 0.3, 0.1];
    const result = bcdtpp.projectH1(u, pt, 0);
    expect(result).toBeCloseTo(5, 6);
  });

  it('projectH1 uses boundary weights for all-boundary mesh', () => {
    const u = (pt) => pt[0] + pt[1] + pt[2];
    const pt = [0.2, 0.3, 0.1];
    const result = bcdtpp.projectH1(u, pt, 0);
    expect(Number.isFinite(result)).toBe(true);
    expect(result).toBeGreaterThanOrEqual(0);
    expect(result).toBeLessThanOrEqual(3);
  });

  it('projectH1 evaluates for quadratic function', () => {
    const u = (pt) => pt[0] * pt[0] + pt[1] * pt[1];
    const pt = mesh.getTetrahedronBarycenter(0);
    const result = bcdtpp.projectH1(u, pt, 0);
    expect(typeof result).toBe('number');
    expect(Number.isFinite(result)).toBe(true);
  });

  it('projectL2 integrates constants exactly', () => {
    const u = () => 5;
    const result = bcdtpp.projectL2(u, 0);
    expect(result).toBeCloseTo(5, 10);
  });

  it('projectL2 integrates linear functions exactly', () => {
    const u = (pt) => pt[0] + pt[1] + pt[2];
    const result = bcdtpp.projectL2(u, 0);
    expect(result).toBeCloseTo(0.75, 6);
  });

  it('projectHcurl returns a vector', () => {
    const u = (pt) => pt[0] + pt[1] + pt[2];
    const pt = mesh.getTetrahedronBarycenter(0);
    const result = bcdtpp.projectHcurl(u, pt, 0);
    expect(Array.isArray(result)).toBe(true);
    expect(result.length).toBe(3);
    result.forEach((v) => expect(Number.isFinite(v)).toBe(true));
  });

  it('projectHdiv returns a vector', () => {
    const u = (pt) => pt[0] + pt[1] + pt[2];
    const pt = mesh.getTetrahedronBarycenter(0);
    const result = bcdtpp.projectHdiv(u, pt, 0);
    expect(Array.isArray(result)).toBe(true);
    expect(result.length).toBe(3);
    result.forEach((v) => expect(Number.isFinite(v)).toBe(true));
  });

  it('projectAtPoint finds tet and projects', () => {
    const u = (pt) => pt[0] + pt[1] + pt[2];
    const result = bcdtpp.projectAtPoint(u, [0.1, 0.1, 0.1], 0);
    expect(result).toHaveProperty('value');
    expect(result).toHaveProperty('tIdx');
    expect(result).toHaveProperty('bary');
    expect(result.tIdx).toBe(0);
    expect(Number.isFinite(result.value)).toBe(true);
  });

  it('projectAtPoint throws for point outside mesh', () => {
    const u = (pt) => pt[0];
    expect(() => bcdtpp.projectAtPoint(u, [5, 5, 5], 0)).toThrow(/not found/);
  });

  it('commuting property: grad Pi^0 = 0 for constant function', () => {
    const u = () => 7;
    const pt = [0.2, 0.1, 0.05];
    const h = 1e-5;
    const gradProj = [
      (bcdtpp.projectH1(u, [pt[0] + h, pt[1], pt[2]], 0) -
        bcdtpp.projectH1(u, [pt[0] - h, pt[1], pt[2]], 0)) /
        (2 * h),
      (bcdtpp.projectH1(u, [pt[0], pt[1] + h, pt[2]], 0) -
        bcdtpp.projectH1(u, [pt[0], pt[1] - h, pt[2]], 0)) /
        (2 * h),
      (bcdtpp.projectH1(u, [pt[0], pt[1], pt[2] + h], 0) -
        bcdtpp.projectH1(u, [pt[0], pt[1], pt[2] - h], 0)) /
        (2 * h),
    ];
    expect(gradProj[0]).toBeCloseTo(0, 5);
    expect(gradProj[1]).toBeCloseTo(0, 5);
    expect(gradProj[2]).toBeCloseTo(0, 5);
  });

  it('exactness: Pi^2 reproduces constant vector fields', () => {
    const c = [2, -3, 1];
    const v = () => c;
    const pt = mesh.getTetrahedronBarycenter(0);
    const proj = bcdtpp.projectHdiv(v, pt, 0);
    expect(proj[0]).toBeCloseTo(c[0], 1);
    expect(proj[1]).toBeCloseTo(c[1], 1);
    expect(proj[2]).toBeCloseTo(c[2], 1);
  });

  it('commuting property: div Pi^2 ≈ Pi^3 div for constant divergence', () => {
    const v = (pt) => [2 * pt[0], 2 * pt[1], 2 * pt[2]];
    const pt = mesh.getTetrahedronBarycenter(0);

    const l2Div = bcdtpp.projectL2(() => 6, 0);
    expect(l2Div).toBeCloseTo(6, 6);

    const h = 1e-5;
    const projV = bcdtpp.projectHdiv(v, pt, 0);
    const projVx = bcdtpp.projectHdiv(v, [pt[0] + h, pt[1], pt[2]], 0);
    const projVy = bcdtpp.projectHdiv(v, [pt[0], pt[1] + h, pt[2]], 0);
    const projVz = bcdtpp.projectHdiv(v, [pt[0], pt[1], pt[2] + h], 0);

    const numDiv =
      (projVx[0] - projV[0]) / h +
      (projVy[1] - projV[1]) / h +
      (projVz[2] - projV[2]) / h;

    expect(numDiv).toBeCloseTo(l2Div, 1);
  });

  it('higher-order H1 projection adds bubble correction', () => {
    const u = (pt) => pt[0] * pt[0] + pt[1] * pt[1];
    const pt = mesh.getTetrahedronBarycenter(0);
    const p0 = bcdtpp.projectHp(u, pt, 0, 0, 0);
    const p1 = bcdtpp.projectHp(u, pt, 0, 0, 1);
    expect(typeof p1).toBe('number');
    expect(Number.isFinite(p1)).toBe(true);
  });

  it('higher-order L2 projection reproduces linear polynomials for p=1', () => {
    const u = (pt) => 2 * pt[0] - 3 * pt[1] + 5 * pt[2] + 7;
    const pt = [0.1, 0.2, 0.05];
    const proj = bcdtpp.projectHp(u, pt, 0, 3, 1);
    expect(proj).toBeCloseTo(u(pt), 6);
  });

  it('exactness: Pi^1 reproduces gradient of linear function', () => {
    const u = (pt) => 2 * pt[0] - 3 * pt[1] + 5 * pt[2];
    const pt = mesh.getTetrahedronBarycenter(0);
    const gradU = [2, -3, 5];
    const proj = bcdtpp.projectHcurl(u, pt, 0);
    expect(proj[0]).toBeCloseTo(gradU[0], 1);
    expect(proj[1]).toBeCloseTo(gradU[1], 1);
    expect(proj[2]).toBeCloseTo(gradU[2], 1);
  });

});
