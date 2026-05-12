import type {Mesh} from '../mesh.js';
import type {Whitney} from '../whitney.js';

export class HdivProjector {
  constructor(mesh: Mesh, whitney: Whitney, quadratureOrder: number);
  project(
    u: (point: number[]) => number | number[],
    point: number[],
    tIdx: number,
    boundaryFaceSet: Set<number>,
  ): number[];
  computeFaceDof(
    u: (point: number[]) => number | number[],
    fIdx: number,
  ): number;
  computeInteriorFaceCoeff(
    u: (point: number[]) => number | number[],
    tIdx: number,
    f: number,
    isScalar: boolean,
  ): number;
  projectRing(
    u: (point: number[]) => number | number[],
    point: number[],
    tIdx: number,
    boundaryFaceSet: Set<number>,
  ): number[];
  extendBoundary(
    boundaryData: Map<number, number>,
    point: number[],
    tIdx: number,
    boundaryFaceSet: Set<number>,
  ): number[];
}
