import type {Mesh} from '../mesh.js';
import type {Whitney} from '../whitney.js';

export class L2Projector {
  constructor(mesh: Mesh, whitney: Whitney, quadratureOrder: number);
  project(u: (point: number[]) => number, tIdx: number): number;
}
