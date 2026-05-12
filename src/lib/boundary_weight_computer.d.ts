import type {Mesh} from './mesh.js';
import type {MeshRefinement} from './mesh_refinement.js';

/**
 * Boundary weight computation for BCDTPP projections.
 *
 * Computes vertex patch weights, edge tangents/lengths, and face normals/areas
 * used by trace-preserving boundary DoFs.
 */
export class BoundaryWeightComputer {
  /**
   * @param mesh - The mesh.
   * @param meshRefinement - Refinement data structure.
   * @param onWarning - Optional warning callback.
   */
  constructor(
    mesh: Mesh,
    meshRefinement: MeshRefinement,
    onWarning?: (ctx: {code: string; severity: 'warn' | 'error'; message: string}) => void,
  );

  /**
   * Computes all boundary weights.
   * @returns Object containing vertexBoundaryData, edgeBoundaryData, and
   *   faceBoundaryData maps.
   */
  compute(): {
    vertexBoundaryData: Map<
      number,
      {nodeMap: number[]; invNodeMap: Map<number, number>; psi: number[]}
    >;
    edgeBoundaryData: Map<
      number,
      {v0: number; v1: number; tangent: number[]; length: number}
    >;
    faceBoundaryData: Map<number, {normal: number[]; area: number}>;
  };
}
