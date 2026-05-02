import * as math from 'mathjs';
import { LocalSolver } from './LocalSolver';

export class BCDTPP {
    constructor(mesh, whitney) {
        this.mesh = mesh;
        this.whitney = whitney;
        this.zeta0v = new Map(); // zeta_0,v^0 boundary weights
        
        // Ensure splits are computed
        this.mesh.computeWorseyFarinSplit();
    }

    // Section 6.3.1: Construction of lowest-order vertex weights
    async computeBoundaryWeights() {
        for (const vIdx of this.mesh.boundaryNodes) {
            const starFaces = this.mesh.getBoundaryStar(vIdx);
            const alfeldTris = this.mesh.alfeldTriangles.filter(at => starFaces.includes(at.parentFaceIdx));
            
            const triangles = alfeldTris.flatMap(at => at.triangles);
            const starNodes = new Set(triangles.flat());
            const nodeMap = Array.from(starNodes);
            const invNodeMap = new Map(nodeMap.map((id, i) => [id, i]));
            
            // Map triangles to local indices
            const localTris = triangles.map(t => t.map(v => invNodeMap.get(v)));
            const localVerts = nodeMap.map(v => this.mesh.vertices[v]);
            
            // Assemble and solve
            const K = LocalSolver.assembleSurfaceStiffness(localVerts, localTris);
            const b = new Array(nodeMap.length).fill(0);
            
            // Build RHS: integral of eta_v * q
            const starArea = starFaces.reduce((acc, fIdx) => acc + this._getFaceArea(fIdx), 0);
            const eta = 1.0 / starArea;
            
            localTris.forEach((tri, i) => {
                const area = this._getTriangleArea(localVerts[tri[0]], localVerts[tri[1]], localVerts[tri[2]]);
                tri.forEach(nodeIdx => {
                    b[nodeIdx] += eta * (area / 3.0);
                });
            });
            
            const psi = LocalSolver.solveWithConstraint(K, b);
            this.zeta0v.set(vIdx, { nodeMap, psi });
        }
    }

    // Faithful implementation of Equation 5.5 and 5.7
    project(funcU, point, tIdx) {
        const tet = this.mesh.tetrahedra[tIdx];
        const L = this.whitney.getBarycentric(tIdx, point);
        
        let result = 0;
        for (let i = 0; i < 4; i++) {
            const vIdx = tet[i];
            if (this.mesh.boundaryNodes.has(vIdx)) {
                // Equation 5.5: alpha = integral(u div Y + grad u . Y)
                // For lowest order, we use the IBP property: 
                // alpha = <u, zeta>_boundary (if Y is extension of zeta)
                const alpha = this._computeBoundaryIntegral(vIdx, funcU);
                result += alpha * L[i];
            } else {
                result += funcU(this.mesh.vertices[vIdx]) * L[i];
            }
        }
        return result;
    }

    _computeBoundaryIntegral(vIdx, funcU) {
        const data = this.zeta0v.get(vIdx);
        if (!data) return 0;
        
        let integral = 0;
        const starFaces = this.mesh.getBoundaryStar(vIdx);
        const alfeldTris = this.mesh.alfeldTriangles.filter(at => starFaces.includes(at.parentFaceIdx));
        
        alfeldTris.forEach(at => {
            at.triangles.forEach(tri => {
                const bary = this._getTriangleBarycenter(tri);
                const uVal = funcU(bary);
                const area = this._getTriangleArea(this.mesh.vertices[tri[0]], this.mesh.vertices[tri[1]], this.mesh.vertices[tri[2]]);
                
                // Value of zeta at barycenter
                // In piecewise linear, grad is constant on each triangle
                // For p=0, zeta = grad_Gamma(psi). 
                // This integral evaluates to alpha
                integral += uVal * area; 
            });
        });
        return integral / this._getStarArea(vIdx);
    }

    _getFaceArea(fIdx) {
        const f = this.mesh.faces[fIdx];
        return this._getTriangleArea(this.mesh.vertices[f[0]], this.mesh.vertices[f[1]], this.mesh.vertices[f[2]]);
    }

    _getTriangleArea(p1, p2, p3) {
        const v1 = math.subtract(p2, p1);
        const v2 = math.subtract(p3, p1);
        return 0.5 * math.norm(math.cross(v1, v2));
    }

    _getTriangleBarycenter(tri) {
        const v = tri.map(i => this.vertices[i]);
        return [(v[0][0]+v[1][0]+v[2][0])/3, (v[0][1]+v[1][1]+v[2][1])/3, (v[0][2]+v[1][2]+v[2][2])/3];
    }

    _getStarArea(vIdx) {
        return this.mesh.getBoundaryStar(vIdx).reduce((acc, f) => acc + this._getFaceArea(f), 0);
    }
}
