import * as math from 'mathjs';

export class Whitney {
    constructor(mesh) {
        this.mesh = mesh;
    }

    getBarycentric(tIdx, point) {
        const tet = this.mesh.tetrahedra[tIdx];
        const v = tet.map(i => this.mesh.vertices[i]);
        
        const T = [
            [v[0][0]-v[3][0], v[1][0]-v[3][0], v[2][0]-v[3][0]],
            [v[0][1]-v[3][1], v[1][1]-v[3][1], v[2][1]-v[3][1]],
            [v[0][2]-v[3][2], v[1][2]-v[3][2], v[2][2]-v[3][2]]
        ];
        
        const b = [point[0]-v[3][0], point[1]-v[3][1], point[2]-v[3][2]];
        try {
            const L = math.lusolve(T, b).map(x => x[0]);
            return [L[0], L[1], L[2], 1 - L[0] - L[1] - L[2]];
        } catch(e) {
            return [0,0,0,0];
        }
    }

    getGradBarycentric(tIdx) {
        const tet = this.mesh.tetrahedra[tIdx];
        const v = tet.map(i => this.mesh.vertices[i]);
        const T = [
            [v[0][0]-v[3][0], v[1][0]-v[3][0], v[2][0]-v[3][0]],
            [v[0][1]-v[3][1], v[1][1]-v[3][1], v[2][1]-v[3][1]],
            [v[0][2]-v[3][2], v[1][2]-v[3][2], v[2][2]-v[3][2]]
        ];
        try {
            const Tinv = math.inv(T);
            const g0 = [Tinv[0][0], Tinv[1][0], Tinv[2][0]];
            const g1 = [Tinv[0][1], Tinv[1][1], Tinv[2][1]];
            const g2 = [Tinv[0][2], Tinv[1][2], Tinv[2][2]];
            const g3 = [-g0[0]-g1[0]-g2[0], -g0[1]-g1[1]-g2[1], -g0[2]-g1[2]-g2[2]];
            return [g0, g1, g2, g3];
        } catch(e) {
            return [[0,0,0],[0,0,0],[0,0,0],[0,0,0]];
        }
    }
}
