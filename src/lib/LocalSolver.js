import * as math from 'mathjs';

export class LocalSolver {
    static assembleSurfaceStiffness(vertices, triangles) {
        const n = vertices.length;
        const K = math.zeros(n, n, 'sparse');
        
        triangles.forEach(tri => {
            const v = tri.map(i => vertices[i]);
            const ke = this._triangleStiffness(v);
            
            for (let i = 0; i < 3; i++) {
                for (let j = 0; j < 3; j++) {
                    const val = K.get([tri[i], tri[j]]);
                    K.set([tri[i], tri[j]], val + ke[i][j]);
                }
            }
        });
        
        return K;
    }

    static _triangleStiffness(v) {
        // Local stiffness matrix for -\Delta_\Gamma on a triangle
        const v1 = math.subtract(v[1], v[0]);
        const v2 = math.subtract(v[2], v[0]);
        
        const cross = math.cross(v1, v2);
        const area = 0.5 * math.norm(cross);
        
        // Cotangent formula or Jacobian-based
        const G = [
            math.subtract(v[2], v[1]),
            math.subtract(v[0], v[2]),
            math.subtract(v[1], v[0])
        ];
        
        const ke = math.zeros(3, 3).toArray();
        for (let i = 0; i < 3; i++) {
            for (let j = 0; j < 3; j++) {
                ke[i][j] = math.dot(G[i], G[j]) / (4 * area);
            }
        }
        return ke;
    }

    static solveWithConstraint(K, b) {
        // Apply mean-zero constraint: sum(x) = 0
        // We replace the last row with ones and set b[last] = 0
        // (This is a common trick for the singular Laplacian)
        const n = K.size()[0];
        const K_dense = K.toArray();
        for (let j = 0; j < n; j++) {
            K_dense[n-1][j] = 1.0;
        }
        const b_mod = [...b];
        b_mod[n-1] = 0;
        
        return math.lusolve(K_dense, b_mod).map(x => x[0]);
    }
}
