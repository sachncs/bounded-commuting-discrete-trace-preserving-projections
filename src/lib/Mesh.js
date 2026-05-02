import * as math from 'mathjs';

export class Mesh {
    constructor(vertices, tetrahedra) {
        this.vertices = vertices; // [[x,y,z], ...]
        this.tetrahedra = tetrahedra; // [[v0,v1,v2,v3], ...]
        this.nodes = vertices.length;
        this.numTets = tetrahedra.length;
        
        this.faces = [];
        this.edges = [];
        this.boundaryFaces = [];
        this.boundaryEdges = [];
        this.boundaryNodes = new Set();
        
        this.faceToTets = [];
        this.edgeToFaces = [];
        
        this._buildTopology();
        
        // Split data
        this.faceBarycenters = []; // Indices in vertices
        this.tetBarycenters = [];  // Indices in vertices
        this.alfeldTriangles = []; // { parentFaceIdx, triangles: [[v0,v1,v_f], ...] }
        this.worseyTetrahedra = []; // { parentTetIdx, tetrahedra: [[v0,v1,v2,v3], ...] }
    }

    _buildTopology() {
        const faceMap = new Map();
        const getFaceKey = (f) => [...f].sort((a, b) => a - b).join(',');

        for (let tIdx = 0; tIdx < this.numTets; tIdx++) {
            const tet = this.tetrahedra[tIdx];
            const localFaces = [
                [tet[1], tet[2], tet[3]], // opposite v0
                [tet[0], tet[2], tet[3]], // opposite v1
                [tet[0], tet[1], tet[3]], // opposite v2
                [tet[0], tet[1], tet[2]]  // opposite v3
            ];

            for (const f of localFaces) {
                const key = getFaceKey(f);
                if (!faceMap.has(key)) {
                    faceMap.set(key, { verts: f, tets: [tIdx] });
                } else {
                    faceMap.get(key).tets.push(tIdx);
                }
            }
        }

        const edgeMap = new Map();
        const getEdgeKey = (e) => [...e].sort((a, b) => a - b).join(',');

        this._faceKeyToIdx = new Map();
        
        let fIdx = 0;
        for (const [key, data] of faceMap.entries()) {
            this.faces.push(data.verts);
            this.faceToTets.push(data.tets);
            this._faceKeyToIdx.set(key, fIdx);
            
            const isBoundary = data.tets.length === 1;
            if (isBoundary) {
                this.boundaryFaces.push(fIdx);
                data.verts.forEach(v => this.boundaryNodes.add(v));
            }

            const fEdges = [
                [data.verts[0], data.verts[1]],
                [data.verts[1], data.verts[2]],
                [data.verts[2], data.verts[0]]
            ];

            for (const e of fEdges) {
                const eKey = getEdgeKey(e);
                if (!edgeMap.has(eKey)) {
                    edgeMap.set(eKey, { verts: e, isBoundary, faces: [fIdx] });
                } else {
                    const eData = edgeMap.get(eKey);
                    eData.faces.push(fIdx);
                    if (isBoundary) eData.isBoundary = true;
                }
            }
            fIdx++;
        }

        let eIdx = 0;
        for (const [key, data] of edgeMap.entries()) {
            this.edges.push(data.verts);
            this.edgeToFaces.push(data.faces);
            if (data.isBoundary) {
                this.boundaryEdges.push(eIdx);
            }
            eIdx++;
        }
    }

    // Section 6.1.3: Alfeld Split of boundary faces
    computeAlfeldSplit() {
        this.faceBarycenters = new Array(this.faces.length).fill(-1);
        this.alfeldTriangles = [];

        this.boundaryFaces.forEach(fIdx => {
            const f = this.faces[fIdx];
            const bary = this.getFaceBarycenter(fIdx);
            const vBaryIdx = this.vertices.length;
            this.vertices.push(bary);
            this.faceBarycenters[fIdx] = vBaryIdx;

            const subTriangles = [
                [f[0], f[1], vBaryIdx],
                [f[1], f[2], vBaryIdx],
                [f[2], f[0], vBaryIdx]
            ];
            this.alfeldTriangles.push({ parentFaceIdx: fIdx, triangles: subTriangles });
        });
    }

    // Section 6.1.4: Worsey-Farin Split of tetrahedra
    computeWorseyFarinSplit() {
        if (this.faceBarycenters.length === 0) this.computeAlfeldSplit();

        this.tetBarycenters = new Array(this.numTets).fill(-1);
        this.worseyTetrahedra = [];

        for (let tIdx = 0; tIdx < this.numTets; tIdx++) {
            const tet = this.tetrahedra[tIdx];
            const bary = this.getTetBarycenter(tIdx);
            const vTetBaryIdx = this.vertices.length;
            this.vertices.push(bary);
            this.tetBarycenters[tIdx] = vTetBaryIdx;

            // Each face of the tet gives 3 sub-triangles (Alfeld split)
            // Each sub-triangle + tet barycenter = sub-tetrahedron
            // Total 4 faces * 3 triangles = 12 sub-tets
            const tetSubTets = [];
            const tFaces = this._getTetFaces(tIdx);
            
            tFaces.forEach(fIdx => {
                const f = this.faces[fIdx];
                let fvBaryIdx = this.faceBarycenters[fIdx];
                if (fvBaryIdx === -1) {
                    // Face was not on boundary, but we need barycenter for WF split
                    const fbary = this.getFaceBarycenter(fIdx);
                    fvBaryIdx = this.vertices.length;
                    this.vertices.push(fbary);
                    this.faceBarycenters[fIdx] = fvBaryIdx;
                }

                const subTris = [
                    [f[0], f[1], fvBaryIdx],
                    [f[1], f[2], fvBaryIdx],
                    [f[2], f[0], fvBaryIdx]
                ];

                subTris.forEach(tri => {
                    tetSubTets.push([...tri, vTetBaryIdx]);
                });
            });
            this.worseyTetrahedra.push({ parentTetIdx: tIdx, tetrahedra: tetSubTets });
        }
    }

    getFaceBarycenter(fIdx) {
        const f = this.faces[fIdx];
        const v = f.map(i => this.vertices[i]);
        return [
            (v[0][0]+v[1][0]+v[2][0])/3,
            (v[0][1]+v[1][1]+v[2][1])/3,
            (v[0][2]+v[1][2]+v[2][2])/3
        ];
    }

    getTetBarycenter(tIdx) {
        const t = this.tetrahedra[tIdx];
        const v = t.map(i => this.vertices[i]);
        return [
            (v[0][0]+v[1][0]+v[2][0]+v[3][0])/4,
            (v[0][1]+v[1][1]+v[2][1]+v[3][1])/4,
            (v[0][2]+v[1][2]+v[2][2]+v[3][2])/4
        ];
    }

    getStar(vIdx) {
        const star = [];
        for (let i = 0; i < this.numTets; i++) {
            if (this.tetrahedra[i].includes(vIdx)) star.push(i);
        }
        return star;
    }

    _getTetFaces(tIdx) {
        const tet = this.tetrahedra[tIdx];
        const local = [[tet[1],tet[2],tet[3]], [tet[0],tet[2],tet[3]], [tet[0],tet[1],tet[3]], [tet[0],tet[1],tet[2]]];
        return local.map(lf => {
            const key = [...lf].sort((a,b)=>a-b).join(',');
            // We need a map from key to faceIdx. I'll add it to buildTopology.
            return this._faceKeyToIdx.get(key);
        });
    }
}
