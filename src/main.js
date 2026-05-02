import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import * as math from 'mathjs';

import { Mesh } from './lib/Mesh';
import { Whitney } from './lib/Whitney';
import { BCDTPP } from './lib/BCDTPP';

// --- Setup Three.js ---
const container = document.getElementById('canvas-container');
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(45, container.clientWidth / container.clientHeight, 0.1, 100);
camera.position.set(3, 2, 4);

const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: true });
renderer.setSize(container.clientWidth, container.clientHeight);
renderer.setPixelRatio(window.devicePixelRatio);
container.appendChild(renderer.domElement);

const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;

const ambientLight = new THREE.AmbientLight(0xffffff, 0.6);
scene.add(ambientLight);
const dirLight = new THREE.DirectionalLight(0xffffff, 0.8);
dirLight.position.set(5, 10, 7);
scene.add(dirLight);

// --- Core state ---
let currentMesh, currentWhitney, currentBCDTPP;
let meshGroup = new THREE.Group();
scene.add(meshGroup);

const defaultMesh = {
    vertices: [[0,0,0], [1,0,0], [0,1,0], [0,0,1]],
    tetrahedra: [[0,1,2,3]]
};

async function loadMesh(data) {
    meshGroup.clear();
    currentMesh = new Mesh(data.vertices, data.tetrahedra);
    currentWhitney = new Whitney(currentMesh);
    currentBCDTPP = new BCDTPP(currentMesh, currentWhitney);

    console.log("Computing splits and weights...");
    await currentBCDTPP.computeBoundaryWeights();

    // Visuals
    const geom = new THREE.BufferGeometry();
    const positions = [];
    currentMesh.tetrahedra.forEach(tet => {
        const v = tet.map(i => currentMesh.vertices[i]);
        const faces = [[1,2,3], [0,2,3], [0,1,3], [0,1,2]];
        faces.forEach(f => {
            f.forEach(vi => positions.push(...v[vi]));
        });
    });
    geom.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    geom.computeVertexNormals();
    
    // Visualize Alfeld split lines
    const alfeldLines = [];
    currentMesh.alfeldTriangles.forEach(at => {
        at.triangles.forEach(tri => {
            const v = tri.map(i => currentMesh.vertices[i]);
            alfeldLines.push(...v[0], ...v[1], ...v[1], ...v[2], ...v[2], ...v[0]);
        });
    });
    const alfeldGeom = new THREE.BufferGeometry();
    alfeldGeom.setAttribute('position', new THREE.Float32BufferAttribute(alfeldLines, 3));
    const alfeldMesh = new THREE.LineSegments(alfeldGeom, new THREE.LineBasicMaterial({ color: 0xfab005 }));
    meshGroup.add(alfeldMesh);

    const mat = new THREE.MeshPhysicalMaterial({ color: 0x6366f1, transparent: true, opacity: 0.4, roughness: 0.3, metalness: 0.1 });
    const mesh = new THREE.Mesh(geom, mat);
    meshGroup.add(mesh);
    
    evaluate();
}

function evaluate() {
    const expr = document.getElementById('func-input').value;
    try {
        const compiled = math.compile(expr);
        const u = (pt) => compiled.evaluate({x: pt[0], y: pt[1], z: pt[2]});
        
        const pt = currentMesh.getTetBarycenter(0);
        
        const h1 = currentBCDTPP.projectH1(u, pt, 0);
        const hcurl = currentBCDTPP.projectHcurl(u, pt, 0);
        const hdiv = currentBCDTPP.projectHdiv(u, pt, 0);
        const l2 = currentBCDTPP.projectL2(u, 0);

        document.getElementById('val-h1').innerText = h1.toFixed(6);
        document.getElementById('val-hcurl').innerText = math.norm(hcurl).toFixed(6);
        document.getElementById('val-hdiv').innerText = math.norm(hdiv).toFixed(6);
        document.getElementById('val-l2').innerText = l2.toFixed(6);
    } catch(e) {
        console.error(e);
    }
}

// --- Events ---
document.getElementById('func-input').addEventListener('input', evaluate);

document.getElementById('mesh-upload').addEventListener('change', (e) => {
    const file = e.target.files[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = (ev) => {
        const data = JSON.parse(ev.target.result);
        loadMesh(data);
        document.getElementById('mesh-status').innerText = `Loaded: ${file.name}`;
    };
    reader.readAsText(file);
});

window.addEventListener('resize', () => {
    camera.aspect = container.clientWidth / container.clientHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(container.clientWidth, container.clientHeight);
});

// Start
loadMesh(defaultMesh);

function animate() {
    requestAnimationFrame(animate);
    controls.update();
    renderer.render(scene, camera);
}
animate();
