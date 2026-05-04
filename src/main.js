import * as THREE from 'three';
import {OrbitControls} from 'three/examples/jsm/controls/OrbitControls.js';

import {Mesh} from './lib/mesh.js';
import {Whitney} from './lib/whitney.js';
import {Bcdtpp} from './lib/bcdtpp.js';
import {norm} from './lib/math_utils.js';

// --- Setup Three.js ---
const container = document.getElementById('canvas-container');
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(
  45,
  container.clientWidth / container.clientHeight,
  0.1,
  100,
);
camera.position.set(3, 2, 4);

const renderer = new THREE.WebGLRenderer({antialias: true, alpha: true});
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
let currentMesh;
let currentWhitney;
let currentBcdtpp;
const meshGroup = new THREE.Group();
scene.add(meshGroup);

const defaultMesh = {
  vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
  tetrahedra: [[0, 1, 2, 3]],
};

/**
 * Compiles a math expression string into an evaluable function.
 * @param {string} expr
 * @return {function(number, number, number): number}
 */
function compileExpression(expr) {
  const mathFunctions = [
    'sin', 'cos', 'tan', 'asin', 'acos', 'atan', 'atan2',
    'exp', 'log', 'log10', 'log2', 'sqrt', 'cbrt',
    'abs', 'ceil', 'floor', 'round', 'trunc',
    'pow', 'min', 'max', 'sign', 'PI', 'E',
  ];
  let wrapped = expr;
  for (const fn of mathFunctions) {
    const re = new RegExp(`\\b${fn}\\b`, 'g');
    wrapped = wrapped.replace(re, `Math.${fn}`);
  }
  return new Function('x', 'y', 'z', `return (${wrapped});`);
}

async function loadMesh(data) {
  meshGroup.traverse((obj) => {
    if (obj.geometry) obj.geometry.dispose();
    if (obj.material) {
      if (Array.isArray(obj.material)) {
        obj.material.forEach((m) => m.dispose());
      } else {
        obj.material.dispose();
      }
    }
  });
  meshGroup.clear();
  try {
    currentMesh = new Mesh(data.vertices, data.tetrahedra);
    currentWhitney = new Whitney(currentMesh);
    currentBcdtpp = new Bcdtpp(currentMesh, currentWhitney, {
      quadratureOrder: 3,
    });

    console.log('Computing splits and weights...');
    currentBcdtpp.computeBoundaryWeights();
    currentBcdtpp.buildPointLocator();

    const geom = new THREE.BufferGeometry();
    const positions = [];
    currentMesh.tetrahedra.forEach((tet) => {
      const v = tet.map((i) => currentMesh.vertices[i]);
      const faces = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
      faces.forEach((f) => {
        f.forEach((vi) => positions.push(...v[vi]));
      });
    });
    geom.setAttribute(
      'position',
      new THREE.Float32BufferAttribute(positions, 3),
    );
    geom.computeVertexNormals();

    const alfeldLines = [];
    currentMesh.alfeldTriangles.forEach((at) => {
      at.triangles.forEach((tri) => {
        const v = tri.map((i) => currentMesh.vertices[i]);
        alfeldLines.push(
          ...v[0], ...v[1],
          ...v[1], ...v[2],
          ...v[2], ...v[0],
        );
      });
    });
    const alfeldGeom = new THREE.BufferGeometry();
    alfeldGeom.setAttribute(
      'position',
      new THREE.Float32BufferAttribute(alfeldLines, 3),
    );
    const alfeldMesh = new THREE.LineSegments(
      alfeldGeom,
      new THREE.LineBasicMaterial({color: 0xfab005}),
    );
    meshGroup.add(alfeldMesh);

    const mat = new THREE.MeshPhysicalMaterial({
      color: 0x6366f1,
      transparent: true,
      opacity: 0.4,
      roughness: 0.3,
      metalness: 0.1,
    });
    const mesh = new THREE.Mesh(geom, mat);
    meshGroup.add(mesh);

    document.getElementById('mesh-status').innerText =
      `Tets: ${currentMesh.tetrahedronCount}, Verts: ${currentMesh.vertexCount}`;
  } catch (e) {
    console.error(e);
    document.getElementById('mesh-status').innerText = `Error: ${e.message}`;
    return;
  }

  evaluate();
}

function evaluate() {
  if (!currentBcdtpp) {
    return;
  }
  const expr = document.getElementById('func-input').value;
  try {
    const u = compileExpression(expr);

    const pt = currentMesh.getTetrahedronBarycenter(0);
    const found = currentBcdtpp.pointLocator.findTetrahedron(pt);
    const tIdx = found ? found.tIdx : 0;

    const h1 = currentBcdtpp.projectH1(u, pt, tIdx);
    const hcurl = currentBcdtpp.projectHcurl(u, pt, tIdx);
    const hdiv = currentBcdtpp.projectHdiv(u, pt, tIdx);
    const l2 = currentBcdtpp.projectL2(u, tIdx);

    document.getElementById('val-h1').innerText = h1.toFixed(6);
    document.getElementById('val-hcurl').innerText = norm(hcurl).toFixed(6);
    document.getElementById('val-hdiv').innerText = norm(hdiv).toFixed(6);
    document.getElementById('val-l2').innerText = l2.toFixed(6);
  } catch (e) {
    console.error(e);
    document.getElementById('val-h1').innerText = 'ERR';
    document.getElementById('val-hcurl').innerText = 'ERR';
    document.getElementById('val-hdiv').innerText = 'ERR';
    document.getElementById('val-l2').innerText = 'ERR';
  }
}

// --- Events ---
document.getElementById('func-input').addEventListener('input', evaluate);

document.getElementById('mesh-upload').addEventListener('change', (e) => {
  const file = e.target.files[0];
  if (!file) {
    return;
  }
  const reader = new FileReader();
  reader.onload = (ev) => {
    try {
      const data = JSON.parse(ev.target.result);
      loadMesh(data);
      document.getElementById('mesh-status').innerText = `Loaded: ${file.name}`;
    } catch (err) {
      document.getElementById('mesh-status').innerText =
        `Parse error: ${err.message}`;
    }
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
