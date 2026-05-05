import * as THREE from 'three';
import {OrbitControls} from 'three/examples/jsm/controls/OrbitControls.js';

/**
 * Initializes the Three.js renderer, scene, camera, and controls.
 * @param {string} containerId
 * @return {{scene: THREE.Scene, camera: THREE.PerspectiveCamera, renderer: THREE.WebGLRenderer, controls: OrbitControls, meshGroup: THREE.Group, container: HTMLElement}}
 */
export function initRenderer(containerId) {
  const container = document.getElementById(containerId);
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

  const meshGroup = new THREE.Group();
  scene.add(meshGroup);

  window.addEventListener('resize', () => {
    camera.aspect = container.clientWidth / container.clientHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(container.clientWidth, container.clientHeight);
  });

  function animate() {
    requestAnimationFrame(animate);
    controls.update();
    renderer.render(scene, camera);
  }
  animate();

  return {scene, camera, renderer, controls, meshGroup, container};
}

/**
 * Clears and disposes all objects in a Three.js group.
 * @param {THREE.Group} meshGroup
 */
export function clearMeshGroup(meshGroup) {
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
}

/**
 * Builds a BufferGeometry from tetrahedron faces.
 * @param {!Object} mesh
 * @return {THREE.BufferGeometry}
 */
function buildTetrahedronGeometry(mesh) {
  const positions = [];
  mesh.tetrahedra.forEach((tet) => {
    const v = tet.map((i) => mesh.vertices[i]);
    const faces = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
    faces.forEach((f) => {
      f.forEach((vi) => positions.push(...v[vi]));
    });
  });
  const geom = new THREE.BufferGeometry();
  geom.setAttribute(
    'position',
    new THREE.Float32BufferAttribute(positions, 3),
  );
  geom.computeVertexNormals();
  return geom;
}

/**
 * Builds a LineSegments geometry from Alfeld-split triangles.
 * @param {!Object} mesh
 * @param {!Object} meshRefinement
 * @return {THREE.BufferGeometry}
 */
function buildAlfeldGeometry(mesh, meshRefinement) {
  const lines = [];
  meshRefinement.alfeldTriangles.forEach((at) => {
    at.triangles.forEach((tri) => {
      const v = tri.map((i) => mesh.vertices[i]);
      lines.push(...v[0], ...v[1], ...v[1], ...v[2], ...v[2], ...v[0]);
    });
  });
  const geom = new THREE.BufferGeometry();
  geom.setAttribute(
    'position',
    new THREE.Float32BufferAttribute(lines, 3),
  );
  return geom;
}

/**
 * Renders a mesh into the Three.js group.
 * @param {!Object} mesh
 * @param {!Object} meshRefinement
 * @param {THREE.Group} meshGroup
 */
export function renderMesh(mesh, meshRefinement, meshGroup) {
  const geom = buildTetrahedronGeometry(mesh);
  const mat = new THREE.MeshPhysicalMaterial({
    color: 0x6366f1,
    transparent: true,
    opacity: 0.4,
    roughness: 0.3,
    metalness: 0.1,
  });
  meshGroup.add(new THREE.Mesh(geom, mat));

  const alfeldGeom = buildAlfeldGeometry(mesh, meshRefinement);
  const alfeldMat = new THREE.LineBasicMaterial({color: 0xfab005});
  meshGroup.add(new THREE.LineSegments(alfeldGeom, alfeldMat));
}
