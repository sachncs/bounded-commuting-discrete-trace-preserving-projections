import {Mesh} from './lib/mesh.js';
import {Whitney} from './lib/whitney.js';
import {Bcdtpp} from './lib/bcdtpp.js';
import * as renderer from './renderer.js';
import * as ui from './ui.js';

// --- Core state ---
let currentMesh;
let currentWhitney;
let currentBcdtpp;
let meshGroup;

const defaultMesh = {
  vertices: [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]],
  tetrahedra: [[0, 1, 2, 3]],
};

function loadMesh(data) {
  ui.validateMeshData(data);

  renderer.clearMeshGroup(meshGroup);

  try {
    currentMesh = new Mesh(data.vertices, data.tetrahedra);
    currentWhitney = new Whitney(currentMesh);
    currentBcdtpp = new Bcdtpp(currentMesh, currentWhitney, {
      quadratureOrder: 3,
    });

    console.log('Computing splits and weights...');
    currentBcdtpp.computeBoundaryWeights();
    currentBcdtpp.buildPointLocator();

    renderer.renderMesh(currentMesh, currentBcdtpp.meshRefinement, meshGroup);

    ui.setStatus(
      `Tets: ${currentMesh.tetrahedronCount}, Verts: ${currentMesh.vertexCount}`,
    );
  } catch (e) {
    console.error(e);
    ui.setStatus(`Error: ${e.message}`);
    return;
  }

  evaluate();
}

function evaluate() {
  if (!currentBcdtpp) {
    return;
  }
  const expr = ui.getExpression();
  const values = ui.evaluateExpression(expr, currentBcdtpp, currentMesh, (e) => {
    console.error('Expression evaluation failed:', e.message);
  });
  ui.updateResults(values);
}

// --- Init ---
const rendererResult = renderer.initRenderer('canvas-container');
if (!rendererResult) {
  ui.setStatus('Error: Canvas container #canvas-container not found');
  console.error('Failed to initialize renderer: #canvas-container not found');
} else {
  meshGroup = rendererResult.meshGroup;
}

ui.setupUI({
  onExpressionChange: evaluate,
  onMeshUpload: (data, filename) => {
    loadMesh(data);
    ui.setStatus(`Loaded: ${filename}`);
  },
});

loadMesh(defaultMesh);
