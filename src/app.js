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

async function loadMesh(data) {
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

    renderer.renderMesh(currentMesh, meshGroup);

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
  const values = ui.evaluateExpression(expr, currentBcdtpp, currentMesh);
  ui.updateResults(values);
}

// --- Init ---
const {meshGroup: mg} = renderer.initRenderer('canvas-container');
meshGroup = mg;

ui.setupUI({
  onExpressionChange: evaluate,
  onMeshUpload: (data, filename) => {
    loadMesh(data);
    ui.setStatus(`Loaded: ${filename}`);
  },
});

loadMesh(defaultMesh);
