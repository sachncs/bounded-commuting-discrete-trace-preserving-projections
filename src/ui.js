import {compileExpression} from './lib/safe_eval.js';
import {norm} from './lib/math_utils.js';

/**
 * Sets up DOM event listeners.
 * @param {{onExpressionChange: function(): void, onMeshUpload: function(!Object, string): void}} callbacks
 */
export function setupUI({onExpressionChange, onMeshUpload}) {
  document.getElementById('func-input').addEventListener('input', onExpressionChange);

  document.getElementById('mesh-upload').addEventListener('change', (e) => {
    const file = e.target.files[0];
    if (!file) {
      return;
    }
    const reader = new FileReader();
    reader.onload = (ev) => {
      try {
        const data = JSON.parse(ev.target.result);
        onMeshUpload(data, file.name);
      } catch (err) {
        setStatus(`Parse error: ${err.message}`);
      }
    };
    reader.readAsText(file);
  });
}

/**
 * @param {string} text
 */
export function setStatus(text) {
  document.getElementById('mesh-status').innerText = text;
}

/**
 * Updates the result display fields.
 * @param {{h1: string, hcurl: string, hdiv: string, l2: string}} values
 */
export function updateResults(values) {
  document.getElementById('val-h1').innerText = values.h1;
  document.getElementById('val-hcurl').innerText = values.hcurl;
  document.getElementById('val-hdiv').innerText = values.hdiv;
  document.getElementById('val-l2').innerText = values.l2;
}

/**
 * @return {string}
 */
export function getExpression() {
  return document.getElementById('func-input').value;
}

/**
 * Evaluates a compiled expression against the current BCDTPP state.
 * @param {string} expr
 * @param {!Object} bcdtpp
 * @param {!Object} mesh
 * @return {{h1: string, hcurl: string, hdiv: string, l2: string}}
 */
export function evaluateExpression(expr, bcdtpp, mesh) {
  try {
    const u = compileExpression(expr);
    const pt = mesh.getTetrahedronBarycenter(0);
    const found = bcdtpp.pointLocator.findTetrahedron(pt);
    const tIdx = found ? found.tIdx : 0;

    const h1 = bcdtpp.projectH1(u, pt, tIdx);
    const hcurl = bcdtpp.projectHcurl(u, pt, tIdx);
    const hdiv = bcdtpp.projectHdiv(u, pt, tIdx);
    const l2 = bcdtpp.projectL2(u, tIdx);

    return {
      h1: h1.toFixed(6),
      hcurl: norm(hcurl).toFixed(6),
      hdiv: norm(hdiv).toFixed(6),
      l2: l2.toFixed(6),
    };
  } catch (e) {
    console.error(e);
    return {h1: 'ERR', hcurl: 'ERR', hdiv: 'ERR', l2: 'ERR'};
  }
}

/**
 * Validates uploaded mesh JSON before construction.
 * @param {!Object} data
 * @throws {Error} If the mesh data is malformed.
 */
export function validateMeshData(data) {
  if (!data || typeof data !== 'object') {
    throw new Error('Mesh data must be an object');
  }
  if (!Array.isArray(data.vertices)) {
    throw new Error('Mesh data must contain a "vertices" array');
  }
  if (!Array.isArray(data.tetrahedra)) {
    throw new Error('Mesh data must contain a "tetrahedra" array');
  }

  for (let i = 0; i < data.vertices.length; i++) {
    const v = data.vertices[i];
    if (!Array.isArray(v) || v.length !== 3 || !v.every((n) => typeof n === 'number' && Number.isFinite(n))) {
      throw new Error(`Vertex ${i} must be an array of 3 finite numbers`);
    }
  }

  const vCount = data.vertices.length;
  for (let i = 0; i < data.tetrahedra.length; i++) {
    const t = data.tetrahedra[i];
    if (!Array.isArray(t) || t.length !== 4 || !t.every((n) => typeof n === 'number' && Number.isInteger(n))) {
      throw new Error(`Tetrahedron ${i} must be an array of 4 integer indices`);
    }
    for (let j = 0; j < 4; j++) {
      if (t[j] < 0 || t[j] >= vCount) {
        throw new Error(
          `Tetrahedron ${i} contains out-of-bounds vertex index ${t[j]} (valid range: 0-${vCount - 1})`,
        );
      }
    }
  }
}
