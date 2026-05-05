/**
 * UI module for the BCDTPP web application.
 *
 * **Coupling note:** This module is tightly coupled to a specific HTML page
 * structure.  It references hardcoded DOM element IDs:
 *   - 'func-input'    – expression input field
 *   - 'mesh-upload'   – mesh file upload input
 *   - 'mesh-status'   – status message display
 *   - 'val-h1'        – H1 result display
 *   - 'val-hcurl'     – H(curl) result display
 *   - 'val-hdiv'      – H(div) result display
 *   - 'val-l2'        – L2 result display
 *
 * These IDs must exist in the DOM; the module is not reusable outside the
 * bundled index.html without modification.
 */

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
 *
 * Errors are swallowed and replaced with 'ERR' strings to avoid leaking
 * internal details to the UI.  If you need diagnostics, pass an onError
 * callback.
 *
 * @param {string} expr
 * @param {!Object} bcdtpp
 * @param {!Object} mesh
 * @param {function(Error)=} onError - Optional error handler.
 * @return {{h1: string, hcurl: string, hdiv: string, l2: string}}
 */
export function evaluateExpression(expr, bcdtpp, mesh, onError) {
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
    if (onError) onError(e);
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

  // Check for coincident vertices (duplicates within 1e-6).
  const vertexKey = (v) =>
    `${Math.round(v[0] * 1e6)},${Math.round(v[1] * 1e6)},${Math.round(v[2] * 1e6)}`;
  const seenVerts = new Set();
  for (let i = 0; i < data.vertices.length; i++) {
    const key = vertexKey(data.vertices[i]);
    if (seenVerts.has(key)) {
      throw new Error(`Vertices ${i} is coincident with an earlier vertex`);
    }
    seenVerts.add(key);
  }

  const vCount = data.vertices.length;
  const seenTets = new Set();
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
    if (new Set(t).size !== 4) {
      throw new Error(`Tetrahedron ${i} contains duplicate vertex indices`);
    }
    const sorted = [...t].sort((a, b) => a - b);
    const tKey = `${sorted[0]},${sorted[1]},${sorted[2]},${sorted[3]}`;
    if (seenTets.has(tKey)) {
      throw new Error(`Tetrahedron ${i} is a duplicate of an earlier tetrahedron`);
    }
    seenTets.add(tKey);
  }
}
