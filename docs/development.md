# Development Workflow

This repository follows a zero-config philosophy: every script is a single `npm run <command>` away.

## Scripts

| Script | Command | Purpose |
|--------|---------|---------|
| Lint | `npm run lint` | Check StandardJS style on source and tests |
| Lint Fix | `npm run lint:fix` | Auto-fix StandardJS violations |
| Test | `npm test` | Run the full Mocha test suite (166 tests) |
| Test Watch | `npm run test:watch` | Run tests in watch mode during development |
| Coverage | `npm run test:coverage` | Run tests with c8; outputs text, HTML, and lcov |
| Build | `npm run build` | Build Rollup bundles: ESM, CJS, UMD |
| Docs | `npm run docs` | Regenerate `docs/api.md` from JSDoc comments |
| Full Build | `npm run build:full` | Build + docs in one step |

## Typical PR Workflow

1. Make changes in `src/lib/` or `tests/`.
2. Run `npm run lint:fix` to auto-format.
3. Run `npm test` to verify behavior.
4. Run `npm run build` to ensure bundles compile.
5. Commit and push.

## Adding a New Test

Create a file in `tests/` ending in `.test.js`. Mocha will pick it up automatically.

```javascript
import { expect } from 'chai'
import { Mesh } from '../src/lib/mesh.js'

describe('My new feature', () => {
  it('does something', () => {
    const mesh = new Mesh(...)
    expect(mesh.tetrahedronCount).to.equal(1)
  })
})
```

## Code Style

We use [StandardJS](https://standardjs.com). There is no `.eslintrc` or `.prettierrc`; the style is enforced entirely by the `standard` CLI.

Key rules:
- 2 spaces indentation
- No semicolons
- Single quotes
- No trailing commas in function parameters
- Camel case identifiers
