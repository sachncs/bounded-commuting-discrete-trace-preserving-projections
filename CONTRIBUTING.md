# Contributing to Bcdtpp

Thank you for your interest in contributing! This document outlines the workflow for proposing changes, reporting issues, and submitting pull requests.

## Getting Started

1. Fork the repository on GitHub.
2. Clone your fork locally.
3. Install dependencies:
   ```bash
   npm install
   ```

## Development Workflow

We use a zero-config toolchain:

```bash
npm run lint        # Check code style with StandardJS
npm run lint:fix    # Auto-fix style issues
npm test            # Run the full Mocha test suite
npm run test:coverage  # Run tests with c8 coverage report
npm run build       # Build Rollup bundles (ESM, CJS, UMD)
npm run docs        # Regenerate API documentation
```

All contributions must pass lint and tests before being merged.

## Commit Conventions

We follow [Conventional Commits](https://www.conventionalcommits.org/) for commit messages. This enables automatic changelog generation and semantic versioning.

### Format

```
<type>(<scope>): <description>

[optional body]

[optional footer(s)]
```

### Types

- **feat**: A new feature
- **fix**: A bug fix
- **docs**: Documentation only changes
- **style**: Code style changes (formatting, missing semi-colons, etc.)
- **refactor**: Code change that neither fixes a bug nor adds a feature
- **perf**: Performance improvement
- **test**: Adding missing tests or correcting existing tests
- **chore**: Build process or auxiliary tool changes
- **ci**: CI configuration changes
- **revert**: Revert a previous commit

### Examples

```bash
git commit -m "feat(mesh): add anisotropic mesh support"
git commit -m "fix(point-locator): correct barycentric coordinate calculation"
git commit -m "docs: update API reference for Bcdtpp class"
git commit -m "test(bcdtpp): add edge cases for degenerate tetrahedra"
git commit -m "chore: update rollup to v4"
```

### Scope

Optional scope should reference the module affected:

- `mesh` - Mesh topology and geometry
- `whitney` - Whitney forms and barycentric utilities
- `bcdtpp` - Main projection class
- `point-locator` - AABB tree point location
- `quadrature` - Gaussian quadrature
- `math-utils` - Linear algebra primitives
- `local-solver` - Boundary-patch assembly
- `higher-order` - Higher-order projection
- `docs` - Documentation

## Pull Request Process

1. Create a feature branch from `master`:
   ```bash
   git checkout -b feat/my-change
   ```
2. Make your changes, add tests, and ensure `npm test` passes.
3. Run `npm run lint:fix` to auto-format.
4. Commit with a conventional commit message describing the change.
5. Push to your fork and open a Pull Request against `master`.
6. Fill out the PR template with a summary and test plan.
7. Ensure CI passes and request a review from a maintainer.

## Reporting Bugs

Please use the [Bug Report issue template](.github/ISSUE_TEMPLATE/bug_report.md) and include:
- A minimal reproduction case or code snippet
- Expected vs actual behavior
- Node.js version and OS

## Requesting Features

Please use the [Feature Request issue template](.github/ISSUE_TEMPLATE/feature_request.md) and describe:
- The problem or limitation you are facing
- The proposed solution or API change
- Any alternatives you have considered

## Code of Conduct

This project adheres to the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code.
