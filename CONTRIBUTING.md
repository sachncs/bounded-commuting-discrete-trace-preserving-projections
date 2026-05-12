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

## Pull Request Process

1. Create a feature branch from `master`:
   ```bash
   git checkout -b feature/my-change
   ```
2. Make your changes, add tests, and ensure `npm test` passes.
3. Run `npm run lint:fix` to auto-format.
4. Commit with a clear message describing the change.
5. Push to your fork and open a Pull Request against `master`.
6. Fill out the PR template with a summary and test plan.

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
