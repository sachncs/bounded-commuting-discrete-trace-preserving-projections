# Publishing

## Versioning

This project follows [Semantic Versioning](https://semver.org/):

- **MAJOR**: Breaking API changes
- **MINOR**: New features, backwards-compatible
- **PATCH**: Bug fixes, backwards-compatible

## Manual Publish

1. Ensure you are on `master` and the working tree is clean.
2. Run the full verification suite:
   ```bash
   npm run lint
   npm test
   npm run build
   npm run docs
   ```
3. Bump the version:
   ```bash
   npm version patch   # or minor / major
   ```
4. Push the tag:
   ```bash
   git push && git push --tags
   ```
5. Publish to npm:
   ```bash
   npm publish --access public
   ```

## Automated Publish (CI)

The repository includes a GitHub Actions workflow (`.github/workflows/publish.yml`) that triggers on release creation.

1. Create a release on GitHub with a new tag (e.g., `v0.1.1`).
2. The workflow will:
   - Run lint and tests
   - Build bundles
   - Generate docs
   - Publish to npm using `NPM_TOKEN` secret

Ensure the `NPM_TOKEN` secret is configured in the repository settings before using automated publish.
