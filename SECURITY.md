# Security Policy

## Supported Versions

| Version | Supported          |
| ------- | ------------------ |
| 0.1.x   | :white_check_mark: |

## Reporting a Vulnerability

If you discover a security vulnerability in this project, please report it responsibly.

1. **Do not** open a public issue describing the vulnerability.
2. Instead, contact the maintainer directly via the email associated with the repository owner.
3. Provide a clear description of the issue, steps to reproduce, and potential impact.
4. Allow reasonable time for the issue to be addressed before disclosing it publicly.

We take all security reports seriously and will respond as quickly as possible.

## Security Best Practices for Consumers

- This library performs numerical computation on user-provided mesh data. Always validate inputs before passing them to `Mesh` or `Bcdtpp` constructors.
- Do not execute untrusted function inputs directly through projection methods without sandboxing.
