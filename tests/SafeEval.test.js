import {describe, it, expect} from 'vitest';
import {compileExpression} from '../src/lib/safe_eval.js';

describe('safe_eval', () => {
  it('evaluates arithmetic', () => {
    const f = compileExpression('2 + 3 * 4');
    expect(f(0, 0, 0)).toBe(14);
  });

  it('evaluates variables', () => {
    const f = compileExpression('x + y * z');
    expect(f(1, 2, 3)).toBe(7);
  });

  it('evaluates math functions', () => {
    const f = compileExpression('sin(x) + cos(y)');
    expect(f(Math.PI / 2, 0, 0)).toBeCloseTo(2, 6);
  });

  it('treats ^ as exponentiation', () => {
    const f = compileExpression('x^2 + y^3');
    expect(f(3, 2, 0)).toBe(17);
  });

  it('evaluates constants', () => {
    const f = compileExpression('PI + E');
    expect(f(0, 0, 0)).toBeCloseTo(Math.PI + Math.E, 6);
  });

  it('evaluates multi-argument functions', () => {
    const f = compileExpression('min(x, y, z)');
    expect(f(3, 1, 2)).toBe(1);
    const g = compileExpression('atan2(y, x)');
    expect(g(1, 1, 0)).toBeCloseTo(Math.PI / 4, 6);
  });

  it('handles unary minus', () => {
    const f = compileExpression('-x + 5');
    expect(f(3, 0, 0)).toBe(2);
  });

  it('rejects unknown identifiers', () => {
    expect(() => compileExpression('foo(x)')).toThrow('Unknown');
  });

  it('rejects empty expressions', () => {
    expect(() => compileExpression('')).toThrow('Empty');
  });

  it('rejects division by zero', () => {
    const f = compileExpression('1 / x');
    expect(() => f(0, 0, 0)).toThrow('Division by zero');
  });

  it('handles scientific notation', () => {
    const f = compileExpression('1e3 + x');
    expect(f(0, 0, 0)).toBe(1000);
    const g = compileExpression('1e-3 + x');
    expect(g(0, 0, 0)).toBe(0.001);
    const h = compileExpression('.5 + x');
    expect(h(0, 0, 0)).toBe(0.5);
    const i = compileExpression('2.5e+2 + x');
    expect(i(0, 0, 0)).toBe(250);
  });

  it('rejects expressions that are too long', () => {
    const longExpr = '1'.repeat(2049);
    expect(() => compileExpression(longExpr)).toThrow('too long');
  });

  it('rejects expressions with too many tokens', () => {
    const manyTokens = Array(600).fill('1').join('+');
    expect(() => compileExpression(manyTokens)).toThrow('Too many tokens');
  });

  it('rejects deeply nested expressions', () => {
    const deep = '1' + '^1'.repeat(70);
    expect(() => compileExpression(deep)).toThrow('too deep');
  });
});
