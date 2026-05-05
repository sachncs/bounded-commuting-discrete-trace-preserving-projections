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
});
