/**
 * Safe mathematical expression evaluator.
 *
 * Parses expressions containing numbers, variables (x, y, z),
 * standard math functions, and operators (+, -, *, /, ^).
 *
 * No use of eval or new Function.
 */

const OPS = {
  PLUS: 'PLUS',
  MINUS: 'MINUS',
  MUL: 'MUL',
  DIV: 'DIV',
  POW: 'POW',
  LPAREN: 'LPAREN',
  RPAREN: 'RPAREN',
  COMMA: 'COMMA',
  NUMBER: 'NUMBER',
  IDENT: 'IDENT',
  EOF: 'EOF',
};

const MATH_FUNCS = new Set([
  'sin', 'cos', 'tan', 'asin', 'acos', 'atan', 'atan2',
  'exp', 'log', 'log10', 'log2', 'sqrt', 'cbrt',
  'abs', 'ceil', 'floor', 'round', 'trunc',
  'pow', 'min', 'max', 'sign',
]);

const MATH_CONSTS = {
  PI: Math.PI,
  E: Math.E,
};

const VARIABLES = new Set(['x', 'y', 'z']);

function tokenize(expr) {
  const tokens = [];
  let i = 0;
  while (i < expr.length) {
    const c = expr[i];
    if (/\s/.test(c)) {
      i++;
      continue;
    }
    if (/[0-9]/.test(c) || (c === '.' && /[0-9]/.test(expr[i + 1]))) {
      let j = i;
      while (j < expr.length && (/[0-9]/.test(expr[j]) || expr[j] === '.')) j++;
      tokens.push({type: OPS.NUMBER, value: parseFloat(expr.slice(i, j))});
      i = j;
      continue;
    }
    if (/[a-zA-Z_]/.test(c)) {
      let j = i;
      while (j < expr.length && /[a-zA-Z0-9_]/.test(expr[j])) j++;
      tokens.push({type: OPS.IDENT, value: expr.slice(i, j)});
      i = j;
      continue;
    }
    switch (c) {
      case '+': tokens.push({type: OPS.PLUS}); break;
      case '-': tokens.push({type: OPS.MINUS}); break;
      case '*': tokens.push({type: OPS.MUL}); break;
      case '/': tokens.push({type: OPS.DIV}); break;
      case '^': tokens.push({type: OPS.POW}); break;
      case '(': tokens.push({type: OPS.LPAREN}); break;
      case ')': tokens.push({type: OPS.RPAREN}); break;
      case ',': tokens.push({type: OPS.COMMA}); break;
      default:
        throw new Error(`Unexpected character '${c}' at position ${i}`);
    }
    i++;
  }
  tokens.push({type: OPS.EOF});
  return tokens;
}

function parse(tokens) {
  let pos = 0;
  const peek = () => tokens[pos];
  const advance = () => tokens[pos++];
  const expect = (type) => {
    if (peek().type !== type) {
      throw new Error(`Expected ${type} but found ${peek().type}`);
    }
    return advance();
  };

  function parseExpression() {
    return parseAddSub();
  }

  function parseAddSub() {
    let node = parseMulDiv();
    while (true) {
      const t = peek();
      if (t.type === OPS.PLUS) {
        advance();
        node = {op: 'add', left: node, right: parseMulDiv()};
      } else if (t.type === OPS.MINUS) {
        advance();
        node = {op: 'sub', left: node, right: parseMulDiv()};
      } else {
        break;
      }
    }
    return node;
  }

  function parseMulDiv() {
    let node = parsePower();
    while (true) {
      const t = peek();
      if (t.type === OPS.MUL) {
        advance();
        node = {op: 'mul', left: node, right: parsePower()};
      } else if (t.type === OPS.DIV) {
        advance();
        node = {op: 'div', left: node, right: parsePower()};
      } else {
        break;
      }
    }
    return node;
  }

  function parsePower() {
    let node = parseUnary();
    const t = peek();
    if (t.type === OPS.POW) {
      advance();
      node = {op: 'pow', left: node, right: parsePower()};
    }
    return node;
  }

  function parseUnary() {
    const t = peek();
    if (t.type === OPS.PLUS) {
      advance();
      return parseUnary();
    }
    if (t.type === OPS.MINUS) {
      advance();
      return {op: 'neg', arg: parseUnary()};
    }
    return parsePrimary();
  }

  function parsePrimary() {
    const t = peek();
    if (t.type === OPS.NUMBER) {
      advance();
      return {op: 'num', value: t.value};
    }
    if (t.type === OPS.IDENT) {
      advance();
      const name = t.value;
      if (peek().type === OPS.LPAREN) {
        advance();
        const args = [];
        if (peek().type !== OPS.RPAREN) {
          args.push(parseExpression());
          while (peek().type === OPS.COMMA) {
            advance();
            args.push(parseExpression());
          }
        }
        expect(OPS.RPAREN);
        if (!MATH_FUNCS.has(name)) {
          throw new Error(`Unknown function '${name}'`);
        }
        return {op: 'call', name, args};
      }
      if (VARIABLES.has(name)) {
        return {op: 'var', name};
      }
      if (name in MATH_CONSTS) {
        return {op: 'num', value: MATH_CONSTS[name]};
      }
      throw new Error(`Unknown identifier '${name}'`);
    }
    if (t.type === OPS.LPAREN) {
      advance();
      const node = parseExpression();
      expect(OPS.RPAREN);
      return node;
    }
    throw new Error(`Unexpected token ${t.type}`);
  }

  const ast = parseExpression();
  if (peek().type !== OPS.EOF) {
    throw new Error(`Unexpected token ${peek().type} after end of expression`);
  }
  return ast;
}

function evaluate(ast, vars) {
  switch (ast.op) {
    case 'num': return ast.value;
    case 'var': {
      const v = vars[ast.name];
      if (v === undefined) throw new Error(`Undefined variable '${ast.name}'`);
      return v;
    }
    case 'neg': return -evaluate(ast.arg, vars);
    case 'add': return evaluate(ast.left, vars) + evaluate(ast.right, vars);
    case 'sub': return evaluate(ast.left, vars) - evaluate(ast.right, vars);
    case 'mul': return evaluate(ast.left, vars) * evaluate(ast.right, vars);
    case 'div': {
      const denom = evaluate(ast.right, vars);
      if (denom === 0) throw new Error('Division by zero');
      return evaluate(ast.left, vars) / denom;
    }
    case 'pow': return Math.pow(evaluate(ast.left, vars), evaluate(ast.right, vars));
    case 'call': {
      const args = ast.args.map((a) => evaluate(a, vars));
      switch (ast.name) {
        case 'sin': return Math.sin(args[0]);
        case 'cos': return Math.cos(args[0]);
        case 'tan': return Math.tan(args[0]);
        case 'asin': return Math.asin(args[0]);
        case 'acos': return Math.acos(args[0]);
        case 'atan': return args.length === 2 ? Math.atan2(args[0], args[1]) : Math.atan(args[0]);
        case 'atan2': return Math.atan2(args[0], args[1]);
        case 'exp': return Math.exp(args[0]);
        case 'log': return Math.log(args[0]);
        case 'log10': return Math.log10(args[0]);
        case 'log2': return Math.log2(args[0]);
        case 'sqrt': return Math.sqrt(args[0]);
        case 'cbrt': return Math.cbrt(args[0]);
        case 'abs': return Math.abs(args[0]);
        case 'ceil': return Math.ceil(args[0]);
        case 'floor': return Math.floor(args[0]);
        case 'round': return Math.round(args[0]);
        case 'trunc': return Math.trunc(args[0]);
        case 'pow': return Math.pow(args[0], args[1]);
        case 'min': return Math.min(...args);
        case 'max': return Math.max(...args);
        case 'sign': return Math.sign(args[0]);
        default:
          throw new Error(`Unknown function '${ast.name}'`);
      }
    }
    default:
      throw new Error(`Unknown AST node ${ast.op}`);
  }
}

/**
 * Compiles a math expression string into an evaluable function.
 * @param {string} expr
 * @return {function(number, number, number): number}
 */
export function compileExpression(expr) {
  if (!expr || expr.trim() === '') {
    throw new Error('Empty expression');
  }
  const tokens = tokenize(expr);
  const ast = parse(tokens);
  return (x, y, z) => evaluate(ast, {x, y, z});
}
