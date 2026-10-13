/* Reproducible practice utilities. No network requests or external dependencies. */
(function () {
  'use strict';
  globalThis.ChemGenerators = Object.create(null);
  function rng(seed) {
    let a = seed >>> 0;
    const random = () => {
      a = (a + 0x6D2B79F5) >>> 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t ^= t + Math.imul(t ^ (t >>> 7), 61 | t);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
    return { random, int: (min, max) => Math.floor(random() * (max - min + 1)) + min, pick: xs => xs[Math.floor(random() * xs.length)] };
  }
  function fmt(x, sig = 4) {
    if (!Number.isFinite(x)) throw new Error('Nonfinite chemistry value');
    if (x === 0) return '0';
    return Number(sigfig(x, sig)).toString();
  }
  function sigfig(x, sig = 3) {
    if (!Number.isFinite(x)) throw new Error('Nonfinite chemistry value');
    if (x === 0) return x.toPrecision(sig);
    // Guard decimal halfway values against a few ulps of binary arithmetic
    // noise, then use half-up rounding on the magnitude and retain zeros.
    const [mantissa, power] = Math.abs(x).toExponential().split('e');
    const scaled = Number(mantissa) * 10 ** (sig - 1);
    const digits = Math.round(scaled + 4 * Number.EPSILON * scaled);
    const rounded = Number(`${digits}e${Number(power) - sig + 1}`) * Math.sign(x);
    return rounded.toPrecision(sig);
  }
  function parseNumber(raw) {
    if (typeof raw !== 'string' || !raw.trim()) return null;
    let text = raw.trim().replace(/−/g, '-');
    const scientific = text.match(/^([+-]?(?:\d+\.?\d*|\.\d+))\s*[×xX*]\s*10\s*\^\s*([+-]?\d+)$/);
    if (scientific) text = scientific[1] + 'e' + scientific[2];
    if (!/^[+-]?(?:\d+\.?\d*|\.\d+)(?:e[+-]?\d+)?$/i.test(text)) return null;
    const value = Number(text);
    // Do not silently turn a nonzero value outside JS's range into zero.
    if (value === 0 && /[1-9]/.test(text.split(/[eE]/)[0])) return null;
    return Number.isFinite(value) ? value : null;
  }
  function grade(field, raw) {
    if (field.type === 'select') {
      const valid = field.options.some(x => x.value === raw);
      return { valid, correct: valid && raw === field.answer, value: raw };
    }
    const value = parseNumber(raw);
    if (value === null) return { valid: false, correct: false, value: null };
    // Relative tolerance matters for tiny Ksp values; no large absolute floor.
    const tolerance = field.tolerance ?? 0.01;
    // Exact counts must not gain a floating-point epsilon allowance.
    const correct = tolerance === 0 || field.answer === 0
      ? value === field.answer
      : Math.abs(value - field.answer) <= Math.abs(field.answer) * tolerance + Number.EPSILON * Math.abs(field.answer);
    return { valid: true, correct, value };
  }
  globalThis.ChemCore = {
    version: '1', rng, fmt, sigfig, parseNumber, grade,
    num: (id, label, answer, unit = '') => ({ id, label, type: 'number', answer, unit, tolerance: 0.01 }),
    choice: (id, label, answer, options) => ({ id, label, type: 'select', answer, options })
  };
})();
