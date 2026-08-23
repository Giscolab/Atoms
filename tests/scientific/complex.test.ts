import { describe, expect, it } from 'vitest';

import {
  addComplex,
  complexExponential,
  complexModulusSquared,
  complexPhase,
  conjugateComplex,
  createComplex,
  multiplyComplex,
  scaleComplex,
  subtractComplex,
  type Complex,
} from '../../src/science/math/complex';
import { ANALYTIC_COMPONENT_TOLERANCE, phaseDistance } from './numericAssertions';

function expectComponent(actual: number, expected: number): void {
  expect(Math.abs(actual - expected)).toBeLessThanOrEqual(ANALYTIC_COMPONENT_TOLERANCE);
}

function expectComplexComponents(actual: Complex, expectedRe: number, expectedIm: number): void {
  expectComponent(actual.re, expectedRe);
  expectComponent(actual.im, expectedIm);
}

describe('arithmétique complexe pure', () => {
  it('crée exactement zéro, un et l’unité imaginaire', () => {
    expect(createComplex(0, 0)).toEqual({ re: 0, im: 0 });
    expect(createComplex(1, 0)).toEqual({ re: 1, im: 0 });
    expect(createComplex(0, 1)).toEqual({ re: 0, im: 1 });
  });

  it('additionne et soustrait les deux composantes', () => {
    const left = createComplex(2, -3);
    const right = createComplex(-5, 7);
    expect(addComplex(left, right)).toEqual({ re: -3, im: 4 });
    expect(subtractComplex(left, right)).toEqual({ re: 7, im: -10 });
  });

  it('multiplie deux complexes et applique un scalaire', () => {
    expect(multiplyComplex(createComplex(0, 1), createComplex(0, 1))).toEqual({ re: -1, im: 0 });
    expect(multiplyComplex(createComplex(2, 3), createComplex(-4, 5))).toEqual({
      re: -23,
      im: -2,
    });
    expect(scaleComplex(createComplex(1.5, -4), -2)).toEqual({ re: -3, im: 8 });
  });

  it('conjugue z et vérifie z fois son conjugué ainsi que son module carré', () => {
    const value = createComplex(3, 4);
    const conjugate = conjugateComplex(value);
    expect(conjugate).toEqual({ re: 3, im: -4 });
    expect(multiplyComplex(value, conjugate)).toEqual({ re: 25, im: 0 });
    expect(complexModulusSquared(value)).toBe(25);
    expect(complexModulusSquared(createComplex(0, 0))).toBe(0);
  });

  it('construit exp(i phi) aux angles cardinaux', () => {
    expectComplexComponents(complexExponential(0), 1, 0);
    expectComplexComponents(complexExponential(Math.PI / 2), 0, 1);
    expectComplexComponents(complexExponential(Math.PI), -1, 0);
    expectComplexComponents(complexExponential(-Math.PI / 2), 0, -1);
  });

  it.each([
    [1, 1, Math.PI / 4],
    [-1, 1, (3 * Math.PI) / 4],
    [-1, -1, (-3 * Math.PI) / 4],
    [1, -1, -Math.PI / 4],
  ])('retourne la phase correcte dans chaque quadrant pour (%s, %s)', (re, im, expected) => {
    const actual = complexPhase(createComplex(re, im));
    expect(actual).not.toBeNull();
    if (actual === null) throw new Error('Une amplitude non nulle doit posséder une phase.');
    expect(phaseDistance(actual, expected)).toBeLessThanOrEqual(ANALYTIC_COMPONENT_TOLERANCE);
  });

  it('représente explicitement comme null la phase de l’amplitude exactement nulle', () => {
    expect(complexPhase(createComplex(0, 0))).toBeNull();
    expect(complexPhase(createComplex(-0, 0))).toBeNull();
  });

  it('rejette toute composante, tout angle ou tout scalaire non fini', () => {
    const zero = createComplex(0, 0);
    const invalid = { re: Number.NaN, im: 0 };

    expect(() => createComplex(Number.NaN, 0)).toThrow(RangeError);
    expect(() => createComplex(0, Number.POSITIVE_INFINITY)).toThrow(RangeError);
    expect(() => addComplex(invalid, zero)).toThrow(RangeError);
    expect(() => subtractComplex(zero, invalid)).toThrow(RangeError);
    expect(() => multiplyComplex(invalid, zero)).toThrow(RangeError);
    expect(() => scaleComplex(zero, Number.NEGATIVE_INFINITY)).toThrow(RangeError);
    expect(() => conjugateComplex(invalid)).toThrow(RangeError);
    expect(() => complexModulusSquared(invalid)).toThrow(RangeError);
    expect(() => complexExponential(Number.NaN)).toThrow(RangeError);
    expect(() => complexPhase(invalid)).toThrow(RangeError);
  });

  it('échoue explicitement si une opération finie produirait un résultat non fini', () => {
    const maximum = createComplex(Number.MAX_VALUE, 0);
    expect(() => scaleComplex(maximum, 2)).toThrow(RangeError);
    expect(() => multiplyComplex(maximum, createComplex(2, 0))).toThrow(RangeError);
    expect(() => complexModulusSquared(maximum)).toThrow(RangeError);
  });
});
