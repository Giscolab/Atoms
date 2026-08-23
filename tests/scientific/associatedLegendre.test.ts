import { describe, expect, it } from 'vitest';

import { associatedLegendre } from '../../src/science/special/associatedLegendre';
import { LOW_DEGREE_POLYNOMIAL_ABSOLUTE_TOLERANCE } from './numericAssertions';

function expectPolynomialValue(actual: number, expected: number): void {
  expect(Math.abs(actual - expected)).toBeLessThanOrEqual(LOW_DEGREE_POLYNOMIAL_ABSOLUTE_TOLERANCE);
}

describe('fonctions de Ferrers avec phase de Condon–Shortley', () => {
  it('respecte P_0^0(x) = 1 et P_1^0(x) = x', () => {
    expect(associatedLegendre(0, 0, 0.25)).toBe(1);
    expect(associatedLegendre(1, 0, 0.25)).toBe(0.25);
  });

  it('inclut la phase de Condon–Shortley dans P_1^1', () => {
    const x = 0.3;
    expectPolynomialValue(associatedLegendre(1, 1, x), -Math.sqrt(1 - x * x));
  });

  it('respecte plusieurs expressions analytiques de degré deux et trois', () => {
    const x = 0.3;
    const root = Math.sqrt(1 - x * x);
    expectPolynomialValue(associatedLegendre(2, 0, x), (3 * x * x - 1) / 2);
    expectPolynomialValue(associatedLegendre(2, 1, x), -3 * x * root);
    expectPolynomialValue(associatedLegendre(2, 2, x), 3 * (1 - x * x));
    expectPolynomialValue(associatedLegendre(3, 3, x), -15 * root ** 3);
  });

  it('respecte les valeurs aux bornes sans clamp silencieux', () => {
    expect(associatedLegendre(4, 0, 1)).toBe(1);
    expect(associatedLegendre(3, 0, -1)).toBe(-1);
    expect(Math.abs(associatedLegendre(4, 2, 1))).toBe(0);
  });

  it.each([
    [-1, 0, 0],
    [1.5, 0, 0],
    [1, -1, 0],
    [1, 2, 0],
    [1, 0, -1.000_000_000_1],
    [1, 0, Number.NaN],
  ])('rejette le triplet hors domaine (%s, %s, %s)', (degree, order, x) => {
    expect(() => associatedLegendre(degree, order, x)).toThrow(RangeError);
  });
});
