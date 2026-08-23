import { describe, expect, it } from 'vitest';

import { generalizedLaguerre } from '../../src/science/special/generalizedLaguerre';
import { LOW_DEGREE_POLYNOMIAL_ABSOLUTE_TOLERANCE } from './numericAssertions';

describe('polynômes de Laguerre généralisés DLMF', () => {
  it('respecte L_0^alpha(x) = 1 et L_1^alpha(x) = 1 + alpha - x', () => {
    expect(generalizedLaguerre(0, 2.5, 7)).toBe(1);
    expect(generalizedLaguerre(1, 2.5, 0.75)).toBe(2.75);
  });

  it('respecte L_2^0(x) = 1 - 2x + x²/2', () => {
    const x = 0.4;
    const expected = 1 - 2 * x + (x * x) / 2;
    expect(Math.abs(generalizedLaguerre(2, 0, x) - expected)).toBeLessThanOrEqual(
      LOW_DEGREE_POLYNOMIAL_ABSOLUTE_TOLERANCE,
    );
  });

  it('respecte des cas analytiques supplémentaires de la convention DLMF', () => {
    const x = 0.3;
    const l3 = 1 - 3 * x + (3 * x * x) / 2 - (x * x * x) / 6;
    expect(Math.abs(generalizedLaguerre(3, 0, x) - l3)).toBeLessThanOrEqual(
      LOW_DEGREE_POLYNOMIAL_ABSOLUTE_TOLERANCE,
    );
    expect(generalizedLaguerre(4, 2, 0)).toBe(15);
  });

  it.each([
    [-1, 0, 0],
    [1.5, 0, 0],
    [1, -1, 0],
    [1, Number.NaN, 0],
    [1, 0, Number.POSITIVE_INFINITY],
  ])('rejette le triplet hors domaine (%s, %s, %s)', (degree, alpha, x) => {
    expect(() => generalizedLaguerre(degree, alpha, x)).toThrow(RangeError);
  });

  it('échoue explicitement lorsqu’une récurrence dépasse le domaine fini de Number', () => {
    expect(() => generalizedLaguerre(2, 0, Number.MAX_VALUE)).toThrow(RangeError);
  });
});
