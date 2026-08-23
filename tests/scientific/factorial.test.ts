import { describe, expect, it } from 'vitest';

import { factorial } from '../../src/science/special/factorial';

describe('factorielle', () => {
  it('respecte les cas de base et des valeurs entières vérifiables exactement', () => {
    expect(factorial(0)).toBe(1);
    expect(factorial(1)).toBe(1);
    expect(factorial(5)).toBe(120);
    expect(factorial(10)).toBe(3_628_800);
  });

  it.each([-1, 0.5, Number.NaN, Number.POSITIVE_INFINITY])(
    'rejette l’entrée invalide %s',
    (value) => {
      expect(() => factorial(value)).toThrow(RangeError);
    },
  );

  it('échoue explicitement lorsque le résultat dépasserait le domaine fini de Number', () => {
    expect(Number.isFinite(factorial(170))).toBe(true);
    expect(() => factorial(171)).toThrow(RangeError);
  });
});
