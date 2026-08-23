import { describe, expect, it } from 'vitest';

import {
  assertValidQuantumNumbers,
  isValidQuantumNumbers,
} from '../../src/science/quantum/quantumNumbers';

const INVALID_QUANTUM_NUMBERS: readonly unknown[] = [
  null,
  {},
  { n: '2', l: 1, m: 0 },
  { n: 0, l: 0, m: 0 },
  { n: 1.5, l: 0, m: 0 },
  { n: 2, l: -1, m: 0 },
  { n: 3, l: 1.5, m: 0 },
  { n: 2, l: 2, m: 0 },
  { n: 3, l: 1, m: -2 },
  { n: 3, l: 1, m: 2 },
  { n: 3, l: 1, m: 0.5 },
  { n: 3, l: Number.NaN, m: 0 },
  { n: Number.MAX_SAFE_INTEGER + 1, l: 0, m: 0 },
  { n: Number.POSITIVE_INFINITY, l: 0, m: 0 },
];

describe('nombres quantiques hydrogénoïdes', () => {
  it('accepte les bornes physiques n >= 1, 0 <= l <= n - 1 et -l <= m <= l', () => {
    expect(isValidQuantumNumbers({ n: 1, l: 0, m: 0 })).toBe(true);
    expect(isValidQuantumNumbers({ n: 4, l: 3, m: -3 })).toBe(true);
    expect(isValidQuantumNumbers({ n: 4, l: 3, m: 3 })).toBe(true);
    expect(() => {
      assertValidQuantumNumbers({ n: 5, l: 2, m: 1 });
    }).not.toThrow();
  });

  it.each(INVALID_QUANTUM_NUMBERS)(
    'rejette explicitement l’entrée invalide %j',
    (quantumNumbers) => {
      expect(isValidQuantumNumbers(quantumNumbers)).toBe(false);
      expect(() => {
        assertValidQuantumNumbers(quantumNumbers);
      }).toThrow(RangeError);
    },
  );
});
