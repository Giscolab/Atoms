import { describe, expect, it } from 'vitest';

import { parseWavefunctionData } from '../../src/data/wavefunctionData';

describe('parseWavefunctionData', () => {
  it('valide de façon déterministe un tableau de triplets numériques finis', () => {
    const input = {
      points: [
        [0, 1, -2],
        [3.5, 0, 4],
      ],
    };

    expect(parseWavefunctionData(input)).toEqual(input);
    expect(parseWavefunctionData(input)).toEqual(input);
  });

  it.each([
    null,
    {},
    { points: 'invalid' },
    { points: [[0, 1]] },
    { points: [[0, 1, Number.POSITIVE_INFINITY]] },
  ])('rejette les données structurellement invalides sans interprétation physique', (input) => {
    expect(() => parseWavefunctionData(input)).toThrow();
  });
});
