import { describe, expect, it } from 'vitest';

import { createSeededRandom, SAMPLER_PRNG_VERSION } from '../../src/sampling/rng';

describe('PRNG scientifique reproductible', () => {
  it('reproduit le vecteur entier xorshift32 de référence pour la seed 1', () => {
    const rng = createSeededRandom(1);
    expect(Array.from({ length: 5 }, () => rng.nextUint32())).toEqual([
      270_369, 67_634_689, 2_647_435_461, 307_599_695, 2_398_689_233,
    ]);
    expect(SAMPLER_PRNG_VERSION).toBe('xorshift32-marsaglia-seed-v1');
  });

  it('documente le remplacement de l’état nul et conserve les seeds uint32', () => {
    expect(createSeededRandom(0)).toMatchObject({
      normalizedSeed: 0,
      initialState: 0x6d2b_79f5,
    });
    expect(createSeededRandom(0xffff_ffff)).toMatchObject({
      normalizedSeed: 0xffff_ffff,
      initialState: 0xffff_ffff,
    });
  });

  it('produit deux suites identiques pour une même seed et distinctes pour deux seeds', () => {
    const first = createSeededRandom(0x1234_5678);
    const second = createSeededRandom(0x1234_5678);
    const other = createSeededRandom(0x1234_5679);
    const firstSequence = Array.from({ length: 1_024 }, () => first.nextUint32());
    const secondSequence = Array.from({ length: 1_024 }, () => second.nextUint32());
    const otherSequence = Array.from({ length: 1_024 }, () => other.nextUint32());
    expect(secondSequence).toEqual(firstSequence);
    expect(otherSequence).not.toEqual(firstSequence);
  });

  it('maintient nextFloat dans [0,1) sans Math.random', () => {
    const rng = createSeededRandom(42);
    for (let index = 0; index < 10_000; index += 1) {
      const value = rng.nextFloat();
      expect(value).toBeGreaterThanOrEqual(0);
      expect(value).toBeLessThan(1);
    }
  });

  it.each([-1, 0.5, Number.NaN, Number.POSITIVE_INFINITY, 0x1_0000_0000])(
    'rejette la seed hors uint32 %s',
    (seed) => {
      expect(() => createSeededRandom(seed)).toThrow(RangeError);
    },
  );
});
