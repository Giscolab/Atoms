import { afterEach, describe, expect, it, vi } from 'vitest';

import { createSeededRandom } from '../../src/sampling/rng';
import { randomUint32 } from '../../src/ui/randomSeed';

afterEach(() => {
  vi.restoreAllMocks();
  vi.unstubAllGlobals();
});

describe('seed aléatoire de l’interface', () => {
  it('conserve une sortie Crypto occupant les 32 bits', () => {
    vi.stubGlobal('crypto', {
      getRandomValues(values: Uint32Array): Uint32Array {
        values[0] = 0xffff_ffff;
        return values;
      },
    });

    expect(randomUint32()).toBe(0xffff_ffff);
  });

  it('fournit une seed uint32 acceptée par le sampler quand Crypto est absent', () => {
    vi.stubGlobal('crypto', undefined);
    vi.spyOn(Date, 'now').mockReturnValue(0x8000_0000);

    const seed = randomUint32();

    expect(Number.isSafeInteger(seed)).toBe(true);
    expect(seed).toBeGreaterThanOrEqual(0);
    expect(seed).toBeLessThanOrEqual(0xffff_ffff);
    expect(createSeededRandom(seed).normalizedSeed).toBe(seed);
  });

  it('reste compatible avec le sampler quand Crypto lève une erreur', () => {
    vi.stubGlobal('crypto', {
      getRandomValues(): never {
        throw new Error('Crypto indisponible');
      },
    });
    vi.spyOn(Date, 'now').mockReturnValue(0xffff_ffff);

    const seed = randomUint32();

    expect(createSeededRandom(seed).normalizedSeed).toBe(seed);
  });
});
