/**
 * xorshift32 de Marsaglia, avec opérations exactement définies sur uint32.
 *
 * Contrat v1 de seed : la seed utilisateur doit être un entier de [0, 2³²-1].
 * Elle devient directement l'état uint32, sauf zéro (état absorbant), remplacé
 * par 0x6d2b79f5. `nextFloat()` divise la sortie uint32 par 2³² et retourne [0,1).
 * Ce générateur vise la reproductibilité numérique, pas la cryptographie.
 *
 * @see George Marsaglia, "Xorshift RNGs", Journal of Statistical Software 8(14), 2003.
 * https://doi.org/10.18637/jss.v008.i14
 */
export const SAMPLER_PRNG_VERSION = 'xorshift32-marsaglia-seed-v1' as const;

const UINT32_RANGE = 0x1_0000_0000;
const ZERO_SEED_REPLACEMENT = 0x6d2b_79f5;
const MAX_UINT32 = 0xffff_ffff;

export interface SeededRandomGenerator {
  readonly initialState: number;
  readonly normalizedSeed: number;
  nextFloat(): number;
  nextUint32(): number;
}

function requireUint32Seed(seed: number): number {
  if (!Number.isSafeInteger(seed) || seed < 0 || seed > MAX_UINT32) {
    throw new RangeError(`La seed doit être un entier non signé sur 32 bits : ${seed}.`);
  }
  return seed >>> 0;
}

export function createSeededRandom(seed: number): SeededRandomGenerator {
  const normalizedSeed = requireUint32Seed(seed);
  const initialState = normalizedSeed === 0 ? ZERO_SEED_REPLACEMENT : normalizedSeed;
  let state = initialState;

  const nextUint32 = (): number => {
    let next = state;
    next ^= next << 13;
    next ^= next >>> 17;
    next ^= next << 5;
    state = next >>> 0;
    return state;
  };

  return {
    initialState,
    normalizedSeed,
    nextFloat(): number {
      return nextUint32() / UINT32_RANGE;
    },
    nextUint32,
  };
}
