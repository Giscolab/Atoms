/** Budget de faux rejet de la campagne statistique déterministe : un sur un million. */
export const SAMPLING_TEST_FAMILY_ALPHA = 1e-6;

function requirePositiveSafeInteger(value: number, quantity: string): void {
  if (!Number.isSafeInteger(value) || value <= 0) {
    throw new RangeError(`${quantity} doit être un entier sûr strictement positif.`);
  }
}

/**
 * Borne uniforme Dvoretzky–Kiefer–Wolfowitz avec constante optimale de Massart,
 * corrigée par Bonferroni pour `comparisonCount` comparaisons :
 * epsilon = sqrt(log(2K/alpha)/(2N)).
 *
 * @see Pascal Massart, Annals of Probability 18(3), 1990.
 * https://doi.org/10.1214/aop/1176990746
 */
export function dkwMassartBound(
  sampleCount: number,
  comparisonCount: number,
  alpha = SAMPLING_TEST_FAMILY_ALPHA,
): number {
  requirePositiveSafeInteger(sampleCount, "Le nombre d'échantillons");
  requirePositiveSafeInteger(comparisonCount, 'Le nombre de comparaisons');
  if (!Number.isFinite(alpha) || alpha <= 0 || alpha >= 1) {
    throw new RangeError('Le risque alpha doit être fini et strictement compris entre zéro et un.');
  }
  return Math.sqrt(Math.log((2 * comparisonCount) / alpha) / (2 * sampleCount));
}

/** Statistique bilatérale exacte de Kolmogorov entre une ECDF et une CDF continue. */
export function kolmogorovDistance(
  samples: readonly number[] | Float32Array | Float64Array,
  theoreticalCdf: (value: number) => number,
): number {
  if (samples.length === 0) throw new RangeError("L'ECDF requiert au moins un échantillon.");
  const sorted = Array.from(samples).sort((left, right) => left - right);
  let maximumDistance = 0;
  for (let index = 0; index < sorted.length; index += 1) {
    const sample = sorted[index];
    if (sample === undefined || !Number.isFinite(sample)) {
      throw new RangeError("Les échantillons de l'ECDF doivent être finis.");
    }
    const probability = theoreticalCdf(sample);
    if (!Number.isFinite(probability) || probability < 0 || probability > 1) {
      throw new RangeError('La CDF théorique doit retourner une probabilité dans [0,1].');
    }
    const lowerEmpirical = index / sorted.length;
    const upperEmpirical = (index + 1) / sorted.length;
    maximumDistance = Math.max(
      maximumDistance,
      Math.abs(probability - lowerEmpirical),
      Math.abs(upperEmpirical - probability),
    );
  }
  return maximumDistance;
}

export function typedArrayValue(
  values: Float32Array | Float64Array,
  index: number,
  quantity: string,
): number {
  const value = values[index];
  if (value === undefined) throw new RangeError(`${quantity} manque à l'indice ${index}.`);
  return value;
}
