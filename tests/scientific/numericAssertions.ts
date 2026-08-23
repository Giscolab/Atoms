/** Marge pour une courte chaîne d'opérations IEEE-754 en double précision. */
export const SHORT_FLOAT64_RELATIVE_TOLERANCE = 16 * Number.EPSILON;

/** Marge absolue pour les polynômes analytiques de faible degré testés dans ce lot. */
export const LOW_DEGREE_POLYNOMIAL_ABSOLUTE_TOLERANCE = 64 * Number.EPSILON;

export function relativeError(actual: number, expected: number): number {
  if (expected === 0) {
    throw new RangeError("L'erreur relative requiert une valeur attendue non nulle.");
  }
  return Math.abs((actual - expected) / expected);
}
