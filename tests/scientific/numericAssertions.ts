/** Marge pour une courte chaîne d'opérations IEEE-754 en double précision. */
export const SHORT_FLOAT64_RELATIVE_TOLERANCE = 16 * Number.EPSILON;

/** Marge absolue pour les polynômes analytiques de faible degré testés dans ce lot. */
export const LOW_DEGREE_POLYNOMIAL_ABSOLUTE_TOLERANCE = 64 * Number.EPSILON;

/**
 * Marge composante par composante pour les formes fermées de faible degré.
 * Le facteur couvre les racines carrées et fonctions trigonométriques qui
 * s'ajoutent aux courtes chaînes d'opérations déjà testées au Lot 2.
 */
export const ANALYTIC_COMPONENT_TOLERANCE = 512 * Number.EPSILON;

/**
 * Marge absolue des quadratures radiales et angulaires composites. Elle est
 * supérieure à l'erreur observée entre les maillages convergés documentés par
 * les tests, tout en restant plusieurs ordres sous les écarts physiques testés.
 */
export const NUMERICAL_NORMALIZATION_ABSOLUTE_TOLERANCE = 1e-8;

/** Une orthogonalité vise zéro et requiert donc une marge absolue distincte. */
export const NUMERICAL_ORTHOGONALITY_ABSOLUTE_TOLERANCE = 1e-8;

/** Marge de la véritable quadrature sphérique 3D, volontairement plus coûteuse. */
export const THREE_DIMENSIONAL_NORMALIZATION_ABSOLUTE_TOLERANCE = 1e-7;

/** Stabilité exigée entre les maillages 3D grossier et fin utilisés en validation. */
export const THREE_DIMENSIONAL_CONVERGENCE_ABSOLUTE_TOLERANCE = 5e-7;

export function relativeError(actual: number, expected: number): number {
  if (expected === 0) {
    throw new RangeError("L'erreur relative requiert une valeur attendue non nulle.");
  }
  return Math.abs((actual - expected) / expected);
}

/** Distance circulaire minimale entre deux phases exprimées en radians. */
export function phaseDistance(actual: number, expected: number): number {
  if (!Number.isFinite(actual) || !Number.isFinite(expected)) {
    throw new RangeError('La comparaison de phase requiert deux angles finis.');
  }
  return Math.abs(Math.atan2(Math.sin(actual - expected), Math.cos(actual - expected)));
}
