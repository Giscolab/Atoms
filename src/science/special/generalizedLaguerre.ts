function validateInputs(degree: number, alpha: number, x: number): void {
  if (!Number.isSafeInteger(degree) || degree < 0) {
    throw new RangeError(`Le degré k doit être un entier sûr non négatif : ${degree}.`);
  }
  if (!Number.isFinite(alpha) || alpha <= -1) {
    throw new RangeError(
      `Le paramètre alpha doit être fini et strictement supérieur à -1 : ${alpha}.`,
    );
  }
  if (!Number.isFinite(x)) {
    throw new RangeError(`L'argument x doit être fini : ${x}.`);
  }
}

function requireFinitePolynomial(value: number, degree: number): number {
  if (!Number.isFinite(value)) {
    throw new RangeError(`L_${degree}^(alpha)(x) dépasse le domaine fini de Number.`);
  }
  return value;
}

/**
 * Polynôme de Laguerre généralisé L_k^(alpha)(x), convention DLMF.
 *
 * Domaine : k entier >= 0, alpha réel > -1 et x réel fini. Le cas hydrogénoïde futur
 * utilisera alpha = 2l + 1. La récurrence est celle de DLMF 18.9.1, table 18.9.1.
 */
export function generalizedLaguerre(degree: number, alpha: number, x: number): number {
  validateInputs(degree, alpha, x);
  if (degree === 0) return 1;

  let previous = 1;
  let current = requireFinitePolynomial(1 + alpha - x, 1);

  for (let currentDegree = 1; currentDegree < degree; currentDegree += 1) {
    const next =
      ((2 * currentDegree + alpha + 1 - x) * current - (currentDegree + alpha) * previous) /
      (currentDegree + 1);
    previous = current;
    current = requireFinitePolynomial(next, currentDegree + 1);
  }

  return current;
}
