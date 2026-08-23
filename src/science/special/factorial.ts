/** Factorielle des entiers non négatifs dans le domaine fini de Number. */
export function factorial(value: number): number {
  if (!Number.isSafeInteger(value) || value < 0) {
    throw new RangeError(`La factorielle requiert un entier sûr non négatif : ${value}.`);
  }

  let result = 1;
  for (let factor = 2; factor <= value; factor += 1) {
    result *= factor;
    if (!Number.isFinite(result)) {
      throw new RangeError(`La factorielle de ${value} dépasse le domaine fini de Number.`);
    }
  }
  return result;
}
