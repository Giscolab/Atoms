function validateInputs(degree: number, order: number, x: number): void {
  if (!Number.isSafeInteger(degree) || degree < 0) {
    throw new RangeError(`Le degré l doit être un entier sûr non négatif : ${degree}.`);
  }
  if (!Number.isSafeInteger(order) || order < 0 || order > degree) {
    throw new RangeError(`L'ordre m doit être un entier sûr tel que 0 <= m <= l : ${order}.`);
  }
  if (!Number.isFinite(x) || x < -1 || x > 1) {
    throw new RangeError(`L'argument x doit appartenir à l'intervalle [-1, 1] : ${x}.`);
  }
}

function requireFiniteFunction(value: number, degree: number, order: number): number {
  if (!Number.isFinite(value)) {
    throw new RangeError(`P_${degree}^${order}(x) dépasse le domaine fini de Number.`);
  }
  return value;
}

/**
 * Fonction de Ferrers P_l^m(x) de première espèce sur [-1, 1].
 *
 * La convention DLMF 14.7.8 est utilisée : la phase de Condon–Shortley (-1)^m est
 * intégrée. Le domaine de ce lot est l entier >= 0 et 0 <= m <= l ; m < 0 est réservé
 * au futur traitement des harmoniques sphériques.
 */
export function associatedLegendre(degree: number, order: number, x: number): number {
  validateInputs(degree, order, x);

  let diagonal = 1;
  if (order > 0) {
    const root = Math.sqrt((1 - x) * (1 + x));
    let oddFactor = 1;
    for (let currentOrder = 1; currentOrder <= order; currentOrder += 1) {
      diagonal = requireFiniteFunction(-oddFactor * root * diagonal, currentOrder, currentOrder);
      oddFactor += 2;
    }
  }

  if (degree === order) return diagonal;

  let previousPrevious = diagonal;
  let previous = requireFiniteFunction(x * (2 * order + 1) * diagonal, order + 1, order);
  if (degree === order + 1) return previous;

  for (let currentDegree = order + 2; currentDegree <= degree; currentDegree += 1) {
    const current =
      ((2 * currentDegree - 1) * x * previous - (currentDegree + order - 1) * previousPrevious) /
      (currentDegree - order);
    previousPrevious = previous;
    previous = requireFiniteFunction(current, currentDegree, order);
  }

  return previous;
}
