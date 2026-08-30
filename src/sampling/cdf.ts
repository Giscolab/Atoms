export interface MonotonicCdfTable {
  readonly coordinates: Float64Array;
  readonly probabilities: Float64Array;
}

function requireArrayValue(values: Float64Array, index: number, quantity: string): number {
  const value = values[index];
  if (value === undefined) {
    throw new RangeError(`${quantity} manque à l'indice ${index}.`);
  }
  return value;
}

/** Valide une CDF tabulée une seule fois, avant son utilisation intensive. */
export function createMonotonicCdfTable(
  coordinates: Float64Array,
  probabilities: Float64Array,
): MonotonicCdfTable {
  if (coordinates.length < 2 || coordinates.length !== probabilities.length) {
    throw new RangeError('Une CDF requiert au moins deux coordonnées et autant de probabilités.');
  }

  let previousCoordinate = requireArrayValue(coordinates, 0, 'La coordonnée de CDF');
  let previousProbability = requireArrayValue(probabilities, 0, 'La probabilité cumulée');
  if (!Number.isFinite(previousCoordinate) || previousProbability !== 0) {
    throw new RangeError(
      'La CDF doit commencer sur une coordonnée finie avec une probabilité nulle.',
    );
  }

  for (let index = 1; index < coordinates.length; index += 1) {
    const coordinate = requireArrayValue(coordinates, index, 'La coordonnée de CDF');
    const probability = requireArrayValue(probabilities, index, 'La probabilité cumulée');
    if (!Number.isFinite(coordinate) || coordinate <= previousCoordinate) {
      throw new RangeError(
        'Les coordonnées de CDF doivent être finies et strictement croissantes.',
      );
    }
    if (!Number.isFinite(probability) || probability < previousProbability || probability > 1) {
      throw new RangeError(
        'Les probabilités de CDF doivent être finies, croissantes et dans [0,1].',
      );
    }
    previousCoordinate = coordinate;
    previousProbability = probability;
  }

  if (previousProbability !== 1) {
    throw new RangeError('La dernière probabilité de la CDF doit valoir exactement un.');
  }
  return { coordinates, probabilities };
}

function requireUnitProbability(probability: number): void {
  if (!Number.isFinite(probability) || probability < 0 || probability > 1) {
    throw new RangeError(`La probabilité doit être finie et appartenir à [0,1] : ${probability}.`);
  }
}

function firstProbabilityAtLeast(table: MonotonicCdfTable, probability: number): number {
  let lower = 0;
  let upper = table.probabilities.length - 1;
  while (lower < upper) {
    const middle = lower + Math.floor((upper - lower) / 2);
    if (requireArrayValue(table.probabilities, middle, 'La probabilité cumulée') < probability) {
      lower = middle + 1;
    } else {
      upper = middle;
    }
  }
  return lower;
}

/** Inversion linéaire d'une CDF monotone ; les plateaux ne sont jamais sélectionnés. */
export function invertMonotonicCdf(table: MonotonicCdfTable, probability: number): number {
  requireUnitProbability(probability);
  const lastIndex = table.probabilities.length - 1;
  if (probability === 0) return requireArrayValue(table.coordinates, 0, 'La coordonnée de CDF');
  if (probability === 1) {
    return requireArrayValue(table.coordinates, lastIndex, 'La coordonnée de CDF');
  }

  const upperIndex = firstProbabilityAtLeast(table, probability);
  const lowerIndex = upperIndex - 1;
  const lowerProbability = requireArrayValue(
    table.probabilities,
    lowerIndex,
    'La probabilité cumulée',
  );
  const upperProbability = requireArrayValue(
    table.probabilities,
    upperIndex,
    'La probabilité cumulée',
  );
  const lowerCoordinate = requireArrayValue(table.coordinates, lowerIndex, 'La coordonnée de CDF');
  const upperCoordinate = requireArrayValue(table.coordinates, upperIndex, 'La coordonnée de CDF');

  const interpolation = (probability - lowerProbability) / (upperProbability - lowerProbability);
  return lowerCoordinate + interpolation * (upperCoordinate - lowerCoordinate);
}

/** Évalue par interpolation linéaire la CDF normalisée sur son domaine tabulé. */
export function evaluateMonotonicCdf(table: MonotonicCdfTable, coordinate: number): number {
  if (!Number.isFinite(coordinate)) {
    throw new RangeError(`La coordonnée de CDF doit être finie : ${coordinate}.`);
  }
  const firstCoordinate = requireArrayValue(table.coordinates, 0, 'La coordonnée de CDF');
  const lastIndex = table.coordinates.length - 1;
  const lastCoordinate = requireArrayValue(table.coordinates, lastIndex, 'La coordonnée de CDF');
  if (coordinate <= firstCoordinate) return 0;
  if (coordinate >= lastCoordinate) return 1;

  let lower = 0;
  let upper = lastIndex;
  while (upper - lower > 1) {
    const middle = lower + Math.floor((upper - lower) / 2);
    if (requireArrayValue(table.coordinates, middle, 'La coordonnée de CDF') <= coordinate) {
      lower = middle;
    } else {
      upper = middle;
    }
  }

  const lowerCoordinate = requireArrayValue(table.coordinates, lower, 'La coordonnée de CDF');
  const upperCoordinate = requireArrayValue(table.coordinates, upper, 'La coordonnée de CDF');
  const lowerProbability = requireArrayValue(table.probabilities, lower, 'La probabilité cumulée');
  const upperProbability = requireArrayValue(table.probabilities, upper, 'La probabilité cumulée');
  const interpolation = (coordinate - lowerCoordinate) / (upperCoordinate - lowerCoordinate);
  return lowerProbability + interpolation * (upperProbability - lowerProbability);
}
