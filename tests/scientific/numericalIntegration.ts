const TWO_PI = 2 * Math.PI;

function requireFiniteValue(value: number, quantity: string): void {
  if (!Number.isFinite(value)) {
    throw new RangeError(`${quantity} doit être un nombre fini.`);
  }
}

function requirePositiveSafeInteger(value: number, quantity: string): void {
  if (!Number.isSafeInteger(value) || value <= 0) {
    throw new RangeError(`${quantity} doit être un entier sûr strictement positif.`);
  }
}

function evaluateFinite(integrand: (value: number) => number, value: number): number {
  const result = integrand(value);
  if (!Number.isFinite(result)) {
    throw new RangeError(`L'intégrande doit rester finie au point ${value}.`);
  }
  return result;
}

/**
 * Règle de Simpson composite déterministe sur un intervalle fini.
 * Le nombre de sous-intervalles doit être pair ; ce helper reste exclusivement
 * un outil de validation et ne constitue ni une CDF ni un intégrateur de production.
 */
export function compositeSimpson(
  integrand: (value: number) => number,
  lowerBound: number,
  upperBound: number,
  intervals: number,
): number {
  requireFiniteValue(lowerBound, 'La borne inférieure');
  requireFiniteValue(upperBound, 'La borne supérieure');
  requirePositiveSafeInteger(intervals, 'Le nombre de sous-intervalles');
  if (upperBound <= lowerBound) {
    throw new RangeError('La borne supérieure doit être strictement supérieure à la borne basse.');
  }
  if (intervals % 2 !== 0) {
    throw new RangeError('La règle de Simpson requiert un nombre pair de sous-intervalles.');
  }

  const step = (upperBound - lowerBound) / intervals;
  let weightedSum = evaluateFinite(integrand, lowerBound) + evaluateFinite(integrand, upperBound);

  for (let index = 1; index < intervals; index += 1) {
    const value = lowerBound + index * step;
    weightedSum += (index % 2 === 0 ? 2 : 4) * evaluateFinite(integrand, value);
  }

  const result = (weightedSum * step) / 3;
  if (!Number.isFinite(result)) {
    throw new RangeError("Le résultat de l'intégration de Simpson n'est pas fini.");
  }
  return result;
}

/**
 * Trapèzes périodiques sur [0, 2π), sans dupliquer le point 2π identique à zéro.
 * Cette règle est adaptée aux faibles modes de Fourier des harmoniques testées.
 */
export function integratePeriodicPhi(
  integrand: (phiRadians: number) => number,
  intervals: number,
): number {
  requirePositiveSafeInteger(intervals, "Le nombre d'intervalles azimutaux");
  const step = TWO_PI / intervals;
  let sum = 0;
  for (let index = 0; index < intervals; index += 1) {
    sum += evaluateFinite(integrand, index * step);
  }
  const result = sum * step;
  if (!Number.isFinite(result)) {
    throw new RangeError("Le résultat de l'intégration azimutale n'est pas fini.");
  }
  return result;
}

/**
 * Intégrale sur l'angle solide. La mesure sin(theta) dtheta dphi est appliquée
 * ici exactement une fois et reste indépendante du code scientifique testé.
 */
export function integrateSolidAngle(
  integrand: (thetaRadians: number, phiRadians: number) => number,
  thetaIntervals: number,
  phiIntervals: number,
): number {
  return compositeSimpson(
    (thetaRadians) =>
      Math.sin(thetaRadians) *
      integratePeriodicPhi((phiRadians) => integrand(thetaRadians, phiRadians), phiIntervals),
    0,
    Math.PI,
    thetaIntervals,
  );
}

/**
 * Intégrale sphérique complète. Le facteur r² et la mesure angulaire sont
 * appliqués par le helper de test, jamais par la densité volumique de production.
 */
export function integrateSphericalVolume(
  integrand: (rBohr: number, thetaRadians: number, phiRadians: number) => number,
  maximumRadiusBohr: number,
  radialIntervals: number,
  thetaIntervals: number,
  phiIntervals: number,
): number {
  requireFiniteValue(maximumRadiusBohr, 'Le rayon maximal');
  if (maximumRadiusBohr <= 0) {
    throw new RangeError('Le rayon maximal doit être strictement positif.');
  }
  return compositeSimpson(
    (rBohr) =>
      rBohr *
      rBohr *
      integrateSolidAngle(
        (thetaRadians, phiRadians) => integrand(rBohr, thetaRadians, phiRadians),
        thetaIntervals,
        phiIntervals,
      ),
    0,
    maximumRadiusBohr,
    radialIntervals,
  );
}
