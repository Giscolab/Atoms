import {
  REAL_ORBITAL_DEFINITIONS,
  realSphericalHarmonic,
  type RealOrbitalName,
} from '../science/hydrogen/realOrbitals';
import { sphericalHarmonic } from '../science/hydrogen/sphericalHarmonics';
import { complexModulusSquared, type Complex } from '../science/math/complex';
import type { SeededRandomGenerator } from './rng';

const TWO_PI = 2 * Math.PI;
const DENSITY_BOUND_ROUNDING_MARGIN = 512 * Number.EPSILON;

export interface SampledAngularDirection {
  readonly amplitude: Complex;
  readonly densityInverseSteradians: number;
  readonly phiRadians: number;
  readonly proposalCount: number;
  readonly thetaRadians: number;
  readonly xUnit: number;
  readonly yUnit: number;
  readonly zUnit: number;
}

/**
 * Théorème d'addition DLMF 14.30.9 :
 * Σ_m |Y_l^m(Ω)|² = (2l+1)/(4π).
 * Il borne chaque harmonique et, par Cauchy–Schwarz, toute combinaison réelle
 * normalisée de même l. L'unité est sr⁻¹.
 *
 * @see https://dlmf.nist.gov/14.30.E9
 */
export function angularDensityUpperBound(degree: number): number {
  if (!Number.isSafeInteger(degree) || degree < 0) {
    throw new RangeError(`Le degré angulaire doit être un entier sûr positif ou nul : ${degree}.`);
  }
  const bound = (2 * degree + 1) / (4 * Math.PI);
  if (!Number.isFinite(bound) || bound <= 0) {
    throw new RangeError("La borne de densité angulaire n'est pas représentable.");
  }
  return bound;
}

function boundedDensity(amplitude: Complex, upperBound: number): number {
  const density = complexModulusSquared(amplitude);
  if (density > upperBound * (1 + DENSITY_BOUND_ROUNDING_MARGIN)) {
    throw new RangeError(
      `La densité angulaire ${density} dépasse sa borne scientifique ${upperBound}.`,
    );
  }
  return Math.min(density, upperBound);
}

function directionFromAngles(
  amplitude: Complex,
  densityInverseSteradians: number,
  proposalCount: number,
  zUnit: number,
  phiRadians: number,
): SampledAngularDirection {
  const thetaRadians = Math.acos(zUnit);
  const transverseRadius = Math.sqrt(Math.max(0, 1 - zUnit * zUnit));
  // Convention sphérique du noyau : x=sin(theta)cos(phi), y=sin(theta)sin(phi), z=cos(theta).
  return {
    amplitude,
    densityInverseSteradians,
    phiRadians,
    proposalCount,
    thetaRadians,
    xUnit: transverseRadius * Math.cos(phiRadians),
    yUnit: transverseRadius * Math.sin(phiRadians),
    zUnit,
  };
}

/**
 * Rejet depuis une direction uniforme sur la sphère : z=cos(theta) est uniforme
 * sur [-1,1] et phi sur [0,2π). La proposition porte donc déjà exactement la
 * mesure dΩ = sin(theta)dtheta dphi ; aucun facteur sin(theta) supplémentaire
 * n'est appliqué à |Y_l^m|².
 */
export function sampleComplexAngularDirection(
  degree: number,
  order: number,
  rng: SeededRandomGenerator,
): SampledAngularDirection {
  const upperBound = angularDensityUpperBound(degree);
  for (let proposalCount = 1; proposalCount <= Number.MAX_SAFE_INTEGER; proposalCount += 1) {
    const zUnit = 2 * rng.nextFloat() - 1;
    const thetaRadians = Math.acos(zUnit);
    // |Y_l^m|² est indépendant de phi : le rejet se fait à phi=0, puis phi est tiré uniformément.
    const proposalAmplitude = sphericalHarmonic(degree, order, thetaRadians, 0);
    const density = boundedDensity(proposalAmplitude, upperBound);
    if (rng.nextFloat() * upperBound < density) {
      const phiRadians = TWO_PI * rng.nextFloat();
      const amplitude = sphericalHarmonic(degree, order, thetaRadians, phiRadians);
      return directionFromAngles(amplitude, density, proposalCount, zUnit, phiRadians);
    }
  }
  throw new RangeError("L'échantillonnage angulaire complexe a épuisé le compteur sûr d'essais.");
}

export function sampleRealAngularDirection(
  name: RealOrbitalName,
  rng: SeededRandomGenerator,
): SampledAngularDirection {
  if (!Object.hasOwn(REAL_ORBITAL_DEFINITIONS, name)) {
    throw new RangeError(`Nom d'orbitale réelle inconnu : ${name}.`);
  }
  const degree = REAL_ORBITAL_DEFINITIONS[name].l;
  const upperBound = angularDensityUpperBound(degree);
  for (let proposalCount = 1; proposalCount <= Number.MAX_SAFE_INTEGER; proposalCount += 1) {
    const zUnit = 2 * rng.nextFloat() - 1;
    const phiRadians = TWO_PI * rng.nextFloat();
    const thetaRadians = Math.acos(zUnit);
    const amplitude = realSphericalHarmonic(name, thetaRadians, phiRadians);
    const density = boundedDensity(amplitude, upperBound);
    if (rng.nextFloat() * upperBound < density) {
      return directionFromAngles(amplitude, density, proposalCount, zUnit, phiRadians);
    }
  }
  throw new RangeError("L'échantillonnage angulaire réel a épuisé le compteur sûr d'essais.");
}
