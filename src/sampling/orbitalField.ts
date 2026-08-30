import {
  REAL_ORBITAL_DEFINITIONS,
  realOrbitalWavefunction,
  type RealOrbitalName,
} from '../science/hydrogen/realOrbitals';
import { hydrogenWavefunction } from '../science/hydrogen/wavefunction';
import { complexModulusSquared, complexPhase, type Complex } from '../science/math/complex';
import { assertValidQuantumNumbers } from '../science/quantum/quantumNumbers';
import type { OrbitalSamplingState } from './contracts';
import {
  ORBITAL_FIELD_VERSION,
  type OrbitalFieldGrid,
  type OrbitalFieldRequest,
} from './fieldContracts';

/**
 * Limites purement numériques du champ N³. À N=64, les trois buffers de sortie
 * représentent 786 432 valeurs Float32, indépendamment des buffers temporaires.
 * Ces limites ne restreignent aucun nombre quantique ni aucune grandeur physique.
 */
export const MIN_ORBITAL_FIELD_RESOLUTION = 8;
export const MAX_ORBITAL_FIELD_RESOLUTION = 64;

function isRealOrbitalName(value: unknown): value is RealOrbitalName {
  return typeof value === 'string' && Object.hasOwn(REAL_ORBITAL_DEFINITIONS, value);
}

/** Valide l'état et retourne son degré angulaire l sans réimplémenter ψ. */
export function orbitalDegree(state: OrbitalSamplingState): number {
  if (state.basis === 'complex') {
    assertValidQuantumNumbers(state);
    return state.l;
  }

  if (!isRealOrbitalName(state.orbital)) {
    throw new RangeError(`Nom d'orbitale réelle inconnu : ${String(state.orbital)}.`);
  }
  const definition = REAL_ORBITAL_DEFINITIONS[state.orbital];
  if (!Number.isSafeInteger(state.n) || state.n < definition.minimumPrincipalQuantumNumber) {
    throw new RangeError(
      `L'orbitale ${state.orbital} requiert un entier sûr n >= ${definition.minimumPrincipalQuantumNumber}.`,
    );
  }
  return definition.l;
}

/**
 * Évalue exclusivement les fonctions d'onde validées en Phase 3. Le rayon est
 * exprimé en a₀ ; theta est l'angle polaire depuis +z et phi l'azimut autour de z.
 */
export function evaluateOrbitalWavefunction(
  state: OrbitalSamplingState,
  radiusBohr: number,
  thetaRadians: number,
  phiRadians: number,
): Complex {
  return state.basis === 'complex'
    ? hydrogenWavefunction(state.n, state.l, state.m, radiusBohr, thetaRadians, phiRadians)
    : realOrbitalWavefunction(state.n, state.orbital, radiusBohr, thetaRadians, phiRadians);
}

function validateRequest(request: OrbitalFieldRequest): void {
  if (
    !Number.isSafeInteger(request.resolution) ||
    request.resolution < MIN_ORBITAL_FIELD_RESOLUTION ||
    request.resolution > MAX_ORBITAL_FIELD_RESOLUTION
  ) {
    throw new RangeError(
      `La résolution numérique doit être un entier sûr dans [${MIN_ORBITAL_FIELD_RESOLUTION}, ${MAX_ORBITAL_FIELD_RESOLUTION}] : ${request.resolution}.`,
    );
  }
  if (!Number.isFinite(request.extentBohr) || request.extentBohr <= 0) {
    throw new RangeError("L'étendue du champ doit être finie et strictement positive en a₀.");
  }
  if (!Number.isFinite(Math.hypot(request.extentBohr, request.extentBohr, request.extentBohr))) {
    throw new RangeError(
      "L'étendue du champ est trop grande pour conserver un rayon cartésien fini.",
    );
  }
}

function gridCoordinate(index: number, resolution: number, extentBohr: number): number {
  return extentBohr * (2 * (index / (resolution - 1)) - 1);
}

function sphericalCoordinates(
  xBohr: number,
  yBohr: number,
  zBohr: number,
): readonly [radiusBohr: number, thetaRadians: number, phiRadians: number] {
  const radiusBohr = Math.hypot(xBohr, yBohr, zBohr);
  if (radiusBohr === 0) return [0, 0, 0];

  // Le bornage absorbe uniquement l'arrondi de z/r aux pôles de la grille.
  const polarCosine = Math.max(-1, Math.min(1, zBohr / radiusBohr));
  return [radiusBohr, Math.acos(polarCosine), Math.atan2(yBohr, xBohr)];
}

/**
 * Calcule ψ sur une grille cartésienne uniforme couvrant [-extent,+extent]³.
 * L'indexation est x + N(y + Nz), donc x varie le plus vite. La densité fournie
 * est la densité volumique |ψ|² en a₀⁻³ avant normalisation par le maximum de la
 * grille ; elle n'inclut aucun jacobien sphérique.
 *
 * @see MIT OpenCourseWare 8.04, Lecture Notes 20–21, coordonnées sphériques.
 * @see MIT OpenCourseWare 8.04, Lecture Notes 22, séparation ψ=R_nl Y_l^m.
 */
export function computeOrbitalField(
  state: OrbitalSamplingState,
  request: OrbitalFieldRequest,
): OrbitalFieldGrid {
  orbitalDegree(state);
  validateRequest(request);

  const { extentBohr, resolution } = request;
  const valueCount = resolution ** 3;
  const rawDensity = new Float64Array(valueCount);
  const rawRealAmplitude = new Float64Array(valueCount);
  const phaseRadians = new Float32Array(valueCount);
  let maximumDensityPerCubicBohr = 0;
  let maximumWavefunctionAmplitude = 0;

  for (let zIndex = 0; zIndex < resolution; zIndex += 1) {
    const zBohr = gridCoordinate(zIndex, resolution, extentBohr);
    for (let yIndex = 0; yIndex < resolution; yIndex += 1) {
      const yBohr = gridCoordinate(yIndex, resolution, extentBohr);
      for (let xIndex = 0; xIndex < resolution; xIndex += 1) {
        const xBohr = gridCoordinate(xIndex, resolution, extentBohr);
        const fieldIndex = xIndex + resolution * (yIndex + resolution * zIndex);
        const [radiusBohr, thetaRadians, phiRadians] = sphericalCoordinates(xBohr, yBohr, zBohr);
        const amplitude = evaluateOrbitalWavefunction(state, radiusBohr, thetaRadians, phiRadians);
        const density = complexModulusSquared(amplitude);
        const amplitudeMagnitude = Math.sqrt(density);
        const phase = complexPhase(amplitude);

        rawDensity[fieldIndex] = density;
        rawRealAmplitude[fieldIndex] = amplitude.re;
        phaseRadians[fieldIndex] = phase ?? Number.NaN;
        maximumDensityPerCubicBohr = Math.max(maximumDensityPerCubicBohr, density);
        maximumWavefunctionAmplitude = Math.max(maximumWavefunctionAmplitude, amplitudeMagnitude);
      }
    }
  }

  if (maximumDensityPerCubicBohr <= 0 || maximumWavefunctionAmplitude <= 0) {
    throw new RangeError(
      'La grille ne contient aucune amplitude représentable ; réduisez son étendue ou augmentez sa résolution.',
    );
  }

  const densityNormalized = new Float32Array(valueCount);
  const signedAmplitudeNormalized = new Float32Array(valueCount);
  for (let index = 0; index < valueCount; index += 1) {
    const density = rawDensity[index];
    const realAmplitude = rawRealAmplitude[index];
    if (density === undefined || realAmplitude === undefined) {
      throw new RangeError(`Valeur de champ absente à l'indice ${index}.`);
    }
    densityNormalized[index] = density / maximumDensityPerCubicBohr;
    signedAmplitudeNormalized[index] = realAmplitude / maximumWavefunctionAmplitude;
  }

  return {
    densityNormalized,
    extentBohr,
    fieldVersion: ORBITAL_FIELD_VERSION,
    maximumDensityPerCubicBohr,
    maximumWavefunctionAmplitude,
    nodesAvailable: state.basis === 'real' || state.m === 0,
    phaseRadians,
    resolution,
    signedAmplitudeNormalized,
  };
}
