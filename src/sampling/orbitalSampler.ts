import {
  REAL_ORBITAL_DEFINITIONS,
  realOrbitalWavefunction,
  type RealOrbitalName,
} from '../science/hydrogen/realOrbitals';
import { hydrogenWavefunction } from '../science/hydrogen/wavefunction';
import { complexPhase } from '../science/math/complex';
import { assertValidQuantumNumbers } from '../science/quantum/quantumNumbers';
import {
  ORBITAL_SAMPLER_VERSION,
  type OrbitalSampleSet,
  type OrbitalSamplingRequest,
  type OrbitalSamplingState,
} from './contracts';
import {
  angularDensityUpperBound,
  sampleComplexAngularDirection,
  sampleRealAngularDirection,
  type SampledAngularDirection,
} from './angularSampler';
import { getRadialCdf, sampleRadiusBohr } from './radialSampler';
import { createSeededRandom, SAMPLER_PRNG_VERSION } from './rng';

/**
 * Plafond de ressources du sampler pur : sept Float32 par point, soit environ
 * 28 Mio de buffers pour un million de points. Ce seuil est une politique
 * mémoire explicite, jamais une limite physique du modèle hydrogène.
 */
export const MAX_ORBITAL_SAMPLE_COUNT = 1_000_000;

function isRealOrbitalName(value: unknown): value is RealOrbitalName {
  return typeof value === 'string' && Object.hasOwn(REAL_ORBITAL_DEFINITIONS, value);
}

function requireSampleCount(sampleCount: number): void {
  if (
    !Number.isSafeInteger(sampleCount) ||
    sampleCount <= 0 ||
    sampleCount > MAX_ORBITAL_SAMPLE_COUNT
  ) {
    throw new RangeError(
      `Le nombre d'échantillons doit être un entier sûr dans [1, ${MAX_ORBITAL_SAMPLE_COUNT}] : ${sampleCount}.`,
    );
  }
}

function validateState(state: OrbitalSamplingState): number {
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

function copyState(state: OrbitalSamplingState): OrbitalSamplingState {
  return state.basis === 'complex'
    ? { basis: 'complex', n: state.n, l: state.l, m: state.m }
    : { basis: 'real', n: state.n, orbital: state.orbital };
}

function angularDirection(
  state: OrbitalSamplingState,
  rng: ReturnType<typeof createSeededRandom>,
): SampledAngularDirection {
  return state.basis === 'complex'
    ? sampleComplexAngularDirection(state.l, state.m, rng)
    : sampleRealAngularDirection(state.orbital, rng);
}

function wavefunctionPhase(
  state: OrbitalSamplingState,
  radiusBohr: number,
  direction: SampledAngularDirection,
): number | null {
  const amplitude =
    state.basis === 'complex'
      ? hydrogenWavefunction(
          state.n,
          state.l,
          state.m,
          radiusBohr,
          direction.thetaRadians,
          direction.phiRadians,
        )
      : realOrbitalWavefunction(
          state.n,
          state.orbital,
          radiusBohr,
          direction.thetaRadians,
          direction.phiRadians,
        );
  return complexPhase(amplitude);
}

/**
 * Pipeline commun Phase 4 : état scientifique → CDF radiale et densité angulaire
 * → positions cartésiennes en a₀. Les points sont des échantillons indépendants
 * de |ψ|²dV ; ils ne constituent ni électrons multiples ni trajectoires.
 */
export function sampleOrbital(request: OrbitalSamplingRequest): OrbitalSampleSet {
  requireSampleCount(request.sampleCount);
  const degree = validateState(request.state);
  const rng = createSeededRandom(request.seed);
  const radialCdf = getRadialCdf(request.state.n, degree);
  const positionsBohr = new Float32Array(request.sampleCount * 3);
  const radiiBohr = new Float32Array(request.sampleCount);
  const thetaRadians = new Float32Array(request.sampleCount);
  const phiRadians = new Float32Array(request.sampleCount);
  const phaseRadians = new Float32Array(request.sampleCount);
  let angularProposalCount = 0;
  let maximumObservedRadiusBohr = 0;
  let undefinedPhaseCount = 0;

  for (let sampleIndex = 0; sampleIndex < request.sampleCount; sampleIndex += 1) {
    const radiusBohr = sampleRadiusBohr(radialCdf, rng.nextFloat());
    const direction = angularDirection(request.state, rng);
    const vectorOffset = sampleIndex * 3;
    positionsBohr[vectorOffset] = radiusBohr * direction.xUnit;
    positionsBohr[vectorOffset + 1] = radiusBohr * direction.yUnit;
    positionsBohr[vectorOffset + 2] = radiusBohr * direction.zUnit;
    radiiBohr[sampleIndex] = radiusBohr;
    thetaRadians[sampleIndex] = direction.thetaRadians;
    phiRadians[sampleIndex] = direction.phiRadians;
    angularProposalCount += direction.proposalCount;
    maximumObservedRadiusBohr = Math.max(maximumObservedRadiusBohr, radiusBohr);

    const phase = wavefunctionPhase(request.state, radiusBohr, direction);
    if (phase === null) {
      phaseRadians[sampleIndex] = Number.NaN;
      undefinedPhaseCount += 1;
    } else {
      phaseRadians[sampleIndex] = phase;
    }
  }

  return {
    positionsBohr,
    radiiBohr,
    thetaRadians,
    phiRadians,
    phaseRadians,
    metadata: {
      angularAcceptanceRate: request.sampleCount / angularProposalCount,
      angularDensityUpperBoundInverseSteradians: angularDensityUpperBound(degree),
      angularProposalCount,
      initialPrngState: rng.initialState,
      maximumObservedRadiusBohr,
      normalizedSeed: rng.normalizedSeed,
      positionUnits: 'bohr',
      prngVersion: SAMPLER_PRNG_VERSION,
      radial: radialCdf.metadata,
      sampleCount: request.sampleCount,
      samplerVersion: ORBITAL_SAMPLER_VERSION,
      seed: request.seed,
      state: copyState(request.state),
      undefinedPhaseCount,
    },
  };
}
