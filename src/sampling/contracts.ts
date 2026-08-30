import type { RealOrbitalName } from '../science/hydrogen/realOrbitals';

/** Version du contrat complet de reproductibilité du sampler scientifique. */
export const ORBITAL_SAMPLER_VERSION = 'hydrogen-orbital-sampler-v1' as const;

export interface ComplexOrbitalState {
  readonly basis: 'complex';
  readonly l: number;
  readonly m: number;
  readonly n: number;
}

export interface RealOrbitalState {
  readonly basis: 'real';
  readonly n: number;
  readonly orbital: RealOrbitalName;
}

/**
 * État stationnaire explicite. Une orbitale réelle est une combinaison normalisée
 * des états complexes ±m et n'est jamais assimilée à un m individuel.
 */
export type OrbitalSamplingState = ComplexOrbitalState | RealOrbitalState;

export interface OrbitalSamplingRequest {
  readonly sampleCount: number;
  /** Entier non signé sur 32 bits ; voir le contrat de transformation dans rng.ts. */
  readonly seed: number;
  readonly state: OrbitalSamplingState;
}

export interface RadialSamplingMetadata {
  readonly cacheVersion: string;
  readonly cdfConvergenceError: number;
  readonly cdfNodeCount: number;
  readonly domainExpansionCount: number;
  readonly estimatedQuadratureMassError: number;
  readonly estimatedTruncatedMassUpperBound: number;
  readonly includedProbabilityMass: number;
  readonly integrationIntervalCount: number;
  readonly maximumRadiusBohr: number;
  readonly normalizationResidual: number;
  readonly resolutionRefinementCount: number;
  readonly targetTruncatedMass: number;
  readonly truncatedProbabilityMass: number;
}

export interface OrbitalSamplingMetadata {
  readonly angularAcceptanceRate: number;
  readonly angularDensityUpperBoundInverseSteradians: number;
  readonly angularProposalCount: number;
  readonly initialPrngState: number;
  readonly maximumObservedRadiusBohr: number;
  readonly normalizedSeed: number;
  readonly positionUnits: 'bohr';
  readonly prngVersion: string;
  readonly radial: RadialSamplingMetadata;
  readonly sampleCount: number;
  readonly samplerVersion: typeof ORBITAL_SAMPLER_VERSION;
  readonly seed: number;
  readonly state: OrbitalSamplingState;
  readonly undefinedPhaseCount: number;
}

export interface OrbitalSampleSet {
  /** Angle azimutal phi dans [0,2π), en radians. */
  readonly phiRadians: Float32Array;
  /** Phase de ψ en radians ; NaN marque exclusivement une amplitude exactement nulle. */
  readonly phaseRadians: Float32Array;
  /** Triplets cartésiens x,y,z en a₀, avec l'axe polaire porté par z. */
  readonly positionsBohr: Float32Array;
  /** Rayon r/a₀ avant conversion cartésienne. */
  readonly radiiBohr: Float32Array;
  /** Angle polaire theta dans [0,π], en radians. */
  readonly thetaRadians: Float32Array;
  readonly metadata: OrbitalSamplingMetadata;
}
