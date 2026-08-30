import { assertValidQuantumNumbers } from '../science/quantum/quantumNumbers';
import {
  HYDROGEN_CHARACTERISTIC_RADIUS_BOHR,
  hydrogenRadialWavefunction,
} from '../science/hydrogen/radialWavefunction';
import type { RadialSamplingMetadata } from './contracts';
import {
  createMonotonicCdfTable,
  evaluateMonotonicCdf,
  invertMonotonicCdf,
  type MonotonicCdfTable,
} from './cdf';

export const RADIAL_CDF_CACHE_VERSION = 'radial-cdf-simpson-adaptive-v1' as const;

export interface RadialCdfOptions {
  readonly cdfConvergenceTolerance?: number;
  readonly initialIntegrationIntervals?: number;
  readonly initialMaximumRadiusInNSquaredAMu?: number;
  readonly maximumDomainExpansions?: number;
  readonly maximumResolutionRefinements?: number;
  readonly targetTruncatedMass?: number;
}

interface ResolvedRadialCdfOptions {
  readonly cdfConvergenceTolerance: number;
  readonly initialIntegrationIntervals: number;
  readonly initialMaximumRadiusInNSquaredAMu: number;
  readonly maximumDomainExpansions: number;
  readonly maximumResolutionRefinements: number;
  readonly targetTruncatedMass: number;
}

interface RawRadialCdf {
  readonly cumulativeMass: Float64Array;
  readonly integrationIntervalCount: number;
  readonly mass: number;
  readonly radiiBohr: Float64Array;
}

export interface RadialCdf {
  readonly metadata: RadialSamplingMetadata;
}

const DEFAULT_OPTIONS: ResolvedRadialCdfOptions = {
  cdfConvergenceTolerance: 1e-8,
  initialIntegrationIntervals: 2_048,
  initialMaximumRadiusInNSquaredAMu: 8,
  maximumDomainExpansions: 6,
  maximumResolutionRefinements: 6,
  targetTruncatedMass: 1e-8,
};

const radialCdfCache = new Map<string, RadialCdf>();
const radialCdfTables = new WeakMap<RadialCdf, MonotonicCdfTable>();

function requirePositiveFinite(value: number, quantity: string): void {
  if (!Number.isFinite(value) || value <= 0) {
    throw new RangeError(`${quantity} doit être un nombre fini strictement positif.`);
  }
}

function requireNonNegativeSafeInteger(value: number, quantity: string): void {
  if (!Number.isSafeInteger(value) || value < 0) {
    throw new RangeError(`${quantity} doit être un entier sûr positif ou nul.`);
  }
}

function resolveOptions(options: RadialCdfOptions | undefined): ResolvedRadialCdfOptions {
  const resolved = { ...DEFAULT_OPTIONS, ...options };
  requirePositiveFinite(resolved.cdfConvergenceTolerance, 'La tolérance de convergence CDF');
  requirePositiveFinite(resolved.initialMaximumRadiusInNSquaredAMu, 'Le domaine radial initial');
  requirePositiveFinite(resolved.targetTruncatedMass, 'La masse tronquée cible');
  if (resolved.targetTruncatedMass >= 1) {
    throw new RangeError('La masse tronquée cible doit être strictement inférieure à un.');
  }
  if (
    !Number.isSafeInteger(resolved.initialIntegrationIntervals) ||
    resolved.initialIntegrationIntervals < 2 ||
    resolved.initialIntegrationIntervals % 2 !== 0
  ) {
    throw new RangeError("Le nombre initial d'intervalles doit être un entier sûr pair >= 2.");
  }
  requireNonNegativeSafeInteger(
    resolved.maximumDomainExpansions,
    "Le nombre maximal d'extensions du domaine",
  );
  requireNonNegativeSafeInteger(
    resolved.maximumResolutionRefinements,
    'Le nombre maximal de raffinements',
  );
  return resolved;
}

function radialProbabilityDensity(n: number, l: number, rBohr: number): number {
  const amplitude = hydrogenRadialWavefunction(n, l, rBohr);
  const density = rBohr * rBohr * amplitude * amplitude;
  if (!Number.isFinite(density) || density < 0) {
    throw new RangeError('La densité radiale r²|R_nl|² doit rester finie et positive ou nulle.');
  }
  return density;
}

/**
 * CDF cumulative par panneaux de Simpson composites. Chaque incrément
 * h(f₀+4f₁+f₂)/3 est positif, ce qui garantit la monotonie avant normalisation.
 * La mesure est P_nl(r)dr = r²|R_nl(r)|²dr, avec r exprimé en a₀.
 *
 * @see MIT OpenCourseWare 8.04, Lecture Notes 22, équations radiales de l'hydrogène.
 */
function integrateRadialCdf(
  n: number,
  l: number,
  maximumRadiusBohr: number,
  integrationIntervals: number,
): RawRadialCdf {
  const panelCount = integrationIntervals / 2;
  const radiiBohr = new Float64Array(panelCount + 1);
  const cumulativeMass = new Float64Array(panelCount + 1);
  const step = maximumRadiusBohr / integrationIntervals;
  let leftDensity = radialProbabilityDensity(n, l, 0);
  let mass = 0;

  for (let panel = 0; panel < panelCount; panel += 1) {
    const leftRadius = 2 * panel * step;
    const middleRadius = leftRadius + step;
    const rightRadius = leftRadius + 2 * step;
    const middleDensity = radialProbabilityDensity(n, l, middleRadius);
    const rightDensity = radialProbabilityDensity(n, l, rightRadius);
    const panelMass = (step * (leftDensity + 4 * middleDensity + rightDensity)) / 3;
    if (!Number.isFinite(panelMass) || panelMass < 0) {
      throw new RangeError("Un panneau d'intégration radiale doit avoir une masse finie positive.");
    }
    mass += panelMass;
    radiiBohr[panel + 1] = rightRadius;
    cumulativeMass[panel + 1] = mass;
    leftDensity = rightDensity;
  }

  if (!Number.isFinite(mass) || mass <= 0) {
    throw new RangeError("La masse de probabilité radiale intégrée n'est pas exploitable.");
  }
  return { cumulativeMass, integrationIntervalCount: integrationIntervals, mass, radiiBohr };
}

function rawArrayValue(values: Float64Array, index: number, quantity: string): number {
  const value = values[index];
  if (value === undefined) throw new RangeError(`${quantity} manque à l'indice ${index}.`);
  return value;
}

function normalizedProbability(raw: RawRadialCdf, index: number): number {
  return rawArrayValue(raw.cumulativeMass, index, 'La masse radiale cumulée') / raw.mass;
}

function cdfConvergenceError(coarse: RawRadialCdf, fine: RawRadialCdf): number {
  let maximumDifference = 0;
  for (let coarseIndex = 0; coarseIndex < coarse.cumulativeMass.length; coarseIndex += 1) {
    const fineIndex = 2 * coarseIndex;
    const difference = Math.abs(
      normalizedProbability(coarse, coarseIndex) - normalizedProbability(fine, fineIndex),
    );
    maximumDifference = Math.max(maximumDifference, difference);
  }
  return maximumDifference;
}

function normalizeCdf(raw: RawRadialCdf): MonotonicCdfTable {
  const probabilities = new Float64Array(raw.cumulativeMass.length);
  const lastIndex = probabilities.length - 1;
  for (let index = 0; index < lastIndex; index += 1) {
    probabilities[index] = Math.min(1, normalizedProbability(raw, index));
  }
  probabilities[lastIndex] = 1;
  return createMonotonicCdfTable(raw.radiiBohr, probabilities);
}

function cacheKey(n: number, l: number, options: ResolvedRadialCdfOptions): string {
  return [
    RADIAL_CDF_CACHE_VERSION,
    n,
    l,
    options.cdfConvergenceTolerance,
    options.initialIntegrationIntervals,
    options.initialMaximumRadiusInNSquaredAMu,
    options.maximumDomainExpansions,
    options.maximumResolutionRefinements,
    options.targetTruncatedMass,
  ].join(':');
}

/**
 * Construit puis met en cache une CDF radiale normalisée. Le domaine est doublé
 * jusqu'à ce que le résidu de la normalisation analytique unitaire, augmenté de
 * l'écart entre deux résolutions Simpson, respecte la masse tronquée cible.
 */
export function getRadialCdf(n: number, l: number, options?: RadialCdfOptions): RadialCdf {
  assertValidQuantumNumbers({ n, l, m: 0 });
  const resolved = resolveOptions(options);
  const key = cacheKey(n, l, resolved);
  const cached = radialCdfCache.get(key);
  if (cached !== undefined) return cached;

  const initialMaximumRadiusBohr =
    resolved.initialMaximumRadiusInNSquaredAMu * n * n * HYDROGEN_CHARACTERISTIC_RADIUS_BOHR;
  requirePositiveFinite(initialMaximumRadiusBohr, 'Le domaine radial initial en a₀');

  let maximumRadiusBohr = initialMaximumRadiusBohr;
  for (
    let domainExpansionCount = 0;
    domainExpansionCount <= resolved.maximumDomainExpansions;
    domainExpansionCount += 1
  ) {
    let coarse = integrateRadialCdf(n, l, maximumRadiusBohr, resolved.initialIntegrationIntervals);
    let convergedFine: RawRadialCdf | undefined;
    let measuredCdfError = Number.POSITIVE_INFINITY;
    let measuredMassError = Number.POSITIVE_INFINITY;
    let resolutionRefinementCount = 0;

    for (let refinement = 0; refinement <= resolved.maximumResolutionRefinements; refinement += 1) {
      const fine = integrateRadialCdf(n, l, maximumRadiusBohr, coarse.integrationIntervalCount * 2);
      measuredCdfError = cdfConvergenceError(coarse, fine);
      measuredMassError = Math.abs(fine.mass - coarse.mass);
      resolutionRefinementCount = refinement + 1;
      if (
        measuredCdfError <= resolved.cdfConvergenceTolerance &&
        measuredMassError <= resolved.cdfConvergenceTolerance
      ) {
        convergedFine = fine;
        break;
      }
      coarse = fine;
    }

    if (convergedFine === undefined) {
      throw new RangeError(
        `La CDF radiale (${n},${l}) ne converge pas en résolution sur [0, ${maximumRadiusBohr}] a₀.`,
      );
    }

    const normalizationResidual = 1 - convergedFine.mass;
    const truncatedProbabilityMass = Math.max(0, normalizationResidual);
    const estimatedTruncatedMassUpperBound = Math.abs(normalizationResidual) + measuredMassError;
    if (estimatedTruncatedMassUpperBound <= resolved.targetTruncatedMass) {
      const radialCdf: RadialCdf = {
        metadata: {
          cacheVersion: RADIAL_CDF_CACHE_VERSION,
          cdfConvergenceError: measuredCdfError,
          cdfNodeCount: convergedFine.cumulativeMass.length,
          domainExpansionCount,
          estimatedQuadratureMassError: measuredMassError,
          estimatedTruncatedMassUpperBound,
          includedProbabilityMass: convergedFine.mass,
          integrationIntervalCount: convergedFine.integrationIntervalCount,
          maximumRadiusBohr,
          normalizationResidual,
          resolutionRefinementCount,
          targetTruncatedMass: resolved.targetTruncatedMass,
          truncatedProbabilityMass,
        },
      };
      radialCdfTables.set(radialCdf, normalizeCdf(convergedFine));
      radialCdfCache.set(key, radialCdf);
      return radialCdf;
    }

    maximumRadiusBohr *= 2;
    requirePositiveFinite(maximumRadiusBohr, 'Le domaine radial étendu en a₀');
  }

  throw new RangeError(
    `La masse radiale tronquée de l'état (${n},${l}) ne converge pas vers ${resolved.targetTruncatedMass}.`,
  );
}

function requireRadialTable(cdf: RadialCdf): MonotonicCdfTable {
  const table = radialCdfTables.get(cdf);
  if (table === undefined) {
    throw new RangeError("La CDF radiale n'a pas été construite par ce module.");
  }
  return table;
}

export function sampleRadiusBohr(cdf: RadialCdf, unitProbability: number): number {
  return invertMonotonicCdf(requireRadialTable(cdf), unitProbability);
}

export function evaluateRadialCdf(cdf: RadialCdf, radiusBohr: number): number {
  return evaluateMonotonicCdf(requireRadialTable(cdf), radiusBohr);
}

/** Réservé aux tests et aux changements explicites de version numérique. */
export function clearRadialCdfCache(): void {
  radialCdfCache.clear();
}
