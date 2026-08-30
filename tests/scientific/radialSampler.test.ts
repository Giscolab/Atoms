import { describe, expect, it } from 'vitest';

import { HYDROGEN_CHARACTERISTIC_RADIUS_BOHR } from '../../src/science/hydrogen/radialWavefunction';
import {
  createMonotonicCdfTable,
  evaluateMonotonicCdf,
  invertMonotonicCdf,
} from '../../src/sampling/cdf';
import { sampleOrbital } from '../../src/sampling/orbitalSampler';
import {
  clearRadialCdfCache,
  evaluateRadialCdf,
  getRadialCdf,
  sampleRadiusBohr,
} from '../../src/sampling/radialSampler';
import { dkwMassartBound, kolmogorovDistance } from './samplingStatistics';

const RADIAL_STATISTICAL_COMPARISONS = 6;

function oneSCdf(radiusBohr: number): number {
  if (radiusBohr <= 0) return 0;
  const x = radiusBohr / HYDROGEN_CHARACTERISTIC_RADIUS_BOHR;
  return 1 - Math.exp(-2 * x) * (1 + 2 * x + 2 * x * x);
}

function twoPCdf(radiusBohr: number): number {
  if (radiusBohr <= 0) return 0;
  const x = radiusBohr / HYDROGEN_CHARACTERISTIC_RADIUS_BOHR;
  return 1 - Math.exp(-x) * (1 + x + x ** 2 / 2 + x ** 3 / 6 + x ** 4 / 24);
}

function maximumCdfError(
  cdf: ReturnType<typeof getRadialCdf>,
  oracle: (radiusBohr: number) => number,
): number {
  let maximumError = 0;
  const maximumRadius = cdf.metadata.maximumRadiusBohr;
  for (let index = 0; index <= 2_000; index += 1) {
    const radius = (index / 2_000) * maximumRadius;
    maximumError = Math.max(
      maximumError,
      Math.abs(evaluateRadialCdf(cdf, radius) - oracle(radius)),
    );
  }
  return maximumError;
}

describe('CDF monotone générique', () => {
  it('inverse sans division nulle une CDF contenant un plateau', () => {
    const table = createMonotonicCdfTable(
      new Float64Array([0, 1, 2, 3]),
      new Float64Array([0, 0.25, 0.25, 1]),
    );
    expect(invertMonotonicCdf(table, 0)).toBe(0);
    expect(invertMonotonicCdf(table, 0.25)).toBe(1);
    expect(invertMonotonicCdf(table, 0.5)).toBeGreaterThan(2);
    expect(invertMonotonicCdf(table, 1)).toBe(3);
    expect(evaluateMonotonicCdf(table, 1.5)).toBe(0.25);
  });

  it('rejette les tables, coordonnées et probabilités invalides', () => {
    expect(() => createMonotonicCdfTable(new Float64Array([0]), new Float64Array([0]))).toThrow(
      RangeError,
    );
    expect(() =>
      createMonotonicCdfTable(new Float64Array([0, 1, 2]), new Float64Array([0, 0.8, 0.7])),
    ).toThrow(RangeError);
    const valid = createMonotonicCdfTable(new Float64Array([0, 1]), new Float64Array([0, 1]));
    expect(() => invertMonotonicCdf(valid, -Number.EPSILON)).toThrow(RangeError);
    expect(() => evaluateMonotonicCdf(valid, Number.NaN)).toThrow(RangeError);
  });
});

describe('CDF et sampling radial hydrogène', () => {
  it.each([
    [1, 0, oneSCdf],
    [2, 1, twoPCdf],
  ] as const)('converge vers la CDF fermée de l’état (%i,%i)', (n, l, oracle) => {
    const cdf = getRadialCdf(n, l);
    const error = maximumCdfError(cdf, oracle);
    expect(error).toBeLessThanOrEqual(2e-5);
    expect(cdf.metadata.cdfConvergenceError).toBeLessThanOrEqual(1e-8);
    expect(cdf.metadata.estimatedTruncatedMassUpperBound).toBeLessThanOrEqual(
      cdf.metadata.targetTruncatedMass,
    );
    expect(Math.abs(cdf.metadata.normalizationResidual)).toBeLessThanOrEqual(1e-8);
    expect(cdf.metadata.includedProbabilityMass).toBeGreaterThan(0.999_999_99);
    expect(cdf.metadata.maximumRadiusBohr).toBeGreaterThan(0);
    expect(cdf.metadata.cdfNodeCount).toBeGreaterThan(2);
  });

  it('normalise, inverse et met en cache par paramètres scientifiques et numériques', () => {
    clearRadialCdfCache();
    const first = getRadialCdf(2, 0);
    const same = getRadialCdf(2, 0);
    const differentResolution = getRadialCdf(2, 0, { initialIntegrationIntervals: 4_096 });
    expect(same).toBe(first);
    expect(differentResolution).not.toBe(first);
    for (const probability of [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1]) {
      const radius = sampleRadiusBohr(first, probability);
      expect(radius).toBeGreaterThanOrEqual(0);
      expect(Math.abs(evaluateRadialCdf(first, radius) - probability)).toBeLessThanOrEqual(1e-12);
    }
  });

  it('passe une comparaison ECDF/CDF 1s avec la borne DKW corrigée', () => {
    const sampleCount = 16_384;
    const result = sampleOrbital({
      sampleCount,
      seed: 0x5eed_1001,
      state: { basis: 'complex', n: 1, l: 0, m: 0 },
    });
    const numericalCdfError = maximumCdfError(getRadialCdf(1, 0), oneSCdf);
    const tailBudget = result.metadata.radial.estimatedTruncatedMassUpperBound;

    for (const prefixLength of [2_048, 8_192, sampleCount]) {
      const prefix = result.radiiBohr.slice(0, prefixLength);
      const distance = kolmogorovDistance(prefix, oneSCdf);
      const statisticalBound = dkwMassartBound(prefixLength, RADIAL_STATISTICAL_COMPARISONS);
      expect(distance).toBeLessThanOrEqual(statisticalBound + numericalCdfError + tailBudget);
    }
  });

  it('rejette les états et politiques numériques invalides', () => {
    expect(() => getRadialCdf(2, 2)).toThrow(RangeError);
    expect(() => getRadialCdf(1, 0, { initialIntegrationIntervals: 3 })).toThrow(RangeError);
    expect(() => getRadialCdf(1, 0, { targetTruncatedMass: 1 })).toThrow(RangeError);
    expect(() => getRadialCdf(1, 0, { cdfConvergenceTolerance: 0 })).toThrow(RangeError);
  });
});
