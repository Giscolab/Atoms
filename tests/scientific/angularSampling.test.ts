import { describe, expect, it } from 'vitest';

import {
  REAL_ORBITAL_DEFINITIONS,
  realSphericalHarmonic,
  type RealOrbitalName,
} from '../../src/science/hydrogen/realOrbitals';
import { sphericalHarmonic } from '../../src/science/hydrogen/sphericalHarmonics';
import { complexModulusSquared } from '../../src/science/math/complex';
import { angularDensityUpperBound } from '../../src/sampling/angularSampler';
import { sampleOrbital } from '../../src/sampling/orbitalSampler';
import { dkwMassartBound, kolmogorovDistance, typedArrayValue } from './samplingStatistics';

const ANGULAR_STATISTICAL_COMPARISONS = 24;
const COMPLEX_SAMPLE_COUNT = 8_192;
const REAL_SAMPLE_COUNT = 4_096;

function complexPPlusMinusOneZCdf(z: number): number {
  if (z <= -1) return 0;
  if (z >= 1) return 1;
  return 0.5 + (3 * z) / 4 - z ** 3 / 4;
}

function axisSquaredCdf(axisCoordinate: number): number {
  if (axisCoordinate <= -1) return 0;
  if (axisCoordinate >= 1) return 1;
  return (axisCoordinate ** 3 + 1) / 2;
}

function uniformUnitCdf(value: number): number {
  return Math.max(0, Math.min(1, value));
}

function directionCoordinate(
  result: ReturnType<typeof sampleOrbital>,
  sampleIndex: number,
  component: 0 | 1 | 2,
): number {
  const radius = typedArrayValue(result.radiiBohr, sampleIndex, 'Le rayon échantillonné');
  return (
    typedArrayValue(
      result.positionsBohr,
      sampleIndex * 3 + component,
      'La position échantillonnée',
    ) / radius
  );
}

function mean(values: readonly number[]): number {
  return values.reduce((sum, value) => sum + value, 0) / values.length;
}

describe('sampling angulaire complexe', () => {
  it('applique la borne DLMF 14.30.9 aux états complexes et réels supportés', () => {
    for (let degree = 0; degree <= 4; degree += 1) {
      const upperBound = angularDensityUpperBound(degree);
      expect(Math.abs(upperBound - (2 * degree + 1) / (4 * Math.PI))).toBeLessThanOrEqual(
        16 * Number.EPSILON,
      );
      for (let order = -degree; order <= degree; order += 1) {
        for (let thetaIndex = 0; thetaIndex <= 32; thetaIndex += 1) {
          const theta = (thetaIndex * Math.PI) / 32;
          for (let phiIndex = 0; phiIndex < 32; phiIndex += 1) {
            const phi = (phiIndex * 2 * Math.PI) / 32;
            const density = complexModulusSquared(sphericalHarmonic(degree, order, theta, phi));
            expect(density).toBeLessThanOrEqual(upperBound * (1 + 512 * Number.EPSILON));
          }
        }
      }
    }

    for (const [name, definition] of Object.entries(REAL_ORBITAL_DEFINITIONS) as ReadonlyArray<
      readonly [RealOrbitalName, (typeof REAL_ORBITAL_DEFINITIONS)[RealOrbitalName]]
    >) {
      const upperBound = angularDensityUpperBound(definition.l);
      for (let thetaIndex = 0; thetaIndex <= 32; thetaIndex += 1) {
        const theta = (thetaIndex * Math.PI) / 32;
        for (let phiIndex = 0; phiIndex < 32; phiIndex += 1) {
          const phi = (phiIndex * 2 * Math.PI) / 32;
          const density = complexModulusSquared(realSphericalHarmonic(name, theta, phi));
          expect(density).toBeLessThanOrEqual(upperBound * (1 + 512 * Number.EPSILON));
        }
      }
    }
  });

  it('échantillonne 1s isotropiquement avec une acceptation unitaire', () => {
    const result = sampleOrbital({
      sampleCount: COMPLEX_SAMPLE_COUNT,
      seed: 0x1000_0001,
      state: { basis: 'complex', n: 1, l: 0, m: 0 },
    });
    const zUnit = Array.from({ length: COMPLEX_SAMPLE_COUNT }, (_, index) =>
      directionCoordinate(result, index, 2),
    );
    const distance = kolmogorovDistance(zUnit, (z) => uniformUnitCdf((z + 1) / 2));
    expect(distance).toBeLessThanOrEqual(
      dkwMassartBound(COMPLEX_SAMPLE_COUNT, ANGULAR_STATISTICAL_COMPARISONS),
    );
    expect(result.metadata.angularProposalCount).toBe(COMPLEX_SAMPLE_COUNT);
    expect(result.metadata.angularAcceptanceRate).toBe(1);
  });

  it('conserve la même densité pour m=+1 et m=-1, avec phi uniforme', () => {
    const request = {
      sampleCount: COMPLEX_SAMPLE_COUNT,
      seed: 0x2000_0002,
    } as const;
    const positive = sampleOrbital({
      ...request,
      state: { basis: 'complex', n: 2, l: 1, m: 1 },
    });
    const negative = sampleOrbital({
      ...request,
      state: { basis: 'complex', n: 2, l: 1, m: -1 },
    });
    expect(negative.positionsBohr).toEqual(positive.positionsBohr);
    expect(negative.radiiBohr).toEqual(positive.radiiBohr);
    expect(negative.thetaRadians).toEqual(positive.thetaRadians);
    expect(negative.phiRadians).toEqual(positive.phiRadians);

    const zUnit = Array.from({ length: COMPLEX_SAMPLE_COUNT }, (_, index) =>
      directionCoordinate(positive, index, 2),
    );
    const phiUnit = Array.from(positive.phiRadians, (phiRadians) => phiRadians / (2 * Math.PI));
    const bound = dkwMassartBound(COMPLEX_SAMPLE_COUNT, ANGULAR_STATISTICAL_COMPARISONS);
    expect(kolmogorovDistance(zUnit, complexPPlusMinusOneZCdf)).toBeLessThanOrEqual(bound);
    expect(kolmogorovDistance(phiUnit, uniformUnitCdf)).toBeLessThanOrEqual(bound);
  });
});

describe('sampling angulaire des vraies combinaisons réelles', () => {
  it.each([
    ['p_x', 0],
    ['p_y', 1],
    ['p_z', 2],
  ] as const)('oriente %s sur son axe cartésien', (orbital, component) => {
    const result = sampleOrbital({
      sampleCount: REAL_SAMPLE_COUNT,
      seed: 0x3000_0000 + component,
      state: { basis: 'real', n: 2, orbital },
    });
    const axisCoordinates = Array.from({ length: REAL_SAMPLE_COUNT }, (_, index) =>
      directionCoordinate(result, index, component),
    );
    expect(kolmogorovDistance(axisCoordinates, axisSquaredCdf)).toBeLessThanOrEqual(
      dkwMassartBound(REAL_SAMPLE_COUNT, ANGULAR_STATISTICAL_COMPARISONS),
    );
  });

  it.each([
    ['d_xy', [3 / 7, 3 / 7, 1 / 7, 1 / 7]],
    ['d_xz', [3 / 7, 1 / 7, 3 / 7, 1 / 7]],
    ['d_yz', [1 / 7, 3 / 7, 3 / 7, 1 / 7]],
    ['d_x2_y2', [3 / 7, 3 / 7, 1 / 7, 1 / 21]],
    ['d_z2', [5 / 21, 5 / 21, 11 / 21, null]],
  ] as const)('retrouve les moments directionnels de %s', (orbital, expected) => {
    const result = sampleOrbital({
      sampleCount: REAL_SAMPLE_COUNT,
      seed: 0x4000_0000 + Object.keys(REAL_ORBITAL_DEFINITIONS).indexOf(orbital),
      state: { basis: 'real', n: 3, orbital },
    });
    const squaredX: number[] = [];
    const squaredY: number[] = [];
    const squaredZ: number[] = [];
    for (let index = 0; index < REAL_SAMPLE_COUNT; index += 1) {
      squaredX.push(directionCoordinate(result, index, 0) ** 2);
      squaredY.push(directionCoordinate(result, index, 1) ** 2);
      squaredZ.push(directionCoordinate(result, index, 2) ** 2);
    }
    const hoeffdingBound = dkwMassartBound(REAL_SAMPLE_COUNT, ANGULAR_STATISTICAL_COMPARISONS);
    expect(Math.abs(mean(squaredX) - expected[0])).toBeLessThanOrEqual(hoeffdingBound);
    expect(Math.abs(mean(squaredY) - expected[1])).toBeLessThanOrEqual(hoeffdingBound);
    expect(Math.abs(mean(squaredZ) - expected[2])).toBeLessThanOrEqual(hoeffdingBound);
    if (expected[3] !== null) {
      const distinctiveMoment =
        orbital === 'd_xz'
          ? mean(squaredX.map((value, index) => value * (squaredZ[index] ?? Number.NaN)))
          : orbital === 'd_yz'
            ? mean(squaredY.map((value, index) => value * (squaredZ[index] ?? Number.NaN)))
            : mean(squaredX.map((value, index) => value * (squaredY[index] ?? Number.NaN)));
      expect(Math.abs(distinctiveMoment - expected[3])).toBeLessThanOrEqual(hoeffdingBound);
    }
  });
});
