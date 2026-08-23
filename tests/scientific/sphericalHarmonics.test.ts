import { describe, expect, it } from 'vitest';

import { sphericalHarmonic } from '../../src/science/hydrogen/sphericalHarmonics';
import type { Complex } from '../../src/science/math/complex';
import { integrateSolidAngle } from './numericalIntegration';
import {
  ANALYTIC_COMPONENT_TOLERANCE,
  NUMERICAL_NORMALIZATION_ABSOLUTE_TOLERANCE,
  NUMERICAL_ORTHOGONALITY_ABSOLUTE_TOLERANCE,
} from './numericAssertions';

const ANGULAR_THETA_INTERVALS = 1_024;
const ANGULAR_PHI_INTERVALS = 64;

function expectComponent(actual: number, expected: number): void {
  expect(Math.abs(actual - expected)).toBeLessThanOrEqual(ANALYTIC_COMPONENT_TOLERANCE);
}

function expectComplexComponents(actual: Complex, expectedRe: number, expectedIm: number): void {
  expectComponent(actual.re, expectedRe);
  expectComponent(actual.im, expectedIm);
}

function modulusSquared(value: Complex): number {
  return value.re * value.re + value.im * value.im;
}

function angularOverlap(
  leftDegree: number,
  leftOrder: number,
  rightDegree: number,
  rightOrder: number,
): Complex {
  const overlapComponent = (imaginary: boolean): number =>
    integrateSolidAngle(
      (thetaRadians, phiRadians) => {
        const left = sphericalHarmonic(leftDegree, leftOrder, thetaRadians, phiRadians);
        const right = sphericalHarmonic(rightDegree, rightOrder, thetaRadians, phiRadians);
        return imaginary
          ? left.re * right.im - left.im * right.re
          : left.re * right.re + left.im * right.im;
      },
      ANGULAR_THETA_INTERVALS,
      ANGULAR_PHI_INTERVALS,
    );

  return { re: overlapComponent(false), im: overlapComponent(true) };
}

describe('harmoniques sphériques complexes normalisées', () => {
  it('retrouve Y_0^0, Y_1^0, Y_1^1 et Y_1^-1 composante par composante', () => {
    const points = [
      { theta: Math.PI / 3, phi: Math.PI / 5 },
      { theta: 0.41, phi: -1.37 },
    ] as const;

    for (const { theta, phi } of points) {
      const sinTheta = Math.sin(theta);
      const y00 = 1 / Math.sqrt(4 * Math.PI);
      const y10 = Math.sqrt(3 / (4 * Math.PI)) * Math.cos(theta);
      const orderOneAmplitude = Math.sqrt(3 / (8 * Math.PI)) * sinTheta;

      expectComplexComponents(sphericalHarmonic(0, 0, theta, phi), y00, 0);
      expectComplexComponents(sphericalHarmonic(1, 0, theta, phi), y10, 0);
      expectComplexComponents(
        sphericalHarmonic(1, 1, theta, phi),
        -orderOneAmplitude * Math.cos(phi),
        -orderOneAmplitude * Math.sin(phi),
      );
      expectComplexComponents(
        sphericalHarmonic(1, -1, theta, phi),
        orderOneAmplitude * Math.cos(phi),
        -orderOneAmplitude * Math.sin(phi),
      );
    }
  });

  it('ne double pas la phase de Condon–Short dans Y_1^1', () => {
    const y11 = sphericalHarmonic(1, 1, Math.PI / 2, 0);
    const y1Minus1 = sphericalHarmonic(1, -1, Math.PI / 2, 0);
    expect(y11.re).toBeLessThan(0);
    expect(y1Minus1.re).toBeGreaterThan(0);
    expectComponent(y11.im, 0);
    expectComponent(y1Minus1.im, 0);
  });

  it('applique la périodicité naturelle en phi sans restreindre les angles finis', () => {
    const reference = sphericalHarmonic(3, -2, 1.13, -0.72);
    for (const turns of [-7, -1, 1, 9]) {
      const periodic = sphericalHarmonic(3, -2, 1.13, -0.72 + turns * 2 * Math.PI);
      expectComplexComponents(periodic, reference.re, reference.im);
    }
  });

  it.each([
    [0, 0],
    [1, -1],
    [1, 0],
    [1, 1],
    [2, -2],
    [2, -1],
    [2, 0],
    [2, 1],
    [2, 2],
    [3, -3],
    [3, 0],
    [3, 2],
  ])('normalise numériquement Y_%i^%i avec la mesure sin(theta)', (degree, order) => {
    const normalization = integrateSolidAngle(
      (thetaRadians, phiRadians) =>
        modulusSquared(sphericalHarmonic(degree, order, thetaRadians, phiRadians)),
      ANGULAR_THETA_INTERVALS,
      ANGULAR_PHI_INTERVALS,
    );
    expect(Math.abs(normalization - 1)).toBeLessThanOrEqual(
      NUMERICAL_NORMALIZATION_ABSOLUTE_TOLERANCE,
    );
  });

  it('vérifie la convergence du maillage angulaire sur un cas défavorable m=0', () => {
    const integrateY30 = (thetaIntervals: number, phiIntervals: number): number =>
      integrateSolidAngle(
        (thetaRadians, phiRadians) =>
          modulusSquared(sphericalHarmonic(3, 0, thetaRadians, phiRadians)),
        thetaIntervals,
        phiIntervals,
      );

    const coarse = integrateY30(512, 32);
    const fine = integrateY30(ANGULAR_THETA_INTERVALS, ANGULAR_PHI_INTERVALS);
    expect(Math.abs(fine - 1)).toBeLessThanOrEqual(NUMERICAL_NORMALIZATION_ABSOLUTE_TOLERANCE);
    expect(Math.abs(fine - coarse)).toBeLessThanOrEqual(NUMERICAL_NORMALIZATION_ABSOLUTE_TOLERANCE);
  });

  it('intègre explicitement la mesure d’angle solide à 4 pi', () => {
    const solidAngle = integrateSolidAngle(() => 1, ANGULAR_THETA_INTERVALS, ANGULAR_PHI_INTERVALS);
    expect(Math.abs(solidAngle - 4 * Math.PI)).toBeLessThanOrEqual(
      NUMERICAL_NORMALIZATION_ABSOLUTE_TOLERANCE,
    );
  });

  it.each([
    [2, -1, 2, -1, 1],
    [2, -2, 2, 1, 0],
    [0, 0, 2, 0, 0],
    [1, 1, 3, 1, 0],
    [3, -2, 2, 2, 0],
  ])(
    'vérifie l’orthogonalité complexe de Y_%i^%i et Y_%i^%i',
    (leftDegree, leftOrder, rightDegree, rightOrder, expectedReal) => {
      const overlap = angularOverlap(leftDegree, leftOrder, rightDegree, rightOrder);
      expect(Math.abs(overlap.re - expectedReal)).toBeLessThanOrEqual(
        expectedReal === 0
          ? NUMERICAL_ORTHOGONALITY_ABSOLUTE_TOLERANCE
          : NUMERICAL_NORMALIZATION_ABSOLUTE_TOLERANCE,
      );
      expect(Math.abs(overlap.im)).toBeLessThanOrEqual(NUMERICAL_ORTHOGONALITY_ABSOLUTE_TOLERANCE);
    },
  );

  it.each([
    [1, 1, 0.67, -0.43],
    [2, 1, 1.31, 0.72],
    [2, 2, 2.19, -2.34],
    [3, 1, 0.92, 4.13],
    [3, 3, 2.42, -1.18],
  ])('respecte Y_%i^(-%i) = (-1)^m conjugate(Y_l^m)', (degree, positiveOrder, theta, phi) => {
    const positive = sphericalHarmonic(degree, positiveOrder, theta, phi);
    const negative = sphericalHarmonic(degree, -positiveOrder, theta, phi);
    const sign = positiveOrder % 2 === 0 ? 1 : -1;

    expectComplexComponents(negative, sign * positive.re, -sign * positive.im);
    expectComponent(modulusSquared(negative), modulusSquared(positive));
  });

  it.each([
    [0, 0, 0.71, -0.38],
    [1, 1, 1.12, 0.44],
    [2, -1, 0.83, -1.27],
    [3, 2, 2.07, 3.41],
  ])('respecte la parité (-1)^l de Y_%i^%i', (degree, order, theta, phi) => {
    const reference = sphericalHarmonic(degree, order, theta, phi);
    const inverted = sphericalHarmonic(degree, order, Math.PI - theta, phi + Math.PI);
    const sign = degree % 2 === 0 ? 1 : -1;
    expectComplexComponents(inverted, sign * reference.re, sign * reference.im);
  });

  it('rejette les degrés, ordres et angles hors domaine', () => {
    const invalidCalls: ReadonlyArray<() => Complex> = [
      () => sphericalHarmonic(-1, 0, 0, 0),
      () => sphericalHarmonic(1.5, 0, 0, 0),
      () => sphericalHarmonic(1, 2, 0, 0),
      () => sphericalHarmonic(2, 0.5, 0, 0),
      () => sphericalHarmonic(1, 0, -1e-12, 0),
      () => sphericalHarmonic(1, 0, Math.PI + 1e-12, 0),
      () => sphericalHarmonic(1, 0, Number.NaN, 0),
      () => sphericalHarmonic(1, 0, 0, Number.POSITIVE_INFINITY),
    ];

    for (const invalidCall of invalidCalls) {
      expect(invalidCall).toThrow(RangeError);
    }
  });
});
