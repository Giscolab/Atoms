import { describe, expect, it } from 'vitest';

import { hydrogenProbabilityDensity } from '../../src/science/hydrogen/probabilityDensity';
import {
  HYDROGEN_CHARACTERISTIC_RADIUS_BOHR,
  hydrogenRadialWavefunction,
} from '../../src/science/hydrogen/radialWavefunction';
import { sphericalHarmonic } from '../../src/science/hydrogen/sphericalHarmonics';
import {
  hydrogenWavefunction,
  hydrogenWavefunctionPhase,
} from '../../src/science/hydrogen/wavefunction';
import type { Complex } from '../../src/science/math/complex';
import {
  compositeSimpson,
  integrateSolidAngle,
  integrateSphericalVolume,
} from './numericalIntegration';
import {
  ANALYTIC_COMPONENT_TOLERANCE,
  NUMERICAL_NORMALIZATION_ABSOLUTE_TOLERANCE,
  NUMERICAL_ORTHOGONALITY_ABSOLUTE_TOLERANCE,
  THREE_DIMENSIONAL_CONVERGENCE_ABSOLUTE_TOLERANCE,
  THREE_DIMENSIONAL_NORMALIZATION_ABSOLUTE_TOLERANCE,
  phaseDistance,
} from './numericAssertions';

interface HydrogenState {
  readonly l: number;
  readonly m: number;
  readonly n: number;
}

const RADIAL_OVERLAP_INTERVALS = 16_384;
const OVERLAP_THETA_INTERVALS = 512;
const OVERLAP_PHI_INTERVALS = 32;

function expectComponent(actual: number, expected: number): void {
  const scale = Math.max(1, Math.abs(expected));
  expect(Math.abs(actual - expected)).toBeLessThanOrEqual(ANALYTIC_COMPONENT_TOLERANCE * scale);
}

function expectComplexComponents(actual: Complex, expectedRe: number, expectedIm: number): void {
  expectComponent(actual.re, expectedRe);
  expectComponent(actual.im, expectedIm);
}

function angularOverlap(left: HydrogenState, right: HydrogenState): Complex {
  const component = (imaginary: boolean): number =>
    integrateSolidAngle(
      (thetaRadians, phiRadians) => {
        const leftAngular = sphericalHarmonic(left.l, left.m, thetaRadians, phiRadians);
        const rightAngular = sphericalHarmonic(right.l, right.m, thetaRadians, phiRadians);
        return imaginary
          ? leftAngular.re * rightAngular.im - leftAngular.im * rightAngular.re
          : leftAngular.re * rightAngular.re + leftAngular.im * rightAngular.im;
      },
      OVERLAP_THETA_INTERVALS,
      OVERLAP_PHI_INTERVALS,
    );
  return { re: component(false), im: component(true) };
}

function radialOverlap(left: HydrogenState, right: HydrogenState): number {
  const maximumN = Math.max(left.n, right.n);
  const maximumRadiusBohr = 20 * maximumN * maximumN * HYDROGEN_CHARACTERISTIC_RADIUS_BOHR;
  return compositeSimpson(
    (rBohr) =>
      rBohr *
      rBohr *
      hydrogenRadialWavefunction(left.n, left.l, rBohr) *
      hydrogenRadialWavefunction(right.n, right.l, rBohr),
    0,
    maximumRadiusBohr,
    RADIAL_OVERLAP_INTERVALS,
  );
}

/** L'intégrale 3D se factorise ici suivant la séparation analytique R_nl Y_l^m. */
function separatedWavefunctionOverlap(left: HydrogenState, right: HydrogenState): Complex {
  const radial = radialOverlap(left, right);
  const angular = angularOverlap(left, right);
  return { re: radial * angular.re, im: radial * angular.im };
}

describe('fonction d’onde complexe et densité volumique', () => {
  it('valide le helper Simpson sur un polynôme cubique indépendant du moteur', () => {
    const integral = compositeSimpson((value) => value ** 3, 0, 1, 16);
    expectComponent(integral, 1 / 4);
  });

  it('retrouve directement la forme fermée de psi_100', () => {
    const a = HYDROGEN_CHARACTERISTIC_RADIUS_BOHR;
    const rBohr = 0.73 * a;
    const expectedAmplitude = Math.exp(-rBohr / a) / (Math.sqrt(Math.PI) * a ** 1.5);

    for (const [theta, phi] of [
      [0, 0],
      [Math.PI / 3, -0.82],
      [Math.PI, 7.1],
    ] as const) {
      expectComplexComponents(
        hydrogenWavefunction(1, 0, 0, rBohr, theta, phi),
        expectedAmplitude,
        0,
      );
    }
  });

  it('retrouve directement les parties réelle et imaginaire de psi_211', () => {
    const a = HYDROGEN_CHARACTERISTIC_RADIUS_BOHR;
    const rBohr = 1.37 * a;
    const theta = 1.11;
    const phi = -0.63;
    const signedAmplitude =
      (-(rBohr / a) * Math.exp(-rBohr / (2 * a)) * Math.sin(theta)) /
      (8 * Math.sqrt(Math.PI) * a ** 1.5);

    expectComplexComponents(
      hydrogenWavefunction(2, 1, 1, rBohr, theta, phi),
      signedAmplitude * Math.cos(phi),
      signedAmplitude * Math.sin(phi),
    );
  });

  it('compose explicitement R_nl et Y_l^m pour un état supérieur', () => {
    const state = { n: 3, l: 1, m: -1 } as const;
    const rBohr = 2.4 * HYDROGEN_CHARACTERISTIC_RADIUS_BOHR;
    const theta = 0.91;
    const phi = 1.37;
    const radial = hydrogenRadialWavefunction(state.n, state.l, rBohr);
    const angular = sphericalHarmonic(state.l, state.m, theta, phi);
    expectComplexComponents(
      hydrogenWavefunction(state.n, state.l, state.m, rBohr, theta, phi),
      radial * angular.re,
      radial * angular.im,
    );
  });

  it('expose |psi|² sans inclure r² ni sin(theta)', () => {
    const a = HYDROGEN_CHARACTERISTIC_RADIUS_BOHR;
    const radii = [0, 0.91 * a] as const;

    for (const rBohr of radii) {
      const expectedDensity = Math.exp((-2 * rBohr) / a) / (Math.PI * a ** 3);
      const polarDensity = hydrogenProbabilityDensity(1, 0, 0, rBohr, 0, -0.4);
      const equatorialDensity = hydrogenProbabilityDensity(1, 0, 0, rBohr, Math.PI / 2, 2.1);
      expectComponent(polarDensity, expectedDensity);
      expectComponent(equatorialDensity, expectedDensity);
      expectComponent(polarDensity, equatorialDensity);
    }

    expect(hydrogenProbabilityDensity(1, 0, 0, 0, 0, 0)).toBeGreaterThan(0);
  });

  it('sépare une densité azimutalement invariante de la phase de m=1', () => {
    const a = HYDROGEN_CHARACTERISTIC_RADIUS_BOHR;
    const rBohr = 1.4 * a;
    const theta = 1.03;
    const firstPhi = 0.21;
    const secondPhi = 1.07;
    const firstDensity = hydrogenProbabilityDensity(2, 1, 1, rBohr, theta, firstPhi);
    const secondDensity = hydrogenProbabilityDensity(2, 1, 1, rBohr, theta, secondPhi);
    const firstPhase = hydrogenWavefunctionPhase(2, 1, 1, rBohr, theta, firstPhi);
    const secondPhase = hydrogenWavefunctionPhase(2, 1, 1, rBohr, theta, secondPhi);

    expectComponent(firstDensity, secondDensity);
    expect(firstPhase).not.toBeNull();
    expect(secondPhase).not.toBeNull();
    if (firstPhase === null || secondPhase === null) {
      throw new Error('Une amplitude non nulle doit posséder une phase.');
    }
    expect(phaseDistance(secondPhase - firstPhase, secondPhi - firstPhi)).toBeLessThanOrEqual(
      ANALYTIC_COMPONENT_TOLERANCE,
    );
  });

  it('retourne null pour la phase sur un nœud exact et conserve la phase des lobes', () => {
    const a = HYDROGEN_CHARACTERISTIC_RADIUS_BOHR;
    expect(hydrogenWavefunctionPhase(2, 1, 1, 0, 0.8, 0.4)).toBeNull();
    expect(hydrogenProbabilityDensity(2, 1, 1, 0, 0.8, 0.4)).toBe(0);

    const phase1s = hydrogenWavefunctionPhase(1, 0, 0, a, 0.8, 0.4);
    const phase2sOuterLobe = hydrogenWavefunctionPhase(2, 0, 0, 3 * a, 0.8, 0.4);
    expect(phase1s).not.toBeNull();
    expect(phase2sOuterLobe).not.toBeNull();
    if (phase1s === null || phase2sOuterLobe === null) {
      throw new Error('Les points hors nœud doivent posséder une phase.');
    }
    expect(phaseDistance(phase1s, 0)).toBeLessThanOrEqual(ANALYTIC_COMPONENT_TOLERANCE);
    expect(phaseDistance(phase2sOuterLobe, Math.PI)).toBeLessThanOrEqual(
      ANALYTIC_COMPONENT_TOLERANCE,
    );
  });

  it('normalise réellement psi_211 en trois dimensions et vérifie la convergence', () => {
    const maximumRadiusBohr = 80 * HYDROGEN_CHARACTERISTIC_RADIUS_BOHR;
    const integrateDensity = (
      radialIntervals: number,
      thetaIntervals: number,
      phiIntervals: number,
    ): number =>
      integrateSphericalVolume(
        (rBohr, thetaRadians, phiRadians) =>
          hydrogenProbabilityDensity(2, 1, 1, rBohr, thetaRadians, phiRadians),
        maximumRadiusBohr,
        radialIntervals,
        thetaIntervals,
        phiIntervals,
      );

    const coarse = integrateDensity(512, 64, 8);
    const fine = integrateDensity(1_024, 128, 16);
    expect(Math.abs(fine - 1)).toBeLessThanOrEqual(
      THREE_DIMENSIONAL_NORMALIZATION_ABSOLUTE_TOLERANCE,
    );
    expect(Math.abs(fine - coarse)).toBeLessThanOrEqual(
      THREE_DIMENSIONAL_CONVERGENCE_ABSOLUTE_TOLERANCE,
    );
  });

  it.each([
    { n: 1, l: 0, m: 0 },
    { n: 2, l: 0, m: 0 },
    { n: 2, l: 1, m: -1 },
    { n: 3, l: 2, m: 2 },
  ])('normalise par séparation l’état (%j)', (state) => {
    const normalization = separatedWavefunctionOverlap(state, state);
    expect(Math.abs(normalization.re - 1)).toBeLessThanOrEqual(
      NUMERICAL_NORMALIZATION_ABSOLUTE_TOLERANCE,
    );
    expect(Math.abs(normalization.im)).toBeLessThanOrEqual(
      NUMERICAL_ORTHOGONALITY_ABSOLUTE_TOLERANCE,
    );
  });

  it.each([
    [
      { n: 1, l: 0, m: 0 },
      { n: 2, l: 0, m: 0 },
    ],
    [
      { n: 2, l: 0, m: 0 },
      { n: 2, l: 1, m: 0 },
    ],
    [
      { n: 2, l: 1, m: 0 },
      { n: 2, l: 1, m: 1 },
    ],
    [
      { n: 3, l: 1, m: 0 },
      { n: 3, l: 2, m: 0 },
    ],
  ])('vérifie une orthogonalité complète séparée entre (%j) et (%j)', (left, right) => {
    const overlap = separatedWavefunctionOverlap(left, right);
    expect(Math.abs(overlap.re)).toBeLessThanOrEqual(NUMERICAL_ORTHOGONALITY_ABSOLUTE_TOLERANCE);
    expect(Math.abs(overlap.im)).toBeLessThanOrEqual(NUMERICAL_ORTHOGONALITY_ABSOLUTE_TOLERANCE);
  });

  it('produit toujours une densité finie et positive ou nulle sur des états admissibles', () => {
    const samples = [
      [1, 0, 0, 0.3, 0.2, -0.7],
      [2, 0, 0, 2.3, 1.1, 0.4],
      [3, 1, -1, 4.2, 2.2, 8.1],
      [3, 2, 2, 1.7, Math.PI, -3.2],
    ] as const;

    for (const [n, l, m, rBohr, theta, phi] of samples) {
      const density = hydrogenProbabilityDensity(n, l, m, rBohr, theta, phi);
      expect(Number.isFinite(density)).toBe(true);
      expect(density).toBeGreaterThanOrEqual(0);
    }
  });

  it('rejette les états et coordonnées spatiales invalides', () => {
    const invalidCalls: ReadonlyArray<() => unknown> = [
      () => hydrogenWavefunction(0, 0, 0, 0, 0, 0),
      () => hydrogenWavefunction(2, 2, 0, 0, 0, 0),
      () => hydrogenWavefunction(2, 1, 2, 0, 0, 0),
      () => hydrogenWavefunction(1, 0, 0, -1e-12, 0, 0),
      () => hydrogenWavefunction(1, 0, 0, Number.NaN, 0, 0),
      () => hydrogenWavefunction(1, 0, 0, 0, -1e-12, 0),
      () => hydrogenWavefunction(1, 0, 0, 0, Math.PI + 1e-12, 0),
      () => hydrogenWavefunction(1, 0, 0, 0, 0, Number.POSITIVE_INFINITY),
      () => hydrogenProbabilityDensity(1, 0, 0, -1, 0, 0),
      () => hydrogenWavefunctionPhase(1, 0, 0, 0, Number.NaN, 0),
    ];

    for (const invalidCall of invalidCalls) {
      expect(invalidCall).toThrow(RangeError);
    }
  });
});
