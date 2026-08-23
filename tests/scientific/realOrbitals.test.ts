import { describe, expect, it } from 'vitest';

import {
  complexModulusSquared,
  conjugateComplex,
  multiplyComplex,
  type Complex,
} from '../../src/science/math/complex';
import {
  hydrogenRadialWavefunction,
  HYDROGEN_CHARACTERISTIC_RADIUS_BOHR,
} from '../../src/science/hydrogen/radialWavefunction';
import {
  REAL_ORBITAL_DEFINITIONS,
  realOrbitalWavefunction,
  realSphericalHarmonic,
  type RealOrbitalName,
} from '../../src/science/hydrogen/realOrbitals';
import { integrateSolidAngle, integrateSphericalVolume } from './numericalIntegration';
import {
  ANALYTIC_COMPONENT_TOLERANCE,
  NUMERICAL_NORMALIZATION_ABSOLUTE_TOLERANCE,
  NUMERICAL_ORTHOGONALITY_ABSOLUTE_TOLERANCE,
  THREE_DIMENSIONAL_NORMALIZATION_ABSOLUTE_TOLERANCE,
} from './numericAssertions';

const P_ORBITALS = ['p_x', 'p_y', 'p_z'] as const satisfies readonly RealOrbitalName[];
const D_ORBITALS = [
  'd_xy',
  'd_xz',
  'd_yz',
  'd_x2_y2',
  'd_z2',
] as const satisfies readonly RealOrbitalName[];
const REAL_ORBITALS = [...P_ORBITALS, ...D_ORBITALS] as const;

const ANGULAR_THETA_INTERVALS = 320;
const ANGULAR_PHI_INTERVALS = 64;
const VOLUME_RADIAL_INTERVALS = 400;
const VOLUME_THETA_INTERVALS = 192;
const VOLUME_PHI_INTERVALS = 12;
const RADIAL_EXTENT_IN_N_SQUARED_A_MU = 16;

function expectComplexComponents(
  actual: Complex,
  expectedRe: number,
  expectedIm = 0,
  tolerance = ANALYTIC_COMPONENT_TOLERANCE,
): void {
  expect(Math.abs(actual.re - expectedRe)).toBeLessThanOrEqual(tolerance);
  expect(Math.abs(actual.im - expectedIm)).toBeLessThanOrEqual(tolerance);
}

function cartesianDirection(
  thetaRadians: number,
  phiRadians: number,
): {
  readonly x: number;
  readonly y: number;
  readonly z: number;
} {
  const sinTheta = Math.sin(thetaRadians);
  return {
    x: sinTheta * Math.cos(phiRadians),
    y: sinTheta * Math.sin(phiRadians),
    z: Math.cos(thetaRadians),
  };
}

function signedCartesianOracle(
  name: RealOrbitalName,
  thetaRadians: number,
  phiRadians: number,
): number {
  const { x, y, z } = cartesianDirection(thetaRadians, phiRadians);
  const pNormalization = Math.sqrt(3 / (4 * Math.PI));
  const dCrossNormalization = Math.sqrt(15 / (4 * Math.PI));

  switch (name) {
    case 'p_x':
      return pNormalization * x;
    case 'p_y':
      return pNormalization * y;
    case 'p_z':
      return pNormalization * z;
    case 'd_xy':
      return dCrossNormalization * x * y;
    case 'd_xz':
      return dCrossNormalization * x * z;
    case 'd_yz':
      return dCrossNormalization * y * z;
    case 'd_x2_y2':
      return Math.sqrt(15 / (16 * Math.PI)) * (x * x - y * y);
    case 'd_z2':
      return Math.sqrt(5 / (16 * Math.PI)) * (3 * z * z - 1);
  }
}

function angularOverlap(leftName: RealOrbitalName, rightName: RealOrbitalName): Complex {
  const component = (part: 're' | 'im'): number =>
    integrateSolidAngle(
      (thetaRadians, phiRadians) => {
        const left = realSphericalHarmonic(leftName, thetaRadians, phiRadians);
        const right = realSphericalHarmonic(rightName, thetaRadians, phiRadians);
        return multiplyComplex(conjugateComplex(left), right)[part];
      },
      ANGULAR_THETA_INTERVALS,
      ANGULAR_PHI_INTERVALS,
    );

  return { re: component('re'), im: component('im') };
}

function distinctPairs(
  names: readonly RealOrbitalName[],
): readonly (readonly [RealOrbitalName, RealOrbitalName])[] {
  const pairs: [RealOrbitalName, RealOrbitalName][] = [];
  for (let leftIndex = 0; leftIndex < names.length; leftIndex += 1) {
    for (let rightIndex = leftIndex + 1; rightIndex < names.length; rightIndex += 1) {
      const left = names[leftIndex];
      const right = names[rightIndex];
      if (left !== undefined && right !== undefined) pairs.push([left, right]);
    }
  }
  return pairs;
}

describe('orbitales réelles standard', () => {
  it('expose les métadonnées de combinaison et les n minimaux attendus', () => {
    expect(REAL_ORBITAL_DEFINITIONS).toEqual({
      p_x: { absoluteOrder: 1, kind: 'cosine', l: 1, minimumPrincipalQuantumNumber: 2 },
      p_y: { absoluteOrder: 1, kind: 'sine', l: 1, minimumPrincipalQuantumNumber: 2 },
      p_z: { absoluteOrder: 0, kind: 'zonal', l: 1, minimumPrincipalQuantumNumber: 2 },
      d_xy: { absoluteOrder: 2, kind: 'sine', l: 2, minimumPrincipalQuantumNumber: 3 },
      d_xz: { absoluteOrder: 1, kind: 'cosine', l: 2, minimumPrincipalQuantumNumber: 3 },
      d_yz: { absoluteOrder: 1, kind: 'sine', l: 2, minimumPrincipalQuantumNumber: 3 },
      d_x2_y2: { absoluteOrder: 2, kind: 'cosine', l: 2, minimumPrincipalQuantumNumber: 3 },
      d_z2: { absoluteOrder: 0, kind: 'zonal', l: 2, minimumPrincipalQuantumNumber: 3 },
    });
  });

  it('rejette les noms inconnus et les nombres principaux sous le minimum orbital', () => {
    const invalidName = 'f_xyz' as RealOrbitalName;
    expect(() => realSphericalHarmonic(invalidName, 1, 1)).toThrow(RangeError);
    expect(() => realOrbitalWavefunction(4, invalidName, 1, 1, 1)).toThrow(RangeError);

    for (const name of P_ORBITALS) {
      expect(() => realOrbitalWavefunction(1, name, 1, 1, 1)).toThrow(RangeError);
    }
    for (const name of D_ORBITALS) {
      expect(() => realOrbitalWavefunction(2, name, 1, 1, 1)).toThrow(RangeError);
    }
    expect(() => realOrbitalWavefunction(2.5, 'p_x', 1, 1, 1)).toThrow(RangeError);
    expect(() => realOrbitalWavefunction(Number.NaN, 'p_x', 1, 1, 1)).toThrow(RangeError);
  });

  it.each(REAL_ORBITALS)('respecte l’oracle cartésien signé de %s', (name) => {
    const thetaRadians = 1.1;
    const phiRadians = 0.7;
    expectComplexComponents(
      realSphericalHarmonic(name, thetaRadians, phiRadians),
      signedCartesianOracle(name, thetaRadians, phiRadians),
    );
  });

  it('oriente p_x, p_y et p_z sur leurs axes et plans nodaux respectifs', () => {
    const xAxis = { theta: Math.PI / 2, phi: 0 };
    const yAxis = { theta: Math.PI / 2, phi: Math.PI / 2 };
    const zAxis = { theta: 0, phi: 0 };

    expect(realSphericalHarmonic('p_x', xAxis.theta, xAxis.phi).re).toBeGreaterThan(0);
    expect(realSphericalHarmonic('p_y', yAxis.theta, yAxis.phi).re).toBeGreaterThan(0);
    expect(realSphericalHarmonic('p_z', zAxis.theta, zAxis.phi).re).toBeGreaterThan(0);
    expectComplexComponents(realSphericalHarmonic('p_x', yAxis.theta, yAxis.phi), 0);
    expectComplexComponents(realSphericalHarmonic('p_y', xAxis.theta, xAxis.phi), 0);
    expectComplexComponents(realSphericalHarmonic('p_z', xAxis.theta, xAxis.phi), 0);
  });

  it('respecte les plans nodaux des orbitales d_xy, d_xz, d_yz et d_x2_y2', () => {
    expect(realSphericalHarmonic('d_xy', Math.PI / 2, Math.PI / 4).re).toBeGreaterThan(0);
    expectComplexComponents(realSphericalHarmonic('d_xy', Math.PI / 2, 0), 0);
    expectComplexComponents(realSphericalHarmonic('d_xy', Math.PI / 2, Math.PI / 2), 0);

    expect(realSphericalHarmonic('d_xz', Math.PI / 4, 0).re).toBeGreaterThan(0);
    expectComplexComponents(realSphericalHarmonic('d_xz', Math.PI / 2, 0), 0);
    expectComplexComponents(realSphericalHarmonic('d_xz', Math.PI / 4, Math.PI / 2), 0);

    expect(realSphericalHarmonic('d_yz', Math.PI / 4, Math.PI / 2).re).toBeGreaterThan(0);
    expectComplexComponents(realSphericalHarmonic('d_yz', Math.PI / 2, Math.PI / 2), 0);
    expectComplexComponents(realSphericalHarmonic('d_yz', Math.PI / 4, 0), 0);

    expect(realSphericalHarmonic('d_x2_y2', Math.PI / 2, 0).re).toBeGreaterThan(0);
    expectComplexComponents(realSphericalHarmonic('d_x2_y2', Math.PI / 2, Math.PI / 4), 0);
    expectComplexComponents(realSphericalHarmonic('d_x2_y2', Math.PI / 2, (3 * Math.PI) / 4), 0);
  });

  it('place le cône nodal de d_z2 à cos²(theta)=1/3', () => {
    const coneTheta = Math.acos(1 / Math.sqrt(3));
    expect(realSphericalHarmonic('d_z2', 0, 0).re).toBeGreaterThan(0);
    expectComplexComponents(realSphericalHarmonic('d_z2', coneTheta, 0.37), 0);
    expectComplexComponents(realSphericalHarmonic('d_z2', Math.PI - coneTheta, 1.2), 0);
  });

  it('maintient une partie imaginaire résiduelle au niveau de l’arrondi', () => {
    for (const name of REAL_ORBITALS) {
      const angular = realSphericalHarmonic(name, 1.234, 2.345);
      expect(Math.abs(angular.im)).toBeLessThanOrEqual(ANALYTIC_COMPONENT_TOLERANCE);
      const n = REAL_ORBITAL_DEFINITIONS[name].minimumPrincipalQuantumNumber;
      const wavefunction = realOrbitalWavefunction(n, name, 2.1, 1.234, 2.345);
      expect(Math.abs(wavefunction.im)).toBeLessThanOrEqual(ANALYTIC_COMPONENT_TOLERANCE);
    }
  });

  it.each(REAL_ORBITALS)('normalise la partie angulaire de %s', (name) => {
    const normalization = integrateSolidAngle(
      (thetaRadians, phiRadians) =>
        complexModulusSquared(realSphericalHarmonic(name, thetaRadians, phiRadians)),
      ANGULAR_THETA_INTERVALS,
      ANGULAR_PHI_INTERVALS,
    );
    expect(Math.abs(normalization - 1)).toBeLessThanOrEqual(
      NUMERICAL_NORMALIZATION_ABSOLUTE_TOLERANCE,
    );
  });

  it.each(REAL_ORBITALS)('normalise en trois dimensions l’orbitale %s', (name) => {
    const n = REAL_ORBITAL_DEFINITIONS[name].minimumPrincipalQuantumNumber;
    const maximumRadiusBohr =
      RADIAL_EXTENT_IN_N_SQUARED_A_MU * n * n * HYDROGEN_CHARACTERISTIC_RADIUS_BOHR;
    const normalization = integrateSphericalVolume(
      (rBohr, thetaRadians, phiRadians) =>
        complexModulusSquared(realOrbitalWavefunction(n, name, rBohr, thetaRadians, phiRadians)),
      maximumRadiusBohr,
      VOLUME_RADIAL_INTERVALS,
      VOLUME_THETA_INTERVALS,
      VOLUME_PHI_INTERVALS,
    );
    expect(Math.abs(normalization - 1)).toBeLessThanOrEqual(
      THREE_DIMENSIONAL_NORMALIZATION_ABSOLUTE_TOLERANCE,
    );
  });

  it.each([...distinctPairs(P_ORBITALS), ...distinctPairs(D_ORBITALS)])(
    'rend %s et %s orthogonales',
    (leftName, rightName) => {
      const overlap = angularOverlap(leftName, rightName);
      expect(Math.abs(overlap.re)).toBeLessThanOrEqual(NUMERICAL_ORTHOGONALITY_ABSOLUTE_TOLERANCE);
      expect(Math.abs(overlap.im)).toBeLessThanOrEqual(NUMERICAL_ORTHOGONALITY_ABSOLUTE_TOLERANCE);
    },
  );

  it.each([
    [3, P_ORBITALS],
    [4, D_ORBITALS],
  ] as const)('partage le même facteur radial pour n=%i au sein de la famille', (n, names) => {
    for (const rBohr of [1.5, 5]) {
      const l = REAL_ORBITAL_DEFINITIONS[names[0]].l;
      const expectedRadialDensity = hydrogenRadialWavefunction(n, l, rBohr) ** 2;
      for (const name of names) {
        const angularlyIntegratedDensity = integrateSolidAngle(
          (thetaRadians, phiRadians) =>
            complexModulusSquared(
              realOrbitalWavefunction(n, name, rBohr, thetaRadians, phiRadians),
            ),
          ANGULAR_THETA_INTERVALS,
          ANGULAR_PHI_INTERVALS,
        );
        expect(Math.abs(angularlyIntegratedDensity - expectedRadialDensity)).toBeLessThanOrEqual(
          NUMERICAL_NORMALIZATION_ABSOLUTE_TOLERANCE,
        );
      }
    }
  });
});
