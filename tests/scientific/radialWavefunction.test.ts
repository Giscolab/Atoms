import { describe, expect, it } from 'vitest';

import {
  HYDROGEN_CHARACTERISTIC_RADIUS_BOHR,
  hydrogenRadialCoordinate,
  hydrogenRadialWavefunction,
} from '../../src/science/hydrogen/radialWavefunction';
import { compositeSimpson } from './numericalIntegration';
import {
  ANALYTIC_COMPONENT_TOLERANCE,
  NUMERICAL_NORMALIZATION_ABSOLUTE_TOLERANCE,
  NUMERICAL_ORTHOGONALITY_ABSOLUTE_TOLERANCE,
} from './numericAssertions';

const RADIAL_INTERVALS = 16_384;
const COARSE_RADIAL_INTERVALS = RADIAL_INTERVALS / 2;
const RADIAL_EXTENT_IN_SCALES = 20;

function radialMaximumBohr(n: number, extent = RADIAL_EXTENT_IN_SCALES): number {
  return extent * n * n * HYDROGEN_CHARACTERISTIC_RADIUS_BOHR;
}

function radialNormalization(n: number, l: number, intervals = RADIAL_INTERVALS): number {
  return compositeSimpson(
    (rBohr) => {
      const radial = hydrogenRadialWavefunction(n, l, rBohr);
      return rBohr * rBohr * radial * radial;
    },
    0,
    radialMaximumBohr(n),
    intervals,
  );
}

function radialOverlap(
  leftN: number,
  rightN: number,
  l: number,
  intervals = RADIAL_INTERVALS,
): number {
  return compositeSimpson(
    (rBohr) =>
      rBohr *
      rBohr *
      hydrogenRadialWavefunction(leftN, l, rBohr) *
      hydrogenRadialWavefunction(rightN, l, rBohr),
    0,
    radialMaximumBohr(Math.max(leftN, rightN)),
    intervals,
  );
}

describe('fonction radiale hydrogène avec masse réduite', () => {
  it('exprime explicitement a_mu et rho dans le contrat spatial en a0', () => {
    const scale = HYDROGEN_CHARACTERISTIC_RADIUS_BOHR;
    expect(scale).toBeGreaterThan(1);
    expect(Math.abs(hydrogenRadialCoordinate(1, scale) - 2)).toBeLessThanOrEqual(
      ANALYTIC_COMPONENT_TOLERANCE,
    );
    expect(Math.abs(hydrogenRadialCoordinate(4, 2 * scale) - 1)).toBeLessThanOrEqual(
      ANALYTIC_COMPONENT_TOLERANCE,
    );
  });

  it('retrouve R_10 à plusieurs rayons', () => {
    const scale = HYDROGEN_CHARACTERISTIC_RADIUS_BOHR;
    for (const rBohr of [0, 0.37 * scale, 2 * scale, 5.5 * scale]) {
      const expected = 2 * scale ** (-3 / 2) * Math.exp(-rBohr / scale);
      expect(Math.abs(hydrogenRadialWavefunction(1, 0, rBohr) - expected)).toBeLessThanOrEqual(
        ANALYTIC_COMPONENT_TOLERANCE,
      );
    }
  });

  it('retrouve R_20 à plusieurs rayons, y compris son nœud', () => {
    const scale = HYDROGEN_CHARACTERISTIC_RADIUS_BOHR;
    for (const rBohr of [0, 0.83 * scale, 2 * scale, 4.2 * scale]) {
      const expected =
        (1 / (2 * Math.sqrt(2))) *
        scale ** (-3 / 2) *
        (2 - rBohr / scale) *
        Math.exp(-rBohr / (2 * scale));
      expect(Math.abs(hydrogenRadialWavefunction(2, 0, rBohr) - expected)).toBeLessThanOrEqual(
        ANALYTIC_COMPONENT_TOLERANCE,
      );
    }
  });

  it('retrouve R_21 à plusieurs rayons', () => {
    const scale = HYDROGEN_CHARACTERISTIC_RADIUS_BOHR;
    for (const rBohr of [0, 0.51 * scale, 2.3 * scale, 6 * scale]) {
      const expected =
        (1 / (2 * Math.sqrt(6))) *
        scale ** (-3 / 2) *
        (rBohr / scale) *
        Math.exp(-rBohr / (2 * scale));
      expect(Math.abs(hydrogenRadialWavefunction(2, 1, rBohr) - expected)).toBeLessThanOrEqual(
        ANALYTIC_COMPONENT_TOLERANCE,
      );
    }
  });

  it.each([
    [1, 0],
    [2, 0],
    [2, 1],
    [3, 0],
    [3, 1],
    [3, 2],
    [4, 3],
  ])('normalise numériquement R_%i%i avec la mesure r² dr', (n, l) => {
    expect(Math.abs(radialNormalization(n, l) - 1)).toBeLessThanOrEqual(
      NUMERICAL_NORMALIZATION_ABSOLUTE_TOLERANCE,
    );
  });

  it('établit la convergence en résolution et en extension du domaine radial', () => {
    const n = 4;
    const l = 0;
    const coarse = radialNormalization(n, l, COARSE_RADIAL_INTERVALS);
    const fine = radialNormalization(n, l);
    const extendedIntervals = 19_660;
    const extended = compositeSimpson(
      (rBohr) => {
        const radial = hydrogenRadialWavefunction(n, l, rBohr);
        return rBohr * rBohr * radial * radial;
      },
      0,
      radialMaximumBohr(n, 24),
      extendedIntervals,
    );

    expect(Math.abs(fine - coarse)).toBeLessThanOrEqual(NUMERICAL_NORMALIZATION_ABSOLUTE_TOLERANCE);
    expect(Math.abs(extended - fine)).toBeLessThanOrEqual(
      NUMERICAL_NORMALIZATION_ABSOLUTE_TOLERANCE,
    );
    expect(Math.abs(extended - 1)).toBeLessThanOrEqual(NUMERICAL_NORMALIZATION_ABSOLUTE_TOLERANCE);
  });

  it('vérifie l’orthogonalité radiale des états de même l', () => {
    expect(Math.abs(radialOverlap(1, 2, 0))).toBeLessThanOrEqual(
      NUMERICAL_ORTHOGONALITY_ABSOLUTE_TOLERANCE,
    );
    expect(Math.abs(radialOverlap(2, 3, 1))).toBeLessThanOrEqual(
      NUMERICAL_ORTHOGONALITY_ABSOLUTE_TOLERANCE,
    );
  });

  it('retrouve numériquement les nœuds radiaux analytiques de 2s, 3s et 3p', () => {
    const scale = HYDROGEN_CHARACTERISTIC_RADIUS_BOHR;
    const roots3s = [(3 * scale * (3 - Math.sqrt(3))) / 2, (3 * scale * (3 + Math.sqrt(3))) / 2];

    expect(Math.abs(hydrogenRadialWavefunction(2, 0, 2 * scale))).toBeLessThanOrEqual(
      ANALYTIC_COMPONENT_TOLERANCE,
    );
    for (const root of roots3s) {
      expect(Math.abs(hydrogenRadialWavefunction(3, 0, root))).toBeLessThanOrEqual(
        ANALYTIC_COMPONENT_TOLERANCE,
      );
    }
    expect(Math.abs(hydrogenRadialWavefunction(3, 1, 6 * scale))).toBeLessThanOrEqual(
      ANALYTIC_COMPONENT_TOLERANCE,
    );
    expect(hydrogenRadialWavefunction(2, 1, 0)).toBe(0);
  });

  it('rejette explicitement les rayons et nombres quantiques hors domaine', () => {
    expect(() => hydrogenRadialCoordinate(0, 0)).toThrow(RangeError);
    expect(() => hydrogenRadialCoordinate(1, -1)).toThrow(RangeError);
    expect(() => hydrogenRadialCoordinate(1, Number.NaN)).toThrow(RangeError);
    expect(() => hydrogenRadialWavefunction(2, 2, 0)).toThrow(RangeError);
    expect(() => hydrogenRadialWavefunction(2, 0.5, 0)).toThrow(RangeError);
    expect(() => hydrogenRadialWavefunction(1, 0, Number.POSITIVE_INFINITY)).toThrow(RangeError);
    expect(() => hydrogenRadialWavefunction(100, 99, 1)).toThrow(RangeError);
  });
});
