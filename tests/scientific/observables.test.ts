import { describe, expect, it } from 'vitest';

import {
  hydrogenExpectedInverseRadiusPerBohr,
  hydrogenExpectedRadiusBohr,
  hydrogenMostProbableRadius1sBohr,
  hydrogenNodeCounts,
} from '../../src/science/hydrogen/observables';
import {
  HYDROGEN_CHARACTERISTIC_RADIUS_BOHR,
  hydrogenRadialWavefunction,
} from '../../src/science/hydrogen/radialWavefunction';
import { compositeSimpson } from './numericalIntegration';
import {
  ANALYTIC_COMPONENT_TOLERANCE,
  NUMERICAL_NORMALIZATION_ABSOLUTE_TOLERANCE,
} from './numericAssertions';

const RADIAL_INTERVALS = 16_384;

function integrateRadialMoment(n: number, l: number, power: number): number {
  const maximumRadius = 20 * n * n * HYDROGEN_CHARACTERISTIC_RADIUS_BOHR;
  return compositeSimpson(
    (rBohr) => {
      const radial = hydrogenRadialWavefunction(n, l, rBohr);
      return rBohr ** (power + 2) * radial * radial;
    },
    0,
    maximumRadius,
    RADIAL_INTERVALS,
  );
}

describe('nœuds et observables analytiques de ¹H', () => {
  it.each([
    [1, 0, 0, 0, 0],
    [2, 0, 1, 0, 1],
    [2, 1, 0, 1, 1],
    [3, 0, 2, 0, 2],
    [3, 1, 1, 1, 2],
    [3, 2, 0, 2, 2],
  ])('centralise les nœuds de l’état (%i,%i)', (n, l, radialNodes, angularNodes, totalNodes) => {
    expect(hydrogenNodeCounts(n, l)).toEqual({ radialNodes, angularNodes, totalNodes });
  });

  it('retrouve les valeurs analytiques simples en unités de a0', () => {
    const scale = HYDROGEN_CHARACTERISTIC_RADIUS_BOHR;
    expect(Math.abs(hydrogenExpectedRadiusBohr(1, 0) - (3 * scale) / 2)).toBeLessThanOrEqual(
      ANALYTIC_COMPONENT_TOLERANCE,
    );
    expect(Math.abs(hydrogenExpectedRadiusBohr(2, 0) - 6 * scale)).toBeLessThanOrEqual(
      ANALYTIC_COMPONENT_TOLERANCE,
    );
    expect(Math.abs(hydrogenExpectedRadiusBohr(2, 1) - 5 * scale)).toBeLessThanOrEqual(
      ANALYTIC_COMPONENT_TOLERANCE,
    );
    expect(Math.abs(hydrogenExpectedInverseRadiusPerBohr(1) - 1 / scale)).toBeLessThanOrEqual(
      ANALYTIC_COMPONENT_TOLERANCE,
    );
    expect(hydrogenMostProbableRadius1sBohr()).toBe(scale);
  });

  it.each([
    [1, 0],
    [2, 0],
    [2, 1],
  ])('compare <r> analytique à l’intégrale radiale de l’état (%i,%i)', (n, l) => {
    const numerical = integrateRadialMoment(n, l, 1);
    expect(Math.abs(numerical - hydrogenExpectedRadiusBohr(n, l))).toBeLessThanOrEqual(
      NUMERICAL_NORMALIZATION_ABSOLUTE_TOLERANCE,
    );
  });

  it.each([
    [1, 0],
    [2, 0],
    [2, 1],
  ])('compare <1/r> analytique à l’intégrale radiale de l’état (%i,%i)', (n, l) => {
    const numerical = integrateRadialMoment(n, l, -1);
    expect(Math.abs(numerical - hydrogenExpectedInverseRadiusPerBohr(n))).toBeLessThanOrEqual(
      NUMERICAL_NORMALIZATION_ABSOLUTE_TOLERANCE,
    );
  });

  it('confirme que la distribution radiale 1s culmine à a_mu', () => {
    const mode = hydrogenMostProbableRadius1sBohr();
    const radialProbability = (rBohr: number): number => {
      const radial = hydrogenRadialWavefunction(1, 0, rBohr);
      return rBohr * rBohr * radial * radial;
    };
    const displacement = mode / 100;
    expect(radialProbability(mode)).toBeGreaterThan(radialProbability(mode - displacement));
    expect(radialProbability(mode)).toBeGreaterThan(radialProbability(mode + displacement));
  });

  it('rejette les nombres quantiques invalides', () => {
    expect(() => hydrogenNodeCounts(0, 0)).toThrow(RangeError);
    expect(() => hydrogenNodeCounts(2, 2)).toThrow(RangeError);
    expect(() => hydrogenExpectedRadiusBohr(2, -1)).toThrow(RangeError);
    expect(() => hydrogenExpectedInverseRadiusPerBohr(Number.NaN)).toThrow(RangeError);
  });
});
