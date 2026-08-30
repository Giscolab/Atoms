import { describe, expect, it } from 'vitest';

import { HYDROGEN_CHARACTERISTIC_RADIUS_BOHR } from '../../src/science/hydrogen/radialWavefunction';
import {
  tabulateAngularDensityCut,
  tabulateRadialDistribution,
} from '../../src/sampling/orbitalCharts';

function maximumIndex(values: Float64Array): number {
  let result = 0;
  for (let index = 1; index < values.length; index += 1) {
    const candidate = values[index];
    const current = values[result];
    if (candidate === undefined || current === undefined) {
      throw new RangeError('La série de test ne doit contenir aucune valeur absente.');
    }
    if (candidate > current) result = index;
  }
  return result;
}

describe('données scientifiques des graphiques orbitaux', () => {
  it('tabule P_10(r), conserve son maximum brut et normalise seulement l’affichage', () => {
    const characteristicRadius = HYDROGEN_CHARACTERISTIC_RADIUS_BOHR;
    const result = tabulateRadialDistribution(
      { basis: 'complex', n: 1, l: 0, m: 0 },
      4 * characteristicRadius,
      401,
    );

    expect(result.radiiBohr).toHaveLength(401);
    expect(result.radiiBohr[0]).toBe(0);
    expect(result.radiiBohr.at(-1)).toBeCloseTo(4 * characteristicRadius, 14);
    expect(result.rawDensityInverseBohr[0]).toBe(0);
    expect(result.maximumRawDensityInverseBohr).toBeGreaterThan(0);
    expect(Math.max(...result.normalizedDensity)).toBe(1);

    const peakIndex = maximumIndex(result.rawDensityInverseBohr);
    expect(result.radiiBohr[peakIndex]).toBeCloseTo(characteristicRadius, 12);
    expect(result.rawDensityInverseBohr[peakIndex]).toBe(result.maximumRawDensityInverseBohr);
  });

  it('retient déterministement xz pour p_z et montre ses deux lobes polaires', () => {
    const result = tabulateAngularDensityCut({ basis: 'real', n: 2, orbital: 'p_z' }, 361, 720);

    expect(result.plane).toBe('xz');
    expect(result.planeLabel).toBe('Plan xz');
    expect(result.normalizedRadius[0]).toBeLessThan(1e-28);
    expect(result.normalizedRadius[90]).toBeCloseTo(1, 14);
    expect(result.normalizedRadius[180]).toBeLessThan(1e-28);
    expect(result.normalizedRadius[270]).toBeCloseTo(1, 14);
    expect(result.anglesRadians[0]).toBe(0);
    expect(result.anglesRadians.at(-1)).toBe(2 * Math.PI);
    expect(result.rawDensityInverseSteradians.at(-1)).toBe(result.rawDensityInverseSteradians[0]);
  });

  it('retient xy pour d_xy, avec des nœuds sur les axes et des lobes diagonaux', () => {
    const result = tabulateAngularDensityCut({ basis: 'real', n: 3, orbital: 'd_xy' }, 361, 720);

    expect(result.plane).toBe('xy');
    expect(result.normalizedRadius[0]).toBeLessThan(1e-28);
    expect(result.normalizedRadius[45]).toBeCloseTo(1, 14);
    expect(result.normalizedRadius[90]).toBeLessThan(1e-28);
    expect(result.normalizedRadius[135]).toBeCloseTo(1, 14);
  });

  it('choisit xy par ordre stable pour la coupe isotrope de 1s', () => {
    const result = tabulateAngularDensityCut({ basis: 'complex', n: 1, l: 0, m: 0 }, 73, 144);

    expect(result.plane).toBe('xy');
    for (const radius of result.normalizedRadius) expect(radius).toBeCloseTo(1, 14);
  });

  it('rend la distribution radiale indépendante de m et de l’orientation réelle', () => {
    const commonArguments = [30, 301] as const;
    const negativeOrder = tabulateRadialDistribution(
      { basis: 'complex', n: 3, l: 2, m: -2 },
      ...commonArguments,
    );
    const positiveOrder = tabulateRadialDistribution(
      { basis: 'complex', n: 3, l: 2, m: 1 },
      ...commonArguments,
    );
    const realOrientation = tabulateRadialDistribution(
      { basis: 'real', n: 3, orbital: 'd_xy' },
      ...commonArguments,
    );

    expect(positiveOrder.rawDensityInverseBohr).toEqual(negativeOrder.rawDensityInverseBohr);
    expect(realOrientation.rawDensityInverseBohr).toEqual(negativeOrder.rawDensityInverseBohr);
    expect(positiveOrder.normalizedDensity).toEqual(negativeOrder.normalizedDensity);
    expect(realOrientation.normalizedDensity).toEqual(negativeOrder.normalizedDensity);
  });

  it('rejette les étendues et résolutions impropres à une tabulation', () => {
    const state = { basis: 'complex', n: 1, l: 0, m: 0 } as const;
    expect(() => tabulateRadialDistribution(state, 0, 101)).toThrow(RangeError);
    expect(() => tabulateRadialDistribution(state, 10, 2)).toThrow(RangeError);
    expect(() => tabulateAngularDensityCut(state, 2, 360)).toThrow(RangeError);
    expect(() => tabulateAngularDensityCut(state, 360, 3)).toThrow(RangeError);
  });
});
