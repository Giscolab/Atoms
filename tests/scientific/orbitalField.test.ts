import { describe, expect, it } from 'vitest';

import {
  realOrbitalWavefunction,
  type RealOrbitalName,
} from '../../src/science/hydrogen/realOrbitals';
import { hydrogenWavefunction } from '../../src/science/hydrogen/wavefunction';
import { complexModulusSquared, complexPhase } from '../../src/science/math/complex';
import { ORBITAL_FIELD_VERSION } from '../../src/sampling/fieldContracts';
import {
  computeOrbitalField,
  evaluateOrbitalWavefunction,
  MAX_ORBITAL_FIELD_RESOLUTION,
  MIN_ORBITAL_FIELD_RESOLUTION,
  orbitalDegree,
} from '../../src/sampling/orbitalField';
import type { OrbitalSamplingState } from '../../src/sampling/contracts';

function gridCoordinate(index: number, resolution: number, extentBohr: number): number {
  return extentBohr * (2 * (index / (resolution - 1)) - 1);
}

function fieldIndex(xIndex: number, yIndex: number, zIndex: number, resolution: number): number {
  return xIndex + resolution * (yIndex + resolution * zIndex);
}

function directCartesianAmplitude(
  state: OrbitalSamplingState,
  xIndex: number,
  yIndex: number,
  zIndex: number,
  resolution: number,
  extentBohr: number,
) {
  const xBohr = gridCoordinate(xIndex, resolution, extentBohr);
  const yBohr = gridCoordinate(yIndex, resolution, extentBohr);
  const zBohr = gridCoordinate(zIndex, resolution, extentBohr);
  const radiusBohr = Math.hypot(xBohr, yBohr, zBohr);
  const thetaRadians =
    radiusBohr === 0 ? 0 : Math.acos(Math.max(-1, Math.min(1, zBohr / radiusBohr)));
  const phiRadians = radiusBohr === 0 ? 0 : Math.atan2(yBohr, xBohr);
  return evaluateOrbitalWavefunction(state, radiusBohr, thetaRadians, phiRadians);
}

describe('champ orbital cartésien pur', () => {
  it('produit trois buffers N³ normalisés et conserve les maxima bruts de ψ_100', () => {
    const resolution = 9;
    const extentBohr = 3;
    const field = computeOrbitalField(
      { basis: 'complex', n: 1, l: 0, m: 0 },
      { extentBohr, resolution },
    );
    const expectedLength = resolution ** 3;
    const centerIndex = fieldIndex(4, 4, 4, resolution);
    const centerAmplitude = hydrogenWavefunction(1, 0, 0, 0, 0, 0);

    expect(field.fieldVersion).toBe(ORBITAL_FIELD_VERSION);
    expect(field.extentBohr).toBe(extentBohr);
    expect(field.resolution).toBe(resolution);
    expect(field.densityNormalized).toHaveLength(expectedLength);
    expect(field.phaseRadians).toHaveLength(expectedLength);
    expect(field.signedAmplitudeNormalized).toHaveLength(expectedLength);
    expect(field.densityNormalized.buffer).not.toBe(field.phaseRadians.buffer);
    expect(field.phaseRadians.buffer).not.toBe(field.signedAmplitudeNormalized.buffer);
    expect(field.maximumDensityPerCubicBohr).toBeCloseTo(
      complexModulusSquared(centerAmplitude),
      13,
    );
    expect(field.maximumWavefunctionAmplitude).toBeCloseTo(
      Math.hypot(centerAmplitude.re, centerAmplitude.im),
      13,
    );
    expect(field.maximumWavefunctionAmplitude ** 2).toBeCloseTo(
      field.maximumDensityPerCubicBohr,
      13,
    );
    expect(field.densityNormalized[centerIndex]).toBe(1);
    expect(field.signedAmplitudeNormalized[centerIndex]).toBe(1);
    expect(field.phaseRadians[centerIndex]).toBe(0);
    expect(field.nodesAvailable).toBe(true);

    let maximumNormalizedDensity = 0;
    for (let index = 0; index < expectedLength; index += 1) {
      const density = field.densityNormalized[index];
      const signedAmplitude = field.signedAmplitudeNormalized[index];
      const phase = field.phaseRadians[index];
      expect(density).toBeDefined();
      expect(signedAmplitude).toBeDefined();
      expect(phase).toBeDefined();
      if (density === undefined || signedAmplitude === undefined || phase === undefined) continue;
      expect(Number.isFinite(density)).toBe(true);
      expect(density).toBeGreaterThanOrEqual(0);
      expect(density).toBeLessThanOrEqual(1);
      expect(Number.isFinite(signedAmplitude)).toBe(true);
      expect(signedAmplitude).toBeGreaterThanOrEqual(-1);
      expect(signedAmplitude).toBeLessThanOrEqual(1);
      expect(Number.isNaN(phase) || (phase >= -Math.PI && phase <= Math.PI)).toBe(true);
      maximumNormalizedDensity = Math.max(maximumNormalizedDensity, density);
    }
    expect(maximumNormalizedDensity).toBe(1);
  });

  it('respecte exactement l’indexation x-fastest et la convention polaire z', () => {
    const state = { basis: 'complex', n: 2, l: 1, m: 1 } as const;
    const resolution = 8;
    const extentBohr = 4;
    const field = computeOrbitalField(state, { extentBohr, resolution });

    for (const [xIndex, yIndex, zIndex] of [
      [6, 2, 5],
      [7, 2, 5],
      [6, 3, 5],
      [6, 2, 6],
    ] as const) {
      const index = fieldIndex(xIndex, yIndex, zIndex, resolution);
      const amplitude = directCartesianAmplitude(
        state,
        xIndex,
        yIndex,
        zIndex,
        resolution,
        extentBohr,
      );
      const expectedPhase = complexPhase(amplitude);
      expect(expectedPhase).not.toBeNull();
      if (expectedPhase === null) throw new Error('Le point de contrôle doit être hors nœud.');

      expect(field.densityNormalized[index]).toBeCloseTo(
        complexModulusSquared(amplitude) / field.maximumDensityPerCubicBohr,
        6,
      );
      expect(field.signedAmplitudeNormalized[index]).toBeCloseTo(
        amplitude.re / field.maximumWavefunctionAmplitude,
        6,
      );
      expect(field.phaseRadians[index]).toBeCloseTo(expectedPhase, 6);
    }
    expect(field.nodesAvailable).toBe(false);
  });

  it('marque la phase indéfinie au nœud exact et préserve les lobes signés réels', () => {
    const resolution = 9;
    const center = 4;
    const field = computeOrbitalField(
      { basis: 'real', n: 2, orbital: 'p_x' },
      { extentBohr: 6, resolution },
    );
    const originIndex = fieldIndex(center, center, center, resolution);
    const negativeXIndex = fieldIndex(center - 1, center, center, resolution);
    const positiveXIndex = fieldIndex(center + 1, center, center, resolution);
    const negativeLobe = field.signedAmplitudeNormalized[negativeXIndex];
    const positiveLobe = field.signedAmplitudeNormalized[positiveXIndex];

    expect(field.nodesAvailable).toBe(true);
    expect(field.densityNormalized[originIndex]).toBe(0);
    expect(field.signedAmplitudeNormalized[originIndex]).toBe(0);
    expect(Number.isNaN(field.phaseRadians[originIndex])).toBe(true);
    expect(negativeLobe).toBeDefined();
    expect(positiveLobe).toBeDefined();
    if (negativeLobe === undefined || positiveLobe === undefined) {
      throw new Error('Les deux lobes doivent appartenir à la grille.');
    }
    expect(negativeLobe * positiveLobe).toBeLessThan(0);
    expect(Math.abs(negativeLobe)).toBeCloseTo(Math.abs(positiveLobe), 6);
  });

  it('n’annonce des nœuds ψ=0 que pour une fonction choisie réelle', () => {
    expect(
      computeOrbitalField({ basis: 'complex', n: 2, l: 1, m: 0 }, { extentBohr: 6, resolution: 8 })
        .nodesAvailable,
    ).toBe(true);
    expect(
      computeOrbitalField({ basis: 'complex', n: 2, l: 1, m: -1 }, { extentBohr: 6, resolution: 8 })
        .nodesAvailable,
    ).toBe(false);
    expect(
      computeOrbitalField(
        { basis: 'real', n: 3, orbital: 'd_xy' },
        { extentBohr: 10, resolution: 8 },
      ).nodesAvailable,
    ).toBe(true);
  });

  it('délègue l’amplitude et le degré aux contrats scientifiques de Phase 3', () => {
    const complexState = { basis: 'complex', n: 3, l: 2, m: -1 } as const;
    const realState = { basis: 'real', n: 3, orbital: 'd_xz' } as const;
    const coordinates = [2.7, 1.1, -0.6] as const;

    expect(orbitalDegree(complexState)).toBe(2);
    expect(orbitalDegree(realState)).toBe(2);
    expect(evaluateOrbitalWavefunction(complexState, ...coordinates)).toEqual(
      hydrogenWavefunction(3, 2, -1, ...coordinates),
    );
    expect(evaluateOrbitalWavefunction(realState, ...coordinates)).toEqual(
      realOrbitalWavefunction(3, 'd_xz', ...coordinates),
    );
  });

  it('rejette les limites numériques, étendues et états invalides', () => {
    const validState = { basis: 'complex', n: 1, l: 0, m: 0 } as const;
    const invalidRealState = {
      basis: 'real',
      n: 2,
      orbital: 'inconnue' as RealOrbitalName,
    } as const;
    const invalidCalls: ReadonlyArray<() => unknown> = [
      () =>
        computeOrbitalField(validState, {
          extentBohr: 2,
          resolution: MIN_ORBITAL_FIELD_RESOLUTION - 1,
        }),
      () =>
        computeOrbitalField(validState, {
          extentBohr: 2,
          resolution: MAX_ORBITAL_FIELD_RESOLUTION + 1,
        }),
      () => computeOrbitalField(validState, { extentBohr: 2, resolution: 8.5 }),
      () => computeOrbitalField(validState, { extentBohr: 0, resolution: 8 }),
      () => computeOrbitalField(validState, { extentBohr: -1, resolution: 8 }),
      () => computeOrbitalField(validState, { extentBohr: Number.NaN, resolution: 8 }),
      () =>
        computeOrbitalField(validState, {
          extentBohr: Number.POSITIVE_INFINITY,
          resolution: 8,
        }),
      () => computeOrbitalField(validState, { extentBohr: Number.MAX_VALUE, resolution: 8 }),
      () =>
        computeOrbitalField(
          { basis: 'complex', n: 2, l: 2, m: 0 },
          { extentBohr: 2, resolution: 8 },
        ),
      () => computeOrbitalField(invalidRealState, { extentBohr: 2, resolution: 8 }),
      () =>
        computeOrbitalField(
          { basis: 'real', n: 1, orbital: 'p_x' },
          { extentBohr: 2, resolution: 8 },
        ),
    ];

    for (const invalidCall of invalidCalls) expect(invalidCall).toThrow(RangeError);
  });
});
