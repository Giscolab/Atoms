import { describe, expect, it } from 'vitest';

import { hydrogenWavefunctionPhase } from '../../src/science/hydrogen/wavefunction';
import type { RealOrbitalName } from '../../src/science/hydrogen/realOrbitals';
import { ORBITAL_SAMPLER_VERSION, type OrbitalSampleSet } from '../../src/sampling/contracts';
import { MAX_ORBITAL_SAMPLE_COUNT, sampleOrbital } from '../../src/sampling/orbitalSampler';
import { SAMPLER_PRNG_VERSION } from '../../src/sampling/rng';
import { phaseDistance } from './numericAssertions';
import { typedArrayValue } from './samplingStatistics';

function expectFiniteSampleBuffers(result: OrbitalSampleSet): void {
  for (const value of result.positionsBohr) expect(Number.isFinite(value)).toBe(true);
  for (const value of result.radiiBohr) expect(Number.isFinite(value)).toBe(true);
  for (const value of result.thetaRadians) expect(Number.isFinite(value)).toBe(true);
  for (const value of result.phiRadians) expect(Number.isFinite(value)).toBe(true);
  for (const value of result.phaseRadians) expect(Number.isFinite(value)).toBe(true);
}

describe('orchestrateur scientifique unique', () => {
  it('produit des buffers explicites en a₀ avec la convention sphérique polaire z', () => {
    const sampleCount = 256;
    const result = sampleOrbital({
      sampleCount,
      seed: 0x1234_abcd,
      state: { basis: 'real', n: 3, orbital: 'd_xy' },
    });
    expect(result.positionsBohr).toHaveLength(3 * sampleCount);
    expect(result.radiiBohr).toHaveLength(sampleCount);
    expect(result.thetaRadians).toHaveLength(sampleCount);
    expect(result.phiRadians).toHaveLength(sampleCount);
    expect(result.phaseRadians).toHaveLength(sampleCount);
    expectFiniteSampleBuffers(result);

    for (let index = 0; index < sampleCount; index += 1) {
      const radius = typedArrayValue(result.radiiBohr, index, 'Le rayon');
      const theta = typedArrayValue(result.thetaRadians, index, "L'angle polaire");
      const phi = typedArrayValue(result.phiRadians, index, "L'angle azimutal");
      const offset = index * 3;
      const tolerance = 2e-5 * Math.max(1, radius);
      expect(
        Math.abs(
          typedArrayValue(result.positionsBohr, offset, 'x') -
            radius * Math.sin(theta) * Math.cos(phi),
        ),
      ).toBeLessThanOrEqual(tolerance);
      expect(
        Math.abs(
          typedArrayValue(result.positionsBohr, offset + 1, 'y') -
            radius * Math.sin(theta) * Math.sin(phi),
        ),
      ).toBeLessThanOrEqual(tolerance);
      expect(
        Math.abs(typedArrayValue(result.positionsBohr, offset + 2, 'z') - radius * Math.cos(theta)),
      ).toBeLessThanOrEqual(tolerance);
    }

    expect(result.metadata).toMatchObject({
      positionUnits: 'bohr',
      prngVersion: SAMPLER_PRNG_VERSION,
      sampleCount,
      samplerVersion: ORBITAL_SAMPLER_VERSION,
      seed: 0x1234_abcd,
      state: { basis: 'real', n: 3, orbital: 'd_xy' },
      undefinedPhaseCount: 0,
    });
  });

  it('reproduit exactement les buffers pour un même état, une même seed et une même version', () => {
    const request = {
      sampleCount: 128,
      seed: 0xcafe_babe,
      state: { basis: 'complex', n: 3, l: 2, m: 1 },
    } as const;
    const first = sampleOrbital(request);
    const second = sampleOrbital(request);
    expect(second.positionsBohr).toEqual(first.positionsBohr);
    expect(second.radiiBohr).toEqual(first.radiiBohr);
    expect(second.thetaRadians).toEqual(first.thetaRadians);
    expect(second.phiRadians).toEqual(first.phiRadians);
    expect(second.phaseRadians).toEqual(first.phaseRadians);
    expect(second.metadata).toEqual(first.metadata);
  });

  it('expose la phase du noyau complexe sans introduire de palette dans le sampler', () => {
    const result = sampleOrbital({
      sampleCount: 128,
      seed: 0xface_0001,
      state: { basis: 'complex', n: 2, l: 1, m: 1 },
    });
    expect('phaseColorsRgb' in result).toBe(false);
    for (let index = 0; index < result.metadata.sampleCount; index += 1) {
      const expected = hydrogenWavefunctionPhase(
        2,
        1,
        1,
        typedArrayValue(result.radiiBohr, index, 'Le rayon'),
        typedArrayValue(result.thetaRadians, index, 'theta'),
        typedArrayValue(result.phiRadians, index, 'phi'),
      );
      expect(expected).not.toBeNull();
      if (expected !== null) {
        expect(
          phaseDistance(typedArrayValue(result.phaseRadians, index, 'La phase'), expected),
        ).toBeLessThanOrEqual(1e-5);
      }
    }
  });

  it('rejette uniformément les requêtes complexes, réelles, seeds et tailles invalides', () => {
    const invalidRealName = 'f_xyz' as RealOrbitalName;
    const invalidCalls: ReadonlyArray<() => unknown> = [
      () =>
        sampleOrbital({ sampleCount: 0, seed: 1, state: { basis: 'complex', n: 1, l: 0, m: 0 } }),
      () =>
        sampleOrbital({ sampleCount: 1.5, seed: 1, state: { basis: 'complex', n: 1, l: 0, m: 0 } }),
      () =>
        sampleOrbital({
          sampleCount: MAX_ORBITAL_SAMPLE_COUNT + 1,
          seed: 1,
          state: { basis: 'complex', n: 1, l: 0, m: 0 },
        }),
      () =>
        sampleOrbital({ sampleCount: 1, seed: -1, state: { basis: 'complex', n: 1, l: 0, m: 0 } }),
      () =>
        sampleOrbital({ sampleCount: 1, seed: 1, state: { basis: 'complex', n: 2, l: 2, m: 0 } }),
      () =>
        sampleOrbital({ sampleCount: 1, seed: 1, state: { basis: 'complex', n: 2, l: 1, m: 2 } }),
      () =>
        sampleOrbital({ sampleCount: 1, seed: 1, state: { basis: 'real', n: 1, orbital: 'p_x' } }),
      () =>
        sampleOrbital({
          sampleCount: 1,
          seed: 1,
          state: { basis: 'real', n: 4, orbital: invalidRealName },
        }),
    ];
    for (const invalidCall of invalidCalls) expect(invalidCall).toThrow(RangeError);
  });
});
