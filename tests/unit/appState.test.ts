import { describe, expect, it } from 'vitest';

import { MAX_ORBITAL_SAMPLE_COUNT } from '../../src/sampling/orbitalSampler';
import {
  createAppState,
  DEFAULT_ORBITAL_SAMPLE_COUNT,
  DEFAULT_ORBITAL_SEED,
  normalizeAppState,
  normalizeOrbitalState,
  normalizeRenderingState,
  normalizeSamplingConfiguration,
} from '../../src/state/appState';

function containsArrayBufferView(value: unknown, visited = new Set<object>()): boolean {
  if (ArrayBuffer.isView(value)) return true;
  if (typeof value !== 'object' || value === null || visited.has(value)) return false;
  visited.add(value);
  return Object.values(value).some((nestedValue) => containsArrayBufferView(nestedValue, visited));
}

describe('état applicatif scientifique de Phase 7', () => {
  it('reproduit intégralement les défauts de la maquette validée sans buffer partagé', () => {
    const state = createAppState();

    expect(state).toEqual({
      legacy: {
        legacy2DInitialized: false,
        showLegacy2D: false,
      },
      orbital: { basis: 'real', n: 3, orbital: 'd_xy' },
      rendering: {
        cameraRotationEnabled: true,
        displayMode: 'hybrid',
        isoDensityFraction: 0.2,
        observable: 'phase',
        pointOpacity: 0.85,
        pointSizePixels: 1.2,
        showAxes: true,
        showNodes: true,
        theme: 'dark',
      },
      sampling: {
        sampleCount: DEFAULT_ORBITAL_SAMPLE_COUNT,
        seed: DEFAULT_ORBITAL_SEED,
      },
    });
    expect(state.sampling.sampleCount).toBe(15_000);
    expect(containsArrayBufferView(state)).toBe(false);
  });

  it('conserve des représentations discriminées réellement distinctes', () => {
    const complex = normalizeOrbitalState({ basis: 'complex', n: 4, l: 2, m: -1 });
    const real = normalizeOrbitalState({ basis: 'real', n: 3, orbital: 'd_z2' });

    expect(complex).toEqual({ basis: 'complex', n: 4, l: 2, m: -1 });
    expect(real).toEqual({ basis: 'real', n: 3, orbital: 'd_z2' });
    expect('orbital' in complex).toBe(false);
    expect('l' in real).toBe(false);
    expect('m' in real).toBe(false);
  });

  it('normalise par copie profonde sans muter ni partager les sous-états', () => {
    const source = createAppState();
    const normalized = normalizeAppState(source);

    expect(normalized).toEqual(source);
    expect(normalized).not.toBe(source);
    expect(normalized.legacy).not.toBe(source.legacy);
    expect(normalized.orbital).not.toBe(source.orbital);
    expect(normalized.rendering).not.toBe(source.rendering);
    expect(normalized.sampling).not.toBe(source.sampling);
  });

  it('accepte chaque variante graphique et le thème clair dans leurs domaines', () => {
    const base = createAppState();
    const normalized = normalizeAppState({
      ...base,
      orbital: { basis: 'complex', n: 2, l: 1, m: 1 },
      rendering: {
        ...base.rendering,
        cameraRotationEnabled: false,
        displayMode: 'isosurface',
        isoDensityFraction: 0.5,
        observable: 'density',
        pointOpacity: 0,
        pointSizePixels: 0.25,
        showAxes: false,
        showNodes: false,
        theme: 'light',
      },
      sampling: { sampleCount: 1, seed: 0xffff_ffff },
    });

    expect(normalized.orbital).toEqual({ basis: 'complex', n: 2, l: 1, m: 1 });
    expect(normalized.rendering).toMatchObject({
      displayMode: 'isosurface',
      observable: 'density',
      theme: 'light',
    });
    expect(normalized.sampling).toEqual({ sampleCount: 1, seed: 0xffff_ffff });
  });

  it('refuse les états quantiques invalides au lieu de les corriger silencieusement', () => {
    const invalidStates: readonly unknown[] = [
      null,
      { basis: 'inconnue', n: 1 },
      { basis: 'complex', n: 0, l: 0, m: 0 },
      { basis: 'complex', n: 2, l: 2, m: 0 },
      { basis: 'complex', n: 2, l: 1, m: 2 },
      { basis: 'real', n: 1, orbital: 'p_x' },
      { basis: 'real', n: 3, orbital: 'inconnue' },
    ];

    for (const invalidState of invalidStates) {
      expect(() => normalizeOrbitalState(invalidState)).toThrow();
    }
  });

  it('valide les domaines techniques du rendu sans imposer de seuil physique', () => {
    const valid = createAppState().rendering;
    const invalidRenderingStates: readonly unknown[] = [
      { ...valid, displayMode: 'volume' },
      { ...valid, observable: 'amplitude' },
      { ...valid, theme: 'sepia' },
      { ...valid, isoDensityFraction: 0 },
      { ...valid, isoDensityFraction: 1 },
      { ...valid, isoDensityFraction: Number.NaN },
      { ...valid, pointOpacity: -0.01 },
      { ...valid, pointOpacity: 1.01 },
      { ...valid, pointSizePixels: 0 },
      { ...valid, pointSizePixels: Number.POSITIVE_INFINITY },
      { ...valid, showAxes: 'oui' },
      { ...valid, cameraRotationEnabled: 1 },
    ];

    for (const invalidState of invalidRenderingStates) {
      expect(() => normalizeRenderingState(invalidState)).toThrow();
    }
  });

  it('valide la configuration reproductible selon les contrats du sampler', () => {
    expect(normalizeSamplingConfiguration({ sampleCount: 32, seed: 0 })).toEqual({
      sampleCount: 32,
      seed: 0,
    });

    const invalidConfigurations: readonly unknown[] = [
      null,
      { sampleCount: 0, seed: 1 },
      { sampleCount: 1.5, seed: 1 },
      { sampleCount: MAX_ORBITAL_SAMPLE_COUNT + 1, seed: 1 },
      { sampleCount: 1, seed: -1 },
      { sampleCount: 1, seed: 0x1_0000_0000 },
      { sampleCount: 1, seed: 1.5 },
    ];
    for (const invalidConfiguration of invalidConfigurations) {
      expect(() => normalizeSamplingConfiguration(invalidConfiguration)).toThrow();
    }
  });

  it('conserve explicitement le module 2D différé et valide ses drapeaux', () => {
    const state = createAppState();
    expect(() =>
      normalizeAppState({
        ...state,
        legacy: { legacy2DInitialized: false, showLegacy2D: 'oui' },
      }),
    ).toThrow(TypeError);
  });
});
