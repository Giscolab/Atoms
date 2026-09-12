import { describe, expect, it } from 'vitest';

import { ORBIT_CAMERA_ELEVATION_LIMIT_RADIANS } from '../../src/rendering/renderingContracts';
import { createAppState } from '../../src/state/appState';
import {
  createScientificSnapshot,
  normalizeScientificSnapshot,
  parseScientificSnapshotJson,
  SCIENTIFIC_SNAPSHOT_FORMAT,
  SCIENTIFIC_SNAPSHOT_SCHEMA_VERSION,
  serializeScientificSnapshot,
} from '../../src/state/scientificSnapshot';

const CAMERA = {
  azimuthRadians: 0.65,
  distanceBohr: 24.5,
  elevationRadians: 0.42,
} as const;

function snapshotObject(): ReturnType<typeof createScientificSnapshot> {
  return createScientificSnapshot(createAppState(), CAMERA);
}

function mutableSnapshot(): Record<string, unknown> {
  return JSON.parse(JSON.stringify(snapshotObject())) as Record<string, unknown>;
}

describe('snapshot scientifique Atoms', () => {
  it('sérialise un état reproductible, versionné et auto-descriptif', () => {
    const snapshot = snapshotObject();

    expect(snapshot.format).toBe(SCIENTIFIC_SNAPSHOT_FORMAT);
    expect(snapshot.schemaVersion).toBe(SCIENTIFIC_SNAPSHOT_SCHEMA_VERSION);
    expect(snapshot.units).toEqual({
      angle: 'radian',
      isoDensity: 'fraction-of-grid-maximum',
      length: 'bohr',
      pointSize: 'css-pixel',
    });
    expect(snapshot.provenance).toMatchObject({
      fieldVersion: 'hydrogen-orbital-field-v1',
      prngVersion: 'xorshift32-marsaglia-seed-v1',
      samplerVersion: 'hydrogen-orbital-sampler-v1',
      scienceEngineVersion: 'hydrogen-wavefunctions-phase3-v1',
      workerProtocolVersion: 'orbital-worker-protocol-v2',
    });
    expect(snapshot.view.camera).toEqual(CAMERA);
  });

  it('produit un JSON déterministe sans timestamp parasite et round-trip exact', () => {
    const snapshot = snapshotObject();
    const first = serializeScientificSnapshot(snapshot);
    const second = serializeScientificSnapshot(snapshot);

    expect(first).toBe(second);
    expect(first).not.toContain('timestamp');
    expect(first.endsWith('\n')).toBe(true);
    expect(parseScientificSnapshotJson(first)).toEqual(snapshot);
  });

  it.each([
    ['format', 'other-format'],
    ['schemaVersion', 2],
  ])('refuse un %s incompatible', (key, value) => {
    const candidate = mutableSnapshot();
    candidate[key] = value;
    expect(() => normalizeScientificSnapshot(candidate)).toThrow(/incompatible/i);
  });

  it('refuse les unités ambiguës ou converties implicitement', () => {
    const candidate = mutableSnapshot();
    const units = candidate.units as Record<string, unknown>;
    units.length = 'meter';

    expect(() => normalizeScientificSnapshot(candidate)).toThrow(/unité de longueur/i);
  });

  it.each([
    ['samplerVersion', 'hydrogen-orbital-sampler-v0'],
    ['prngVersion', 'random-v2'],
    ['fieldVersion', 'hydrogen-orbital-field-v2'],
    ['scienceEngineVersion', 'future-engine'],
    ['workerProtocolVersion', 'orbital-worker-protocol-v99'],
  ])('refuse une provenance %s différente', (key, value) => {
    const candidate = mutableSnapshot();
    const provenance = candidate.provenance as Record<string, unknown>;
    provenance[key] = value;

    expect(() => normalizeScientificSnapshot(candidate)).toThrow(/version/i);
  });

  it('refuse une dérive des paramètres numériques internes du Worker', () => {
    const candidate = mutableSnapshot();
    const provenance = candidate.provenance as Record<string, unknown>;
    const workerOptions = provenance.workerOptions as Record<string, unknown>;
    workerOptions.fieldResolution = 48;

    expect(() => normalizeScientificSnapshot(candidate)).toThrow(/fieldResolution/i);
  });

  it('refuse les champs inconnus au lieu de les ignorer silencieusement', () => {
    const candidate = mutableSnapshot();
    candidate.comment = 'champ arbitraire';

    expect(() => normalizeScientificSnapshot(candidate)).toThrow(/exactement les champs/i);
  });

  it('refuse un état orbital qui mélange les bases réelle et complexe', () => {
    const candidate = mutableSnapshot();
    const state = candidate.state as Record<string, unknown>;
    state.orbital = { basis: 'real', l: 2, m: 0, n: 3, orbital: 'd_z2' };

    expect(() => normalizeScientificSnapshot(candidate)).toThrow(/exactement les champs/i);
  });

  it('refuse un état hors du domaine n = 1…9 exposé par l’interface', () => {
    const candidate = mutableSnapshot();
    const state = candidate.state as Record<string, unknown>;
    state.orbital = { basis: 'complex', l: 2, m: 0, n: 10 };

    expect(() => normalizeScientificSnapshot(candidate)).toThrow(/n dans \[1, 9\]/i);
  });

  it('refuse des paramètres valides pour le moteur mais non représentables par l’interface', () => {
    const tooFewSamples = mutableSnapshot();
    const sampling = (tooFewSamples.state as Record<string, unknown>).sampling as Record<
      string,
      unknown
    >;
    sampling.sampleCount = 1_000;
    expect(() => normalizeScientificSnapshot(tooFewSamples)).toThrow(/nombre d'échantillons/i);

    const offStepOpacity = mutableSnapshot();
    const rendering = (offStepOpacity.state as Record<string, unknown>).rendering as Record<
      string,
      unknown
    >;
    rendering.pointOpacity = 0.83;
    expect(() => normalizeScientificSnapshot(offStepOpacity)).toThrow(/pas de 0.05/i);
  });

  it('refuse une caméra non reproductible', () => {
    const candidate = mutableSnapshot();
    const view = candidate.view as Record<string, unknown>;
    view.camera = {
      ...CAMERA,
      elevationRadians: ORBIT_CAMERA_ELEVATION_LIMIT_RADIANS + 0.001,
    };

    expect(() => normalizeScientificSnapshot(candidate)).toThrow(/élévation caméra/i);
  });

  it('signale distinctement un JSON malformé', () => {
    expect(() => parseScientificSnapshotJson('{"format":')).toThrow(/JSON scientifique invalide/i);
  });
});
