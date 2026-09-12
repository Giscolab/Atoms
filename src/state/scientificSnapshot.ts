import { ORBITAL_SAMPLER_VERSION } from '../sampling/contracts';
import { ORBITAL_FIELD_VERSION } from '../sampling/fieldContracts';
import { SAMPLER_PRNG_VERSION } from '../sampling/rng';
import { normalizeAppState, type AppState } from './appState';
import {
  DEFAULT_ORBITAL_WORKER_OPTIONS,
  ORBITAL_WORKER_PROTOCOL_VERSION,
  SCIENCE_ENGINE_VERSION,
  type OrbitalWorkerOptions,
} from '../workers/orbitalSamplingProtocol';
import {
  ORBIT_CAMERA_ELEVATION_LIMIT_RADIANS,
  type OrbitCameraState,
} from '../rendering/renderingContracts';

export const SCIENTIFIC_SNAPSHOT_FORMAT = 'atoms-scientific-snapshot' as const;
export const SCIENTIFIC_SNAPSHOT_SCHEMA_VERSION = 1 as const;
export const SCIENTIFIC_SNAPSHOT_MAX_PRINCIPAL_QUANTUM_NUMBER = 9 as const;

export interface ScientificSnapshotUnits {
  readonly angle: 'radian';
  readonly isoDensity: 'fraction-of-grid-maximum';
  readonly length: 'bohr';
  readonly pointSize: 'css-pixel';
}

export interface ScientificSnapshotProvenance {
  readonly fieldVersion: typeof ORBITAL_FIELD_VERSION;
  readonly prngVersion: typeof SAMPLER_PRNG_VERSION;
  readonly samplerVersion: typeof ORBITAL_SAMPLER_VERSION;
  readonly scienceEngineVersion: typeof SCIENCE_ENGINE_VERSION;
  readonly workerOptions: OrbitalWorkerOptions;
  readonly workerProtocolVersion: typeof ORBITAL_WORKER_PROTOCOL_VERSION;
}

export interface ScientificSnapshotView {
  readonly camera: OrbitCameraState;
}

export interface ScientificSnapshot {
  readonly format: typeof SCIENTIFIC_SNAPSHOT_FORMAT;
  readonly provenance: ScientificSnapshotProvenance;
  readonly schemaVersion: typeof SCIENTIFIC_SNAPSHOT_SCHEMA_VERSION;
  readonly state: AppState;
  readonly units: ScientificSnapshotUnits;
  readonly view: ScientificSnapshotView;
}

const SNAPSHOT_UNITS: ScientificSnapshotUnits = {
  angle: 'radian',
  isoDensity: 'fraction-of-grid-maximum',
  length: 'bohr',
  pointSize: 'css-pixel',
};

const SNAPSHOT_UI_LIMITS = {
  isoDensityFraction: { max: 0.8, min: 0.05, step: 0.01 },
  pointOpacity: { max: 1, min: 0.2, step: 0.05 },
  pointSizePixels: { max: 4, min: 0.6, step: 0.1 },
  sampleCount: { max: 60_000, min: 2_000, step: 1_000 },
} as const;

function requireRecord(value: unknown, label: string): Record<string, unknown> {
  if (typeof value !== 'object' || value === null || Array.isArray(value)) {
    throw new TypeError(`${label} doit être un objet JSON.`);
  }
  return value as Record<string, unknown>;
}

function requireExactKeys(
  record: Record<string, unknown>,
  expectedKeys: readonly string[],
  label: string,
): void {
  const actual = Object.keys(record).sort();
  const expected = [...expectedKeys].sort();
  if (actual.length !== expected.length || actual.some((key, index) => key !== expected[index])) {
    throw new RangeError(
      `${label} doit contenir exactement les champs : ${expectedKeys.join(', ')}.`,
    );
  }
}

function requireLiteral<T extends string | number>(value: unknown, expected: T, label: string): T {
  if (value !== expected) {
    throw new RangeError(
      `${label} incompatible : ${String(value)} (attendu : ${String(expected)}).`,
    );
  }
  return expected;
}

function requireUiRepresentableNumber(
  value: unknown,
  limits: { readonly max: number; readonly min: number; readonly step: number },
  label: string,
): number {
  if (
    typeof value !== 'number' ||
    !Number.isFinite(value) ||
    value < limits.min ||
    value > limits.max
  ) {
    throw new RangeError(`${label} doit appartenir à [${limits.min}, ${limits.max}].`);
  }
  const stepIndex = (value - limits.min) / limits.step;
  if (Math.abs(stepIndex - Math.round(stepIndex)) > 1e-9) {
    throw new RangeError(`${label} doit respecter un pas de ${limits.step}.`);
  }
  return value;
}

function normalizeCameraState(value: unknown): OrbitCameraState {
  const camera = requireRecord(value, 'La caméra');
  requireExactKeys(camera, ['azimuthRadians', 'distanceBohr', 'elevationRadians'], 'La caméra');
  const azimuthRadians = camera.azimuthRadians;
  const distanceBohr = camera.distanceBohr;
  const elevationRadians = camera.elevationRadians;
  if (
    typeof azimuthRadians !== 'number' ||
    !Number.isFinite(azimuthRadians) ||
    typeof elevationRadians !== 'number' ||
    !Number.isFinite(elevationRadians)
  ) {
    throw new RangeError('Les angles de caméra doivent être des nombres finis en radians.');
  }
  if (Math.abs(elevationRadians) > ORBIT_CAMERA_ELEVATION_LIMIT_RADIANS) {
    throw new RangeError('L’élévation caméra dépasse la limite de navigation du renderer.');
  }
  if (typeof distanceBohr !== 'number' || !Number.isFinite(distanceBohr) || distanceBohr <= 0) {
    throw new RangeError('La distance caméra doit être finie et strictement positive en a₀.');
  }
  return { azimuthRadians, distanceBohr, elevationRadians };
}

function assertStrictAppStateShape(value: unknown): void {
  const state = requireRecord(value, "L'état applicatif");
  requireExactKeys(state, ['orbital', 'rendering', 'sampling'], "L'état applicatif");

  const orbital = requireRecord(state.orbital, "L'état orbital");
  if (orbital.basis === 'complex') {
    requireExactKeys(orbital, ['basis', 'l', 'm', 'n'], "L'état orbital complexe");
  } else if (orbital.basis === 'real') {
    requireExactKeys(orbital, ['basis', 'n', 'orbital'], "L'orbitale réelle");
  } else {
    throw new RangeError(`Base orbitale inconnue : ${String(orbital.basis)}.`);
  }
  if (
    typeof orbital.n !== 'number' ||
    !Number.isSafeInteger(orbital.n) ||
    orbital.n < 1 ||
    orbital.n > SCIENTIFIC_SNAPSHOT_MAX_PRINCIPAL_QUANTUM_NUMBER
  ) {
    throw new RangeError(
      `Le snapshot Atoms prend en charge n dans [1, ${SCIENTIFIC_SNAPSHOT_MAX_PRINCIPAL_QUANTUM_NUMBER}].`,
    );
  }

  const rendering = requireRecord(state.rendering, "L'état de rendu");
  requireExactKeys(
    rendering,
    [
      'cameraRotationEnabled',
      'displayMode',
      'isoDensityFraction',
      'observable',
      'pointOpacity',
      'pointSizePixels',
      'showAxes',
      'showNodes',
      'theme',
    ],
    "L'état de rendu",
  );
  requireUiRepresentableNumber(
    rendering.isoDensityFraction,
    SNAPSHOT_UI_LIMITS.isoDensityFraction,
    "Le seuil d'isodensité du snapshot",
  );
  requireUiRepresentableNumber(
    rendering.pointOpacity,
    SNAPSHOT_UI_LIMITS.pointOpacity,
    "L'opacité des points du snapshot",
  );
  requireUiRepresentableNumber(
    rendering.pointSizePixels,
    SNAPSHOT_UI_LIMITS.pointSizePixels,
    'La taille des points du snapshot',
  );

  const sampling = requireRecord(state.sampling, "La configuration d'échantillonnage");
  requireExactKeys(sampling, ['sampleCount', 'seed'], "La configuration d'échantillonnage");
  requireUiRepresentableNumber(
    sampling.sampleCount,
    SNAPSHOT_UI_LIMITS.sampleCount,
    "Le nombre d'échantillons du snapshot",
  );
}

function normalizeWorkerOptions(value: unknown): OrbitalWorkerOptions {
  const options = requireRecord(value, 'Les paramètres numériques du Worker');
  requireExactKeys(
    options,
    [
      'angularChartPointCount',
      'fieldResolution',
      'planeExplorationPointCount',
      'radialChartPointCount',
      'radialCoverageProbability',
    ],
    'Les paramètres numériques du Worker',
  );

  for (const [key, expected] of Object.entries(DEFAULT_ORBITAL_WORKER_OPTIONS)) {
    if (options[key] !== expected) {
      throw new RangeError(
        `Paramètre numérique incompatible ${key}=${String(options[key])} (attendu : ${String(expected)}).`,
      );
    }
  }
  return { ...DEFAULT_ORBITAL_WORKER_OPTIONS };
}

function normalizeUnits(value: unknown): ScientificSnapshotUnits {
  const units = requireRecord(value, 'Les unités');
  requireExactKeys(units, ['angle', 'isoDensity', 'length', 'pointSize'], 'Les unités');
  return {
    angle: requireLiteral(units.angle, SNAPSHOT_UNITS.angle, "L'unité angulaire"),
    isoDensity: requireLiteral(
      units.isoDensity,
      SNAPSHOT_UNITS.isoDensity,
      "L'unité du seuil d'isodensité",
    ),
    length: requireLiteral(units.length, SNAPSHOT_UNITS.length, "L'unité de longueur"),
    pointSize: requireLiteral(
      units.pointSize,
      SNAPSHOT_UNITS.pointSize,
      "L'unité de taille des points",
    ),
  };
}

function normalizeProvenance(value: unknown): ScientificSnapshotProvenance {
  const provenance = requireRecord(value, 'La provenance scientifique');
  requireExactKeys(
    provenance,
    [
      'fieldVersion',
      'prngVersion',
      'samplerVersion',
      'scienceEngineVersion',
      'workerOptions',
      'workerProtocolVersion',
    ],
    'La provenance scientifique',
  );
  return {
    fieldVersion: requireLiteral(
      provenance.fieldVersion,
      ORBITAL_FIELD_VERSION,
      'La version du champ orbital',
    ),
    prngVersion: requireLiteral(provenance.prngVersion, SAMPLER_PRNG_VERSION, 'La version du PRNG'),
    samplerVersion: requireLiteral(
      provenance.samplerVersion,
      ORBITAL_SAMPLER_VERSION,
      'La version du sampler',
    ),
    scienceEngineVersion: requireLiteral(
      provenance.scienceEngineVersion,
      SCIENCE_ENGINE_VERSION,
      'La version du noyau scientifique',
    ),
    workerOptions: normalizeWorkerOptions(provenance.workerOptions),
    workerProtocolVersion: requireLiteral(
      provenance.workerProtocolVersion,
      ORBITAL_WORKER_PROTOCOL_VERSION,
      'La version du protocole Worker',
    ),
  };
}

export function createScientificSnapshot(
  state: AppState,
  camera: OrbitCameraState,
): ScientificSnapshot {
  return {
    format: SCIENTIFIC_SNAPSHOT_FORMAT,
    provenance: {
      fieldVersion: ORBITAL_FIELD_VERSION,
      prngVersion: SAMPLER_PRNG_VERSION,
      samplerVersion: ORBITAL_SAMPLER_VERSION,
      scienceEngineVersion: SCIENCE_ENGINE_VERSION,
      workerOptions: { ...DEFAULT_ORBITAL_WORKER_OPTIONS },
      workerProtocolVersion: ORBITAL_WORKER_PROTOCOL_VERSION,
    },
    schemaVersion: SCIENTIFIC_SNAPSHOT_SCHEMA_VERSION,
    state: normalizeAppState(state),
    units: { ...SNAPSHOT_UNITS },
    view: { camera: normalizeCameraState(camera) },
  };
}

export function normalizeScientificSnapshot(value: unknown): ScientificSnapshot {
  const snapshot = requireRecord(value, 'Le snapshot scientifique');
  requireExactKeys(
    snapshot,
    ['format', 'provenance', 'schemaVersion', 'state', 'units', 'view'],
    'Le snapshot scientifique',
  );
  requireLiteral(snapshot.format, SCIENTIFIC_SNAPSHOT_FORMAT, 'Le format');
  requireLiteral(
    snapshot.schemaVersion,
    SCIENTIFIC_SNAPSHOT_SCHEMA_VERSION,
    'La version du schéma',
  );
  assertStrictAppStateShape(snapshot.state);
  const view = requireRecord(snapshot.view, 'La vue');
  requireExactKeys(view, ['camera'], 'La vue');

  return {
    format: SCIENTIFIC_SNAPSHOT_FORMAT,
    provenance: normalizeProvenance(snapshot.provenance),
    schemaVersion: SCIENTIFIC_SNAPSHOT_SCHEMA_VERSION,
    state: normalizeAppState(snapshot.state),
    units: normalizeUnits(snapshot.units),
    view: { camera: normalizeCameraState(view.camera) },
  };
}

export function parseScientificSnapshotJson(text: string): ScientificSnapshot {
  let parsed: unknown;
  try {
    parsed = JSON.parse(text) as unknown;
  } catch (error) {
    const detail = error instanceof Error ? error.message : String(error);
    throw new SyntaxError(`JSON scientifique invalide : ${detail}`, { cause: error });
  }
  return normalizeScientificSnapshot(parsed);
}

export function serializeScientificSnapshot(snapshot: ScientificSnapshot): string {
  return `${JSON.stringify(normalizeScientificSnapshot(snapshot), null, 2)}
`;
}
