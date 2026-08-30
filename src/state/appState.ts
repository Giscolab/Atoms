import type {
  OrbitalDisplayMode,
  OrbitalObservable,
  RenderTheme,
} from '../rendering/renderingContracts';
import { REAL_ORBITAL_DEFINITIONS, type RealOrbitalName } from '../science/hydrogen/realOrbitals';
import { assertValidQuantumNumbers } from '../science/quantum/quantumNumbers';
import {
  type ComplexOrbitalState,
  type OrbitalSamplingState,
  type RealOrbitalState,
} from '../sampling/contracts';
import { MAX_ORBITAL_SAMPLE_COUNT } from '../sampling/orbitalSampler';
import { createSeededRandom } from '../sampling/rng';

export type { ComplexOrbitalState, OrbitalSamplingState, RealOrbitalState };

/** Seed reproductible par défaut ; 0x41544f4d encode ASCII « ATOM ». */
export const DEFAULT_ORBITAL_SEED = 0x4154_4f4d;
export const DEFAULT_ORBITAL_SAMPLE_COUNT = 15_000;

export interface OrbitalSamplingConfiguration {
  readonly sampleCount: number;
  readonly seed: number;
}

export interface OrbitalRenderingState {
  readonly cameraRotationEnabled: boolean;
  readonly displayMode: OrbitalDisplayMode;
  readonly isoDensityFraction: number;
  readonly observable: OrbitalObservable;
  readonly pointOpacity: number;
  readonly pointSizePixels: number;
  readonly showAxes: boolean;
  readonly showNodes: boolean;
  readonly theme: RenderTheme;
}

export interface Legacy2DState {
  readonly legacy2DInitialized: boolean;
  readonly showLegacy2D: boolean;
}

/**
 * État applicatif sérialisable de Phase 7. Les grands buffers transférables du
 * sampler et du champ orbital restent la propriété du Worker et du renderer.
 */
export interface AppState {
  readonly legacy: Legacy2DState;
  readonly orbital: OrbitalSamplingState;
  readonly rendering: OrbitalRenderingState;
  readonly sampling: OrbitalSamplingConfiguration;
}

function requireRecord(value: unknown, quantity: string): Record<string, unknown> {
  if (typeof value !== 'object' || value === null) {
    throw new TypeError(`${quantity} doit être un objet.`);
  }
  return value as Record<string, unknown>;
}

function isRealOrbitalName(value: unknown): value is RealOrbitalName {
  return typeof value === 'string' && Object.hasOwn(REAL_ORBITAL_DEFINITIONS, value);
}

/** Valide puis copie l'union scientifique sans la convertir dans une autre base. */
export function normalizeOrbitalState(value: unknown): OrbitalSamplingState {
  const state = requireRecord(value, "L'état orbital");
  if (state.basis === 'complex') {
    assertValidQuantumNumbers(state);
    return { basis: 'complex', n: state.n, l: state.l, m: state.m };
  }

  if (state.basis === 'real') {
    if (!isRealOrbitalName(state.orbital)) {
      throw new RangeError(`Nom d'orbitale réelle inconnu : ${String(state.orbital)}.`);
    }
    const definition = REAL_ORBITAL_DEFINITIONS[state.orbital];
    const principalQuantumNumber = state.n;
    if (
      typeof principalQuantumNumber !== 'number' ||
      !Number.isSafeInteger(principalQuantumNumber) ||
      principalQuantumNumber < definition.minimumPrincipalQuantumNumber
    ) {
      throw new RangeError(
        `L'orbitale ${state.orbital} requiert un entier sûr n >= ${definition.minimumPrincipalQuantumNumber}.`,
      );
    }
    return { basis: 'real', n: principalQuantumNumber, orbital: state.orbital };
  }

  throw new RangeError(`Base orbitale inconnue : ${String(state.basis)}.`);
}

function requireBoolean(value: unknown, quantity: string): boolean {
  if (typeof value !== 'boolean') throw new TypeError(`${quantity} doit être un booléen.`);
  return value;
}

function requireFiniteNumber(value: unknown, quantity: string): number {
  if (typeof value !== 'number' || !Number.isFinite(value)) {
    throw new RangeError(`${quantity} doit être un nombre fini.`);
  }
  return value;
}

function requireDisplayMode(value: unknown): OrbitalDisplayMode {
  if (value !== 'cloud' && value !== 'hybrid' && value !== 'isosurface') {
    throw new RangeError(`Mode d'affichage inconnu : ${String(value)}.`);
  }
  return value;
}

function requireObservable(value: unknown): OrbitalObservable {
  if (value !== 'density' && value !== 'phase') {
    throw new RangeError(`Observable de rendu inconnue : ${String(value)}.`);
  }
  return value;
}

function requireTheme(value: unknown): RenderTheme {
  if (value !== 'dark' && value !== 'light') {
    throw new RangeError(`Thème inconnu : ${String(value)}.`);
  }
  return value;
}

/** Valide puis copie les paramètres visuels, sans altérer l'état reçu. */
export function normalizeRenderingState(value: unknown): OrbitalRenderingState {
  const rendering = requireRecord(value, "L'état de rendu");
  const isoDensityFraction = requireFiniteNumber(
    rendering.isoDensityFraction,
    "Le seuil d'isodensité",
  );
  const pointOpacity = requireFiniteNumber(rendering.pointOpacity, "L'opacité des points");
  const pointSizePixels = requireFiniteNumber(
    rendering.pointSizePixels,
    'La taille des points en pixels',
  );

  if (isoDensityFraction <= 0 || isoDensityFraction >= 1) {
    throw new RangeError("Le seuil d'isodensité doit appartenir à ]0,1[.");
  }
  if (pointOpacity < 0 || pointOpacity > 1) {
    throw new RangeError("L'opacité des points doit appartenir à [0,1].");
  }
  if (pointSizePixels <= 0) {
    throw new RangeError('La taille des points en pixels doit être strictement positive.');
  }

  return {
    cameraRotationEnabled: requireBoolean(
      rendering.cameraRotationEnabled,
      'Le réglage de rotation caméra',
    ),
    displayMode: requireDisplayMode(rendering.displayMode),
    isoDensityFraction,
    observable: requireObservable(rendering.observable),
    pointOpacity,
    pointSizePixels,
    showAxes: requireBoolean(rendering.showAxes, "L'affichage des axes"),
    showNodes: requireBoolean(rendering.showNodes, "L'affichage des nœuds"),
    theme: requireTheme(rendering.theme),
  };
}

/** Valide puis copie la configuration reproductible du sampler. */
export function normalizeSamplingConfiguration(value: unknown): OrbitalSamplingConfiguration {
  const sampling = requireRecord(value, 'La configuration du sampler');
  if (
    typeof sampling.sampleCount !== 'number' ||
    !Number.isSafeInteger(sampling.sampleCount) ||
    sampling.sampleCount <= 0 ||
    sampling.sampleCount > MAX_ORBITAL_SAMPLE_COUNT
  ) {
    throw new RangeError(
      `Le nombre d'échantillons doit être un entier sûr dans [1, ${MAX_ORBITAL_SAMPLE_COUNT}].`,
    );
  }
  if (typeof sampling.seed !== 'number') {
    throw new RangeError('La seed doit être un entier non signé sur 32 bits.');
  }

  return {
    sampleCount: sampling.sampleCount,
    seed: createSeededRandom(sampling.seed).normalizedSeed,
  };
}

function normalizeLegacy2DState(value: unknown): Legacy2DState {
  const legacy = requireRecord(value, "L'état du module 2D legacy");
  return {
    legacy2DInitialized: requireBoolean(
      legacy.legacy2DInitialized,
      "L'initialisation du module 2D legacy",
    ),
    showLegacy2D: requireBoolean(legacy.showLegacy2D, "L'affichage du module 2D legacy"),
  };
}

/** Normalisation pure : validation et copie profonde des quatre sous-états. */
export function normalizeAppState(value: unknown): AppState {
  const state = requireRecord(value, "L'état applicatif");
  return {
    legacy: normalizeLegacy2DState(state.legacy),
    orbital: normalizeOrbitalState(state.orbital),
    rendering: normalizeRenderingState(state.rendering),
    sampling: normalizeSamplingConfiguration(state.sampling),
  };
}

/** État initial correspondant à la maquette UI validée. */
export function createAppState(): AppState {
  return normalizeAppState({
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
}
