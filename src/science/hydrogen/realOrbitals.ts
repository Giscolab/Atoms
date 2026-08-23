import {
  addComplex,
  createComplex,
  multiplyComplex,
  scaleComplex,
  subtractComplex,
  type Complex,
} from '../math/complex';
import { hydrogenRadialWavefunction } from './radialWavefunction';
import { sphericalHarmonic } from './sphericalHarmonics';

export type RealOrbitalName = 'p_x' | 'p_y' | 'p_z' | 'd_xy' | 'd_xz' | 'd_yz' | 'd_x2_y2' | 'd_z2';

export type RealOrbitalCombinationKind = 'zonal' | 'cosine' | 'sine';

export interface RealOrbitalDefinition {
  readonly absoluteOrder: number;
  readonly kind: RealOrbitalCombinationKind;
  readonly l: number;
  readonly minimumPrincipalQuantumNumber: number;
}

/**
 * Métadonnées des orbitales réelles standard de Phase 3.
 *
 * Les familles p et d sont des combinaisons des états complexes ±m, jamais des
 * renommages d'un état m individuel.
 */
export const REAL_ORBITAL_DEFINITIONS = {
  p_x: { absoluteOrder: 1, kind: 'cosine', l: 1, minimumPrincipalQuantumNumber: 2 },
  p_y: { absoluteOrder: 1, kind: 'sine', l: 1, minimumPrincipalQuantumNumber: 2 },
  p_z: { absoluteOrder: 0, kind: 'zonal', l: 1, minimumPrincipalQuantumNumber: 2 },
  d_xy: { absoluteOrder: 2, kind: 'sine', l: 2, minimumPrincipalQuantumNumber: 3 },
  d_xz: { absoluteOrder: 1, kind: 'cosine', l: 2, minimumPrincipalQuantumNumber: 3 },
  d_yz: { absoluteOrder: 1, kind: 'sine', l: 2, minimumPrincipalQuantumNumber: 3 },
  d_x2_y2: { absoluteOrder: 2, kind: 'cosine', l: 2, minimumPrincipalQuantumNumber: 3 },
  d_z2: { absoluteOrder: 0, kind: 'zonal', l: 2, minimumPrincipalQuantumNumber: 3 },
} as const satisfies Readonly<Record<RealOrbitalName, RealOrbitalDefinition>>;

const IMAGINARY_UNIT = createComplex(0, 1);

function isRealOrbitalName(value: unknown): value is RealOrbitalName {
  return typeof value === 'string' && Object.hasOwn(REAL_ORBITAL_DEFINITIONS, value);
}

function requireRealOrbitalDefinition(name: unknown): RealOrbitalDefinition {
  if (!isRealOrbitalName(name)) {
    throw new RangeError(`Nom d'orbitale réelle inconnu : ${String(name)}.`);
  }
  return REAL_ORBITAL_DEFINITIONS[name];
}

function evaluateRealSphericalHarmonic(
  definition: RealOrbitalDefinition,
  thetaRadians: number,
  phiRadians: number,
): Complex {
  if (definition.kind === 'zonal') {
    return sphericalHarmonic(definition.l, 0, thetaRadians, phiRadians);
  }

  const order = definition.absoluteOrder;
  const negativeOrder = sphericalHarmonic(definition.l, -order, thetaRadians, phiRadians);
  const positiveOrder = sphericalHarmonic(definition.l, order, thetaRadians, phiRadians);
  const condonShortleySign = order % 2 === 0 ? 1 : -1;
  const signedPositiveOrder = scaleComplex(positiveOrder, condonShortleySign);

  if (definition.kind === 'cosine') {
    // C_lm = [Y_l^(-m) + (-1)^m Y_l^m] / sqrt(2).
    return scaleComplex(addComplex(negativeOrder, signedPositiveOrder), Math.SQRT1_2);
  }

  // S_lm = i [Y_l^(-m) - (-1)^m Y_l^m] / sqrt(2).
  return scaleComplex(
    multiplyComplex(IMAGINARY_UNIT, subtractComplex(negativeOrder, signedPositiveOrder)),
    Math.SQRT1_2,
  );
}

/**
 * Partie angulaire réelle normalisée, avec les signes cartésiens usuels :
 * p_x ∝ x/r, p_y ∝ y/r et p_z ∝ z/r.
 */
export function realSphericalHarmonic(
  name: RealOrbitalName,
  thetaRadians: number,
  phiRadians: number,
): Complex {
  return evaluateRealSphericalHarmonic(
    requireRealOrbitalDefinition(name),
    thetaRadians,
    phiRadians,
  );
}

/** Fonction d'onde d'une orbitale réelle, normalisée dans le contrat radial en a0. */
export function realOrbitalWavefunction(
  n: number,
  name: RealOrbitalName,
  rBohr: number,
  thetaRadians: number,
  phiRadians: number,
): Complex {
  const definition = requireRealOrbitalDefinition(name);
  if (!Number.isSafeInteger(n) || n < definition.minimumPrincipalQuantumNumber) {
    throw new RangeError(
      `L'orbitale ${name} requiert un entier sûr n >= ${definition.minimumPrincipalQuantumNumber} : ${n}.`,
    );
  }

  const radialAmplitude = hydrogenRadialWavefunction(n, definition.l, rBohr);
  const angularAmplitude = evaluateRealSphericalHarmonic(definition, thetaRadians, phiRadians);
  return scaleComplex(angularAmplitude, radialAmplitude);
}
