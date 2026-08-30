import { hydrogenRadialWavefunction } from '../science/hydrogen/radialWavefunction';
import {
  REAL_ORBITAL_DEFINITIONS,
  realSphericalHarmonic,
  type RealOrbitalName,
} from '../science/hydrogen/realOrbitals';
import { sphericalHarmonic } from '../science/hydrogen/sphericalHarmonics';
import { complexModulusSquared, type Complex } from '../science/math/complex';
import { assertValidQuantumNumbers } from '../science/quantum/quantumNumbers';
import type { OrbitalSamplingState } from './contracts';

const TWO_PI = 2 * Math.PI;
const PLANE_SCORE_ROUNDING_MARGIN = 2048 * Number.EPSILON;

export type PrincipalPlane = 'xy' | 'xz' | 'yz';

export interface RadialDistributionChartData {
  /** Maximum de P_nl(r), exprimé en a₀⁻¹, avant normalisation graphique. */
  readonly maximumRawDensityInverseBohr: number;
  /** P_nl(r) / max(P_nl), grandeur sans dimension dans [0,1]. */
  readonly normalizedDensity: Float64Array;
  /** Valeurs de P_nl(r) = r² |R_nl(r)|², exprimées en a₀⁻¹. */
  readonly rawDensityInverseBohr: Float64Array;
  /** Abscisses r/a₀, de zéro à l'étendue demandée, bornes incluses. */
  readonly radiiBohr: Float64Array;
}

export interface AngularDensityCutChartData {
  /** Angles du plan, de 0 à 2π, bornes incluses. */
  readonly anglesRadians: Float64Array;
  /** Maximum brut de |Y(Ω)|² sur la coupe retenue, en sr⁻¹. */
  readonly maximumRawDensityInverseSteradians: number;
  /** Rayon polaire graphique |Y(Ω)|² / max_coupe(|Y(Ω)|²). */
  readonly normalizedRadius: Float64Array;
  readonly plane: PrincipalPlane;
  readonly planeLabel: string;
  /** Moyenne exploratoire de chaque coupe ; ce score n'est pas une probabilité. */
  readonly planeMeanDensityScoresInverseSteradians: Readonly<Record<PrincipalPlane, number>>;
  /** Valeurs brutes de |Y(Ω)|² sur le grand cercle retenu, en sr⁻¹. */
  readonly rawDensityInverseSteradians: Float64Array;
}

export interface OrbitalChartDataOptions {
  readonly angularPointCount: number;
  readonly radialExtentBohr: number;
  readonly radialPointCount: number;
  readonly planeExplorationPointCount: number;
}

export interface OrbitalChartData {
  readonly angularCut: AngularDensityCutChartData;
  readonly radialDistribution: RadialDistributionChartData;
}

interface ValidatedAngularState {
  readonly degree: number;
  readonly evaluate: (thetaRadians: number, phiRadians: number) => Complex;
}

interface UnitDirection {
  readonly x: number;
  readonly y: number;
  readonly z: number;
}

const PRINCIPAL_PLANES = ['xy', 'xz', 'yz'] as const satisfies readonly PrincipalPlane[];

const PLANE_LABELS: Readonly<Record<PrincipalPlane, string>> = {
  xy: 'Plan xy',
  xz: 'Plan xz',
  yz: 'Plan yz',
};

function requirePointCount(pointCount: number, quantity: string, minimum: number): void {
  if (!Number.isSafeInteger(pointCount) || pointCount < minimum) {
    throw new RangeError(`${quantity} doit être un entier sûr supérieur ou égal à ${minimum}.`);
  }
}

function requirePositiveFinite(value: number, quantity: string): void {
  if (!Number.isFinite(value) || value <= 0) {
    throw new RangeError(`${quantity} doit être strictement positive et finie : ${value}.`);
  }
}

function typedArrayValue(values: Float64Array, index: number, quantity: string): number {
  const value = values[index];
  if (value === undefined) {
    throw new RangeError(`${quantity} est absente à l'index ${index}.`);
  }
  return value;
}

function isRealOrbitalName(value: unknown): value is RealOrbitalName {
  return typeof value === 'string' && Object.hasOwn(REAL_ORBITAL_DEFINITIONS, value);
}

function validateAngularState(state: OrbitalSamplingState): ValidatedAngularState {
  if (state.basis === 'complex') {
    assertValidQuantumNumbers(state);
    return {
      degree: state.l,
      evaluate: (thetaRadians, phiRadians) =>
        sphericalHarmonic(state.l, state.m, thetaRadians, phiRadians),
    };
  }

  if (!isRealOrbitalName(state.orbital)) {
    throw new RangeError(`Nom d'orbitale réelle inconnu : ${String(state.orbital)}.`);
  }
  const definition = REAL_ORBITAL_DEFINITIONS[state.orbital];
  if (!Number.isSafeInteger(state.n) || state.n < definition.minimumPrincipalQuantumNumber) {
    throw new RangeError(
      `L'orbitale ${state.orbital} requiert un entier sûr n >= ${definition.minimumPrincipalQuantumNumber}.`,
    );
  }
  return {
    degree: definition.l,
    evaluate: (thetaRadians, phiRadians) =>
      realSphericalHarmonic(state.orbital, thetaRadians, phiRadians),
  };
}

function radialDegree(state: OrbitalSamplingState): number {
  return validateAngularState(state).degree;
}

/**
 * Tabule la densité radiale P_nl(r) = r² |R_nl(r)|².
 *
 * `R_nl` est normalisée par rapport à `rBohr` dans le noyau Phase 3 ; P_nl est
 * donc le coefficient numérique d'une densité en a₀⁻¹. Cette courbe contient le
 * facteur radial r², mais aucune partie angulaire.
 *
 * @see MIT OpenCourseWare 8.04, Lecture Notes 22, sections 1 et 2.
 */
export function tabulateRadialDistribution(
  state: OrbitalSamplingState,
  extentBohr: number,
  pointCount: number,
): RadialDistributionChartData {
  requirePositiveFinite(extentBohr, "L'étendue radiale en a₀");
  requirePointCount(pointCount, 'Le nombre de points radiaux', 3);
  const degree = radialDegree(state);
  const radiiBohr = new Float64Array(pointCount);
  const rawDensityInverseBohr = new Float64Array(pointCount);
  let maximumRawDensityInverseBohr = 0;

  for (let index = 0; index < pointCount; index += 1) {
    const radiusBohr = (extentBohr * index) / (pointCount - 1);
    const radialAmplitude = hydrogenRadialWavefunction(state.n, degree, radiusBohr);
    const rawDensity = radiusBohr * radiusBohr * radialAmplitude * radialAmplitude;
    if (!Number.isFinite(rawDensity) || rawDensity < 0) {
      throw new RangeError("La densité radiale tabulée n'est pas représentable finiment.");
    }
    radiiBohr[index] = radiusBohr;
    rawDensityInverseBohr[index] = rawDensity;
    maximumRawDensityInverseBohr = Math.max(maximumRawDensityInverseBohr, rawDensity);
  }

  if (!(maximumRawDensityInverseBohr > 0)) {
    throw new RangeError(
      'La grille radiale ne contient aucune densité positive ; augmentez son étendue ou sa résolution.',
    );
  }

  const normalizedDensity = new Float64Array(pointCount);
  for (let index = 0; index < pointCount; index += 1) {
    normalizedDensity[index] =
      typedArrayValue(rawDensityInverseBohr, index, 'La densité radiale') /
      maximumRawDensityInverseBohr;
  }

  return {
    maximumRawDensityInverseBohr,
    normalizedDensity,
    rawDensityInverseBohr,
    radiiBohr,
  };
}

function directionInPlane(plane: PrincipalPlane, angleRadians: number): UnitDirection {
  const firstAxis = Math.cos(angleRadians);
  const secondAxis = Math.sin(angleRadians);
  switch (plane) {
    case 'xy':
      return { x: firstAxis, y: secondAxis, z: 0 };
    case 'xz':
      return { x: firstAxis, y: 0, z: secondAxis };
    case 'yz':
      return { x: 0, y: firstAxis, z: secondAxis };
  }
}

function sphericalAngles(direction: UnitDirection): readonly [number, number] {
  const thetaRadians = Math.acos(Math.max(-1, Math.min(1, direction.z)));
  const azimuth = Math.atan2(direction.y, direction.x);
  return [thetaRadians, azimuth < 0 ? azimuth + TWO_PI : azimuth];
}

function densityInPlane(
  evaluate: ValidatedAngularState['evaluate'],
  plane: PrincipalPlane,
  angleRadians: number,
): number {
  const [thetaRadians, phiRadians] = sphericalAngles(directionInPlane(plane, angleRadians));
  const density = complexModulusSquared(evaluate(thetaRadians, phiRadians));
  if (!Number.isFinite(density) || density < 0) {
    throw new RangeError("La densité angulaire de la coupe n'est pas représentable finiment.");
  }
  return density;
}

function meanPlaneDensity(
  evaluate: ValidatedAngularState['evaluate'],
  plane: PrincipalPlane,
  explorationPointCount: number,
): number {
  let sum = 0;
  for (let index = 0; index < explorationPointCount; index += 1) {
    sum += densityInPlane(evaluate, plane, (TWO_PI * index) / explorationPointCount);
  }
  return sum / explorationPointCount;
}

function exceedsWithRoundingMargin(candidate: number, current: number): boolean {
  const scale = Math.max(1, Math.abs(candidate), Math.abs(current));
  return candidate > current + PLANE_SCORE_ROUNDING_MARGIN * scale;
}

/**
 * Construit une coupe de la densité angulaire |Y(Ω)|² sur un grand cercle.
 *
 * Il ne s'agit pas d'une distribution probabiliste en angle : ni la mesure de
 * volume ni le facteur sin(theta) ne sont intégrés. Les trois plans principaux
 * sont explorés avec le même maillage périodique ; celui dont la densité moyenne
 * de coupe est la plus forte est retenu. Les égalités numériques sont départagées
 * dans l'ordre stable xy, xz, yz.
 *
 * @see NIST DLMF 14.30.1 pour la convention des harmoniques sphériques.
 */
export function tabulateAngularDensityCut(
  state: OrbitalSamplingState,
  pointCount: number,
  planeExplorationPointCount: number,
): AngularDensityCutChartData {
  requirePointCount(pointCount, 'Le nombre de points de la coupe angulaire', 3);
  requirePointCount(planeExplorationPointCount, "Le nombre de points d'exploration des plans", 4);
  const { evaluate } = validateAngularState(state);
  const scoreEntries = PRINCIPAL_PLANES.map(
    (plane) => [plane, meanPlaneDensity(evaluate, plane, planeExplorationPointCount)] as const,
  );
  const planeMeanDensityScoresInverseSteradians = Object.fromEntries(scoreEntries) as Record<
    PrincipalPlane,
    number
  >;
  let plane: PrincipalPlane = 'xy';
  let bestScore = planeMeanDensityScoresInverseSteradians[plane];
  for (const candidate of PRINCIPAL_PLANES.slice(1)) {
    const candidateScore = planeMeanDensityScoresInverseSteradians[candidate];
    if (exceedsWithRoundingMargin(candidateScore, bestScore)) {
      plane = candidate;
      bestScore = candidateScore;
    }
  }

  const anglesRadians = new Float64Array(pointCount);
  const rawDensityInverseSteradians = new Float64Array(pointCount);
  let maximumRawDensityInverseSteradians = 0;
  for (let index = 0; index < pointCount; index += 1) {
    const angleRadians = (TWO_PI * index) / (pointCount - 1);
    const periodicAngle = index === pointCount - 1 ? 0 : angleRadians;
    const density = densityInPlane(evaluate, plane, periodicAngle);
    anglesRadians[index] = angleRadians;
    rawDensityInverseSteradians[index] = density;
    maximumRawDensityInverseSteradians = Math.max(maximumRawDensityInverseSteradians, density);
  }

  if (!(maximumRawDensityInverseSteradians > 0)) {
    throw new RangeError('La coupe angulaire retenue ne contient aucune densité positive.');
  }

  const normalizedRadius = new Float64Array(pointCount);
  for (let index = 0; index < pointCount; index += 1) {
    normalizedRadius[index] =
      typedArrayValue(rawDensityInverseSteradians, index, 'La densité angulaire') /
      maximumRawDensityInverseSteradians;
  }

  return {
    anglesRadians,
    maximumRawDensityInverseSteradians,
    normalizedRadius,
    plane,
    planeLabel: PLANE_LABELS[plane],
    planeMeanDensityScoresInverseSteradians,
    rawDensityInverseSteradians,
  };
}

export function buildOrbitalChartData(
  state: OrbitalSamplingState,
  options: OrbitalChartDataOptions,
): OrbitalChartData {
  return {
    angularCut: tabulateAngularDensityCut(
      state,
      options.angularPointCount,
      options.planeExplorationPointCount,
    ),
    radialDistribution: tabulateRadialDistribution(
      state,
      options.radialExtentBohr,
      options.radialPointCount,
    ),
  };
}
