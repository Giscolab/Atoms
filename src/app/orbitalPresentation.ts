import {
  hydrogenSchrodingerEnergyElectronVolts,
  HYDROGEN_REDUCED_MASS_RATIO,
} from '../science/hydrogen/energy';
import {
  hydrogenExpectedRadiusBohr,
  hydrogenNodeCounts,
  type HydrogenNodeCounts,
} from '../science/hydrogen/observables';
import {
  REAL_ORBITAL_DEFINITIONS,
  type RealOrbitalDefinition,
  type RealOrbitalName,
} from '../science/hydrogen/realOrbitals';
import { assertValidQuantumNumbers } from '../science/quantum/quantumNumbers';
import { bohrToMeters, bohrToPicometers, metersToNanometers } from '../science/units/atomicUnits';
import type { OrbitalSamplingState } from '../sampling/contracts';

const SPECTROSCOPIC_ORBITAL_LETTERS = ['s', 'p', 'd', 'f', 'g', 'h', 'i'] as const;

export interface OrbitalEnergyPresentation {
  readonly electronVolts: number;
  readonly modelLabel: string;
  readonly reducedMassRatio: number;
  readonly unit: 'eV';
}

export interface OrbitalExpectedRadiusPresentation {
  readonly bohr: number;
  readonly nanometers: number;
  readonly picometers: number;
  readonly unitLabels: {
    readonly bohr: 'a₀';
    readonly nanometers: 'nm';
    readonly picometers: 'pm';
  };
}

export interface OrbitalPresentation {
  readonly angularDegree: number;
  readonly basisLabel: string;
  readonly energy: OrbitalEnergyPresentation;
  readonly expectedRadius: OrbitalExpectedRadiusPresentation;
  readonly nodes: HydrogenNodeCounts;
  readonly notation: string;
  readonly realCombinationLabel: string | null;
}

interface ResolvedOrbitalState {
  readonly angularDegree: number;
  readonly realDefinition: RealOrbitalDefinition | null;
}

function isRealOrbitalName(value: unknown): value is RealOrbitalName {
  return typeof value === 'string' && Object.hasOwn(REAL_ORBITAL_DEFINITIONS, value);
}

function resolveOrbitalState(state: OrbitalSamplingState): ResolvedOrbitalState {
  if (state.basis === 'complex') {
    assertValidQuantumNumbers(state);
    return { angularDegree: state.l, realDefinition: null };
  }

  if (!isRealOrbitalName(state.orbital)) {
    throw new RangeError(`Nom d'orbitale réelle inconnu : ${String(state.orbital)}.`);
  }

  const realDefinition = REAL_ORBITAL_DEFINITIONS[state.orbital];
  assertValidQuantumNumbers({ n: state.n, l: realDefinition.l, m: 0 });
  return { angularDegree: realDefinition.l, realDefinition };
}

function magneticOrderLabel(order: number): string {
  return order > 0 ? `+${order}` : String(order);
}

function spectroscopicFamilyLabel(n: number, angularDegree: number): string {
  const orbitalLetter = SPECTROSCOPIC_ORBITAL_LETTERS[angularDegree];
  return orbitalLetter === undefined ? `n=${n}, l=${angularDegree}` : `${n}${orbitalLetter}`;
}

function orbitalNotation(state: OrbitalSamplingState, angularDegree: number): string {
  if (state.basis === 'real') return `${state.n}${state.orbital}`;

  const family = spectroscopicFamilyLabel(state.n, angularDegree);
  return `${family} (m = ${magneticOrderLabel(state.m)})`;
}

function realCombinationLabel(definition: RealOrbitalDefinition | null): string | null {
  if (definition === null) return null;
  if (definition.absoluteOrder === 0) return 'État réel m = 0 (harmonique zonale)';
  return `Combinaison réelle normalisée des états m = ±${definition.absoluteOrder}`;
}

/**
 * Prépare les grandeurs et libellés scientifiques destinés à l'interface.
 *
 * Cette frontière applicative délègue tous les calculs physiques et conversions
 * au noyau Phase 3. L'UI peut ainsi afficher ces valeurs sans recopier de formule.
 */
export function buildOrbitalPresentation(state: OrbitalSamplingState): OrbitalPresentation {
  const { angularDegree, realDefinition } = resolveOrbitalState(state);
  const expectedRadiusBohr = hydrogenExpectedRadiusBohr(state.n, angularDegree);
  const expectedRadiusMeters = bohrToMeters(expectedRadiusBohr);

  return {
    angularDegree,
    basisLabel: state.basis === 'complex' ? 'Base complexe |n,l,m⟩' : 'Base réelle normalisée',
    energy: {
      electronVolts: hydrogenSchrodingerEnergyElectronVolts(state.n),
      modelLabel: '¹H non relativiste — masse réduite électron-proton',
      reducedMassRatio: HYDROGEN_REDUCED_MASS_RATIO,
      unit: 'eV',
    },
    expectedRadius: {
      bohr: expectedRadiusBohr,
      nanometers: metersToNanometers(expectedRadiusMeters),
      picometers: bohrToPicometers(expectedRadiusBohr),
      unitLabels: {
        bohr: 'a₀',
        nanometers: 'nm',
        picometers: 'pm',
      },
    },
    nodes: hydrogenNodeCounts(state.n, angularDegree),
    notation: orbitalNotation(state, angularDegree),
    realCombinationLabel: realCombinationLabel(realDefinition),
  };
}
