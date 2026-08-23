import { BOHR_RADIUS, ELECTRON_MASS, ELECTRON_PROTON_MASS_RATIO } from '../constants/codata2022';
import { assertValidPrincipalQuantumNumber } from '../quantum/quantumNumbers';
import { hartreeToElectronVolts, hartreeToJoules } from '../units/atomicUnits';

function requirePositiveFiniteMass(massKilograms: number, name: string): void {
  if (!Number.isFinite(massKilograms) || massKilograms <= 0) {
    throw new RangeError(
      `${name} doit être une masse strictement positive et finie en kilogrammes.`,
    );
  }
}

/** Calcule μ = m1 m2 / (m1 + m2) pour deux masses positives exprimées en kilogrammes. */
export function reducedMassKilograms(mass1Kilograms: number, mass2Kilograms: number): number {
  requirePositiveFiniteMass(mass1Kilograms, 'La première masse');
  requirePositiveFiniteMass(mass2Kilograms, 'La seconde masse');

  const smallerMass = Math.min(mass1Kilograms, mass2Kilograms);
  const largerMass = Math.max(mass1Kilograms, mass2Kilograms);
  // Forme équivalente à m1*m2/(m1+m2), sans produit ni somme intermédiaire
  // susceptibles de déborder ou de sous-déborder pour des masses finies extrêmes.
  const reducedMass = smallerMass / (1 + smallerMass / largerMass);
  if (!Number.isFinite(reducedMass) || reducedMass <= 0) {
    throw new RangeError(
      'La masse réduite ne peut pas être représentée comme un nombre fini positif.',
    );
  }
  return reducedMass;
}

const electronProtonMassRatio = ELECTRON_PROTON_MASS_RATIO.value;

/** Rapport μ/m_e = 1/(1 + m_e/m_p) fondé sur le rapport direct CODATA 2022. */
export const HYDROGEN_REDUCED_MASS_RATIO = 1 / (1 + electronProtonMassRatio);

/** Masse réduite centrale électron-proton du modèle non relativiste de ¹H. */
export const HYDROGEN_REDUCED_MASS_KILOGRAMS = ELECTRON_MASS.value * HYDROGEN_REDUCED_MASS_RATIO;

/** Rayon caractéristique a_μ = a_0 m_e/μ du problème relatif électron-proton. */
export const HYDROGEN_CHARACTERISTIC_RADIUS_METERS =
  BOHR_RADIUS.value / HYDROGEN_REDUCED_MASS_RATIO;

/**
 * Énergie coulombienne non relativiste de ¹H en Hartree conventionnels :
 * E_n/E_h = -(μ/m_e)/(2 n²).
 */
export function hydrogenSchrodingerEnergyHartree(n: number): number {
  assertValidPrincipalQuantumNumber(n);
  return -HYDROGEN_REDUCED_MASS_RATIO / (2 * n * n);
}

/** Énergie coulombienne non relativiste de ¹H en joules. */
export function hydrogenSchrodingerEnergyJoules(n: number): number {
  return hartreeToJoules(hydrogenSchrodingerEnergyHartree(n));
}

/** Énergie coulombienne non relativiste de ¹H en électronvolts. */
export function hydrogenSchrodingerEnergyElectronVolts(n: number): number {
  return hartreeToElectronVolts(hydrogenSchrodingerEnergyHartree(n));
}
