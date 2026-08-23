import { BOHR_RADIUS, ELEMENTARY_CHARGE, HARTREE_ENERGY } from '../constants/codata2022';

const PICOMETERS_PER_METER = 1e12;
const NANOMETERS_PER_METER = 1e9;
const PICOMETERS_PER_BOHR = BOHR_RADIUS.value * PICOMETERS_PER_METER;
const ELECTRON_VOLTS_PER_HARTREE = HARTREE_ENERGY.value / ELEMENTARY_CHARGE.value;

function requireFinite(value: number, quantity: string): void {
  if (!Number.isFinite(value)) {
    throw new RangeError(`${quantity} doit être un nombre fini.`);
  }
}

function requireFiniteResult(value: number, quantity: string): number {
  if (!Number.isFinite(value)) {
    throw new RangeError(`La conversion de ${quantity} dépasse le domaine des nombres finis.`);
  }
  return value;
}

/** Convertit une longueur exprimée en rayons de Bohr conventionnels vers les mètres. */
export function bohrToMeters(valueInBohr: number): number {
  requireFinite(valueInBohr, 'La longueur en rayons de Bohr');
  return requireFiniteResult(valueInBohr * BOHR_RADIUS.value, 'la longueur en mètres');
}

/** Convertit une longueur exprimée en mètres vers les rayons de Bohr conventionnels. */
export function metersToBohr(valueInMeters: number): number {
  requireFinite(valueInMeters, 'La longueur en mètres');
  return requireFiniteResult(valueInMeters / BOHR_RADIUS.value, 'la longueur en rayons de Bohr');
}

/** Convertit une longueur exprimée en rayons de Bohr conventionnels vers les picomètres. */
export function bohrToPicometers(valueInBohr: number): number {
  requireFinite(valueInBohr, 'La longueur en rayons de Bohr');
  return requireFiniteResult(valueInBohr * PICOMETERS_PER_BOHR, 'la longueur en picomètres');
}

/** Convertit une longueur exprimée en picomètres vers les rayons de Bohr conventionnels. */
export function picometersToBohr(valueInPicometers: number): number {
  requireFinite(valueInPicometers, 'La longueur en picomètres');
  return requireFiniteResult(
    valueInPicometers / PICOMETERS_PER_BOHR,
    'la longueur en rayons de Bohr',
  );
}

/** Convertit une longueur exprimée en mètres vers les nanomètres. */
export function metersToNanometers(valueInMeters: number): number {
  requireFinite(valueInMeters, 'La longueur en mètres');
  return requireFiniteResult(valueInMeters * NANOMETERS_PER_METER, 'la longueur en nanomètres');
}

/** Convertit une longueur exprimée en nanomètres vers les mètres. */
export function nanometersToMeters(valueInNanometers: number): number {
  requireFinite(valueInNanometers, 'La longueur en nanomètres');
  return requireFiniteResult(valueInNanometers / NANOMETERS_PER_METER, 'la longueur en mètres');
}

/** Convertit une énergie exprimée en Hartree conventionnels vers les joules. */
export function hartreeToJoules(valueInHartree: number): number {
  requireFinite(valueInHartree, "L'énergie en Hartree");
  return requireFiniteResult(valueInHartree * HARTREE_ENERGY.value, "l'énergie en joules");
}

/** Convertit une énergie exprimée en joules vers les Hartree conventionnels. */
export function joulesToHartree(valueInJoules: number): number {
  requireFinite(valueInJoules, "L'énergie en joules");
  return requireFiniteResult(valueInJoules / HARTREE_ENERGY.value, "l'énergie en Hartree");
}

/** Convertit une énergie exprimée en Hartree vers les électronvolts. */
export function hartreeToElectronVolts(valueInHartree: number): number {
  requireFinite(valueInHartree, "L'énergie en Hartree");
  return requireFiniteResult(
    valueInHartree * ELECTRON_VOLTS_PER_HARTREE,
    "l'énergie en électronvolts",
  );
}

/** Convertit une énergie exprimée en électronvolts vers les Hartree. */
export function electronVoltsToHartree(valueInElectronVolts: number): number {
  requireFinite(valueInElectronVolts, "L'énergie en électronvolts");
  return requireFiniteResult(
    valueInElectronVolts / ELECTRON_VOLTS_PER_HARTREE,
    "l'énergie en Hartree",
  );
}
