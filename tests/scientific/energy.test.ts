import { describe, expect, it } from 'vitest';

import {
  BOHR_RADIUS,
  ELECTRON_MASS,
  ELECTRON_PROTON_MASS_RATIO,
} from '../../src/science/constants/codata2022';
import {
  HYDROGEN_CHARACTERISTIC_RADIUS_METERS,
  HYDROGEN_REDUCED_MASS_KILOGRAMS,
  HYDROGEN_REDUCED_MASS_RATIO,
  hydrogenSchrodingerEnergyElectronVolts,
  hydrogenSchrodingerEnergyHartree,
  hydrogenSchrodingerEnergyJoules,
  reducedMassKilograms,
} from '../../src/science/hydrogen/energy';
import { hartreeToElectronVolts, hartreeToJoules } from '../../src/science/units/atomicUnits';
import { relativeError, SHORT_FLOAT64_RELATIVE_TOLERANCE } from './numericAssertions';

describe('modèle énergétique non relativiste de ¹H', () => {
  it('utilise une masse réduite positive et strictement inférieure à m_e', () => {
    expect(HYDROGEN_REDUCED_MASS_KILOGRAMS).toBeGreaterThan(0);
    expect(HYDROGEN_REDUCED_MASS_KILOGRAMS).toBeLessThan(ELECTRON_MASS.value);
    expect(HYDROGEN_REDUCED_MASS_RATIO).toBeGreaterThan(0);
    expect(HYDROGEN_REDUCED_MASS_RATIO).toBeLessThan(1);
    expect(HYDROGEN_REDUCED_MASS_RATIO).toBe(1 / (1 + ELECTRON_PROTON_MASS_RATIO.value));
  });

  it('distingue a_mu du rayon de Bohr conventionnel a0', () => {
    expect(HYDROGEN_CHARACTERISTIC_RADIUS_METERS).toBeGreaterThan(BOHR_RADIUS.value);
    expect(
      relativeError(HYDROGEN_CHARACTERISTIC_RADIUS_METERS, 5.294_654_094_602_462_5e-11),
    ).toBeLessThanOrEqual(SHORT_FLOAT64_RELATIVE_TOLERANCE);
  });

  it('produit des niveaux liés négatifs avec la dépendance exacte en 1/n²', () => {
    const energy1 = hydrogenSchrodingerEnergyHartree(1);
    const energy2 = hydrogenSchrodingerEnergyHartree(2);
    const energy3 = hydrogenSchrodingerEnergyHartree(3);

    expect(energy1).toBeLessThan(0);
    expect(energy2).toBeLessThan(0);
    expect(energy3).toBeLessThan(0);
    expect(relativeError(Math.abs(energy2), Math.abs(energy1) / 4)).toBeLessThanOrEqual(
      SHORT_FLOAT64_RELATIVE_TOLERANCE,
    );
    expect(relativeError(Math.abs(energy3), Math.abs(energy1) / 9)).toBeLessThanOrEqual(
      SHORT_FLOAT64_RELATIVE_TOLERANCE,
    );
  });

  it('retrouve la valeur centrale du modèle, sans la confondre avec une valeur spectroscopique', () => {
    expect(relativeError(HYDROGEN_REDUCED_MASS_RATIO, 0.999_455_679_424_761_5)).toBeLessThanOrEqual(
      SHORT_FLOAT64_RELATIVE_TOLERANCE,
    );
    expect(
      relativeError(HYDROGEN_REDUCED_MASS_KILOGRAMS, 9.104_425_288_916_783e-31),
    ).toBeLessThanOrEqual(SHORT_FLOAT64_RELATIVE_TOLERANCE);
    expect(
      relativeError(hydrogenSchrodingerEnergyHartree(1), -0.499_727_839_712_380_8),
    ).toBeLessThanOrEqual(SHORT_FLOAT64_RELATIVE_TOLERANCE);
    expect(
      relativeError(hydrogenSchrodingerEnergyJoules(1), -2.178_685_811_725_458e-18),
    ).toBeLessThanOrEqual(SHORT_FLOAT64_RELATIVE_TOLERANCE);
    expect(
      relativeError(hydrogenSchrodingerEnergyElectronVolts(1), -13.598_287_264_283_359),
    ).toBeLessThanOrEqual(SHORT_FLOAT64_RELATIVE_TOLERANCE);
  });

  it('exprime la même énergie en Hartree, joules et électronvolts', () => {
    const energyHartree = hydrogenSchrodingerEnergyHartree(3);
    expect(
      relativeError(hydrogenSchrodingerEnergyJoules(3), hartreeToJoules(energyHartree)),
    ).toBeLessThanOrEqual(SHORT_FLOAT64_RELATIVE_TOLERANCE);
    expect(
      relativeError(
        hydrogenSchrodingerEnergyElectronVolts(3),
        hartreeToElectronVolts(energyHartree),
      ),
    ).toBeLessThanOrEqual(SHORT_FLOAT64_RELATIVE_TOLERANCE);
  });

  it('rejette les masses et nombres quantiques invalides', () => {
    expect(() => reducedMassKilograms(0, 1)).toThrow(RangeError);
    expect(() => reducedMassKilograms(1, Number.NaN)).toThrow(RangeError);
    expect(() => hydrogenSchrodingerEnergyHartree(0)).toThrow(RangeError);
    expect(() => hydrogenSchrodingerEnergyJoules(1.5)).toThrow(RangeError);
  });

  it('calcule une masse réduite représentable sans overflow ou underflow intermédiaire', () => {
    expect(reducedMassKilograms(1e308, 1e308)).toBe(5e307);
    expect(reducedMassKilograms(1e-300, 1e-300)).toBe(5e-301);
  });
});
