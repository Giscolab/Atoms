import { describe, expect, it } from 'vitest';

import {
  bohrToMeters,
  bohrToPicometers,
  electronVoltsToHartree,
  hartreeToElectronVolts,
  hartreeToJoules,
  joulesToHartree,
  metersToBohr,
  metersToNanometers,
  nanometersToMeters,
  picometersToBohr,
} from '../../src/science/units/atomicUnits';
import { relativeError, SHORT_FLOAT64_RELATIVE_TOLERANCE } from './numericAssertions';

describe('conversions des unités atomiques', () => {
  it('préserve une longueur lors des allers-retours a0 ↔ m et a0 ↔ pm', () => {
    const valueInBohr = 2.75;
    expect(relativeError(metersToBohr(bohrToMeters(valueInBohr)), valueInBohr)).toBeLessThanOrEqual(
      SHORT_FLOAT64_RELATIVE_TOLERANCE,
    );
    expect(
      relativeError(picometersToBohr(bohrToPicometers(valueInBohr)), valueInBohr),
    ).toBeLessThanOrEqual(SHORT_FLOAT64_RELATIVE_TOLERANCE);
  });

  it('établit avec la valeur CODATA que 1 a0 vaut environ 52,9 pm', () => {
    const codataPicometers = 52.917_721_054_4;
    expect(relativeError(bohrToPicometers(1), codataPicometers)).toBeLessThanOrEqual(
      SHORT_FLOAT64_RELATIVE_TOLERANCE,
    );
  });

  it('préserve une longueur lors de la conversion m → nm → m', () => {
    const valueInMeters = 2.5e-9;
    expect(
      relativeError(nanometersToMeters(metersToNanometers(valueInMeters)), valueInMeters),
    ).toBeLessThanOrEqual(SHORT_FLOAT64_RELATIVE_TOLERANCE);
  });

  it('préserve une énergie lors des allers-retours Hartree ↔ J et Hartree ↔ eV', () => {
    const valueInHartree = -3.25;
    expect(
      relativeError(joulesToHartree(hartreeToJoules(valueInHartree)), valueInHartree),
    ).toBeLessThanOrEqual(SHORT_FLOAT64_RELATIVE_TOLERANCE);
    expect(
      relativeError(electronVoltsToHartree(hartreeToElectronVolts(valueInHartree)), valueInHartree),
    ).toBeLessThanOrEqual(SHORT_FLOAT64_RELATIVE_TOLERANCE);
  });

  it('retrouve la conversion CODATA du Hartree vers les électronvolts', () => {
    const codataElectronVolts = 27.211_386_245_981;
    // Le tableau XXXIII donne 27.211 386 245 981(30) eV : l'incertitude-type
    // absolue de cette représentation publiée vaut donc 3,0 × 10^-11 eV.
    const codataStandardUncertaintyElectronVolts = 3e-11;
    expect(Math.abs(hartreeToElectronVolts(1) - codataElectronVolts)).toBeLessThanOrEqual(
      codataStandardUncertaintyElectronVolts,
    );
  });

  it('rejette les entrées non finies au lieu de produire silencieusement NaN', () => {
    expect(() => bohrToMeters(Number.NaN)).toThrow(RangeError);
    expect(() => metersToNanometers(Number.POSITIVE_INFINITY)).toThrow(RangeError);
    expect(() => electronVoltsToHartree(Number.NEGATIVE_INFINITY)).toThrow(RangeError);
  });

  it('échoue explicitement lorsqu’une conversion dépasse le domaine fini de Number', () => {
    expect(() => bohrToPicometers(Number.MAX_VALUE)).toThrow(RangeError);
  });

  it('préserve les résultats subnormaux encore représentables sans intermédiaire sous-débordé', () => {
    expect(bohrToPicometers(1e-314)).toBeGreaterThan(0);
    expect(picometersToBohr(1e-312)).toBeGreaterThan(0);
    expect(hartreeToElectronVolts(1e-307)).toBeGreaterThan(0);
    expect(electronVoltsToHartree(1e-307)).toBeGreaterThan(0);
  });
});
