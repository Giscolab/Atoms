import { describe, expect, it } from 'vitest';

import {
  BOHR_RADIUS,
  CODATA_2022_REFERENCE,
  ELECTRON_MASS,
  ELECTRON_PROTON_MASS_RATIO,
  ELEMENTARY_CHARGE,
  HARTREE_ENERGY,
  PLANCK_CONSTANT,
  PROTON_MASS,
  REDUCED_PLANCK_CONSTANT,
  SPEED_OF_LIGHT_IN_VACUUM,
} from '../../src/science/constants/codata2022';

describe('constantes CODATA 2022', () => {
  it('identifie explicitement la version et la publication de référence', () => {
    expect(CODATA_2022_REFERENCE.adjustment).toBe('CODATA 2022');
    expect(CODATA_2022_REFERENCE.doi).toBe('10.1063/5.0279860');
    expect(CODATA_2022_REFERENCE.url).toBe('https://physics.nist.gov/cuu/pdf/JPCRD2022CODATA.pdf');
  });

  it('représente les constantes définissant le SI comme exactes', () => {
    expect(SPEED_OF_LIGHT_IN_VACUUM).toMatchObject({
      exactInSi: true,
      relativeStandardUncertainty: 0,
      standardUncertainty: 0,
      unit: 'm s^-1',
      value: 299_792_458,
    });
    expect(PLANCK_CONSTANT).toMatchObject({
      exactInSi: true,
      relativeStandardUncertainty: 0,
      standardUncertainty: 0,
      unit: 'J s',
      value: 6.626_070_15e-34,
    });
    expect(ELEMENTARY_CHARGE).toMatchObject({
      exactInSi: true,
      relativeStandardUncertainty: 0,
      standardUncertainty: 0,
      unit: 'C',
      value: 1.602_176_634e-19,
    });
    expect(REDUCED_PLANCK_CONSTANT).toMatchObject({
      exactInSi: true,
      relativeStandardUncertainty: 0,
      standardUncertainty: 0,
      unit: 'J s',
      value: PLANCK_CONSTANT.value / (2 * Math.PI),
    });
  });

  it('conserve les valeurs centrales et incertitudes publiées des constantes non exactes', () => {
    expect(ELECTRON_MASS).toMatchObject({
      exactInSi: false,
      relativeStandardUncertainty: 3.1e-10,
      standardUncertainty: 2.8e-40,
      unit: 'kg',
      value: 9.109_383_713_9e-31,
    });
    expect(PROTON_MASS).toMatchObject({
      exactInSi: false,
      relativeStandardUncertainty: 3.1e-10,
      standardUncertainty: 5.2e-37,
      unit: 'kg',
      value: 1.672_621_925_95e-27,
    });
    expect(ELECTRON_PROTON_MASS_RATIO).toMatchObject({
      exactInSi: false,
      relativeStandardUncertainty: 1.7e-11,
      standardUncertainty: 9.4e-15,
      unit: '1',
      value: 5.446_170_214_889e-4,
    });
    expect(BOHR_RADIUS).toMatchObject({
      exactInSi: false,
      relativeStandardUncertainty: 1.6e-10,
      standardUncertainty: 8.2e-21,
      unit: 'm',
      value: 5.291_772_105_44e-11,
    });
    expect(HARTREE_ENERGY).toMatchObject({
      exactInSi: false,
      relativeStandardUncertainty: 1.1e-12,
      standardUncertainty: 4.8e-30,
      unit: 'J',
      value: 4.359_744_722_206e-18,
    });
  });

  it('rattache chaque constante à la même publication CODATA 2022', () => {
    const constants = [
      SPEED_OF_LIGHT_IN_VACUUM,
      PLANCK_CONSTANT,
      REDUCED_PLANCK_CONSTANT,
      ELEMENTARY_CHARGE,
      ELECTRON_MASS,
      PROTON_MASS,
      ELECTRON_PROTON_MASS_RATIO,
      BOHR_RADIUS,
      HARTREE_ENERGY,
    ];

    for (const constant of constants) {
      expect(constant.source).toBe(CODATA_2022_REFERENCE);
    }
  });
});
