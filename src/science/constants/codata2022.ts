export interface CodataReference {
  readonly adjustment: 'CODATA 2022';
  readonly doi: string;
  readonly title: string;
  readonly url: string;
}

export interface CodataConstant {
  readonly exactInSi: boolean;
  readonly relativeStandardUncertainty: number;
  readonly source: CodataReference;
  readonly standardUncertainty: number;
  readonly symbol: string;
  readonly unit: string;
  readonly value: number;
}

/** Valeurs recommandées issues de l'ajustement CODATA 2022, tableaux XXXII–XXXIII. */
export const CODATA_2022_REFERENCE: CodataReference = {
  adjustment: 'CODATA 2022',
  doi: '10.1063/5.0279860',
  title: 'CODATA recommended values of the fundamental physical constants: 2022',
  url: 'https://physics.nist.gov/cuu/pdf/JPCRD2022CODATA.pdf',
};

/** Vitesse de la lumière dans le vide, constante définissant le SI. */
export const SPEED_OF_LIGHT_IN_VACUUM = {
  exactInSi: true,
  relativeStandardUncertainty: 0,
  source: CODATA_2022_REFERENCE,
  standardUncertainty: 0,
  symbol: 'c',
  unit: 'm s^-1',
  value: 299_792_458,
} as const satisfies CodataConstant;

/** Constante de Planck, constante définissant le SI. */
export const PLANCK_CONSTANT = {
  exactInSi: true,
  relativeStandardUncertainty: 0,
  source: CODATA_2022_REFERENCE,
  standardUncertainty: 0,
  symbol: 'h',
  unit: 'J s',
  value: 6.626_070_15e-34,
} as const satisfies CodataConstant;

/** Constante de Planck réduite, définie exactement par ħ = h/(2π). */
export const REDUCED_PLANCK_CONSTANT = {
  exactInSi: true,
  relativeStandardUncertainty: 0,
  source: CODATA_2022_REFERENCE,
  standardUncertainty: 0,
  symbol: 'ħ',
  unit: 'J s',
  value: PLANCK_CONSTANT.value / (2 * Math.PI),
} as const satisfies CodataConstant;

/** Valeur absolue de la charge élémentaire, constante définissant le SI. */
export const ELEMENTARY_CHARGE = {
  exactInSi: true,
  relativeStandardUncertainty: 0,
  source: CODATA_2022_REFERENCE,
  standardUncertainty: 0,
  symbol: 'e',
  unit: 'C',
  value: 1.602_176_634e-19,
} as const satisfies CodataConstant;

/** Masse de l'électron, valeur recommandée et ajustée CODATA 2022. */
export const ELECTRON_MASS = {
  exactInSi: false,
  relativeStandardUncertainty: 3.1e-10,
  source: CODATA_2022_REFERENCE,
  standardUncertainty: 2.8e-40,
  symbol: 'm_e',
  unit: 'kg',
  value: 9.109_383_713_9e-31,
} as const satisfies CodataConstant;

/** Masse du proton, valeur recommandée et ajustée CODATA 2022. */
export const PROTON_MASS = {
  exactInSi: false,
  relativeStandardUncertainty: 3.1e-10,
  source: CODATA_2022_REFERENCE,
  standardUncertainty: 5.2e-37,
  symbol: 'm_p',
  unit: 'kg',
  value: 1.672_621_925_95e-27,
} as const satisfies CodataConstant;

/** Rapport de masses électron-proton directement recommandé par CODATA 2022. */
export const ELECTRON_PROTON_MASS_RATIO = {
  exactInSi: false,
  relativeStandardUncertainty: 1.7e-11,
  source: CODATA_2022_REFERENCE,
  standardUncertainty: 9.4e-15,
  symbol: 'm_e/m_p',
  unit: '1',
  value: 5.446_170_214_889e-4,
} as const satisfies CodataConstant;

/** Rayon de Bohr conventionnel, défini avec la masse de l'électron m_e. */
export const BOHR_RADIUS = {
  exactInSi: false,
  relativeStandardUncertainty: 1.6e-10,
  source: CODATA_2022_REFERENCE,
  standardUncertainty: 8.2e-21,
  symbol: 'a_0',
  unit: 'm',
  value: 5.291_772_105_44e-11,
} as const satisfies CodataConstant;

/** Énergie de Hartree conventionnelle, valeur recommandée CODATA 2022. */
export const HARTREE_ENERGY = {
  exactInSi: false,
  relativeStandardUncertainty: 1.1e-12,
  source: CODATA_2022_REFERENCE,
  standardUncertainty: 4.8e-30,
  symbol: 'E_h',
  unit: 'J',
  value: 4.359_744_722_206e-18,
} as const satisfies CodataConstant;
