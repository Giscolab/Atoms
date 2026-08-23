import { assertValidQuantumNumbers } from '../quantum/quantumNumbers';
import { HYDROGEN_CHARACTERISTIC_RADIUS_BOHR } from './radialWavefunction';

export interface HydrogenNodeCounts {
  readonly radialNodes: number;
  readonly angularNodes: number;
  readonly totalNodes: number;
}

function requirePositiveFiniteResult(value: number, quantity: string): number {
  if (!Number.isFinite(value) || value <= 0) {
    throw new RangeError(`${quantity} doit être représentable par un nombre fini positif.`);
  }
  return value;
}

/**
 * Nombres de nœuds d'un état propre hydrogénoïde : n-l-1 radiaux, l angulaires
 * et n-1 au total. Le point r=0 dû au facteur r^l n'est pas compté comme un
 * nœud radial physique.
 *
 * @see https://dlmf.nist.gov/18.39.E37
 * @see MIT OpenCourseWare 8.04, Lecture Notes 22, section 2.
 */
export function hydrogenNodeCounts(n: number, l: number): HydrogenNodeCounts {
  assertValidQuantumNumbers({ n, l, m: 0 });

  return {
    radialNodes: n - l - 1,
    angularNodes: l,
    totalNodes: n - 1,
  };
}

/**
 * Valeur moyenne analytique ⟨r⟩/a₀ pour ¹H non relativiste :
 * (a_μ/a₀)[3n²-l(l+1)]/2.
 *
 * L'échelle a_μ du problème à masse réduite est celle exposée dans MIT 8.04,
 * Lecture Notes 22, section 1.
 *
 * @see https://ocw.mit.edu/courses/5-80-small-molecule-spectroscopy-and-dynamics-fall-2008/a6931596f907971b77a92bf01588bb6b_02pset_ans_sp94.pdf
 */
export function hydrogenExpectedRadiusBohr(n: number, l: number): number {
  assertValidQuantumNumbers({ n, l, m: 0 });

  const expectedRadius = (HYDROGEN_CHARACTERISTIC_RADIUS_BOHR / 2) * (3 * n * n - l * (l + 1));
  return requirePositiveFiniteResult(expectedRadius, 'La valeur moyenne ⟨r⟩/a₀');
}

/**
 * Valeur moyenne analytique a₀⟨1/r⟩ pour ¹H non relativiste :
 * 1/[(a_μ/a₀)n²]. Le résultat est le coefficient numérique en a₀⁻¹.
 *
 * @see https://ocw.mit.edu/courses/5-80-small-molecule-spectroscopy-and-dynamics-fall-2008/a6931596f907971b77a92bf01588bb6b_02pset_ans_sp94.pdf
 * @see MIT OpenCourseWare 8.04, Lecture Notes 22, section 1, pour l'échelle a_μ.
 */
export function hydrogenExpectedInverseRadiusPerBohr(n: number): number {
  assertValidQuantumNumbers({ n, l: 0, m: 0 });

  const expectedInverseRadius = 1 / (HYDROGEN_CHARACTERISTIC_RADIUS_BOHR * n * n);
  return requirePositiveFiniteResult(expectedInverseRadius, 'La valeur moyenne a₀⟨1/r⟩');
}

/**
 * Rayon le plus probable de la distribution radiale 1s, exprimé en a₀ :
 * r_max/a₀ = a_μ/a₀.
 *
 * Ce maximum résulte de r²|R_10(r)|² avec R_10 donné par DLMF 18.39.37 et
 * l'échelle de masse réduite de MIT 8.04, Lecture Notes 22, section 1.
 */
export function hydrogenMostProbableRadius1sBohr(): number {
  return requirePositiveFiniteResult(
    HYDROGEN_CHARACTERISTIC_RADIUS_BOHR,
    'Le rayon radial le plus probable de 1s',
  );
}
