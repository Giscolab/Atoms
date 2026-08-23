import { assertValidQuantumNumbers } from '../quantum/quantumNumbers';
import { factorial } from '../special/factorial';
import { generalizedLaguerre } from '../special/generalizedLaguerre';
import { HYDROGEN_REDUCED_MASS_RATIO } from './energy';

/**
 * Longueur caractéristique du problème relatif électron-proton exprimée en a₀ :
 * a_μ/a₀ = m_e/μ = 1/(μ/m_e).
 */
export const HYDROGEN_CHARACTERISTIC_RADIUS_BOHR = 1 / HYDROGEN_REDUCED_MASS_RATIO;

function requireFiniteNonNegativeRadius(rBohr: number): void {
  if (!Number.isFinite(rBohr) || rBohr < 0) {
    throw new RangeError(`Le rayon doit être fini et positif ou nul en unités de a₀ : ${rBohr}.`);
  }
}

function requireFiniteResult(value: number, quantity: string): number {
  if (!Number.isFinite(value)) {
    throw new RangeError(`${quantity} dépasse le domaine fini de Number.`);
  }
  return value;
}

/**
 * Coordonnée radiale sans dimension ρ = 2r/(n a_μ).
 *
 * `rBohr` est la valeur numérique de r/a₀ et `a_μ/a₀` vaut
 * {@link HYDROGEN_CHARACTERISTIC_RADIUS_BOHR}.
 */
export function hydrogenRadialCoordinate(n: number, rBohr: number): number {
  assertValidQuantumNumbers({ n, l: 0, m: 0 });
  requireFiniteNonNegativeRadius(rBohr);

  const rho = (rBohr / n) * (2 / HYDROGEN_CHARACTERISTIC_RADIUS_BOHR);
  return requireFiniteResult(rho, 'La coordonnée radiale ρ');
}

/**
 * Fonction radiale hydrogène normalisée R_nl(r), avec masse réduite.
 *
 * Le rayon d'entrée `rBohr` est exprimé en a₀. La valeur retournée est donc le
 * coefficient numérique de R_nl en a₀⁻³ᐟ² et vérifie
 * ∫₀∞ rBohr² |R_nl(rBohr)|² drBohr = 1.
 *
 * Convention : ρ = 2r/(n a_μ), k = n-l-1 et α = 2l+1. La phase globale est
 * celle de la forme radiale DLMF, avec les Laguerre généralisés de convention DLMF.
 *
 * @see https://dlmf.nist.gov/18.39.E37
 * @see MIT OpenCourseWare 8.04, Lecture Notes 22, sections 1 et 2.
 */
export function hydrogenRadialWavefunction(n: number, l: number, rBohr: number): number {
  assertValidQuantumNumbers({ n, l, m: 0 });
  const rho = hydrogenRadialCoordinate(n, rBohr);
  const degree = n - l - 1;
  const alpha = 2 * l + 1;

  const factorialRatio = factorial(degree) / factorial(n + l);
  const inverseScaledRadius = 2 / (n * HYDROGEN_CHARACTERISTIC_RADIUS_BOHR);
  const normalizationSquared = (inverseScaledRadius ** 3 * factorialRatio) / (2 * n);

  if (!Number.isFinite(normalizationSquared) || normalizationSquared <= 0) {
    throw new RangeError('La normalisation radiale ne peut pas être représentée finiment.');
  }

  const normalization = Math.sqrt(normalizationSquared);
  const exponential = Math.exp(-rho / 2);
  const radialPower = rho ** l;
  const laguerre = generalizedLaguerre(degree, alpha, rho);

  requireFiniteResult(radialPower, 'Le facteur radial ρ^l');
  return requireFiniteResult(
    normalization * exponential * radialPower * laguerre,
    'La fonction radiale R_nl',
  );
}
