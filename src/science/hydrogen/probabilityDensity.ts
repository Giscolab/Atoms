import { complexModulusSquared } from '../math/complex';
import { hydrogenWavefunction } from './wavefunction';

/**
 * Densité volumique |psi_nlm|² en a0^(-3).
 *
 * Cette grandeur n'inclut ni r² ni sin(theta) : ces facteurs appartiennent à
 * la mesure sphérique r² sin(theta) dr dtheta dphi, pas à la densité elle-même.
 */
export function hydrogenProbabilityDensity(
  n: number,
  l: number,
  m: number,
  rBohr: number,
  thetaRadians: number,
  phiRadians: number,
): number {
  return complexModulusSquared(hydrogenWavefunction(n, l, m, rBohr, thetaRadians, phiRadians));
}
