import { complexPhase, scaleComplex, type Complex } from '../math/complex';
import { assertValidQuantumNumbers } from '../quantum/quantumNumbers';
import { hydrogenRadialWavefunction } from './radialWavefunction';
import { sphericalHarmonic } from './sphericalHarmonics';

/**
 * Fonction d'onde stationnaire normalisée psi_nlm = R_nl Y_l^m du modèle ¹H.
 *
 * rBohr est exprimé en rayons de Bohr conventionnels a0 ; les angles sont en
 * radians. Le résultat complexe est exprimé en a0^(-3/2).
 */
export function hydrogenWavefunction(
  n: number,
  l: number,
  m: number,
  rBohr: number,
  thetaRadians: number,
  phiRadians: number,
): Complex {
  assertValidQuantumNumbers({ n, l, m });
  const radialAmplitude = hydrogenRadialWavefunction(n, l, rBohr);
  const angularAmplitude = sphericalHarmonic(l, m, thetaRadians, phiRadians);
  return scaleComplex(angularAmplitude, radialAmplitude);
}

/**
 * Phase de psi en radians dans [-pi, pi], ou null lorsque psi est exactement nulle.
 */
export function hydrogenWavefunctionPhase(
  n: number,
  l: number,
  m: number,
  rBohr: number,
  thetaRadians: number,
  phiRadians: number,
): number | null {
  return complexPhase(hydrogenWavefunction(n, l, m, rBohr, thetaRadians, phiRadians));
}
