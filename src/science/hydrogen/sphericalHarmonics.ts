import { complexExponential, conjugateComplex, scaleComplex, type Complex } from '../math/complex';
import { assertValidQuantumNumbers } from '../quantum/quantumNumbers';
import { associatedLegendre } from '../special/associatedLegendre';
import { factorial } from '../special/factorial';

const TWO_PI = 2 * Math.PI;

function validateAngularArguments(
  degree: number,
  order: number,
  thetaRadians: number,
  phiRadians: number,
): void {
  // n=l+1 permet de réutiliser exactement les contraintes centrales sur l et m.
  assertValidQuantumNumbers({ n: degree + 1, l: degree, m: order });
  if (!Number.isFinite(thetaRadians) || thetaRadians < 0 || thetaRadians > Math.PI) {
    throw new RangeError(`theta doit être fini et appartenir à [0, pi] : ${thetaRadians}.`);
  }
  if (!Number.isFinite(phiRadians)) {
    throw new RangeError(`phi doit être un angle fini en radians : ${phiRadians}.`);
  }
}

function positiveOrderSphericalHarmonic(
  degree: number,
  order: number,
  thetaRadians: number,
  phiRadians: number,
): Complex {
  const normalization = Math.sqrt(
    ((2 * degree + 1) / (4 * Math.PI)) * (factorial(degree - order) / factorial(degree + order)),
  );
  if (!Number.isFinite(normalization) || normalization <= 0) {
    throw new RangeError(`La normalisation de Y_${degree}^${order} n'est pas représentable.`);
  }

  const ferrers = associatedLegendre(degree, order, Math.cos(thetaRadians));
  // Réduction périodique avant la multiplication par m : phi reste accepté sur
  // tout le domaine fini sans risque d'overflow pour m*phi.
  const periodicPhi = phiRadians % TWO_PI;
  const azimuthalPhase = complexExponential(order * periodicPhi);
  return scaleComplex(azimuthalPhase, normalization * ferrers);
}

/**
 * Harmonique sphérique complexe normalisée Y_l^m(theta, phi), convention DLMF 14.30.1.
 *
 * `associatedLegendre` contient déjà la phase de Condon–Shortley (-1)^m : aucun
 * facteur supplémentaire n'est appliqué ici. theta est dans [0, pi], phi est tout
 * réel fini en radians, et le résultat est sans dimension.
 */
export function sphericalHarmonic(
  degree: number,
  order: number,
  thetaRadians: number,
  phiRadians: number,
): Complex {
  validateAngularArguments(degree, order, thetaRadians, phiRadians);
  if (order >= 0) {
    return positiveOrderSphericalHarmonic(degree, order, thetaRadians, phiRadians);
  }

  const positiveOrder = -order;
  const positiveHarmonic = positiveOrderSphericalHarmonic(
    degree,
    positiveOrder,
    thetaRadians,
    phiRadians,
  );
  const negativeOrderSign = positiveOrder % 2 === 0 ? 1 : -1;
  // DLMF 14.30.6 : Y_l^(-m) = (-1)^m conjugate(Y_l^m), m > 0.
  return scaleComplex(conjugateComplex(positiveHarmonic), negativeOrderSign);
}
