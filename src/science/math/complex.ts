/** Nombre complexe minimal utilisé par le noyau scientifique pur. */
export interface Complex {
  readonly im: number;
  readonly re: number;
}

function requireFinite(value: number, quantity: string): void {
  if (!Number.isFinite(value)) {
    throw new RangeError(`${quantity} doit être un nombre fini.`);
  }
}

function requireFiniteComplex(value: Complex, quantity: string): void {
  requireFinite(value.re, `La partie réelle de ${quantity}`);
  requireFinite(value.im, `La partie imaginaire de ${quantity}`);
}

function finiteComplexResult(re: number, im: number, operation: string): Complex {
  if (!Number.isFinite(re) || !Number.isFinite(im)) {
    throw new RangeError(`Le résultat de ${operation} dépasse le domaine fini de Number.`);
  }
  return { re, im };
}

/** Crée un nombre complexe dont les deux composantes sont finies. */
export function createComplex(re: number, im: number): Complex {
  requireFinite(re, 'La partie réelle');
  requireFinite(im, 'La partie imaginaire');
  return { re, im };
}

export function addComplex(left: Complex, right: Complex): Complex {
  requireFiniteComplex(left, 'l’opérande gauche');
  requireFiniteComplex(right, 'l’opérande droit');
  return finiteComplexResult(left.re + right.re, left.im + right.im, "l'addition complexe");
}

export function subtractComplex(left: Complex, right: Complex): Complex {
  requireFiniteComplex(left, 'l’opérande gauche');
  requireFiniteComplex(right, 'l’opérande droit');
  return finiteComplexResult(left.re - right.re, left.im - right.im, 'la soustraction complexe');
}

export function multiplyComplex(left: Complex, right: Complex): Complex {
  requireFiniteComplex(left, 'l’opérande gauche');
  requireFiniteComplex(right, 'l’opérande droit');
  return finiteComplexResult(
    left.re * right.re - left.im * right.im,
    left.re * right.im + left.im * right.re,
    'la multiplication complexe',
  );
}

export function scaleComplex(value: Complex, scalar: number): Complex {
  requireFiniteComplex(value, 'la valeur complexe');
  requireFinite(scalar, 'Le scalaire');
  return finiteComplexResult(
    value.re * scalar,
    value.im * scalar,
    'la multiplication complexe par un scalaire',
  );
}

export function conjugateComplex(value: Complex): Complex {
  requireFiniteComplex(value, 'la valeur complexe');
  return { re: value.re, im: -value.im };
}

/** Retourne |z|², réel et positif ou nul. */
export function complexModulusSquared(value: Complex): number {
  requireFiniteComplex(value, 'la valeur complexe');
  const result = value.re * value.re + value.im * value.im;
  if (!Number.isFinite(result)) {
    throw new RangeError('Le module carré complexe dépasse le domaine fini de Number.');
  }
  return result;
}

/** Retourne exp(i angle) pour un angle fini exprimé en radians. */
export function complexExponential(angleRadians: number): Complex {
  requireFinite(angleRadians, "L'angle");
  return { re: Math.cos(angleRadians), im: Math.sin(angleRadians) };
}

/**
 * Retourne arg(z) dans [-π, π] en radians. La phase d'une amplitude exactement
 * nulle n'existe pas physiquement et est donc représentée par null.
 */
export function complexPhase(value: Complex): number | null {
  requireFiniteComplex(value, 'la valeur complexe');
  if (value.re === 0 && value.im === 0) return null;
  return Math.atan2(value.im, value.re);
}
