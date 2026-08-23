export interface QuantumNumbers {
  readonly l: number;
  readonly m: number;
  readonly n: number;
}

export function isValidPrincipalQuantumNumber(n: number): boolean {
  return Number.isSafeInteger(n) && n >= 1;
}

export function assertValidPrincipalQuantumNumber(n: number): void {
  if (!isValidPrincipalQuantumNumber(n)) {
    throw new RangeError(
      `Le nombre quantique principal n doit être un entier sûr supérieur à 0 : ${n}.`,
    );
  }
}

/** Vérifie n >= 1, 0 <= l <= n - 1 et -l <= m <= l pour des entiers sûrs. */
export function isValidQuantumNumbers(value: unknown): value is QuantumNumbers {
  if (typeof value !== 'object' || value === null) return false;

  const candidate = value as Partial<Record<keyof QuantumNumbers, unknown>>;
  const { n, l, m } = candidate;
  if (typeof n !== 'number' || typeof l !== 'number' || typeof m !== 'number') return false;

  return (
    isValidPrincipalQuantumNumber(n) &&
    Number.isSafeInteger(l) &&
    Number.isSafeInteger(m) &&
    l >= 0 &&
    l <= n - 1 &&
    m >= -l &&
    m <= l
  );
}

export function assertValidQuantumNumbers(value: unknown): asserts value is QuantumNumbers {
  if (!isValidQuantumNumbers(value)) {
    throw new RangeError(
      'Les nombres quantiques doivent respecter n >= 1, 0 <= l < n et |m| <= l.',
    );
  }
}
