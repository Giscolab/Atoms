export type Theme = 'dark' | 'light';

export type ParsedOrbitalOption = readonly [number, number, number];

const INTEGER_TEXT_PATTERN = /^-?\d+$/u;

export function chooseTheme(storedTheme: string | null, prefersLight: boolean): Theme {
  if (storedTheme === 'dark' || storedTheme === 'light') {
    return storedTheme;
  }

  return prefersLight ? 'light' : 'dark';
}

export function parseOrbitalOption(value: string): ParsedOrbitalOption | null {
  const [first, second, third, ...additionalParts] = value.split('_');

  if (
    first === undefined ||
    second === undefined ||
    third === undefined ||
    additionalParts.length > 0 ||
    !INTEGER_TEXT_PATTERN.test(first) ||
    !INTEGER_TEXT_PATTERN.test(second) ||
    !INTEGER_TEXT_PATTERN.test(third)
  ) {
    return null;
  }

  const parsed = [Number(first), Number(second), Number(third)] as const;
  return parsed.every(Number.isSafeInteger) ? parsed : null;
}

export function inclusiveIntegerRange(min: number, max: number): number[] {
  if (!Number.isSafeInteger(min) || !Number.isSafeInteger(max)) {
    throw new RangeError('Range bounds must be safe integers.');
  }

  const values: number[] = [];
  for (let value = min; value <= max; value += 1) {
    values.push(value);
  }

  return values;
}
