import { describe, expect, it } from 'vitest';

import { chooseTheme, inclusiveIntegerRange, parseOrbitalOption } from '../../src/uiState';

describe('chooseTheme', () => {
  it('gives a valid stored theme priority over the system preference', () => {
    expect(chooseTheme('dark', true)).toBe('dark');
    expect(chooseTheme('light', false)).toBe('light');
  });

  it('uses the system preference when no valid theme is stored', () => {
    expect(chooseTheme(null, true)).toBe('light');
    expect(chooseTheme(null, false)).toBe('dark');
    expect(chooseTheme('unknown', true)).toBe('light');
  });
});

describe('parseOrbitalOption', () => {
  it('parses exactly three underscore-separated integers', () => {
    expect(parseOrbitalOption('2_1_-1')).toEqual([2, 1, -1]);
    expect(parseOrbitalOption('0_99_-123')).toEqual([0, 99, -123]);
  });

  it.each(['', '2_1', '2_1_0_extra', '2__0', '2_1.5_0', '2_1e0_0', '2_1_Infinity'])(
    'rejects the malformed option %j',
    (value) => {
      expect(parseOrbitalOption(value)).toBeNull();
    },
  );

  it('rejects integer text that cannot be represented safely', () => {
    expect(parseOrbitalOption('9007199254740992_0_0')).toBeNull();
  });
});

describe('inclusiveIntegerRange', () => {
  it('includes both bounds and crosses zero without gaps', () => {
    expect(inclusiveIntegerRange(-2, 2)).toEqual([-2, -1, 0, 1, 2]);
  });

  it('returns one value for equal bounds and no values for reversed bounds', () => {
    expect(inclusiveIntegerRange(4, 4)).toEqual([4]);
    expect(inclusiveIntegerRange(4, 3)).toEqual([]);
  });

  it('rejects non-integer or non-finite bounds', () => {
    expect(() => inclusiveIntegerRange(0.5, 2)).toThrow(RangeError);
    expect(() => inclusiveIntegerRange(0, Number.POSITIVE_INFINITY)).toThrow(RangeError);
  });
});
