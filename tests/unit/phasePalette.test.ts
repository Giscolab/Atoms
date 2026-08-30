import { describe, expect, it } from 'vitest';
import {
  normalizePhaseRadians,
  phaseColorSrgb,
  UNDEFINED_PHASE_COLOR_SRGB,
} from '../../src/rendering/phasePalette';

describe('palette cyclique de phase', () => {
  it('normalise les tours positifs et négatifs dans [0,2π)', () => {
    expect(normalizePhaseRadians(2 * Math.PI)).toBe(0);
    expect(normalizePhaseRadians(-Math.PI / 2)).toBeCloseTo((3 * Math.PI) / 2, 14);
  });

  it('produit exactement la même couleur après un tour complet', () => {
    const phase = 0.731;
    const reference = phaseColorSrgb(phase);
    const afterTurn = phaseColorSrgb(phase + 2 * Math.PI);
    const beforeTurn = phaseColorSrgb(phase - 2 * Math.PI);
    for (let channel = 0; channel < 3; channel += 1) {
      expect(afterTurn[channel]).toBeCloseTo(reference[channel] ?? Number.NaN, 14);
      expect(beforeTurn[channel]).toBeCloseTo(reference[channel] ?? Number.NaN, 14);
    }
  });

  it('distingue les phases opposées tout en restant dans le gamut sRGB', () => {
    const positive = phaseColorSrgb(0);
    const negative = phaseColorSrgb(Math.PI);
    expect(negative).not.toEqual(positive);
    for (const channel of [...positive, ...negative]) {
      expect(channel).toBeGreaterThanOrEqual(0);
      expect(channel).toBeLessThanOrEqual(1);
    }
  });

  it('réserve une couleur neutre aux phases non finies', () => {
    expect(phaseColorSrgb(Number.NaN)).toBe(UNDEFINED_PHASE_COLOR_SRGB);
    expect(phaseColorSrgb(Number.POSITIVE_INFINITY)).toBe(UNDEFINED_PHASE_COLOR_SRGB);
  });
});
