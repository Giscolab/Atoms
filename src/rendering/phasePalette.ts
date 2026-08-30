export type SrgbColor = readonly [red: number, green: number, blue: number];

/** Cyan utilisé lorsque seule la densité |ψ|² est affichée. */
export const DENSITY_COLOR_SRGB: SrgbColor = [0.337, 0.706, 0.843];

/** Couleur neutre réservée aux points où la phase de ψ est indéfinie. */
export const UNDEFINED_PHASE_COLOR_SRGB: SrgbColor = [0.62, 0.65, 0.68];

const FULL_TURN_RADIANS = 2 * Math.PI;

/*
 * Ancres cycliques inspirées de couleurs distinguables en déficience
 * rouge-vert : bleu ciel, vert bleuté, vermillon et mauve. La première ancre
 * est répétée implicitement à 2π, ce qui évite une rupture à la couture.
 */
const PHASE_ANCHORS: readonly SrgbColor[] = [
  [0.337, 0.706, 0.843],
  [0, 0.62, 0.451],
  [0.835, 0.369, 0],
  [0.8, 0.475, 0.655],
];

function interpolate(start: number, end: number, amount: number): number {
  return start + (end - start) * amount;
}

/** Ramène une phase finie dans [0,2π). */
export function normalizePhaseRadians(phaseRadians: number): number {
  if (!Number.isFinite(phaseRadians)) return Number.NaN;
  const remainder = phaseRadians % FULL_TURN_RADIANS;
  return remainder < 0 ? remainder + FULL_TURN_RADIANS : remainder;
}

/**
 * Convertit une phase en couleur sRGB cyclique. Cette fonction est uniquement
 * une convention de visualisation : elle ne calcule aucune grandeur physique.
 */
export function phaseColorSrgb(phaseRadians: number): SrgbColor {
  const normalizedPhase = normalizePhaseRadians(phaseRadians);
  if (Number.isNaN(normalizedPhase)) return UNDEFINED_PHASE_COLOR_SRGB;

  const scaled = (normalizedPhase / FULL_TURN_RADIANS) * PHASE_ANCHORS.length;
  const anchorIndex = Math.floor(scaled) % PHASE_ANCHORS.length;
  const nextAnchorIndex = (anchorIndex + 1) % PHASE_ANCHORS.length;
  const start = PHASE_ANCHORS[anchorIndex];
  const end = PHASE_ANCHORS[nextAnchorIndex];
  if (!start || !end) throw new RangeError('Intervalle de palette de phase invalide.');

  const localAmount = scaled - Math.floor(scaled);
  const smoothAmount = localAmount * localAmount * (3 - 2 * localAmount);
  return [
    interpolate(start[0], end[0], smoothAmount),
    interpolate(start[1], end[1], smoothAmount),
    interpolate(start[2], end[2], smoothAmount),
  ];
}
