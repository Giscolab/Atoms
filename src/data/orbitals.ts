import type { NamedSamplerKey } from '../science/legacyScience';

export const SAMPLER_MAP: Partial<Record<string, NamedSamplerKey>> = {
  '2_1_0': '2p_z',
  '2_1_1': '2p_x',
  '2_1_-1': '2p_y',
  '3_1_0': '3p_z',
};

const ORBITAL_LETTERS: string[] = ['s', 'p', 'd', 'f', 'g', 'h', 'i'];

const ORBITAL_FAMILY: Record<string, string> = {
  s: 'Distribution sphérique autour du noyau.',
  p: 'Deux lobes principaux avec un plan nodal.',
  d: 'Forme complexe a plusieurs lobes (type trèfle).',
  f: 'Forme multi-lobes avec structure plus fine.',
  g: 'Orbitale de haut ordre avec nombreux noeuds angulaires.',
  h: 'Orbitale de haut ordre, très structurée.',
  i: 'Orbitale de tres haut ordre et forte complexité angulaire.',
};

export const PANEL_MIN_N = 1;
export const PANEL_MAX_N = 9;

export interface OrbitalDescription {
  letter: string;
  familySummary: string;
  magneticProjectionSummary: string;
  radialNodes: number;
  angularNodes: number;
}

export function describeOrbital(n: number, l: number, m: number): OrbitalDescription {
  const letter = ORBITAL_LETTERS[l] || `l=${l}`;
  const familySummary = ORBITAL_FAMILY[letter] || `Famille orbitale definie par l=${l}.`;
  let magneticProjectionSummary: string;
  if (m === 0) magneticProjectionSummary = 'Symetrie axiale (m=0).';
  else if (m > 0) magneticProjectionSummary = `Projection positive (m=+${m}).`;
  else magneticProjectionSummary = `Projection negative (m=${m}).`;

  return {
    letter,
    familySummary,
    magneticProjectionSummary,
    radialNodes: Math.max(0, n - l - 1),
    angularNodes: l,
  };
}
