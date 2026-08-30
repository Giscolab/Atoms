/** Version numérique du champ cartésien destiné aux représentations scientifiques. */
export const ORBITAL_FIELD_VERSION = 'hydrogen-orbital-field-v1' as const;

export interface OrbitalFieldGrid {
  /** Densité |ψ|² normalisée par son maximum sur la grille, dans [0,1]. */
  readonly densityNormalized: Float32Array;
  /** Demi-largeur physique de la grille cubique, en a₀. */
  readonly extentBohr: number;
  readonly fieldVersion: typeof ORBITAL_FIELD_VERSION;
  /** Maximum volumique brut de |ψ|², coefficient numérique en a₀⁻³. */
  readonly maximumDensityPerCubicBohr: number;
  /** Maximum brut de |ψ|, coefficient numérique en a₀⁻³ᐟ². */
  readonly maximumWavefunctionAmplitude: number;
  /** Une isosurface ψ=0 est interprétable si la fonction choisie est réelle. */
  readonly nodesAvailable: boolean;
  /** Phase de ψ sur la grille ; NaN marque uniquement une amplitude nulle. */
  readonly phaseRadians: Float32Array;
  /** Nombre de cellules par axe ; indexation x + N(y + Nz). */
  readonly resolution: number;
  /** Partie réelle de ψ normalisée par le maximum de |ψ|, dans [-1,1]. */
  readonly signedAmplitudeNormalized: Float32Array;
}

export interface OrbitalFieldRequest {
  readonly extentBohr: number;
  readonly resolution: number;
}
