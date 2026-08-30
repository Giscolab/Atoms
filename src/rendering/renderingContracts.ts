import type { OrbitalSampleSet } from '../sampling/contracts';
import type { OrbitalFieldGrid } from '../sampling/fieldContracts';

export type RenderTheme = 'dark' | 'light';
export type OrbitalObservable = 'density' | 'phase';
export type OrbitalDisplayMode = 'cloud' | 'hybrid' | 'isosurface';

export interface OrbitalRenderDataset {
  readonly field: OrbitalFieldGrid;
  readonly samples: Pick<OrbitalSampleSet, 'phaseRadians' | 'positionsBohr'>;
}

export interface OrbitalAppearance {
  readonly displayMode: OrbitalDisplayMode;
  /** Seuil d'isodensité comme fraction du maximum de la grille, dans (0,1). */
  readonly isoDensityFraction: number;
  readonly observable: OrbitalObservable;
  readonly pointOpacity: number;
  readonly pointSizePixels: number;
  readonly showAxes: boolean;
  readonly showNodes: boolean;
  readonly theme: RenderTheme;
}

export interface SceneDiagnostics {
  readonly geometries: number;
  readonly textures: number;
  readonly triangles: number;
}

export interface SceneRenderer {
  dispose(): void;
  fitCameraToOrbital(): void;
  getCameraDistance(): number;
  getDiagnostics(): SceneDiagnostics;
  hasOrbital(): boolean;
  renderFrame(): void;
  resize(viewport: HTMLElement): void;
  rotateCamera(azimuthDelta: number, elevationDelta: number): void;
  rotateCameraAutomatically(deltaRadians: number): void;
  setAppearance(appearance: OrbitalAppearance): void;
  setOrbital(dataset: OrbitalRenderDataset): void;
  zoomCamera(distanceDelta: number): void;
}
