import type { OrbitalSampleSet } from '../sampling/contracts';
import type { OrbitalFieldGrid } from '../sampling/fieldContracts';

export type RenderTheme = 'dark' | 'light';
export type OrbitalObservable = 'density' | 'phase';
export type OrbitalDisplayMode = 'cloud' | 'hybrid' | 'isosurface';

export const ORBIT_CAMERA_ELEVATION_LIMIT_RADIANS = Math.PI / 2 - 0.02;

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

export interface OrbitCameraState {
  readonly azimuthRadians: number;
  readonly distanceBohr: number;
  readonly elevationRadians: number;
}

export interface SceneDiagnostics {
  readonly cloudPoints: number;
  readonly geometries: number;
  readonly materials: number;
  readonly programs: number;
  /** Position-buffer fingerprint only: independent of phase, theme and material. */
  readonly surfaceFingerprint: string;
  readonly surfaceTriangles: number;
  readonly surfaceVertices: number;
  readonly textures: number;
  readonly triangles: number;
}

export interface SceneRenderer {
  capturePng(): Promise<Blob>;
  dispose(): void;
  fitCameraToOrbital(): void;
  getCameraDistance(): number;
  getCameraState(): OrbitCameraState;
  getDiagnostics(): SceneDiagnostics;
  hasOrbital(): boolean;
  renderFrame(): void;
  resize(viewport: HTMLElement): void;
  rotateCamera(azimuthDelta: number, elevationDelta: number): void;
  rotateCameraAutomatically(deltaRadians: number): void;
  setAppearance(appearance: OrbitalAppearance): void;
  setCameraState(camera: OrbitCameraState): void;
  setOrbital(dataset: OrbitalRenderDataset): void;
  zoomCamera(distanceDelta: number): void;
}
