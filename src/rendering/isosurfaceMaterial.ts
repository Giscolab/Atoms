import * as THREE from 'three';
import type { OrbitalDisplayMode, RenderTheme } from './renderingContracts';

/** Artistic surface settings, never properties of the atom or of its charge. */
export function createIsosurfaceMaterial(
  mode: OrbitalDisplayMode,
  theme: RenderTheme,
): THREE.MeshPhysicalMaterial {
  const hybrid = mode === 'hybrid';
  return new THREE.MeshPhysicalMaterial({
    vertexColors: true,
    // Only the outward-facing skin contributes: no stacked front/back veils.
    side: THREE.FrontSide,
    transparent: hybrid,
    opacity: hybrid ? 0.58 : 1,
    depthWrite: !hybrid,
    metalness: 0,
    // A broad matte highlight keeps the scientific phase colors dominant.
    roughness: theme === 'dark' ? 0.5 : 0.46,
    clearcoat: 0.1,
    clearcoatRoughness: 0.62,
    // Neutral highlights preserve the meaning of the cyclic phase palette.
    specularIntensity: theme === 'dark' ? 0.42 : 0.34,
    emissive: 0x000000,
    transmission: 0,
    flatShading: false,
  });
}

/** A small camera-oriented studio rig; no textures, shadows or postprocess passes. */
export function createIsosurfaceLighting(): {
  readonly rig: THREE.Group;
  setTheme(theme: RenderTheme): void;
} {
  const rig = new THREE.Group();
  const ambient = new THREE.HemisphereLight(0xffffff, 0x24323b, 0.48);
  const key = new THREE.DirectionalLight(0xffffff, 2.1);
  key.position.set(-3, 4, 5);
  const fill = new THREE.DirectionalLight(0xffffff, 0.45);
  fill.position.set(4, -1, 2);
  const rim = new THREE.DirectionalLight(0xffffff, 1.7);
  rim.position.set(2, 3, -4);
  rig.add(ambient, key, fill, rim);
  return {
    rig,
    setTheme(theme): void {
      ambient.intensity = theme === 'dark' ? 0.48 : 0.72;
      key.intensity = theme === 'dark' ? 2.1 : 1.8;
      fill.intensity = theme === 'dark' ? 0.45 : 0.58;
      rim.intensity = theme === 'dark' ? 1.7 : 1.3;
    },
  };
}
