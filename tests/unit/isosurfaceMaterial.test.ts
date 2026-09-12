import { describe, expect, it } from 'vitest';
import { FrontSide, Light } from 'three';
import {
  createIsosurfaceLighting,
  createIsosurfaceMaterial,
} from '../../src/rendering/isosurfaceMaterial';

describe('isodensity presentation', () => {
  it.each(['dark', 'light'] as const)('has an opaque, depth-tested skin in %s', (theme) => {
    const material = createIsosurfaceMaterial('isosurface', theme);
    expect(material.transparent).toBe(false);
    expect(material.opacity).toBe(1);
    expect(material.depthWrite).toBe(true);
    expect(material.depthTest).toBe(true);
    expect(material.side).toBe(FrontSide);
    expect(material.vertexColors).toBe(true);
    expect(material.flatShading).toBe(false);
    expect(material.transmission).toBe(0);
    expect(material.emissive.getHex()).toBe(0);
    expect(material.metalness).toBe(0);
    expect(material.roughness).toBeGreaterThanOrEqual(0.46);
    expect(material.clearcoat).toBeLessThanOrEqual(0.1);
    expect(material.specularIntensity).toBeLessThanOrEqual(0.42);
    material.dispose();
  });

  it('leaves cloud samples visible through a single hybrid skin', () => {
    const material = createIsosurfaceMaterial('hybrid', 'dark');
    expect(material.transparent).toBe(true);
    expect(material.opacity).toBeGreaterThan(0.5);
    expect(material.opacity).toBeLessThan(0.85);
    expect(material.side).toBe(FrontSide);
    expect(material.depthWrite).toBe(false);
    material.dispose();
  });

  it('keeps bounded, neutral studio lights and reuses them on theme changes', () => {
    const lighting = createIsosurfaceLighting();
    const children = [...lighting.rig.children];
    for (const theme of ['light', 'dark', 'light'] as const) {
      lighting.setTheme(theme);
      expect(lighting.rig.children).toEqual(children);
      for (const light of lighting.rig.children) {
        expect(light).toBeInstanceOf(Light);
        if (!(light instanceof Light)) throw new Error('Expected a light');
        expect(light.intensity).toBeGreaterThan(0);
        expect(light.intensity).toBeLessThan(4);
        expect(light.color.getHex()).toBe(0xffffff);
      }
    }
  });
});
