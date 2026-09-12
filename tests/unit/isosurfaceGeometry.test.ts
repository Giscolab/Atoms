import * as THREE from 'three';
import { describe, expect, it } from 'vitest';
import { ORBITAL_FIELD_VERSION, type OrbitalFieldGrid } from '../../src/sampling/fieldContracts';
import { computeOrbitalField } from '../../src/sampling/orbitalField';
import { createIsosurfaceGeometry } from '../../src/rendering/isosurfaceGeometry';
import { phaseColorSrgb, UNDEFINED_PHASE_COLOR_SRGB } from '../../src/rendering/phasePalette';

function fieldFrom(
  n: number,
  density: (x: number, y: number, z: number) => number,
): OrbitalFieldGrid {
  const values = new Float32Array(n ** 3);
  for (let z = 0; z < n; z += 1)
    for (let y = 0; y < n; y += 1)
      for (let x = 0; x < n; x += 1) {
        values[x + n * (y + n * z)] = density(
          (4 * x) / (n - 1) - 2,
          (4 * y) / (n - 1) - 2,
          (4 * z) / (n - 1) - 2,
        );
      }
  return {
    densityNormalized: values,
    extentBohr: 2,
    fieldVersion: ORBITAL_FIELD_VERSION,
    maximumDensityPerCubicBohr: 1,
    maximumWavefunctionAmplitude: 1,
    nodesAvailable: false,
    phaseRadians: new Float32Array(values.length),
    resolution: n,
    signedAmplitudeNormalized: Float32Array.from(values, Math.sqrt),
  };
}

function assertMesh(geometry: THREE.BufferGeometry): void {
  const positions = geometry.getAttribute('position');
  const normals = geometry.getAttribute('normal');
  expect(positions.count).toBeGreaterThan(0);
  expect(positions.count % 3).toBe(0);
  const a = new THREE.Vector3();
  const b = new THREE.Vector3();
  const c = new THREE.Vector3();
  const face = new THREE.Vector3();
  const edge = new THREE.Vector3();
  const normal = new THREE.Vector3();
  for (let index = 0; index < positions.count; index += 3) {
    a.fromBufferAttribute(positions, index);
    b.fromBufferAttribute(positions, index + 1);
    c.fromBufferAttribute(positions, index + 2);
    face.subVectors(b, a).cross(edge.subVectors(c, a));
    expect(face.lengthSq()).toBeGreaterThan(0);
    expect(Number.isFinite(face.lengthSq())).toBe(true);
    normal.fromBufferAttribute(normals, index);
    expect(face.dot(normal)).toBeGreaterThan(0);
    for (let vertex = index; vertex < index + 3; vertex += 1) {
      normal.fromBufferAttribute(normals, vertex);
      expect(normal.length()).toBeCloseTo(1, 5);
    }
  }
}

function closedComponentCount(geometry: THREE.BufferGeometry): number {
  const position = geometry.getAttribute('position');
  const vertices = new Map<string, number>();
  const parents: number[] = [];
  const edges = new Map<string, number>();
  const root = (value: number): number => {
    let result = value;
    while (parents[result] !== result) result = parents[result] ?? result;
    return result;
  };
  const id = (index: number): number => {
    const key = [position.getX(index), position.getY(index), position.getZ(index)]
      .map((value) => Math.round(value * 1e5))
      .join(',');
    let value = vertices.get(key);
    if (value === undefined) {
      value = vertices.size;
      vertices.set(key, value);
      parents.push(value);
    }
    return value;
  };
  for (let index = 0; index < position.count; index += 3) {
    const triangle = [id(index), id(index + 1), id(index + 2)] as const;
    for (const [a, b] of [
      [triangle[0], triangle[1]],
      [triangle[1], triangle[2]],
      [triangle[2], triangle[0]],
    ] as const) {
      parents[root(b)] = root(a);
      const key = `${Math.min(a, b)},${Math.max(a, b)}`;
      edges.set(key, (edges.get(key) ?? 0) + 1);
    }
  }
  expect([...edges.values()].every((count) => count === 2)).toBe(true);
  return new Set(parents.map((_, index) => root(index))).size;
}

describe('density isosurface tessellation', () => {
  it('includes the outermost cells of an affine threshold plane in physical units', () => {
    const field = fieldFrom(8, (x) => (x + 2) / 4);
    const geometry = createIsosurfaceGeometry(field, 0.2, 'density');
    assertMesh(geometry);
    const positions = geometry.getAttribute('position');
    const normals = geometry.getAttribute('normal');
    for (let index = 0; index < positions.count; index += 1) {
      expect(positions.getX(index)).toBeCloseTo(-1.2, 5);
      expect(normals.getX(index)).toBeCloseTo(-1, 6);
    }
    expect(geometry.boundingBox?.min.y).toBeCloseTo(-2, 5);
    expect(geometry.boundingBox?.max.y).toBeCloseTo(2, 5);
    expect(geometry.boundingBox?.min.z).toBeCloseTo(-2, 5);
    expect(geometry.boundingBox?.max.z).toBeCloseTo(2, 5);
    geometry.dispose();
  });

  it('keeps sphere vertices on the supplied interpolated level with outward unit normals', () => {
    const n = 17;
    const field = fieldFrom(n, (x, y, z) => 1 - (x * x + y * y + z * z) / 12);
    const geometry = createIsosurfaceGeometry(field, 0.94, 'density');
    assertMesh(geometry);
    expect(closedComponentCount(geometry)).toBe(1);
    const positions = geometry.getAttribute('position');
    const normals = geometry.getAttribute('normal');
    const interpolatedSquare = (value: number): number => {
      const low = (Math.floor(((value + 2) * (n - 1)) / 4) * 4) / (n - 1) - 2;
      const high = low + 4 / (n - 1);
      return low * low + ((high * high - low * low) * (value - low)) / (high - low);
    };
    const position = new THREE.Vector3();
    const normal = new THREE.Vector3();
    for (let index = 0; index < positions.count; index += 1) {
      position.fromBufferAttribute(positions, index);
      normal.fromBufferAttribute(normals, index);
      const density =
        1 -
        (interpolatedSquare(position.x) +
          interpolatedSquare(position.y) +
          interpolatedSquare(position.z)) /
          12;
      expect(density).toBeCloseTo(0.94, 6);
      expect(position.normalize().dot(normal)).toBeGreaterThan(0.9999);
    }
    geometry.dispose();
  });

  it.each([
    ['p_z', 2, 10, 2],
    ['d_xy', 3, 16, 4],
  ] as const)('preserves closed, separate %s lobes', (orbital, n, extentBohr, count) => {
    const field = computeOrbitalField(
      { basis: 'real', n, orbital },
      { extentBohr, resolution: 32 },
    );
    const geometry = createIsosurfaceGeometry(field, 0.2, 'phase');
    assertMesh(geometry);
    expect(closedComponentCount(geometry)).toBe(count);
    geometry.dispose();
  });

  it('interpolates complex phase across the ±π seam without a false zero-phase stripe', () => {
    const field = fieldFrom(8, (x) => (x + 2) / 4);
    for (let index = 0; index < field.phaseRadians.length; index += 1)
      field.phaseRadians[index] = index % 8 < 4 ? Math.PI - 0.001 : -Math.PI + 0.001;
    const before = [
      field.densityNormalized.slice(),
      field.phaseRadians.slice(),
      field.signedAmplitudeNormalized.slice(),
    ];
    const geometry = createIsosurfaceGeometry(field, 0.5, 'phase');
    const colors = geometry.getAttribute('color');
    const target = phaseColorSrgb(Math.PI);
    const linear = new THREE.Color().setRGB(target[0], target[1], target[2], THREE.SRGBColorSpace);
    for (let index = 0; index < colors.count; index += 1) {
      expect(colors.getX(index)).toBeCloseTo(linear.r, 5);
      expect(colors.getY(index)).toBeCloseTo(linear.g, 5);
      expect(colors.getZ(index)).toBeCloseTo(linear.b, 5);
    }
    expect(field.densityNormalized).toEqual(before[0]);
    expect(field.phaseRadians).toEqual(before[1]);
    expect(field.signedAmplitudeNormalized).toEqual(before[2]);
    geometry.dispose();
  });

  it('uses a neutral color when interpolated complex amplitudes cancel', () => {
    const field = fieldFrom(8, (_x, _y, z) => (z + 2) / 4);
    // Opposite phases across x; select vertices on the exact x=0 midpoint.
    for (let index = 0; index < field.phaseRadians.length; index += 1)
      field.phaseRadians[index] = index % 8 < 4 ? 0 : Math.PI;
    const geometry = createIsosurfaceGeometry(field, 0.5, 'phase');
    const positions = geometry.getAttribute('position');
    const colors = geometry.getAttribute('color');
    const target = new THREE.Color().setRGB(...UNDEFINED_PHASE_COLOR_SRGB, THREE.SRGBColorSpace);
    let tested = 0;
    for (let index = 0; index < positions.count; index += 1)
      if (Math.abs(positions.getX(index)) < 1e-6) {
        expect(colors.getX(index)).toBeCloseTo(target.r, 5);
        expect(colors.getY(index)).toBeCloseTo(target.g, 5);
        expect(colors.getZ(index)).toBeCloseTo(target.b, 5);
        tested += 1;
      }
    expect(tested).toBeGreaterThan(0);
    geometry.dispose();
  });

  it('returns an empty compact geometry when no cell crosses the threshold and rejects invalid grids', () => {
    const field = fieldFrom(8, () => 0.2);
    const geometry = createIsosurfaceGeometry(field, 0.4, 'density');
    expect(geometry.getAttribute('position').count).toBe(0);
    geometry.dispose();
    expect(() => createIsosurfaceGeometry({ ...field, resolution: 1e6 }, 0.4, 'density')).toThrow(
      RangeError,
    );
    expect(() => createIsosurfaceGeometry(field, Number.NaN, 'density')).toThrow(RangeError);
  });
});
