import * as THREE from 'three';
import { MarchingCubes } from 'three/examples/jsm/objects/MarchingCubes.js';
import type { OrbitalFieldGrid } from '../sampling/fieldContracts';
import { DENSITY_COLOR_SRGB, phaseColorSrgb } from './phasePalette';
import type { OrbitalObservable } from './renderingContracts';

/** Numerical tessellation refinement only; this never resamples the scientific engine. */
const REFINEMENT = 2;

function sample(values: Float32Array, n: number, x: number, y: number, z: number): number {
  const cx = Math.max(0, Math.min(n - 1, x));
  const cy = Math.max(0, Math.min(n - 1, y));
  const cz = Math.max(0, Math.min(n - 1, z));
  const ix = Math.min(n - 2, Math.floor(cx));
  const iy = Math.min(n - 2, Math.floor(cy));
  const iz = Math.min(n - 2, Math.floor(cz));
  const tx = cx - ix;
  const ty = cy - iy;
  const tz = cz - iz;
  let value = 0;
  for (let dz = 0; dz <= 1; dz += 1) {
    for (let dy = 0; dy <= 1; dy += 1) {
      for (let dx = 0; dx <= 1; dx += 1) {
        const weight = (dx ? tx : 1 - tx) * (dy ? ty : 1 - ty) * (dz ? tz : 1 - tz);
        if (weight > 0) value += weight * (values[ix + dx + n * (iy + dy + n * (iz + dz))] ?? 0);
      }
    }
  }
  return value;
}

function gradientBuffers(
  field: OrbitalFieldGrid,
): readonly [Float32Array, Float32Array, Float32Array] {
  const n = field.resolution;
  const result = [
    new Float32Array(n ** 3),
    new Float32Array(n ** 3),
    new Float32Array(n ** 3),
  ] as const;
  const spacing = (2 * field.extentBohr) / (n - 1);
  for (let z = 0; z < n; z += 1) {
    for (let y = 0; y < n; y += 1) {
      for (let x = 0; x < n; x += 1) {
        const index = x + n * (y + n * z);
        for (const [axis, coordinate, stride] of [
          [0, x, 1],
          [1, y, n],
          [2, z, n * n],
        ] as const) {
          const low = coordinate > 0 ? index - stride : index;
          const high = coordinate < n - 1 ? index + stride : index;
          const intervals = coordinate > 0 && coordinate < n - 1 ? 2 : 1;
          // Outward from the high-density region: minus the density gradient.
          result[axis][index] =
            ((field.densityNormalized[low] ?? 0) - (field.densityNormalized[high] ?? 0)) /
            (intervals * spacing);
        }
      }
    }
  }
  return result;
}

function phaseComponents(field: OrbitalFieldGrid): readonly [Float32Array, Float32Array] {
  const real = new Float32Array(field.densityNormalized.length);
  const imaginary = new Float32Array(real.length);
  for (let index = 0; index < real.length; index += 1) {
    const phase = field.phaseRadians[index] ?? Number.NaN;
    if (!Number.isFinite(phase)) continue;
    const amplitude = Math.sqrt(field.densityNormalized[index] ?? 0);
    real[index] = amplitude * Math.cos(phase);
    imaginary[index] = amplitude * Math.sin(phase);
  }
  return [real, imaginary];
}

/**
 * Extracts the same normalized density threshold in physical a₀ coordinates.
 * Trilinear subdivision improves tessellation, not the resolution of the field.
 * A clamped border includes every original boundary cell in Three's extraction
 * loop (indices 1 through size-3). Source buffers and their normalization stay intact.
 */
export function createIsosurfaceGeometry(
  field: OrbitalFieldGrid,
  threshold: number,
  observable: OrbitalObservable,
): THREE.BufferGeometry {
  const n = field.resolution;
  if (
    !Number.isSafeInteger(n) ||
    n < 2 ||
    n > 64 ||
    !Number.isFinite(field.extentBohr) ||
    field.extentBohr <= 0
  ) {
    throw new RangeError(
      'The rendering grid requires 2–64 samples per axis and a positive finite extent.',
    );
  }
  if (!Number.isFinite(threshold) || threshold <= 0 || threshold >= 1)
    throw new RangeError('The density threshold must be in (0,1).');
  if (
    [field.densityNormalized, field.phaseRadians, field.signedAmplitudeNormalized].some(
      (values) => values.length !== n ** 3,
    )
  )
    throw new RangeError('Inconsistent rendering field buffer lengths.');
  if (
    field.densityNormalized.some((value) => !Number.isFinite(value) || value < 0 || value > 1) ||
    field.signedAmplitudeNormalized.some((value) => !Number.isFinite(value))
  )
    throw new RangeError('The rendering field contains invalid values.');
  const size = REFINEMENT * (n - 1) + 3;
  const density = new Float32Array(size ** 3);
  for (let z = 0; z < size; z += 1) {
    for (let y = 0; y < size; y += 1) {
      for (let x = 0; x < size; x += 1) {
        density[x + size * (y + size * z)] = sample(
          field.densityNormalized,
          n,
          (x - 1) / REFINEMENT,
          (y - 1) / REFINEMENT,
          (z - 1) / REFINEMENT,
        );
      }
    }
  }
  let activeCells = 0;
  const offsets = [
    0,
    1,
    size,
    size + 1,
    size * size,
    size * size + 1,
    size * size + size,
    size * size + size + 1,
  ];
  for (let z = 1; z < size - 2; z += 1) {
    for (let y = 1; y < size - 2; y += 1) {
      for (let x = 1; x < size - 2; x += 1) {
        const index = x + size * (y + size * z);
        let below = false;
        let above = false;
        for (const offset of offsets) {
          if ((density[index + offset] ?? 0) < threshold) below = true;
          else above = true;
        }
        if (below && above) activeCells += 1;
      }
    }
  }
  const geometry = new THREE.BufferGeometry();
  if (!activeCells) {
    for (const name of ['position', 'normal', 'color'])
      geometry.setAttribute(name, new THREE.BufferAttribute(new Float32Array(), 3));
    return geometry;
  }
  const material = new THREE.MeshBasicMaterial();
  const extractor = new MarchingCubes(size, material, false, false, 5 * activeCells);
  try {
    extractor.field.set(density);
    extractor.isolation = threshold;
    extractor.update();
    const gradients = gradientBuffers(field);
    const phases = observable === 'phase' && !field.nodesAvailable ? phaseComponents(field) : null;
    const positions = new Float32Array(extractor.count * 3);
    const normals = new Float32Array(positions.length);
    const colors = new Float32Array(positions.length);
    const color = new THREE.Color();
    const a = new THREE.Vector3();
    const b = new THREE.Vector3();
    const c = new THREE.Vector3();
    const face = new THREE.Vector3();
    const edge = new THREE.Vector3();
    const normal = new THREE.Vector3();
    let output = 0;
    const originalIndex = (value: number): number => (((value + 1) * size) / 2 - 1) / REFINEMENT;
    const physical = (value: number): number =>
      field.extentBohr * ((2 * originalIndex(value)) / (n - 1) - 1);
    const readVertex = (offset: number, target: THREE.Vector3): THREE.Vector3 =>
      target.set(
        physical(extractor.positionArray[offset] ?? Number.NaN),
        physical(extractor.positionArray[offset + 1] ?? Number.NaN),
        physical(extractor.positionArray[offset + 2] ?? Number.NaN),
      );
    for (let offset = 0; offset < extractor.count * 3; offset += 9) {
      readVertex(offset, a);
      readVertex(offset + 3, b);
      readVertex(offset + 6, c);
      face.subVectors(b, a).cross(edge.subVectors(c, a));
      if (
        !Number.isFinite(face.lengthSq()) ||
        face.lengthSq() <= Number.EPSILON ** 2 * field.extentBohr ** 4
      )
        continue;
      face.normalize();
      const cx = (((a.x + b.x + c.x) / (3 * field.extentBohr) + 1) * (n - 1)) / 2;
      const cy = (((a.y + b.y + c.y) / (3 * field.extentBohr) + 1) * (n - 1)) / 2;
      const cz = (((a.z + b.z + c.z) / (3 * field.extentBohr) + 1) * (n - 1)) / 2;
      normal.set(
        sample(gradients[0], n, cx, cy, cz),
        sample(gradients[1], n, cx, cy, cz),
        sample(gradients[2], n, cx, cy, cz),
      );
      const reversed = face.dot(normal) < 0;
      if (reversed) face.negate();
      for (const vertex of reversed ? [a, c, b] : [a, b, c]) {
        const x = ((vertex.x / field.extentBohr + 1) * (n - 1)) / 2;
        const y = ((vertex.y / field.extentBohr + 1) * (n - 1)) / 2;
        const z = ((vertex.z / field.extentBohr + 1) * (n - 1)) / 2;
        normal.set(
          sample(gradients[0], n, x, y, z),
          sample(gradients[1], n, x, y, z),
          sample(gradients[2], n, x, y, z),
        );
        if (normal.lengthSq() === 0) normal.copy(face);
        else normal.normalize();
        let phase = Number.NaN;
        if (observable === 'phase') {
          if (field.nodesAvailable) {
            const signed = sample(field.signedAmplitudeNormalized, n, x, y, z);
            phase = signed === 0 ? Number.NaN : signed > 0 ? 0 : Math.PI;
          } else if (phases) {
            const real = sample(phases[0], n, x, y, z);
            const imaginary = sample(phases[1], n, x, y, z);
            // Relative roundoff at exact cancellation is not a defined phase.
            // Display interpolation of supplied normalized amplitudes only, not a physics evaluation.
            const amplitude = Math.sqrt(Math.max(0, sample(field.densityNormalized, n, x, y, z)));
            if (Math.hypot(real, imaginary) > 8 * 2 ** -23 * amplitude)
              phase = Math.atan2(imaginary, real);
          }
        }
        const srgb = observable === 'phase' ? phaseColorSrgb(phase) : DENSITY_COLOR_SRGB;
        color.setRGB(srgb[0], srgb[1], srgb[2], THREE.SRGBColorSpace);
        vertex.toArray(positions, output);
        normal.toArray(normals, output);
        color.toArray(colors, output);
        output += 3;
      }
    }
    geometry.setAttribute('position', new THREE.BufferAttribute(positions.slice(0, output), 3));
    geometry.setAttribute('normal', new THREE.BufferAttribute(normals.slice(0, output), 3));
    geometry.setAttribute('color', new THREE.BufferAttribute(colors.slice(0, output), 3));
    geometry.computeBoundingBox();
    geometry.computeBoundingSphere();
    return geometry;
  } finally {
    extractor.geometry.dispose();
    material.dispose();
  }
}
