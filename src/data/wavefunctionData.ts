export type Point3 = [number, number, number];

export interface WavefunctionData {
  points: Point3[];
}

function isPoint3(value: unknown): value is Point3 {
  return (
    Array.isArray(value) &&
    value.length === 3 &&
    value.every((coordinate) => typeof coordinate === 'number' && Number.isFinite(coordinate))
  );
}

export function parseWavefunctionData(value: unknown): WavefunctionData {
  if (typeof value !== 'object' || value === null || !('points' in value)) {
    throw new Error('Missing "points" array');
  }
  const { points } = value;
  if (!Array.isArray(points) || !points.every(isPoint3)) {
    throw new Error('Invalid "points" array');
  }
  return { points };
}
