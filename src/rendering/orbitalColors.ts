import { legacyOrbitalIntensity } from '../science/legacyScience';

export type Color3 = [number, number, number];

function heatFire(v: number): Color3 {
  v = Math.max(0, Math.min(1, v));
  const colors: Color3[] = [
    [0, 0, 0],
    [0.5, 0, 0.99],
    [0.8, 0, 0],
    [1, 0.5, 0],
    [1, 1, 0],
    [1, 1, 1],
  ];
  const scaledValue = v * 5,
    index = Math.min(~~scaledValue, 4),
    interpolation = scaledValue - index;
  const start = colors[index];
  const end = colors[index + 1];
  if (!start || !end) throw new RangeError('Invalid heat-map interval');
  return [
    start[0] + interpolation * (end[0] - start[0]),
    start[1] + interpolation * (end[1] - start[1]),
    start[2] + interpolation * (end[2] - start[2]),
  ];
}

function densityColor(v: number): Color3 {
  v = Math.max(0, Math.min(1, v));
  let r, g, b;
  if (v < 0.5) {
    const f = v / 0.5;
    r = 1;
    g = 1 - 0.5 * f;
    b = 0;
  } else {
    const f = (v - 0.5) / 0.5;
    r = 1 - 0.5 * f;
    g = 0.5 * (1 - f);
    b = 0.5 * f;
  }
  return [r, g, b];
}

export function legacyOrbitalColor(
  r: number,
  theta: number,
  phi: number,
  n: number,
  l: number,
  m: number,
  useAltColor: boolean,
): Color3 {
  const intensity = legacyOrbitalIntensity(r, theta, phi, n, l, m);
  const t = intensity * 1.5 * Math.pow(5, n);
  return useAltColor ? densityColor(Math.min(t, 1)) : heatFire(t);
}
