import type { Color3 } from './orbitalColors';

export interface Vector2 {
  x: number;
  y: number;
}

export interface LegacyWavePoint {
  lp: Vector2;
  dir: Vector2;
}

export interface LegacyWave {
  energy: number;
  sigma: number;
  k: number;
  phase: number;
  a: number;
  pos: Vector2;
  dir: Vector2;
  col: Color3;
  points: LegacyWavePoint[];
}

export interface LegacyAtom2D {
  pos: Vector2;
  n: number;
  angle: number;
}

export function createLegacyWave(
  energy: number,
  pos: Vector2,
  dir: Vector2,
  col: Color3 = [1, 0.6, 0.1],
): LegacyWave {
  const sigma = 30;
  const length = Math.sqrt(dir.x * dir.x + dir.y * dir.y);
  const normalizedDirection = { x: dir.x / length, y: dir.y / length };
  const points: LegacyWavePoint[] = [];
  for (let x = -sigma; x <= sigma; x += 0.5) {
    points.push({
      lp: { x: pos.x + x * normalizedDirection.x, y: pos.y + x * normalizedDirection.y },
      dir: { x: normalizedDirection.x * 120, y: normalizedDirection.y * 120 },
    });
  }
  return {
    energy,
    sigma,
    k: 0.4,
    phase: 0,
    a: 6,
    pos: { ...pos },
    dir: normalizedDirection,
    col,
    points,
  };
}

export function getLegacyWaveDisplayPoint(wave: LegacyWave, point: LegacyWavePoint): Vector2 {
  const perp = { x: -point.dir.y, y: point.dir.x };
  const perpendicularLength = Math.sqrt(perp.x * perp.x + perp.y * perp.y) || 1;
  const displacement =
    wave.a *
    Math.sin(wave.k * Math.sqrt(point.lp.x * point.lp.x + point.lp.y * point.lp.y) - wave.phase);
  return {
    x: point.lp.x + (perp.x / perpendicularLength) * displacement,
    y: point.lp.y + (perp.y / perpendicularLength) * displacement,
  };
}

export function updateLegacyWave(
  wave: LegacyWave,
  dt: number,
  width: number,
  height: number,
): boolean {
  wave.phase += 20 * dt;
  for (const point of wave.points) {
    point.lp.x += point.dir.x * dt;
    point.lp.y += point.dir.y * dt;
    if (Math.abs(point.lp.x) > width || Math.abs(point.lp.y) > height) return true;
  }
  return false;
}

export function createLegacyAtom(pos: Vector2): LegacyAtom2D {
  return {
    pos: { ...pos },
    n: 1,
    angle: Math.random() * Math.PI * 2,
  };
}

export function createLegacyAtomRing(count: number, radius: number): LegacyAtom2D[] {
  const atoms: LegacyAtom2D[] = [];
  for (let i = 0; i < count; i++) {
    const angle = (2 * Math.PI * i) / count;
    atoms.push(createLegacyAtom({ x: Math.cos(angle) * radius, y: Math.sin(angle) * radius }));
  }
  return atoms;
}

export function getLegacyElectronPosition(atom: LegacyAtom2D, orbitDistance: number): Vector2 {
  return {
    x: Math.cos(atom.angle) * atom.n * orbitDistance + atom.pos.x,
    y: Math.sin(atom.angle) * atom.n * orbitDistance + atom.pos.y,
  };
}

export function updateLegacyAtom(atom: LegacyAtom2D): void {
  atom.angle += 0.06;
}

export function createLegacyIncidentWaves(
  energy: number,
  count: number,
  startX: number,
  startY: number,
  verticalSpacing: number,
): LegacyWave[] {
  const waves: LegacyWave[] = [];
  for (let i = 0; i < count; i++) {
    waves.push(
      createLegacyWave(energy, { x: startX, y: i * verticalSpacing + startY }, { x: -1, y: 0 }),
    );
  }
  return waves;
}

export function createLegacyRadialWaves(
  energy: number,
  origin: Vector2,
  count: number,
): LegacyWave[] {
  const waves: LegacyWave[] = [];
  for (let i = 0; i < count; i++) {
    const angle = Math.random() * 2 * Math.PI;
    waves.push(createLegacyWave(energy, origin, { x: Math.cos(angle), y: Math.sin(angle) }));
  }
  return waves;
}
