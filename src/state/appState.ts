export interface QuantumState {
  n: number;
  l: number;
  m: number;
  N: number;
}

export interface CloudBuffers {
  positions: Float32Array;
  radii: Float32Array;
}

export interface AppState {
  quantum: QuantumState;
  cloud: CloudBuffers;
  animating: boolean;
  pointSize: number;
  useAltColor: boolean;
  show2D: boolean;
  simulation2DInitialized: boolean;
}

export function createAppState(): AppState {
  return {
    quantum: { n: 2, l: 1, m: 0, N: 15000 },
    cloud: {
      positions: new Float32Array(0),
      radii: new Float32Array(0),
    },
    animating: true,
    pointSize: 0.24,
    useAltColor: false,
    show2D: false,
    simulation2DInitialized: false,
  };
}

export function clampValue(value: number, min: number, max: number): number {
  return Math.max(min, Math.min(max, value));
}

export function applyAtomValues(
  quantum: QuantumState,
  nValue: number,
  lValue: number,
  mValue: number,
): void {
  quantum.n = Math.max(1, nValue || 1);
  quantum.l = Math.max(0, Math.min(quantum.n - 1, lValue || 0));
  quantum.m = clampValue(mValue || 0, -quantum.l, quantum.l);
}

export function normalizeQuantumState(quantum: QuantumState): void {
  if (quantum.l > quantum.n - 1) quantum.l = quantum.n - 1;
  if (quantum.l < 0) quantum.l = 0;
  if (quantum.m > quantum.l) quantum.m = quantum.l;
  if (quantum.m < -quantum.l) quantum.m = -quantum.l;
}
