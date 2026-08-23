export type Point3 = [number, number, number];

interface RadialCdf {
  cdf: Float64Array;
  rMax: number;
  M: number;
}

interface AngularCdf {
  cdf: Float64Array;
  M: number;
}

export interface SampledOrbitalPoint {
  position: Point3;
  r: number;
  theta: number;
  phi: number;
}

function requireNumber(values: ArrayLike<number>, index: number): number {
  const value = values[index];
  if (value === undefined) {
    throw new RangeError(`Indice numérique hors limites : ${index}.`);
  }
  return value;
}

function gamma(z: number): number {
  if (z < 0.5) return Math.PI / (Math.sin(Math.PI * z) * gamma(1 - z));
  z -= 1;
  const C = [
    0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
    -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
    1.5056327351493116e-7,
  ];
  let x = requireNumber(C, 0);
  for (let i = 1; i < 9; i++) x += requireNumber(C, i) / (z + i);
  const t = z + 7.5;
  return Math.sqrt(2 * Math.PI) * Math.pow(t, z + 0.5) * Math.exp(-t) * x;
}

const _rCDF: Record<string, RadialCdf | undefined> = {};
const _tCDF: Record<string, AngularCdf | undefined> = {};

function buildRCDF(n: number, l: number): RadialCdf {
  const key = `${n}_${l}`;
  const cached = _rCDF[key];
  if (cached) return cached;
  const M = 4096,
    rMax = 10 * n * n;
  const dr = rMax / (M - 1);
  const cdf = new Float64Array(M);
  let sum = 0;
  const k = n - l - 1,
    alpha = 2 * l + 1;
  for (let i = 0; i < M; i++) {
    const r = i * dr,
      rho = (2 * r) / n;
    let L = 1;
    if (k === 1) L = 1 + alpha - rho;
    else if (k > 1) {
      let Lm2 = 1,
        Lm1 = 1 + alpha - rho;
      for (let j = 2; j <= k; j++) {
        L = ((2 * j - 1 + alpha - rho) * Lm1 - (j - 1 + alpha) * Lm2) / j;
        Lm2 = Lm1;
        Lm1 = L;
      }
    }
    const norm = (Math.pow(2 / n, 3) * gamma(n - l)) / (2 * n * gamma(n + l + 1));
    const R = Math.sqrt(norm) * Math.exp(-rho / 2) * Math.pow(rho || 0, l) * L;
    const pdf = r * r * R * R;
    if (isFinite(pdf)) sum += pdf;
    cdf[i] = sum;
  }
  for (let i = 0; i < M; i++) cdf[i] = requireNumber(cdf, i) / sum;
  return (_rCDF[key] = { cdf, rMax, M });
}

function buildTCDF(l: number, absM: number): AngularCdf {
  const key = `${l}_${absM}`;
  const cached = _tCDF[key];
  if (cached) return cached;
  const M = 2048;
  const dt = Math.PI / (M - 1);
  const cdf = new Float64Array(M);
  let sum = 0;
  for (let i = 0; i < M; i++) {
    const theta = i * dt,
      x = Math.cos(theta);
    let Pmm = 1;
    if (absM > 0) {
      const s = Math.sqrt((1 - x) * (1 + x));
      let fact = 1;
      for (let j = 1; j <= absM; j++) {
        Pmm *= -fact * s;
        fact += 2;
      }
    }
    let Plm = Pmm;
    if (l > absM) {
      let Pm1m = x * (2 * absM + 1) * Pmm;
      if (l === absM + 1) {
        Plm = Pm1m;
      } else {
        let pp = Pmm;
        for (let ll = absM + 2; ll <= l; ll++) {
          const Pll = ((2 * ll - 1) * x * Pm1m - (ll + absM - 1) * pp) / (ll - absM);
          pp = Pm1m;
          Pm1m = Pll;
        }
        Plm = Pm1m;
      }
    }
    const pdf = Math.sin(theta) * Plm * Plm;
    if (isFinite(pdf)) sum += pdf;
    cdf[i] = sum;
  }
  for (let i = 0; i < M; i++) cdf[i] = requireNumber(cdf, i) / (sum || 1);
  return (_tCDF[key] = { cdf, M });
}

function bsearch(cdf: Float64Array, u: number): number {
  let lo = 0,
    hi = cdf.length - 1;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    if (requireNumber(cdf, mid) < u) lo = mid + 1;
    else hi = mid;
  }
  return lo;
}

function sampleR(n: number, l: number): number {
  const { cdf, rMax, M } = buildRCDF(n, l);
  return bsearch(cdf, Math.random()) * (rMax / (M - 1));
}

function sampleTheta(l: number, m: number): number {
  const { cdf, M } = buildTCDF(l, Math.abs(m));
  return bsearch(cdf, Math.random()) * (Math.PI / (M - 1));
}

function samplePhi(): number {
  return 2 * Math.PI * Math.random();
}

const _a0 = 52.9; // Bohr radius pm

function _rProb2p(r: number): number {
  const x = r / _a0;
  return (x ** 4 / 24) * Math.exp(-x);
}

function _rProb3p(r: number): number {
  const x = r / _a0;
  const R = x * (1 - x / 6) * Math.exp(-x / 3);
  return R * R * r * r;
}

function _sampleR2p(): number {
  const mx = 15 * _a0,
    pm = _rProb2p(4 * _a0);
  for (;;) {
    const r = Math.random() * mx;
    if (Math.random() <= _rProb2p(r) / pm) return r;
  }
}

function _sampleR3p(): number {
  const mx = 25 * _a0,
    pm = _rProb3p(8 * _a0);
  for (;;) {
    const r = Math.random() * mx;
    if (Math.random() <= _rProb3p(r) / pm) return r;
  }
}

const _rT = (): number => Math.acos(1 - 2 * Math.random()),
  _rP = (): number => 2 * Math.PI * Math.random();

export function sphericalToCartesian(r: number, th: number, ph: number): Point3 {
  const s = Math.sin(th);
  return [r * s * Math.cos(ph), r * Math.cos(th), r * s * Math.sin(ph)];
}

export const NAMED_SAMPLERS = {
  '2p_x': () => {
    const r = _sampleR2p();
    let t: number, p: number;
    do {
      t = _rT();
    } while (Math.random() > Math.sin(t) ** 3);
    do {
      p = _rP();
    } while (Math.random() > Math.cos(p) ** 2);
    return sphericalToCartesian(r, t, p);
  },
  '2p_y': () => {
    const r = _sampleR2p();
    let t: number, p: number;
    do {
      t = _rT();
    } while (Math.random() > Math.sin(t) ** 3);
    do {
      p = _rP();
    } while (Math.random() > Math.sin(p) ** 2);
    return sphericalToCartesian(r, t, p);
  },
  '2p_z': () => {
    const r = _sampleR2p();
    let t: number;
    do {
      t = _rT();
    } while (Math.random() > Math.cos(t) ** 2);
    return sphericalToCartesian(r, t, _rP());
  },
  '3p_z': () => {
    const r = _sampleR3p();
    let t: number;
    do {
      t = _rT();
    } while (Math.random() > Math.cos(t) ** 2);
    return sphericalToCartesian(r, t, _rP());
  },
} satisfies Record<string, () => Point3>;

export type NamedSamplerKey = keyof typeof NAMED_SAMPLERS;

export function sampleLegacyOrbitalPoint(
  n: number,
  l: number,
  m: number,
  namedKey: NamedSamplerKey | undefined,
): SampledOrbitalPoint {
  const namedFn = namedKey ? NAMED_SAMPLERS[namedKey] : undefined;
  if (namedFn) {
    const position = namedFn();
    const [x, y, z] = position;
    const r = Math.sqrt(x * x + y * y + z * z);
    const theta = Math.acos(Math.max(-1, Math.min(1, y / Math.max(r, 1e-9))));
    const phi = Math.atan2(z, x);
    return { position, r, theta, phi };
  }

  const r = sampleR(n, l);
  const theta = sampleTheta(l, m);
  const phi = samplePhi();
  return { position: sphericalToCartesian(r, theta, phi), r, theta, phi };
}

export function legacyOrbitalIntensity(
  r: number,
  theta: number,
  _phi: number,
  n: number,
  l: number,
  m: number,
): number {
  const rho = (2 * r) / n,
    k = n - l - 1,
    alpha = 2 * l + 1;
  let L = 1;
  if (k === 1) L = 1 + alpha - rho;
  else if (k > 1) {
    let Lm2 = 1,
      Lm1 = 1 + alpha - rho;
    for (let j = 2; j <= k; j++) {
      L = ((2 * j - 1 + alpha - rho) * Lm1 - (j - 1 + alpha) * Lm2) / j;
      Lm2 = Lm1;
      Lm1 = L;
    }
  }
  const norm = (Math.pow(2 / n, 3) * gamma(n - l)) / (2 * n * gamma(n + l + 1));
  const R = Math.sqrt(norm) * Math.exp(-rho / 2) * Math.pow(rho || 0, l) * L;
  const absM = Math.abs(m),
    x = Math.cos(theta);
  let Pmm = 1;
  if (absM > 0) {
    const s = Math.sqrt((1 - x) * (1 + x));
    let f = 1;
    for (let j = 1; j <= absM; j++) {
      Pmm *= -f * s;
      f += 2;
    }
  }
  let Plm = Pmm;
  if (l > absM) {
    let Pm1m = x * (2 * absM + 1) * Pmm;
    if (l === absM + 1) {
      Plm = Pm1m;
    } else {
      let pp = Pmm;
      for (let ll = absM + 2; ll <= l; ll++) {
        const Pll = ((2 * ll - 1) * x * Pm1m - (ll + absM - 1) * pp) / (ll - absM);
        pp = Pm1m;
        Pm1m = Pll;
      }
      Plm = Pm1m;
    }
  }
  return R * R * Plm * Plm;
}

export function legacyProbabilityFlow(px: number, py: number, pz: number, m: number): Point3 {
  const r = Math.sqrt(px * px + py * py + pz * pz);
  if (r < 1e-6) return [0, 0, 0];
  const theta = Math.acos(Math.max(-1, Math.min(1, py / r)));
  const phi = Math.atan2(pz, px);
  const st = Math.sin(theta);
  if (Math.abs(st) < 1e-4) return [0, 0, 0];
  const vm = m / (r * st);
  return [-vm * Math.sin(phi), 0, vm * Math.cos(phi)];
}

export const LEGACY_IMPORT_SCALE = 5.29; // Bohr → pm

export function convertLegacyImportedPoint([pointX, pointY, pointZ]: Point3): SampledOrbitalPoint {
  const x = pointX * LEGACY_IMPORT_SCALE,
    y = pointY * LEGACY_IMPORT_SCALE,
    z = pointZ * LEGACY_IMPORT_SCALE;
  const r = Math.sqrt(x * x + y * y + z * z);
  const theta = Math.acos(Math.max(-1, Math.min(1, y / Math.max(r, 1e-6))));
  const phi = Math.atan2(z, x);
  return { position: [x, y, z], r, theta, phi };
}

export const LEGACY_FLOW_DT = 0.4;

export function advanceLegacyProbabilityFlowPoint(
  px: number,
  py: number,
  pz: number,
  targetRadius: number,
  m: number,
): Point3 {
  const [vx, , vz] = legacyProbabilityFlow(px, py, pz, m);
  const nx = px + vx * LEGACY_FLOW_DT,
    nz = pz + vz * LEGACY_FLOW_DT;
  const nr = Math.sqrt(nx * nx + py * py + nz * nz);
  if (nr > 0.01) {
    const scale = targetRadius / nr;
    return [nx * scale, py, nz * scale];
  }
  return [px, py, pz];
}

export function legacyEnergyEv(n: number): number {
  return -13.6 / (n * n);
}

export const LEGACY_PHOTOELECTRIC_ENERGY = -13.6 / 4 - -13.6;
