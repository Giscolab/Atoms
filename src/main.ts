import * as THREE from 'three';

import { chooseTheme, inclusiveIntegerRange, parseOrbitalOption, type Theme } from './uiState';

// Phase 1 preserves the legacy r128 color pipeline for visual equivalence.
THREE.ColorManagement.enabled = false;

type Point3 = [number, number, number];
type Color3 = Point3;

interface RadialCdf {
  cdf: Float64Array;
  rMax: number;
  M: number;
}

interface AngularCdf {
  cdf: Float64Array;
  M: number;
}

interface QuantumState {
  n: number;
  l: number;
  m: number;
  N: number;
}

type ElementConstructor<T extends HTMLElement> = new () => T;

function requireElement<T extends HTMLElement>(id: string, constructor: ElementConstructor<T>): T {
  const element = document.getElementById(id);
  if (!(element instanceof constructor)) {
    throw new Error(`Élément #${id} introuvable ou de type inattendu.`);
  }
  return element;
}

function requireCanvas2dContext(canvas: HTMLCanvasElement): CanvasRenderingContext2D {
  const context = canvas.getContext('2d');
  if (!context) throw new Error('Contexte Canvas 2D indisponible.');
  return context;
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

const NAMED_SAMPLERS = {
  '2p_x': () => {
    const r = _sampleR2p();
    let t: number, p: number;
    do {
      t = _rT();
    } while (Math.random() > Math.sin(t) ** 3);
    do {
      p = _rP();
    } while (Math.random() > Math.cos(p) ** 2);
    return s2c(r, t, p);
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
    return s2c(r, t, p);
  },
  '2p_z': () => {
    const r = _sampleR2p();
    let t: number;
    do {
      t = _rT();
    } while (Math.random() > Math.cos(t) ** 2);
    return s2c(r, t, _rP());
  },
  '3p_z': () => {
    const r = _sampleR3p();
    let t: number;
    do {
      t = _rT();
    } while (Math.random() > Math.cos(t) ** 2);
    return s2c(r, t, _rP());
  },
} satisfies Record<string, () => Point3>;

type NamedSamplerKey = keyof typeof NAMED_SAMPLERS;

const SAMPLER_MAP: Partial<Record<string, NamedSamplerKey>> = {
  '2_1_0': '2p_z',
  '2_1_1': '2p_x',
  '2_1_-1': '2p_y',
  '3_1_0': '3p_z',
};

interface WavefunctionData {
  points: Point3[];
}

function isPoint3(value: unknown): value is Point3 {
  return (
    Array.isArray(value) &&
    value.length === 3 &&
    value.every((coordinate) => typeof coordinate === 'number' && Number.isFinite(coordinate))
  );
}

function parseWavefunctionData(value: unknown): WavefunctionData {
  if (typeof value !== 'object' || value === null || !('points' in value)) {
    throw new Error('Missing "points" array');
  }
  const { points } = value;
  if (!Array.isArray(points) || !points.every(isPoint3)) {
    throw new Error('Invalid "points" array');
  }
  return { points };
}

async function loadWavefunction(file: File): Promise<void> {
  return new Promise<void>((resolve, reject) => {
    const reader = new FileReader();
    reader.onload = () => {
      try {
        if (typeof reader.result !== 'string') {
          throw new Error('Invalid JSON file contents');
        }
        const data = parseWavefunctionData(JSON.parse(reader.result) as unknown);
        const BPM = 5.29; // Bohr → pm
        const N = data.points.length;
        posArr = new Float32Array(N * 3);
        rArr = new Float32Array(N);
        const colArr = new Float32Array(N * 3);
        data.points.forEach(([pointX, pointY, pointZ], i) => {
          const x = pointX * BPM,
            y = pointY * BPM,
            z = pointZ * BPM;
          posArr[i * 3] = x;
          posArr[i * 3 + 1] = y;
          posArr[i * 3 + 2] = z;
          const radius = Math.sqrt(x * x + y * y + z * z);
          rArr[i] = radius;
          const th = Math.acos(Math.max(-1, Math.min(1, y / Math.max(radius, 1e-6))));
          const ph = Math.atan2(z, x);
          const [cr, cg, cb] = orbitalColor(radius, th, ph, qn.n, qn.l, qn.m);
          colArr[i * 3] = cr;
          colArr[i * 3 + 1] = cg;
          colArr[i * 3 + 2] = cb;
        });
        if (cloudMesh) {
          scene.remove(cloudMesh);
          cloudMesh.geometry.dispose();
          cloudMesh.material.dispose();
          cloudMesh = null;
        }
        const geo = new THREE.BufferGeometry();
        geo.setAttribute('position', new THREE.BufferAttribute(posArr.slice(), 3));
        geo.setAttribute('color', new THREE.BufferAttribute(colArr, 3));
        cloudMesh = new THREE.Points(
          geo,
          new THREE.PointsMaterial({
            size: ptSizeVal,
            vertexColors: true,
            transparent: true,
            opacity: 0.88,
            sizeAttenuation: true,
            depthWrite: false,
          }),
        );
        scene.add(cloudMesh);
        qn.N = N;
        requireElement('nPartsVal', HTMLOutputElement).textContent = N.toLocaleString();
        updateHUD();
        hideLoading();
        resolve();
      } catch (error) {
        reject(error instanceof Error ? error : new Error(String(error)));
      }
    };
    reader.onerror = () => {
      reject(reader.error ?? new Error('Unable to read JSON file'));
    };
    reader.readAsText(file);
  });
}

function s2c(r: number, th: number, ph: number): Point3 {
  const s = Math.sin(th);
  return [r * s * Math.cos(ph), r * Math.cos(th), r * s * Math.sin(ph)];
}

function heatFire(v: number): Color3 {
  v = Math.max(0, Math.min(1, v));
  const s: Color3[] = [
    [0, 0, 0],
    [0.5, 0, 0.99],
    [0.8, 0, 0],
    [1, 0.5, 0],
    [1, 1, 0],
    [1, 1, 1],
  ];
  const sv = v * 5,
    i = Math.min(~~sv, 4),
    t = sv - i;
  const start = s[i];
  const end = s[i + 1];
  if (!start || !end) throw new RangeError('Invalid heat-map interval');
  return [
    start[0] + t * (end[0] - start[0]),
    start[1] + t * (end[1] - start[1]),
    start[2] + t * (end[2] - start[2]),
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
let useAltColor = false;

function orbitalColor(
  r: number,
  theta: number,
  _phi: number,
  n: number,
  l: number,
  m: number,
): Color3 {
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
  const intensity = R * R * Plm * Plm;
  const t = intensity * 1.5 * Math.pow(5, n);
  return useAltColor ? densityColor(Math.min(t, 1)) : heatFire(t);
}

function probFlow(px: number, py: number, pz: number, m: number): Point3 {
  const r = Math.sqrt(px * px + py * py + pz * pz);
  if (r < 1e-6) return [0, 0, 0];
  const theta = Math.acos(Math.max(-1, Math.min(1, py / r)));
  const phi = Math.atan2(pz, px);
  const st = Math.sin(theta);
  if (Math.abs(st) < 1e-4) return [0, 0, 0];
  const vm = m / (r * st);
  return [-vm * Math.sin(phi), 0, vm * Math.cos(phi)];
}

const canvas3d = requireElement('atomSimCanvas', HTMLCanvasElement);

const renderer = new THREE.WebGLRenderer({ canvas: canvas3d, antialias: true, alpha: false });
renderer.outputColorSpace = THREE.LinearSRGBColorSpace;
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 1.5));
renderer.setClearColor(0x050302, 1);

const scene = new THREE.Scene();
const cam = new THREE.PerspectiveCamera(60, 800 / 600, 0.01, 2000);

function resizeCanvas(): void {
  const vp = requireElement('viewport', HTMLElement);
  const w = vp.clientWidth - 40;
  const h = vp.clientHeight - 40;
  canvas3d.width = w;
  canvas3d.height = h;
  cam.aspect = w / h;
  cam.updateProjectionMatrix();
  renderer.setSize(w, h);
}

const orb = { az: 0.3, el: Math.PI / 2.5, r: 60 };
let dragging = false,
  lastX = 0,
  lastY = 0;

function updateCamera(): void {
  const el = Math.max(0.02, Math.min(Math.PI - 0.02, orb.el));
  cam.position.set(
    orb.r * Math.sin(el) * Math.cos(orb.az),
    orb.r * Math.cos(el),
    orb.r * Math.sin(el) * Math.sin(orb.az),
  );
  cam.lookAt(0, 0, 0);
}
updateCamera();

canvas3d.addEventListener('mousedown', (e) => {
  dragging = true;
  lastX = e.clientX;
  lastY = e.clientY;
});
window.addEventListener('mouseup', () => (dragging = false));
window.addEventListener('mousemove', (e) => {
  if (!dragging) return;
  orb.az += (e.clientX - lastX) * 0.008;
  orb.el -= (e.clientY - lastY) * 0.008;
  lastX = e.clientX;
  lastY = e.clientY;
  updateCamera();
  requireElement('iDist', HTMLElement).textContent = orb.r.toFixed(1);
});
canvas3d.addEventListener(
  'wheel',
  (e) => {
    orb.r = Math.max(5, Math.min(400, orb.r + e.deltaY * 0.05));
    updateCamera();
    requireElement('iDist', HTMLElement).textContent = orb.r.toFixed(1);
    e.preventDefault();
  },
  { passive: false },
);

let tLast: Touch | null = null;
canvas3d.addEventListener('touchstart', (e) => {
  const touch = e.touches[0];
  if (touch) tLast = touch;
});
canvas3d.addEventListener(
  'touchmove',
  (e) => {
    if (!tLast) return;
    const t = e.touches[0];
    if (!t) return;
    orb.az += (t.clientX - tLast.clientX) * 0.01;
    orb.el -= (t.clientY - tLast.clientY) * 0.01;
    tLast = t;
    updateCamera();
    e.preventDefault();
  },
  { passive: false },
);

(function () {
  const N = 1200,
    pos = new Float32Array(N * 3),
    col = new Float32Array(N * 3);
  for (let i = 0; i < N; i++) {
    const th = Math.acos(2 * Math.random() - 1),
      ph = 2 * Math.PI * Math.random(),
      r = 400 + Math.random() * 200;
    pos[i * 3] = r * Math.sin(th) * Math.cos(ph);
    pos[i * 3 + 1] = r * Math.cos(th);
    pos[i * 3 + 2] = r * Math.sin(th) * Math.sin(ph);
    const b = 0.08 + Math.random() * 0.25;
    col[i * 3] = b * 0.9;
    col[i * 3 + 1] = b * 0.7;
    col[i * 3 + 2] = b * 0.3;
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  g.setAttribute('color', new THREE.BufferAttribute(col, 3));
  scene.add(
    new THREE.Points(
      g,
      new THREE.PointsMaterial({ size: 0.4, vertexColors: true, sizeAttenuation: true }),
    ),
  );
})();

const nucleus = new THREE.Mesh(
  new THREE.SphereGeometry(0.5, 16, 16),
  new THREE.MeshBasicMaterial({ color: 0xff6030 }),
);
scene.add(nucleus);
const ring = new THREE.Mesh(
  new THREE.RingGeometry(0.7, 1.0, 32),
  new THREE.MeshBasicMaterial({
    color: 0xff4020,
    side: THREE.DoubleSide,
    transparent: true,
    opacity: 0.3,
  }),
);
scene.add(ring);

type CloudMesh = THREE.Points<THREE.BufferGeometry, THREE.PointsMaterial>;

let cloudMesh: CloudMesh | null = null;
let posArr = new Float32Array(0),
  rArr = new Float32Array(0);
const qn: QuantumState = { n: 2, l: 1, m: 0, N: 15000 };
let animating = true;
let ptSizeVal = 0.24;
const THEME_STORAGE_KEY = 'atom-theme';
const THEME_CLEAR_COLOR: Record<Theme, number> = { dark: 0x050302, light: 0xf3ece4 };

function getPreferredTheme(): Theme {
  let storedTheme: string | null = null;
  try {
    storedTheme = localStorage.getItem(THEME_STORAGE_KEY);
  } catch {}
  return chooseTheme(storedTheme, window.matchMedia('(prefers-color-scheme: light)').matches);
}

function applyTheme(theme: Theme, persist = true): void {
  const normalizedTheme: Theme = theme === 'light' ? 'light' : 'dark';
  document.documentElement.setAttribute('data-theme', normalizedTheme);
  const themeToggle = requireElement('themeToggle', HTMLInputElement);
  themeToggle.checked = normalizedTheme === 'light';
  renderer.setClearColor(THEME_CLEAR_COLOR[normalizedTheme], 1);
  if (persist) {
    try {
      localStorage.setItem(THEME_STORAGE_KEY, normalizedTheme);
    } catch {}
  }
}

const ORBITAL_LETTERS: string[] = ['s', 'p', 'd', 'f', 'g', 'h', 'i'];
const ORBITAL_FAMILY: Record<string, string> = {
  s: 'Distribution sphérique autour du noyau.',
  p: 'Deux lobes principaux avec un plan nodal.',
  d: 'Forme complexe a plusieurs lobes (type trèfle).',
  f: 'Forme multi-lobes avec structure plus fine.',
  g: 'Orbitale de haut ordre avec nombreux noeuds angulaires.',
  h: 'Orbitale de haut ordre, très structurée.',
  i: 'Orbitale de tres haut ordre et forte complexité angulaire.',
};
const PANEL_MIN_N = 1;
const PANEL_MAX_N = 9;

function clampValue(v: number, min: number, max: number): number {
  return Math.max(min, Math.min(max, v));
}

function orbitalLetter(l: number): string {
  return ORBITAL_LETTERS[l] || `l=${l}`;
}

function orbitalFamilySummary(l: number): string {
  const letter = orbitalLetter(l);
  return ORBITAL_FAMILY[letter] || `Famille orbitale definie par l=${l}.`;
}

function magneticProjectionSummary(m: number): string {
  if (m === 0) return 'Symetrie axiale (m=0).';
  if (m > 0) return `Projection positive (m=+${m}).`;
  return `Projection negative (m=${m}).`;
}

function optionText(option: HTMLOptionElement | undefined): string {
  if (!option) return '';
  return option.textContent.replace(/\s+/g, ' ').trim();
}

function setSelectOptions(
  select: HTMLSelectElement,
  values: number[],
  selectedValue: number,
): void {
  select.innerHTML = '';
  for (const value of values) {
    const option = document.createElement('option');
    option.value = String(value);
    option.textContent = String(value);
    select.appendChild(option);
  }
  select.value = String(selectedValue);
}

function refreshAtomForm(fromState = true): void {
  const nSel = requireElement('atomN', HTMLSelectElement);
  const lSel = requireElement('atomL', HTMLSelectElement);
  const mSel = requireElement('atomM', HTMLSelectElement);

  const requestedN = fromState ? qn.n : Number(nSel.value || qn.n);
  const maxN = Math.max(PANEL_MAX_N, qn.n);
  const nVal = clampValue(requestedN, PANEL_MIN_N, maxN);
  setSelectOptions(nSel, inclusiveIntegerRange(PANEL_MIN_N, maxN), nVal);

  const requestedL = fromState ? qn.l : Number(lSel.value || qn.l);
  const lMax = Math.max(0, nVal - 1);
  const lVal = clampValue(requestedL, 0, lMax);
  setSelectOptions(lSel, inclusiveIntegerRange(0, lMax), lVal);

  const requestedM = fromState ? qn.m : Number(mSel.value || qn.m);
  const mVal = clampValue(requestedM, -lVal, lVal);
  setSelectOptions(mSel, inclusiveIntegerRange(-lVal, lVal), mVal);
}

function updatePresetInfo(): void {
  const { n, l, m } = qn;
  const letter = orbitalLetter(l);
  const radialNodes = Math.max(0, n - l - 1);
  const angularNodes = l;

  const sel = requireElement('orbitalSel', HTMLSelectElement);
  const matchValue = `${n}_${l}_${m}`;
  const selected = [...sel.options].find((o) => o.value === matchValue);
  const label = selected ? optionText(selected) : `${n}${letter} (m=${m})`;
  const sourceTag = selected ? 'Prereglage du catalogue' : "Etat defini via l'onglet Atome";

  requireElement('presetName', HTMLElement).textContent = label;
  requireElement('presetSummary', HTMLElement).textContent =
    `${orbitalFamilySummary(l)} ${sourceTag}.`;
  requireElement('presetFamily', HTMLElement).textContent = letter;
  requireElement('presetNumbers', HTMLElement).textContent = `n=${n}, l=${l}, m=${m}`;
  requireElement('presetNodes', HTMLElement).textContent =
    `${radialNodes} radial, ${angularNodes} angulaire`;
  requireElement('presetProjection', HTMLElement).textContent = magneticProjectionSummary(m);
  requireElement('presetLock', HTMLElement).textContent =
    'Contraintes fixes: n>=1, 0<=l<=n-1, -l<=m<=l.';
}

function setPanelTab(tabId: 'presets' | 'atom'): void {
  const showPresets = tabId === 'presets';

  const presets = requireElement('panel-view-presets', HTMLElement);
  const atom = requireElement('panel-view-atom', HTMLElement);

  presets.classList.toggle('active', showPresets);
  atom.classList.toggle('active', !showPresets);

  presets.hidden = !showPresets;
  atom.hidden = showPresets;

  requireElement('tabBtnPresets', HTMLButtonElement).setAttribute(
    'aria-selected',
    showPresets ? 'true' : 'false',
  );
  requireElement('tabBtnAtom', HTMLButtonElement).setAttribute(
    'aria-selected',
    showPresets ? 'false' : 'true',
  );
}

function applyAtomForm(): void {
  const nVal = Number(requireElement('atomN', HTMLSelectElement).value);
  const lVal = Number(requireElement('atomL', HTMLSelectElement).value);
  const mVal = Number(requireElement('atomM', HTMLSelectElement).value);

  qn.n = Math.max(1, nVal || 1);
  qn.l = Math.max(0, Math.min(qn.n - 1, lVal || 0));
  qn.m = clampValue(mVal || 0, -qn.l, qn.l);

  updateHUD();
  generateCloud();
}

function showLoading(msg?: string): void {
  const ld = requireElement('loading', HTMLElement);
  ld.style.display = 'block';
  ld.textContent = msg || "Calcul de l'orbitale...";
}
function hideLoading(): void {
  requireElement('loading', HTMLElement).style.display = 'none';
  requireElement('progressBar', HTMLElement).textContent = '';
}

const _dotTex = (() => {
  const c = document.createElement('canvas');
  c.width = c.height = 16;
  const x = requireCanvas2dContext(c);
  const g = x.createRadialGradient(8, 8, 0, 8, 8, 8);
  g.addColorStop(0, 'rgba(255,255,255,1)');
  g.addColorStop(0.5, 'rgba(255,255,255,0.8)');
  g.addColorStop(1, 'rgba(255,255,255,0)');
  x.fillStyle = g;
  x.fillRect(0, 0, 16, 16);
  return new THREE.CanvasTexture(c);
})();

function generateCloud(callback?: () => void): void {
  const { n, l, m, N } = qn;
  showLoading("Calcul de l'orbitale...");

  setTimeout(() => {
    posArr = new Float32Array(N * 3);
    rArr = new Float32Array(N);
    const colArr = new Float32Array(N * 3);
    let done = 0;
    const chunk = 2000;

    function doChunk(): void {
      const end = Math.min(done + chunk, N);
      for (let i = done; i < end; i++) {
        let position: Point3;
        let r: number, theta: number, phi: number;
        const namedKey = SAMPLER_MAP[`${n}_${l}_${m}`];
        const namedFn = namedKey ? NAMED_SAMPLERS[namedKey] : undefined;
        if (namedFn) {
          position = namedFn();
          const [x, y, z] = position;
          r = Math.sqrt(x * x + y * y + z * z);
          theta = Math.acos(Math.max(-1, Math.min(1, y / Math.max(r, 1e-9))));
          phi = Math.atan2(z, x);
        } else {
          r = sampleR(n, l);
          theta = sampleTheta(l, m);
          phi = samplePhi();
          position = s2c(r, theta, phi);
        }
        const [x, y, z] = position;
        posArr[i * 3] = x;
        posArr[i * 3 + 1] = y;
        posArr[i * 3 + 2] = z;
        rArr[i] = r;
        const [cr, cg, cb] = orbitalColor(r, theta, phi, n, l, m);
        colArr[i * 3] = cr;
        colArr[i * 3 + 1] = cg;
        colArr[i * 3 + 2] = cb;
      }
      done = end;
      requireElement('progressBar', HTMLElement).textContent =
        `${Math.round((done / N) * 100)}% - ${done.toLocaleString()} / ${N.toLocaleString()}`;
      if (done < N) {
        setTimeout(doChunk, 0);
        return;
      }

      // Build Three.js geometry
      setTimeout(() => {
        if (cloudMesh) {
          scene.remove(cloudMesh);
          cloudMesh.geometry.dispose();
          cloudMesh.material.dispose();
          cloudMesh = null;
        }
        const geo = new THREE.BufferGeometry();
        geo.setAttribute('position', new THREE.BufferAttribute(posArr.slice(), 3));
        geo.setAttribute('color', new THREE.BufferAttribute(colArr, 3));
        cloudMesh = new THREE.Points(
          geo,
          new THREE.PointsMaterial({
            size: ptSizeVal,
            vertexColors: true,
            transparent: true,
            opacity: 0.88,
            sizeAttenuation: true,
            depthWrite: false,
            map: _dotTex,
            alphaTest: 0.01,
          }),
        );
        scene.add(cloudMesh);
        hideLoading();
        updateHUD();
        if (callback) callback();
      }, 0);
    }
    doChunk();
  }, 50);
}

function updateHUD(): void {
  const { n, l, m, N } = qn;
  requireElement('qn-n', HTMLElement).textContent = String(n);
  requireElement('qn-l', HTMLElement).textContent = String(l);
  requireElement('qn-m', HTMLElement).textContent = String(m);
  requireElement('energy-val', HTMLElement).textContent = `${(-13.6 / (n * n)).toFixed(3)} eV`;
  requireElement('iOrb', HTMLElement).textContent = `${n}, ${l}, ${m}`;
  requireElement('iPts', HTMLElement).textContent = N.toLocaleString();
  requireElement('iDist', HTMLElement).textContent = orb.r.toFixed(1);
  const sel = requireElement('orbitalSel', HTMLSelectElement);
  const match = `${n}_${l}_${m}`;
  const hasMatch = [...sel.options].some((o) => o.value === match);
  const customOption = sel.querySelector('option[data-custom-state="true"]');
  if (customOption) customOption.remove();
  if (hasMatch) {
    sel.value = match;
  } else {
    const custom = document.createElement('option');
    custom.value = match;
    custom.textContent = `Etat perso (${n}, ${l}, ${m})`;
    custom.dataset.customState = 'true';
    sel.prepend(custom);
    sel.value = match;
  }
  refreshAtomForm(true);
  updatePresetInfo();
}

function applyOrbitalOption(value: string): void {
  const parsed = parseOrbitalOption(value);
  if (!parsed) throw new Error(`Option orbitale invalide : ${value}`);
  [qn.n, qn.l, qn.m] = parsed;
}

const orbitalSelect = requireElement('orbitalSel', HTMLSelectElement);
const particleInput = requireElement('nParts', HTMLInputElement);
const pointSizeInput = requireElement('ptSize', HTMLInputElement);
const animationToggle = requireElement('animToggle', HTMLInputElement);
const themeToggle = requireElement('themeToggle', HTMLInputElement);
const jsonInput = requireElement('jsonInput', HTMLInputElement);
const colorToggle = requireElement('colorToggle', HTMLInputElement);
const button2d = requireElement('btn2d', HTMLButtonElement);

orbitalSelect.addEventListener('change', () => {
  applyOrbitalOption(orbitalSelect.value);
  updateHUD();
  generateCloud();
});

requireElement('btnGen', HTMLButtonElement).addEventListener('click', () => {
  generateCloud();
});
requireElement('tabBtnPresets', HTMLButtonElement).addEventListener('click', () => {
  setPanelTab('presets');
});
requireElement('tabBtnAtom', HTMLButtonElement).addEventListener('click', () => {
  setPanelTab('atom');
});
requireElement('atomN', HTMLSelectElement).addEventListener('change', () => {
  refreshAtomForm(false);
});
requireElement('atomL', HTMLSelectElement).addEventListener('change', () => {
  refreshAtomForm(false);
});
requireElement('atomM', HTMLSelectElement).addEventListener('change', () => {
  refreshAtomForm(false);
});
requireElement('btnApplyAtom', HTMLButtonElement).addEventListener('click', applyAtomForm);
setPanelTab('presets');
applyTheme(getPreferredTheme(), false);

requireElement('btnRandom', HTMLButtonElement).addEventListener('click', () => {
  const opts = [...orbitalSelect.options];
  const pick = opts[Math.floor(Math.random() * opts.length)];
  if (!pick) return;
  orbitalSelect.value = pick.value;
  applyOrbitalOption(pick.value);
  updateHUD();
  generateCloud();
});

particleInput.addEventListener('change', () => {
  qn.N = Number(particleInput.value);
  requireElement('nPartsVal', HTMLOutputElement).textContent = qn.N.toLocaleString();
  generateCloud();
});
particleInput.addEventListener('input', () => {
  qn.N = Number(particleInput.value);
  requireElement('nPartsVal', HTMLOutputElement).textContent = qn.N.toLocaleString();
});

pointSizeInput.addEventListener('input', () => {
  ptSizeVal = Number(pointSizeInput.value) * 0.08;
  requireElement('ptSizeVal', HTMLOutputElement).textContent = pointSizeInput.value;
  if (cloudMesh) cloudMesh.material.size = ptSizeVal;
});

animationToggle.addEventListener('change', () => {
  animating = animationToggle.checked;
});

themeToggle.addEventListener('change', () => {
  applyTheme(themeToggle.checked ? 'light' : 'dark');
});

jsonInput.addEventListener('change', () => {
  const file = jsonInput.files?.[0];
  if (!file) return;
  showLoading('Chargement JSON...');
  void loadWavefunction(file).catch((error: unknown) => {
    hideLoading();
    const message = error instanceof Error ? error.message : String(error);
    alert(`Erreur JSON : ${message}`);
  });
  jsonInput.value = '';
});

colorToggle.addEventListener('change', () => {
  useAltColor = colorToggle.checked;
  generateCloud();
});

let show2D = false,
  sim2DInit = false;
button2d.addEventListener('click', () => {
  show2D = !show2D;
  requireElement('panel2d', HTMLElement).classList.toggle('visible', show2D);
  button2d.classList.toggle('active', show2D);
  if (show2D && !sim2DInit) init2D();
});

window.addEventListener('keydown', (e) => {
  if (e.repeat) return;
  const k = e.key.toLowerCase();
  let dirty = false;
  if (k === 'w') {
    qn.n++;
    dirty = true;
  } else if (k === 's') {
    if (qn.n > 1) qn.n--;
    dirty = true;
  } else if (k === 'e') {
    qn.l++;
    dirty = true;
  } else if (k === 'd') {
    if (qn.l > 0) qn.l--;
    dirty = true;
  } else if (k === 'r') {
    qn.m++;
    dirty = true;
  } else if (k === 'f') {
    qn.m--;
    dirty = true;
  } else if (k === 't') {
    qn.N = Math.min(qn.N * 2, 80000);
    particleInput.value = String(qn.N);
    requireElement('nPartsVal', HTMLOutputElement).textContent = qn.N.toLocaleString();
    dirty = true;
  } else if (k === 'g') {
    qn.N = Math.max(Math.round(qn.N / 2), 2000);
    particleInput.value = String(qn.N);
    requireElement('nPartsVal', HTMLOutputElement).textContent = qn.N.toLocaleString();
    dirty = true;
  } else if (k === 'a') {
    animating = !animating;
    animationToggle.checked = animating;
  } else if (k === 'q') {
    button2d.click();
  }
  if (dirty) {
    if (qn.l > qn.n - 1) qn.l = qn.n - 1;
    if (qn.l < 0) qn.l = 0;
    if (qn.m > qn.l) qn.m = qn.l;
    if (qn.m < -qn.l) qn.m = -qn.l;
    generateCloud();
  }
});

const DT = 0.4;
function animateFlow(): void {
  if (!cloudMesh || !animating || qn.m === 0) return;
  const pos = cloudMesh.geometry.getAttribute('position');
  if (!(pos instanceof THREE.BufferAttribute)) {
    throw new Error('Attribut de position du nuage indisponible.');
  }
  const m = qn.m;
  const availablePointCount = Math.min(qn.N, pos.count, rArr.length, Math.floor(posArr.length / 3));
  for (let i = 0; i < availablePointCount; i++) {
    const px = requireNumber(posArr, i * 3),
      py = requireNumber(posArr, i * 3 + 1),
      pz = requireNumber(posArr, i * 3 + 2);
    const [vx, , vz] = probFlow(px, py, pz, m);
    const nx = px + vx * DT,
      nz = pz + vz * DT;
    const nr = Math.sqrt(nx * nx + py * py + nz * nz);
    if (nr > 0.01) {
      const sc = requireNumber(rArr, i) / nr;
      posArr[i * 3] = nx * sc;
      posArr[i * 3 + 2] = nz * sc;
    }
    pos.setXYZ(
      i,
      requireNumber(posArr, i * 3),
      requireNumber(posArr, i * 3 + 1),
      requireNumber(posArr, i * 3 + 2),
    );
  }
  pos.needsUpdate = true;
}

let fpsTimes: number[] = [],
  fpsLast = performance.now();
function tickFPS(): void {
  const now = performance.now();
  fpsTimes.push(now);
  fpsTimes = fpsTimes.filter((t) => now - t < 1000);
  if (now - fpsLast > 400) {
    requireElement('iFps', HTMLElement).textContent = String(fpsTimes.length);
    fpsLast = now;
  }
}

let frame = 0;
function loop3D(): void {
  requestAnimationFrame(loop3D);
  frame++;
  tickFPS();
  if (cloudMesh) {
    if (animating && qn.m !== 0 && frame % 2 === 0) animateFlow();
    if (qn.m === 0 && animating) cloudMesh.rotation.y += 0.002;
  }
  nucleus.rotation.y += 0.01;
  renderer.render(scene, cam);
}

interface Vector2 {
  x: number;
  y: number;
}

function init2D(): void {
  sim2DInit = true;
  const cv = requireElement('c2d', HTMLCanvasElement);
  const ctx = requireCanvas2dContext(cv);
  const W2 = 300,
    H2 = 195,
    orbitDist = 15;

  class WavePoint {
    lp: Vector2;
    dir: Vector2;

    constructor(lp: Vector2, dir: Vector2) {
      this.lp = lp;
      this.dir = dir;
    }
  }

  class Wave {
    energy: number;
    sigma: number;
    k: number;
    phase: number;
    a: number;
    pos: Vector2;
    dir: Vector2;
    col: Color3;
    points: WavePoint[];

    constructor(energy: number, pos: Vector2, dir: Vector2, col: Color3 = [1, 0.6, 0.1]) {
      this.energy = energy;
      this.sigma = 30;
      this.k = 0.4;
      this.phase = 0;
      this.a = 6;
      this.pos = { ...pos };
      const len = Math.sqrt(dir.x * dir.x + dir.y * dir.y);
      this.dir = { x: dir.x / len, y: dir.y / len };
      this.col = col;
      this.points = [];
      for (let x = -this.sigma; x <= this.sigma; x += 0.5) {
        this.points.push(
          new WavePoint(
            { x: pos.x + x * this.dir.x, y: pos.y + x * this.dir.y },
            { x: this.dir.x * 120, y: this.dir.y * 120 },
          ),
        );
      }
    }
    draw(): void {
      ctx.strokeStyle = `rgba(${~~(this.col[0] * 255)},${~~(this.col[1] * 255)},${~~(this.col[2] * 255)},.9)`;
      ctx.lineWidth = 0.7;
      ctx.beginPath();
      let first = true;
      for (const p of this.points) {
        const perp = { x: -p.dir.y, y: p.dir.x };
        const pl = Math.sqrt(perp.x * perp.x + perp.y * perp.y) || 1;
        const yd =
          this.a * Math.sin(this.k * Math.sqrt(p.lp.x * p.lp.x + p.lp.y * p.lp.y) - this.phase);
        const dx = p.lp.x + (perp.x / pl) * yd,
          dy = p.lp.y + (perp.y / pl) * yd;
        if (first) {
          ctx.moveTo(dx + W2 / 2, dy + H2 / 2);
          first = false;
        } else ctx.lineTo(dx + W2 / 2, dy + H2 / 2);
      }
      ctx.stroke();
    }
    update(dt: number): boolean {
      this.phase += 20 * dt;
      for (const p of this.points) {
        p.lp.x += p.dir.x * dt;
        p.lp.y += p.dir.y * dt;
        if (Math.abs(p.lp.x) > W2 || Math.abs(p.lp.y) > H2) return true;
      }
      return false;
    }
  }

  class Atom2D {
    pos: Vector2;
    n: number;
    angle: number;

    constructor(pos: Vector2) {
      this.pos = { ...pos };
      this.n = 1;
      this.angle = Math.random() * Math.PI * 2;
    }
    draw(): void {
      ctx.strokeStyle = 'rgba(255,154,60,.18)';
      ctx.lineWidth = 0.5;
      ctx.beginPath();
      ctx.arc(this.pos.x + W2 / 2, this.pos.y + H2 / 2, this.n * orbitDist, 0, Math.PI * 2);
      ctx.stroke();
      const ex = Math.cos(this.angle) * this.n * orbitDist + this.pos.x + W2 / 2;
      const ey = Math.sin(this.angle) * this.n * orbitDist + this.pos.y + H2 / 2;
      ctx.fillStyle = '#ff9a3c';
      ctx.beginPath();
      ctx.arc(ex, ey, 1.5, 0, Math.PI * 2);
      ctx.fill();
      ctx.fillStyle = '#ff4020';
      ctx.beginPath();
      ctx.arc(this.pos.x + W2 / 2, this.pos.y + H2 / 2, 3, 0, Math.PI * 2);
      ctx.fill();
    }
    update(): void {
      this.angle += 0.06;
    }
  }

  const atoms2d: Atom2D[] = [];
  for (let i = 0; i < 12; i++) {
    const a = (2 * Math.PI * i) / 12;
    atoms2d.push(new Atom2D({ x: Math.cos(a) * 55, y: Math.sin(a) * 55 }));
  }

  const waves2d: Wave[] = [];
  const E12 = -13.6 / 4 - -13.6;
  for (let i = 0; i < 8; i++)
    waves2d.push(new Wave(E12, { x: 80, y: i * 22 - 76 }, { x: -1, y: 0 }));

  cv.addEventListener('click', (e) => {
    const rx = e.offsetX - W2 / 2,
      ry = e.offsetY - H2 / 2;
    for (let i = 0; i < 18; i++) {
      const a = Math.random() * 2 * Math.PI;
      waves2d.push(new Wave(E12, { x: rx, y: ry }, { x: Math.cos(a), y: Math.sin(a) }));
    }
  });

  function loop2D(): void {
    if (!show2D) {
      requestAnimationFrame(loop2D);
      return;
    }
    ctx.fillStyle = 'rgba(5,3,2,.85)';
    ctx.fillRect(0, 0, W2, H2);
    for (const a of atoms2d) {
      a.draw();
      a.update();
    }
    for (let i = 0; i < waves2d.length;) {
      const wave = waves2d[i];
      if (!wave) break;
      if (wave.energy === 0) {
        i++;
        continue;
      }
      wave.draw();
      if (wave.update(0.016)) waves2d.splice(i, 1);
      else i++;
    }
    requestAnimationFrame(loop2D);
  }
  loop2D();
}

window.addEventListener('resize', resizeCanvas);

resizeCanvas();
updateHUD();
generateCloud(() => {
  loop3D();
});
