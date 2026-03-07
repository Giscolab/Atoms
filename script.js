'use strict';

// â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
//  PHYSICS ENGINE
// â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

function gamma(z) {
  if (z < 0.5) return Math.PI / (Math.sin(Math.PI * z) * gamma(1 - z));
  z -= 1;
  const C = [0.99999999999980993,676.5203681218851,-1259.1392167224028,
    771.32342877765313,-176.61502916214059,12.507343278686905,
    -0.13857109526572012,9.9843695780195716e-6,1.5056327351493116e-7];
  let x = C[0];
  for (let i = 1; i < 9; i++) x += C[i] / (z + i);
  const t = z + 7.5;
  return Math.sqrt(2 * Math.PI) * Math.pow(t, z + 0.5) * Math.exp(-t) * x;
}

// â”€â”€â”€ Cached CDF tables â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
const _rCDF = {}, _tCDF = {};

function buildRCDF(n, l) {
  const key = `${n}_${l}`;
  if (_rCDF[key]) return _rCDF[key];
  const M = 4096, rMax = 10 * n * n;
  const dr = rMax / (M - 1);
  const cdf = new Float64Array(M);
  let sum = 0;
  const k = n - l - 1, alpha = 2 * l + 1;
  for (let i = 0; i < M; i++) {
    const r = i * dr, rho = 2 * r / n;
    let L = 1;
    if (k === 1) L = 1 + alpha - rho;
    else if (k > 1) {
      let Lm2 = 1, Lm1 = 1 + alpha - rho;
      for (let j = 2; j <= k; j++) {
        L = ((2*j-1+alpha-rho)*Lm1 - (j-1+alpha)*Lm2) / j;
        Lm2 = Lm1; Lm1 = L;
      }
    }
    const norm = Math.pow(2/n, 3) * gamma(n-l) / (2*n*gamma(n+l+1));
    const R = Math.sqrt(norm) * Math.exp(-rho/2) * Math.pow(rho||0, l) * L;
    const pdf = r * r * R * R;
    if (isFinite(pdf)) sum += pdf;
    cdf[i] = sum;
  }
  for (let i = 0; i < M; i++) cdf[i] /= sum;
  return (_rCDF[key] = { cdf, rMax, M });
}

function buildTCDF(l, absM) {
  const key = `${l}_${absM}`;
  if (_tCDF[key]) return _tCDF[key];
  const M = 2048;
  const dt = Math.PI / (M - 1);
  const cdf = new Float64Array(M);
  let sum = 0;
  for (let i = 0; i < M; i++) {
    const theta = i * dt, x = Math.cos(theta);
    let Pmm = 1;
    if (absM > 0) {
      const s = Math.sqrt((1-x)*(1+x));
      let fact = 1;
      for (let j = 1; j <= absM; j++) { Pmm *= -fact * s; fact += 2; }
    }
    let Plm = Pmm;
    if (l > absM) {
      let Pm1m = x * (2*absM+1) * Pmm;
      if (l === absM+1) { Plm = Pm1m; }
      else {
        let pp = Pmm;
        for (let ll = absM+2; ll <= l; ll++) {
          const Pll = ((2*ll-1)*x*Pm1m - (ll+absM-1)*pp) / (ll-absM);
          pp = Pm1m; Pm1m = Pll;
        }
        Plm = Pm1m;
      }
    }
    const pdf = Math.sin(theta) * Plm * Plm;
    if (isFinite(pdf)) sum += pdf;
    cdf[i] = sum;
  }
  for (let i = 0; i < M; i++) cdf[i] /= sum || 1;
  return (_tCDF[key] = { cdf, M });
}

function bsearch(cdf, u) {
  let lo = 0, hi = cdf.length - 1;
  while (lo < hi) { const mid = (lo+hi)>>1; if (cdf[mid] < u) lo=mid+1; else hi=mid; }
  return lo;
}

function sampleR(n, l) {
  const { cdf, rMax, M } = buildRCDF(n, l);
  return bsearch(cdf, Math.random()) * (rMax / (M-1));
}
function sampleTheta(l, m) {
  const { cdf, M } = buildTCDF(l, Math.abs(m));
  return bsearch(cdf, Math.random()) * (Math.PI / (M-1));
}
function samplePhi() { return 2 * Math.PI * Math.random(); }

// â”€â”€â”€ Rejection samplers analytiques (par orbital nommÃ©e) â”€
const _a0 = 52.9; // Bohr radius pm
function _rProb2p(r){ const x=r/_a0; return (x**4/24)*Math.exp(-x); }
function _rProb3p(r){ const x=r/_a0; const R=x*(1-x/6)*Math.exp(-x/3); return R*R*r*r; }
function _sampleR2p(){ const mx=15*_a0,pm=_rProb2p(4*_a0); for(;;){const r=Math.random()*mx; if(Math.random()<=_rProb2p(r)/pm)return r;} }
function _sampleR3p(){ const mx=25*_a0,pm=_rProb3p(8*_a0); for(;;){const r=Math.random()*mx; if(Math.random()<=_rProb3p(r)/pm)return r;} }
const _rT=()=>Math.acos(1-2*Math.random()), _rP=()=>2*Math.PI*Math.random();

const NAMED_SAMPLERS = {
  '2p_x': ()=>{ const r=_sampleR2p(); let t,p; for(;;){t=_rT();if(Math.random()<=Math.sin(t)**3)break;} for(;;){p=_rP();if(Math.random()<=Math.cos(p)**2)break;} return s2c(r,t,p); },
  '2p_y': ()=>{ const r=_sampleR2p(); let t,p; for(;;){t=_rT();if(Math.random()<=Math.sin(t)**3)break;} for(;;){p=_rP();if(Math.random()<=Math.sin(p)**2)break;} return s2c(r,t,p); },
  '2p_z': ()=>{ const r=_sampleR2p(); let t; for(;;){t=_rT();if(Math.random()<=Math.cos(t)**2)break;} return s2c(r,t,_rP()); },
  '3p_z': ()=>{ const r=_sampleR3p(); let t; for(;;){t=_rT();if(Math.random()<=Math.cos(t)**2)break;} return s2c(r,t,_rP()); }
};

// â”€â”€â”€ Mapping select value â†’ NAMED_SAMPLER key â”€â”€â”€â”€
const SAMPLER_MAP = {
  '2_1_0': '2p_z', '2_1_1': '2p_x', '2_1_-1': '2p_y',
  '3_1_0': '3p_z'
};

// â”€â”€â”€ JSON wavefunction loader (format: {points:[[x,y,z],â€¦]} en Bohr) â”€â”€
async function loadWavefunction(file) {
  return new Promise((resolve, reject) => {
    const reader = new FileReader();
    reader.onload = ev => {
      try {
        const data = JSON.parse(ev.target.result);
        if (!data.points) throw new Error('Missing "points" array');
        const BPM = 5.29; // Bohr â†’ pm
        const N = data.points.length;
        posArr = new Float32Array(N*3);
        rArr   = new Float32Array(N);
        const colArr = new Float32Array(N*3);
        data.points.forEach((p,i) => {
          const x=p[0]*BPM, y=p[1]*BPM, z=p[2]*BPM;
          posArr[i*3]=x; posArr[i*3+1]=y; posArr[i*3+2]=z;
          rArr[i]=Math.sqrt(x*x+y*y+z*z);
          const th=Math.acos(Math.max(-1,Math.min(1,y/Math.max(rArr[i],1e-6))));
          const ph=Math.atan2(z,x);
          const [cr,cg,cb]=orbitalColor(rArr[i],th,ph,qn.n,qn.l,qn.m);
          colArr[i*3]=cr; colArr[i*3+1]=cg; colArr[i*3+2]=cb;
        });
        if (cloudMesh) scene.remove(cloudMesh);
        const geo = new THREE.BufferGeometry();
        geo.setAttribute('position', new THREE.BufferAttribute(posArr.slice(),3));
        geo.setAttribute('color',    new THREE.BufferAttribute(colArr,3));
        cloudMesh = new THREE.Points(geo, new THREE.PointsMaterial({
          size:ptSizeVal, vertexColors:true, transparent:true,
          opacity:.88, sizeAttenuation:true, depthWrite:false
        }));
        scene.add(cloudMesh);
        qn.N = N;
        document.getElementById('nPartsVal').textContent = N.toLocaleString();
        updateHUD(); hideLoading();
        resolve();
      } catch(e){ reject(e); }
    };
    reader.onerror = reject;
    reader.readAsText(file);
  });
}

function s2c(r, th, ph) {
  const s = Math.sin(th);
  return [r*s*Math.cos(ph), r*Math.cos(th), r*s*Math.sin(ph)];
}

function heatFire(v) {
  v = Math.max(0, Math.min(1, v));
  const s = [[0,0,0],[0.5,0,0.99],[0.8,0,0],[1,.5,0],[1,1,0],[1,1,1]];
  const sv = v * 5, i = Math.min(~~sv, 4), t = sv - i;
  return [s[i][0]+t*(s[i+1][0]-s[i][0]),s[i][1]+t*(s[i+1][1]-s[i][1]),s[i][2]+t*(s[i+1][2]-s[i][2])];
}

// â”€â”€â”€ Alternative colormap: yellow â†’ orange/red â†’ purple â”€â”€
function densityColor(v) {
  v = Math.max(0, Math.min(1, v));
  let r, g, b;
  if (v < 0.5) { const f=v/0.5;      r=1;       g=1-0.5*f; b=0; }
  else         { const f=(v-0.5)/0.5; r=1-0.5*f; g=0.5*(1-f); b=0.5*f; }
  return [r, g, b];
}
let useAltColor = false; // toggled by UI

function orbitalColor(r, theta, phi, n, l, m) {
  const rho = 2*r/n, k = n-l-1, alpha = 2*l+1;
  let L = 1;
  if (k===1) L = 1+alpha-rho;
  else if (k>1) {
    let Lm2=1, Lm1=1+alpha-rho;
    for (let j=2;j<=k;j++){L=((2*j-1+alpha-rho)*Lm1-(j-1+alpha)*Lm2)/j;Lm2=Lm1;Lm1=L;}
  }
  const norm = Math.pow(2/n,3)*gamma(n-l)/(2*n*gamma(n+l+1));
  const R = Math.sqrt(norm)*Math.exp(-rho/2)*Math.pow(rho||0,l)*L;
  const absM = Math.abs(m), x = Math.cos(theta);
  let Pmm=1;
  if (absM>0){const s=Math.sqrt((1-x)*(1+x));let f=1;for(let j=1;j<=absM;j++){Pmm*=-f*s;f+=2;}}
  let Plm=Pmm;
  if (l>absM){
    let Pm1m=x*(2*absM+1)*Pmm;
    if(l===absM+1){Plm=Pm1m;}
    else{let pp=Pmm;for(let ll=absM+2;ll<=l;ll++){const Pll=((2*ll-1)*x*Pm1m-(ll+absM-1)*pp)/(ll-absM);pp=Pm1m;Pm1m=Pll;}Plm=Pm1m;}
  }
  const intensity = R*R*Plm*Plm;
  const t = intensity * 1.5 * Math.pow(5, n);
  return useAltColor ? densityColor(Math.min(t,1)) : heatFire(t);
}

function probFlow(px, py, pz, m) {
  const r = Math.sqrt(px*px+py*py+pz*pz);
  if (r < 1e-6) return [0,0,0];
  const theta = Math.acos(Math.max(-1,Math.min(1,py/r)));
  const phi   = Math.atan2(pz, px);
  const st    = Math.sin(theta);
  if (Math.abs(st) < 1e-4) return [0,0,0];
  const vm = m / (r * st);
  return [-vm*Math.sin(phi), 0, vm*Math.cos(phi)];
}

// â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
//  THREE.JS 3D RENDERER
// â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•

const canvas3d = document.getElementById('atomSimCanvas');

const renderer = new THREE.WebGLRenderer({ canvas: canvas3d, antialias: true, alpha: false });
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 1.5));
renderer.setClearColor(0x050302, 1);

const scene = new THREE.Scene();
const cam   = new THREE.PerspectiveCamera(60, 800/600, 0.01, 2000);

// â”€â”€â”€ Fit canvas to viewport â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
function resizeCanvas() {
  const vp = document.getElementById('viewport');
  const w  = vp.clientWidth  - 40;
  const h  = vp.clientHeight - 40;
  canvas3d.width  = w;
  canvas3d.height = h;
  cam.aspect = w / h;
  cam.updateProjectionMatrix();
  renderer.setSize(w, h);
}

// â”€â”€â”€ Orbit state â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
let orb = { az: 0.3, el: Math.PI/2.5, r: 60 };
let dragging = false, lastX = 0, lastY = 0;

function updateCamera() {
  const el = Math.max(0.02, Math.min(Math.PI-0.02, orb.el));
  cam.position.set(
    orb.r*Math.sin(el)*Math.cos(orb.az),
    orb.r*Math.cos(el),
    orb.r*Math.sin(el)*Math.sin(orb.az)
  );
  cam.lookAt(0, 0, 0);
}
updateCamera();

canvas3d.addEventListener('mousedown', e => { dragging=true; lastX=e.clientX; lastY=e.clientY; });
window.addEventListener('mouseup',  () => dragging=false);
window.addEventListener('mousemove', e => {
  if (!dragging) return;
  orb.az += (e.clientX-lastX)*0.008;
  orb.el -= (e.clientY-lastY)*0.008;
  lastX=e.clientX; lastY=e.clientY;
  updateCamera();
  document.getElementById('iDist').textContent = orb.r.toFixed(1);
});
canvas3d.addEventListener('wheel', e => {
  orb.r = Math.max(5, Math.min(400, orb.r + e.deltaY*0.05));
  updateCamera();
  document.getElementById('iDist').textContent = orb.r.toFixed(1);
  e.preventDefault();
}, {passive:false});

let tLast = null;
canvas3d.addEventListener('touchstart', e => { tLast = e.touches[0]; });
canvas3d.addEventListener('touchmove',  e => {
  if (!tLast) return;
  const t = e.touches[0];
  orb.az += (t.clientX-tLast.clientX)*0.01;
  orb.el -= (t.clientY-tLast.clientY)*0.01;
  tLast = t; updateCamera(); e.preventDefault();
}, {passive:false});

// â”€â”€â”€ Stars background â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
(function(){
  const N=1200, pos=new Float32Array(N*3), col=new Float32Array(N*3);
  for(let i=0;i<N;i++){
    const th=Math.acos(2*Math.random()-1), ph=2*Math.PI*Math.random(), r=400+Math.random()*200;
    pos[i*3]=r*Math.sin(th)*Math.cos(ph); pos[i*3+1]=r*Math.cos(th); pos[i*3+2]=r*Math.sin(th)*Math.sin(ph);
    const b=0.08+Math.random()*0.25; col[i*3]=b*.9; col[i*3+1]=b*.7; col[i*3+2]=b*.3;
  }
  const g=new THREE.BufferGeometry();
  g.setAttribute('position',new THREE.BufferAttribute(pos,3));
  g.setAttribute('color',   new THREE.BufferAttribute(col,3));
  scene.add(new THREE.Points(g,new THREE.PointsMaterial({size:.4,vertexColors:true,sizeAttenuation:true})));
})();

// â”€â”€â”€ Nucleus â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
const nucleus = new THREE.Mesh(
  new THREE.SphereGeometry(0.5,16,16),
  new THREE.MeshBasicMaterial({color:0xff6030})
);
scene.add(nucleus);
const ring = new THREE.Mesh(
  new THREE.RingGeometry(0.7,1.0,32),
  new THREE.MeshBasicMaterial({color:0xff4020,side:THREE.DoubleSide,transparent:true,opacity:.3})
);
scene.add(ring);

// â”€â”€â”€ Cloud state â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
let cloudMesh = null;
let posArr = null, rArr = null;
let qn = {n:2, l:1, m:0, N:15000};
let animating  = true;
let ptSizeVal  = 0.24;
const THEME_STORAGE_KEY = 'atom-theme';
const THEME_CLEAR_COLOR = { dark: 0x050302, light: 0xf3ece4 };

function getPreferredTheme() {
  try {
    const storedTheme = localStorage.getItem(THEME_STORAGE_KEY);
    if (storedTheme === 'light' || storedTheme === 'dark') return storedTheme;
  } catch (e) {}
  if (window.matchMedia && window.matchMedia('(prefers-color-scheme: light)').matches) {
    return 'light';
  }
  return 'dark';
}

function applyTheme(theme, persist = true) {
  const normalizedTheme = theme === 'light' ? 'light' : 'dark';
  document.documentElement.setAttribute('data-theme', normalizedTheme);
  const themeToggle = document.getElementById('themeToggle');
  if (themeToggle) themeToggle.checked = normalizedTheme === 'light';
  renderer.setClearColor(THEME_CLEAR_COLOR[normalizedTheme], 1);
  if (persist) {
    try {
      localStorage.setItem(THEME_STORAGE_KEY, normalizedTheme);
    } catch (e) {}
  }
}

const ORBITAL_LETTERS = ['s','p','d','f','g','h','i'];
const ORBITAL_FAMILY = {
  s: 'Distribution spherique autour du noyau.',
  p: 'Deux lobes principaux avec un plan nodal.',
  d: 'Forme complexe a plusieurs lobes (type trefle).',
  f: 'Forme multi-lobes avec structure plus fine.',
  g: 'Orbitale de haut ordre avec nombreux noeuds angulaires.',
  h: 'Orbitale de haut ordre, tres structuree.',
  i: 'Orbitale de tres haut ordre et forte complexite angulaire.'
};
const PANEL_MIN_N = 1;
const PANEL_MAX_N = 9;

function clampValue(v, min, max) {
  return Math.max(min, Math.min(max, v));
}

function orbitalLetter(l) {
  return ORBITAL_LETTERS[l] || `l=${l}`;
}

function orbitalFamilySummary(l) {
  const letter = orbitalLetter(l);
  return ORBITAL_FAMILY[letter] || `Famille orbitale definie par l=${l}.`;
}

function magneticProjectionSummary(m) {
  if (m === 0) return 'Symetrie axiale (m=0).';
  if (m > 0) return `Projection positive (m=+${m}).`;
  return `Projection negative (m=${m}).`;
}

function optionText(option) {
  if (!option) return '';
  return option.textContent.replace(/\s+/g, ' ').trim();
}

function intRange(min, max) {
  const values = [];
  for (let i = min; i <= max; i++) values.push(i);
  return values;
}

function setSelectOptions(select, values, selectedValue) {
  select.innerHTML = '';
  for (const value of values) {
    const option = document.createElement('option');
    option.value = String(value);
    option.textContent = String(value);
    select.appendChild(option);
  }
  select.value = String(selectedValue);
}

function refreshAtomForm(fromState = true) {
  const nSel = document.getElementById('atomN');
  const lSel = document.getElementById('atomL');
  const mSel = document.getElementById('atomM');
  if (!nSel || !lSel || !mSel) return;

  const requestedN = fromState ? qn.n : Number(nSel.value || qn.n);
  const maxN = Math.max(PANEL_MAX_N, qn.n);
  const nVal = clampValue(requestedN, PANEL_MIN_N, maxN);
  setSelectOptions(nSel, intRange(PANEL_MIN_N, maxN), nVal);

  const requestedL = fromState ? qn.l : Number(lSel.value || qn.l);
  const lMax = Math.max(0, nVal - 1);
  const lVal = clampValue(requestedL, 0, lMax);
  setSelectOptions(lSel, intRange(0, lMax), lVal);

  const requestedM = fromState ? qn.m : Number(mSel.value || qn.m);
  const mVal = clampValue(requestedM, -lVal, lVal);
  setSelectOptions(mSel, intRange(-lVal, lVal), mVal);
}

function updatePresetInfo() {
  const { n, l, m } = qn;
  const letter = orbitalLetter(l);
  const radialNodes = Math.max(0, n - l - 1);
  const angularNodes = l;

  const sel = document.getElementById('orbitalSel');
  const matchValue = `${n}_${l}_${m}`;
  const selected = sel ? [...sel.options].find(o => o.value === matchValue) : null;
  const label = selected ? optionText(selected) : `${n}${letter} (m=${m})`;
  const sourceTag = selected ? 'Prereglage du catalogue' : "Etat defini via l'onglet Atome";

  document.getElementById('presetName').textContent = label;
  document.getElementById('presetSummary').textContent = `${orbitalFamilySummary(l)} ${sourceTag}.`;
  document.getElementById('presetFamily').textContent = letter;
  document.getElementById('presetNumbers').textContent = `n=${n}, l=${l}, m=${m}`;
  document.getElementById('presetNodes').textContent = `${radialNodes} radial, ${angularNodes} angulaire`;
  document.getElementById('presetProjection').textContent = magneticProjectionSummary(m);
  document.getElementById('presetLock').textContent = 'Contraintes fixes: n>=1, 0<=l<=n-1, -l<=m<=l.';
}

function setPanelTab(tabId) {
  const showPresets = tabId === 'presets';
  document.getElementById('panel-view-presets').classList.toggle('active', showPresets);
  document.getElementById('panel-view-atom').classList.toggle('active', !showPresets);
  document.getElementById('tabBtnPresets').setAttribute('aria-selected', showPresets ? 'true' : 'false');
  document.getElementById('tabBtnAtom').setAttribute('aria-selected', showPresets ? 'false' : 'true');
}

function applyAtomForm() {
  const nVal = Number(document.getElementById('atomN').value);
  const lVal = Number(document.getElementById('atomL').value);
  const mVal = Number(document.getElementById('atomM').value);
  qn.n = Math.max(1, nVal || 1);
  qn.l = Math.max(0, Math.min(qn.n - 1, lVal || 0));
  qn.m = clampValue(mVal || 0, -qn.l, qn.l);
  updateHUD();
  generateCloud();
}

// â”€â”€â”€ Loading helpers â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
function showLoading(msg) {
  const ld = document.getElementById('loading');
  ld.style.display = 'block';
  ld.textContent   = msg || 'Calcul de l\'orbitale...';
}
function hideLoading() {
  document.getElementById('loading').style.display    = 'none';
  document.getElementById('progressBar').textContent  = '';
}

// â”€â”€â”€ Circular dot texture (anti square-pixel) â”€â”€â”€â”€
const _dotTex = (() => {
  const c = document.createElement('canvas'); c.width = c.height = 16;
  const x = c.getContext('2d');
  const g = x.createRadialGradient(8,8,0, 8,8,8);
  g.addColorStop(0,   'rgba(255,255,255,1)');
  g.addColorStop(0.5, 'rgba(255,255,255,0.8)');
  g.addColorStop(1,   'rgba(255,255,255,0)');
  x.fillStyle = g; x.fillRect(0,0,16,16);
  return new THREE.CanvasTexture(c);
})();

// â”€â”€â”€ Generate cloud (chunked, non-blocking) â”€â”€â”€â”€â”€â”€
function generateCloud(callback) {
  const { n, l, m, N } = qn;
  showLoading('Calcul de l\'orbitale...');

  setTimeout(() => {
    posArr = new Float32Array(N*3);
    rArr   = new Float32Array(N);
    const colArr = new Float32Array(N*3);
    let done = 0;
    const chunk = 2000;

    function doChunk() {
      const end = Math.min(done + chunk, N);
      for (let i = done; i < end; i++) {
        let x, y, z, r, theta, phi;
        const namedFn = NAMED_SAMPLERS[SAMPLER_MAP[`${n}_${l}_${m}`]];
        if (namedFn) {
          [x,y,z] = namedFn();
          r = Math.sqrt(x*x+y*y+z*z);
          theta = Math.acos(Math.max(-1,Math.min(1,y/Math.max(r,1e-9))));
          phi   = Math.atan2(z,x);
        } else {
          r     = sampleR(n, l);
          theta = sampleTheta(l, m);
          phi   = samplePhi();
          [x,y,z] = s2c(r, theta, phi);
        }
        posArr[i*3]=x; posArr[i*3+1]=y; posArr[i*3+2]=z;
        rArr[i] = r;
        const [cr,cg,cb] = orbitalColor(r,theta,phi,n,l,m);
        colArr[i*3]=cr; colArr[i*3+1]=cg; colArr[i*3+2]=cb;
      }
      done = end;
      document.getElementById('progressBar').textContent =
        `${Math.round(done/N*100)}%  â€”  ${done.toLocaleString()} / ${N.toLocaleString()}`;
      if (done < N) { setTimeout(doChunk, 0); return; }

      // Build Three.js geometry
      setTimeout(() => {
        if (cloudMesh) scene.remove(cloudMesh);
        const geo = new THREE.BufferGeometry();
        geo.setAttribute('position', new THREE.BufferAttribute(posArr.slice(), 3));
        geo.setAttribute('color',    new THREE.BufferAttribute(colArr, 3));
        cloudMesh = new THREE.Points(geo, new THREE.PointsMaterial({
          size: ptSizeVal, vertexColors: true,
          transparent: true, opacity: 0.88,
          sizeAttenuation: true, depthWrite: false,
          map: _dotTex, alphaTest: 0.01
        }));
        scene.add(cloudMesh);
        hideLoading();
        updateHUD();
        if (callback) callback();
      }, 0);
    }
    doChunk();
  }, 50);
}

// â”€â”€â”€ HUD / Panel sync â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
function updateHUD() {
  const { n, l, m, N } = qn;
  document.getElementById('qn-n').textContent      = n;
  document.getElementById('qn-l').textContent      = l;
  document.getElementById('qn-m').textContent      = m;
  document.getElementById('energy-val').textContent = `${(-13.6/(n*n)).toFixed(3)} eV`;
  document.getElementById('iOrb').textContent      = `${n}, ${l}, ${m}`;
  document.getElementById('iPts').textContent      = N.toLocaleString();
  document.getElementById('iDist').textContent     = orb.r.toFixed(1);
  const sel = document.getElementById('orbitalSel');
  const match = `${n}_${l}_${m}`;
  const hasMatch = [...sel.options].some(o => o.value === match);
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

// â”€â”€â”€ Panel controls wiring â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
document.getElementById('orbitalSel').addEventListener('change', function() {
  const p = this.value.split('_').map(Number);
  qn.n = p[0]; qn.l = p[1]; qn.m = p[2];
  updateHUD();
});

document.getElementById('btnGen').addEventListener('click', () => generateCloud());
document.getElementById('tabBtnPresets').addEventListener('click', () => setPanelTab('presets'));
document.getElementById('tabBtnAtom').addEventListener('click', () => setPanelTab('atom'));
document.getElementById('atomN').addEventListener('change', () => refreshAtomForm(false));
document.getElementById('atomL').addEventListener('change', () => refreshAtomForm(false));
document.getElementById('atomM').addEventListener('change', () => refreshAtomForm(false));
document.getElementById('btnApplyAtom').addEventListener('click', applyAtomForm);
setPanelTab('presets');
applyTheme(getPreferredTheme(), false);

// â”€â”€â”€ Random orbital button â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
document.getElementById('btnRandom').addEventListener('click', () => {
  const opts = [...document.getElementById('orbitalSel').options];
  const pick = opts[Math.floor(Math.random() * opts.length)];
  document.getElementById('orbitalSel').value = pick.value;
  const p = pick.value.split('_').map(Number);
  qn.n = p[0]; qn.l = p[1]; qn.m = p[2];
  generateCloud();
});

document.getElementById('nParts').addEventListener('change', function() {
  qn.N = +this.value;
  document.getElementById('nPartsVal').textContent = qn.N.toLocaleString();
  generateCloud();
});
document.getElementById('nParts').addEventListener('input', function() {
  qn.N = +this.value;
  document.getElementById('nPartsVal').textContent = qn.N.toLocaleString();
});

document.getElementById('ptSize').addEventListener('input', function() {
  ptSizeVal = +this.value * 0.08;
  document.getElementById('ptSizeVal').textContent = this.value;
  if (cloudMesh) cloudMesh.material.size = ptSizeVal;
});

document.getElementById('animToggle').addEventListener('change', function() {
  animating = this.checked;
});

document.getElementById('themeToggle').addEventListener('change', function() {
  applyTheme(this.checked ? 'light' : 'dark');
});

// â”€â”€â”€ JSON import wiring â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
document.getElementById('jsonInput').addEventListener('change', function() {
  const file = this.files[0]; if (!file) return;
  showLoading('Chargement JSON...');
  loadWavefunction(file).catch(e => { hideLoading(); alert('Erreur JSON : ' + e.message); });
  this.value = '';
});

// â”€â”€â”€ Alt colormap toggle â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
document.getElementById('colorToggle').addEventListener('change', function() {
  useAltColor = this.checked;
  generateCloud();
});

// â”€â”€â”€ 2D panel toggle â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
let show2D = false, sim2DInit = false;
document.getElementById('btn2d').addEventListener('click', () => {
  show2D = !show2D;
  document.getElementById('panel2d').classList.toggle('visible', show2D);
  document.getElementById('btn2d').classList.toggle('active',   show2D);
  if (show2D && !sim2DInit) init2D();
});

// â”€â”€â”€ Keyboard controls â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
window.addEventListener('keydown', e => {
  if (e.repeat) return;
  const k = e.key.toLowerCase();
  let dirty = false;
  if      (k==='w'){qn.n++;dirty=true;}
  else if (k==='s'){if(qn.n>1)qn.n--;dirty=true;}
  else if (k==='e'){qn.l++;dirty=true;}
  else if (k==='d'){if(qn.l>0)qn.l--;dirty=true;}
  else if (k==='r'){qn.m++;dirty=true;}
  else if (k==='f'){qn.m--;dirty=true;}
  else if (k==='t'){
    qn.N=Math.min(qn.N*2,80000);
    document.getElementById('nParts').value=qn.N;
    document.getElementById('nPartsVal').textContent=qn.N.toLocaleString();
    dirty=true;
  } else if (k==='g'){
    qn.N=Math.max(Math.round(qn.N/2),2000);
    document.getElementById('nParts').value=qn.N;
    document.getElementById('nPartsVal').textContent=qn.N.toLocaleString();
    dirty=true;
  } else if (k==='a'){
    animating=!animating;
    document.getElementById('animToggle').checked=animating;
  } else if (k==='q'){
    document.getElementById('btn2d').click();
  }
  if (dirty) {
    if(qn.l>qn.n-1) qn.l=qn.n-1;
    if(qn.l<0)      qn.l=0;
    if(qn.m>qn.l)   qn.m=qn.l;
    if(qn.m<-qn.l)  qn.m=-qn.l;
    generateCloud();
  }
});

// â”€â”€â”€ Animate probability current â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
const DT = 0.4;
function animateFlow() {
  if (!cloudMesh || !animating || qn.m === 0) return;
  const pos = cloudMesh.geometry.attributes.position;
  const m   = qn.m;
  for (let i = 0; i < qn.N; i++) {
    const px=posArr[i*3], py=posArr[i*3+1], pz=posArr[i*3+2];
    const [vx,,vz] = probFlow(px,py,pz,m);
    const nx=px+vx*DT, nz=pz+vz*DT;
    const nr=Math.sqrt(nx*nx+py*py+nz*nz);
    if(nr>0.01){
      const sc=rArr[i]/nr;
      posArr[i*3]=nx*sc; posArr[i*3+2]=nz*sc;
    }
    pos.array[i*3]=posArr[i*3];
    pos.array[i*3+1]=posArr[i*3+1];
    pos.array[i*3+2]=posArr[i*3+2];
  }
  pos.needsUpdate = true;
}

// â”€â”€â”€ FPS tracker â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
let fpsTimes = [], fpsLast = performance.now();
function tickFPS() {
  const now = performance.now();
  fpsTimes.push(now);
  fpsTimes = fpsTimes.filter(t => now - t < 1000);
  if (now - fpsLast > 400) {
    document.getElementById('iFps').textContent = fpsTimes.length;
    fpsLast = now;
  }
}

// â”€â”€â”€ 3D render loop â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
let frame = 0;
function loop3D() {
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

// â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
//  2D PHOTOELECTRIC SIMULATION
// â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•â•
function init2D() {
  sim2DInit = true;
  const cv  = document.getElementById('c2d');
  const ctx = cv.getContext('2d');
  const W2=300, H2=195, orbitDist=15;

  class WavePoint { constructor(lp,dir){this.lp=lp;this.dir=dir;} }

  class Wave {
    constructor(energy,pos,dir,col=[1,.6,.1]){
      this.energy=energy;this.sigma=30;this.k=0.4;
      this.phase=0;this.a=6;this.pos={...pos};
      const len=Math.sqrt(dir.x*dir.x+dir.y*dir.y);
      this.dir={x:dir.x/len,y:dir.y/len};this.col=col;this.points=[];
      for(let x=-this.sigma;x<=this.sigma;x+=0.5){
        this.points.push(new WavePoint(
          {x:pos.x+x*this.dir.x,y:pos.y+x*this.dir.y},
          {x:this.dir.x*120,y:this.dir.y*120}
        ));
      }
    }
    draw(){
      ctx.strokeStyle=`rgba(${~~(this.col[0]*255)},${~~(this.col[1]*255)},${~~(this.col[2]*255)},.9)`;
      ctx.lineWidth=0.7;ctx.beginPath();let first=true;
      for(const p of this.points){
        const perp={x:-p.dir.y,y:p.dir.x};
        const pl=Math.sqrt(perp.x*perp.x+perp.y*perp.y)||1;
        const yd=this.a*Math.sin(this.k*Math.sqrt(p.lp.x*p.lp.x+p.lp.y*p.lp.y)-this.phase);
        const dx=p.lp.x+perp.x/pl*yd, dy=p.lp.y+perp.y/pl*yd;
        if(first){ctx.moveTo(dx+W2/2,dy+H2/2);first=false;}
        else ctx.lineTo(dx+W2/2,dy+H2/2);
      }
      ctx.stroke();
    }
    update(dt){
      this.phase+=20*dt;
      for(const p of this.points){
        p.lp.x+=p.dir.x*dt; p.lp.y+=p.dir.y*dt;
        if(Math.abs(p.lp.x)>W2||Math.abs(p.lp.y)>H2) return true;
      }
      return false;
    }
  }

  class Atom2D {
    constructor(pos){this.pos={...pos};this.n=1;this.angle=Math.random()*Math.PI*2;}
    draw(){
      ctx.strokeStyle='rgba(255,154,60,.18)';ctx.lineWidth=0.5;
      ctx.beginPath();ctx.arc(this.pos.x+W2/2,this.pos.y+H2/2,this.n*orbitDist,0,Math.PI*2);ctx.stroke();
      const ex=Math.cos(this.angle)*this.n*orbitDist+this.pos.x+W2/2;
      const ey=Math.sin(this.angle)*this.n*orbitDist+this.pos.y+H2/2;
      ctx.fillStyle='#ff9a3c';ctx.beginPath();ctx.arc(ex,ey,1.5,0,Math.PI*2);ctx.fill();
      ctx.fillStyle='#ff4020';ctx.beginPath();ctx.arc(this.pos.x+W2/2,this.pos.y+H2/2,3,0,Math.PI*2);ctx.fill();
    }
    update(){this.angle+=0.06;}
  }

  const atoms2d=[];
  for(let i=0;i<12;i++){
    const a=2*Math.PI*i/12;
    atoms2d.push(new Atom2D({x:Math.cos(a)*55, y:Math.sin(a)*55}));
  }

  const waves2d=[];
  const E12=-13.6/4-(-13.6);
  for(let i=0;i<8;i++) waves2d.push(new Wave(E12,{x:80,y:i*22-76},{x:-1,y:0}));

  cv.addEventListener('click', e => {
    const rx=e.offsetX-W2/2, ry=e.offsetY-H2/2;
    for(let i=0;i<18;i++){
      const a=Math.random()*2*Math.PI;
      waves2d.push(new Wave(E12,{x:rx,y:ry},{x:Math.cos(a),y:Math.sin(a)}));
    }
  });

  function loop2D(){
    if(!show2D){requestAnimationFrame(loop2D);return;}
    ctx.fillStyle='rgba(5,3,2,.85)';ctx.fillRect(0,0,W2,H2);
    for(const a of atoms2d){a.draw();a.update();}
    for(let i=0;i<waves2d.length;){
      if(waves2d[i].energy===0){i++;continue;}
      waves2d[i].draw();
      if(waves2d[i].update(0.016)) waves2d.splice(i,1); else i++;
    }
    requestAnimationFrame(loop2D);
  }
  loop2D();
}

// â”€â”€â”€ Resize â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
window.addEventListener('resize', resizeCanvas);

// â”€â”€â”€ Boot â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€â”€
resizeCanvas();
updateHUD();
generateCloud(() => { loop3D(); });
