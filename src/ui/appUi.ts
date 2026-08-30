import type { OrbitalChartData } from '../sampling/orbitalCharts';
import { REAL_ORBITAL_DEFINITIONS, type RealOrbitalName } from '../science/hydrogen/realOrbitals';
import type { OrbitalRenderingState, AppState } from '../state/appState';
import type { OrbitalPresentation } from '../app/orbitalPresentation';
import type { OrbitalSamplingState } from '../sampling/contracts';
import { renderAngularChart, renderRadialChart } from './chartRenderer';
import { requireElement } from './dom';

const THEME_STORAGE_KEY = 'atoms-theme';
const UI_MAX_N = 9;

export interface UiCallbacks {
  generate(): void;
  resetCamera(): void;
  setOrbital(orbital: OrbitalSamplingState): void;
  setSampling(sampling: AppState['sampling']): void;
  setRendering(patch: Partial<OrbitalRenderingState>): void;
  setLegacy2DVisible(visible: boolean): void;
}

export interface AppUi {
  bind(callbacks: UiCallbacks): void;
  finishGeneration(message?: string): void;
  renderCharts(charts: OrbitalChartData): void;
  renderState(state: AppState, presentation: OrbitalPresentation, nodesAvailable?: boolean): void;
  setEngineStatus(message: string, state?: 'ready' | 'busy' | 'error'): void;
  showGeneration(message?: string): void;
  updateFps(fps: number): void;
  updateHud(state: AppState, cameraDistance: number): void;
  updateProgress(completed: number, total: number): void;
}

function integerOptions(
  select: HTMLSelectElement,
  min: number,
  max: number,
  selected: number,
): void {
  select.replaceChildren();
  for (let value = min; value <= max; value += 1) {
    const option = document.createElement('option');
    option.value = String(value);
    option.textContent = String(value);
    option.selected = value === selected;
    select.appendChild(option);
  }
}

function finiteInteger(value: string, fallback: number, min: number, max: number): number {
  const parsed = Number(value);
  if (!Number.isSafeInteger(parsed)) return fallback;
  return Math.max(min, Math.min(max, parsed));
}

function formatInteger(value: number): string {
  return value.toLocaleString('fr-FR');
}

function formatDecimal(value: number, digits: number): string {
  return value.toLocaleString('fr-FR', {
    maximumFractionDigits: digits,
    minimumFractionDigits: digits,
  });
}

function selectedRealOrbital(select: HTMLSelectElement): RealOrbitalName {
  const value = select.value;
  if (Object.hasOwn(REAL_ORBITAL_DEFINITIONS, value)) return value as RealOrbitalName;
  return 'd_xy';
}

function randomUint32(): number {
  const values = new Uint32Array(1);
  try {
    globalThis.crypto.getRandomValues(values);
    return values[0] ?? 0;
  } catch {
    // Un contexte local sans Crypto API conserve une seed non nulle bornée.
  }
  return (Date.now() >>> 0) ^ 0x41544f4d;
}

function themeFromDocument(): 'dark' | 'light' {
  return document.documentElement.getAttribute('data-theme') === 'light' ? 'light' : 'dark';
}

export function createAppUi(initialState: AppState): AppUi {
  let currentState = initialState;
  let callbacks: UiCallbacks | null = null;

  const quantumN = requireElement('quantumN', HTMLSelectElement);
  const quantumL = requireElement('quantumL', HTMLSelectElement);
  const quantumM = requireElement('quantumM', HTMLSelectElement);
  const realOrbital = requireElement('realOrbital', HTMLSelectElement);
  const basisComplex = requireElement('basisComplex', HTMLInputElement);
  const basisReal = requireElement('basisReal', HTMLInputElement);
  const complexMField = requireElement('complexMField', HTMLElement);
  const realOrbitalField = requireElement('realOrbitalField', HTMLElement);
  const sampleCount = requireElement('sampleCount', HTMLInputElement);
  const seedInput = requireElement('seedInput', HTMLInputElement);
  const isoThreshold = requireElement('isoThreshold', HTMLInputElement);
  const pointOpacity = requireElement('pointOpacity', HTMLInputElement);
  const pointSize = requireElement('pointSize', HTMLInputElement);
  const nodesToggle = requireElement('nodesToggle', HTMLInputElement);
  const generationStatus = requireElement('generationStatus', HTMLElement);
  const progressTrack = generationStatus.parentElement?.querySelector('.progress-track');

  function requireCallbacks(): UiCallbacks {
    if (!callbacks) throw new Error('L’interface Atoms doit être liée avant utilisation.');
    return callbacks;
  }

  function applyTheme(theme: 'dark' | 'light'): void {
    document.documentElement.setAttribute('data-theme', theme);
    const darkButton = requireElement('themeDark', HTMLButtonElement);
    const lightButton = requireElement('themeLight', HTMLButtonElement);
    darkButton.classList.toggle('active', theme === 'dark');
    lightButton.classList.toggle('active', theme === 'light');
    darkButton.setAttribute('aria-pressed', String(theme === 'dark'));
    lightButton.setAttribute('aria-pressed', String(theme === 'light'));
    try {
      localStorage.setItem(THEME_STORAGE_KEY, theme);
    } catch {
      /* stockage indisponible : le thème reste local */
    }
  }

  function updateMOptions(l: number, selected: number): void {
    integerOptions(quantumM, -l, l, Math.max(-l, Math.min(l, selected)));
  }

  function renderQuantumControls(state: AppState): void {
    const { orbital } = state;
    const n = orbital.n;
    const l = orbital.basis === 'complex' ? orbital.l : REAL_ORBITAL_DEFINITIONS[orbital.orbital].l;
    integerOptions(quantumN, 1, UI_MAX_N, n);
    integerOptions(quantumL, 0, n - 1, l);
    updateMOptions(l, orbital.basis === 'complex' ? orbital.m : 0);
    quantumL.disabled = orbital.basis === 'real';
    basisComplex.checked = orbital.basis === 'complex';
    basisReal.checked = orbital.basis === 'real';
    complexMField.hidden = orbital.basis !== 'complex';
    realOrbitalField.hidden = orbital.basis !== 'real';
    if (orbital.basis === 'real') realOrbital.value = orbital.orbital;
  }

  function renderRenderingControls(rendering: OrbitalRenderingState): void {
    requireElement('observableDensity', HTMLInputElement).checked =
      rendering.observable === 'density';
    requireElement('observablePhase', HTMLInputElement).checked = rendering.observable === 'phase';
    requireElement('displayMode', HTMLSelectElement).value = rendering.displayMode;
    isoThreshold.value = String(rendering.isoDensityFraction);
    requireElement('isoThresholdValue', HTMLOutputElement).textContent =
      `${Math.round(rendering.isoDensityFraction * 100)} %`;
    sampleCount.value = String(currentState.sampling.sampleCount);
    requireElement('sampleCountValue', HTMLOutputElement).textContent = formatInteger(
      currentState.sampling.sampleCount,
    );
    seedInput.value = String(currentState.sampling.seed >>> 0);
    pointOpacity.value = String(rendering.pointOpacity);
    requireElement('pointOpacityValue', HTMLOutputElement).textContent =
      `${Math.round(rendering.pointOpacity * 100)} %`;
    pointSize.value = String(rendering.pointSizePixels);
    requireElement('pointSizeValue', HTMLOutputElement).textContent =
      `${formatDecimal(rendering.pointSizePixels, 1)} px`;
    requireElement('axesToggle', HTMLInputElement).checked = rendering.showAxes;
    nodesToggle.checked = rendering.showNodes;
    requireElement('motionToggle', HTMLInputElement).checked = rendering.cameraRotationEnabled;
  }

  function emitComplexState(): void {
    const n = finiteInteger(quantumN.value, currentState.orbital.n, 1, UI_MAX_N);
    const l = finiteInteger(quantumL.value, 0, 0, n - 1);
    const m = finiteInteger(quantumM.value, 0, -l, l);
    requireCallbacks().setOrbital({ basis: 'complex', n, l, m });
  }

  function emitRealState(): void {
    const orbital = selectedRealOrbital(realOrbital);
    const definition = REAL_ORBITAL_DEFINITIONS[orbital];
    const n = Math.max(
      finiteInteger(quantumN.value, currentState.orbital.n, 1, UI_MAX_N),
      definition.minimumPrincipalQuantumNumber,
    );
    requireCallbacks().setOrbital({ basis: 'real', n, orbital });
  }

  function bindEvents(): void {
    const bound = requireCallbacks();
    basisComplex.addEventListener('change', () => {
      if (!basisComplex.checked) return;
      const n = finiteInteger(quantumN.value, 3, 1, UI_MAX_N);
      const l = finiteInteger(quantumL.value, 0, 0, n - 1);
      bound.setOrbital({ basis: 'complex', n, l, m: 0 });
    });
    basisReal.addEventListener('change', () => {
      if (basisReal.checked) emitRealState();
    });
    quantumN.addEventListener('change', () => {
      if (basisComplex.checked) emitComplexState();
      else emitRealState();
    });
    quantumL.addEventListener('change', () => {
      if (basisComplex.checked) emitComplexState();
    });
    quantumM.addEventListener('change', () => {
      if (basisComplex.checked) emitComplexState();
    });
    realOrbital.addEventListener('change', emitRealState);

    requireElement('observableDensity', HTMLInputElement).addEventListener('change', () => {
      bound.setRendering({ observable: 'density' });
    });
    requireElement('observablePhase', HTMLInputElement).addEventListener('change', () => {
      bound.setRendering({ observable: 'phase' });
    });
    requireElement('displayMode', HTMLSelectElement).addEventListener('change', (event) => {
      bound.setRendering({
        displayMode: (event.currentTarget as HTMLSelectElement)
          .value as OrbitalRenderingState['displayMode'],
      });
    });
    isoThreshold.addEventListener('input', () => {
      const value = Number(isoThreshold.value);
      requireElement('isoThresholdValue', HTMLOutputElement).textContent =
        `${Math.round(value * 100)} %`;
      bound.setRendering({ isoDensityFraction: value });
    });
    pointOpacity.addEventListener('input', () => {
      const value = Number(pointOpacity.value);
      requireElement('pointOpacityValue', HTMLOutputElement).textContent =
        `${Math.round(value * 100)} %`;
      bound.setRendering({ pointOpacity: value });
    });
    pointSize.addEventListener('input', () => {
      const value = Number(pointSize.value);
      requireElement('pointSizeValue', HTMLOutputElement).textContent =
        `${formatDecimal(value, 1)} px`;
      bound.setRendering({ pointSizePixels: value });
    });
    requireElement('axesToggle', HTMLInputElement).addEventListener('change', (event) => {
      bound.setRendering({ showAxes: (event.currentTarget as HTMLInputElement).checked });
    });
    nodesToggle.addEventListener('change', (event) => {
      bound.setRendering({ showNodes: (event.currentTarget as HTMLInputElement).checked });
    });
    requireElement('motionToggle', HTMLInputElement).addEventListener('change', (event) => {
      bound.setRendering({
        cameraRotationEnabled: (event.currentTarget as HTMLInputElement).checked,
      });
    });

    sampleCount.addEventListener('input', () => {
      requireElement('sampleCountValue', HTMLOutputElement).textContent = formatInteger(
        Number(sampleCount.value),
      );
    });
    sampleCount.addEventListener('change', () => {
      const value = finiteInteger(
        sampleCount.value,
        currentState.sampling.sampleCount,
        2000,
        60000,
      );
      bound.setSampling({ ...currentState.sampling, sampleCount: value });
      bound.generate();
    });
    seedInput.addEventListener('change', () => {
      const value = finiteInteger(seedInput.value, currentState.sampling.seed, 0, 0xffff_ffff);
      seedInput.value = String(value);
      bound.setSampling({ ...currentState.sampling, seed: value });
      bound.generate();
    });
    requireElement('randomizeSeed', HTMLButtonElement).addEventListener('click', () => {
      const value = randomUint32();
      seedInput.value = String(value);
      bound.setSampling({ ...currentState.sampling, seed: value });
      bound.generate();
    });
    requireElement('generateButton', HTMLButtonElement).addEventListener('click', () => {
      bound.generate();
    });
    requireElement('resetCamera', HTMLButtonElement).addEventListener('click', () => {
      bound.resetCamera();
    });
    requireElement('themeDark', HTMLButtonElement).addEventListener('click', () => {
      bound.setRendering({ theme: 'dark' });
    });
    requireElement('themeLight', HTMLButtonElement).addEventListener('click', () => {
      bound.setRendering({ theme: 'light' });
    });
    requireElement('settingsButton', HTMLButtonElement).addEventListener('click', () => {
      quantumN.focus();
    });

    const button2d = requireElement('btn2d', HTMLButtonElement);
    const toggle2d = (): void => {
      const next = !currentState.legacy.showLegacy2D;
      bound.setLegacy2DVisible(next);
      button2d.setAttribute('aria-expanded', String(next));
    };
    button2d.addEventListener('click', toggle2d);
    requireElement('close2d', HTMLButtonElement).addEventListener('click', toggle2d);
    window.addEventListener('keydown', (event) => {
      if (event.defaultPrevented || event.repeat) return;
      const target = event.target;
      if (
        target instanceof HTMLInputElement ||
        target instanceof HTMLSelectElement ||
        target instanceof HTMLTextAreaElement ||
        target instanceof HTMLButtonElement
      )
        return;
      if (event.key.toLowerCase() === 'q') toggle2d();
    });
  }

  const api: AppUi = {
    bind(nextCallbacks): void {
      callbacks = nextCallbacks;
      let storedTheme: string | null = null;
      try {
        storedTheme = localStorage.getItem(THEME_STORAGE_KEY);
      } catch {
        /* stockage indisponible */
      }
      const theme =
        storedTheme === 'light' || storedTheme === 'dark' ? storedTheme : themeFromDocument();
      if (theme !== currentState.rendering.theme) {
        // Réhydrate le modèle applicatif avant le premier rendu : le thème
        // persiste ainsi réellement au rechargement, pas uniquement dans le DOM.
        nextCallbacks.setRendering({ theme });
      } else {
        applyTheme(theme);
      }
      bindEvents();
    },
    finishGeneration(message = 'État prêt'): void {
      generationStatus.textContent = message;
      generationStatus.dataset.visible = 'false';
      if (progressTrack instanceof HTMLElement) progressTrack.dataset.visible = 'false';
      api.setEngineStatus('Moteur prêt', 'ready');
    },
    renderCharts(charts): void {
      renderRadialChart(requireElement('radialChart', SVGSVGElement), charts.radialDistribution);
      renderAngularChart(requireElement('angularChart', SVGSVGElement), charts.angularCut);
      requireElement('radialChartLabel', HTMLElement).textContent =
        `${formatDecimal(charts.radialDistribution.radiiBohr.at(-1) ?? 0, 1)} a₀`;
      requireElement('angularChartLabel', HTMLElement).textContent =
        `${charts.angularCut.planeLabel} · coupe géométrique de |Y(Ω)|²`;
      requireElement('angularPlaneLabel', HTMLElement).textContent = charts.angularCut.planeLabel;
    },
    renderState(state, presentation, nodesAvailable = true): void {
      currentState = state;
      applyTheme(state.rendering.theme);
      renderQuantumControls(state);
      renderRenderingControls(state.rendering);
      requireElement('quantumNValue', HTMLOutputElement).textContent = String(state.orbital.n);
      const degree = presentation.angularDegree;
      requireElement('quantumLValue', HTMLOutputElement).textContent = String(degree);
      requireElement('quantumMValue', HTMLOutputElement).textContent =
        state.orbital.basis === 'complex'
          ? String(state.orbital.m)
          : (presentation.realCombinationLabel?.match(/±\d+/u)?.[0] ?? '0');
      requireElement('orbitalNotation', HTMLElement).textContent = presentation.notation;
      requireElement('basisLabel', HTMLElement).textContent =
        state.orbital.basis === 'complex' ? 'complexe' : 'réelle';
      requireElement('orbitalSummary', HTMLElement).textContent =
        presentation.realCombinationLabel ??
        `${presentation.basisLabel} · projection m = ${state.orbital.basis === 'complex' ? state.orbital.m : 0}`;
      requireElement('energyValue', HTMLElement).innerHTML =
        `${formatDecimal(presentation.energy.electronVolts, 3)} <small>eV</small>`;
      requireElement('radiusValue', HTMLElement).innerHTML =
        `${formatDecimal(presentation.expectedRadius.bohr, 1)} <small>a₀</small>`;
      requireElement('scaleValue', HTMLElement).innerHTML =
        `${formatDecimal(presentation.expectedRadius.picometers, 1)} <small>pm</small>`;
      requireElement('modelLabel', HTMLElement).textContent = presentation.energy.modelLabel;
      requireElement('orbitalCaption', HTMLElement).textContent =
        `${presentation.notation} · ${state.rendering.observable === 'phase' ? 'phase de ψ' : '|ψ|² dV'}`;
      requireElement('radialNodesValue', HTMLElement).textContent = String(
        presentation.nodes.radialNodes,
      );
      requireElement('angularNodesValue', HTMLElement).textContent = String(
        presentation.nodes.angularNodes,
      );
      requireElement('totalNodesValue', HTMLElement).textContent = String(
        presentation.nodes.totalNodes,
      );
      nodesToggle.disabled = !nodesAvailable;
      nodesToggle.title = nodesAvailable
        ? 'Afficher les surfaces nodales ψ = 0'
        : 'Les nœuds de phase complexe non réelle ne sont pas affichés';
      const panel2d = requireElement('panel2d', HTMLElement);
      panel2d.hidden = !state.legacy.showLegacy2D;
      panel2d.classList.toggle('visible', state.legacy.showLegacy2D);
      requireElement('btn2d', HTMLButtonElement).setAttribute(
        'aria-expanded',
        String(state.legacy.showLegacy2D),
      );
    },
    setEngineStatus(message, status = 'ready'): void {
      const element = requireElement('engineStatus', HTMLElement);
      const textNode = element.lastChild;
      if (textNode) textNode.textContent = ` ${message}`;
      element.dataset.state = status;
    },
    showGeneration(message = 'Calcul de l’état orbital…'): void {
      generationStatus.textContent = message;
      generationStatus.dataset.visible = 'true';
      if (progressTrack instanceof HTMLElement) progressTrack.dataset.visible = 'true';
      api.setEngineStatus('Calcul…', 'busy');
    },
    updateFps(fps): void {
      requireElement('iFps', HTMLElement).textContent = String(fps);
    },
    updateHud(state, cameraDistance): void {
      requireElement('iPts', HTMLElement).textContent = formatInteger(state.sampling.sampleCount);
      requireElement('iOrb', HTMLElement).textContent =
        state.orbital.basis === 'real'
          ? `${state.orbital.n}${state.orbital.orbital}`
          : `${state.orbital.n},${state.orbital.l},${state.orbital.m}`;
      requireElement('iDist', HTMLElement).innerHTML =
        `${formatDecimal(cameraDistance, 1)} <small>a₀</small>`;
      requireElement('iSeed', HTMLElement).textContent =
        state.sampling.seed === 0x4154_4f4d ? 'ATOM' : String(state.sampling.seed >>> 0);
    },
    updateProgress(completed, total): void {
      const percent = total > 0 ? Math.round((completed / total) * 100) : 0;
      generationStatus.dataset.progress = String(percent);
      const bar = document.getElementById('generationProgress');
      if (bar instanceof HTMLElement) bar.style.width = `${percent}%`;
    },
  };

  return api;
}
