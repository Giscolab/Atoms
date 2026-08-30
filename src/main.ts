import { createAnimationLoop } from './app/animationLoop';
import { buildOrbitalPresentation } from './app/orbitalPresentation';
import {
  createOrbitalWorkerClient,
  OrbitalWorkerCancellationError,
} from './app/orbitalWorkerClient';
import { createPhotoelectric2D } from './rendering/photoelectric2d';
import { createSceneRenderer } from './rendering/sceneRenderer';
import { LEGACY_PHOTOELECTRIC_ENERGY } from './science/legacyScience';
import { createAppState, normalizeAppState, type AppState } from './state/appState';
import { createAppUi } from './ui/appUi';
import { requireElement } from './ui/dom';
import { bindViewportControls } from './ui/viewportControls';

/**
 * Gestion contextuelle du scroll.
 * - Mobile (≤ 900px) : scroll global autorisé, restauration automatique conservée.
 * - Desktop (> 900px) : scroll forcé à 0, y compris lors du passage mobile → desktop.
 */
const compactLayout = window.matchMedia('(max-width: 900px)');

function resetDocumentScrollForDesktop(): void {
  if (compactLayout.matches) {
    // En mode mobile, on restaure le comportement natif du navigateur.
    if ('scrollRestoration' in history) {
      history.scrollRestoration = 'auto';
    }
    return;
  }

  if ('scrollRestoration' in history) {
    history.scrollRestoration = 'manual';
  }

  window.scrollTo({ top: 0, left: 0, behavior: 'auto' });
}

resetDocumentScrollForDesktop();
compactLayout.addEventListener('change', resetDocumentScrollForDesktop);

let state: AppState = createAppState();
if (window.matchMedia('(prefers-reduced-motion: reduce)').matches) {
  state = normalizeAppState({
    ...state,
    rendering: { ...state.rendering, cameraRotationEnabled: false },
  });
}

const canvas3d = requireElement('atomSimCanvas', HTMLCanvasElement);
const viewport = requireElement('viewport', HTMLElement);
const renderer = createSceneRenderer(canvas3d);
const ui = createAppUi(state);
const workerClient = createOrbitalWorkerClient();

// État runtime impératif : le canvas 2D a-t-il été réellement créé ?
let photoelectricInitialized = false;
// Variables partagées (contexte unique)
let lastNodesAvailable = state.orbital.basis === 'real' || state.orbital.m === 0;
let generationSequence = 0;

function sameOrbitalState(left: AppState['orbital'], right: AppState['orbital']): boolean {
  if (left.basis !== right.basis || left.n !== right.n) return false;
  if (left.basis === 'complex' && right.basis === 'complex') {
    return left.l === right.l && left.m === right.m;
  }
  return left.basis === 'real' && right.basis === 'real' && left.orbital === right.orbital;
}

function renderState(): void {
  ui.renderState(state, buildOrbitalPresentation(state.orbital), lastNodesAvailable);
  ui.updateHud(state, renderer.getCameraDistance());
}

function updateState(next: AppState): void {
  state = normalizeAppState(next);
  lastNodesAvailable = state.orbital.basis === 'real' || state.orbital.m === 0;
  renderer.setAppearance(state.rendering);
  renderState(); // Inclut déjà ui.updateHud()
}

function initializeLegacy2DIfNeeded(): void {
  if (photoelectricInitialized) return;

  createPhotoelectric2D(
    requireElement('c2d', HTMLCanvasElement),
    () => state.legacy.showLegacy2D,
    LEGACY_PHOTOELECTRIC_ENERGY,
  ).init();

  photoelectricInitialized = true;

  if (!state.legacy.legacy2DInitialized) {
    updateState(
      normalizeAppState({
        ...state,
        legacy: {
          ...state.legacy,
          legacy2DInitialized: true,
        },
      }),
    );
  }
}

function setLegacy2DVisible(visible: boolean): void {
  updateState(
    normalizeAppState({
      ...state,
      legacy: {
        ...state.legacy,
        showLegacy2D: visible,
      },
    }),
  );

  if (visible) {
    initializeLegacy2DIfNeeded();
  }
}

function showGenerationError(message: string): void {
  ui.setEngineStatus('Erreur de calcul', 'error');
  ui.showGeneration(`Erreur : ${message}`);
  const status = requireElement('generationStatus', HTMLElement);
  status.dataset.visible = 'true';
}

async function generateCurrentOrbital(): Promise<void> {
  const requestState = state;
  const jobId = `orbital-${++generationSequence}`;
  ui.showGeneration('Échantillonnage de |ψ|²…');

  try {
    const payload = await workerClient.generate(
      {
        jobId,
        sampleCount: requestState.sampling.sampleCount,
        seed: requestState.sampling.seed,
        state: requestState.orbital,
      },
      (progress) => {
        const labels: Record<string, string> = {
          charts: 'Calcul des distributions…',
          field: 'Calcul du champ volumique…',
          sampling: 'Échantillonnage de |ψ|²…',
          transfer: 'Transfert vers le rendu…',
        };
        ui.showGeneration(labels[progress.stage] ?? 'Calcul orbital…');
        ui.updateProgress(progress.completed, progress.total);
      },
    );

    // Résultat obsolète si l'état a changé pendant le calcul.
    if (
      !sameOrbitalState(state.orbital, requestState.orbital) ||
      state.sampling.sampleCount !== requestState.sampling.sampleCount ||
      state.sampling.seed !== requestState.sampling.seed
    ) {
      return;
    }

    renderer.setOrbital({
      field: payload.field,
      samples: {
        phaseRadians: payload.sampleSet.phaseRadians,
        positionsBohr: payload.sampleSet.positionsBohr,
      },
    });

    lastNodesAvailable = payload.field.nodesAvailable;
    renderer.setAppearance(state.rendering);
    ui.renderState(state, buildOrbitalPresentation(state.orbital), lastNodesAvailable);
    ui.renderCharts(payload.charts);
    ui.updateHud(state, renderer.getCameraDistance());
    ui.finishGeneration('État prêt');
  } catch (error) {
    if (error instanceof OrbitalWorkerCancellationError) return;
    showGenerationError(error instanceof Error ? error.message : String(error));
  }
}

const animationLoop = createAnimationLoop(renderer, {
  autoRotate: () => state.rendering.cameraRotationEnabled,
  onFps: (fps) => ui.updateFps(fps),
});

ui.bind({
  generate: () => {
    void generateCurrentOrbital();
  },
  resetCamera: () => {
    renderer.fitCameraToOrbital();
    ui.updateHud(state, renderer.getCameraDistance());
  },
  setLegacy2DVisible,
  setOrbital: (orbital) => {
    updateState(normalizeAppState({ ...state, orbital }));
    void generateCurrentOrbital();
  },
  setRendering: (patch) => {
    updateState(normalizeAppState({ ...state, rendering: { ...state.rendering, ...patch } }));
  },
  setSampling: (sampling) => {
    updateState(normalizeAppState({ ...state, sampling }));
  },
});

bindViewportControls(canvas3d, renderer, (distance) => {
  ui.updateHud(state, distance);
});

function handleResize(): void {
  renderer.resize(viewport);
}
renderer.resize(viewport);
window.addEventListener('resize', handleResize);
if ('ResizeObserver' in window) {
  new ResizeObserver(handleResize).observe(viewport);
}

window.addEventListener('beforeunload', () => {
  animationLoop.stop();
  workerClient.dispose();
  renderer.dispose();
});

renderState();
animationLoop.start();
void generateCurrentOrbital();