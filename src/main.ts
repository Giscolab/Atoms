import { createAnimationLoop } from './app/animationLoop';
import { createCloudController } from './app/cloudController';
import { createPhotoelectric2D } from './rendering/photoelectric2d';
import { createSceneRenderer } from './rendering/sceneRenderer';
import { LEGACY_PHOTOELECTRIC_ENERGY, legacyEnergyEv } from './science/legacyScience';
import { createAppState } from './state/appState';
import { createAppUi } from './ui/appUi';
import { requireElement } from './ui/dom';
import { bindViewportControls } from './ui/viewportControls';
import { readWavefunctionFile } from './ui/wavefunctionFile';

const state = createAppState();
const canvas3d = requireElement('atomSimCanvas', HTMLCanvasElement);
const renderer = createSceneRenderer(canvas3d);
const ui = createAppUi(state);

const refreshHud = (): void => {
  ui.updateHud(legacyEnergyEv(state.quantum.n), renderer.getCameraDistance());
};

const cloudController = createCloudController(state, renderer, ui, refreshHud);

function loadWavefunction(file: File): Promise<void> {
  return new Promise<void>((resolve, reject) => {
    readWavefunctionFile(
      file,
      (value) => {
        cloudController.importWavefunctionData(value);
        resolve();
      },
      reject,
    );
  });
}

bindViewportControls(canvas3d, renderer, (distance) => {
  ui.updateDistance(distance);
});

ui.bind({
  generateCloud: () => {
    cloudController.generate();
  },
  refreshHud,
  setPointSize: (pointSize) => {
    renderer.setPointSize(pointSize);
  },
  setThemeColor: (color) => {
    renderer.setClearColor(color);
  },
  loadWavefunction,
  initialize2D: () => {
    createPhotoelectric2D(
      requireElement('c2d', HTMLCanvasElement),
      () => state.show2D,
      LEGACY_PHOTOELECTRIC_ENERGY,
    ).init();
  },
});

const animationLoop = createAnimationLoop(state, renderer, ui);

window.addEventListener('resize', () => {
  renderer.resize(requireElement('viewport', HTMLElement));
});

renderer.resize(requireElement('viewport', HTMLElement));
refreshHud();
cloudController.generate(() => {
  animationLoop.start();
});
