import type { SceneRenderer } from '../rendering/sceneRenderer';
import { advanceLegacyProbabilityFlowPoint } from '../science/legacyScience';
import type { AppState } from '../state/appState';
import type { AppUi } from '../ui/appUi';

export interface AnimationLoop {
  start(): void;
}

function requireNumber(values: ArrayLike<number>, index: number): number {
  const value = values[index];
  if (value === undefined) {
    throw new RangeError(`Indice numérique hors limites : ${index}.`);
  }
  return value;
}

export function createAnimationLoop(
  state: AppState,
  renderer: SceneRenderer,
  ui: AppUi,
): AnimationLoop {
  function animateFlow(): void {
    if (!renderer.hasCloud() || !state.animating || state.quantum.m === 0) return;
    const positionWriter = renderer.getCloudPositionWriter();
    if (!positionWriter) return;
    const m = state.quantum.m;
    const availablePointCount = Math.min(
      state.quantum.N,
      positionWriter.count,
      state.cloud.radii.length,
      Math.floor(state.cloud.positions.length / 3),
    );
    for (let i = 0; i < availablePointCount; i++) {
      const px = requireNumber(state.cloud.positions, i * 3),
        py = requireNumber(state.cloud.positions, i * 3 + 1),
        pz = requireNumber(state.cloud.positions, i * 3 + 2);
      const [x, y, z] = advanceLegacyProbabilityFlowPoint(
        px,
        py,
        pz,
        requireNumber(state.cloud.radii, i),
        m,
      );
      state.cloud.positions[i * 3] = x;
      state.cloud.positions[i * 3 + 1] = y;
      state.cloud.positions[i * 3 + 2] = z;
      positionWriter.setXYZ(
        i,
        requireNumber(state.cloud.positions, i * 3),
        requireNumber(state.cloud.positions, i * 3 + 1),
        requireNumber(state.cloud.positions, i * 3 + 2),
      );
    }
    positionWriter.commit();
  }

  let fpsTimes: number[] = [],
    fpsLast = performance.now();

  function tickFPS(): void {
    const now = performance.now();
    fpsTimes.push(now);
    fpsTimes = fpsTimes.filter((time) => now - time < 1000);
    if (now - fpsLast > 400) {
      ui.updateFps(fpsTimes.length);
      fpsLast = now;
    }
  }

  let frame = 0;

  function loop3D(): void {
    requestAnimationFrame(loop3D);
    frame++;
    tickFPS();
    if (renderer.hasCloud()) {
      if (state.animating && state.quantum.m !== 0 && frame % 2 === 0) animateFlow();
      if (state.quantum.m === 0 && state.animating) renderer.rotateCloudY(0.002);
    }
    renderer.renderFrame();
  }

  return { start: loop3D };
}
