import type { SceneRenderer } from '../rendering/renderingContracts';

export interface AnimationLoop {
  start(): void;
  stop(): void;
}

export interface AnimationLoopOptions {
  readonly autoRotate: () => boolean;
  readonly onFps?: (fps: number) => void;
}

/**
 * Boucle de présentation : elle ne modifie jamais les échantillons. Le seul
 * mouvement automatique est une rotation de caméra, clairement distincte de
 * toute dynamique physique de l'état stationnaire.
 */
export function createAnimationLoop(
  renderer: SceneRenderer,
  options: AnimationLoopOptions,
): AnimationLoop {
  let frameHandle: number | null = null;
  let running = false;
  let lastTime = 0;
  let recentFrames: number[] = [];
  let lastFpsReport = 0;

  function tick(time: number): void {
    if (!running) return;
    frameHandle = requestAnimationFrame(tick);
    const deltaMilliseconds = lastTime === 0 ? 16.7 : Math.min(100, Math.max(0, time - lastTime));
    lastTime = time;
    if (options.autoRotate()) renderer.rotateCameraAutomatically(deltaMilliseconds * 0.000045);
    renderer.renderFrame();

    recentFrames.push(time);
    recentFrames = recentFrames.filter((frameTime) => time - frameTime < 1000);
    if (time - lastFpsReport >= 500) {
      options.onFps?.(recentFrames.length);
      lastFpsReport = time;
    }
  }

  return {
    start(): void {
      if (running) return;
      running = true;
      lastTime = 0;
      lastFpsReport = performance.now();
      recentFrames = [];
      frameHandle = requestAnimationFrame(tick);
    },
    stop(): void {
      running = false;
      if (frameHandle !== null) cancelAnimationFrame(frameHandle);
      frameHandle = null;
    },
  };
}
