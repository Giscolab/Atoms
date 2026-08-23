import { SAMPLER_MAP } from '../data/orbitals';
import { parseWavefunctionData } from '../data/wavefunctionData';
import { legacyOrbitalColor } from '../rendering/orbitalColors';
import type { SceneRenderer } from '../rendering/sceneRenderer';
import { convertLegacyImportedPoint, sampleLegacyOrbitalPoint } from '../science/legacyScience';
import type { AppState } from '../state/appState';
import type { AppUi } from '../ui/appUi';

export interface CloudController {
  generate(callback?: () => void): void;
  importWavefunctionData(value: unknown): void;
}

export function createCloudController(
  state: AppState,
  renderer: SceneRenderer,
  ui: AppUi,
  refreshHud: () => void,
): CloudController {
  return {
    generate(callback): void {
      const { n, l, m, N } = state.quantum;
      ui.showLoading("Calcul de l'orbitale...");

      setTimeout(() => {
        state.cloud.positions = new Float32Array(N * 3);
        state.cloud.radii = new Float32Array(N);
        const colors = new Float32Array(N * 3);
        let done = 0;
        const chunk = 2000;

        function doChunk(): void {
          const end = Math.min(done + chunk, N);
          for (let i = done; i < end; i++) {
            const namedKey = SAMPLER_MAP[`${n}_${l}_${m}`];
            const { position, r, theta, phi } = sampleLegacyOrbitalPoint(n, l, m, namedKey);
            const [x, y, z] = position;
            state.cloud.positions[i * 3] = x;
            state.cloud.positions[i * 3 + 1] = y;
            state.cloud.positions[i * 3 + 2] = z;
            state.cloud.radii[i] = r;
            const [red, green, blue] = legacyOrbitalColor(
              r,
              theta,
              phi,
              n,
              l,
              m,
              state.useAltColor,
            );
            colors[i * 3] = red;
            colors[i * 3 + 1] = green;
            colors[i * 3 + 2] = blue;
          }
          done = end;
          ui.updateProgress(done, N);
          if (done < N) {
            setTimeout(doChunk, 0);
            return;
          }

          setTimeout(() => {
            renderer.replaceGeneratedCloud(state.cloud.positions.slice(), colors, state.pointSize);
            ui.hideLoading();
            refreshHud();
            if (callback) callback();
          }, 0);
        }
        doChunk();
      }, 50);
    },
    importWavefunctionData(value): void {
      const data = parseWavefunctionData(value);
      const N = data.points.length;
      state.cloud.positions = new Float32Array(N * 3);
      state.cloud.radii = new Float32Array(N);
      const colors = new Float32Array(N * 3);
      data.points.forEach((point, i) => {
        const { position, r, theta, phi } = convertLegacyImportedPoint(point);
        const [x, y, z] = position;
        state.cloud.positions[i * 3] = x;
        state.cloud.positions[i * 3 + 1] = y;
        state.cloud.positions[i * 3 + 2] = z;
        state.cloud.radii[i] = r;
        const [red, green, blue] = legacyOrbitalColor(
          r,
          theta,
          phi,
          state.quantum.n,
          state.quantum.l,
          state.quantum.m,
          state.useAltColor,
        );
        colors[i * 3] = red;
        colors[i * 3 + 1] = green;
        colors[i * 3 + 2] = blue;
      });
      renderer.replaceImportedCloud(state.cloud.positions.slice(), colors, state.pointSize);
      state.quantum.N = N;
      ui.updateParticleCount(N);
      refreshHud();
      ui.hideLoading();
    },
  };
}
