import {
  createLegacyAtomRing,
  createLegacyIncidentWaves,
  createLegacyRadialWaves,
  getLegacyElectronPosition,
  getLegacyWaveDisplayPoint,
  updateLegacyAtom,
  updateLegacyWave,
  type LegacyAtom2D,
  type LegacyWave,
} from './legacyPhotoelectricModel';

export interface Photoelectric2D {
  init(): void;
}

function requireCanvas2dContext(canvas: HTMLCanvasElement): CanvasRenderingContext2D {
  const context = canvas.getContext('2d');
  if (!context) throw new Error('Contexte Canvas 2D indisponible.');
  return context;
}

export function createPhotoelectric2D(
  canvas: HTMLCanvasElement,
  isVisible: () => boolean,
  energy: number,
): Photoelectric2D {
  return {
    init(): void {
      const context = requireCanvas2dContext(canvas);
      const width = 300,
        height = 195,
        orbitDistance = 15;

      function drawWave(wave: LegacyWave): void {
        context.strokeStyle = `rgba(${~~(wave.col[0] * 255)},${~~(wave.col[1] * 255)},${~~(wave.col[2] * 255)},.9)`;
        context.lineWidth = 0.7;
        context.beginPath();
        let first = true;
        for (const point of wave.points) {
          const displayPoint = getLegacyWaveDisplayPoint(wave, point);
          if (first) {
            context.moveTo(displayPoint.x + width / 2, displayPoint.y + height / 2);
            first = false;
          } else context.lineTo(displayPoint.x + width / 2, displayPoint.y + height / 2);
        }
        context.stroke();
      }

      function drawAtom(atom: LegacyAtom2D): void {
        context.strokeStyle = 'rgba(255,154,60,.18)';
        context.lineWidth = 0.5;
        context.beginPath();
        context.arc(
          atom.pos.x + width / 2,
          atom.pos.y + height / 2,
          atom.n * orbitDistance,
          0,
          Math.PI * 2,
        );
        context.stroke();
        const electron = getLegacyElectronPosition(atom, orbitDistance);
        context.fillStyle = '#ff9a3c';
        context.beginPath();
        context.arc(electron.x + width / 2, electron.y + height / 2, 1.5, 0, Math.PI * 2);
        context.fill();
        context.fillStyle = '#ff4020';
        context.beginPath();
        context.arc(atom.pos.x + width / 2, atom.pos.y + height / 2, 3, 0, Math.PI * 2);
        context.fill();
      }

      const atoms2d = createLegacyAtomRing(12, 55);
      const waves2d = createLegacyIncidentWaves(energy, 8, 80, -76, 22);

      canvas.addEventListener('click', (event) => {
        const rx = event.offsetX - width / 2,
          ry = event.offsetY - height / 2;
        waves2d.push(...createLegacyRadialWaves(energy, { x: rx, y: ry }, 18));
      });

      function loop2D(): void {
        if (!isVisible()) {
          requestAnimationFrame(loop2D);
          return;
        }
        context.fillStyle = 'rgba(5,3,2,.85)';
        context.fillRect(0, 0, width, height);
        for (const atom of atoms2d) {
          drawAtom(atom);
          updateLegacyAtom(atom);
        }
        for (let i = 0; i < waves2d.length;) {
          const wave = waves2d[i];
          if (!wave) break;
          if (wave.energy === 0) {
            i++;
            continue;
          }
          drawWave(wave);
          if (updateLegacyWave(wave, 0.016, width, height)) waves2d.splice(i, 1);
          else i++;
        }
        requestAnimationFrame(loop2D);
      }
      loop2D();
    },
  };
}
