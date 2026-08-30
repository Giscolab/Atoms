export interface ViewportCamera {
  fitCameraToOrbital?(): void;
  getCameraDistance(): number;
  rotateCamera(azimuthDelta: number, elevationDelta: number): void;
  zoomCamera(distanceDelta: number): void;
}

export function bindViewportControls(
  canvas: HTMLCanvasElement,
  camera: ViewportCamera,
  updateDistance: (distance: number) => void,
  onReset?: () => void,
): void {
  let pointerId: number | null = null;
  let lastX = 0;
  let lastY = 0;

  const syncDistance = (): void => {
    updateDistance(camera.getCameraDistance());
  };

  canvas.addEventListener('pointerdown', (event) => {
    if (pointerId !== null) return;
    pointerId = event.pointerId;
    lastX = event.clientX;
    lastY = event.clientY;
    canvas.setPointerCapture(event.pointerId);
    canvas.focus({ preventScroll: true });
  });
  canvas.addEventListener('pointermove', (event) => {
    if (pointerId !== event.pointerId) return;
    camera.rotateCamera((event.clientX - lastX) * 0.008, (event.clientY - lastY) * 0.008);
    lastX = event.clientX;
    lastY = event.clientY;
    syncDistance();
  });
  const releasePointer = (event: PointerEvent): void => {
    if (pointerId !== event.pointerId) return;
    pointerId = null;
    if (canvas.hasPointerCapture(event.pointerId)) canvas.releasePointerCapture(event.pointerId);
  };
  canvas.addEventListener('pointerup', releasePointer);
  canvas.addEventListener('pointercancel', releasePointer);
  canvas.addEventListener('lostpointercapture', () => {
    pointerId = null;
  });

  canvas.addEventListener(
    'wheel',
    (event) => {
      camera.zoomCamera(event.deltaY * 0.045);
      syncDistance();
      event.preventDefault();
    },
    { passive: false },
  );

  canvas.addEventListener('keydown', (event) => {
    const step = event.shiftKey ? 0.14 : 0.07;
    switch (event.key) {
      case 'ArrowLeft':
        camera.rotateCamera(-step, 0);
        break;
      case 'ArrowRight':
        camera.rotateCamera(step, 0);
        break;
      case 'ArrowUp':
        camera.rotateCamera(0, -step);
        break;
      case 'ArrowDown':
        camera.rotateCamera(0, step);
        break;
      case '+':
      case '=':
        camera.zoomCamera(-1);
        break;
      case '-':
      case '_':
        camera.zoomCamera(1);
        break;
      case '0':
        camera.fitCameraToOrbital?.();
        onReset?.();
        break;
      default:
        return;
    }
    syncDistance();
    event.preventDefault();
  });
}
