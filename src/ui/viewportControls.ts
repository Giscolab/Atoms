export interface ViewportCamera {
  getCameraDistance(): number;
  rotateCamera(azimuthDelta: number, elevationDelta: number): void;
  zoomCamera(distanceDelta: number): void;
}

export function bindViewportControls(
  canvas: HTMLCanvasElement,
  camera: ViewportCamera,
  updateDistance: (distance: number) => void,
): void {
  let dragging = false,
    lastX = 0,
    lastY = 0;

  canvas.addEventListener('mousedown', (event) => {
    dragging = true;
    lastX = event.clientX;
    lastY = event.clientY;
  });
  window.addEventListener('mouseup', () => (dragging = false));
  window.addEventListener('mousemove', (event) => {
    if (!dragging) return;
    camera.rotateCamera((event.clientX - lastX) * 0.008, (event.clientY - lastY) * 0.008);
    lastX = event.clientX;
    lastY = event.clientY;
    updateDistance(camera.getCameraDistance());
  });
  canvas.addEventListener(
    'wheel',
    (event) => {
      camera.zoomCamera(event.deltaY * 0.05);
      updateDistance(camera.getCameraDistance());
      event.preventDefault();
    },
    { passive: false },
  );

  let lastTouch: Touch | null = null;
  canvas.addEventListener('touchstart', (event) => {
    const touch = event.touches[0];
    if (touch) lastTouch = touch;
  });
  canvas.addEventListener(
    'touchmove',
    (event) => {
      if (!lastTouch) return;
      const touch = event.touches[0];
      if (!touch) return;
      camera.rotateCamera(
        (touch.clientX - lastTouch.clientX) * 0.01,
        (touch.clientY - lastTouch.clientY) * 0.01,
      );
      lastTouch = touch;
      event.preventDefault();
    },
    { passive: false },
  );
}
