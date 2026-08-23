import * as THREE from 'three';

// Phase 1 preserves the legacy r128 color pipeline for visual equivalence.
THREE.ColorManagement.enabled = false;

type CloudMesh = THREE.Points<THREE.BufferGeometry, THREE.PointsMaterial>;

export interface CloudPositionWriter {
  count: number;
  setXYZ(index: number, x: number, y: number, z: number): void;
  commit(): void;
}

export interface SceneRenderer {
  resize(viewport: HTMLElement): void;
  setClearColor(color: number): void;
  getCameraDistance(): number;
  rotateCamera(azimuthDelta: number, elevationDelta: number): void;
  zoomCamera(distanceDelta: number): void;
  replaceGeneratedCloud(positions: Float32Array, colors: Float32Array, pointSize: number): void;
  replaceImportedCloud(positions: Float32Array, colors: Float32Array, pointSize: number): void;
  setPointSize(pointSize: number): void;
  hasCloud(): boolean;
  getCloudPositionWriter(): CloudPositionWriter | null;
  rotateCloudY(delta: number): void;
  renderFrame(): void;
}

function requireCanvas2dContext(canvas: HTMLCanvasElement): CanvasRenderingContext2D {
  const context = canvas.getContext('2d');
  if (!context) throw new Error('Contexte Canvas 2D indisponible.');
  return context;
}

export function createSceneRenderer(canvas: HTMLCanvasElement): SceneRenderer {
  const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: false });
  renderer.outputColorSpace = THREE.LinearSRGBColorSpace;
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 1.5));
  renderer.setClearColor(0x050302, 1);

  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(60, 800 / 600, 0.01, 2000);
  const orbit = { az: 0.3, el: Math.PI / 2.5, r: 60 };

  function updateCamera(): void {
    const el = Math.max(0.02, Math.min(Math.PI - 0.02, orbit.el));
    camera.position.set(
      orbit.r * Math.sin(el) * Math.cos(orbit.az),
      orbit.r * Math.cos(el),
      orbit.r * Math.sin(el) * Math.sin(orbit.az),
    );
    camera.lookAt(0, 0, 0);
  }
  updateCamera();

  (function () {
    const N = 1200,
      pos = new Float32Array(N * 3),
      col = new Float32Array(N * 3);
    for (let i = 0; i < N; i++) {
      const th = Math.acos(2 * Math.random() - 1),
        ph = 2 * Math.PI * Math.random(),
        r = 400 + Math.random() * 200;
      pos[i * 3] = r * Math.sin(th) * Math.cos(ph);
      pos[i * 3 + 1] = r * Math.cos(th);
      pos[i * 3 + 2] = r * Math.sin(th) * Math.sin(ph);
      const b = 0.08 + Math.random() * 0.25;
      col[i * 3] = b * 0.9;
      col[i * 3 + 1] = b * 0.7;
      col[i * 3 + 2] = b * 0.3;
    }
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    geometry.setAttribute('color', new THREE.BufferAttribute(col, 3));
    scene.add(
      new THREE.Points(
        geometry,
        new THREE.PointsMaterial({ size: 0.4, vertexColors: true, sizeAttenuation: true }),
      ),
    );
  })();

  const nucleus = new THREE.Mesh(
    new THREE.SphereGeometry(0.5, 16, 16),
    new THREE.MeshBasicMaterial({ color: 0xff6030 }),
  );
  scene.add(nucleus);
  const ring = new THREE.Mesh(
    new THREE.RingGeometry(0.7, 1.0, 32),
    new THREE.MeshBasicMaterial({
      color: 0xff4020,
      side: THREE.DoubleSide,
      transparent: true,
      opacity: 0.3,
    }),
  );
  scene.add(ring);

  const dotTexture = (() => {
    const dotCanvas = document.createElement('canvas');
    dotCanvas.width = dotCanvas.height = 16;
    const context = requireCanvas2dContext(dotCanvas);
    const gradient = context.createRadialGradient(8, 8, 0, 8, 8, 8);
    gradient.addColorStop(0, 'rgba(255,255,255,1)');
    gradient.addColorStop(0.5, 'rgba(255,255,255,0.8)');
    gradient.addColorStop(1, 'rgba(255,255,255,0)');
    context.fillStyle = gradient;
    context.fillRect(0, 0, 16, 16);
    return new THREE.CanvasTexture(dotCanvas);
  })();

  let cloudMesh: CloudMesh | null = null;

  function removeCloud(): void {
    if (!cloudMesh) return;
    scene.remove(cloudMesh);
    cloudMesh.geometry.dispose();
    cloudMesh.material.dispose();
    cloudMesh = null;
  }

  function replaceGeneratedCloud(
    positions: Float32Array,
    colors: Float32Array,
    pointSize: number,
  ): void {
    removeCloud();
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
    cloudMesh = new THREE.Points(
      geometry,
      new THREE.PointsMaterial({
        size: pointSize,
        vertexColors: true,
        transparent: true,
        opacity: 0.88,
        sizeAttenuation: true,
        depthWrite: false,
        map: dotTexture,
        alphaTest: 0.01,
      }),
    );
    scene.add(cloudMesh);
  }

  function replaceImportedCloud(
    positions: Float32Array,
    colors: Float32Array,
    pointSize: number,
  ): void {
    removeCloud();
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
    cloudMesh = new THREE.Points(
      geometry,
      new THREE.PointsMaterial({
        size: pointSize,
        vertexColors: true,
        transparent: true,
        opacity: 0.88,
        sizeAttenuation: true,
        depthWrite: false,
      }),
    );
    scene.add(cloudMesh);
  }

  return {
    resize(viewport): void {
      const width = viewport.clientWidth - 40;
      const height = viewport.clientHeight - 40;
      canvas.width = width;
      canvas.height = height;
      camera.aspect = width / height;
      camera.updateProjectionMatrix();
      renderer.setSize(width, height);
    },
    setClearColor(color): void {
      renderer.setClearColor(color, 1);
    },
    getCameraDistance(): number {
      return orbit.r;
    },
    rotateCamera(azimuthDelta, elevationDelta): void {
      orbit.az += azimuthDelta;
      orbit.el -= elevationDelta;
      updateCamera();
    },
    zoomCamera(distanceDelta): void {
      orbit.r = Math.max(5, Math.min(400, orbit.r + distanceDelta));
      updateCamera();
    },
    replaceGeneratedCloud,
    replaceImportedCloud,
    setPointSize(pointSize): void {
      if (cloudMesh) cloudMesh.material.size = pointSize;
    },
    hasCloud(): boolean {
      return cloudMesh !== null;
    },
    getCloudPositionWriter(): CloudPositionWriter | null {
      if (!cloudMesh) return null;
      const position = cloudMesh.geometry.getAttribute('position');
      if (!(position instanceof THREE.BufferAttribute)) {
        throw new Error('Attribut de position du nuage indisponible.');
      }
      return {
        count: position.count,
        setXYZ(index, x, y, z): void {
          position.setXYZ(index, x, y, z);
        },
        commit(): void {
          position.needsUpdate = true;
        },
      };
    },
    rotateCloudY(delta): void {
      if (cloudMesh) cloudMesh.rotation.y += delta;
    },
    renderFrame(): void {
      nucleus.rotation.y += 0.01;
      renderer.render(scene, camera);
    },
  };
}
