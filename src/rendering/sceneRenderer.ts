import * as THREE from 'three';
import { MarchingCubes } from 'three/examples/jsm/objects/MarchingCubes.js';
import { DENSITY_COLOR_SRGB, phaseColorSrgb, type SrgbColor } from './phasePalette';
import type {
  OrbitalAppearance,
  OrbitalRenderDataset,
  RenderTheme,
  SceneDiagnostics,
  SceneRenderer,
} from './renderingContracts';

const DEFAULT_APPEARANCE: OrbitalAppearance = {
  displayMode: 'hybrid',
  isoDensityFraction: 0.2,
  observable: 'density',
  pointOpacity: 0.82,
  pointSizePixels: 1.5,
  showAxes: true,
  showNodes: false,
  theme: 'dark',
};

const CAMERA_FIELD_OF_VIEW_DEGREES = 48;
const CAMERA_ELEVATION_LIMIT_RADIANS = Math.PI / 2 - 0.02;
const CAMERA_FIT_MARGIN = 1.16;
const MINIMUM_RENDER_DIMENSION_PIXELS = 1;
const MAXIMUM_DEVICE_PIXEL_RATIO = 2;
const GRID_DIVISION_COUNT = 12;

type OrbitalCloud = THREE.Points<THREE.BufferGeometry, THREE.PointsMaterial>;

interface OrbitCameraState {
  azimuthRadians: number;
  distanceBohr: number;
  elevationRadians: number;
}

function assertFinitePositive(value: number, label: string): void {
  if (!Number.isFinite(value) || value <= 0) {
    throw new RangeError(`${label} doit être fini et strictement positif.`);
  }
}

function validateAppearance(appearance: OrbitalAppearance): void {
  if (
    !Number.isFinite(appearance.isoDensityFraction) ||
    appearance.isoDensityFraction <= 0 ||
    appearance.isoDensityFraction >= 1
  ) {
    throw new RangeError("Le seuil d'isodensité doit appartenir à l'intervalle ]0,1[.");
  }
  assertFinitePositive(appearance.pointSizePixels, 'La taille des points en pixels');
  if (
    !Number.isFinite(appearance.pointOpacity) ||
    appearance.pointOpacity < 0 ||
    appearance.pointOpacity > 1
  ) {
    throw new RangeError("L'opacité des points doit appartenir à l'intervalle [0,1].");
  }
}

function validateDataset(dataset: OrbitalRenderDataset): void {
  const { field, samples } = dataset;
  assertFinitePositive(field.extentBohr, "L'étendue de la grille en a₀");
  if (!Number.isSafeInteger(field.resolution) || field.resolution < 4) {
    throw new RangeError(
      'La résolution de la grille doit être un entier sûr supérieur ou égal à 4.',
    );
  }

  const fieldValueCount = field.resolution ** 3;
  if (
    field.densityNormalized.length !== fieldValueCount ||
    field.phaseRadians.length !== fieldValueCount ||
    field.signedAmplitudeNormalized.length !== fieldValueCount
  ) {
    throw new RangeError('Les buffers du champ ne correspondent pas à sa résolution cubique.');
  }
  if (samples.positionsBohr.length % 3 !== 0) {
    throw new RangeError('Le buffer de positions doit contenir des triplets cartésiens.');
  }
  if (samples.phaseRadians.length * 3 !== samples.positionsBohr.length) {
    throw new RangeError('Le buffer de phase doit contenir une valeur par position cartésienne.');
  }
}

function materialList(material: THREE.Material | THREE.Material[]): THREE.Material[] {
  return material instanceof THREE.Material ? [material] : material;
}

function disposeMaterials(material: THREE.Material | THREE.Material[]): void {
  for (const item of materialList(material)) item.dispose();
}

function createDotTexture(): THREE.CanvasTexture {
  const dotCanvas = document.createElement('canvas');
  dotCanvas.width = 32;
  dotCanvas.height = 32;
  const context = dotCanvas.getContext('2d');
  if (!context) throw new Error('Contexte Canvas 2D indisponible pour la texture des points.');

  const gradient = context.createRadialGradient(16, 16, 0, 16, 16, 16);
  gradient.addColorStop(0, 'rgba(255,255,255,1)');
  gradient.addColorStop(0.5, 'rgba(255,255,255,0.88)');
  gradient.addColorStop(1, 'rgba(255,255,255,0)');
  context.fillStyle = gradient;
  context.fillRect(0, 0, 32, 32);

  const texture = new THREE.CanvasTexture(dotCanvas);
  texture.colorSpace = THREE.SRGBColorSpace;
  texture.needsUpdate = true;
  return texture;
}

function writeWorkingColor(
  target: Float32Array,
  offset: number,
  colorSrgb: SrgbColor,
  scratchColor: THREE.Color,
): void {
  scratchColor.setRGB(colorSrgb[0], colorSrgb[1], colorSrgb[2], THREE.SRGBColorSpace);
  target[offset] = scratchColor.r;
  target[offset + 1] = scratchColor.g;
  target[offset + 2] = scratchColor.b;
}

function colorBufferForPhases(
  phasesRadians: Float32Array,
  observable: OrbitalAppearance['observable'],
): Float32Array {
  const colors = new Float32Array(phasesRadians.length * 3);
  const scratchColor = new THREE.Color();
  for (let index = 0; index < phasesRadians.length; index += 1) {
    const phase = phasesRadians[index];
    const color = observable === 'phase' ? phaseColorSrgb(phase ?? Number.NaN) : DENSITY_COLOR_SRGB;
    writeWorkingColor(colors, index * 3, color, scratchColor);
  }
  return colors;
}

function createPointsMaterial(
  appearance: OrbitalAppearance,
  dotTexture: THREE.Texture,
): THREE.PointsMaterial {
  return new THREE.PointsMaterial({
    alphaTest: 0.015,
    depthWrite: false,
    map: dotTexture,
    opacity: appearance.pointOpacity,
    size: appearance.pointSizePixels,
    // Le contrôle est explicitement exprimé en pixels, indépendamment de la distance caméra.
    sizeAttenuation: false,
    transparent: appearance.pointOpacity < 1,
    vertexColors: true,
  });
}

function maximumMarchingCubeTriangleCount(resolution: number): number {
  // MarchingCubes parcourt (N-3)^3 cellules et une cellule produit au plus 5 triangles.
  return Math.max(1, 5 * (resolution - 3) ** 3);
}

function applyPhysicalGridTransform(surface: MarchingCubes, extentBohr: number): void {
  const resolution = surface.resolution;
  /*
   * MarchingCubes place l'indice i en 2i/N-1, tandis que le champ échantillonne
   * -extent + 2*extent*i/(N-1). Cette transformation affine fait coïncider
   * exactement les sommets avec les coordonnées cartésiennes en a₀.
   */
  const scaleBohr = (extentBohr * resolution) / (resolution - 1);
  const offsetBohr = extentBohr / (resolution - 1);
  surface.scale.setScalar(scaleBohr);
  surface.position.setScalar(offsetBohr);
}

function orbitalRadiusBohr(dataset: OrbitalRenderDataset): number {
  const positions = dataset.samples.positionsBohr;
  let maximumRadius = dataset.field.extentBohr;
  for (let index = 0; index < positions.length; index += 3) {
    const x = positions[index];
    const y = positions[index + 1];
    const z = positions[index + 2];
    if (x === undefined || y === undefined || z === undefined) continue;
    maximumRadius = Math.max(maximumRadius, Math.hypot(x, y, z));
  }
  return maximumRadius;
}

function themeColors(theme: RenderTheme): {
  readonly background: THREE.ColorRepresentation;
  readonly gridCenter: THREE.ColorRepresentation;
  readonly gridLines: THREE.ColorRepresentation;
  readonly node: THREE.ColorRepresentation;
  readonly nucleus: THREE.ColorRepresentation;
} {
  return theme === 'dark'
    ? {
        background: 0x07131a,
        gridCenter: 0x4f6973,
        gridLines: 0x233943,
        node: 0xe8edf0,
        nucleus: 0xf4f7f8,
      }
    : {
        background: 0xf2f7f8,
        gridCenter: 0x71868d,
        gridLines: 0xb9c8cc,
        node: 0x33464d,
        nucleus: 0x24373e,
      };
}

export function createSceneRenderer(canvas: HTMLCanvasElement): SceneRenderer {
  const renderer = new THREE.WebGLRenderer({
    antialias: true,
    canvas,
    powerPreference: 'high-performance',
  });
  renderer.outputColorSpace = THREE.SRGBColorSpace;
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, MAXIMUM_DEVICE_PIXEL_RATIO));

  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(CAMERA_FIELD_OF_VIEW_DEGREES, 1, 0.01, 1000);
  // Le repère scientifique conserve z comme axe polaire vertical.
  camera.up.set(0, 0, 1);
  const orbit: OrbitCameraState = {
    azimuthRadians: 0.65,
    distanceBohr: 20,
    elevationRadians: 0.42,
  };

  const ambientLight = new THREE.HemisphereLight(0xe8f7ff, 0x18252b, 1.35);
  const directionalLight = new THREE.DirectionalLight(0xffffff, 1.7);
  directionalLight.position.set(4, -3, 6);
  scene.add(ambientLight, directionalLight);

  const nucleusGeometry = new THREE.SphereGeometry(1, 20, 14);
  const nucleusMaterial = new THREE.MeshStandardMaterial({ roughness: 0.25 });
  const schematicNucleus = new THREE.Mesh(nucleusGeometry, nucleusMaterial);
  schematicNucleus.name = 'Noyau schématique — non représenté à l’échelle';
  scene.add(schematicNucleus);

  const dotTexture = createDotTexture();
  let appearance: OrbitalAppearance = { ...DEFAULT_APPEARANCE };
  let dataset: OrbitalRenderDataset | null = null;
  let cloud: OrbitalCloud | null = null;
  let densitySurface: MarchingCubes | null = null;
  let nodeSurface: MarchingCubes | null = null;
  let guides: THREE.Group | null = null;
  let guideGeometries: THREE.BufferGeometry[] = [];
  let guideMaterials: THREE.Material[] = [];
  let currentOrbitalRadiusBohr = 10;

  function updateCamera(): void {
    orbit.elevationRadians = THREE.MathUtils.clamp(
      orbit.elevationRadians,
      -CAMERA_ELEVATION_LIMIT_RADIANS,
      CAMERA_ELEVATION_LIMIT_RADIANS,
    );
    const equatorialDistance = orbit.distanceBohr * Math.cos(orbit.elevationRadians);
    camera.position.set(
      equatorialDistance * Math.cos(orbit.azimuthRadians),
      equatorialDistance * Math.sin(orbit.azimuthRadians),
      orbit.distanceBohr * Math.sin(orbit.elevationRadians),
    );
    camera.near = Math.max(0.001, currentOrbitalRadiusBohr / 1000);
    camera.far = Math.max(100, orbit.distanceBohr + currentOrbitalRadiusBohr * 6);
    camera.lookAt(0, 0, 0);
    camera.updateProjectionMatrix();
  }

  function applyTheme(): void {
    const colors = themeColors(appearance.theme);
    renderer.setClearColor(colors.background, 1);
    nucleusMaterial.color.set(colors.nucleus);
    ambientLight.intensity = appearance.theme === 'dark' ? 1.35 : 1.65;
    directionalLight.intensity = appearance.theme === 'dark' ? 1.7 : 1.35;
  }

  function removeCloud(): void {
    if (!cloud) return;
    scene.remove(cloud);
    cloud.geometry.dispose();
    cloud.material.dispose();
    cloud = null;
  }

  function removeDensitySurface(): void {
    if (!densitySurface) return;
    scene.remove(densitySurface);
    densitySurface.geometry.dispose();
    disposeMaterials(densitySurface.material);
    densitySurface = null;
  }

  function removeNodeSurface(): void {
    if (!nodeSurface) return;
    scene.remove(nodeSurface);
    nodeSurface.geometry.dispose();
    disposeMaterials(nodeSurface.material);
    nodeSurface = null;
  }

  function removeGuides(): void {
    if (!guides) return;
    scene.remove(guides);
    for (const geometry of guideGeometries) geometry.dispose();
    for (const material of guideMaterials) material.dispose();
    guideGeometries = [];
    guideMaterials = [];
    guides = null;
  }

  function rebuildGuides(): void {
    removeGuides();
    if (!dataset) return;
    const colors = themeColors(appearance.theme);
    const size = dataset.field.extentBohr * 2;
    const nextGuides = new THREE.Group();
    const grid = new THREE.GridHelper(
      size,
      GRID_DIVISION_COUNT,
      colors.gridCenter,
      colors.gridLines,
    );
    // GridHelper est construit dans xz ; le plan scientifique équatorial est xy.
    grid.rotation.x = Math.PI / 2;
    grid.material.transparent = true;
    grid.material.opacity = appearance.theme === 'dark' ? 0.34 : 0.48;
    const axes = new THREE.AxesHelper(dataset.field.extentBohr * 1.08);
    guideGeometries = [grid.geometry, axes.geometry];
    guideMaterials = [...materialList(grid.material), ...materialList(axes.material)];
    nextGuides.add(grid, axes);
    nextGuides.visible = appearance.showAxes;
    guides = nextGuides;
    scene.add(nextGuides);
  }

  function rebuildCloud(): void {
    removeCloud();
    if (!dataset || appearance.displayMode === 'isosurface') return;
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(dataset.samples.positionsBohr, 3));
    geometry.setAttribute(
      'color',
      new THREE.BufferAttribute(
        colorBufferForPhases(dataset.samples.phaseRadians, appearance.observable),
        3,
      ),
    );
    geometry.computeBoundingSphere();
    cloud = new THREE.Points(geometry, createPointsMaterial(appearance, dotTexture));
    cloud.frustumCulled = true;
    scene.add(cloud);
  }

  function createDensitySurface(): MarchingCubes {
    if (!dataset) throw new Error('Aucun champ orbital à trianguler.');
    const { field } = dataset;
    const material = new THREE.MeshStandardMaterial({
      depthWrite: false,
      metalness: 0,
      opacity: 0.48,
      roughness: 0.58,
      side: THREE.DoubleSide,
      transparent: true,
      vertexColors: true,
    });
    const surface = new MarchingCubes(
      field.resolution,
      material,
      false,
      true,
      maximumMarchingCubeTriangleCount(field.resolution),
    );
    surface.field.set(field.densityNormalized);
    surface.palette.set(colorBufferForPhases(field.phaseRadians, appearance.observable));
    surface.isolation = appearance.isoDensityFraction;
    applyPhysicalGridTransform(surface, field.extentBohr);
    surface.renderOrder = 1;
    surface.update();
    return surface;
  }

  function createNodeSurface(): MarchingCubes {
    if (!dataset) throw new Error('Aucun champ orbital à trianguler.');
    const { field } = dataset;
    const material = new THREE.MeshBasicMaterial({
      color: themeColors(appearance.theme).node,
      depthWrite: false,
      // Les surfaces nodales sont un guide discret, jamais une seconde opaque.
      opacity: appearance.theme === 'dark' ? 0.11 : 0.18,
      side: THREE.DoubleSide,
      transparent: true,
    });
    const surface = new MarchingCubes(
      field.resolution,
      material,
      false,
      false,
      maximumMarchingCubeTriangleCount(field.resolution),
    );
    surface.field.set(field.signedAmplitudeNormalized);
    surface.isolation = 0;
    applyPhysicalGridTransform(surface, field.extentBohr);
    surface.renderOrder = 2;
    surface.update();
    return surface;
  }

  function rebuildSurfaces(): void {
    removeDensitySurface();
    removeNodeSurface();
    if (!dataset) return;

    if (appearance.displayMode !== 'cloud') {
      densitySurface = createDensitySurface();
      scene.add(densitySurface);
    }
    if (appearance.showNodes && dataset.field.nodesAvailable) {
      nodeSurface = createNodeSurface();
      scene.add(nodeSurface);
    }
  }

  function fitCameraToOrbital(): void {
    if (!dataset) return;
    currentOrbitalRadiusBohr = Math.max(Number.EPSILON, orbitalRadiusBohr(dataset));
    const verticalFieldOfViewRadians = THREE.MathUtils.degToRad(camera.fov);
    orbit.distanceBohr =
      (currentOrbitalRadiusBohr / Math.sin(verticalFieldOfViewRadians / 2)) * CAMERA_FIT_MARGIN;
    updateCamera();
  }

  applyTheme();
  updateCamera();

  return {
    dispose(): void {
      removeCloud();
      removeDensitySurface();
      removeNodeSurface();
      removeGuides();
      dotTexture.dispose();
      nucleusGeometry.dispose();
      nucleusMaterial.dispose();
      renderer.dispose();
    },
    fitCameraToOrbital,
    getCameraDistance(): number {
      return orbit.distanceBohr;
    },
    getDiagnostics(): SceneDiagnostics {
      return {
        geometries: renderer.info.memory.geometries,
        textures: renderer.info.memory.textures,
        triangles: renderer.info.render.triangles,
      };
    },
    hasOrbital(): boolean {
      return dataset !== null;
    },
    renderFrame(): void {
      renderer.render(scene, camera);
    },
    resize(viewport): void {
      const width = Math.max(MINIMUM_RENDER_DIMENSION_PIXELS, Math.floor(viewport.clientWidth));
      const height = Math.max(MINIMUM_RENDER_DIMENSION_PIXELS, Math.floor(viewport.clientHeight));
      camera.aspect = width / height;
      camera.updateProjectionMatrix();
      renderer.setSize(width, height, false);
    },
    rotateCamera(azimuthDelta, elevationDelta): void {
      if (!Number.isFinite(azimuthDelta) || !Number.isFinite(elevationDelta)) return;
      orbit.azimuthRadians += azimuthDelta;
      orbit.elevationRadians -= elevationDelta;
      updateCamera();
    },
    rotateCameraAutomatically(deltaRadians): void {
      if (!Number.isFinite(deltaRadians)) return;
      orbit.azimuthRadians += deltaRadians;
      updateCamera();
    },
    setAppearance(nextAppearance): void {
      validateAppearance(nextAppearance);
      const previousAppearance = appearance;
      appearance = { ...nextAppearance };

      const themeChanged = previousAppearance.theme !== appearance.theme;
      const cloudGeometryChanged =
        previousAppearance.observable !== appearance.observable ||
        previousAppearance.displayMode !== appearance.displayMode;
      const cloudMaterialChanged =
        previousAppearance.pointOpacity !== appearance.pointOpacity ||
        previousAppearance.pointSizePixels !== appearance.pointSizePixels;
      const surfacesChanged =
        previousAppearance.displayMode !== appearance.displayMode ||
        previousAppearance.isoDensityFraction !== appearance.isoDensityFraction ||
        previousAppearance.observable !== appearance.observable ||
        previousAppearance.showNodes !== appearance.showNodes ||
        themeChanged;

      if (themeChanged) {
        applyTheme();
        rebuildGuides();
      }
      if (cloudGeometryChanged) {
        rebuildCloud();
      } else if (cloudMaterialChanged && cloud) {
        cloud.material.size = appearance.pointSizePixels;
        cloud.material.opacity = appearance.pointOpacity;
        cloud.material.transparent = appearance.pointOpacity < 1;
        cloud.material.needsUpdate = true;
      }
      if (surfacesChanged) rebuildSurfaces();
      if (previousAppearance.showAxes !== appearance.showAxes && guides) {
        guides.visible = appearance.showAxes;
      }
    },
    setOrbital(nextDataset): void {
      validateDataset(nextDataset);
      removeCloud();
      removeDensitySurface();
      removeNodeSurface();
      dataset = nextDataset;
      currentOrbitalRadiusBohr = orbitalRadiusBohr(nextDataset);

      // Taille volontairement schématique et explicitement non physique.
      schematicNucleus.scale.setScalar(Math.max(0.08, nextDataset.field.extentBohr * 0.006));
      rebuildGuides();
      rebuildCloud();
      rebuildSurfaces();
      fitCameraToOrbital();
    },
    zoomCamera(distanceDelta): void {
      if (!Number.isFinite(distanceDelta)) return;
      const minimumDistance = Math.max(0.1, currentOrbitalRadiusBohr * 0.2);
      const maximumDistance = Math.max(100, currentOrbitalRadiusBohr * 40);
      orbit.distanceBohr = THREE.MathUtils.clamp(
        orbit.distanceBohr + distanceDelta,
        minimumDistance,
        maximumDistance,
      );
      updateCamera();
    },
  };
}

export type { SceneRenderer } from './renderingContracts';
