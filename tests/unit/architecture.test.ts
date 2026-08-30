import { existsSync, readFileSync, readdirSync } from 'node:fs';
import { join, relative } from 'node:path';
import { fileURLToPath } from 'node:url';

import { describe, expect, it } from 'vitest';

const sourceRoot = fileURLToPath(new URL('../../src/', import.meta.url));

function collectTypeScriptFiles(directory: string): string[] {
  const files: string[] = [];
  for (const entry of readdirSync(directory, { withFileTypes: true })) {
    const path = join(directory, entry.name);
    if (entry.isDirectory()) files.push(...collectTypeScriptFiles(path));
    else if (entry.isFile() && entry.name.endsWith('.ts')) files.push(path);
  }
  return files.sort();
}

function sourceName(path: string): string {
  return relative(sourceRoot, path).replaceAll('\\', '/');
}

function readSource(path: string): string {
  return readFileSync(path, 'utf8');
}

describe('frontières architecturales', () => {
  it('garde la couche scientifique indépendante de Three.js et du navigateur', async () => {
    const scienceDirectory = join(sourceRoot, 'science');
    const scienceFiles = collectTypeScriptFiles(scienceDirectory);
    expect(scienceFiles.length).toBeGreaterThan(0);

    const forbiddenPatterns = [
      /(?:from\s+|import\s*\()\s*['"]three(?:\/[^'"]*)?['"]/u,
      /\b(?:document|window|HTMLElement|HTMLCanvasElement|CanvasRenderingContext2D|OffscreenCanvas|FileReader|WebGLRenderer|Worker|SharedWorker)\b/u,
      /(?:from\s+|import\s*\()\s*['"]node:/u,
      /(?:from\s+|import\s*\()\s*['"](?:\.\.\/)+(?:app|data|rendering|sampling|state|ui|workers)(?:\/|['"])/u,
    ];

    for (const path of scienceFiles) {
      const source = readSource(path);
      for (const pattern of forbiddenPatterns) {
        expect(source, `${sourceName(path)} contient ${String(pattern)}`).not.toMatch(pattern);
      }
    }

    await expect(import('../../src/science/legacyScience')).resolves.toBeDefined();
  });

  it('confine Three.js au rendu et les formules legacy hors du rendu et de l’UI', () => {
    const sourceFiles = collectTypeScriptFiles(sourceRoot);
    const threeImporters = sourceFiles
      .filter((path) =>
        /(?:from\s+|import\s*\()\s*['"]three(?:\/[^'"]*)?['"]/u.test(readSource(path)),
      )
      .map(sourceName);
    expect(threeImporters.length).toBeGreaterThan(0);
    expect(threeImporters.every((path) => path.startsWith('rendering/'))).toBe(true);

    const formulaDefinitionPattern =
      /-13\.6|52\.9|5\.29|NAMED_SAMPLERS|\bfunction\s+(?:gamma|buildRCDF|buildTCDF|legacyOrbitalIntensity|legacyProbabilityFlow|sampleLegacyOrbitalPoint)\b/u;
    for (const directoryName of ['rendering', 'ui']) {
      const files = collectTypeScriptFiles(join(sourceRoot, directoryName));
      for (const path of files) {
        expect(readSource(path), sourceName(path)).not.toMatch(formulaDefinitionPattern);
      }
    }

    for (const path of collectTypeScriptFiles(join(sourceRoot, 'ui'))) {
      expect(readSource(path), sourceName(path)).not.toMatch(
        /legacyOrbitalColor|legacyOrbitalIntensity|sampleLegacyOrbitalPoint/u,
      );
    }
  });

  it('laisse main comme point de composition et isole le sampler scientifique', () => {
    const mainSource = readSource(join(sourceRoot, 'main.ts'));
    expect(mainSource).not.toMatch(
      /\bTHREE\b|from\s+['"]three|Math\.random|setTimeout|requestAnimationFrame|-13\.6|52\.9|5\.29/u,
    );

    const samplingDirectory = join(sourceRoot, 'sampling');
    expect(existsSync(samplingDirectory)).toBe(true);
    const forbiddenSamplingPatterns = [
      /(?:from\s+|import\s*\()\s*['"]three(?:\/[^'"]*)?['"]/u,
      /\b(?:document|window|HTMLElement|HTMLCanvasElement|CanvasRenderingContext2D|OffscreenCanvas|FileReader|WebGLRenderer|Worker|SharedWorker)\b/u,
      /(?:from\s+|import\s*\()\s*['"]node:/u,
      /(?:from\s+|import\s*\()\s*['"](?:\.\.\/)+(?:app|data|rendering|state|ui|workers)(?:\/|['"])/u,
      /\bMath\.random\b/u,
      /\b(?:SAMPLER_MAP|NAMED_SAMPLERS|legacyScience)\b/u,
    ];
    for (const path of collectTypeScriptFiles(samplingDirectory)) {
      const source = readSource(path);
      for (const pattern of forbiddenSamplingPatterns) {
        expect(source, `${sourceName(path)} contient ${String(pattern)}`).not.toMatch(pattern);
      }
    }
  });

  it('branche le nouveau noyau dans les frontières applicatives prévues', () => {
    const runtimeFiles = collectTypeScriptFiles(sourceRoot).filter((path) => {
      const name = sourceName(path);
      return !name.startsWith('science/') && !name.startsWith('sampling/');
    });
    const legacyBoundaryFiles = new Set([
      'main.ts',
      'rendering/photoelectric2d.ts',
      'rendering/legacyPhotoelectricModel.ts',
    ]);
    for (const path of runtimeFiles) {
      const name = sourceName(path);
      const source = readSource(path);
      if (!legacyBoundaryFiles.has(name)) {
        expect(source, `${name} ne doit pas importer la couche legacy`).not.toMatch(
          /(?:from\s+|import\s*\()\s*['"][^'"]*\/science\/legacyScience['"]/u,
        );
      }
    }
    expect(readSource(join(sourceRoot, 'app', 'orbitalPresentation.ts'))).toMatch(
      /\/science\/hydrogen\//u,
    );
    expect(readSource(join(sourceRoot, 'app', 'orbitalWorkerClient.ts'))).toMatch(/\/sampling\//u);
  });

  it('confine le Worker à l’orchestration du sampler et garde Three.js hors de sa couche', () => {
    const workersDirectory = join(sourceRoot, 'workers');
    expect(existsSync(workersDirectory)).toBe(true);
    const forbiddenPatterns = [
      /(?:from\s+|import\s*\()\s*['"]three(?:\/[^'"]*)?['"]/u,
      /(?:from\s+|import\s*\()\s*['"](?:\.\.\/)+(?:data|rendering|state|ui)(?:\/|['"])/u,
      /\b(?:document|window|HTMLElement|HTMLCanvasElement|WebGLRenderer|Math\.random)\b/u,
      /\b(?:SAMPLER_MAP|NAMED_SAMPLERS|legacyScience)\b/u,
    ];
    for (const path of collectTypeScriptFiles(workersDirectory)) {
      const source = readSource(path);
      for (const pattern of forbiddenPatterns) {
        expect(source, `${sourceName(path)} contient ${String(pattern)}`).not.toMatch(pattern);
      }
    }
  });

  it('n’introduit ni sampling ni dynamique dans le noyau scientifique de Phase 3', () => {
    const phaseThreeDirectories = ['hydrogen', 'math'];
    const forbiddenPatterns = [
      /\bMath\.random\b/u,
      /\b(?:CDF|RNG)\b/iu,
      /\b(?:random|sample|sampler|sampling)\b/iu,
      /\b(?:velocity|trajectory|Worker|SharedWorker)\b/iu,
    ];

    for (const directoryName of phaseThreeDirectories) {
      for (const path of collectTypeScriptFiles(join(sourceRoot, 'science', directoryName))) {
        const source = readSource(path);
        for (const pattern of forbiddenPatterns) {
          expect(source, `${sourceName(path)} contient ${String(pattern)}`).not.toMatch(pattern);
        }
      }
    }
  });
});
