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

describe('frontières du Lot 1.2', () => {
  it('garde la couche scientifique indépendante de Three.js et du navigateur', async () => {
    const scienceDirectory = join(sourceRoot, 'science');
    const scienceFiles = collectTypeScriptFiles(scienceDirectory);
    expect(scienceFiles.length).toBeGreaterThan(0);

    const forbiddenPatterns = [
      /(?:from\s+|import\s*\()\s*['"]three(?:\/[^'"]*)?['"]/u,
      /\b(?:document|window|HTMLElement|HTMLCanvasElement|CanvasRenderingContext2D|WebGLRenderer)\b/u,
      /(?:from\s+|import\s*\()\s*['"]\.\.\/(?:app|data|rendering|state|ui)(?:\/|['"])/u,
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

  it('laisse main comme point de composition et ne crée aucun Worker prématuré', () => {
    const mainSource = readSource(join(sourceRoot, 'main.ts'));
    expect(mainSource).not.toMatch(
      /\bTHREE\b|from\s+['"]three|Math\.random|setTimeout|requestAnimationFrame|-13\.6|52\.9|5\.29/u,
    );
    expect(existsSync(join(sourceRoot, 'workers'))).toBe(false);
  });
});
