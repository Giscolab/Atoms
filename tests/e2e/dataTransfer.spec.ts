import { readFile } from 'node:fs/promises';

import { expect, test, type Page } from '@playwright/test';

const GENERATION_TIMEOUT = 45_000;

async function waitForGeneration(page: Page): Promise<void> {
  await expect(page.locator('#generationStatus')).toHaveAttribute('data-visible', 'false', {
    timeout: GENERATION_TIMEOUT,
  });
  await expect(page.locator('#engineStatus')).toContainText('prêt');
}

async function exportSnapshot(page: Page): Promise<{ filename: string; text: string }> {
  const downloadPromise = page.waitForEvent('download');
  await page.locator('#exportSnapshotButton').click();
  const download = await downloadPromise;
  const path = await download.path();
  if (!path) throw new Error('Le snapshot exporté doit être disponible localement pour le test.');
  return {
    filename: download.suggestedFilename(),
    text: await readFile(path, 'utf8'),
  };
}

test('exporte et réimporte un snapshot scientifique reproductible', async ({ page }) => {
  await page.goto('/');
  await waitForGeneration(page);

  const exported = await exportSnapshot(page);
  expect(exported.filename).toMatch(/^atoms-.*-seed-1096044365\.atoms\.json$/);
  const snapshot = JSON.parse(exported.text) as Record<string, unknown>;
  expect(snapshot).toMatchObject({
    format: 'atoms-scientific-snapshot',
    schemaVersion: 1,
  });
  expect(exported.text).not.toContain('timestamp');

  await page.locator('label[for="basisComplex"]').click();
  await page.locator('#quantumN').selectOption('2');
  await page.locator('#quantumL').selectOption('1');
  await page.locator('#quantumM').selectOption('-1');
  await waitForGeneration(page);
  await expect(page.locator('#iOrb')).toHaveText('2,1,-1');

  await page.locator('#snapshotImportInput').setInputFiles({
    buffer: Buffer.from(exported.text, 'utf8'),
    mimeType: 'application/json',
    name: 'restauration.atoms.json',
  });
  await waitForGeneration(page);

  await expect(page.locator('#iOrb')).toHaveText('3d_xy');
  await expect(page.locator('#iSeed')).toHaveText('ATOM');
  await expect(page.locator('#transferStatus')).toContainText('Snapshot importé');
  await expect(page.locator('#transferStatus')).toHaveAttribute('data-state', 'success');
});

test('refuse un snapshot de version incompatible sans altérer l’état courant', async ({ page }) => {
  await page.goto('/');
  await waitForGeneration(page);
  const exported = await exportSnapshot(page);
  const snapshot = JSON.parse(exported.text) as { schemaVersion: number };
  snapshot.schemaVersion = 99;

  await page.locator('#snapshotImportInput').setInputFiles({
    buffer: Buffer.from(`${JSON.stringify(snapshot)}\n`, 'utf8'),
    mimeType: 'application/json',
    name: 'incompatible.atoms.json',
  });

  await expect(page.locator('#transferStatus')).toContainText('Import refusé');
  await expect(page.locator('#transferStatus')).toHaveAttribute('data-state', 'error');
  await expect(page.locator('#iOrb')).toHaveText('3d_xy');
  await expect(page.locator('#iSeed')).toHaveText('ATOM');
});

test('exporte une vraie capture PNG de la vue 3D', async ({ page }) => {
  await page.goto('/');
  await waitForGeneration(page);

  const downloadPromise = page.waitForEvent('download');
  await page.locator('#capturePngButton').click();
  const download = await downloadPromise;
  expect(download.suggestedFilename()).toMatch(/^atoms-.*\.png$/);
  const path = await download.path();
  if (!path) throw new Error('La capture PNG doit être disponible localement pour le test.');
  const png = await readFile(path);

  expect([...png.subarray(0, 8)]).toEqual([137, 80, 78, 71, 13, 10, 26, 10]);
  expect(png.byteLength).toBeGreaterThan(1000);
  await expect(page.locator('#transferStatus')).toContainText('Capture PNG exportée');
});
