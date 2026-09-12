import { expect, test, type Page } from '@playwright/test';
import type { SceneDiagnostics } from '../../src/rendering/renderingContracts';

test.use({ contextOptions: { reducedMotion: 'reduce' } });

async function diagnostics(page: Page): Promise<SceneDiagnostics> {
  return page.evaluate(
    () =>
      new Promise<SceneDiagnostics>((resolve, reject) => {
        requestAnimationFrame(() => {
          requestAnimationFrame(() => {
            const host = window as unknown as { __atomsDiagnostics?: () => SceneDiagnostics };
            if (!host.__atomsDiagnostics) reject(new Error('Development diagnostics required'));
            else resolve(host.__atomsDiagnostics());
          });
        });
      }),
  );
}

test('préserve la géométrie entre apparences, change le seuil et stabilise les ressources', async ({
  page,
}, testInfo) => {
  test.setTimeout(120_000);
  const errors: string[] = [];
  const readbackWarnings: string[] = [];
  page.on('pageerror', (error) => errors.push(error.message));
  page.on('console', (message) => {
    if (message.type() === 'warning' && /GPU stall due to ReadPixels/.test(message.text())) {
      readbackWarnings.push(message.text());
      return;
    }
    if (
      message.type() === 'error' ||
      (message.type() === 'warning' && /WebGL|THREE|shader/i.test(message.text()))
    ) {
      errors.push(message.text());
    }
  });
  await page.addInitScript(() => {
    localStorage.setItem('atoms-theme', 'dark');
  });
  expect((await page.goto('/Atoms/?diagnostics=1'))?.ok()).toBe(true);
  await expect(page.locator('#engineStatus')).toHaveAttribute('data-state', 'ready', {
    timeout: 45_000,
  });
  await expect(page.locator('#generationStatus')).toHaveAttribute('data-visible', 'false');
  await expect(page.locator('#motionToggle')).not.toBeChecked();
  await expect(page.locator('#iOrb')).toHaveText('3d_xy');
  const hybrid = await diagnostics(page);
  expect(hybrid.cloudPoints).toBe(15000);
  expect(hybrid.surfaceTriangles).toBeGreaterThan(0);
  expect(hybrid.surfaceVertices).toBeGreaterThan(0);
  expect(hybrid.surfaceFingerprint).not.toBe('');

  await page.locator('#displayMode').selectOption('cloud');
  const cloud = await diagnostics(page);
  expect(cloud.cloudPoints).toBe(15000);
  expect(cloud.surfaceVertices).toBe(0);
  expect(cloud.surfaceTriangles).toBe(0);
  await page.locator('#displayMode').selectOption('isosurface');
  const isolated = await diagnostics(page);
  expect(isolated.cloudPoints).toBe(0);
  expect(isolated.surfaceFingerprint).toBe(hybrid.surfaceFingerprint);
  await page.locator('#isoThreshold').fill('0.4');
  const raised = await diagnostics(page);
  expect(raised.surfaceTriangles).toBeGreaterThan(0);
  expect(raised.surfaceFingerprint).not.toBe(hybrid.surfaceFingerprint);
  await page.locator('#isoThreshold').fill('0.2');
  expect((await diagnostics(page)).surfaceFingerprint).toBe(hybrid.surfaceFingerprint);
  await page.locator('label[for="observablePhase"]').click();
  expect((await diagnostics(page)).surfaceFingerprint).toBe(hybrid.surfaceFingerprint);
  await page.locator('#themeLight').click();
  await expect(page.locator('html')).toHaveAttribute('data-theme', 'light');
  expect((await diagnostics(page)).surfaceFingerprint).toBe(hybrid.surfaceFingerprint);

  // Each comparison returns to exactly the same state; the first cycle warms lazy programs.
  let stable: SceneDiagnostics | undefined;
  for (let cycle = 0; cycle < 5; cycle++) {
    await page.locator('#displayMode').selectOption('cloud');
    await diagnostics(page);
    await page.locator('#displayMode').selectOption('isosurface');
    await page.locator('#isoThreshold').fill('0.4');
    await diagnostics(page);
    await page.locator('#themeLight').click();
    await page.locator('label[for="observablePhase"]').click();
    await diagnostics(page);
    await page.locator('#isoThreshold').fill('0.2');
    await page.locator('#displayMode').selectOption('hybrid');
    await page.locator('#themeDark').click();
    await page.locator('label[for="observableDensity"]').click();
    const current = await diagnostics(page);
    expect(current.surfaceFingerprint).toBe(hybrid.surfaceFingerprint);
    expect(current.cloudPoints).toBe(15000);
    if (stable) expect(current).toEqual(stable);
    if (cycle === 0) stable = current;
  }
  await testInfo.attach('readback-performance-warnings.json', {
    body: JSON.stringify(readbackWarnings),
    contentType: 'application/json',
  });
  expect(errors).toEqual([]);
});
