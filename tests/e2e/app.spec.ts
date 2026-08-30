import { expect, test, type Page } from '@playwright/test';

const GENERATION_TIMEOUT = 45_000;

function collectRuntimeErrors(page: Page): string[] {
  const runtimeErrors: string[] = [];
  page.on('pageerror', (error) => runtimeErrors.push(`Erreur de page : ${error.message}`));
  page.on('console', (message) => {
    if (message.type() === 'error') runtimeErrors.push(`Erreur console : ${message.text()}`);
  });
  return runtimeErrors;
}

async function waitForGeneration(page: Page): Promise<void> {
  await expect(page.locator('#generationStatus')).toHaveAttribute('data-visible', 'false', {
    timeout: GENERATION_TIMEOUT,
  });
  await expect(page.locator('#engineStatus')).toContainText('prêt');
}

function digitsOnly(value: string | null): string {
  return value?.replace(/\D/g, '') ?? '';
}

test('initialise le rendu WebGL2, le champ et les analyses sans erreur', async ({ page }) => {
  const runtimeErrors = collectRuntimeErrors(page);
  const response = await page.goto('/');
  expect(response?.ok()).toBe(true);
  await waitForGeneration(page);

  const canvas = page.locator('#atomSimCanvas');
  await expect(canvas).toBeVisible();
  const canvasState = await canvas.evaluate((element) => {
    const canvasElement = element as HTMLCanvasElement;
    return {
      clientHeight: canvasElement.clientHeight,
      clientWidth: canvasElement.clientWidth,
      hasWebGl2: canvasElement.getContext('webgl2') !== null,
      height: canvasElement.height,
      width: canvasElement.width,
    };
  });
  expect(canvasState.hasWebGl2).toBe(true);
  expect(canvasState.width).toBeGreaterThan(0);
  expect(canvasState.height).toBeGreaterThan(0);
  expect(canvasState.clientWidth).toBeGreaterThan(0);
  expect(canvasState.clientHeight).toBeGreaterThan(0);

  await expect(page.locator('#quantumNValue')).toHaveText('3');
  await expect(page.locator('#quantumLValue')).toHaveText('2');
  await expect(page.locator('#quantumMValue')).toHaveText('±2');
  await expect(page.locator('#iOrb')).toHaveText('3d_xy');
  await expect
    .poll(async () => digitsOnly(await page.locator('#iPts').textContent()))
    .toBe('15000');
  await expect(page.locator('#radialChart path.radial-series')).toHaveAttribute('d', /M/);
  await expect(page.locator('#angularChart path.angular-series')).toHaveAttribute('d', /M/);
  await expect(page.locator('#basisReal')).toBeChecked();
  await expect(page.locator('#themeDark')).toHaveAttribute('aria-pressed', 'true');
  await expect(page.locator('#iFps')).toHaveText(/^\d+$/, { timeout: GENERATION_TIMEOUT });
  expect(Number(await page.locator('#iFps').textContent())).toBeGreaterThan(0);
  expect(runtimeErrors, runtimeErrors.join('\n')).toEqual([]);
});

test('sépare base complexe/réelle, seed, thème, nœuds et outil 2D', async ({ page }) => {
  const runtimeErrors = collectRuntimeErrors(page);
  await page.goto('/');
  await waitForGeneration(page);

  await page.locator('#themeLight').click();
  await expect(page.locator('html')).toHaveAttribute('data-theme', 'light');
  await page.reload();
  await waitForGeneration(page);
  await expect(page.locator('html')).toHaveAttribute('data-theme', 'light');
  await page.locator('#themeDark').click();
  await expect(page.locator('html')).toHaveAttribute('data-theme', 'dark');

  await page.locator('label[for="basisComplex"]').click();
  await expect(page.locator('#complexMField')).toBeVisible();
  await expect(page.locator('#realOrbitalField')).toBeHidden();
  await page.locator('#quantumN').selectOption('2');
  await page.locator('#quantumL').selectOption('1');
  await page.locator('#quantumM').selectOption('-1');
  await waitForGeneration(page);
  await expect(page.locator('#iOrb')).toHaveText('2,1,-1');
  await expect(page.locator('#nodesToggle')).toBeDisabled();

  await page.locator('label[for="basisReal"]').click();
  await page.locator('#realOrbital').selectOption('p_z');
  await page.locator('#quantumN').selectOption('2');
  await waitForGeneration(page);
  await expect(page.locator('#realOrbitalField')).toBeVisible();
  await expect(page.locator('#complexMField')).toBeHidden();
  await expect(page.locator('#nodesToggle')).toBeEnabled();

  await page.locator('#sampleCount').fill('2000');
  await waitForGeneration(page);
  await expect.poll(async () => digitsOnly(await page.locator('#iPts').textContent())).toBe('2000');
  await page.locator('#seedInput').fill('123');
  await page.locator('#seedInput').press('Tab');
  await waitForGeneration(page);
  await expect(page.locator('#iSeed')).toHaveText('123');

  await page.locator('#btn2d').click();
  await expect(page.locator('#panel2d')).toBeVisible();
  await expect(page.locator('#btn2d')).toHaveAttribute('aria-expanded', 'true');
  await page.locator('#close2d').click();
  await expect(page.locator('#panel2d')).toBeHidden();
  await page.locator('#viewport').click({ position: { x: 12, y: 84 } });
  await page.keyboard.press('q');
  await expect(page.locator('#panel2d')).toBeVisible();
  await page.keyboard.press('q');
  await expect(page.locator('#panel2d')).toBeHidden();

  await page.locator('#atomSimCanvas').focus();
  await page.keyboard.press('ArrowLeft');
  await page.keyboard.press('0');
  await page.locator('#resetCamera').click();
  expect(runtimeErrors, runtimeErrors.join('\n')).toEqual([]);
});

test('reste utilisable sur mobile sans débordement horizontal', async ({ page }) => {
  await page.setViewportSize({ width: 390, height: 844 });
  await page.goto('/');
  await waitForGeneration(page);
  const layout = await page.evaluate(() => ({
    innerWidth: window.innerWidth,
    scrollWidth: document.documentElement.scrollWidth,
    canvasWidth: document.querySelector<HTMLCanvasElement>('#atomSimCanvas')?.clientWidth ?? 0,
    canvasHeight: document.querySelector<HTMLCanvasElement>('#atomSimCanvas')?.clientHeight ?? 0,
  }));
  expect(layout.scrollWidth).toBeLessThanOrEqual(layout.innerWidth + 1);
  expect(layout.canvasWidth).toBeGreaterThan(0);
  expect(layout.canvasHeight).toBeGreaterThan(0);
  await expect(page.locator('#generateButton')).toBeVisible();
});
