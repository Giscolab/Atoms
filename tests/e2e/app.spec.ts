import { expect, test, type Page } from '@playwright/test';

const GENERATION_TIMEOUT = 30_000;

function collectRuntimeErrors(page: Page): string[] {
  const runtimeErrors: string[] = [];

  page.on('pageerror', (error) => {
    runtimeErrors.push(`Erreur de page : ${error.message}`);
  });

  page.on('console', (message) => {
    if (message.type() === 'error') {
      runtimeErrors.push(`Erreur console : ${message.text()}`);
    }
  });

  return runtimeErrors;
}

async function waitForGeneration(page: Page): Promise<void> {
  await expect(page.locator('#loading')).toBeHidden({
    timeout: GENERATION_TIMEOUT,
  });
  await expect(page.locator('#progressBar')).toHaveText('');
}

function digitsOnly(value: string | null): string {
  return value?.replace(/\D/g, '') ?? '';
}

test('initialise le rendu WebGL2 et le HUD sans erreur', async ({ page }) => {
  const runtimeErrors = collectRuntimeErrors(page);

  const response = await page.goto('/');
  expect(response?.ok()).toBe(true);

  const canvas = page.locator('#atomSimCanvas');
  await expect(canvas).toBeVisible();
  await waitForGeneration(page);

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

  await expect(page.locator('#qn-n')).toHaveText('2');
  await expect(page.locator('#qn-l')).toHaveText('1');
  await expect(page.locator('#qn-m')).toHaveText('0');
  await expect(page.locator('#iOrb')).toHaveText('2, 1, 0');
  await expect
    .poll(async () => digitsOnly(await page.locator('#iPts').textContent()))
    .toBe('15000');

  const fps = page.locator('#iFps');
  await expect(fps).toHaveText(/^\d+$/, { timeout: GENERATION_TIMEOUT });
  expect(Number(await fps.textContent())).toBeGreaterThan(0);

  expect(runtimeErrors, runtimeErrors.join('\n')).toEqual([]);
});

test('préserve les contrôles essentiels à 2 000 points', async ({ page }) => {
  const runtimeErrors = collectRuntimeErrors(page);

  await page.goto('/');
  await waitForGeneration(page);

  const atomTab = page.getByRole('tab', { name: 'Définir atome' });
  const presetsTab = page.getByRole('tab', { name: 'Préréglages' });
  const atomPanel = page.locator('#panel-view-atom');
  const presetsPanel = page.locator('#panel-view-presets');

  await atomTab.click();
  await expect(atomTab).toHaveAttribute('aria-selected', 'true');
  await expect(atomPanel).toBeVisible();
  await expect(presetsPanel).toBeHidden();

  await presetsTab.click();
  await expect(presetsTab).toHaveAttribute('aria-selected', 'true');
  await expect(presetsPanel).toBeVisible();
  await expect(atomPanel).toBeHidden();

  const themeToggle = page.locator('#themeToggle');
  const initialThemeIsLight = await themeToggle.isChecked();
  const expectedTheme = initialThemeIsLight ? 'dark' : 'light';

  await themeToggle.setChecked(!initialThemeIsLight);
  await expect(page.locator('html')).toHaveAttribute('data-theme', expectedTheme);

  const particleInput = page.locator('#nParts');
  await page.locator('#orbitalSel').selectOption('2_1_1');
  await waitForGeneration(page);
  await particleInput.evaluate((element) => {
    const input = element as HTMLInputElement;
    input.value = '20000';
    input.dispatchEvent(new Event('input', { bubbles: true }));
  });
  await page.evaluate(async () => {
    await new Promise<void>((resolve) => {
      requestAnimationFrame(() => {
        requestAnimationFrame(() => {
          resolve();
        });
      });
    });
  });
  expect(runtimeErrors, runtimeErrors.join('\n')).toEqual([]);

  await particleInput.fill('2000');
  await waitForGeneration(page);
  await expect.poll(async () => digitsOnly(await page.locator('#iPts').textContent())).toBe('2000');

  await page.locator('#orbitalSel').selectOption('1_0_0');
  await waitForGeneration(page);
  await expect(page.locator('#qn-n')).toHaveText('1');
  await expect(page.locator('#qn-l')).toHaveText('0');
  await expect(page.locator('#qn-m')).toHaveText('0');
  await expect(page.locator('#iOrb')).toHaveText('1, 0, 0');

  const panel2d = page.locator('#panel2d');
  const button2d = page.locator('#btn2d');

  await expect(panel2d).toBeHidden();

  await button2d.click();
  await expect(panel2d).toBeVisible();
  expect(await panel2d.evaluate((element) => element.hasAttribute('hidden'))).toBe(false);
  await expect(button2d).toHaveClass(/active/);

  await page.keyboard.press('q');
  await expect(panel2d).toBeHidden();
  expect(await panel2d.evaluate((element) => element.hasAttribute('hidden'))).toBe(true);
  await expect(button2d).not.toHaveClass(/active/);

  expect(runtimeErrors, runtimeErrors.join('\n')).toEqual([]);
});
