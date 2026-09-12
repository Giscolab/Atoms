import AxeBuilder from '@axe-core/playwright';
import { expect, test, type Page } from '@playwright/test';

const GENERATION_TIMEOUT = 45_000;
const WCAG_TAGS = ['wcag2a', 'wcag2aa', 'wcag21a', 'wcag21aa'];

async function loadReadyApp(page: Page): Promise {
  const response = await page.goto('/', { waitUntil: 'domcontentloaded' });
  expect(response?.ok()).toBe(true);
  await expect(page.locator('#generationStatus')).toHaveAttribute('data-visible', 'false', {
    timeout: GENERATION_TIMEOUT,
  });
  await expect(page.locator('#engineStatus')).toContainText('prêt', {
    timeout: GENERATION_TIMEOUT,
  });
}

async function expectNoAutomatedViolations(page: Page): Promise {
  const scan = await new AxeBuilder({ page }).withTags(WCAG_TAGS).analyze();
  const violations = scan.violations.map(({ id, impact, nodes }) => ({
    id,
    impact,
    targets: nodes.map((node) => node.target),
  }));
  expect(violations).toEqual([]);
}

test.beforeEach(async ({ page }, testInfo) => {
  // Fix 1 & 2: Augmenter le timeout avant d'exécuter loadReadyApp
  test.setTimeout(90_000);

  // Fix 1: Ignorer immédiatement sous Firefox AVANT de charger l'application
  if (
    testInfo.title.includes('sans violation WCAG') &&
    testInfo.project.name !== 'chromium'
  ) {
    test.skip(true, 'Audit axe exécuté une fois sous Chromium.');
    return;
  }

  await loadReadyApp(page);
});

test('thème sombre sans violation WCAG A/AA automatisable', async ({ page }) => {
  await expect(page.locator('#themeDark')).toHaveAttribute('aria-pressed', 'true');
  await expectNoAutomatedViolations(page);
});

test('thème clair sans violation WCAG A/AA automatisable', async ({ page }) => {
  await page.locator('#themeLight').click();
  await expect(page.locator('#themeLight')).toHaveAttribute('aria-pressed', 'true');
  await expectNoAutomatedViolations(page);
});

/** Inspect the rendered focus indicator, including labels for visually hidden inputs.
 * A nonzero CSS outline alone is insufficient: overflow can clip its entire perimeter.
 */
async function expectVisibleKeyboardFocus(page: Page): Promise {
  const focused = page.locator(':focus');
  await expect(focused).toHaveCount(1);

  const readFocus = async () =>
    focused.evaluate((element) => {
      const indicator = element.matches('.toggle-row input')
        ? element.closest('.toggle-row')
        : element.matches('.segmented input')
          ? element.nextElementSibling
          : element;
      if (!(indicator instanceof HTMLElement)) throw new Error('Indicateur de focus absent.');
      const style = getComputedStyle(indicator);
      const rect = indicator.getBoundingClientRect();
      return {
        id:
          element.id ||
          (element.matches('.insight-panel')
            ? 'analyses'
            : element.matches('.panel-scroll')
              ? 'controls'
              : element.tagName),
        keyboard: element.matches(':focus-visible'),
        visible:
          style.visibility === 'visible' &&
          Number(style.opacity) > 0 &&
          rect.width > 0 &&
          rect.height > 0 &&
          rect.right > 0 &&
          rect.bottom > 0 &&
          rect.left < innerWidth &&
          rect.top < innerHeight,
        outline:
          style.outlineStyle !== 'none' &&
          Number.parseFloat(style.outlineWidth) >= 2 &&
          style.outlineColor !== 'transparent' &&
          style.outlineColor !== 'rgba(0, 0, 0, 0)',
      };
    });

  const initial = await readFocus();
  await expect
    .poll(readFocus, {
      timeout: 3000,
      message: `Focus visible pour ${initial.id}`,
    })
    .toMatchObject({
      id: initial.id,
      keyboard: true,
      visible: true,
      outline: true,
    });
  return initial.id;
}

for (const theme of ['dark', 'light'] as const) {
  for (const basis of ['real', 'complex'] as const) {
    test(`parcours Tab et focus visible : \({theme}, base\){basis}`, async ({ page }) => {
      const visited: string[] = [];
      // No locator.focus(), clicks or assigned tabindex: Tab reaches every group.
      // Native radio groups expose one Tab stop; their choices use arrow keys.
      for (let step = 0; step < 80; step += 1) {
        await page.keyboard.press('Tab');
        let id = await expectVisibleKeyboardFocus(page);
        if (id === (theme === 'dark' ? 'themeDark' : 'themeLight')) {
          await page.keyboard.press('Enter');
          await expect(page.locator('html')).toHaveAttribute('data-theme', theme);
          await expectVisibleKeyboardFocus(page);
        }
        if (id === 'basisReal' || id === 'basisComplex') {
          await page.keyboard.press('ArrowLeft');
          await expect(page.locator('#basisComplex')).toBeChecked();
          if (basis === 'real') {
            await page.keyboard.press('ArrowRight');
            await expect(page.locator('#basisReal')).toBeChecked();
          }
          id = basis === 'real' ? 'basisReal' : 'basisComplex';
        }
        if (id === 'observablePhase') {
          await page.keyboard.press('ArrowLeft');
          await expect(page.locator('#observableDensity')).toBeChecked();
          await page.keyboard.press('ArrowRight');
          await expect(page.locator('#observablePhase')).toBeChecked();
        }
        visited.push(id);
        if (id === 'generateButton') {
          await page.keyboard.press('Enter');
          await expect(page.locator('#generationStatus')).toHaveAttribute('data-visible', 'false', {
            timeout: GENERATION_TIMEOUT,
          });
          await expect(page.locator('#engineStatus')).toContainText('prêt');
        }
        if (id === 'resetCamera') await page.keyboard.press('Enter');
        if (id === 'atomSimCanvas') {
          const distance = await page.locator('#iDist').textContent();
          await page.keyboard.press('+');
          await expect(page.locator('#iDist')).not.toHaveText(distance ?? '');
          await page.keyboard.press('0');
          await expect(page.locator('#iDist')).toHaveText(distance ?? '');
        }
        if (id === 'analyses') break;
      }
      const groups = [
        ['themeDark', 'themeLight', 'settingsButton'],
        [
          basis === 'real' ? 'basisReal' : 'basisComplex',
          'quantumN',
          ...(basis === 'real' ? ['realOrbital'] : ['quantumL', 'quantumM']),
        ],
        ['observablePhase', 'displayMode', 'isoThreshold'],
        ['sampleCount', 'seedInput', 'randomizeSeed'],
        ['pointOpacity', 'pointSize', 'axesToggle', 'nodesToggle', 'motionToggle'],
        ['generateButton'],
        ['resetCamera', 'atomSimCanvas'],
        ['analyses'],
      ];
      let previousGroupEnd = -1;
      for (const group of groups) {
        for (const id of group) expect(visited, `${id} accessible par Tab`).toContain(id);
        const positions = group.map((id) => visited.indexOf(id));
        expect(
          Math.min(...positions),
          `Ordre logique du groupe ${group.join(', ')}`,
        ).toBeGreaterThan(previousGroupEnd);
        previousGroupEnd = Math.max(...positions);
      }
      expect(visited).not.toContain(basis === 'real' ? 'quantumM' : 'realOrbital');
      if (basis === 'real') expect(visited).not.toContain('quantumL'); // Derived from real orbital.
      const analyses = page.locator('.insight-panel');
      const initialScroll = await analyses.evaluate((element) => element.scrollTop);
      await page.keyboard.press('ArrowDown');
      await expect
        .poll(() => analyses.evaluate((element) => element.scrollTop))
        .toBeGreaterThan(initialScroll);
      await page.keyboard.press('Shift+Tab');
      await expect(page.locator('#atomSimCanvas')).toBeFocused();
      await expectVisibleKeyboardFocus(page);
    });
  }
}