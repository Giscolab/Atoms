import { expect, test, type Page } from '@playwright/test';

const GENERATION_TIMEOUT = 45_000;
const SEED = '20260912';
const SAMPLE_COUNT = '2000';

interface VisualCase {
  name: string;
  n: number;
  l: number;
  m: number;
  orbital?: 'p_z' | 'd_z2' | 'd_xy';
  plane: 'xy' | 'xz';
  radialNodes: number;
}

// Node counts and orbital distinctions follow docs/SCIENCE.md and DLMF 14.30.
// The SVG baselines qualify the actual Worker -> UI charts, not the WebGL raster.
const CASES: readonly VisualCase[] = [
  { name: '1s-complex-m0', n: 1, l: 0, m: 0, plane: 'xy', radialNodes: 0 },
  { name: '2s-complex-m0', n: 2, l: 0, m: 0, plane: 'xy', radialNodes: 1 },
  { name: '2p-z-real', n: 2, l: 1, m: 0, orbital: 'p_z', plane: 'xz', radialNodes: 0 },
  { name: '2p-complex-m-plus1', n: 2, l: 1, m: 1, plane: 'xy', radialNodes: 0 },
  { name: '3d-z2-real', n: 3, l: 2, m: 0, orbital: 'd_z2', plane: 'xz', radialNodes: 0 },
  { name: '3d-xy-real', n: 3, l: 2, m: 0, orbital: 'd_xy', plane: 'xy', radialNodes: 0 },
  { name: '4f-complex-m0', n: 4, l: 3, m: 0, plane: 'xz', radialNodes: 0 },
];

async function waitForGeneration(page: Page): Promise<void> {
  await expect(page.locator('#generationStatus')).toHaveAttribute('data-visible', 'false', {
    timeout: GENERATION_TIMEOUT,
  });
  await expect(page.locator('#engineStatus')).toHaveAttribute('data-state', 'ready');
}

async function selectState(page: Page, state: VisualCase): Promise<void> {
  if (state.orbital !== undefined) {
    await page.locator('label[for="basisReal"]').click();
    await page.locator('#realOrbital').selectOption(state.orbital);
    await page.locator('#quantumN').selectOption(String(state.n));
  } else {
    await page.locator('label[for="basisComplex"]').click();
    // First admit every l used by the cases; no reliance on the previous state.
    await page.locator('#quantumN').selectOption('4');
    await page.locator('#quantumL').selectOption(String(state.l));
    await page.locator('#quantumM').selectOption(String(state.m));
    await page.locator('#quantumN').selectOption(String(state.n));
  }
  await waitForGeneration(page);
}

/** Capture the live SVG geometry and its resolved paint, without font rasterization. */
async function captureScientificSvg(page: Page, title: string): Promise<string> {
  return page.evaluate((captureTitle) => {
    const namespace = 'http://www.w3.org/2000/svg';
    const documentSvg = document.createElementNS(namespace, 'svg');
    documentSvg.setAttribute('xmlns', namespace);
    documentSvg.setAttribute('viewBox', '0 0 640 280');
    documentSvg.setAttribute('width', '640');
    documentSvg.setAttribute('height', '280');
    documentSvg.setAttribute('role', 'img');
    const heading = document.createElementNS(namespace, 'title');
    heading.textContent = captureTitle;
    documentSvg.appendChild(heading);
    const backdrop = document.createElementNS(namespace, 'rect');
    backdrop.setAttribute('width', '640');
    backdrop.setAttribute('height', '280');
    backdrop.setAttribute('fill', getComputedStyle(document.body).backgroundColor);
    documentSvg.appendChild(backdrop);

    const properties = [
      'fill',
      'stroke',
      'stroke-width',
      'stroke-linecap',
      'stroke-linejoin',
      'stroke-dasharray',
      'opacity',
      'font-family',
      'font-size',
      'font-weight',
    ];
    for (const [index, id] of ['radialChart', 'angularChart'].entries()) {
      const source = document.getElementById(id);
      if (!(source instanceof SVGSVGElement)) throw new Error(`Missing live chart ${id}`);
      const copy = source.cloneNode(true) as SVGSVGElement;
      copy.setAttribute('x', String(index * 320));
      copy.setAttribute('y', '30');
      copy.setAttribute('width', '320');
      copy.setAttribute('height', id === 'radialChart' ? '190' : '220');
      const originals = source.querySelectorAll('*');
      const copies = copy.querySelectorAll('*');
      originals.forEach((element, elementIndex) => {
        const target = copies[elementIndex];
        if (target === undefined) throw new Error('SVG clone lost an element');
        const computed = getComputedStyle(element);
        target.setAttribute(
          'style',
          properties
            .map((property) => `${property}:${computed.getPropertyValue(property)}`)
            .join(';'),
        );
      });
      documentSvg.appendChild(copy);
    }
    const annotation = document.createElementNS(namespace, 'text');
    annotation.setAttribute('x', '16');
    annotation.setAttribute('y', '18');
    annotation.setAttribute('fill', getComputedStyle(document.body).color);
    annotation.setAttribute('font-family', 'monospace');
    annotation.setAttribute('font-size', '11');
    annotation.textContent = captureTitle.split(' | radial')[0] ?? captureTitle;
    documentSvg.appendChild(annotation);
    const footer = annotation.cloneNode(false) as SVGTextElement;
    footer.setAttribute('y', '272');
    footer.setAttribute('font-size', '10');
    footer.textContent = 'Radial: P(r)/max, r in a0. Angular: |Y|2/max, geometric cut.';
    documentSvg.appendChild(footer);
    // XML is directly viewable and diffable, with no OS-specific image baseline.
    return `${new XMLSerializer().serializeToString(documentSvg)}\n`;
  }, title);
}

test.use({
  viewport: { width: 1440, height: 1000 },
  deviceScaleFactor: 1,
  colorScheme: 'dark',
  contextOptions: { reducedMotion: 'reduce' },
  locale: 'fr-FR',
  timezoneId: 'Europe/Paris',
});

for (const state of CASES) {
  test(`capture scientifique déterministe : ${state.name}`, async ({ page }, testInfo) => {
    test.setTimeout(90_000);
    const runtimeErrors: string[] = [];
    page.on('pageerror', (error) => {
      runtimeErrors.push(error.message);
    });
    await page.addInitScript(() => localStorage.setItem('atoms-theme', 'dark'));
    expect((await page.goto('/'))?.ok()).toBe(true);
    await waitForGeneration(page);
    await expect(page.locator('#motionToggle')).not.toBeChecked();
    await page.locator('#sampleCount').fill(SAMPLE_COUNT);
    await waitForGeneration(page);
    await page.locator('#seedInput').fill(SEED);
    await page.locator('#seedInput').press('Tab');
    await waitForGeneration(page);
    await selectState(page, state);
    await page.locator('#displayMode').selectOption('cloud');
    await page.locator('label[for="observableDensity"]').click();
    await page.locator('#resetCamera').click();

    await expect(page.locator('#themeDark')).toHaveAttribute('aria-pressed', 'true');
    await expect(page.locator('#sampleCount')).toHaveValue(SAMPLE_COUNT);
    await expect(page.locator('#iSeed')).toHaveText(SEED);
    await expect(page.locator('#basisLabel')).toHaveText(
      state.orbital === undefined ? 'complexe' : 'réelle',
    );
    await expect(page.locator('#iOrb')).toHaveText(
      state.orbital === undefined
        ? `${state.n},${state.l},${state.m}`
        : `${state.n}${state.orbital}`,
    );
    await expect(page.locator('#radialNodesValue')).toHaveText(String(state.radialNodes));
    await expect(page.locator('#angularNodesValue')).toHaveText(String(state.l));
    await expect(page.locator('#totalNodesValue')).toHaveText(String(state.n - 1));
    await expect(page.locator('#angularPlaneLabel')).toHaveText(`Plan ${state.plane}`);
    await expect(page.locator('#radialChartLabel')).toContainText('a₀');
    await expect(page.locator('#radialChart .radial-series')).toHaveAttribute('d', /^M .+ L /);
    await expect(page.locator('#angularChart .angular-series')).toHaveAttribute('d', /^M .+ Z$/);
    expect(
      await page
        .locator('#atomSimCanvas')
        .evaluate((element) => (element as HTMLCanvasElement).getContext('webgl2') !== null),
    ).toBe(true);

    const title = `${state.name} | seed ${SEED} | ${SAMPLE_COUNT} points | radial P(r)/max, r in a0 | angular |Y|2/max, ${state.plane}`;
    const first = await captureScientificSvg(page, title);
    const cameraDistance = await page.locator('#iDist').textContent();
    expect(cameraDistance).toMatch(/a₀/u);
    expect(first).toMatchSnapshot(`${state.name}.svg`);

    // A new real Worker generation must leave the scientific capture unchanged.
    await page.locator('#generateButton').click();
    await waitForGeneration(page);
    await page.locator('#resetCamera').click();
    await expect(page.locator('#iDist')).toHaveText(cameraDistance ?? '');
    expect(await captureScientificSvg(page, title)).toBe(first);
    await testInfo.attach(`${state.name}-scientific.svg`, {
      body: first,
      contentType: 'image/svg+xml',
    });
    for (const chart of ['radialChart', 'angularChart']) {
      await testInfo.attach(`${state.name}-${chart}.png`, {
        body: await page.locator(`#${chart}`).screenshot({ animations: 'disabled', caret: 'hide' }),
        contentType: 'image/png',
      });
    }
    expect(runtimeErrors).toEqual([]);
  });
}
