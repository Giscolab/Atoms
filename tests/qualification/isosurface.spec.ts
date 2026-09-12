import { expect, test, type Page } from '@playwright/test';
import { cpus, platform, release } from 'node:os';
import { writeFile } from 'node:fs/promises';
import type { SceneDiagnostics } from '../../src/rendering/renderingContracts';

async function settle(page: Page): Promise<void> {
  await expect(page.locator('#engineStatus')).toHaveAttribute('data-state', 'ready');
  await expect(page.locator('#generationStatus')).toHaveAttribute('data-visible', 'false');
  await page.evaluate(
    () =>
      new Promise<void>((resolve) => {
        requestAnimationFrame(() => {
          requestAnimationFrame(() => {
            resolve();
          });
        });
      }),
  );
}

function summarize(values: number[]): { median: number; p95: number; mean: number } {
  const sorted = [...values].sort((a, b) => a - b);
  if (!sorted.length) throw new Error('No timing samples');
  return {
    median: sorted[Math.floor(sorted.length / 2)] ?? 0,
    p95: sorted[Math.ceil(sorted.length * 0.95) - 1] ?? 0,
    mean: values.reduce((sum, value) => sum + value, 0) / values.length,
  };
}

test('qualifie les trois modes à 15 000 points et capture les surfaces', async ({
  page,
  browser,
}, testInfo) => {
  test.setTimeout(180_000);
  const errors: string[] = [];
  const readbackWarnings: string[] = [];
  page.on('pageerror', (error) => errors.push(error.message));
  page.on('console', (message) => {
    // Canvas screenshots intentionally read back pixels; keep the driver's performance notice.
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
  await settle(page);
  await expect(page.locator('#iOrb')).toHaveText('3d_xy');
  await expect(page.locator('#sampleCount')).toHaveValue('15000');
  await expect(page.locator('#seedInput')).toHaveValue('1096044365');
  await expect(page.locator('#isoThreshold')).toHaveValue('0.2');
  await expect(page.locator('#motionToggle')).not.toBeChecked();
  await expect(page.locator('#observablePhase')).toBeChecked();
  await expect(page.locator('#nodesToggle')).toBeChecked();
  await page.locator('#resetCamera').click();
  const environment = await page.evaluate(() => {
    const canvas = document.querySelector<HTMLCanvasElement>('#atomSimCanvas');
    const gl = canvas?.getContext('webgl2');
    if (!gl) throw new Error('WebGL2 required');
    const debug = gl.getExtension('WEBGL_debug_renderer_info');
    return {
      userAgent: navigator.userAgent,
      gpu: String(gl.getParameter(debug ? debug.UNMASKED_RENDERER_WEBGL : gl.RENDERER)),
      devicePixelRatio,
      viewport: { width: innerWidth, height: innerHeight },
    };
  });
  const modes = [];
  for (const mode of ['cloud', 'hybrid', 'isosurface'] as const) {
    await page.locator('#displayMode').selectOption(mode);
    await settle(page);
    const intervals = await page.evaluate(
      () =>
        new Promise<number[]>((resolve) => {
          const intervals: number[] = [];
          let previous: number | null = null;
          let warmup = 12;
          const frame = (now: number): void => {
            if (warmup > 0) warmup--;
            else if (previous !== null) intervals.push(now - previous);
            previous = now;
            if (intervals.length === 120) resolve(intervals);
            else requestAnimationFrame(frame);
          };
          requestAnimationFrame(frame);
        }),
    );
    const resources = await page.evaluate(() => {
      const host = window as unknown as { __atomsDiagnostics?: () => SceneDiagnostics };
      if (!host.__atomsDiagnostics) throw new Error('Diagnostics required');
      return host.__atomsDiagnostics();
    });
    modes.push({ mode, resources, frameIntervalsMs: summarize(intervals), intervals });
    if (mode !== 'cloud') {
      const path = testInfo.outputPath(`dark-${mode}.png`);
      await page.locator('#atomSimCanvas').screenshot({ path });
      await testInfo.attach(`dark-${mode}.png`, {
        path,
        contentType: 'image/png',
      });
    }
  }
  const thresholdDispatchMs = await page.evaluate(() => {
    const input = document.querySelector<HTMLInputElement>('#isoThreshold');
    if (!input) throw new Error('Threshold control required');
    const timings: number[] = [];
    for (let repeat = 0; repeat < 12; repeat++) {
      input.value = repeat % 2 === 0 ? '0.4' : '0.2';
      const start = performance.now();
      input.dispatchEvent(new Event('input', { bubbles: true }));
      timings.push(performance.now() - start);
    }
    return timings.slice(2);
  });
  await page.locator('#themeLight').click();
  await settle(page);
  const lightCapture = testInfo.outputPath('light-isosurface.png');
  await page.locator('#atomSimCanvas').screenshot({ path: lightCapture });
  await testInfo.attach('light-isosurface.png', {
    path: lightCapture,
    contentType: 'image/png',
  });
  const report = {
    recordedAt: new Date().toISOString(),
    environment: {
      ...environment,
      browser: browser.version(),
      node: process.version,
      platform: platform(),
      osRelease: release(),
      cpu: cpus()[0]?.model,
    },
    orbital: '3d_xy',
    seed: 1096044365,
    sampleCount: 15000,
    threshold: 0.2,
    observable: 'phase',
    showNodes: true,
    modes,
    thresholdDispatchMs: { ...summarize(thresholdDispatchMs), samples: thresholdDispatchMs },
    readbackWarnings,
    notes: [
      '120 requestAnimationFrame intervals per mode after warm-up; browser/compositor cadence, not isolated GPU time.',
      'Threshold input dispatch measures synchronous main-thread UI and extraction work, excluding subsequent GPU rendering.',
      'Fixed development runtime and local backend; before/after comparisons require matching environments.',
    ],
  };
  const reportPath = testInfo.outputPath('isosurface-performance.json');
  await writeFile(reportPath, JSON.stringify(report, null, 2));
  await testInfo.attach('isosurface-performance.json', {
    path: reportPath,
    contentType: 'application/json',
  });
  console.log(
    `ISOSURFACE_QUALIFICATION=${JSON.stringify({ environment: report.environment, modes: modes.map(({ mode, resources, frameIntervalsMs }) => ({ mode, resources, frameIntervalsMs })), thresholdDispatchMs: report.thresholdDispatchMs })}`,
  );
  expect(errors).toEqual([]);
});
