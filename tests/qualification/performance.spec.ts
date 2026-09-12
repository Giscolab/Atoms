import { expect, test, type Page } from '@playwright/test';
import { cpus, platform, release, totalmem } from 'node:os';
import type { OrbitalWorkerJob, OrbitalWorkerResponse } from '../../src/workers/orbitalSamplingProtocol';
import type { SceneDiagnostics } from '../../src/rendering/renderingContracts';

interface ObservedJob {
  jobId: string;
  sampleCount: number;
  seed: number;
  state: OrbitalWorkerJob['state'];
  startedMs: number;
  stages: Partial<Record<'sampling' | 'field' | 'charts' | 'transfer', number>>;
  finishedMs: number | null;
  status: 'pending' | 'result' | 'cancelled' | 'error';
  transferBytes: number;
}

interface QualificationProbe {
  jobs: ObservedJob[];
  liveWorkers: number;
  peakWorkers: number;
  buffersCreated: number;
  buffersDeleted: number;
  texturesCreated: number;
  texturesDeleted: number;
  readyAtMs: number | null;
}

declare global {
  interface Window {
    __qualificationProbe: QualificationProbe;
  }
}

// This observer retains scalar metadata only, never Worker payloads or WebGL handles.
function installProbe(): void {
  const probe: QualificationProbe = {
    jobs: [], liveWorkers: 0, peakWorkers: 0,
    buffersCreated: 0, buffersDeleted: 0, texturesCreated: 0, texturesDeleted: 0,
    readyAtMs: null,
  };
  window.__qualificationProbe = probe;
  const NativeWorker = window.Worker;
  window.Worker = class extends NativeWorker {
    private terminated = false;
    private observedJob: ObservedJob | undefined;

    constructor(url: string | URL, options?: WorkerOptions) {
      super(url, options);
      probe.liveWorkers++;
      probe.peakWorkers = Math.max(probe.peakWorkers, probe.liveWorkers);
      this.addEventListener('message', (event: MessageEvent<OrbitalWorkerResponse>) => {
        const job = this.observedJob;
        const response = event.data;
        if (!job || response.jobId !== job.jobId) return;
        if (response.kind === 'progress') {
          job.stages[response.stage] = performance.now();
        } else {
          job.finishedMs = performance.now();
          job.status = response.kind;
          if (response.kind === 'result') {
            const seen = new Set<ArrayBufferLike>();
            const countBuffers = (value: unknown): void => {
              if (ArrayBuffer.isView(value)) seen.add(value.buffer);
              else if (value !== null && typeof value === 'object') {
                for (const child of Object.values(value)) countBuffers(child);
              }
            };
            countBuffers(response);
            job.transferBytes = [...seen].reduce((sum, buffer) => sum + buffer.byteLength, 0);
          }
        }
      });
    }

    override postMessage(message: unknown, options?: Transferable[] | StructuredSerializeOptions): void {
      // The app owns this Worker and sends only its explicit generation protocol.
      const job = message as OrbitalWorkerJob;
      this.observedJob = {
        jobId: job.jobId, sampleCount: job.sampleCount, seed: job.seed, state: job.state,
        startedMs: performance.now(), stages: {}, finishedMs: null, status: 'pending', transferBytes: 0,
      };
      probe.jobs.push(this.observedJob);
      if (Array.isArray(options)) super.postMessage(message, options);
      else super.postMessage(message, options);
    }

    override terminate(): void {
      if (!this.terminated) {
        this.terminated = true;
        probe.liveWorkers--;
        if (this.observedJob?.status === 'pending') {
          this.observedJob.status = 'cancelled';
          this.observedJob.finishedMs = performance.now();
        }
      }
      super.terminate();
    }
  };

  const gl = WebGL2RenderingContext.prototype;
  const createBuffer = gl.createBuffer;
  const deleteBuffer = gl.deleteBuffer;
  const createTexture = gl.createTexture;
  const deleteTexture = gl.deleteTexture;
  gl.createBuffer = function (): WebGLBuffer | null {
    const buffer = createBuffer.call(this);
    if (buffer) probe.buffersCreated++;
    return buffer;
  };
  gl.deleteBuffer = function (buffer): void {
    if (buffer) probe.buffersDeleted++;
    deleteBuffer.call(this, buffer);
  };
  gl.createTexture = function (): WebGLTexture | null {
    const texture = createTexture.call(this);
    if (texture) probe.texturesCreated++;
    return texture;
  };
  gl.deleteTexture = function (texture): void {
    if (texture) probe.texturesDeleted++;
    deleteTexture.call(this, texture);
  };
  new MutationObserver(() => {
    if (document.querySelector('#engineStatus')?.textContent?.includes('prêt') &&
        document.querySelector<HTMLElement>('#generationStatus')?.dataset.visible === 'false') {
      probe.readyAtMs = performance.now();
    }
  }).observe(document, { subtree: true, childList: true, attributes: true, characterData: true });
}

async function ready(page: Page): Promise<void> {
  await expect(page.locator('#generationStatus')).toHaveAttribute('data-visible', 'false', { timeout: 60_000 });
  await expect(page.locator('#engineStatus')).toContainText('prêt');
  await page.evaluate(() => new Promise<void>((resolve) => {
    requestAnimationFrame(() => requestAnimationFrame(() => resolve()));
  }));
  const active = await page.evaluate(() => ({
    workers: window.__qualificationProbe.liveWorkers,
    pending: window.__qualificationProbe.jobs.filter((job) => job.status === 'pending').length,
  }));
  expect(active).toEqual({ workers: 1, pending: 0 });
}

async function regenerate(page: Page): Promise<ObservedJob & { readyMs: number; receiveToReadyMs: number }> {
  await page.locator('#generateButton').click();
  await ready(page);
  return page.evaluate(() => {
    const probe = window.__qualificationProbe;
    const job = probe.jobs.at(-1);
    if (!job || job.status !== 'result' || job.finishedMs === null || probe.readyAtMs === null) {
      throw new Error('Missing completed generation measurement');
    }
    return { ...job, readyMs: probe.readyAtMs - job.startedMs, receiveToReadyMs: probe.readyAtMs - job.finishedMs };
  });
}

function quantiles(values: number[]): { median: number; p95: number; minimum: number; maximum: number } {
  const sorted = [...values].sort((a, b) => a - b);
  return {
    median: sorted[Math.floor(sorted.length / 2)] ?? 0,
    p95: sorted[Math.min(sorted.length - 1, Math.floor(sorted.length * 0.95))] ?? 0,
    minimum: sorted[0] ?? 0,
    maximum: sorted.at(-1) ?? 0,
  };
}

test('qualifie génération, CPU, cadence et ressources sur des cycles réels', async ({ page, browser, context }, testInfo) => {
  test.setTimeout(360_000);
  const errors: string[] = [];
  page.on('pageerror', (error) => errors.push(error.message));
  page.on('console', (message) => { if (message.type() === 'error') errors.push(message.text()); });
  await page.addInitScript(installProbe);
  const cdp = await context.newCDPSession(page);
  await cdp.send('HeapProfiler.enable');
  const response = await page.goto('/?diagnostics=1');
  expect(response?.ok()).toBe(true);
  await ready(page);
  const initialReadyMs = await page.evaluate(() => window.__qualificationProbe.readyAtMs);
  const environment = await page.evaluate(() => {
    const canvas = document.querySelector<HTMLCanvasElement>('#atomSimCanvas');
    const gl = canvas?.getContext('webgl2');
    if (!gl) throw new Error('WebGL2 required for qualification');
    const debug = gl.getExtension('WEBGL_debug_renderer_info');
    return {
      userAgent: navigator.userAgent,
      hardwareConcurrency: navigator.hardwareConcurrency,
      devicePixelRatio: window.devicePixelRatio,
      viewport: { width: innerWidth, height: innerHeight },
      gpu: debug ? String(gl.getParameter(debug.UNMASKED_RENDERER_WEBGL)) : String(gl.getParameter(gl.RENDERER)),
      gpuTimerExtension: gl.getExtension('EXT_disjoint_timer_query_webgl2') !== null,
    };
  });
  const sampleMaximum = Number(await page.locator('#sampleCount').getAttribute('max'));
  expect(sampleMaximum).toBe(60_000);
  const generations = [];
  await cdp.send('Profiler.enable');
  await cdp.send('Profiler.start');
  for (const sampleCount of [2000, 15000, sampleMaximum]) {
    await page.locator('#sampleCount').fill(String(sampleCount));
    await ready(page);
    for (let repetition = 0; repetition < 3; repetition++) {
      generations.push({ repetition, ...await regenerate(page) });
    }
  }
  const { profile } = await cdp.send('Profiler.stop');
  await testInfo.attach('main-thread-cpu.cpuprofile', { body: JSON.stringify(profile), contentType: 'application/json' });

  // Source doubles are synthetic; the runtime writes doubles directly to Float32 sampling arrays.
  // This isolated allocation+conversion experiment is NOT a separately timed production stage.
  const conversion = await page.evaluate(() => [2000, 15000, 60000].map((count) => {
    const source = Float64Array.from({ length: count * 3 }, (_, index) => 100 * Math.sin(index * 0.71));
    const elapsedMs: number[] = [];
    let maximumAbsoluteError = 0;
    for (let repeat = 0; repeat < 20; repeat++) {
      const start = performance.now();
      const converted = new Float32Array(source);
      elapsedMs.push(performance.now() - start);
      for (let index = 0; index < source.length; index++) {
        maximumAbsoluteError = Math.max(maximumAbsoluteError, Math.abs((source[index] ?? 0) - (converted[index] ?? 0)));
      }
    }
    return { count, elapsedMs, maximumAbsoluteError, sourceBytes: source.byteLength, destinationBytes: source.length * 4 };
  }));
  const frameIntervals = await page.evaluate(() => new Promise<number[]>((resolve) => {
    const intervals: number[] = [];
    let previous = performance.now();
    const frame = (now: number): void => {
      intervals.push(now - previous);
      previous = now;
      if (intervals.length >= 61) resolve(intervals.slice(1));
      else requestAnimationFrame(frame);
    };
    requestAnimationFrame(frame);
  }));

  const resourceSeries: ({ cycle: number; heapBytes: number; resources: SceneDiagnostics; buffers: number; textures: number })[] = [];
  // Return to precisely the same canonical state, count, seed, and hybrid renderer each cycle.
  // Warm-up precedes comparisons so lazy shader/program creation does not mimic a leak.
  for (let cycle = 0; cycle < 5; cycle++) {
    await page.locator('label[for="basisComplex"]').click();
    await page.locator('#quantumN').selectOption('2');
    await page.locator('#quantumL').selectOption('1');
    await page.locator('#quantumM').selectOption('1');
    await page.locator('#sampleCount').fill(cycle % 2 === 0 ? '2000' : '60000');
    await page.locator('#seedInput').fill(String(700 + cycle));
    await page.locator('#seedInput').press('Tab');
    await ready(page);
    await page.locator('#displayMode').selectOption('isosurface');
    await page.locator('#displayMode').selectOption('cloud');
    await page.locator('label[for="basisReal"]').click();
    await page.locator('#quantumN').selectOption('3');
    await page.locator('#realOrbital').selectOption('d_xy');
    await page.locator('#sampleCount').fill('15000');
    await page.locator('#seedInput').fill('1096044365');
    await page.locator('#seedInput').press('Tab');
    await page.locator('#displayMode').selectOption('hybrid');
    await regenerate(page);
    await cdp.send('HeapProfiler.collectGarbage');
    const heap = await cdp.send('Runtime.getHeapUsage');
    const snapshot = await page.evaluate(() => {
      const host = window as unknown as { __atomsDiagnostics?: () => SceneDiagnostics };
      if (!host.__atomsDiagnostics) throw new Error('Development diagnostics hook absent');
      const probe = window.__qualificationProbe;
      return {
        resources: host.__atomsDiagnostics(),
        buffers: probe.buffersCreated - probe.buffersDeleted,
        textures: probe.texturesCreated - probe.texturesDeleted,
      };
    });
    resourceSeries.push({ cycle, heapBytes: heap.usedSize, ...snapshot });
  }
  // Force overlapping real UI generations within one task, so supersession is actually exercised.
  await page.evaluate(() => {
    const button = document.querySelector<HTMLButtonElement>('#generateButton');
    if (!button) throw new Error('Generation control missing');
    button.click(); button.click(); button.click();
  });
  await ready(page);
  const probe = await page.evaluate(() => window.__qualificationProbe);
  const stableResources = resourceSeries.slice(1);
  const first = stableResources[0];
  if (!first) throw new Error('Missing resource cycles');
  for (const point of stableResources) {
    expect(point.resources).toEqual(first.resources);
    expect(point.buffers).toBe(first.buffers);
    expect(point.textures).toBe(first.textures);
  }
  expect(probe.jobs.filter((job) => job.status === 'result').length).toBeGreaterThanOrEqual(20);
  expect(probe.jobs.filter((job) => job.status === 'cancelled').length).toBeGreaterThanOrEqual(2);
  expect(probe.peakWorkers).toBe(1);
  expect(probe.buffersDeleted).toBeGreaterThan(0);
  expect(probe.jobs.filter((job) => job.status === 'error')).toEqual([]);
  expect(errors).toEqual([]);
  const report = {
    recordedAt: new Date().toISOString(), initialReadyMs,
    environment: { ...environment, browserVersion: browser.version(), node: process.version, platform: platform(), osRelease: release(), cpu: cpus()[0]?.model, logicalCpus: cpus().length, physicalMemoryBytes: totalmem() },
    notes: [
      'Vite development runtime, instrumented Chromium, fixed viewport; wall-clock timings are environment-specific.',
      'Worker stages are main-thread receipt timestamps, not isolated Worker CPU execution time.',
      'CPU profile covers the main thread during the load sweep; Worker execution is represented by stage wall times only.',
      'Frame cadence is requestAnimationFrame scheduling including browser/compositor overhead, not GPU execution time.',
      'Float64 conversion is an isolated synthetic allocation experiment; production sampling already produces Float32 arrays.',
      'Heap values are main-page V8 usedSize after explicit GC; exclude Worker heaps and GPU driver memory.',
      'Plateau asserts equal live Three.js and GL resource counts at the same canonical state after warm-up; bounded stress cannot prove universal absence of leaks.',
    ],
    generations, conversion: conversion.map((entry) => ({ ...entry, summaryMs: quantiles(entry.elapsedMs) })),
    frameCadence: { sampleCount: frameIntervals.length, intervalsMs: quantiles(frameIntervals), observedFps: 1000 / (frameIntervals.reduce((sum, value) => sum + value, 0) / frameIntervals.length) },
    resourceSeries, probe,
  };
  await testInfo.attach('performance-resources.json', { body: JSON.stringify(report, null, 2), contentType: 'application/json' });
});
