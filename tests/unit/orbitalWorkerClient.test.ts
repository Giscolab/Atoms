import { describe, expect, it } from 'vitest';

import {
  createOrbitalWorkerClient,
  OrbitalWorkerCancellationError,
  type OrbitalWorkerLike,
} from '../../src/app/orbitalWorkerClient';
import { buildOrbitalChartData } from '../../src/sampling/orbitalCharts';
import { computeOrbitalField } from '../../src/sampling/orbitalField';
import { sampleOrbital } from '../../src/sampling/orbitalSampler';
import {
  orbitalGenerationTransferables,
  ORBITAL_FIELD_VERSION,
  ORBITAL_WORKER_PROTOCOL_VERSION,
  SCIENCE_ENGINE_VERSION,
  type OrbitalGenerationPayload,
  type OrbitalWorkerJob,
  type OrbitalWorkerResponse,
} from '../../src/workers/orbitalSamplingProtocol';

class FakeWorker implements OrbitalWorkerLike {
  onerror: OrbitalWorkerLike['onerror'] = null;
  onmessage: OrbitalWorkerLike['onmessage'] = null;
  readonly messages: OrbitalWorkerJob[] = [];
  terminated = false;

  emit(response: OrbitalWorkerResponse): void {
    this.onmessage?.({ data: response });
  }

  postMessage(message: OrbitalWorkerJob): void {
    this.messages.push(message);
  }

  terminate(): void {
    this.terminated = true;
  }
}

function onePointPayload(seed = 7): OrbitalGenerationPayload {
  const state = { basis: 'complex', n: 1, l: 0, m: 0 } as const;
  const extentBohr = 4;
  return {
    charts: buildOrbitalChartData(state, {
      angularPointCount: 9,
      planeExplorationPointCount: 8,
      radialExtentBohr: extentBohr,
      radialPointCount: 9,
    }),
    field: computeOrbitalField(state, { extentBohr, resolution: 8 }),
    radialCoverageProbability: 0.999,
    sampleSet: sampleOrbital({ sampleCount: 1, seed, state }),
  };
}

function result(jobId: string, seed = 7): OrbitalWorkerResponse {
  const payload = onePointPayload(seed);
  return {
    ...payload,
    fieldVersion: ORBITAL_FIELD_VERSION,
    jobId,
    kind: 'result',
    protocolVersion: ORBITAL_WORKER_PROTOCOL_VERSION,
    samplerVersion: payload.sampleSet.metadata.samplerVersion,
    scienceVersion: SCIENCE_ENGINE_VERSION,
  };
}

describe('client Worker orbital', () => {
  it('transmet le contrat versionné et résout le résultat complet courant', async () => {
    const workers: FakeWorker[] = [];
    const client = createOrbitalWorkerClient(() => {
      const worker = new FakeWorker();
      workers.push(worker);
      return worker;
    });

    const pending = client.generate({
      jobId: 'job-1',
      options: { fieldResolution: 8 },
      sampleCount: 1,
      seed: 7,
      state: { basis: 'complex', n: 1, l: 0, m: 0 },
    });
    const worker = workers[0];
    expect(worker).toBeDefined();
    expect(worker?.messages[0]).toMatchObject({
      jobId: 'job-1',
      kind: 'generate',
      options: { fieldResolution: 8, radialCoverageProbability: 0.999 },
      protocolVersion: ORBITAL_WORKER_PROTOCOL_VERSION,
      scienceVersion: SCIENCE_ENGINE_VERSION,
    });
    worker?.emit(result('job-1'));
    await expect(pending).resolves.toMatchObject({
      field: { resolution: 8 },
      sampleSet: { metadata: { sampleCount: 1 } },
    });
    client.dispose();
  });

  it('termine le job remplacé et ignore son résultat obsolète', async () => {
    const workers: FakeWorker[] = [];
    const client = createOrbitalWorkerClient(() => {
      const worker = new FakeWorker();
      workers.push(worker);
      return worker;
    });
    const first = client.generate({
      jobId: 'ancien',
      options: { fieldResolution: 8 },
      sampleCount: 1,
      seed: 1,
      state: { basis: 'complex', n: 1, l: 0, m: 0 },
    });
    const firstWorker = workers[0];
    const second = client.generate({
      jobId: 'courant',
      options: { fieldResolution: 8 },
      sampleCount: 1,
      seed: 2,
      state: { basis: 'complex', n: 1, l: 0, m: 0 },
    });
    await expect(first).rejects.toMatchObject({
      code: 'superseded',
      name: OrbitalWorkerCancellationError.name,
    });
    expect(firstWorker?.terminated).toBe(true);
    firstWorker?.emit(result('ancien'));

    const secondWorker = workers[1];
    secondWorker?.emit(result('courant', 2));
    await expect(second).resolves.toMatchObject({ sampleSet: { metadata: { seed: 2 } } });
    client.dispose();
  });

  it('filtre les jobId étrangers, relaie la progression et réutilise le Worker inactif', async () => {
    const worker = new FakeWorker();
    const progress: string[] = [];
    const client = createOrbitalWorkerClient(() => worker);
    const first = client.generate(
      {
        jobId: 'job-a',
        options: { fieldResolution: 8 },
        sampleCount: 1,
        seed: 1,
        state: { basis: 'complex', n: 1, l: 0, m: 0 },
      },
      (message) => progress.push(message.stage),
    );
    worker.emit({ completed: 3, jobId: 'étranger', kind: 'progress', stage: 'transfer', total: 4 });
    worker.emit({ completed: 0, jobId: 'job-a', kind: 'progress', stage: 'sampling', total: 4 });
    worker.emit(result('job-a', 1));
    await first;
    expect(progress).toEqual(['sampling']);

    const second = client.generate({
      jobId: 'job-b',
      options: { fieldResolution: 8 },
      sampleCount: 1,
      seed: 2,
      state: { basis: 'complex', n: 1, l: 0, m: 0 },
    });
    expect(worker.messages).toHaveLength(2);
    worker.emit(result('job-b', 2));
    await second;
    client.dispose();
  });

  it('expose tous les quatorze buffers transférables sans doublon', () => {
    const payload = onePointPayload();
    const transferables = orbitalGenerationTransferables(payload);
    expect(transferables).toHaveLength(14);
    expect(new Set(transferables).size).toBe(14);
    expect(transferables).toContain(payload.sampleSet.positionsBohr.buffer);
    expect(transferables).toContain(payload.field.densityNormalized.buffer);
    expect(transferables).toContain(payload.charts.angularCut.normalizedRadius.buffer);
  });

  it('rejette un résultat dont la version de protocole ne correspond pas', async () => {
    const worker = new FakeWorker();
    const client = createOrbitalWorkerClient(() => worker);
    const pending = client.generate({
      jobId: 'version-incompatible',
      options: { fieldResolution: 8 },
      sampleCount: 1,
      seed: 7,
      state: { basis: 'complex', n: 1, l: 0, m: 0 },
    });
    const incompatible = result('version-incompatible');
    if (incompatible.kind !== 'result') throw new Error('Résultat de test inattendu.');
    worker.emit({ ...incompatible, protocolVersion: 'version-inconnue' });
    await expect(pending).rejects.toThrow(/versions ou métadonnées de résultat incompatibles/u);
    expect(worker.terminated).toBe(true);
    client.dispose();
  });
});
