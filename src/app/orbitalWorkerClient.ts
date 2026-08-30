import {
  ORBITAL_SAMPLER_VERSION,
  type OrbitalSamplingRequest,
  type OrbitalSamplingState,
} from '../sampling/contracts';
import { ORBITAL_FIELD_VERSION } from '../sampling/fieldContracts';
import {
  createOrbitalWorkerJob,
  ORBITAL_WORKER_PROTOCOL_VERSION,
  SCIENCE_ENGINE_VERSION,
  type OrbitalGenerationPayload,
  type OrbitalWorkerOptions,
  type OrbitalWorkerProgress,
  type OrbitalWorkerResponse,
} from '../workers/orbitalSamplingProtocol';

export interface WorkerMessageEventLike<T> {
  readonly data: T;
}

export interface WorkerErrorEventLike {
  readonly message: string;
}

export interface OrbitalWorkerLike {
  onerror: ((event: WorkerErrorEventLike) => void) | null;
  onmessage: ((event: WorkerMessageEventLike<OrbitalWorkerResponse>) => void) | null;
  postMessage(message: ReturnType<typeof createOrbitalWorkerJob>): void;
  terminate(): void;
}

export type OrbitalWorkerFactory = () => OrbitalWorkerLike;

export interface OrbitalWorkerClientRequest extends OrbitalSamplingRequest {
  readonly jobId: string;
  readonly options?: Partial<OrbitalWorkerOptions>;
}

export interface OrbitalWorkerClient {
  dispose(): void;
  generate(
    request: OrbitalWorkerClientRequest,
    onProgress?: (progress: OrbitalWorkerProgress) => void,
  ): Promise<OrbitalGenerationPayload>;
}

export class OrbitalWorkerCancellationError extends Error {
  readonly code: 'disposed' | 'superseded';

  constructor(code: 'disposed' | 'superseded') {
    super(
      code === 'superseded'
        ? 'Le job orbital a été remplacé par une génération plus récente.'
        : 'Le client Worker orbital a été arrêté.',
    );
    this.name = 'OrbitalWorkerCancellationError';
    this.code = code;
  }
}

interface ActiveJob {
  readonly expectedFieldResolution: number;
  readonly expectedRadialCoverageProbability: number;
  readonly expectedSampleCount: number;
  readonly expectedSeed: number;
  readonly expectedState: OrbitalSamplingState;
  readonly jobId: string;
  readonly reject: (reason: Error) => void;
  readonly resolve: (payload: OrbitalGenerationPayload) => void;
  readonly worker: OrbitalWorkerLike;
}

function createBrowserOrbitalWorker(): OrbitalWorkerLike {
  const worker = new Worker(new URL('../workers/orbitalSampler.worker.ts', import.meta.url), {
    name: 'atoms-orbital-sampler',
    type: 'module',
  });
  // L'adaptateur applicatif resserre MessageEvent<unknown> vers le protocole validé ci-dessous.
  return worker as unknown as OrbitalWorkerLike;
}

function workerFailure(message: string): Error {
  return new Error(`Échec du Worker orbital : ${message}`);
}

function sameOrbitalState(left: OrbitalSamplingState, right: OrbitalSamplingState): boolean {
  if (left.basis !== right.basis || left.n !== right.n) return false;
  if (left.basis === 'complex' && right.basis === 'complex') {
    return left.l === right.l && left.m === right.m;
  }
  return left.basis === 'real' && right.basis === 'real' && left.orbital === right.orbital;
}

function hasCompatiblePayload(
  response: Extract<OrbitalWorkerResponse, { kind: 'result' }>,
  active: ActiveJob,
): boolean {
  return (
    response.protocolVersion === ORBITAL_WORKER_PROTOCOL_VERSION &&
    response.scienceVersion === SCIENCE_ENGINE_VERSION &&
    response.samplerVersion === ORBITAL_SAMPLER_VERSION &&
    response.fieldVersion === ORBITAL_FIELD_VERSION &&
    response.radialCoverageProbability === active.expectedRadialCoverageProbability &&
    response.field.resolution === active.expectedFieldResolution &&
    response.sampleSet.metadata.seed === active.expectedSeed &&
    response.sampleSet.metadata.sampleCount === active.expectedSampleCount &&
    sameOrbitalState(response.sampleSet.metadata.state, active.expectedState)
  );
}

export function createOrbitalWorkerClient(
  workerFactory: OrbitalWorkerFactory = createBrowserOrbitalWorker,
): OrbitalWorkerClient {
  let active: ActiveJob | null = null;
  let idleWorker: OrbitalWorkerLike | null = null;
  let disposed = false;

  function terminateActive(code: 'disposed' | 'superseded'): void {
    if (!active) return;
    const previous = active;
    active = null;
    previous.worker.terminate();
    previous.reject(new OrbitalWorkerCancellationError(code));
  }

  function failActive(worker: OrbitalWorkerLike, error: Error): void {
    if (active?.worker !== worker) return;
    const previous = active;
    active = null;
    idleWorker = null;
    worker.terminate();
    previous.reject(error);
  }

  return {
    dispose(): void {
      if (disposed) return;
      disposed = true;
      terminateActive('disposed');
      if (idleWorker) {
        idleWorker.terminate();
        idleWorker = null;
      }
    },

    generate(request, onProgress): Promise<OrbitalGenerationPayload> {
      if (disposed) return Promise.reject(new OrbitalWorkerCancellationError('disposed'));
      if (active) terminateActive('superseded');

      const worker = idleWorker ?? workerFactory();
      idleWorker = null;
      const job = createOrbitalWorkerJob(request);

      return new Promise<OrbitalGenerationPayload>((resolve, reject) => {
        active = {
          expectedFieldResolution: job.options.fieldResolution,
          expectedRadialCoverageProbability: job.options.radialCoverageProbability,
          expectedSampleCount: job.sampleCount,
          expectedSeed: job.seed,
          expectedState: job.state,
          jobId: job.jobId,
          reject,
          resolve,
          worker,
        };

        worker.onmessage = (event): void => {
          if (active?.worker !== worker) return;
          const response = event.data;
          if (response.jobId !== active.jobId) return;

          if (response.kind === 'progress') {
            onProgress?.(response);
            return;
          }
          if (response.kind === 'error') {
            failActive(worker, workerFailure(response.message));
            return;
          }
          if (!hasCompatiblePayload(response, active)) {
            failActive(worker, workerFailure('versions ou métadonnées de résultat incompatibles'));
            return;
          }

          const completed = active;
          active = null;
          idleWorker = worker;
          completed.resolve({
            charts: response.charts,
            field: response.field,
            radialCoverageProbability: response.radialCoverageProbability,
            sampleSet: response.sampleSet,
          });
        };

        worker.onerror = (event): void => {
          failActive(worker, workerFailure(event.message));
        };

        worker.postMessage(job);
      });
    },
  };
}
