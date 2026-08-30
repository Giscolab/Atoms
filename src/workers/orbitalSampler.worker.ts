import { buildOrbitalChartData } from '../sampling/orbitalCharts';
import { computeOrbitalField, orbitalDegree } from '../sampling/orbitalField';
import { sampleOrbital } from '../sampling/orbitalSampler';
import { getRadialCdf, sampleRadiusBohr } from '../sampling/radialSampler';
import {
  assertCompatibleOrbitalWorkerJob,
  orbitalGenerationTransferables,
  ORBITAL_FIELD_VERSION,
  ORBITAL_WORKER_PROTOCOL_VERSION,
  SCIENCE_ENGINE_VERSION,
  type OrbitalWorkerJob,
  type OrbitalWorkerProgressStage,
  type OrbitalWorkerResponse,
} from './orbitalSamplingProtocol';

interface WorkerMessageEvent<T> {
  readonly data: T;
}

interface OrbitalWorkerScope {
  onmessage: ((event: WorkerMessageEvent<OrbitalWorkerJob>) => void) | null;
  postMessage(message: OrbitalWorkerResponse, transfer?: readonly ArrayBuffer[]): void;
}

const workerScope = globalThis as unknown as OrbitalWorkerScope;

function errorMessage(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}

function reportProgress(
  job: OrbitalWorkerJob,
  stage: OrbitalWorkerProgressStage,
  completed: number,
): void {
  workerScope.postMessage({
    completed,
    jobId: job.jobId,
    kind: 'progress',
    stage,
    total: 4,
  });
}

workerScope.onmessage = (event): void => {
  const job = event.data;
  try {
    assertCompatibleOrbitalWorkerJob(job);
    reportProgress(job, 'sampling', 0);
    const sampleSet = sampleOrbital({
      sampleCount: job.sampleCount,
      seed: job.seed,
      state: job.state,
    });

    reportProgress(job, 'field', 1);
    const degree = orbitalDegree(job.state);
    const radialCdf = getRadialCdf(job.state.n, degree);
    const extentBohr = sampleRadiusBohr(radialCdf, job.options.radialCoverageProbability);
    const field = computeOrbitalField(job.state, {
      extentBohr,
      resolution: job.options.fieldResolution,
    });

    reportProgress(job, 'charts', 2);
    const charts = buildOrbitalChartData(job.state, {
      angularPointCount: job.options.angularChartPointCount,
      planeExplorationPointCount: job.options.planeExplorationPointCount,
      radialExtentBohr: extentBohr,
      radialPointCount: job.options.radialChartPointCount,
    });

    reportProgress(job, 'transfer', 4);
    const payload = {
      charts,
      field,
      radialCoverageProbability: job.options.radialCoverageProbability,
      sampleSet,
    };
    workerScope.postMessage(
      {
        ...payload,
        fieldVersion: ORBITAL_FIELD_VERSION,
        jobId: job.jobId,
        kind: 'result',
        protocolVersion: ORBITAL_WORKER_PROTOCOL_VERSION,
        samplerVersion: sampleSet.metadata.samplerVersion,
        scienceVersion: SCIENCE_ENGINE_VERSION,
      },
      orbitalGenerationTransferables(payload),
    );
  } catch (error) {
    workerScope.postMessage({
      jobId: job.jobId,
      kind: 'error',
      message: errorMessage(error),
    });
  }
};
