import {
  ORBITAL_SAMPLER_VERSION,
  type OrbitalSampleSet,
  type OrbitalSamplingState,
} from '../sampling/contracts';
import { ORBITAL_FIELD_VERSION, type OrbitalFieldGrid } from '../sampling/fieldContracts';
import type { OrbitalChartData } from '../sampling/orbitalCharts';

/** Version explicite du noyau de fonctions d'onde consommé par le Worker. */
export const SCIENCE_ENGINE_VERSION = 'hydrogen-wavefunctions-phase3-v1' as const;

/** Version du protocole de génération complet ; indépendante du sampler scientifique. */
export const ORBITAL_WORKER_PROTOCOL_VERSION = 'orbital-worker-protocol-v2' as const;

export interface OrbitalWorkerOptions {
  /** Nombre de sommets par axe de la grille volumique. Paramètre numérique, pas physique. */
  readonly fieldResolution: number;
  /** Quantile de la CDF radiale utilisé comme demi-étendue du cube de rendu. */
  readonly radialCoverageProbability: number;
  readonly angularChartPointCount: number;
  readonly planeExplorationPointCount: number;
  readonly radialChartPointCount: number;
}

export const DEFAULT_ORBITAL_WORKER_OPTIONS: OrbitalWorkerOptions = {
  angularChartPointCount: 361,
  fieldResolution: 32,
  planeExplorationPointCount: 720,
  radialChartPointCount: 301,
  radialCoverageProbability: 0.999,
};

export interface OrbitalWorkerJob {
  readonly jobId: string;
  readonly kind: 'generate';
  readonly options: OrbitalWorkerOptions;
  readonly protocolVersion: string;
  readonly sampleCount: number;
  readonly samplerVersion: string;
  readonly scienceVersion: string;
  readonly seed: number;
  readonly state: OrbitalSamplingState;
}

export type OrbitalWorkerProgressStage = 'sampling' | 'field' | 'charts' | 'transfer';

export interface OrbitalWorkerProgress {
  readonly completed: number;
  readonly jobId: string;
  readonly kind: 'progress';
  readonly stage: OrbitalWorkerProgressStage;
  readonly total: number;
}

export interface OrbitalGenerationPayload {
  readonly charts: OrbitalChartData;
  readonly field: OrbitalFieldGrid;
  readonly radialCoverageProbability: number;
  readonly sampleSet: OrbitalSampleSet;
}

export interface OrbitalWorkerResult extends OrbitalGenerationPayload {
  readonly fieldVersion: string;
  readonly jobId: string;
  readonly kind: 'result';
  readonly protocolVersion: string;
  readonly samplerVersion: string;
  readonly scienceVersion: string;
}

export interface OrbitalWorkerFailure {
  readonly jobId: string;
  readonly kind: 'error';
  readonly message: string;
}

export type OrbitalWorkerResponse =
  OrbitalWorkerFailure | OrbitalWorkerProgress | OrbitalWorkerResult;

function requireJobId(jobId: string): void {
  if (jobId.trim().length === 0) {
    throw new RangeError("L'identifiant du job Worker ne peut pas être vide.");
  }
}

function requireSafePointCount(value: number, quantity: string, minimum: number): void {
  if (!Number.isSafeInteger(value) || value < minimum) {
    throw new RangeError(`${quantity} doit être un entier sûr supérieur ou égal à ${minimum}.`);
  }
}

export function createOrbitalWorkerJob(
  request: Pick<OrbitalWorkerJob, 'jobId' | 'sampleCount' | 'seed' | 'state'> & {
    readonly options?: Partial<OrbitalWorkerOptions>;
  },
): OrbitalWorkerJob {
  requireJobId(request.jobId);
  return {
    jobId: request.jobId,
    kind: 'generate',
    options: { ...DEFAULT_ORBITAL_WORKER_OPTIONS, ...request.options },
    protocolVersion: ORBITAL_WORKER_PROTOCOL_VERSION,
    sampleCount: request.sampleCount,
    samplerVersion: ORBITAL_SAMPLER_VERSION,
    scienceVersion: SCIENCE_ENGINE_VERSION,
    seed: request.seed,
    state: request.state,
  };
}

export function assertCompatibleOrbitalWorkerJob(job: OrbitalWorkerJob): void {
  requireJobId(job.jobId);
  if (job.protocolVersion !== ORBITAL_WORKER_PROTOCOL_VERSION) {
    throw new RangeError(`Version de protocole Worker incompatible : ${job.protocolVersion}.`);
  }
  if (job.scienceVersion !== SCIENCE_ENGINE_VERSION) {
    throw new RangeError(`Version du noyau scientifique incompatible : ${job.scienceVersion}.`);
  }
  if (job.samplerVersion !== ORBITAL_SAMPLER_VERSION) {
    throw new RangeError(`Version du sampler incompatible : ${job.samplerVersion}.`);
  }
  requireSafePointCount(job.options.fieldResolution, 'La résolution du champ', 8);
  if (job.options.fieldResolution > 64) {
    throw new RangeError('La résolution du champ doit rester inférieure ou égale à 64.');
  }
  requireSafePointCount(job.options.angularChartPointCount, 'La résolution angulaire', 3);
  requireSafePointCount(job.options.radialChartPointCount, 'La résolution radiale', 3);
  requireSafePointCount(job.options.planeExplorationPointCount, "La résolution d'exploration", 4);
  if (
    !Number.isFinite(job.options.radialCoverageProbability) ||
    job.options.radialCoverageProbability <= 0 ||
    job.options.radialCoverageProbability >= 1
  ) {
    throw new RangeError('La couverture radiale du rendu doit appartenir strictement à ]0,1[.');
  }
}

function requireArrayBuffer(buffer: ArrayBufferLike): ArrayBuffer {
  if (!(buffer instanceof ArrayBuffer)) {
    throw new TypeError('Le protocole Worker requiert des ArrayBuffer transférables.');
  }
  return buffer;
}

/** Tous les buffers scientifiques sont transférés, jamais partagés entre deux jobs. */
export function orbitalGenerationTransferables(payload: OrbitalGenerationPayload): ArrayBuffer[] {
  const { angularCut, radialDistribution } = payload.charts;
  return [
    payload.sampleSet.positionsBohr.buffer,
    payload.sampleSet.radiiBohr.buffer,
    payload.sampleSet.thetaRadians.buffer,
    payload.sampleSet.phiRadians.buffer,
    payload.sampleSet.phaseRadians.buffer,
    payload.field.densityNormalized.buffer,
    payload.field.phaseRadians.buffer,
    payload.field.signedAmplitudeNormalized.buffer,
    radialDistribution.radiiBohr.buffer,
    radialDistribution.rawDensityInverseBohr.buffer,
    radialDistribution.normalizedDensity.buffer,
    angularCut.anglesRadians.buffer,
    angularCut.rawDensityInverseSteradians.buffer,
    angularCut.normalizedRadius.buffer,
  ].map(requireArrayBuffer);
}

export { ORBITAL_FIELD_VERSION };
