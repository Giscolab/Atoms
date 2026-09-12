/** Seed d'interface compatible avec le contrat uint32 du sampler. */
export function randomUint32(): number {
  const values = new Uint32Array(1);
  try {
    globalThis.crypto.getRandomValues(values);
    return values[0] ?? 0;
  } catch {
    // Un contexte local sans Crypto API conserve une seed bornée sur uint32.
  }
  return ((Date.now() >>> 0) ^ 0x41544f4d) >>> 0;
}
