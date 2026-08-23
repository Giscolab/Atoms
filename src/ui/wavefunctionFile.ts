export type WavefunctionFileSuccess = (value: unknown) => void;
export type WavefunctionFileFailure = (error: Error) => void;

export function readWavefunctionFile(
  file: File,
  onSuccess: WavefunctionFileSuccess,
  onFailure: WavefunctionFileFailure,
): void {
  const reader = new FileReader();
  reader.onload = () => {
    try {
      if (typeof reader.result !== 'string') {
        throw new Error('Invalid JSON file contents');
      }
      onSuccess(JSON.parse(reader.result) as unknown);
    } catch (error) {
      onFailure(error instanceof Error ? error : new Error(String(error)));
    }
  };
  reader.onerror = () => {
    onFailure(reader.error ?? new Error('Unable to read JSON file'));
  };
  reader.readAsText(file);
}
