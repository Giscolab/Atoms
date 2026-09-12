function triggerDownload(url: string, filename: string): void {
  const anchor = document.createElement('a');
  anchor.href = url;
  anchor.download = filename;
  anchor.rel = 'noopener';
  document.body.appendChild(anchor);
  anchor.click();
  anchor.remove();
}

export function downloadBlob(blob: Blob, filename: string): void {
  const url = URL.createObjectURL(blob);
  try {
    triggerDownload(url, filename);
  } finally {
    window.setTimeout(() => {
      URL.revokeObjectURL(url);
    }, 0);
  }
}

export function downloadText(text: string, filename: string): void {
  downloadBlob(new Blob([text], { type: 'application/json;charset=utf-8' }), filename);
}

export function safeFileStem(value: string): string {
  const normalized = value
    .normalize('NFKD')
    .replace(/[^a-zA-Z0-9_-]+/g, '-')
    .replace(/^-+|-+$/g, '')
    .toLowerCase();
  return normalized || 'orbital';
}
