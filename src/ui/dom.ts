export type ElementConstructor<T extends HTMLElement> = new () => T;

export function requireElement<T extends HTMLElement>(
  id: string,
  constructor: ElementConstructor<T>,
): T {
  const element = document.getElementById(id);
  if (!(element instanceof constructor)) {
    throw new Error(`Élément #${id} introuvable ou de type inattendu.`);
  }
  return element;
}
