import type {
  AngularDensityCutChartData,
  RadialDistributionChartData,
} from '../sampling/orbitalCharts';

const SVG_NS = 'http://www.w3.org/2000/svg';

function svgElement<K extends keyof SVGElementTagNameMap>(tag: K): SVGElementTagNameMap[K] {
  return document.createElementNS(SVG_NS, tag);
}

function finiteValue(values: ArrayLike<number>, index: number): number {
  const value = values[index];
  if (value === undefined || !Number.isFinite(value)) return 0;
  return value;
}

/** Convertit une série normalisée en chemin SVG ; aucune donnée scientifique n'est calculée ici. */
export function buildLinePath(
  values: ArrayLike<number>,
  width: number,
  height: number,
  padding: number,
): string {
  if (values.length < 2) return '';
  const plotWidth = Math.max(1, width - 2 * padding);
  const plotHeight = Math.max(1, height - 2 * padding);
  const commands: string[] = [];
  for (let index = 0; index < values.length; index += 1) {
    const x = padding + (plotWidth * index) / (values.length - 1);
    const y = padding + plotHeight * (1 - Math.max(0, Math.min(1, finiteValue(values, index))));
    commands.push(`${index === 0 ? 'M' : 'L'} ${x.toFixed(2)} ${y.toFixed(2)}`);
  }
  return commands.join(' ');
}

/** Convertit une série polaire normalisée en contour SVG fermé. */
export function buildPolarPath(
  radii: ArrayLike<number>,
  centerX: number,
  centerY: number,
  radius: number,
): string {
  if (radii.length < 2) return '';
  const commands: string[] = [];
  for (let index = 0; index < radii.length; index += 1) {
    const angle = (2 * Math.PI * index) / (radii.length - 1) - Math.PI / 2;
    const distance = radius * Math.max(0, Math.min(1, finiteValue(radii, index)));
    const x = centerX + Math.cos(angle) * distance;
    const y = centerY + Math.sin(angle) * distance;
    commands.push(`${index === 0 ? 'M' : 'L'} ${x.toFixed(2)} ${y.toFixed(2)}`);
  }
  return `${commands.join(' ')} Z`;
}

function appendText(
  svg: SVGSVGElement,
  text: string,
  x: number,
  y: number,
  anchor: 'start' | 'middle' | 'end' = 'middle',
): void {
  const label = svgElement('text');
  label.textContent = text;
  label.setAttribute('x', String(x));
  label.setAttribute('y', String(y));
  label.setAttribute('text-anchor', anchor);
  label.setAttribute('class', 'chart-label');
  svg.appendChild(label);
}

function appendLine(
  svg: SVGSVGElement,
  x1: number,
  y1: number,
  x2: number,
  y2: number,
  className: string,
): void {
  const line = svgElement('line');
  line.setAttribute('x1', String(x1));
  line.setAttribute('y1', String(y1));
  line.setAttribute('x2', String(x2));
  line.setAttribute('y2', String(y2));
  line.setAttribute('class', className);
  svg.appendChild(line);
}

function clearSvg(svg: SVGSVGElement): void {
  while (svg.firstChild) svg.removeChild(svg.firstChild);
}

export function renderRadialChart(svg: SVGSVGElement, data: RadialDistributionChartData): void {
  clearSvg(svg);
  const width = 320;
  const height = 190;
  const padding = 24;
  for (let tick = 0; tick <= 4; tick += 1) {
    const y = padding + ((height - 2 * padding) * tick) / 4;
    appendLine(svg, padding, y, width - padding, y, 'chart-grid');
    appendText(svg, (1 - tick / 4).toFixed(1), padding - 6, y + 3, 'end');
  }
  appendLine(svg, padding, padding, padding, height - padding, 'chart-axis');
  appendLine(svg, padding, height - padding, width - padding, height - padding, 'chart-axis');

  const path = svgElement('path');
  path.setAttribute('d', buildLinePath(data.normalizedDensity, width, height, padding));
  path.setAttribute('class', 'chart-series radial-series');
  path.setAttribute('vector-effect', 'non-scaling-stroke');
  svg.appendChild(path);

  const area = svgElement('path');
  area.setAttribute(
    'd',
    `${buildLinePath(data.normalizedDensity, width, height, padding)} L ${width - padding} ${height - padding} L ${padding} ${height - padding} Z`,
  );
  area.setAttribute('class', 'chart-area radial-area');
  svg.insertBefore(area, path);

  const firstRadius = data.radiiBohr[0] ?? 0;
  const lastRadius = data.radiiBohr[data.radiiBohr.length - 1] ?? 0;
  appendText(svg, firstRadius.toFixed(0), padding, height - 8);
  appendText(svg, lastRadius.toFixed(0), width - padding, height - 8);
}

export function renderAngularChart(svg: SVGSVGElement, data: AngularDensityCutChartData): void {
  clearSvg(svg);
  const width = 320;
  const centerX = width / 2;
  const centerY = 106;
  const radius = 78;

  for (const fraction of [0.33, 0.66, 1]) {
    const circle = svgElement('circle');
    circle.setAttribute('cx', String(centerX));
    circle.setAttribute('cy', String(centerY));
    circle.setAttribute('r', String(radius * fraction));
    circle.setAttribute('class', 'chart-grid polar-grid');
    svg.appendChild(circle);
  }
  appendLine(svg, centerX - radius, centerY, centerX + radius, centerY, 'chart-grid');
  appendLine(svg, centerX, centerY - radius, centerX, centerY + radius, 'chart-grid');

  const path = svgElement('path');
  path.setAttribute('d', buildPolarPath(data.normalizedRadius, centerX, centerY, radius));
  path.setAttribute('class', 'chart-series angular-series');
  path.setAttribute('vector-effect', 'non-scaling-stroke');
  svg.appendChild(path);

  const labels: Array<[string, number, number]> = [
    ['90°', centerX, centerY - radius - 8],
    ['0°', centerX + radius + 10, centerY + 3],
    ['180°', centerX - radius - 10, centerY + 3],
    ['270°', centerX, centerY + radius + 17],
  ];
  for (const [label, x, y] of labels) appendText(svg, label, x, y);
}