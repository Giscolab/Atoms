import { describeOrbital, PANEL_MAX_N, PANEL_MIN_N } from '../data/orbitals';
import {
  applyAtomValues,
  clampValue,
  normalizeQuantumState,
  type AppState,
} from '../state/appState';
import { requireElement } from './dom';
import { chooseTheme, inclusiveIntegerRange, parseOrbitalOption, type Theme } from './uiState';

const THEME_STORAGE_KEY = 'atom-theme';
const THEME_CLEAR_COLOR: Record<Theme, number> = { dark: 0x050302, light: 0xf3ece4 };

export interface UiCallbacks {
  generateCloud(): void;
  refreshHud(): void;
  setPointSize(pointSize: number): void;
  setThemeColor(color: number): void;
  loadWavefunction(file: File): Promise<void>;
  initialize2D(): void;
}

export interface AppUi {
  bind(callbacks: UiCallbacks): void;
  showLoading(message?: string): void;
  hideLoading(): void;
  updateProgress(done: number, total: number): void;
  updateHud(energyEv: number, cameraDistance: number): void;
  updateParticleCount(count: number): void;
  updateDistance(distance: number): void;
  updateFps(fps: number): void;
}

function optionText(option: HTMLOptionElement | undefined): string {
  if (!option) return '';
  return option.textContent.replace(/\s+/g, ' ').trim();
}

function setSelectOptions(
  select: HTMLSelectElement,
  values: number[],
  selectedValue: number,
): void {
  select.innerHTML = '';
  for (const value of values) {
    const option = document.createElement('option');
    option.value = String(value);
    option.textContent = String(value);
    select.appendChild(option);
  }
  select.value = String(selectedValue);
}

export function createAppUi(state: AppState): AppUi {
  const orbitalSelect = requireElement('orbitalSel', HTMLSelectElement);
  const particleInput = requireElement('nParts', HTMLInputElement);
  const pointSizeInput = requireElement('ptSize', HTMLInputElement);
  const animationToggle = requireElement('animToggle', HTMLInputElement);
  const themeToggle = requireElement('themeToggle', HTMLInputElement);
  const jsonInput = requireElement('jsonInput', HTMLInputElement);
  const colorToggle = requireElement('colorToggle', HTMLInputElement);
  const button2d = requireElement('btn2d', HTMLButtonElement);

  function getPreferredTheme(): Theme {
    let storedTheme: string | null = null;
    try {
      storedTheme = localStorage.getItem(THEME_STORAGE_KEY);
    } catch {}
    return chooseTheme(storedTheme, window.matchMedia('(prefers-color-scheme: light)').matches);
  }

  function applyTheme(theme: Theme, callbacks: UiCallbacks, persist = true): void {
    const normalizedTheme: Theme = theme === 'light' ? 'light' : 'dark';
    document.documentElement.setAttribute('data-theme', normalizedTheme);
    requireElement('themeToggle', HTMLInputElement).checked = normalizedTheme === 'light';
    callbacks.setThemeColor(THEME_CLEAR_COLOR[normalizedTheme]);
    if (persist) {
      try {
        localStorage.setItem(THEME_STORAGE_KEY, normalizedTheme);
      } catch {}
    }
  }

  function refreshAtomForm(fromState = true): void {
    const nSelect = requireElement('atomN', HTMLSelectElement);
    const lSelect = requireElement('atomL', HTMLSelectElement);
    const mSelect = requireElement('atomM', HTMLSelectElement);
    const quantum = state.quantum;

    const requestedN = fromState ? quantum.n : Number(nSelect.value || quantum.n);
    const maxN = Math.max(PANEL_MAX_N, quantum.n);
    const nValue = clampValue(requestedN, PANEL_MIN_N, maxN);
    setSelectOptions(nSelect, inclusiveIntegerRange(PANEL_MIN_N, maxN), nValue);

    const requestedL = fromState ? quantum.l : Number(lSelect.value || quantum.l);
    const lMax = Math.max(0, nValue - 1);
    const lValue = clampValue(requestedL, 0, lMax);
    setSelectOptions(lSelect, inclusiveIntegerRange(0, lMax), lValue);

    const requestedM = fromState ? quantum.m : Number(mSelect.value || quantum.m);
    const mValue = clampValue(requestedM, -lValue, lValue);
    setSelectOptions(mSelect, inclusiveIntegerRange(-lValue, lValue), mValue);
  }

  function updatePresetInfo(): void {
    const { n, l, m } = state.quantum;
    const description = describeOrbital(n, l, m);

    const select = requireElement('orbitalSel', HTMLSelectElement);
    const matchValue = `${n}_${l}_${m}`;
    const selected = [...select.options].find((option) => option.value === matchValue);
    const label = selected ? optionText(selected) : `${n}${description.letter} (m=${m})`;
    const sourceTag = selected ? 'Prereglage du catalogue' : "Etat defini via l'onglet Atome";

    requireElement('presetName', HTMLElement).textContent = label;
    requireElement('presetSummary', HTMLElement).textContent =
      `${description.familySummary} ${sourceTag}.`;
    requireElement('presetFamily', HTMLElement).textContent = description.letter;
    requireElement('presetNumbers', HTMLElement).textContent = `n=${n}, l=${l}, m=${m}`;
    requireElement('presetNodes', HTMLElement).textContent =
      `${description.radialNodes} radial, ${description.angularNodes} angulaire`;
    requireElement('presetProjection', HTMLElement).textContent =
      description.magneticProjectionSummary;
    requireElement('presetLock', HTMLElement).textContent =
      'Contraintes fixes: n>=1, 0<=l<=n-1, -l<=m<=l.';
  }

  function setPanelTab(tabId: 'presets' | 'atom'): void {
    const showPresets = tabId === 'presets';
    const presets = requireElement('panel-view-presets', HTMLElement);
    const atom = requireElement('panel-view-atom', HTMLElement);

    presets.classList.toggle('active', showPresets);
    atom.classList.toggle('active', !showPresets);
    presets.hidden = !showPresets;
    atom.hidden = showPresets;

    requireElement('tabBtnPresets', HTMLButtonElement).setAttribute(
      'aria-selected',
      showPresets ? 'true' : 'false',
    );
    requireElement('tabBtnAtom', HTMLButtonElement).setAttribute(
      'aria-selected',
      showPresets ? 'false' : 'true',
    );
  }

  function applyOrbitalOption(value: string): void {
    const parsed = parseOrbitalOption(value);
    if (!parsed) throw new Error(`Option orbitale invalide : ${value}`);
    [state.quantum.n, state.quantum.l, state.quantum.m] = parsed;
  }

  const api: AppUi = {
    bind(callbacks): void {
      orbitalSelect.addEventListener('change', () => {
        applyOrbitalOption(orbitalSelect.value);
        callbacks.refreshHud();
        callbacks.generateCloud();
      });

      requireElement('btnGen', HTMLButtonElement).addEventListener('click', () => {
        callbacks.generateCloud();
      });
      requireElement('tabBtnPresets', HTMLButtonElement).addEventListener('click', () => {
        setPanelTab('presets');
      });
      requireElement('tabBtnAtom', HTMLButtonElement).addEventListener('click', () => {
        setPanelTab('atom');
      });
      requireElement('atomN', HTMLSelectElement).addEventListener('change', () => {
        refreshAtomForm(false);
      });
      requireElement('atomL', HTMLSelectElement).addEventListener('change', () => {
        refreshAtomForm(false);
      });
      requireElement('atomM', HTMLSelectElement).addEventListener('change', () => {
        refreshAtomForm(false);
      });
      requireElement('btnApplyAtom', HTMLButtonElement).addEventListener('click', () => {
        const nValue = Number(requireElement('atomN', HTMLSelectElement).value);
        const lValue = Number(requireElement('atomL', HTMLSelectElement).value);
        const mValue = Number(requireElement('atomM', HTMLSelectElement).value);
        applyAtomValues(state.quantum, nValue, lValue, mValue);
        callbacks.refreshHud();
        callbacks.generateCloud();
      });
      setPanelTab('presets');
      applyTheme(getPreferredTheme(), callbacks, false);

      requireElement('btnRandom', HTMLButtonElement).addEventListener('click', () => {
        const options = [...orbitalSelect.options];
        const pick = options[Math.floor(Math.random() * options.length)];
        if (!pick) return;
        orbitalSelect.value = pick.value;
        applyOrbitalOption(pick.value);
        callbacks.refreshHud();
        callbacks.generateCloud();
      });

      particleInput.addEventListener('change', () => {
        state.quantum.N = Number(particleInput.value);
        requireElement('nPartsVal', HTMLOutputElement).textContent =
          state.quantum.N.toLocaleString();
        callbacks.generateCloud();
      });
      particleInput.addEventListener('input', () => {
        state.quantum.N = Number(particleInput.value);
        requireElement('nPartsVal', HTMLOutputElement).textContent =
          state.quantum.N.toLocaleString();
      });

      pointSizeInput.addEventListener('input', () => {
        state.pointSize = Number(pointSizeInput.value) * 0.08;
        requireElement('ptSizeVal', HTMLOutputElement).textContent = pointSizeInput.value;
        callbacks.setPointSize(state.pointSize);
      });

      animationToggle.addEventListener('change', () => {
        state.animating = animationToggle.checked;
      });

      themeToggle.addEventListener('change', () => {
        applyTheme(themeToggle.checked ? 'light' : 'dark', callbacks);
      });

      jsonInput.addEventListener('change', () => {
        const file = jsonInput.files?.[0];
        if (!file) return;
        api.showLoading('Chargement JSON...');
        void callbacks.loadWavefunction(file).catch((error: unknown) => {
          api.hideLoading();
          const message = error instanceof Error ? error.message : String(error);
          alert(`Erreur JSON : ${message}`);
        });
        jsonInput.value = '';
      });

      colorToggle.addEventListener('change', () => {
        state.useAltColor = colorToggle.checked;
        callbacks.generateCloud();
      });

      button2d.addEventListener('click', () => {
        state.show2D = !state.show2D;
        requireElement('panel2d', HTMLElement).classList.toggle('visible', state.show2D);
        button2d.classList.toggle('active', state.show2D);
        if (state.show2D && !state.simulation2DInitialized) {
          state.simulation2DInitialized = true;
          callbacks.initialize2D();
        }
      });

      window.addEventListener('keydown', (event) => {
        if (event.repeat) return;
        const key = event.key.toLowerCase();
        let dirty = false;
        if (key === 'w') {
          state.quantum.n++;
          dirty = true;
        } else if (key === 's') {
          if (state.quantum.n > 1) state.quantum.n--;
          dirty = true;
        } else if (key === 'e') {
          state.quantum.l++;
          dirty = true;
        } else if (key === 'd') {
          if (state.quantum.l > 0) state.quantum.l--;
          dirty = true;
        } else if (key === 'r') {
          state.quantum.m++;
          dirty = true;
        } else if (key === 'f') {
          state.quantum.m--;
          dirty = true;
        } else if (key === 't') {
          state.quantum.N = Math.min(state.quantum.N * 2, 80000);
          particleInput.value = String(state.quantum.N);
          requireElement('nPartsVal', HTMLOutputElement).textContent =
            state.quantum.N.toLocaleString();
          dirty = true;
        } else if (key === 'g') {
          state.quantum.N = Math.max(Math.round(state.quantum.N / 2), 2000);
          particleInput.value = String(state.quantum.N);
          requireElement('nPartsVal', HTMLOutputElement).textContent =
            state.quantum.N.toLocaleString();
          dirty = true;
        } else if (key === 'a') {
          state.animating = !state.animating;
          animationToggle.checked = state.animating;
        } else if (key === 'q') {
          button2d.click();
        }
        if (dirty) {
          normalizeQuantumState(state.quantum);
          callbacks.generateCloud();
        }
      });
    },
    showLoading(message): void {
      const loading = requireElement('loading', HTMLElement);
      loading.style.display = 'block';
      loading.textContent = message || "Calcul de l'orbitale...";
    },
    hideLoading(): void {
      requireElement('loading', HTMLElement).style.display = 'none';
      requireElement('progressBar', HTMLElement).textContent = '';
    },
    updateProgress(done, total): void {
      requireElement('progressBar', HTMLElement).textContent =
        `${Math.round((done / total) * 100)}% - ${done.toLocaleString()} / ${total.toLocaleString()}`;
    },
    updateHud(energyEv, cameraDistance): void {
      const { n, l, m, N } = state.quantum;
      requireElement('qn-n', HTMLElement).textContent = String(n);
      requireElement('qn-l', HTMLElement).textContent = String(l);
      requireElement('qn-m', HTMLElement).textContent = String(m);
      requireElement('energy-val', HTMLElement).textContent = `${energyEv.toFixed(3)} eV`;
      requireElement('iOrb', HTMLElement).textContent = `${n}, ${l}, ${m}`;
      requireElement('iPts', HTMLElement).textContent = N.toLocaleString();
      requireElement('iDist', HTMLElement).textContent = cameraDistance.toFixed(1);
      const select = requireElement('orbitalSel', HTMLSelectElement);
      const match = `${n}_${l}_${m}`;
      const hasMatch = [...select.options].some((option) => option.value === match);
      const customOption = select.querySelector('option[data-custom-state="true"]');
      if (customOption) customOption.remove();
      if (hasMatch) {
        select.value = match;
      } else {
        const custom = document.createElement('option');
        custom.value = match;
        custom.textContent = `Etat perso (${n}, ${l}, ${m})`;
        custom.dataset.customState = 'true';
        select.prepend(custom);
        select.value = match;
      }
      refreshAtomForm(true);
      updatePresetInfo();
    },
    updateParticleCount(count): void {
      requireElement('nPartsVal', HTMLOutputElement).textContent = count.toLocaleString();
    },
    updateDistance(distance): void {
      requireElement('iDist', HTMLElement).textContent = distance.toFixed(1);
    },
    updateFps(fps): void {
      requireElement('iFps', HTMLElement).textContent = String(fps);
    },
  };

  return api;
}
