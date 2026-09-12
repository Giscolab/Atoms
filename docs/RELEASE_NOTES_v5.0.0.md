# Atoms v5.0.0 — Hydrogen Quantum Orbital Visualizer

Atoms v5.0.0 is the first release of the rebuilt scientific 3D orbital pipeline.
It provides an interactive browser laboratory for hydrogen quantum orbitals, backed by a tested TypeScript scientific core and a reproducible rendering pipeline.

## Highlights

- Interactive 3D visualization of hydrogen quantum orbitals with Three.js/WebGL2.
- Complex `|n,l,m⟩` states and explicit normalized real `p`/`d` orbital combinations.
- Probability-density, wavefunction-phase and hybrid display modes.
- Deterministic Monte Carlo sampling from `|ψ|² dV` with an explicit `uint32` seed.
- Web Worker generation with transferable buffers and clean supersession of stale jobs.
- Density isosurfaces and `ψ = 0` nodal surfaces where a real field is physically interpretable.
- Radial and angular analysis panels, reduced-mass hydrogen energy and node counts.
- Responsive dark/light interface with keyboard navigation and reduced-motion support.

## Scientific model

The implemented model is neutral `¹H` in the non-relativistic Coulomb Schrödinger approximation using the electron–proton reduced mass.
CODATA 2022 / NIST constants and NIST DLMF conventions are used for the numerical core.

The release deliberately excludes fine/hyperfine structure, spin–orbit effects, relativistic recoil, finite proton size, Lamb shift, QED corrections, external fields and experimental spectroscopy.

## Validation

The release is qualified by 293 Vitest tests across 27 files, covering quantum-number validation, CODATA constants, units, energy, special functions, radial functions, spherical harmonics, real orbitals, wavefunctions, observables, sampling and Worker behavior.

Scientific regression coverage includes seven deterministic SVG baselines: `1s`, `2s`, real `2p_z`, complex `2p m=+1`, real `3d_z²`, real `3d_xy` and complex `4f`.

Browser qualification covers Chromium, Firefox and WebKit. Accessibility checks include keyboard focus behavior in both themes and bases, plus automated WCAG A/AA axe-core scans under Chromium.

A dedicated qualification suite also stresses generation supersession, Worker concurrency, JavaScript heap behavior and Three.js/WebGL resource stability. Performance figures are environment-dependent measurements, not hardware guarantees.

## Engineering and delivery

- TypeScript strict mode with additional unchecked-index, exact-optional and unused-symbol checks.
- Scientific core typechecked independently from DOM libraries.
- Vite production build and GitHub Pages deployment.
- GitHub Actions gates deployment on dependency audit, formatting, lint, typechecks, unit/scientific tests, production build and multi-browser E2E qualification.
- Reproducible dependency installation from `package-lock.json`.
- Public scientific contract, institutional references and reproducible validation report in `docs/`.

## Known limits

The point cloud is a finite Monte Carlo representation of `|ψ|² dV`; its points are not individual electrons or trajectories. Isosurfaces are finite-grid visualizations and do not establish exact enclosed probability or convergence for every state. Phase colors encode `arg(ψ)` and do not represent charge or multiple particles.

## Requirements

- Node.js `>=24.19.0` for local development.
- A modern browser with WebGL2 support.

Online demo: https://giscolab.github.io/Atoms/

Repository: https://github.com/Giscolab/Atoms

Scientific scope: `docs/SCIENCE.md`

Validation protocol and limits: `docs/VALIDATION.md`

Institutional sources: `docs/REFERENCES.md`

## License

Atoms v5.0.0 is released under the [MIT License](../LICENSE).
