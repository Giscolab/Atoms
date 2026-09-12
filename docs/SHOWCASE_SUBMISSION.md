# OpenAI Showcase — dossier de soumission

Ce document prépare les champs projet du formulaire OpenAI Showcase. Il ne remplace pas les champs personnels, l’attestation de droits ni l’acceptation finale du propriétaire.

## Projet

**Project type**

App — interactive scientific visualization

**Did you use Codex to build this?**

Yes

**Did you use another coding agent to build this?**

Yes — GPT-6 Astra

**Tech stack**

TypeScript, Vite, Three.js/WebGL2, Web Workers, Vitest, Playwright, ESLint, GitHub Actions, GitHub Pages.

**Use cases**

Education, data visualization, interactive experiences, scientific visualization, quantum mechanics.

## Modèles et processus de construction

**Capability showcased**

N/A — the deployed app is a standalone scientific browser experience. OpenAI coding agents were used to build, review and validate it; no OpenAI capability is invoked at runtime.

**OpenAI models and APIs**

Codex and GPT-6 Astra were used during development. The deployed application does not call an OpenAI API at runtime.

**Other models or APIs**

No other AI model or model API is used by the deployed application.

**Building process**

Atoms was developed iteratively with Codex and GPT-6 Astra as engineering agents. They assisted with architecture, TypeScript migration, quantum-model implementation, deterministic sampling, Web Workers, Three.js rendering, accessibility, testing and scientific documentation. I reviewed each stage, corrected results when needed, and required reproducible unit, scientific and multi-browser validation before moving on.

## URLs et exécution

**Public GitHub repository**

https://github.com/Giscolab/Atoms

**Hosted URL**

https://giscolab.github.io/Atoms/

**Setup steps**

Try it directly at https://giscolab.github.io/Atoms/. For local use: `git clone https://github.com/Giscolab/Atoms.git`, `cd Atoms`, `npm ci`, then `npm run dev` and open the Vite URL. Requires Node.js 24.19.0+ and a WebGL2-capable browser.

## Présentation publique

**Title**

Hydrogen Quantum Orbital Visualizer

**Tagline**

Explore hydrogen quantum orbitals in an interactive 3D laboratory with probability density, wavefunction phase, nodal structures and reproducible sampling.

**Description**

Atoms is an interactive browser laboratory for exploring hydrogen quantum orbitals from a tested TypeScript scientific core. It visualizes probability density, wavefunction phase, real-orbital combinations, nodal structures and radial/angular analyses in a responsive Three.js scene.

Scientific sampling runs in Web Workers and is reproducible from an explicit seed. The project includes numerical and statistical validation, deterministic scientific regression tests, accessibility checks and functional qualification across Chromium, Firefox and WebKit.

**Author**

Giscolab

**Public cover image**

https://raw.githubusercontent.com/Giscolab/Atoms/main/docs/captures/Atoms_3d_xy_phase_hybride_n3_l2_mpm2_16000.png

Source locale : `docs/captures/Atoms_3d_xy_phase_hybride_n3_l2_mpm2_16000.png` — 1907 × 914 px.

## À compléter personnellement dans le formulaire

- first name, last name et email ;
- website/social profile si souhaité ;
- confirmer que `Giscolab` est bien le nom d’auteur à afficher ;
- relire l’attestation de droits et le Showcase Gallery Program Agreement ;
- cocher l’attestation et soumettre uniquement après cette vérification personnelle.

Le formulaire officiel accepte explicitement `N/A` pour une capability ou une API absente du runtime. Atoms peut donc décrire honnêtement l’usage de Codex/GPT-6 Astra pendant le développement sans prétendre utiliser une API OpenAI dans l’application déployée.

Source du formulaire vérifiée le 12 septembre 2026 : https://openai.com/form/showcase-submission/
