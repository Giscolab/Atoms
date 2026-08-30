# Atoms — Plan de refonte scientifique, technique et graphique

**Document de référence — 23 août 2026**  
**Projet :** `Giscolab/Atoms`  
**Origine :** réimplémentation web issue du travail pédagogique original de Kavan (`kavan010/Atoms`)  
**Cible de refonte :** `v5.0.0`  
**État du plan :** mis à jour après validation des Phases 4 à 7 ; la Phase 8 reste différée

---

## Légende d'avancement

- [x] terminé et validé ;
- [ ] à faire ;
- lorsqu'un élément est partiellement réalisé, le détail est indiqué explicitement dans son texte.

Le présent document décrit l'état réel du projet et sa trajectoire cible. Le code validé fait foi sur l'implémentation actuelle ; ce plan fixe les prochaines décisions et les critères de sortie.

---

# 1. But de la refonte

Atoms doit devenir un visualiseur scientifique de l'atome d'hydrogène qui soit à la fois :

- [ ] scientifiquement défendable de bout en bout ;
- [ ] graphiquement lisible et moderne ;
- [ ] déterministe et reproductible ;
- [ ] testable scientifiquement et logiciellement ;
- [ ] performant sur navigateur moderne ;
- [ ] maintenable ;
- [ ] transparent sur ses hypothèses physiques ;
- [ ] explicite sur ce qui est une représentation pédagogique et ce qui est une grandeur physique ;
- [x] respectueux du travail de l'auteur original par une attribution visible et durable.

La priorité absolue est la **justesse scientifique**. Aucun effet visuel ne doit suggérer un phénomène physique faux.

---

# 2. Principes non négociables

## 2.1 La physique est la source de vérité

- [x] Le nouveau noyau scientifique est séparé du moteur legacy.
- [x] Une visualisation ne doit jamais être conservée uniquement parce qu'elle est jolie.
- [x] Toute grandeur affichée doit avoir une unité, une définition et une provenance.
- [x] Toute approximation importante doit être nommée explicitement dans l'interface ou la documentation scientifique.

## 2.2 Séparer calcul, échantillonnage, rendu et interface

- [x] Le noyau scientifique actuel ne dépend ni du DOM ni de Three.js.
- [x] Le renderer Three.js est isolé du nouveau noyau scientifique.
- [x] L'UI est séparée du renderer et du noyau scientifique.
- [x] Le sampler scientifique doit devenir un module indépendant du renderer.
- [x] La génération asynchrone doit devenir un service Worker indépendant.

Règle cible :

```text
science → sampling → worker → application → rendering/UI
```

Aucune dépendance inverse ne doit introduire de formule scientifique dans le renderer ou l'UI.

## 2.3 Aucune trajectoire classique fictive de l'électron

- [x] Un nuage de points doit représenter des échantillons de `|ψ|²`, pas des électrons individuels en mouvement.
- [x] Un état stationnaire ne doit jamais être représenté comme une planète orbitant autour du noyau.
- [x] Le mouvement legacy actuel des points doit disparaître du rendu scientifique final.
- [ ] Un éventuel courant de probabilité doit être présenté comme **champ de courant**, jamais comme trajectoires individuelles.

## 2.4 États complexes et orbitales réelles restent distincts

Le logiciel doit distinguer :

```text
Base complexe : |n,l,m⟩
Base réelle    : pₓ, pᵧ, dxy, ...
```

- [x] `m=+1` ne doit pas être assimilé à `p_x`.
- [x] `m=-1` ne doit pas être assimilé à `p_y`.
- [x] Les orbitales réelles doivent être construites comme combinaisons linéaires correctement normalisées des états complexes `±m`.
- [x] La phase et la densité doivent rester deux grandeurs distinctes.

## 2.5 Unités cohérentes de bout en bout

- [x] Les conversions physiques du nouveau noyau sont centralisées et testées.
- [x] `a₀ ↔ m ↔ pm ↔ nm` est défini dans le nouveau noyau.
- [x] `E_h ↔ J ↔ eV` est défini dans le nouveau noyau.
- [ ] Le runtime legacy doit être remplacé afin qu'aucune ancienne échelle incohérente ne subsiste.
- [ ] Les imports/exports doivent déclarer explicitement leurs unités.

## 2.6 Toute assertion scientifique importante doit être testée

- [x] Les constantes, unités, nombres quantiques, énergie et fonctions spéciales de base sont testés.
- [ ] Les fonctions d'onde doivent être testées analytiquement et numériquement.
- [ ] Le sampler doit être testé statistiquement.
- [ ] Le renderer scientifique doit disposer de régressions visuelles déterministes.

## 2.7 Politique Git

- [x] `main` est la seule branche permanente de travail.
- [x] Le point legacy est conservé par tag Git.
- [ ] Les futurs jalons `v5` doivent être conservés par tags, pas par branches historiques permanentes.
- [ ] Aucun historique ne doit être réécrit uniquement pour masquer la provenance du projet.

## 2.8 Attribution

- [x] Kavan / `kavan010` reste crédité clairement.
- [x] Le dépôt original est référencé dans `CREDITS.md`.
- [x] Le README reprend cette attribution de manière visible.

---

# 3. État du projet

## 3.1 État de départ legacy

Le point de départ historique était une application web compacte :

```text
Atoms/
├── index.html
├── styles.css
├── script.js
├── screenshot.png
└── README.md
```

Le moteur scientifique, le rendu Three.js, l'UI, la simulation 2D et les entrées utilisateur étaient regroupés principalement dans `script.js`.

## 3.2 État actuel après les Phases 4 à 7

L'application a été migrée vers Vite + TypeScript et les responsabilités ont été séparées.

Structure actuelle pertinente :

```text
src/
├── app/
│   ├── animationLoop.ts
│   ├── orbitalPresentation.ts
│   └── orbitalWorkerClient.ts
├── data/
│   ├── orbitals.ts
│   └── wavefunctionData.ts
├── rendering/
│   ├── legacyPhotoelectricModel.ts
│   ├── photoelectric2d.ts
│   ├── phasePalette.ts
│   ├── renderingContracts.ts
│   └── sceneRenderer.ts
├── sampling/
│   ├── angularSampler.ts
│   ├── cdf.ts
│   ├── contracts.ts
│   ├── fieldContracts.ts
│   ├── orbitalCharts.ts
│   ├── orbitalField.ts
│   ├── orbitalSampler.ts
│   ├── radialSampler.ts
│   └── rng.ts
├── science/
│   ├── legacyScience.ts
│   ├── constants/
│   │   └── codata2022.ts
│   ├── hydrogen/
│   │   └── energy.ts
│   ├── quantum/
│   │   └── quantumNumbers.ts
│   ├── special/
│   │   ├── associatedLegendre.ts
│   │   ├── factorial.ts
│   │   └── generalizedLaguerre.ts
│   └── units/
│       └── atomicUnits.ts
├── state/
│   └── appState.ts
├── ui/
│   ├── appUi.ts
│   ├── chartRenderer.ts
│   ├── dom.ts
│   ├── uiState.ts
│   ├── viewportControls.ts
│   └── wavefunctionFile.ts
├── workers/
│   ├── orbitalSampler.worker.ts
│   └── orbitalSamplingProtocol.ts
└── main.ts
```

Tests actuels pertinents :

```text
tests/
├── e2e/
├── scientific/
└── unit/
```

Le pipeline orbital 3D visible est désormais branché sur `src/science`, `src/sampling`, le Worker et
le renderer. `legacyScience.ts` est limité au panneau photon–hydrogène 2D historique, conservé
séparément en Phase 8.

## 3.3 Défauts legacy encore ouverts

- [ ] mélange historique d'unités dans le module 2D legacy ;
- [ ] constante d'import JSON historique `5.29` incorrecte pour `a₀ → pm` dans le module 2D ;
- [x] samplers spéciaux `2p_x`, `2p_y`, `2p_z`, `3p_z` retirés du pipeline 3D ;
- [x] confusion historique entre `m=±1` et `p_x/p_y` retirée de l'UI 3D ;
- [x] concurrence possible entre plusieurs générations asynchrones résolue pour le pipeline 3D ;
- [x] panneau 2D rendu réellement visible/cachable et couvert par E2E ;
- [x] animation du nuage retirée du rendu 3D ;
- [x] énergie visible du panneau 3D issue du modèle à masse réduite ;
- [ ] module 2D encore nommé/interprété comme « photoélectrique » sans modèle physique correct ;
- [x] README mis à jour pour le périmètre actuellement livré.

---

# 4. Point juridique, attribution et références

## 4.1 Crédit

Le README final devra conserver une section visible de ce type :

> Atoms est une réimplémentation web indépendante inspirée du projet pédagogique original **Hydrogen Quantum Orbital Visualizer** de **Kavan** (`kavan010/Atoms`). Le dépôt original a fourni le point de départ historique et pédagogique du projet.

Lien historique :

- https://github.com/kavan010/Atoms

État :

- [x] `CREDITS.md` créé ;
- [x] provenance historique explicitée ;
- [ ] reprendre cette attribution dans le README final.

## 4.2 Licence

Le dépôt original observé ne comporte pas de fichier `LICENSE` à sa racine.

Conséquences :

- [ ] ne jamais affirmer arbitrairement « MIT » ;
- [ ] clarifier le statut juridique final avant `v5.0.0` ;
- [ ] distinguer crédit, provenance et licence ;
- [ ] vérifier les droits de redistribution de tout document scientifique ou média conservé dans le dépôt public ;
- [ ] préférer des références officielles et des liens publics lorsqu'une redistribution locale n'est pas nécessaire.

## 4.3 Références scientifiques

Sources prioritaires :

1. NIST / CODATA pour les constantes physiques ;
2. NIST DLMF pour les fonctions spéciales et conventions ;
3. MIT OpenCourseWare pour la formulation pédagogique et certaines dérivations ;
4. NIST ASD pour les comparaisons spectroscopiques lorsque pertinentes.

- [x] `docs/SCIENCE.md` initial existe ;
- [ ] créer `docs/REFERENCES.md` final ;
- [ ] créer `docs/VALIDATION.md` final.

---

# 5. Stratégie Git : une seule branche `main`

## 5.1 Sauvegarde legacy

- [x] créer le tag `legacy-web-2026-08-23` ;
- [x] conserver le point de départ récupérable ;
- [ ] pousser/archiver les tags de jalon nécessaires selon la politique de release.

## 5.2 Politique cible

```text
main        ← seule branche permanente
v5.0.0-*    ← tags de versions
```

- [x] `main` est la branche permanente de travail ;
- [ ] ne créer aucune branche historique permanente pour les phases ;
- [ ] utiliser des tags aux jalons importants.

## 5.3 Remote historique

Les historiques `main` et `upstream/main` n'ont pas d'ancêtre commun et ne doivent pas être greffés artificiellement.

- [x] provenance préservée par documentation ;
- [ ] traiter le remote `upstream` lors de la phase finale d'autonomie GitHub ;
- [ ] ne jamais supprimer ni modifier le dépôt original de l'auteur.

## 5.4 Autonomie GitHub

Avant toute éventuelle sortie d'un réseau de forks :

- [ ] vérifier issues ;
- [ ] vérifier pull requests ;
- [ ] vérifier stars/watchers ;
- [ ] vérifier forks enfants ;
- [ ] sauvegarder localement ;
- [ ] documenter les conséquences ;
- [ ] décider explicitement si l'opération est nécessaire.

## 5.5 Pas de réécriture globale

- [x] aucune réécriture globale de l'historique n'a été nécessaire pour la refonte ;
- [ ] conserver ce principe jusqu'à `v5.0.0` sauf raison juridique ou technique exceptionnelle documentée.

---

# 6. Stack cible moderne

Versions de référence retenues pour la refonte :

| Composant | Version cible | Rôle | État |
|---|---:|---|---|
| Node.js | 24.19.0 LTS | baseline CI | [x] |
| Node.js | 26.x Current | développement local autorisé | [x] |
| TypeScript | 7.x natif / compatibilité outillage | langage principal | [x] |
| Vite | 8.2.2 | dev + build | [x] |
| Three.js | 0.185.1 / r185 | rendu 3D | [x] |
| Vitest | 4.1.10 | tests unitaires/scientifiques | [x] |
| Playwright Test | 1.62.1 | E2E | [x] |
| ESLint | 10.7.0 | analyse statique | [x] |
| Prettier | 3.9.0 | formatage | [x] |

## 6.1 Baseline Node

- [x] Node 24 LTS utilisé comme référence CI ;
- [x] Node Current autorisé localement ;
- [ ] ne pas relever la baseline CI sans décision explicite.

## 6.2 UI sans framework lourd

- [x] DOM natif + TypeScript conservé ;
- [ ] ne pas introduire React ou autre framework sans besoin réel démontré.

## 6.3 WebGL / WebGPU

- [x] WebGL2 via Three.js est la baseline actuelle ;
- [ ] conserver une isolation suffisante pour envisager WebGPU plus tard ;
- [ ] ne pas rendre WebGPU obligatoire pour `v5.0.0` sans validation multi-navigateurs.

---

# 7. Architecture cible mise à jour

Le plan initial prévoyait `src/physics/`. L'architecture validée utilise désormais `src/science/`. **Aucun renommage vers `physics/` ne doit être effectué uniquement pour coller à une ancienne version du plan.**

Architecture cible :

```text
Atoms/
├── public/
│   └── assets/
├── src/
│   ├── main.ts
│   ├── app/
│   ├── data/
│   ├── science/
│   │   ├── legacyScience.ts
│   │   ├── constants/
│   │   │   └── codata2022.ts
│   │   ├── units/
│   │   │   └── atomicUnits.ts
│   │   ├── quantum/
│   │   │   └── quantumNumbers.ts
│   │   ├── special/
│   │   │   ├── factorial.ts
│   │   │   ├── generalizedLaguerre.ts
│   │   │   └── associatedLegendre.ts
│   │   ├── math/
│   │   │   └── complex.ts
│   │   └── hydrogen/
│   │       ├── energy.ts
│   │       ├── radialWavefunction.ts
│   │       ├── sphericalHarmonics.ts
│   │       ├── wavefunction.ts
│   │       ├── probabilityDensity.ts
│   │       ├── realOrbitals.ts
│   │       └── observables.ts
│   ├── sampling/
│   │   ├── rng.ts
│   │   ├── cdf.ts
│   │   ├── radialSampler.ts
│   │   ├── angularSampler.ts
│   │   └── orbitalSampler.ts
│   ├── workers/
│   │   ├── orbital.worker.ts
│   │   └── protocol.ts
│   ├── rendering/
│   │   ├── sceneRenderer.ts
│   │   ├── camera.ts
│   │   ├── pointCloud.ts
│   │   ├── isosurface.ts
│   │   ├── axes.ts
│   │   ├── phaseColors.ts
│   │   └── dispose.ts
│   ├── features/
│   │   └── photonHydrogen2d/
│   │       ├── model.ts
│   │       ├── transitions.ts
│   │       ├── photoionization.ts
│   │       └── renderer2d.ts
│   ├── state/
│   └── ui/
│       ├── controls.ts
│       ├── presets.ts
│       ├── scientificInfo.ts
│       ├── loading.ts
│       └── accessibility.ts
├── tests/
│   ├── unit/
│   ├── scientific/
│   ├── sampling/
│   ├── regression/
│   └── e2e/
├── docs/
│   ├── SCIENCE.md
│   ├── VALIDATION.md
│   └── REFERENCES.md
├── CREDITS.md
├── PLAN.md
├── README.md
└── ...
```

L'arborescence cible est une direction, pas une obligation de renommer des modules déjà correctement isolés.

---

# 8. Contrat scientifique du moteur quantique

## 8.1 Domaine initial

`v5.0.0` modélise en priorité :

- [x] atome d'hydrogène neutre `¹H` comme domaine scientifique choisi ;
- [x] un électron ;
- [x] potentiel coulombien central ;
- [x] approximation non relativiste de Schrödinger ;
- [x] correction de masse réduite électron-proton dans le nouveau noyau ;
- [x] fonctions d'onde analytiques complètes dans le périmètre complexe et réel `p`/`d` ;
- [x] échantillonnage de `|ψ|²` ;
- [x] rendu orbital 3D utilisant exclusivement le nouveau moteur.

Le modèle `v5.0.0` ne doit pas être présenté comme incluant par défaut :

- structure fine ;
- Lamb shift ;
- spin-orbite ;
- QED ;
- champ magnétique externe ;
- effet Stark ;
- effet Zeeman ;
- interactions multi-électroniques.

## 8.2 Constantes et métrologie

Source : CODATA 2022 / NIST.

Déjà implémenté :

- [x] `c` ;
- [x] `h` ;
- [x] `ħ` dérivé de `h` ;
- [x] `e` ;
- [x] `m_e` ;
- [x] `m_p` ;
- [x] `m_e/m_p` direct CODATA ;
- [x] `a₀` ;
- [x] `E_h`.

Pour `¹H`, le modèle utilise directement :

```text
r = m_e / m_p
μ / m_e = 1 / (1 + r)
```

et conserve la relation physique générale :

```text
μ = m_e m_p / (m_e + m_p)
```

## 8.3 Rayon caractéristique

Le rayon de Bohr conventionnel et le rayon du problème relatif restent distincts :

```text
a₀   = rayon de Bohr conventionnel
a_μ  = a₀ m_e / μ
a_μ  = a₀ (1 + m_e/m_p)
```

- [x] distinction `a₀` / `a_μ` implémentée ;
- [x] relation testée ;
- [ ] utiliser `a_μ` explicitement dans la future fonction radiale de `¹H`.

## 8.4 Énergie

Le nouveau moteur utilise :

```text
E_n / E_h = -(μ/m_e) / (2 n²)
```

- [x] résultat disponible en Hartree ;
- [x] conversion en joules ;
- [x] conversion en eV ;
- [x] dépendance `1/n²` testée ;
- [x] distinction modèle de Schrödinger / valeur spectroscopique documentée ;
- [ ] remplacer l'énergie visible legacy par cette grandeur lors de l'intégration UI.

## 8.5 Fonction d'onde

Cible :

```text
ψₙₗₘ(r, θ, φ) = Rₙₗ(r) Yₗᵐ(θ, φ)
```

- [ ] implémenter `Rₙₗ` ;
- [ ] implémenter `Yₗᵐ` ;
- [ ] implémenter `ψₙₗₘ` complexe ;
- [ ] implémenter `|ψₙₗₘ|²` ;
- [ ] implémenter la phase ;
- [ ] valider la normalisation complète.

## 8.6 Coordonnée radiale et masse réduite

Le futur moteur radial de `¹H` doit utiliser une variable adimensionnelle cohérente avec la masse réduite :

```text
ρ = 2r / (n a_μ)
```

Les API peuvent exposer les longueurs en `a₀`, mais la formule physique doit convertir explicitement vers `a_μ` lorsque nécessaire.

- [ ] aucun usage silencieux de `a₀` à la place de `a_μ` ;
- [ ] tester les formes fermées `1s`, `2s`, `2p` avec cette convention.

## 8.7 Partie radiale

Le sampler radial futur doit être basé sur :

```text
Pₙₗ(r) dr = r² |Rₙₗ(r)|² dr
```

- [ ] ne jamais oublier le facteur `r²` ;
- [ ] ne jamais l'appliquer deux fois ;
- [ ] normaliser numériquement `∫₀∞ r²|Rₙₗ|²dr = 1`.

## 8.8 Partie angulaire

Pour un état complexe :

```text
Yₗᵐ(θ, φ) ∝ Pₗᵐ(cos θ) e^(imφ)
```

Convention retenue :

- [x] fonction de Ferrers / Legendre associé DLMF ;
- [x] phase de Condon–Shortley déjà incluse dans `Pₗᵐ` ;
- [x] ne jamais appliquer une seconde phase de Condon–Shortley dans `Yₗᵐ` ;
- [x] traiter correctement `m < 0` ;
- [x] vérifier `Y_l^{-m}` via la relation de conjugaison appropriée ;
- [x] vérifier que `|Y_l^{+m}|² = |Y_l^{-m}|²`.

## 8.9 Orbitales réelles

- [x] `p_x`, `p_y`, `p_z` ;
- [x] `d_xy`, `d_xz`, `d_yz` ;
- [x] `d_x²-y²`, `d_z²` ;
- [ ] familles supérieures nécessaires à l'UI.

Elles doivent être dérivées de combinaisons normalisées d'états complexes `±m`, avec convention de signe documentée et tests associés.

## 8.10 Nœuds

Relations attendues :

```text
nœuds radiaux    = n - l - 1
nœuds angulaires = l
nœuds totaux     = n - 1
```

- [x] comptage simple disponible dans les données/UI legacy ;
- [x] le calcul scientifique est devenu source de vérité ;
- [x] tester les nœuds radiaux ;
- [x] tester les nœuds angulaires ;
- [x] afficher les surfaces nodales lors de la phase renderer/UI.

## 8.11 Unités

Convention cible :

### Interne scientifique

- longueur d'API : `a₀` lorsque pertinent ;
- longueur caractéristique du problème `¹H` : `a_μ` explicitement dérivée ;
- énergie : Hartree `E_h` ;
- angles : radians ;
- fonctions d'onde : normalisation cohérente avec les unités choisies.

### Affichage

- `a₀` ;
- pm ;
- nm ;
- eV.

- [x] conversions du nouveau noyau centralisées ;
- [ ] éliminer toutes les constantes magiques restantes du module 2D historique lors de son remplacement en Phase 8.

---

# 9. Échantillonnage scientifique : contrat cible

Le sampler scientifique ne doit être construit qu'après validation des fonctions d'onde.

## 9.1 Suppression des exceptions ad hoc

Éléments supprimés lors de l'intégration du sampler Phase 4 :

- [x] `_rProb2p` ;
- [x] `_rProb3p` ;
- [x] `_sampleR2p` ;
- [x] `_sampleR3p` ;
- [x] `NAMED_SAMPLERS` ;
- [x] `SAMPLER_MAP` ;
- [x] toute correspondance `m=±1 → p_x/p_y`.

## 9.2 Pipeline commun

Le sampler doit :

- [x] accepter une définition d'état scientifique explicite ;
- [x] accepter soit une base complexe, soit une orbitale réelle correctement définie ;
- [x] produire des positions suivant la densité cible ;
- [x] utiliser les mêmes unités internes pour toutes les orbitales ;
- [x] être indépendant de Three.js ;
- [x] être indépendant du DOM ;
- [x] être déterministe pour une seed et une version données ;
- [x] exposer des métadonnées de génération utiles aux tests.

## 9.3 PRNG seedable

- [x] choisir un algorithme PRNG explicite et documenté ;
- [x] versionner son contrat de reproductibilité ;
- [x] définir précisément la transformation d'une seed utilisateur vers l'état interne ;
- [x] interdire `Math.random()` dans le nouveau sampler ;
- [x] tester des vecteurs déterministes du PRNG ;
- [x] garantir qu'une même seed produit la même suite dans une même version du moteur.

## 9.4 CDF radiale

La CDF radiale doit être construite depuis :

```text
Pₙₗ(r) = r² |Rₙₗ(r)|²
```

À faire :

- [x] intégration numérique documentée ;
- [x] choix adaptatif/justifié du domaine radial ;
- [x] contrôle de la masse de probabilité tronquée ;
- [x] résolution adaptée à l'état ;
- [x] interpolation lors de l'inversion de CDF ;
- [x] monotonie explicitement vérifiée ;
- [x] normalisation explicite ;
- [x] tests de convergence quand domaine/résolution augmentent ;
- [x] cache indexé par paramètres scientifiques ;
- [x] invalidation/versionnement explicite du cache.

## 9.5 Sampling angulaire

Pour les états complexes :

- [x] dériver la distribution angulaire de `|Yₗᵐ|²` ;
- [x] tenir compte de la mesure `sin θ dθ dφ` ;
- [x] exploiter l'uniformité en `φ` des densités complexes lorsque mathématiquement applicable ;
- [x] tester la symétrie `+m/-m`.

Pour les orbitales réelles :

- [x] dériver la distribution depuis la vraie combinaison réelle ;
- [x] ne jamais réintroduire de sampler nommé ad hoc ;
- [x] vérifier l'orientation attendue sans changer la distribution radiale.

## 9.6 Tests statistiques

Les tests doivent être déterministes grâce à une seed fixe.

- [x] comparer l'ECDF radiale à la CDF théorique ;
- [x] comparer les distributions angulaires aux densités théoriques ;
- [x] tester moyennes et quantiles pertinents ;
- [x] tester la convergence avec `N` croissant ;
- [x] utiliser des seuils statistiques justifiés plutôt que des tolérances arbitraires ;
- [x] éviter les tests flakys dépendant d'un RNG non déterministe.

## 9.7 Critère d'uniformité du pipeline

Toutes les orbitales doivent utiliser :

```text
état scientifique
    ↓
fonction d'onde / densité
    ↓
sampler commun
    ↓
positions en unités scientifiques
```

Aucun nom d'orbitale ne doit sélectionner une formule de sampling codée à part.

---

# 10. Concurrence et Web Worker

## 10.1 Problème historique et état livré

Le runtime legacy pouvait lancer plusieurs générations asynchrones qui manipulaient des états
partagés. Le pipeline 3D livré associe désormais chaque génération à un `jobId`, termine le Worker
obsolète et filtre les résultats avant toute mise à jour de la scène ou de l'interface.

- [x] supprimer toute possibilité de mélange entre deux générations ;
- [x] ne plus réutiliser des buffers globaux comme contrat implicite entre jobs.

## 10.2 Contrat de job

Chaque génération doit posséder au minimum :

```text
jobId
state
basis
seed
sampleCount
options
scienceVersion
samplerVersion
```

Le résultat doit contenir :

```text
jobId
positions
metadata
phase / informations optionnelles si demandées
```

## 10.3 Buffers

- [x] utiliser des `Float32Array`/buffers transférables lorsque pertinent ;
- [x] ne pas partager un tableau mutable entre deux jobs ;
- [ ] mesurer le coût conversion Float64 → Float32 avant rendu.

## 10.4 Annulation

- [x] une nouvelle génération rend l'ancienne obsolète ;
- [x] un résultat obsolète ne modifie jamais scène, HUD ou progression ;
- [x] le Worker obsolète est interrompu par terminaison dès que raisonnablement possible ;
- [x] tester les changements rapides d'état.

---

# 11. Rendu 3D scientifique

## 11.1 Nuage probabiliste

- [x] les points doivent être échantillonnés depuis `|ψ|²` ;
- [x] nombre de points réglable ;
- [x] transparence contrôlée ;
- [x] légende expliquant que les points ne sont pas des électrons individuels ;
- [x] seed affichable pour reproduction.

## 11.2 Isosurfaces

- [x] densité choisie explicitement ;
- [ ] possibilité de plusieurs niveaux ;
- [x] conventions documentées ;
- [ ] validation visuelle sur cas connus.

## 11.3 Surfaces nodales

- [x] mode distinct ;
- [x] surfaces nodales dérivées du moteur scientifique ;
- [x] ne pas les confondre avec une isosurface de densité.

## 11.4 Phase / signe

- [x] palette divergente ;
- [x] légende ;
- [x] phase séparée de densité ;
- [ ] accessible en déficience de perception des couleurs autant que possible.

## 11.5 Échelle

- [x] afficher `a₀` ;
- [x] conversion pm ;
- [ ] conversion nm si utile ;
- [x] axes optionnels ;
- [x] repère d'échelle ;
- [x] cadrage automatique basé sur la distribution radiale réelle.

## 11.6 Noyau

- [x] si le noyau est agrandi pour être visible, afficher :

```text
Noyau représenté de manière schématique — taille non à l'échelle.
```

## 11.7 Animation

- [x] supprimer l'animation legacy des points comme pseudo-orbites ;
- [x] autoriser rotation de caméra ;
- [x] autoriser rotation volontaire de la scène comme outil d'observation ;
- [ ] autoriser animation de phase uniquement si clairement étiquetée ;
- [ ] représenter un éventuel courant de probabilité comme champ, jamais comme trajectoires de points.

## 11.8 Ressources Three.js

- [x] renderer isolé ;
- [x] auditer disposal des géométries, matériaux et textures ;
- [ ] tester absence de fuite lors de changements répétés d'état.

---

# 12. UI/UX scientifique

## 12.1 Direction visuelle

Objectif : interface d'instrument scientifique moderne.

- [x] surface 3D prioritaire ;
- [x] panneau compact ;
- [x] hiérarchie typographique claire ;
- [x] effets uniquement fonctionnels ;
- [x] informations scientifiques proches du contrôle concerné ;
- [x] thèmes sombre/clair équivalents ;
- [x] `prefers-reduced-motion`.

## 12.2 Panneau orbital

Structure cible :

```text
État quantique
  n
  l
  m

Représentation
  base complexe / base réelle
  orbitale réelle éventuelle

Rendu
  nuage
  isosurface
  nœuds
  phase

Échelle
  a₀ / pm / nm
```

- [x] séparer base complexe et base réelle ;
- [x] reconstruire les presets ;
- [x] supprimer les labels scientifiques faux ;
- [x] empêcher les combinaisons invalides.

## 12.3 Informations scientifiques

Afficher :

- [x] état choisi ;
- [x] notation spectroscopique ;
- [x] nombre de nœuds ;
- [x] énergie du modèle ;
- [x] hypothèse de masse réduite ;
- [x] statistiques radiales utiles ;
- [x] signification de la couleur ;
- [x] signification du nuage ;
- [x] seed si mode reproductible.

## 12.4 Accessibilité

- [ ] navigation clavier complète ;
- [x] focus et contrôles clavier de base présents ;
- [ ] audit final focus visible ;
- [ ] labels explicites ;
- [ ] `aria-live` limité aux informations pertinentes ;
- [ ] contrastes vérifiés ;
- [ ] aucun contrôle identifié uniquement par couleur.

---

# 13. Module 2D : interactions photon–hydrogène

## 13.1 État actuel

- [x] panneau 2D ouvrable/fermable ;
- [x] raccourci `Q` fonctionnel ;
- [x] E2E couvre l'état `hidden` et la visibilité ;
- [ ] modèle physique actuel à remplacer.

## 13.2 Direction

Renommer le module final en :

**Interactions photon–hydrogène 2D**

Il doit distinguer :

- [ ] transitions liées-liées ;
- [ ] absorption ;
- [ ] émission lorsque correctement modélisée ;
- [ ] photoionisation ;
- [ ] seuil énergétique selon l'état initial.

## 13.3 Photoionisation

Approximation simple :

```text
hν > |E_bound|
K ≈ hν - |E_bound|
```

- [ ] documenter hypothèses ;
- [ ] utiliser les énergies du nouveau moteur ;
- [ ] tester les seuils.

## 13.4 Effet photoélectrique métallique

Si ce sujet est conservé :

```text
K_max = hν - Φ
```

- [ ] le placer dans un module séparé ;
- [ ] ne jamais confondre fonction travail d'un matériau et énergie d'ionisation de `¹H`.

---

# 14. Tests scientifiques obligatoires

## 14.1 Fondations déjà validées

- [x] constantes CODATA ;
- [x] unités ;
- [x] nombres quantiques ;
- [x] masse réduite ;
- [x] énergie ;
- [x] factorielle ;
- [x] Laguerre généralisé ;
- [x] Legendre associé/Ferrers ;
- [x] phase de Condon–Shortley documentée.

## 14.2 Normalisation radiale

Pour plusieurs `(n,l)` :

```text
∫₀∞ r² |Rₙₗ(r)|² dr ≈ 1
```

- [ ] `1s` ;
- [ ] `2s` ;
- [ ] `2p` ;
- [ ] états supérieurs représentatifs.

## 14.3 Normalisation angulaire

```text
∫ |Yₗᵐ(θ,φ)|² dΩ ≈ 1
```

- [ ] plusieurs `l` ;
- [ ] `m=0` ;
- [ ] `m>0` ;
- [ ] `m<0`.

## 14.4 Normalisation complète

```text
∫ |ψₙₗₘ|² dV ≈ 1
```

- [ ] états complexes ;
- [ ] orbitales réelles.

## 14.5 Orthogonalité

- [ ] états différents en `n` ;
- [ ] états différents en `l` ;
- [ ] états différents en `m` ;
- [ ] orbitales réelles orthogonales pertinentes.

## 14.6 Symétries

- [ ] `1s` isotrope ;
- [ ] `+m` / `-m` : même densité en base complexe ;
- [ ] `p_x`, `p_y`, `p_z` : même distribution radiale ;
- [ ] orientations réelles correctes ;
- [ ] parité `(-1)^l`.

## 14.7 Nœuds

- [ ] nœuds radiaux `n-l-1` ;
- [ ] nœuds angulaires `l` ;
- [ ] nœuds totaux `n-1`.

## 14.8 Observables analytiques

Ajouter progressivement :

- [ ] `⟨r⟩` ;
- [ ] `⟨1/r⟩` ;
- [ ] rayon le plus probable de `1s` ;
- [x] énergie en fonction de `n` dans le modèle choisi ;
- [ ] autres références analytiques utiles et bien sourcées.

## 14.9 Tests statistiques du sampler

- [ ] histogramme/ECDF radial vs théorie ;
- [ ] distribution angulaire vs théorie ;
- [ ] moyenne/quantiles ;
- [ ] convergence avec `N` ;
- [ ] seed fixe ;
- [ ] seuils statistiques justifiés.

## 14.10 Régressions visuelles

Créer des captures reproductibles au minimum pour :

- [ ] `1s` ;
- [ ] `2s` ;
- [ ] `2p_z` réel ;
- [ ] `2p` complexe `m=+1` ;
- [ ] `3d_z²` ;
- [ ] `3d_xy` ;
- [ ] un état `4f` représentatif.

---

# 15. Références scientifiques de validation

## 15.1 NIST / CODATA

- Fundamental Physical Constants: https://physics.nist.gov/constants
- CODATA version history: https://physics.nist.gov/cuu/Reference/versioncon.shtml
- Atomic Spectra Database: https://physics.nist.gov/asd
- Hydrogen energy levels: https://physics.nist.gov/PhysRefData/Handbook/Tables/hydrogentable5.htm

## 15.2 NIST DLMF

- https://dlmf.nist.gov/

Utilisation principale :

- fonctions de Laguerre ;
- Legendre/Ferrers ;
- futures harmoniques sphériques ;
- conventions mathématiques.

## 15.3 MIT OpenCourseWare

Utilisation :

- formulation pédagogique du problème de l'hydrogène ;
- masse réduite ;
- moment angulaire ;
- structure des solutions analytiques.

## 15.4 Documentation interne

- [x] `SCIENCE.md` initial ;
- [ ] chaque nouvelle formule doit indiquer convention, source, unité, domaine et tests associés ;
- [ ] `REFERENCES.md` final ;
- [ ] `VALIDATION.md` final.

---

# 16. Tests logiciels et qualité

## 16.1 Vitest

Catégories :

```text
unit
scientific
sampling
regression
```

- [x] `unit` ;
- [x] `scientific` fondations ;
- [x] `sampling` ;
- [ ] `regression` scientifique/visuelle hors E2E.

## 16.2 Playwright

- [x] démarrage application ;
- [x] WebGL2/HUD sans erreur ;
- [x] contrôles essentiels ;
- [x] panneau 2D visible/caché ;
- [ ] génération concurrente ;
- [x] état scientifique complexe/réel ;
- [x] responsive final ;
- [x] clavier final ;
- [x] thèmes finalisés ;
- [x] absence d'erreurs console sur les scénarios couverts ;
- [ ] captures scientifiques stables.

## 16.3 ESLint

- [x] configuration active ;
- [x] lint vert sur l'état actuel ;
- [ ] conserver l'interdiction des globals implicites ;
- [ ] surveiller promesses oubliées ;
- [ ] surveiller conversions douteuses ;
- [ ] éviter `any` non justifié.

## 16.4 TypeScript

- [x] typecheck principal vert ;
- [x] typecheck scientifique sans DOM vert ;
- [ ] auditer avant release la présence effective des options finales :

```text
strict: true
noUncheckedIndexedAccess: true
exactOptionalPropertyTypes: true
noImplicitOverride: true
noFallthroughCasesInSwitch: true
```

Les types doivent progressivement différencier :

- [ ] nombres quantiques validés ;
- [ ] unités/grandeurs ;
- [ ] état complexe ;
- [ ] orbitale réelle ;
- [ ] seed ;
- [ ] résultat d'échantillonnage ;
- [ ] résultat Worker.

---

# 17. Performance

## 17.1 Objectifs

- [ ] UI réactive pendant génération ;
- [ ] annulation logique immédiate d'un job obsolète ;
- [ ] 60 FPS lorsque la charge GPU le permet ;
- [ ] dégradation propre sur matériel faible ;
- [ ] aucune fuite Three.js ;
- [ ] mesures reproductibles avant toute optimisation.

## 17.2 Progressive rendering

Possibilité future :

```text
5 000
15 000
30 000
60 000 points
```

- [ ] ne pas mélanger des lots de deux états différents ;
- [ ] préserver la seed et le jobId ;
- [ ] mesurer le bénéfice réel avant adoption.

## 17.3 Qualité adaptative

Peut adapter :

- pixel ratio ;
- nombre de points ;
- taille de points ;
- fréquence de mise à jour.

Ne doit jamais adapter :

- la fonction d'onde ;
- la distribution de probabilité ;
- les constantes physiques ;
- le modèle scientifique.

---

# 18. Import / export

## 18.1 Import JSON

Schéma cible :

```json
{
  "schemaVersion": 1,
  "units": "bohr",
  "state": {
    "basis": "complex",
    "n": 2,
    "l": 1,
    "m": 0
  },
  "seed": 123456,
  "points": []
}
```

- [ ] unité obligatoire ;
- [ ] version de schéma obligatoire ;
- [ ] état quantique validé ;
- [ ] base explicite ;
- [ ] seed optionnelle/explicite selon le type d'import ;
- [ ] refuser les unités ambiguës ;
- [ ] supprimer la conversion legacy `5.29` lors de la migration.

## 18.2 Export

Priorités :

- [ ] capture PNG ;
- [ ] état JSON reproductible ;
- [ ] paramètres + seed ;
- [ ] CSV/JSON scientifique si utile ;
- [ ] n'annoncer OBJ/glTF que si réellement implémentés et pertinents.

---

# 19. Documentation finale

## 19.1 README.md

Doit contenir uniquement ce qui existe réellement :

- [ ] capture actuelle ;
- [ ] URL réelle ;
- [ ] installation réelle ;
- [ ] scripts npm réels ;
- [ ] stack réellement utilisée ;
- [ ] limites scientifiques ;
- [ ] crédit Kavan / `kavan010` ;
- [ ] statut de licence exact ;
- [ ] lien vers `SCIENCE.md` ;
- [ ] lien vers `VALIDATION.md` ;
- [ ] aucune fonction inexistante.

## 19.2 CREDITS.md

- [x] Kavan / `kavan010` ;
- [x] dépôt original ;
- [x] rôle historique/pédagogique ;
- [ ] audit final avant release.

## 19.3 SCIENCE.md

- [x] contrat initial ;
- [x] CODATA 2022 ;
- [x] masse réduite ;
- [x] rapport direct `m_e/m_p` ;
- [x] énergie ;
- [x] conventions Laguerre/Legendre ;
- [ ] fonctions d'onde ;
- [ ] orbitales réelles ;
- [ ] sampling ;
- [ ] renderer scientifique ;
- [ ] photon–hydrogène 2D.

## 19.4 VALIDATION.md

- [ ] tolérances ;
- [ ] valeurs analytiques ;
- [ ] normalisations ;
- [ ] orthogonalités ;
- [ ] tests statistiques ;
- [ ] captures/régressions ;
- [ ] navigateurs ;
- [ ] performance ;
- [ ] limites connues.

---

# 20. CI

Pipeline cible sur `main` :

```text
npm ci
npm run typecheck
npm run typecheck:science
npm run lint
npm run test
npm run build
npm run test:e2e
```

- [x] CI présente ;
- [x] baseline Node 24 LTS ;
- [x] typecheck scientifique indépendant ajouté ;
- [ ] étendre aux futurs tests sampling/régression ;
- [ ] aucun déploiement si un test scientifique critique échoue.

---

# 21. Séquence de refonte mise à jour

# Phase 0 — Gel et sécurisation ✅

## Actions

- [x] créer le tag `legacy-web-2026-08-23` ;
- [x] normaliser les fins de ligne avec `.gitattributes` ;
- [x] confirmer `main` comme seule branche permanente ;
- [x] conserver l'historique ;
- [x] créer `PLAN.md` ;
- [x] créer `CREDITS.md` minimal ;
- [ ] retirer toutes les déclarations de licence non vérifiées lors de la révision finale du README.

## Critère de sortie

- [x] le point de départ est récupérable ;
- [x] la provenance est documentée.

---

# Phase 1 — Toolchain moderne et séparation des responsabilités ✅

## Actions

- [x] initialiser npm ;
- [x] installer Vite / TypeScript / Three.js ;
- [x] convertir Three.js CDN vers npm/ESM ;
- [x] migrer le legacy vers TypeScript sans correction scientifique majeure ;
- [x] ajouter ESLint ;
- [x] ajouter Prettier ;
- [x] ajouter Vitest ;
- [x] ajouter Playwright ;
- [x] ajouter CI ;
- [x] séparer app / data / rendering / science / state / ui ;
- [x] isoler Three.js dans le renderer ;
- [x] réduire `main.ts` au rôle de composition root ;
- [x] ajouter le typecheck scientifique sans DOM ;
- [x] réparer le bug de visibilité du panneau 2D ;
- [x] ajouter un test E2E de non-régression du panneau 2D.

## Critère de sortie

- [x] l'application tourne sous Vite + TypeScript ;
- [x] le comportement legacy est conservé hors correctif 2D explicitement validé ;
- [x] les responsabilités sont séparées.

---

# Phase 2 — Fondations scientifiques pures ✅

## Actions

### Lot 2

- [x] centraliser CODATA 2022 ;
- [x] créer les conversions d'unités ;
- [x] valider les nombres quantiques ;
- [x] définir la masse réduite ;
- [x] définir `a_μ` ;
- [x] implémenter l'énergie de Schrödinger non relativiste avec masse réduite ;
- [x] implémenter factorielle ;
- [x] implémenter Laguerre généralisé ;
- [x] implémenter Legendre/Ferrers associé ;
- [x] fixer la convention de Condon–Shortley ;
- [x] créer les tests scientifiques ;
- [x] créer `docs/SCIENCE.md` ;
- [x] maintenir l'isolation DOM/Three.js ;
- [x] laisser `legacyScience.ts` intact ;
- [x] ne pas connecter le nouveau moteur au runtime.

### Micro-Lot 2.1 CODATA

- [x] ajouter la constante directe `m_e/m_p` CODATA 2022 ;
- [x] utiliser `μ/m_e = 1/(1+m_e/m_p)` pour `¹H` ;
- [x] conserver les masses SI générales ;
- [x] conserver la fonction générique de masse réduite ;
- [x] mettre à jour les tests métrologiques ;
- [x] mettre à jour `SCIENCE.md` ;
- [x] maintenir le runtime visible inchangé.

## Critère de sortie

- [x] constantes, unités, quantum numbers, énergie et fonctions spéciales de base validés ;
- [x] noyau scientifique compilable sans DOM ;
- [x] aucune dépendance Three.js ;
- [x] tests scientifiques verts ;
- [x] fondations figées avant fonctions d'onde.

---

# Phase 3 — Fonctions d'onde et bases orbitales ✅ VALIDÉE

## Actions

### Mathématique complexe

- [x] créer un type/utilitaires complexes minimaux et purs ;
- [x] tester addition, multiplication, conjugaison, module² et phase ;
- [x] éviter une dépendance externe si elle n'apporte pas de bénéfice réel.

### Fonction radiale

- [x] implémenter `Rₙₗ(r)` ;
- [x] utiliser `a_μ` de manière explicite pour le modèle `¹H` ;
- [x] utiliser `ρ = 2r/(n a_μ)` ;
- [x] documenter la normalisation ;
- [x] valider `1s`, `2s`, `2p` contre formes analytiques ;
- [x] tester normalisation radiale ;
- [x] tester nœuds radiaux.

### Harmoniques sphériques complexes

- [x] implémenter `Yₗᵐ(θ,φ)` ;
- [x] utiliser la convention DLMF déjà choisie ;
- [x] ne pas doubler la phase de Condon–Shortley ;
- [x] gérer `m<0` ;
- [x] tester normalisation angulaire ;
- [x] tester orthogonalité ;
- [x] tester relation `+m/-m` ;
- [x] tester parité.

### Fonction d'onde complète

- [x] implémenter `ψₙₗₘ = RₙₗYₗᵐ` ;
- [x] exposer amplitude complexe ;
- [x] exposer phase ;
- [x] exposer `|ψ|²` ;
- [x] tester normalisation 3D ;
- [x] tester plusieurs états connus.

### Orbitales réelles

- [x] définir une représentation typée des orbitales réelles ;
- [x] implémenter combinaisons `±m` normalisées ;
- [x] `p_x` ;
- [x] `p_y` ;
- [x] `p_z` ;
- [x] famille `d` standard ;
- [x] documenter le périmètre volontairement limité aux familles `p` et `d` ;
- [x] documenter conventions de signe ;
- [x] tester orientation et orthogonalité ;
- [x] vérifier que distribution radiale ne dépend pas de l'orientation réelle.

### Observables et nœuds de base

- [x] centraliser calcul des nombres de nœuds ;
- [x] tester `n-l-1`, `l`, `n-1` ;
- [x] ajouter les premières observables analytiques simples nécessaires à la validation.

### Isolation

- [x] aucun DOM ;
- [x] aucun Three.js ;
- [x] aucun sampler ;
- [x] aucun Worker ;
- [x] aucun branchement runtime pendant la construction initiale.

## Critère de sortie

- [x] le moteur calcule correctement `Rₙₗ` ;
- [x] le moteur calcule correctement `Yₗᵐ` ;
- [x] le moteur calcule `ψₙₗₘ`, sa phase et `|ψ|²` ;
- [x] les orbitales réelles sont des combinaisons scientifiques correctes ;
- [x] normalisation et orthogonalité passent ;
- [x] le tout fonctionne sans DOM, Three.js ni sampling.

---

# Phase 4 — Sampler scientifique

## Actions

### Reproductibilité

- [x] PRNG seedable ;
- [x] contrat de seed documenté ;
- [x] vecteurs de test déterministes ;
- [x] aucun `Math.random()` dans le nouveau sampler.

### Radial

- [x] CDF radiale depuis `r²|Rₙₗ|²` ;
- [x] intégration numérique documentée ;
- [x] domaine radial convergé ;
- [x] contrôle de masse tronquée ;
- [x] inversion CDF interpolée ;
- [x] cache scientifique versionné ;
- [x] tests de convergence.

### Angulaire

- [x] sampling `|Yₗᵐ|²` pour base complexe ;
- [x] sampling des orbitales réelles depuis leur densité réelle ;
- [x] mesure `sinθ` correctement prise en compte ;
- [x] tests de symétrie et orientation.

### Pipeline unique

- [x] supprimer samplers ad hoc ;
- [x] supprimer `NAMED_SAMPLERS` ;
- [x] supprimer `SAMPLER_MAP` ;
- [x] supprimer mapping `m=±1 → p_x/p_y` ;
- [x] unités uniques ;
- [x] même orchestrateur de sampling pour toutes les représentations.

### Validation statistique

- [x] tests ECDF/CDF ;
- [x] tests angulaires ;
- [x] moyennes/quantiles ;
- [x] convergence avec `N` ;
- [x] seuils statistiques justifiés ;
- [x] aucun test flaky.

## Critère de sortie

- [x] toutes les orbitales utilisent le même pipeline scientifique ;
- [x] une seed fixe reproduit la même génération à version identique du moteur ;
- [x] les distributions générées passent les tests statistiques.

---

# Phase 5 — Worker et concurrence

## Actions

- [x] déplacer la génération dans un Web Worker ;
- [x] protocole `jobId` ;
- [x] seed dans le protocole ;
- [x] état scientifique immuable par job ;
- [x] buffers transférables ;
- [x] progression associée au bon job ;
- [x] annulation logique des jobs obsolètes ;
- [x] résultat obsolète ignoré ;
- [x] tests de changements rapides d'état ;
- [x] supprimer les buffers globaux de génération legacy.

## Critère de sortie

- [x] impossible qu'un résultat ancien écrase un nouvel état ;
- [x] l'UI reste réactive pendant la génération.

---

# Phase 6 — Renderer 3D scientifique

## Actions

- [x] brancher le renderer sur le nouveau sampler ;
- [x] nuage probabiliste explicite ;
- [x] échelle physique ;
- [x] cadrage automatique ;
- [x] axes ;
- [x] légendes ;
- [x] mode densité ;
- [x] mode phase ;
- [x] isosurfaces ;
- [x] surfaces nodales ;
- [x] gestion propre des ressources ;
- [x] noyau signalé non à l'échelle ;
- [x] supprimer pseudo-trajectoires des points ;
- [x] rotation caméra/scène clairement distinguée de la dynamique physique ;
- [ ] courant de probabilité éventuel représenté comme champ.

## Critère de sortie

- [x] le rendu est physiquement interprétable ;
- [x] aucun effet visuel n'implique une orbite classique fictive ;
- [x] l'échelle et les légendes sont explicites.

---

# Phase 7 — UI/UX scientifique complète

## Actions

- [x] séparer base complexe et base réelle ;
- [x] reconstruire presets ;
- [x] corriger labels `m` ;
- [x] afficher énergie du nouveau moteur ;
- [x] afficher nœuds ;
- [x] afficher unités ;
- [x] expliquer nuage probabiliste ;
- [x] expliquer phase ;
- [x] expliquer modèle physique actif ;
- [x] moderniser panneau sans surcharge décorative ;
- [x] finaliser sombre/clair ;
- [x] accessibilité clavier ;
- [x] responsive ;
- [x] `prefers-reduced-motion`.

## Critère de sortie

- [x] aucun label UI ne contredit le moteur scientifique ;
- [x] aucune représentation n'est ambiguë sur sa signification physique.

---

# Phase 8 — Interactions photon–hydrogène 2D

## Actions

- [x] panneau visible/cachable ;
- [x] raccourci `Q` ;
- [x] test E2E de visibilité ;
- [ ] supprimer le modèle de Bohr comme représentation physique principale ;
- [ ] renommer le module ;
- [ ] implémenter transitions énergétiques ;
- [ ] distinguer excitation et ionisation ;
- [ ] calculer énergie des photons ;
- [ ] photoionisation ;
- [ ] documenter hypothèses ;
- [ ] ajouter tests.

## Critère de sortie

- [ ] le module 2D possède une signification physique définie et vérifiable.

---

# Phase 9 — Validation globale

## Actions

- [x] infrastructure Vitest ;
- [x] infrastructure Playwright ;
- [x] tests des fondations scientifiques ;
- [x] campagne fonctions d'onde ;
- [x] campagne sampler ;
- [x] campagne Worker/concurrence ;
- [ ] régressions visuelles déterministes ;
- [x] Chrome ;
- [ ] Firefox ;
- [ ] WebKit ;
- [x] responsive ;
- [ ] profilage CPU ;
- [ ] profilage GPU ;
- [ ] chasse aux fuites mémoire ;
- [ ] remplir `VALIDATION.md`.

## Critère de sortie

- [ ] `VALIDATION.md` contient des résultats reproductibles ;
- [ ] toutes les validations critiques passent.

---

# Phase 10 — Documentation, attribution et autonomie GitHub

## Actions

- [ ] réécrire README ;
- [x] `CREDITS.md` initial ;
- [ ] finaliser `CREDITS.md` ;
- [x] `SCIENCE.md` initial ;
- [ ] finaliser `SCIENCE.md` ;
- [ ] créer `REFERENCES.md` ;
- [ ] finaliser `VALIDATION.md` ;
- [ ] clarifier licence ;
- [ ] vérifier droits des documents redistribués ;
- [ ] traiter le remote `upstream` ;
- [ ] vérifier qu'`origin` ne possède que les branches voulues ;
- [ ] décider autonomie/fork network ;
- [ ] configurer déploiement réel ;
- [ ] ajouter badges réels uniquement.

## Critère de sortie

- [ ] le dépôt public décrit exactement ce qu'il contient, ce qu'il calcule et d'où il vient.

---

# Phase 11 — Release 5.0

## Versionnement proposé

```text
v5.0.0-alpha.1  moteur scientifique + sampler opérationnels
v5.0.0-beta.1   rendu + UI stabilisés
v5.0.0-rc.1     validation scientifique complète
v5.0.0          release publique
```

## Aucun `v5.0.0` si

- [ ] toutes les normalisations scientifiques critiques passent ;
- [ ] toutes les unités sont non ambiguës ;
- [ ] `m` et orbitales réelles sont séparés ;
- [ ] les orbitales réelles sont correctes ;
- [ ] le sampler est unique, déterministe et validé ;
- [ ] la concurrence de génération est résolue ;
- [ ] aucune animation ne suggère une orbite électronique classique ;
- [ ] le module 2D est physiquement défini ;
- [ ] le README ne contient aucune fonction inexistante ;
- [ ] le statut de licence est présenté correctement ;
- [ ] CI critique verte.

Les cases ci-dessus doivent toutes être cochées avant la release finale.

---

# 22. Definition of Done `v5.0.0`

## Git / provenance

- [x] `main` est la seule branche permanente de travail.
- [x] Le point legacy est protégé par tag.
- [x] L'origine Kavan / `kavan010` est clairement créditée.
- [ ] Le statut de licence est exact dans toute la documentation publique.

## Toolchain / qualité

- [x] TypeScript compile sans erreur sur l'état actuel.
- [x] Le noyau scientifique compile sans DOM.
- [ ] Toutes les options TypeScript strictes finales ont été auditées.
- [x] ESLint passe sans erreur bloquante.
- [x] Prettier passe sur l'état actuel.
- [x] Tous les tests Vitest actuels passent.
- [x] Tous les tests scientifiques actuellement implémentés passent.
- [x] Tous les tests Playwright actuels passent.
- [x] CI est configurée sur Node 24 LTS.

## Science

- [x] CODATA 2022 centralisé.
- [x] `m_e/m_p` direct CODATA utilisé pour `¹H`.
- [x] `a₀ ↔ pm` centralisé et testé dans le nouveau noyau.
- [x] masse réduite documentée et testée.
- [x] énergie de Schrödinger avec masse réduite implémentée.
- [x] Laguerre et Legendre/Ferrers de base validés.
- [x] `Rₙₗ` validé.
- [x] `Yₗᵐ` validé.
- [x] `ψₙₗₘ` et `|ψ|²` validés.
- [x] orbitales réelles correctement construites dans le périmètre `p`/`d` documenté.
- [x] normalisation complète validée.
- [x] orthogonalité validée.
- [x] nœuds validés scientifiquement.

## Sampling

- [x] toutes les longueurs du runtime orbital utilisent une convention unique.
- [x] aucun sampler orbital spécial incohérent ne subsiste.
- [x] aucun `SAMPLER_MAP` legacy ne subsiste.
- [x] `m=±1` n'est plus étiqueté `p_x/p_y`.
- [x] PRNG seedable.
- [x] CDF radiale validée.
- [x] sampling angulaire validé.
- [x] tests statistiques verts.
- [x] même pipeline pour toutes les représentations.

## Runtime / concurrence

- [x] le nuage visible provient du nouveau `|ψ|²`.
- [x] les jobs concurrents sont sûrs et annulables.
- [x] aucun résultat obsolète ne peut modifier la scène.
- [x] buffers transférables utilisés correctement.

## Rendu

- [x] aucune animation ne représente implicitement une trajectoire classique inexistante.
- [x] le rendu affiche une échelle physique.
- [x] le noyau schématique est signalé comme non à l'échelle.
- [x] densité et phase sont visuellement distinctes.
- [x] légendes scientifiques présentes.
- [ ] régressions visuelles déterministes disponibles.

## UI / 2D

- [x] le panneau 2D fonctionne réellement.
- [ ] le module 2D possède un modèle physique explicite.
- [x] base complexe et base réelle séparées dans l'UI.
- [x] aucune notation UI ne contredit le moteur.
- [ ] accessibilité finale validée.

## Documentation

- [x] README ne décrit que des fonctions présentes dans le périmètre livré.
- [x] `SCIENCE.md` décrit le socle et le pipeline 3D intégrés.
- [ ] `REFERENCES.md` est complet.
- [ ] `VALIDATION.md` documente tolérances et références.
- [ ] licence/provenance/documentation sont juridiquement cohérentes.

---

# 23. Ordre d'exécution actuel

État :

1. [x] gel et tag legacy ;
2. [x] toolchain moderne ;
3. [x] migration TypeScript / Three.js npm ;
4. [x] séparation des responsabilités ;
5. [x] réparation du panneau 2D ;
6. [x] fondations scientifiques ;
7. [x] Micro-Lot 2.1 CODATA ;
8. [x] Phase 3 : fonctions d'onde et bases orbitales ;
9. [x] **Phase 4 : sampler scientifique** ;
10. [x] Phase 5 : Worker/concurrence ;
11. [x] Phase 6 : renderer scientifique ;
12. [x] Phase 7 : UI/UX scientifique ;
13. [ ] Phase 8 : photon–hydrogène 2D ;
14. [ ] Phase 9 : validation globale ;
15. [ ] Phase 10 : documentation/autonomie GitHub ;
16. [ ] Phase 11 : release `v5.0.0`.

Règle : **ne pas commencer une phase dépendante tant que le critère de sortie scientifique de la phase précédente n'est pas atteint.**

---

# 24. Sources techniques du socle

Références de maintenance :

- Node.js releases: https://nodejs.org/en/about/previous-releases
- TypeScript: https://www.typescriptlang.org/
- Vite: https://vite.dev/
- Three.js: https://threejs.org/
- Three.js migration guide: https://github.com/mrdoob/three.js/wiki/Migration-Guide
- Vitest: https://vitest.dev/
- Playwright: https://playwright.dev/
- ESLint: https://eslint.org/
- Prettier: https://prettier.io/

Les numéros de version exacts doivent être vérifiés lors de chaque future montée de dépendance. Le plan ne doit pas forcer une mise à jour automatique simplement parce qu'une nouvelle version existe.

---

# 25. Décision d'architecture actuelle

La refonte `v5` repose désormais sur :

```text
TypeScript
Vite
Three.js / WebGL2 baseline
DOM natif
src/science comme noyau scientifique pur
unités atomiques et conversions centralisées
modèle ¹H non relativiste avec masse réduite
data CODATA 2022 / NIST
fonctions spéciales selon NIST DLMF
Vitest + Playwright
sampler seedable validé
Web Worker validé après le sampler
main uniquement + tags de jalons
```

Décision importante :

```text
Phase 3 = fonctions d'onde et bases orbitales
Phase 4 = sampler scientifique
```

Le sampler ne doit pas être construit avant que `Rₙₗ`, `Yₗᵐ`, `ψₙₗₘ`, `|ψ|²` et les orbitales réelles soient validés.

---

# 26. Condition de changement du plan

Le présent document fait foi tant qu'une décision explicite ne le remplace pas.

Toute modification importante doit :

- [ ] être motivée par une contrainte scientifique, technique, juridique ou de validation ;
- [ ] être inscrite dans `PLAN.md` ;
- [ ] mettre à jour les critères de sortie concernés ;
- [ ] ne pas réécrire rétroactivement l'historique des décisions ;
- [ ] conserver la distinction entre dette legacy et nouveau moteur validé.

---

# 27. Prochaine action officielle

## Phases 4–7 — intégration scientifique et UI

Les Phases 4, 5, 6 et 7 sont implémentées et validées par les tests scientifiques,
le build et les scénarios E2E Chromium. Elles fournissent le pipeline suivant :

```text
état scientifique explicite
PRNG seedable et versionné
CDF radiale issue de r²|Rₙₗ|²
sampling angulaire avec la mesure sinθ dθ dφ
pipeline commun aux bases complexe et réelle
métadonnées de reproductibilité
tests statistiques déterministes
```

Le prochain lot officiel est la **Phase 9 — validation globale**, puis la
documentation finale. La démonstration photon–hydrogène 2D reste explicitement
isolée en Phase 8 et n'est pas présentée comme le modèle orbital 3D.
