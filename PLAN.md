# Atoms — Plan de refonte scientifique, technique et graphique

**Document de référence — 23 août 2026**  
**Projet :** `Giscolab/Atoms`  
**Origine :** réimplémentation web issue du travail pédagogique original de Kavan (`kavan010/Atoms`)  
**Cible de refonte :** `v5.0.0`

---

## 1. But de la refonte

Atoms doit devenir un visualiseur scientifique de l'atome d'hydrogène qui soit à la fois :

- scientifiquement défendable ;
- graphiquement lisible et moderne ;
- déterministe et testable ;
- performant sur navigateur moderne ;
- maintenable ;
- transparent sur ses hypothèses physiques ;
- explicite sur ce qui est une représentation pédagogique et ce qui est une grandeur physique ;
- respectueux du travail de l'auteur original.

La priorité absolue est la **justesse scientifique**. Aucun effet visuel ne doit suggérer un phénomène physique faux.

---

## 2. Principes non négociables

1. **La physique est la source de vérité.**
   - Une visualisation ne doit jamais être conservée uniquement parce qu'elle est jolie.
   - Toute grandeur affichée doit avoir une unité, une définition et une provenance.

2. **Séparer calcul, échantillonnage, rendu et interface.**
   - Le moteur quantique ne connaît pas Three.js.
   - Three.js ne connaît pas les formules physiques.
   - L'UI ne contient aucune formule de physique.

3. **Aucune trajectoire classique fictive de l'électron.**
   - Un nuage de points représente des échantillons de `|ψ|²`, pas des électrons individuels en mouvement.
   - Un état stationnaire ne doit pas être représenté comme une planète orbitant autour du noyau.

4. **Deux notions différentes doivent rester différentes.**
   - Les états propres complexes `|n,l,m⟩`.
   - Les orbitales réelles usuelles `p_x`, `p_y`, `d_xy`, etc., qui sont des combinaisons linéaires d'états `m`.

5. **Unités cohérentes de bout en bout.**
   - Calcul interne en unités atomiques lorsque cela simplifie et sécurise les équations.
   - Conversion vers `a₀`, pm, nm, Hartree et eV uniquement aux frontières d'affichage/import/export.

6. **Toute assertion scientifique importante doit être testée.**

7. **Le dépôt ne conserve qu'une branche de travail : `main`.**
   - Les jalons sont conservés par tags Git, pas par branches permanentes.

8. **L'auteur original doit rester crédité clairement.**
   - Le projet devient autonome techniquement sans effacer sa provenance historique.

---

## 3. État de départ constaté

### 3.1 Architecture actuelle

Le dépôt actuel est une application web sans build :

```text
Atoms/
├── index.html
├── styles.css
├── script.js
├── screenshot.png
└── README.md
```

Le moteur scientifique, le rendu Three.js, l'UI, la simulation 2D et les entrées utilisateur sont regroupés principalement dans `script.js`.

### 3.2 Git

État observé dans l'archive du 23 août 2026 :

- branche locale : `main` uniquement ;
- `origin/main` : `Giscolab/Atoms` ;
- remote `upstream` : `kavan010/Atoms` ;
- `main` et `upstream/main` n'ont actuellement aucun ancêtre commun ;
- l'origine du projet Kavan est préservée explicitement par la documentation et les crédits, indépendamment de la topologie Git ;
- les quatre fichiers texte signalés modifiés dans l'archive ne présentent pas de différence sémantique avec `origin/main` : il s'agit de fins de ligne.

### 3.3 Défauts déjà identifiés

- mélange d'unités entre rayons de Bohr et picomètres ;
- constante d'import JSON `5.29` incorrecte pour une conversion `a₀ → pm` ;
- samplers spéciaux `2p_x`, `2p_y`, `2p_z`, `3p_z` utilisant une échelle différente du moteur générique ;
- confusion entre `m=±1` et les orbitales réelles `p_x` / `p_y` ;
- risque de concurrence entre plusieurs appels de `generateCloud()` à cause des tableaux globaux partagés ;
- panneau 2D masqué par l'attribut HTML `hidden` alors que le code ne le retire pas ;
- animation du nuage susceptible d'être interprétée à tort comme une trajectoire électronique ;
- énergie affichée avec l'approximation `-13.6/n²` sans préciser le modèle ;
- simulation dite « photoélectrique » ne modélisant pas actuellement un véritable effet photoélectrique ni une photoionisation ;
- README annonçant des fonctions, fichiers, commandes, versions et une licence qui ne correspondent pas au dépôt réel.

---

## 4. Point juridique et attribution

### 4.1 Crédit

Le futur README devra conserver une section visible du type :

> Atoms est une réimplémentation web indépendante inspirée du projet pédagogique original **Hydrogen Quantum Orbital Visualizer** de **Kavan** (`kavan010/Atoms`). Le dépôt original a fourni le point de départ historique et pédagogique du projet.

Le lien vers le dépôt original doit être conservé :

- https://github.com/kavan010/Atoms

### 4.2 Licence : point à corriger avant publication finale

Le dépôt original observé ne contient pas de fichier `LICENSE` à sa racine.

Par conséquent :

- le README actuel ne doit pas déclarer arbitrairement « licence MIT » ;
- « code visible publiquement sur GitHub » ne signifie pas automatiquement « code sous licence open source » ;
- le crédit à l'auteur ne remplace pas une licence ;
- la licence finale du projet doit être décidée uniquement après clarification de la provenance et des droits applicables aux éléments conservés.

**Action de refonte :** retirer l'affirmation MIT du README tant qu'elle n'est pas juridiquement fondée.

---

## 5. Stratégie Git : une seule branche `main`

### 5.1 Avant toute opération destructrice

Créer un tag immuable du point de départ :

```text
legacy-web-2026-08-23
```

Ce tag constitue le retour arrière de référence.

### 5.2 Politique cible

```text
main        ← seule branche permanente
v5.0.0-*    ← tags de versions
```

Pas de branches historiques conservées pour la refonte.

### 5.3 Remote `upstream`

Après création de `CREDITS.md` et correction du README :

- conserver le lien historique vers `kavan010/Atoms` dans la documentation ;
- supprimer le remote Git `upstream` du clone de travail ;
- ne jamais supprimer ni modifier le dépôt de l'auteur original.

Les historiques `main` et `upstream/main` n'ont actuellement aucun ancêtre commun. Ils ne doivent pas être greffés artificiellement.

L'origine du projet Kavan est préservée explicitement par la documentation et les crédits, indépendamment de la topologie Git.

### 5.4 Détacher le fork GitHub

Si le dépôt doit devenir officiellement autonome sur GitHub, utiliser la fonction **Leave fork network** lorsque les conditions GitHub sont satisfaites.

Avant de le faire :

- vérifier les issues ;
- vérifier les pull requests ;
- vérifier les stars/watchers ;
- vérifier l'absence de forks enfants ;
- effectuer une sauvegarde locale complète.

GitHub indique que la sortie du réseau de forks conserve les métadonnées des commits mais peut faire perdre les métadonnées GitHub liées au fork (issues, PR, wiki, stars, watchers, commentaires, etc.). L'opération est permanente.

Référence : https://docs.github.com/en/pull-requests/how-tos/work-with-forks/detaching-a-fork

### 5.5 Ne pas réécrire l'histoire

La refonte ne doit pas faire de `git filter-repo`, squash global ou réécriture complète de `main` uniquement pour masquer l'origine du projet.

L'historique existant est utile :

- techniquement ;
- historiquement ;
- pour l'attribution.

---

## 6. Stack cible moderne

Versions stables vérifiées le 23 août 2026.

| Composant | Version cible | Rôle |
|---|---:|---|
| Node.js | **24.19.0 LTS** | environnement CI et référence de développement |
| Node.js | **26.7.0 Current** | autorisé localement, non utilisé comme baseline CI |
| TypeScript | **7.0.2** | langage principal |
| Vite | **8.2.2** | serveur de développement et build |
| Three.js | **0.185.1 / r185** | rendu 3D |
| Vitest | **4.1.10** | tests unitaires et scientifiques |
| Playwright Test | **1.62.1** | tests navigateur end-to-end |
| ESLint | **10.7.0** | analyse statique |
| Prettier | **3.9.0** | formatage déterministe |

### 6.1 Pourquoi Node 24 LTS plutôt que Node 26 comme baseline

Node 26 est la version Current au 23 août 2026. Node recommande les versions LTS pour les applications de production. La CI doit donc être reproductible sur Node 24 LTS, tout en autorisant Node 26 pour le développement local.

Référence : https://nodejs.org/en/about/previous-releases

### 6.2 Pourquoi pas React

La refonte ne doit pas ajouter un framework UI sans nécessité.

Atoms a besoin :

- d'un moteur scientifique robuste ;
- d'un renderer 3D ;
- d'un panneau de contrôle ;
- d'un état applicatif limité.

TypeScript + DOM natif + composants simples suffisent et réduisent :

- la surface de dépendances ;
- les abstractions ;
- le coût de maintenance ;
- les interactions indirectes entre UI et moteur scientifique.

React pourra être reconsidéré uniquement si la complexité réelle de l'interface le justifie plus tard.

### 6.3 WebGL / WebGPU

Baseline `v5.0.0` :

- WebGL2 stable via Three.js ;
- architecture de renderer suffisamment isolée pour permettre un backend WebGPU ultérieur ;
- WebGPU ne doit pas devenir une dépendance fonctionnelle obligatoire avant validation multi-navigateurs.

---

## 7. Architecture cible

```text
Atoms/
├── public/
│   └── assets/
├── src/
│   ├── main.ts
│   │
│   ├── app/
│   │   ├── AtomsApp.ts
│   │   ├── state.ts
│   │   └── commands.ts
│   │
│   ├── physics/
│   │   ├── constants.ts
│   │   ├── units.ts
│   │   ├── quantumNumbers.ts
│   │   ├── energy.ts
│   │   ├── hydrogen/
│   │   │   ├── radialWavefunction.ts
│   │   │   ├── sphericalHarmonics.ts
│   │   │   ├── realOrbitals.ts
│   │   │   ├── wavefunction.ts
│   │   │   ├── probabilityDensity.ts
│   │   │   └── observables.ts
│   │   └── special/
│   │       ├── gamma.ts
│   │       ├── factorial.ts
│   │       ├── laguerre.ts
│   │       └── legendre.ts
│   │
│   ├── sampling/
│   │   ├── rng.ts
│   │   ├── radialSampler.ts
│   │   ├── angularSampler.ts
│   │   ├── orbitalSampler.ts
│   │   └── cdf.ts
│   │
│   ├── workers/
│   │   ├── orbital.worker.ts
│   │   └── protocol.ts
│   │
│   ├── render/
│   │   ├── SceneRenderer.ts
│   │   ├── camera.ts
│   │   ├── pointCloud.ts
│   │   ├── isosurface.ts
│   │   ├── axes.ts
│   │   ├── phaseColors.ts
│   │   └── dispose.ts
│   │
│   ├── features/
│   │   └── photonHydrogen2d/
│   │       ├── model.ts
│   │       ├── transitions.ts
│   │       ├── photoionization.ts
│   │       └── renderer2d.ts
│   │
│   ├── ui/
│   │   ├── controls.ts
│   │   ├── presets.ts
│   │   ├── scientificInfo.ts
│   │   ├── loading.ts
│   │   └── accessibility.ts
│   │
│   └── styles/
│       ├── tokens.css
│       ├── base.css
│       ├── layout.css
│       ├── controls.css
│       └── viewport.css
│
├── tests/
│   ├── unit/
│   ├── scientific/
│   ├── sampling/
│   ├── regression/
│   └── e2e/
│
├── docs/
│   ├── SCIENCE.md
│   ├── VALIDATION.md
│   └── REFERENCES.md
│
├── CREDITS.md
├── PLAN.md
├── README.md
├── index.html
├── package.json
├── package-lock.json
├── tsconfig.json
├── vite.config.ts
├── vitest.config.ts
├── playwright.config.ts
├── eslint.config.js
├── .prettierrc
├── .editorconfig
└── .gitattributes
```

---

## 8. Contrat scientifique du moteur quantique

### 8.1 Domaine initial

`v5.0.0` modélise en priorité :

- l'atome d'hydrogène neutre `¹H` ;
- un électron ;
- un potentiel coulombien central ;
- les solutions analytiques de l'équation de Schrödinger non relativiste.

Les limites du modèle doivent être affichées clairement.

Ne pas présenter le modèle comme incluant par défaut :

- structure fine ;
- Lamb shift ;
- spin-orbite ;
- QED ;
- champ magnétique externe ;
- effet Stark ;
- effet Zeeman ;
- interactions multi-électroniques.

Ces phénomènes pourront être ajoutés plus tard comme modules distincts.

### 8.2 Fonction d'onde

Le cœur doit représenter explicitement :

```text
ψₙₗₘ(r, θ, φ) = Rₙₗ(r) Yₗᵐ(θ, φ)
```

avec normalisation vérifiée numériquement.

### 8.3 Partie radiale

Implémenter les polynômes de Laguerre associés avec une convention documentée et testée.

Le sampler radial utilise la densité radiale :

```text
P(r) dr = r² |Rₙₗ(r)|² dr
```

Le facteur `r²` ne doit jamais être oublié ou appliqué deux fois.

### 8.4 Partie angulaire

Pour un état complexe `|n,l,m⟩` :

```text
Yₗᵐ(θ, φ) ∝ Pₗᵐ(cos θ) e^(imφ)
```

Conséquences à respecter dans l'interface :

- `m=+1` n'est pas `p_x` ;
- `m=-1` n'est pas `p_y` ;
- les densités de probabilité des états `+m` et `-m` sont identiques dans le modèle sans perturbation externe ;
- leur phase / courant de probabilité peut différer.

### 8.5 Orbitales réelles

Les orbitales réelles usuelles doivent être une catégorie UI séparée :

- `p_x`, `p_y`, `p_z` ;
- `d_xy`, `d_xz`, `d_yz` ;
- `d_x²-y²`, `d_z²` ;
- etc.

Elles sont construites comme combinaisons réelles des harmoniques sphériques complexes.

Le moteur doit donc exposer deux modes explicites :

```text
Base complexe : |n,l,m⟩
Base réelle    : pₓ, pᵧ, dxy, ...
```

### 8.6 Nœuds

Le moteur doit vérifier et afficher :

```text
nœuds radiaux   = n - l - 1
nœuds angulaires = l
nœuds totaux    = n - 1
```

La visualisation doit permettre de voir les surfaces nodales lorsque cela est pertinent.

### 8.7 Unités

Convention proposée :

#### Interne

- longueur : rayon de Bohr `a₀` ;
- énergie : Hartree `E_h` ;
- angles : radians.

#### Affichage

- `a₀` ;
- pm ;
- nm ;
- eV.

Les conversions sont centralisées dans `physics/units.ts`.

Aucune constante numérique physique ne doit être copiée à plusieurs endroits.

### 8.8 Constantes physiques

La source de référence pour les constantes est **CODATA 2022 / NIST**, qui est encore le jeu de valeurs recommandé le plus récent au 23 août 2026.

Références :

- https://physics.nist.gov/constants
- https://physics.nist.gov/cuu/Reference/versioncon.shtml

La prochaine révision CODATA est annoncée pour 2026 mais n'est pas encore le jeu recommandé publié au moment de ce plan.

### 8.9 Énergie

Ne plus afficher simplement :

```text
-13.6 / n²
```

sans contexte.

L'application doit distinguer :

1. **énergie du modèle de Schrödinger idéal** ;
2. éventuellement **valeur de référence spectroscopique NIST** lorsque la comparaison a un sens ;
3. les corrections physiques non incluses.

Le NIST fournit des valeurs critiques pour les niveaux de l'hydrogène, y compris des données plus fines que le modèle de Schrödinger simple.

Références :

- https://physics.nist.gov/PhysRefData/ASD/levels_form.html
- https://physics.nist.gov/PhysRefData/Handbook/Tables/hydrogentable5.htm

---

## 9. Échantillonnage : refonte complète

### 9.1 Supprimer les exceptions ad hoc

Supprimer à terme :

- `_rProb2p` ;
- `_rProb3p` ;
- `_sampleR2p` ;
- `_sampleR3p` ;
- `NAMED_SAMPLERS` ;
- `SAMPLER_MAP`.

Toutes les orbitales doivent passer par un moteur cohérent commun.

### 9.2 Sampler générique

Le sampler doit :

- accepter une définition d'état quantique ;
- produire des positions selon `|ψ|²` ;
- utiliser les mêmes unités pour toutes les orbitales ;
- être indépendant du renderer ;
- permettre une graine pseudo-aléatoire reproductible ;
- être testable statistiquement.

### 9.3 CDF

Les CDF ne doivent plus être de simples sommes de points sans contrôle d'erreur.

À faire :

- intégration numérique documentée ;
- interpolation dans l'inversion de CDF ;
- résolution adaptée aux nombres quantiques ;
- tests de convergence ;
- cache indexé par paramètres physiques ;
- invalidation explicite.

### 9.4 Reproductibilité

Introduire un PRNG seedable.

Exemple fonctionnel attendu :

```text
state = 3d, seed = 123456, N = 50 000
→ même tableau d'échantillons à version identique du moteur
```

Cela permet :

- tests de régression ;
- comparaison de renderer ;
- captures reproductibles ;
- débogage scientifique.

---

## 10. Concurrence et Web Worker

### 10.1 Problème actuel

`generateCloud()` écrit dans des tableaux globaux pendant plusieurs tâches asynchrones.

Deux générations rapprochées peuvent donc se mélanger.

### 10.2 Architecture cible

Chaque génération possède :

```text
jobId
state
seed
sampleCount
options
```

Le worker retourne :

```text
jobId
positions
colors / phase éventuelle
metadata
```

Le thread principal n'accepte un résultat que si `jobId` correspond encore à la requête active.

### 10.3 Buffers

Utiliser des `Float32Array` transférables entre worker et thread principal.

Ne pas garder des tableaux globaux réutilisés entre deux jobs.

### 10.4 Annulation

Une nouvelle génération annule logiquement l'ancienne.

Une génération obsolète :

- ne modifie pas la scène ;
- ne modifie pas le HUD ;
- ne modifie pas la progression ;
- peut s'arrêter dès que possible.

---

## 11. Rendu 3D scientifique

### 11.1 Modes de rendu

Prévoir au minimum :

1. **Nuage probabiliste**
   - échantillons de `|ψ|²` ;
   - nombre de points réglable ;
   - transparence contrôlée.

2. **Isosurface de densité**
   - surface d'une densité choisie ;
   - option de plusieurs niveaux ;
   - utile pour comparer aux représentations classiques d'orbitales.

3. **Surfaces nodales**
   - visualisation pédagogique distincte.

4. **Phase / signe**
   - palette divergente cohérente ;
   - légende visible ;
   - ne pas confondre couleur de phase et densité.

### 11.2 Échelle

Le renderer doit afficher une échelle réelle du modèle :

- valeur en `a₀` ;
- conversion pm ;
- axes optionnels ;
- repère d'échelle.

La caméra doit se cadrer automatiquement à partir de la distribution radiale de l'état, pas d'une distance arbitraire commune à toutes les orbitales.

### 11.3 Noyau

Si le noyau est représenté avec une taille visible, indiquer explicitement :

```text
Noyau représenté de manière schématique — taille non à l'échelle.
```

Une représentation à l'échelle réelle rendrait le proton pratiquement invisible à l'échelle atomique.

### 11.4 Animation

Supprimer toute animation pouvant laisser croire que les points sont des électrons classiques en orbite.

Modes autorisés :

- rotation de caméra ;
- rotation volontaire de la scène comme outil d'observation ;
- animation de phase clairement étiquetée ;
- champ de courant de probabilité correctement dérivé et présenté comme champ, pas comme trajectoires individuelles.

### 11.5 Couleurs

Séparer :

- couleur de densité ;
- couleur de phase ;
- sélection UI ;
- thème graphique.

Une palette scientifique doit :

- être lisible en mode sombre et clair ;
- rester interprétable en déficience de perception des couleurs ;
- disposer d'une légende.

---

## 12. Refonte graphique et UX

### 12.1 Direction visuelle

Objectif : interface d'instrument scientifique moderne, pas interface décorative.

Principes :

- surface 3D prioritaire ;
- panneau compact ;
- hiérarchie typographique forte ;
- peu d'effets gratuits ;
- informations scientifiques proches du contrôle concerné ;
- transitions courtes et fonctionnelles ;
- thème sombre et clair réellement équivalents.

### 12.2 Panneau orbital

Séparer clairement :

```text
État quantique
  n
  l
  m

Représentation
  base complexe / base réelle
  orbitales réelles disponibles

Rendu
  points
  isosurface
  nœuds
  phase

Échelle
  a₀ / pm / nm
```

### 12.3 Informations scientifiques

Afficher en temps réel :

- état choisi ;
- notation spectroscopique ;
- nombre de nœuds ;
- énergie du modèle ;
- rayon caractéristique ou statistiques utiles ;
- hypothèse physique active ;
- signification de la couleur ;
- signification du nuage.

### 12.4 Accessibilité

Minimum :

- navigation clavier complète ;
- focus visible ;
- labels explicites ;
- `aria-live` limité aux informations utiles ;
- contrastes vérifiés ;
- `prefers-reduced-motion` ;
- aucun contrôle uniquement identifiable par couleur.

---

## 13. Module 2D : remplacer « Photoélectrique » par une physique cohérente

### 13.1 Diagnostic

Le module actuel :

- dessine des atomes en modèle de Bohr ;
- utilise une différence d'énergie proche de la transition `1 → 2` ;
- dessine des ondes incidentes ;
- ne calcule pas une véritable émission photoélectrique ;
- ne modélise pas une photoionisation réelle.

### 13.2 Direction recommandée

Pour rester cohérent avec Atoms, renommer le module :

**Interactions photon–hydrogène 2D**

et y distinguer :

1. transitions liées-liées ;
2. absorption d'un photon ;
3. émission stimulée / spontanée uniquement si modélisée correctement ;
4. photoionisation lorsque l'énergie du photon dépasse le seuil de l'état considéré.

### 13.3 Photoionisation

Dans l'approximation simple :

```text
hν > |E_bound|
```

puis :

```text
K ≈ hν - |E_bound|
```

avec hypothèses et limites documentées.

### 13.4 Effet photoélectrique classique

Si l'on veut également enseigner l'effet photoélectrique d'un métal :

```text
K_max = hν - Φ
```

ce doit être un **module séparé**, car la fonction travail `Φ` d'un matériau n'est pas le seuil d'ionisation de l'atome d'hydrogène.

---

## 14. Tests scientifiques obligatoires

### 14.1 Nombres quantiques

Tester :

```text
n >= 1
0 <= l <= n - 1
-l <= m <= l
```

Aucun état invalide ne doit atteindre le moteur.

### 14.2 Normalisation radiale

Pour plusieurs couples `(n,l)` :

```text
∫₀∞ r² |Rₙₗ(r)|² dr ≈ 1
```

### 14.3 Normalisation angulaire

```text
∫ |Yₗᵐ(θ,φ)|² dΩ ≈ 1
```

### 14.4 Normalisation complète

```text
∫ |ψₙₗₘ|² dV ≈ 1
```

### 14.5 Orthogonalité

Tester numériquement des paires d'états distincts.

### 14.6 Symétries

Exemples :

- `1s` isotrope ;
- `m` et `-m` : même densité dans la base complexe ;
- `p_x`, `p_y`, `p_z` : mêmes distributions radiales, orientations différentes ;
- symétries de parité attendues.

### 14.7 Nœuds

Vérifier le nombre de nœuds radiaux et angulaires connu analytiquement.

### 14.8 Valeurs moyennes analytiques

Ajouter des tests sur des observables connues, par exemple :

- `⟨r⟩` ;
- `⟨1/r⟩` ;
- rayon le plus probable de `1s` ;
- énergie en fonction de `n` dans le modèle choisi.

### 14.9 Tests statistiques du sampler

Pour une seed fixe et un grand `N` :

- histogramme radial comparé à la distribution théorique ;
- distribution angulaire comparée à la théorie ;
- moyenne et quantiles dans les tolérances ;
- tests de convergence avec `N` croissant.

### 14.10 Régressions visuelles

Playwright doit produire des captures reproductibles pour un ensemble minimal :

```text
1s
2s
2p_z réel
2p état complexe m=+1
3d_z²
3d_xy
4f exemple
```

---

## 15. Références scientifiques de validation

Sources primaires / institutionnelles à privilégier :

### NIST / CODATA

- Fundamental Physical Constants: https://physics.nist.gov/constants
- CODATA version history: https://physics.nist.gov/cuu/Reference/versioncon.shtml
- Atomic Spectra Database: https://physics.nist.gov/asd
- Hydrogen energy levels: https://physics.nist.gov/PhysRefData/Handbook/Tables/hydrogentable5.htm

### NIST DLMF

Utiliser la Digital Library of Mathematical Functions comme référence pour :

- fonctions spéciales ;
- polynômes de Laguerre ;
- fonctions de Legendre associées ;
- conventions mathématiques.

- https://dlmf.nist.gov/

### Documentation interne

`docs/SCIENCE.md` devra préciser pour chaque formule :

- convention utilisée ;
- source ;
- unités ;
- domaine de validité ;
- test associé.

---

## 16. Tests logiciels et qualité

### 16.1 Vitest

Catégories :

```text
unit
scientific
sampling
regression
```

### 16.2 Playwright

Tester :

- démarrage application ;
- changement d'état ;
- génération concurrente ;
- responsive ;
- clavier ;
- thèmes ;
- module 2D ;
- absence d'erreurs console ;
- capture scientifique stable.

### 16.3 ESLint

Configuration stricte.

Interdire notamment :

- globals implicites ;
- promesses oubliées ;
- variables inutilisées ;
- conversions douteuses ;
- `any` TypeScript non justifié.

### 16.4 TypeScript

Activer un niveau strict :

```text
strict: true
noUncheckedIndexedAccess: true
exactOptionalPropertyTypes: true
noImplicitOverride: true
noFallthroughCasesInSwitch: true
```

Les types doivent différencier autant que possible :

- nombre quantique principal ;
- nombre azimutal ;
- nombre magnétique ;
- unités / grandeurs ;
- état complexe ;
- orbitale réelle ;
- résultat d'échantillonnage.

---

## 17. Performance

### 17.1 Objectifs

Cibles à mesurer, pas à supposer :

- UI toujours réactive pendant la génération ;
- annulation immédiate d'un job obsolète ;
- 60 FPS si la charge GPU le permet ;
- dégradation propre sur matériel plus faible ;
- aucune fuite de géométrie ou matériau Three.js.

### 17.2 Progressive rendering

Possibilité de retourner des lots progressifs :

```text
5 000
15 000
30 000
60 000 points
```

sans mélanger plusieurs états.

### 17.3 Qualité adaptative

Adapter éventuellement :

- pixel ratio ;
- nombre de points ;
- taille de points ;
- fréquence de mise à jour ;

mais ne jamais changer la physique pour gagner des FPS.

---

## 18. Import / export

### 18.1 Import JSON

Définir un schéma versionné.

Exemple :

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
  "points": []
}
```

Refuser un fichier dont l'unité est absente ou ambiguë.

### 18.2 Export

Priorité :

- capture PNG ;
- état JSON reproductible ;
- paramètres + seed ;
- données scientifiques CSV/JSON si utile.

OBJ / glTF ne doivent être annoncés que s'ils sont réellement implémentés et scientifiquement pertinents.

---

## 19. Documentation finale

### README.md

Doit contenir uniquement ce qui existe réellement :

- capture actuelle ;
- URL réelle ;
- installation réelle ;
- scripts npm réels ;
- version Three.js réelle ;
- limites scientifiques ;
- crédit de Kavan / `kavan010` ;
- statut de licence exact ;
- lien vers `SCIENCE.md` et `VALIDATION.md`.

### CREDITS.md

Doit conserver :

- Kavan / `kavan010` ;
- dépôt original ;
- rôle du travail original ;
- distinction entre origine pédagogique et réimplémentation web actuelle.

### SCIENCE.md

Contrat scientifique complet.

### VALIDATION.md

Résultats des tests numériques :

- tolérances ;
- valeurs analytiques ;
- comparaisons ;
- limites.

---

## 20. CI

Créer une GitHub Action sur `main` :

```text
npm ci
npm run typecheck
npm run lint
npm run test
npm run build
npm run test:e2e
```

Baseline : Node 24 LTS.

Aucun déploiement si un test scientifique critique échoue.

---

## 21. Séquence de refonte

## Phase 0 — Gel et sécurisation

### Actions

- créer le tag `legacy-web-2026-08-23` ;
- normaliser les fins de ligne avec `.gitattributes` ;
- confirmer `main` comme seule branche ;
- conserver l'historique ;
- créer `PLAN.md` ;
- créer `CREDITS.md` minimal ;
- retirer les déclarations de licence non vérifiées lors de la révision README.

### Critère de sortie

Le point de départ est récupérable exactement.

---

## Phase 1 — Toolchain moderne

### Actions

- initialiser npm ;
- installer Vite / TypeScript / Three.js ;
- créer la structure `src/` ;
- convertir le chargement CDN Three.js en import ESM ;
- ajouter ESLint / Prettier ;
- ajouter Vitest ;
- ajouter Playwright ;
- ajouter CI ;
- conserver le rendu fonctionnel avant modification physique majeure.

### Critère de sortie

L'application actuelle tourne sous Vite + TypeScript avec comportement visuel équivalent.

---

## Phase 2 — Moteur scientifique pur

### Actions

- créer `physics/constants.ts` ;
- créer `units.ts` ;
- écrire les fonctions spéciales ;
- écrire `Rₙₗ` ;
- écrire `Yₗᵐ` ;
- écrire les orbitales réelles ;
- supprimer progressivement les formules du renderer/UI ;
- écrire les premiers tests de normalisation.

### Critère de sortie

Le moteur quantique fonctionne sans DOM ni Three.js et passe ses tests analytiques.

---

## Phase 3 — Sampler scientifique

### Actions

- seedable RNG ;
- CDF radiale validée ;
- sampling angulaire validé ;
- suppression des samplers ad hoc ;
- tests statistiques ;
- unités uniques.

### Critère de sortie

Toutes les orbitales utilisent le même pipeline scientifique.

---

## Phase 4 — Worker et concurrence

### Actions

- déplacer la génération dans un Web Worker ;
- protocole `jobId` ;
- annulation des jobs obsolètes ;
- buffers transférables ;
- progression sûre ;
- test de changements rapides d'orbitale.

### Critère de sortie

Impossible qu'un résultat ancien écrase un nouvel état.

---

## Phase 5 — Renderer 3D

### Actions

- isoler Three.js ;
- gestion propre des ressources ;
- cadrage automatique ;
- échelle physique ;
- axes ;
- légendes ;
- mode densité ;
- mode phase ;
- préparation isosurfaces ;
- suppression des animations trompeuses.

### Critère de sortie

Le rendu est à la fois plus propre et plus fidèle à la signification physique.

---

## Phase 6 — UI/UX complète

### Actions

- séparer états complexes et orbitales réelles ;
- reconstruire les presets ;
- corriger toutes les notations ;
- moderniser le panneau ;
- corriger thème clair/sombre ;
- accessibilité clavier ;
- responsive ;
- panneaux pédagogiques courts mais précis.

### Critère de sortie

Aucun label UI ne contredit le moteur scientifique.

---

## Phase 7 — Interactions photon–hydrogène 2D

### Actions

- supprimer le modèle de Bohr utilisé comme représentation principale ;
- implémenter transitions énergétiques ;
- distinguer excitation et ionisation ;
- calculer l'énergie des photons ;
- documenter les approximations ;
- ajouter tests.

### Critère de sortie

Le module 2D possède une signification physique définie et vérifiable.

---

## Phase 8 — Validation globale

### Actions

- campagne de tests scientifiques ;
- comparaison aux références ;
- tests Playwright ;
- tests Chrome / Firefox / WebKit ;
- tests responsive ;
- profilage CPU/GPU ;
- chasse aux fuites mémoire.

### Critère de sortie

`VALIDATION.md` est rempli avec des résultats reproductibles.

---

## Phase 9 — Documentation, attribution et autonomie GitHub

### Actions

- réécrire README ;
- finaliser CREDITS ;
- clarifier licence ;
- supprimer le remote local `upstream` ;
- vérifier qu'`origin` ne possède que `main` ;
- décider du détachement du réseau de forks ;
- configurer GitHub Pages / déploiement réel ;
- ajouter badges réels uniquement.

### Critère de sortie

Le dépôt public décrit exactement ce qu'il contient et d'où il vient.

---

## Phase 10 — Release 5.0

### Versionnement proposé

```text
v5.0.0-alpha.1  moteur scientifique opérationnel
v5.0.0-beta.1   rendu + UI stabilisés
v5.0.0-rc.1     validation scientifique complète
v5.0.0          release publique
```

### Aucun `v5.0.0` si :

- une normalisation scientifique critique échoue ;
- les unités sont ambiguës ;
- `m` et orbitales réelles sont encore confondus ;
- la concurrence de génération n'est pas résolue ;
- le README contient encore des fonctions inexistantes ;
- le statut de licence reste présenté de manière fausse.

---

## 22. Definition of Done v5.0.0

- [ ] `main` est la seule branche permanente.
- [ ] Le point legacy est protégé par tag.
- [ ] L'origine Kavan / `kavan010` est clairement créditée.
- [ ] Le statut de licence est exact.
- [ ] TypeScript strict passe sans erreur.
- [ ] ESLint passe sans erreur bloquante.
- [ ] Tous les tests Vitest passent.
- [ ] Tous les tests scientifiques critiques passent.
- [ ] Tous les tests Playwright critiques passent.
- [ ] Toutes les longueurs utilisent une convention d'unité unique.
- [ ] La conversion `a₀ ↔ pm` est centralisée et testée.
- [ ] Aucun sampler orbital spécial incohérent ne subsiste.
- [ ] `m=±1` n'est plus étiqueté `p_x/p_y`.
- [ ] Les orbitales réelles sont générées comme combinaisons réelles correctes.
- [ ] Le nuage est explicitement décrit comme échantillonnage de probabilité.
- [ ] Aucune animation ne représente implicitement une trajectoire classique inexistante.
- [ ] Les jobs de génération concurrents sont sûrs et annulables.
- [ ] Le panneau 2D fonctionne réellement.
- [ ] Le module 2D possède un modèle physique explicite.
- [ ] Le rendu affiche une échelle physique.
- [ ] Le noyau schématique est signalé comme non à l'échelle.
- [ ] Le README ne décrit que des fonctions présentes.
- [ ] CI verte sur Node 24 LTS.
- [ ] `VALIDATION.md` documente les tolérances et les références.

---

## 23. Premier lot d'exécution recommandé

Ordre exact pour commencer :

1. créer le tag de sauvegarde ;
2. ajouter `.gitattributes` et stabiliser les fins de ligne ;
3. créer `package.json` ;
4. installer les versions modernes du socle ;
5. passer Three.js du CDN à un module npm ;
6. migrer `script.js` vers TypeScript **sans changer la physique** ;
7. ajouter Vitest et les premiers tests qui reproduisent l'état actuel ;
8. extraire le moteur mathématique dans `src/physics` ;
9. corriger les unités ;
10. supprimer les samplers spéciaux ;
11. séparer base complexe et orbitales réelles ;
12. déplacer l'échantillonnage dans un worker ;
13. seulement ensuite refaire le rendu et l'UI.

Cette séquence évite de modifier simultanément architecture, physique et apparence, ce qui rendrait les régressions difficiles à identifier.

---

## 24. Sources techniques des versions du socle

Vérifiées le 23 août 2026 :

- Node.js releases: https://nodejs.org/en/about/previous-releases
- TypeScript npm: https://www.npmjs.com/package/typescript
- Vite releases: https://vite.dev/releases
- Vite npm: https://www.npmjs.com/package/vite
- Three.js npm: https://www.npmjs.com/package/three
- Three.js migration guide: https://github.com/mrdoob/three.js/wiki/Migration-Guide
- Vitest npm: https://www.npmjs.com/package/vitest
- Playwright Test npm: https://www.npmjs.com/package/@playwright/test
- ESLint: https://eslint.org/
- Prettier 3.9: https://prettier.io/blog/2026/06/27/3.9.0.html

---

## 25. Décision d'architecture

La refonte `v5` part donc sur :

```text
TypeScript 7
Vite 8
Three.js r185
WebGL2 baseline
Web Worker pour le calcul
Vitest + Playwright
DOM natif, sans framework UI
moteur scientifique indépendant du renderer
unités atomiques internes
validation NIST/CODATA/DLMF
main uniquement + tags de versions
```

Ce document fait foi pour la refonte tant qu'une décision explicite ne le remplace pas.
