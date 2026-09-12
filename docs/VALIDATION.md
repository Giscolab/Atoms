# Validation d'Atoms

Ce document décrit les validations reproductibles du moteur scientifique et de l'application web.
Il complète [SCIENCE.md](SCIENCE.md) et [REFERENCES.md](REFERENCES.md). Les résultats de performance
sont des mesures locales datées, pas des garanties universelles.

## Commandes de qualification

```bash
npm audit
npm run lint
npm run typecheck
npm run typecheck:science
npm test
npm run build
npm run test:e2e
npm run test:qualification
npm run format:check
git diff --check
```

La CI GitHub audite les dépendances, vérifie le formatage et les contrôles de qualité, puis exécute les E2E sous Chromium, Firefox et WebKit.
La qualification performance reste séparée car ses temps dépendent du matériel et du navigateur.

## Validation scientifique

La suite Vitest contient actuellement 310 tests dans 28 fichiers. Elle couvre notamment les nombres
quantiques, constantes CODATA, unités, énergie, fonctions spéciales, partie radiale, harmoniques
sphériques complexes, orbitales réelles, fonction d'onde, observables, sampling, Worker et le
contrat des snapshots scientifiques versionnés.

### Tolérances numériques

Les tolérances partagées sont définies dans `tests/scientific/numericAssertions.ts` :

| Validation                                      |              Tolérance |
| ----------------------------------------------- | ---------------------: |
| courte chaîne IEEE-754 relative                 |  `16 × Number.EPSILON` |
| polynômes analytiques de faible degré           |  `64 × Number.EPSILON` |
| identités analytiques composante par composante | `512 × Number.EPSILON` |
| normalisation numérique                         |                 `1e-8` |
| orthogonalité numérique                         |                 `1e-8` |
| normalisation par quadrature 3D                 |                 `1e-7` |
| convergence entre deux maillages 3D             |                 `5e-7` |

Les quadratures de test utilisent Simpson composite en rayon et en `theta`, une règle périodique en
`phi`, et appliquent explicitement les jacobiens `r²` et `sin(theta)`. Les tests comparent aussi des
maillages de résolutions différentes afin de distinguer précision ponctuelle et convergence.

### Sampling statistique

Les tests statistiques utilisent une borne de Dvoretzky–Kiefer–Wolfowitz avec la constante de
Massart et correction de Bonferroni. Le budget de faux rejet de la famille de tests est `1e-6`.
La seed `uint32` rend les campagnes déterministes à version identique du moteur.

Les tests du sampler vérifient les distributions radiales et angulaires, les moments attendus, les
valeurs finies, les domaines autorisés et la reproductibilité.

### Snapshots scientifiques

Le format `atoms-scientific-snapshot` v1 est couvert par 17 tests unitaires dédiés. Ils vérifient le
round-trip déterministe, les unités, la provenance numérique, les paramètres internes du Worker,
les champs inconnus, les bases orbitales, le domaine `n = 1…9`, les plages/pas représentables par
l'interface, la caméra et les erreurs JSON.

Le flux navigateur est qualifié sous Chromium, Firefox et WebKit : export JSON puis réimport réel,
refus d'une version incompatible sans altération de l'état courant, et export d'un PNG dont la
signature binaire est contrôlée. Le schéma public est documenté dans
[`SCIENTIFIC_SNAPSHOTS.md`](SCIENTIFIC_SNAPSHOTS.md) et
[`schemas/atoms-scientific-snapshot-v1.schema.json`](schemas/atoms-scientific-snapshot-v1.schema.json).

## Validation navigateur et accessibilité

Playwright qualifie l'application sous Chromium, Firefox et WebKit. Les scénarios couvrent le
chargement WebGL2, la génération scientifique, les bases complexe/réelle, la seed, les thèmes,
les nœuds, le responsive mobile, le défilement compact et les interactions principales.

Sur les runners Linux GitHub Actions dépourvus de GPU matériel, le projet Firefox de Playwright
active `webgl.force-enabled=true` afin d’autoriser le backend logiciel WebGL2. Ce réglage appartient
uniquement à l’environnement de qualification ; il ne modifie ni le code produit ni le site déployé.

La navigation clavier est testée dans les deux thèmes et les deux bases. Le parcours vérifie les
contrôles quantiques, l'observable, le mode d'affichage, l'échantillonnage, les toggles, la
génération, la caméra, le canvas et le panneau d'analyses. Les éléments atteints par `Tab` doivent
présenter un indicateur de focus visible et rester dans le viewport.

L'audit axe-core est exécuté sous Chromium pour les thèmes sombre et clair avec les tags WCAG A/AA.
Au checkpoint de validation, aucune violation automatisable n'est remontée. Un contrôle automatisé
ne remplace pas une revue humaine exhaustive de l'ergonomie et de la compréhension visuelle.

## Régressions scientifiques déterministes

Sept baselines SVG versionnées couvrent `1s`, `2s`, `2p_z` réel, `2p` complexe `m=+1`, `3d_z²`
réel, `3d_xy` réel et `4f` complexe. Elles utilisent une seed fixe, 2 000 points, un thème et une
caméra stabilisés. Les captures comparent les courbes scientifiques SVG, pas le raster WebGL.

Les baselines sont qualifiées sous Chromium pour éviter les différences de sérialisation CSS sans
valeur scientifique entre moteurs. Les scénarios fonctionnels restent, eux, multi-navigateurs.
Une seconde génération Worker du même état doit reproduire exactement la même capture scientifique.

## Qualification performance et ressources

Mesures locales du 12 septembre 2026 : Windows NT `10.0.26200`, Intel Core i7-14700F, 31,8 Gio de
RAM, Node.js `26.8.1`, Playwright `1.62.1`. Chromium headless utilise ANGLE/SwiftShader dans ce run ;
les temps et la cadence ne représentent donc pas une promesse de performance GPU matérielle.

| Points | Worker → résultat, médiane | État prêt, médiane | Buffers transférés |
| -----: | -------------------------: | -----------------: | -----------------: |
|  2 000 |                    17,0 ms |            23,5 ms |     465 104 octets |
| 15 000 |                    47,4 ms |            54,2 ms |     829 104 octets |
| 60 000 |                 1 022,7 ms |         1 047,0 ms |   2 089 104 octets |

Le stress a observé 76 jobs : 58 résultats, 18 annulations par supersession, zéro erreur et un seul
Worker actif au maximum. Après warm-up, quatre cycles revenus au même état conservent exactement
6 géométries, 6 matériaux, 8 programmes et 2 textures Three.js, ainsi que 15 buffers et 6 textures
WebGL comptés par le probe.

Le heap V8 principal mesuré après collecte explicite passe d'environ 10,58 à 10,78 Mio sur cinq
cycles. Cette variation bornée n'est ni une preuve de fuite ni une preuve universelle d'absence de
fuite. Aucun accroissement des compteurs de ressources Three.js/WebGL n'a été observé au même état.

La cadence `requestAnimationFrame` observée est d'environ 30 FPS sous SwiftShader. Elle inclut le
navigateur et le compositeur et ne mesure pas directement le temps GPU.

## Qualité TypeScript et dépendances

`tsconfig.json` active notamment `strict`, `noUncheckedIndexedAccess`, `exactOptionalPropertyTypes`,
`noImplicitOverride`, `noUnusedLocals` et `noUnusedParameters`. Le noyau scientifique est aussi
compilé séparément sans bibliothèques DOM via `tsconfig.science.json`.
L'audit statique ciblé ne trouve aucun `any`, `@ts-ignore` ni `@ts-expect-error` dans `src/` et
`tests/`. Les rares assertions doubles restantes sont cantonnées aux frontières Worker ou aux
fixtures de test et restent soumises au typecheck et à ESLint.

`npm audit` doit rester à zéro vulnérabilité connue au checkpoint de release. Les versions de
production et de développement sont verrouillées par `package-lock.json` et la CI utilise `npm ci`.

## Limites de la validation

Ces tests établissent des invariants sur le domaine couvert ; ils ne constituent pas une preuve
formelle du modèle quantique ni de toutes les configurations possibles de l'interface.

- le modèle physique est l'hydrogène non relativiste à masse réduite décrit dans `SCIENCE.md` ;
- les quadratures numériques couvrent des états et maillages représentatifs, pas un continuum infini ;
- le rendu utilise une grille finie et des buffers `Float32` et ne prouve pas la convergence de toute
  géométrie d'isosurface ;
- les points Monte-Carlo ne représentent ni des électrons individuels ni des trajectoires ;
- l'audit axe-core ne remplace pas une évaluation humaine complète de l'accessibilité ;
- les mesures de temps, heap, cadence et ressources sont bornées à l'environnement de qualification ;
- les baselines SVG valident les courbes et conventions scientifiques déterministes, pas chaque pixel
  du raster WebGL sur toutes les piles GPU.

Pour les formules, constantes et sources institutionnelles, voir [SCIENCE.md](SCIENCE.md) et
[REFERENCES.md](REFERENCES.md).
