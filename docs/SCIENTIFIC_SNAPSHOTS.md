# Snapshots scientifiques Atoms

Atoms peut exporter et réimporter un état de travail dans un fichier JSON versionné. Ce format vise la reproductibilité d'une configuration scientifique et de sa présentation, pas l'archivage de buffers de calcul.

## Contrat v1

Le format courant est `atoms-scientific-snapshot` avec `schemaVersion: 1`. Son schéma JSON formel est publié dans [`schemas/atoms-scientific-snapshot-v1.schema.json`](schemas/atoms-scientific-snapshot-v1.schema.json).

Un snapshot contient six blocs de premier niveau :

- `format` et `schemaVersion`, qui identifient sans ambiguïté le contrat ;
- `units`, qui fixe explicitement les unités et conventions numériques ;
- `provenance`, qui identifie les versions des algorithmes utilisés ;
- `state.orbital`, qui décrit l'état quantique et sa base de représentation ;
- `state.sampling` et `state.rendering`, qui décrivent l'échantillonnage et le rendu ;
- `view.camera`, qui conserve l'orientation et la distance de la caméra.

Le JSON est sérialisé de façon déterministe et ne contient volontairement aucun timestamp. Deux exports du même état produisent donc le même contenu textuel.

## Rôle de la seed

La seed `uint32` n'est pas un paramètre physique de l'atome. C'est une clé de reproductibilité du tirage pseudo-aléatoire utilisé pour échantillonner `|ψ|² dV`.
À version identique du noyau scientifique, du sampler, du PRNG, du champ et du Worker, un même état orbital, un même nombre d'échantillons et une même seed reproduisent le même nuage Monte-Carlo. Cela permet de comparer des rendus, reproduire une figure et diagnostiquer une régression sans introduire un nouveau tirage aléatoire.

La garantie n'est pas portée par la seed seule. C'est pourquoi le snapshot enregistre aussi :

- `scienceEngineVersion` ;
- `samplerVersion` ;
- `prngVersion` ;
- `fieldVersion` ;
- `workerProtocolVersion` ;
- les paramètres numériques internes du Worker.

Un fichier provenant d'une provenance incompatible est refusé plutôt que converti silencieusement.

## Validation stricte

L'import accepte uniquement les champs prévus par le schéma. Les propriétés supplémentaires, unités inconnues, versions incompatibles, états quantiques invalides ou paramètres hors du domaine de l'interface provoquent une erreur explicite.

Le domaine importable correspond à ce qu'Atoms sait réellement représenter dans l'interface : `n` de 1 à 9, 2 000 à 60 000 échantillons et les plages/pas publiés pour les contrôles de rendu.

L'import applique d'abord la validation complète, régénère ensuite l'orbitale, puis restaure la caméra. Un import refusé ne doit pas être présenté comme une approximation de l'état demandé.

## Unités v1

- longueur : `bohr` (`a₀`) ;
- angle et phase : `radian` ;
- taille des points : `css-pixel` ;
- seuil d'isodensité : `fraction-of-grid-maximum`.

## Capture PNG

La commande de capture exporte la frame 3D courante en PNG depuis le canvas WebGL. Elle est destinée à la communication, aux rapports et aux comparaisons visuelles ; elle ne remplace pas le snapshot JSON, qui reste la source reproductible de l'état.

## Pourquoi pas de CSV

Atoms n'exporte pas de CSV dans ce lot. Un état orbital, sa provenance, ses unités, son rendu et sa caméra forment une structure hiérarchique pour laquelle JSON est le format adapté.

Un CSV ne deviendrait pertinent que pour un jeu de données tabulaire clairement défini, par exemple une série de mesures ou de points destinée à une analyse externe. Aucun export CSV n'est ajouté sans ce besoin scientifique explicite.

## Compatibilité future

Toute évolution incompatible du format doit créer une nouvelle version de schéma. Une future version pourra ajouter une migration explicite, mais la v1 n'effectue aucune conversion implicite d'unités, de bases orbitales ou de provenance numérique.

Le fichier TypeScript de référence est `src/state/scientificSnapshot.ts`. Le JSON Schema documente le format public ; les deux doivent évoluer ensemble et rester couverts par les tests.
