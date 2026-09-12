# Rendu des isodensités

La surface représente un niveau de `|ψ|²`, normalisé par le maximum du champ cartésien
reçu du Worker. Son matériau, ses reflets et sa transparence sont des conventions
graphiques ; ils ne représentent aucune propriété matérielle ou électrique de l’atome.

## Extraction

`src/rendering/isosurfaceGeometry.ts` raffine chaque intervalle de la grille en deux
par interpolation trilinéaire de la **densité**. Les 32 sommets scientifiques par axe
donnent ainsi 63 sommets de tessellation. Une bordure technique recopiée permet à
MarchingCubes de traiter aussi les cellules extérieures du domaine original ; aucune
surface artificielle ne ferme le cube lorsque le niveau atteint son bord.

L’extraction utilise le [MarchingCubes de Three.js](https://threejs.org/docs/pages/MarchingCubes.html),
sans billes, floutage ou relaxation des sommets. Les positions sont converties en `a₀`.
Les normales unitaires interpolent le gradient négatif de la grille originale,
calculé par différences centrées à l’intérieur et unilatérales au bord. Elles servent
à l’éclairage, sans déplacer les sommets. Les triangles dégénérés sont exclus ; les
buffers finaux ne contiennent que les sommets utilisés.

Cette subdivision approxime plus finement le niveau de l’interpolant trilinéaire.
Elle ne récupère pas les structures absentes de la grille scientifique et ne démontre
pas une convergence vers la fonction d’onde analytique. Les ambiguïtés topologiques
du MarchingCubes classique restent une limite pour des champs arbitraires. Les
cas de référence testés ne constituent pas une garantie pour tous les états et seuils.

## Phase et lumière

La géométrie reste identique lorsque l’observable ou le thème change. Pour un champ
réel, le signe de l’amplitude interpolée sélectionne les phases opposées `0` et `π`.
Pour un champ complexe, la phase est interpolée dans le plan complexe, avec une
pondération par l’amplitude ; une moyenne directe des angles serait incorrecte à
la couture `−π/π`. Les couleurs sont converties de sRGB vers l’espace linéaire de
travail avant leur utilisation par Three.js. La palette et sa légende sont conservées.

Le [MeshPhysicalMaterial](https://threejs.org/docs/pages/MeshPhysicalMaterial.html) dédié
emploie une rugosité modérée et un vernis discret, sans métal, transmission ou émission.
Trois lumières directionnelles et une lumière hémisphérique accompagnent l’orientation
de la caméra. Le tone mapping neutre contient les reflets ; le nuage conserve ses
couleurs sans tone mapping. Le mode isodensité utilise une peau opaque avec profondeur.
Le mode hybride conserve une peau transparente à face unique au-dessus du nuage,
évitant l’accumulation des faces avant et arrière de chaque lobe.

Aucun bloom, flou de profondeur, déplacement de surface ou animation supplémentaire
n’est ajouté. Les commandes existantes, le thème clair et la réduction des mouvements
restent applicables.

## Compatibilité

Le noyau scientifique, le sampler, les trois buffers du champ, la résolution Worker
`32³`, ses autres options et toutes les versions numériques sont inchangés. Les
snapshots scientifiques v1 restent compatibles. Ils reproduisent les données et
l’état de présentation ; ils ne figent pas les pixels d’une ancienne version du renderer.

Les preuves numériques, navigateur et performance figurent dans [VALIDATION.md](VALIDATION.md).
