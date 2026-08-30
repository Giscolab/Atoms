# Contrat scientifique d'Atoms

Ce document décrit le noyau scientifique pur construit pendant les Lots 2, 2.1 et les Phases 3 à 7.
Le pipeline 3D visible consomme maintenant ce noyau via le sampler, le Web Worker et le renderer
scientifique. Le module photon–hydrogène 2D reste une démonstration historique séparée, explicitement
maintenue en Phase 8, et ne constitue pas la source du modèle orbital 3D.

## Domaine du modèle

Le modèle cible représente l'atome d'hydrogène neutre `¹H` comme un système électron-proton dans
l'approximation coulombienne non relativiste de Schrödinger. Le passage aux coordonnées du centre
de masse et de mouvement relatif remplace le problème à deux corps par un mouvement libre du
centre de masse et un problème coulombien relatif utilisant la masse réduite.

Ne sont pas inclus : structure fine ou hyperfine, spin-orbite, recul relativiste, taille finie du
proton, Lamb shift, QED, champs externes, effet Stark, effet Zeeman et données spectroscopiques
expérimentales. L'énergie calculée ne doit donc pas être présentée comme l'énergie exacte de
l'hydrogène réel.

Référence physique : MIT OpenCourseWare, _Quantum Physics I_, [Lectures 21–22: Hydrogen
Atom](references/MIT8_04S16_LecNotes22.pdf), notamment les équations (1.13), (1.17)–(1.24) et
(2.31). La [ressource officielle MIT](https://ocw.mit.edu/courses/8-04-quantum-physics-i-spring-2016/resources/mit8_04s16_lecnotes22/)
reste la référence institutionnelle en ligne.

## Constantes et provenance

Atoms utilise les valeurs recommandées de l'ajustement **CODATA 2022** publiées par le NIST dans
les tableaux XXXII et XXXIII de _CODATA recommended values of the fundamental physical constants:
2022_, DOI [`10.1063/5.0279860`](https://doi.org/10.1063/5.0279860). Le
[rapport NIST/CODATA](https://physics.nist.gov/cuu/pdf/JPCRD2022CODATA.pdf) et
[l'historique officiel des versions](https://physics.nist.gov/cuu/Reference/versioncon.shtml) sont
normatifs.

| Grandeur                           |   Symbole |              Valeur utilisée |   Unité |                   Incertitude-type | Statut SI           |
| ---------------------------------- | --------: | ---------------------------: | ------: | ---------------------------------: | ------------------- |
| Vitesse de la lumière dans le vide |       `c` |                `299 792 458` | `m s⁻¹` |                                `0` | exacte              |
| Constante de Planck                |       `h` |       `6.626 070 15 × 10⁻³⁴` |   `J s` |                                `0` | exacte              |
| Constante de Planck réduite        |       `ħ` |                     `h/(2π)` |   `J s` |                                `0` | exacte, dérivée     |
| Charge élémentaire                 |       `e` |      `1.602 176 634 × 10⁻¹⁹` |     `C` |                                `0` | exacte              |
| Masse de l'électron                |     `m_e` |     `9.109 383 7139 × 10⁻³¹` |    `kg` |                   `2.8 × 10⁻⁴⁰ kg` | non exacte, ajustée |
| Masse du proton                    |     `m_p` |   `1.672 621 925 95 × 10⁻²⁷` |    `kg` |                   `5.2 × 10⁻³⁷ kg` | non exacte, ajustée |
| Rapport électron-proton            | `m_e/m_p` |   `5.446 170 214 889 × 10⁻⁴` |     `1` | `9.4 × 10⁻¹⁵` (`uᵣ = 1.7 × 10⁻¹¹`) | non exact, ajusté   |
| Rayon de Bohr conventionnel        |      `a₀` |   `5.291 772 105 44 × 10⁻¹¹` |     `m` |                    `8.2 × 10⁻²¹ m` | non exact, dérivé   |
| Énergie de Hartree                 |     `E_h` | `4.359 744 722 2060 × 10⁻¹⁸` |     `J` |                    `4.8 × 10⁻³⁰ J` | non exacte, dérivée |

« Exacte dans le SI » décrit le statut métrologique de la grandeur ; cela ne signifie pas que sa
valeur est toujours représentable exactement par un `number` IEEE-754. `ħ` est donc calculée depuis
`h/(2π)` au lieu de recopier une décimale tronquée.

La constante de Rydberg n'est pas introduite dans ce lot : `E_h` suffit au modèle énergétique et
évite une constante redondante. Aucune incertitude n'est dérivée pour `μ`, `a_μ` ou les niveaux
d'énergie sans traitement des corrélations de l'ajustement CODATA.

Le rapport sans dimension `m_e/m_p` provient directement de la table XXXIII de CODATA 2022. Pour
le modèle numérique de `¹H`, cette valeur ajustée est la source de vérité du rapport de masses.

## Masse réduite et échelle de longueur

Le modèle `¹H` utilise explicitement la masse réduite électron-proton :

```text
μ = m_e m_p / (m_e + m_p)
```

Pour le modèle numérique particulier de `¹H`, Atoms utilise directement le rapport ajusté CODATA :

```text
r = m_e/m_p
μ/m_e = 1 / (1 + r)
```

Cela évite de reconstruire une grandeur mieux déterminée à partir de deux masses SI arrondies et
corrélées. La fonction générique de masse réduite conserve néanmoins la relation physique générale
ci-dessus.

Le rayon de Bohr CODATA reste la grandeur conventionnelle définie avec `m_e`. Le rayon
caractéristique du problème relatif est une grandeur distincte :

```text
a_μ = a₀ m_e / μ
a_μ = a₀ (1 + m_e/m_p)
```

Il en résulte `0 < μ < m_e` et `a_μ > a₀`. Avec les valeurs centrales retenues,
`a_μ ≈ 52.946540946024626 pm`, tandis que `1 a₀ = 52.9177210544 pm`.

## Politique d'unités

Les calculs analytiques emploient les unités atomiques lorsque cela est pertinent : longueur en
`a₀`, énergie en `E_h`, angles en radians. Les conversions sont explicites aux frontières :

- `a₀ ↔ m` et `a₀ ↔ pm` ;
- `m ↔ nm` ;
- `E_h ↔ J` et `E_h ↔ eV`.

La conversion électronvolt–joule découle de la charge élémentaire exacte :

```text
1 eV = e joules
```

Les anciennes valeurs magiques `5.29`, `52.9`, `27.2` et `13.6` n'appartiennent pas au nouveau
socle.

## Nombres quantiques

Un état quantique admissible respecte :

```text
n entier, n >= 1
l entier, 0 <= l <= n - 1
m entier, -l <= m <= l
```

Les entrées hors de ce domaine sont rejetées explicitement avant tout calcul.

## Énergie de Schrödinger de ¹H

L'énergie coulombienne non relativiste utilise le Hartree conventionnel et la correction explicite
de masse réduite :

```text
E_n / E_h = -(μ/m_e) / (2 n²)
μ/m_e = 1 / (1 + m_e/m_p)
```

Elle est exposée en Hartree, joules et électronvolts. Elle suit `E_n ∝ -1/n²`, mais reste distincte
des valeurs expérimentales ou recommandées des niveaux et de l'énergie d'ionisation, qui incluent
des phénomènes absents du modèle.

## Polynômes de Laguerre généralisés

Atoms adopte les polynômes `L_k^(α)(x)` de la convention NIST DLMF. Le domaine exposé est `k`
entier non négatif, `α > -1` réel et fini, et `x` réel fini. Les orbitales hydrogénoïdes du lot
suivant utiliseront `α = 2l + 1`.

Initialisation et récurrence :

```text
L_0^(α)(x) = 1
L_1^(α)(x) = 1 + α - x
(k + 1)L_(k+1)^(α)(x)
  = (2k + α + 1 - x)L_k^(α)(x) - (k + α)L_(k-1)^(α)(x)
```

Références : NIST DLMF [§18.5](https://dlmf.nist.gov/18.5), cas analytiques, et
[§18.9, table 18.9.1](https://dlmf.nist.gov/18.9#T1), récurrence.

## Fonctions de Legendre associées

Sur `x ∈ [-1, 1]`, le module calcule la fonction de Ferrers de première espèce `P_l^m(x)` pour
`l >= 0` et `0 <= m <= l`. La convention inclut explicitement la phase de Condon–Shortley :

```text
P_l^m(x) = (-1)^m (1 - x²)^(m/2) d^m P_l(x)/dx^m
```

En particulier :

```text
P_1^1(x) = -sqrt(1 - x²)
```

Références : NIST DLMF [§14.7.8](https://dlmf.nist.gov/14.7.E8), définition avec phase, et
[§14.10.3](https://dlmf.nist.gov/14.10.E3), récurrence en degré. Les notes MIT
[Lectures 20–21](references/MIT8_04S16_LecNotes20_21.pdf) placent leur facteur `(-1)^m` au niveau
des harmoniques sphériques ; Atoms choisit au contraire de l'intégrer dès `P_l^m`, conformément à
la convention DLMF retenue. Un futur module d'harmoniques ne devra donc pas l'appliquer une seconde
fois.

## Fonction radiale normalisée

Le calcul spatial reçoit `rBohr = r/a₀`, une distance réelle positive ou nulle exprimée en rayons de
Bohr conventionnels. L'échelle du problème relatif de `¹H` est représentée dans la même unité par :

```text
a = a_μ/a₀ = 1/(μ/m_e)
ρ = 2 rBohr/(n a)
k = n-l-1
α = 2l+1
```

La fonction numérique renvoyée représente `a₀^(3/2) R_nl^SI` et porte donc le contrat dimensionnel
`a₀^(-3/2)` :

```text
R_nl(rBohr) = sqrt[
  (2/(n a))³ (n-l-1)! / (2n (n+l)!)
] exp(-ρ/2) ρ^l L_(n-l-1)^(2l+1)(ρ)
```

Sa normalisation est définie dans la coordonnée numérique du moteur :

```text
∫₀∞ rBohr² |R_nl(rBohr)|² drBohr = 1
```

La formule et sa mesure sont celles de NIST DLMF
[§18.39(ii), équations 18.39.35 et 18.39.37](https://dlmf.nist.gov/18.39.E37), avec `Z=1` et
l'échelle conventionnelle remplacée explicitement par `a_μ`. Les notes MIT
[Hydrogen Atom](references/MIT8_04S16_LecNotes22.pdf), équations (2.28), (2.33)–(2.36), confirment
la séparation, le degré polynomial et l'exponentielle. Leur variable radiale vaut `r/(na)`, soit
`ρ/2` dans le présent contrat ; cette différence de notation ne change pas la fonction.

Les formes analytiques ponctuelles utilisées comme oracles sont :

```text
R_10 = 2 a^(-3/2) exp(-rBohr/a)
R_20 = [1/(2 sqrt(2))] a^(-3/2)
       (2-rBohr/a) exp[-rBohr/(2a)]
R_21 = [1/(2 sqrt(6))] a^(-3/2)
       (rBohr/a) exp[-rBohr/(2a)]
```

## Harmoniques sphériques complexes

Atoms suit NIST DLMF [§14.30.1](https://dlmf.nist.gov/14.30.E1). Pour `m >= 0`, `theta` dans
`[0, π]` et `phi` réel fini en radians :

```text
N_lm = sqrt[(2l+1)/(4π) (l-m)!/(l+m)!]
Y_l^m(theta,phi) = N_lm P_l^m(cos theta) exp(i m phi)
```

`P_l^m` est la fonction de Ferrers du module précédent et contient déjà `(-1)^m`. Le module
`Y_l^m` **n'ajoute donc aucune seconde phase de Condon–Shortley**. Cette décision suit directement
NIST DLMF [§14.7.8](https://dlmf.nist.gov/14.7.E8). Les notes MIT
[Quantum Mechanics in 3D / Angular Momentum](references/MIT8_04S16_LecNotes20_21.pdf), équations
(3.26)–(3.27), utilisent une autre répartition du même signe : leur `P_l^m` n'inclut pas la phase et
leur définition de `Y_l^m` l'ajoute ensuite.

Les ordres négatifs ne sont jamais transmis au module de Ferrers. Pour `m>0`, ils sont construits
par la relation DLMF [14.30.6](https://dlmf.nist.gov/14.30.E6) :

```text
Y_l^(-m) = (-1)^m conjugate(Y_l^m)
```

Les harmoniques sont sans dimension et orthonormées avec la mesure solide explicite :

```text
∫₀²π ∫₀π conjugate(Y_l^m) Y_l'^m' sin(theta) dtheta dphi
  = δ_ll' δ_mm'
```

La parité adoptée est celle de DLMF [14.30.7](https://dlmf.nist.gov/14.30.E7) :

```text
Y_l^m(π-theta, phi+π) = (-1)^l Y_l^m(theta,phi)
```

`phi` est accepté sur tout le domaine réel fini. Sa réduction périodique interne avant le produit
`m phi` préserve exactement la convention tout en évitant un débordement numérique inutile.

## Fonction d'onde, densité et phase

La fonction d'onde stationnaire séparée est :

```text
ψ_nlm(rBohr,theta,phi) = R_nl(rBohr) Y_l^m(theta,phi)
```

Elle est complexe, normalisée dans la mesure `rBohr² sin(theta) drBohr dtheta dphi` et possède le
contrat `a₀^(-3/2)`. La densité volumique est strictement :

```text
|ψ_nlm|² = ψ_nlm conjugate(ψ_nlm)
```

Elle est réelle, positive ou nulle, exprimée en `a₀^(-3)` et ne contient ni `rBohr²` ni
`sin(theta)`. Ces facteurs sont la jacobienne de la mesure sphérique et ne doivent pas être appliqués
deux fois.

La phase est `atan2(Im ψ, Re ψ)` en radians dans `[-π, π]`. Lorsque les deux composantes de `ψ`
sont exactement nulles, la phase est physiquement indéfinie et l'API renvoie `null` plutôt qu'un
angle arbitraire.

## Orbitales réelles standard

Une orbitale réelle n'est pas un renommage d'une valeur particulière de `m`. À partir des
harmoniques complexes normalisées et pour `m>0`, Atoms fixe les combinaisons :

```text
C_lm = [Y_l^(-m) + (-1)^m Y_l^m]/sqrt(2)
S_lm = i [Y_l^(-m) - (-1)^m Y_l^m]/sqrt(2)
```

Les signes globaux sont déterministes et choisis pour donner des coefficients cartésiens positifs :

```text
p_x = C_11          p_y = S_11          p_z = Y_1^0
d_xz = C_21         d_yz = S_21         d_z2 = Y_2^0
d_x2_y2 = C_22      d_xy = S_22
```

Il en résulte notamment `p_x ∝ x/r`, `p_y ∝ y/r`, `p_z ∝ z/r`, ainsi que les formes usuelles
positives en `xy`, `xz`, `yz`, `x²-y²` et `3z²-r²`. Chaque combinaison conserve la normalisation et
utilise exactement la même fonction radiale pour un couple `(n,l)` donné. La représentation
complexe générique accepte tout `l` numériquement représentable ; seuls les noms réels explicites
des familles `p` et `d`, nécessaires à cette phase, sont introduits. Aucune famille réelle `f/g/h/i`
n'est ajoutée sans besoin actuel.

## Nœuds et observables analytiques

Pour un état admissible `(n,l)`, la source de vérité scientifique est :

```text
radialNodes = n-l-1
angularNodes = l
totalNodes = n-1
```

Le zéro de `R_nl(0)` dû au facteur `r^l` pour `l>0` n'est pas compté comme un nœud radial
supplémentaire. Les notes MIT _Hydrogen Atom_, figure 2 et équation (2.34), identifient le degré
`n-l-1` du polynôme radial et son nombre de nœuds.

Les seules observables ajoutées en Phase 3 sont exprimées dans le contrat interne :

```text
<r>/a₀ = (a/2) [3n²-l(l+1)]
<1/r> a₀ = 1/(a n²)
r_most_probable(1s)/a₀ = a
```

Les deux espérances sont les formules hydrogénoïdes usuelles, documentées notamment dans le
[corrigé institutionnel MIT 5.80, problème 7](https://ocw.mit.edu/courses/5-80-small-molecule-spectroscopy-and-dynamics-fall-2008/a6931596f907971b77a92bf01588bb6b_02pset_ans_sp94.pdf),
adaptées de `a₀` à `a_μ`. Le maximum radial `1s` découle directement de la maximisation de
`r²|R_10|²` ; il ne s'agit ni du maximum de `|ψ|²` ni d'une trajectoire.

## Validation et domaine numériques

Les quadratures de validation sont déterministes et restent exclusivement dans
`tests/scientific/numericalIntegration.ts`. Elles utilisent Simpson composite en rayon et en
`theta`, ainsi qu'une règle périodique en `phi`. Les jacobiennes `rBohr²` et `sin(theta)` y sont
appliquées explicitement, jamais dans la densité de production. Les seuils partagés sont :

```text
identités analytiques composante par composante : 512 Number.EPSILON
normalisation numérique                         : 1e-8
orthogonalité numérique                         : 1e-8
normalisation par quadrature 3D                 : 1e-7
stabilité entre deux maillages 3D               : 5e-7
```

Les tests radiaux vérifient séparément la stabilité lors du raffinement du maillage et de
l'extension du domaine fini. Les tests angulaires et 3D comparent également deux résolutions. Ces
seuils sont des bornes absolues adaptées à des intégrales visant respectivement un ou zéro ; ils ne
sont pas réutilisés pour les valeurs CODATA.

Les fonctions rejettent les nombres non finis, les rayons négatifs, les états quantiques invalides
et les angles polaires hors de `[0, π]`. `phi` accepte tout réel fini et est traité modulo `2π`. Les
normalisations réutilisent la factorielle IEEE-754 du socle : une fonction radiale exige
`n+l <= 170`, et une harmonique exige `l+|m| <= 170`. Au-delà, l'API échoue explicitement plutôt que
de retourner une valeur sous-évaluée ou non finie.

## Pipeline 3D intégré (Phases 4 à 7)

La génération visible suit un pipeline unique et versionné :

```text
état complexe ou orbitale réelle
        ↓
ψ, |ψ|² et phase du noyau Phase 3
        ↓
CDF radiale r²|Rₙₗ|² + sampling angulaire avec sin(θ) dθ dφ
        ↓
Web Worker (seed, jobId, versions, buffers transférables)
        ↓
nuage probabiliste, champ cartésien et courbes analytiques
        ↓
renderer Three.js / UI DOM
```

Le nuage affiché contient des positions tirées selon la mesure volumique
`|ψ|² dV` ; ses points ne sont ni des électrons individuels ni des trajectoires. La seed
`uint32` est déterministe à version identique du moteur. Les anciennes animations de points ont été
retirées du rendu 3D ; seule la rotation de caméra est proposée comme outil d'observation.

Le Worker utilise `radialCoverageProbability = 0,999` uniquement pour choisir l'étendue finie du
cube de visualisation à partir de la CDF radiale. Cette convention de cadrage ne tronque ni ne
renormalise la fonction d'onde analytique. Le champ cartésien est une grille numérique de résolution
`32³` par défaut, indexée `x + N(y + Nz)`, et ses densités sont normalisées par le maximum de la
grille pour l'affichage ; les maxima bruts en `a₀⁻³` et `a₀⁻³ᐟ²` restent fournis dans le contrat.

Les courbes latérales séparent les grandeurs : la distribution radiale trace `r²|Rₙₗ|²` en
`a₀⁻¹`, tandis que la coupe angulaire trace `|Yₗᵐ|²` sur le meilleur des grands cercles `xy`, `xz`
ou `yz`. Cette dernière est une coupe géométrique, pas une probabilité angulaire intégrée.

Le renderer propose trois modes d'affichage : densité, phase et hybride. La densité détermine la
géométrie d'une isosurface positive ; la phase est colorée par une palette divergente cyan/corail
avec la convention `arg(ψ) ∈ [-π, π]`. Les surfaces `ψ = 0` ne sont activées que lorsque le champ
réel est interprétable (`m = 0` complexe ou orbitale réelle) et restent distinctes de l'isosurface
de densité. Le noyau est schématique et explicitement indiqué comme non à l'échelle.

## Interprétation des représentations graphiques

- **Points du nuage** : chaque point est un échantillon de la distribution `|ψ|² dV`. Un nuage
  fini est une représentation de Monte-Carlo ; il ne représente ni plusieurs électrons, ni la
  position certaine d'un électron, ni une trajectoire temporelle.
- **Isosurface** : la surface relie les points où la densité volumique normalisée sur la grille
  atteint le seuil d'affichage choisi. Elle ne constitue ni une frontière matérielle de l'atome,
  ni une orbite, ni une surface de charge.
- **Couleur de phase** : les pôles cyan et corail codent un écart de phase de `ψ`. Pour une
  orbitale réelle, ils permettent notamment de distinguer les deux signes, équivalents à des phases
  séparées de `π` dans la convention globale choisie. Pour un état complexe général, la palette
  parcourt cycliquement les phases intermédiaires. Ces couleurs ne désignent donc ni deux charges,
  ni deux électrons différents.
- **Orbitales réelles `d`** : `d_xy` et `d_x²−y²` sont des combinaisons linéaires réelles,
  normalisées et distinctes, construites à partir des composantes complexes `m = ±2`. Elles ne sont
  pas le renommage d'un état complexe possédant une valeur unique de `m`.

Les [captures documentaires](captures/) illustrent ces conventions avec le même moteur orbital 3D.
Dans la capture où le panneau 2D historique est ouvert, seule la visualisation orbitale 3D relève de
ce pipeline : le panneau 2D reste un module indépendant maintenu pour la Phase 8. Les légendes et le
rôle de chaque image sont détaillés dans la section
[Captures d'écran du README](../README.md#captures-décran).

`legacyScience.ts` est désormais limité à la constante du panneau photon–hydrogène 2D historique.
Il n'est pas utilisé pour la génération orbitale 3D. Toute évolution de ce module 2D doit rester dans
la Phase 8 et documenter séparément ses hypothèses physiques.
