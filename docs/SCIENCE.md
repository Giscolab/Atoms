# Contrat scientifique d'Atoms

Ce document décrit uniquement le nouveau socle scientifique introduit au Lot 2. Le moteur
historique reste temporairement responsable de l'application visible et n'est pas une source de
vérité scientifique.

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

## Frontière du Lot 2

Ce lot n'implémente pas les fonctions radiales `R_nl`, les harmoniques sphériques `Y_l^m`, les
fonctions d'onde `ψ_nlm`, la densité `|ψ|²`, les orbitales réelles ou l'échantillonnage. Aucun module
visible n'importe encore ce nouveau socle. L'application continue d'utiliser exclusivement
`legacyScience.ts` jusqu'à un lot de migration explicitement autorisé.
