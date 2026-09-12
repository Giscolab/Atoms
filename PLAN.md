# Atoms — Plan restant

Ce document ne recense plus les étapes déjà réalisées. Le travail validé appartient à l’historique Git, au README et à la documentation scientifique.

**Projet :** `Giscolab/Atoms`
**Cible :** finalisation `v5.0.0` et préparation Showcase
**Règle :** une tâche disparaît de ce fichier dès qu’elle est implémentée et validée.

---

# 1. Validation globale

## 1.1 Régressions visuelles

- [ ] ajouter des régressions visuelles déterministes ;
- [ ] couvrir au minimum `1s`, `2s`, `2p_z`, `2p` complexe `m=+1`, `3d_z²`, `3d_xy` et un état `4f` représentatif ;
- [ ] conserver des captures scientifiques stables et reproductibles.

## 1.2 Navigateurs

- [ ] exécuter et qualifier les scénarios E2E sous Firefox ;
- [ ] exécuter et qualifier les scénarios E2E sous WebKit.

## 1.3 Performance et ressources

- [ ] profiler CPU ;
- [ ] profiler GPU ;
- [ ] mesurer le coût de conversion `Float64 → Float32` avant rendu ;
- [ ] vérifier l’absence de fuite mémoire lors de changements répétés d’état ;
- [ ] vérifier l’absence de fuite de ressources Three.js ;
- [ ] documenter les mesures reproductibles avant toute optimisation.

## 1.4 Accessibilité finale

- [ ] valider la navigation clavier complète ;
- [ ] auditer le focus visible ;
- [ ] vérifier les labels explicites et l’usage pertinent de `aria-live` ;
- [ ] vérifier les contrastes ;
- [ ] vérifier qu’aucun contrôle n’est identifié uniquement par la couleur ;
- [ ] confirmer que la palette reste exploitable autant que possible en cas de déficience de perception des couleurs.

## 1.5 Validation scientifique documentée

- [ ] créer puis remplir `docs/VALIDATION.md` ;
- [ ] y documenter les tolérances numériques utilisées ;
- [ ] y documenter les valeurs analytiques de référence ;
- [ ] y documenter normalisations, orthogonalités et tests statistiques ;
- [ ] y documenter les navigateurs qualifiés, les performances mesurées et les limites connues ;
- [ ] vérifier que les résultats publiés sont reproductibles.

---

# 2. Documentation et statut juridique

- [ ] effectuer l’audit final de `docs/SCIENCE.md` ;
- [ ] auditer `docs/REFERENCES.md` avant release ;
- [ ] clarifier le statut de licence exact du dépôt ;
- [ ] vérifier que README, `SCIENCE.md`, `REFERENCES.md` et futur `VALIDATION.md` décrivent uniquement des fonctions réellement présentes ;
- [ ] vérifier que les limites scientifiques sont explicites ;
- [ ] vérifier que provenance, attribution et licence restent formulées sans ambiguïté.

---

# 3. Import / export scientifique

## 3.1 Import JSON

- [ ] définir une version de schéma obligatoire ;
- [ ] déclarer explicitement les unités ;
- [ ] valider l’état quantique importé ;
- [ ] déclarer explicitement la base complexe ou réelle ;
- [ ] rendre la seed explicite lorsque le format l’exige ;
- [ ] refuser les unités ambiguës.

## 3.2 Export

- [ ] permettre une capture PNG ;
- [ ] permettre un état JSON reproductible ;
- [ ] inclure paramètres et seed ;
- [ ] ajouter un export scientifique CSV/JSON seulement s’il apporte une utilité réelle.

---

# 4. Audits de qualité restants

- [ ] auditer les options TypeScript strictes finales ;
- [ ] effectuer un audit final TypeScript/ESLint ciblant notamment promesses oubliées, conversions douteuses et `any` non justifiés ;
- [ ] étendre la CI aux futures régressions visuelles/scientifiques critiques.

---

# 5. Évolutions facultatives hors blocage `v5.0.0`

Ces éléments restent des pistes et ne doivent pas retarder la release s’ils ne sont pas retenus.

- [ ] représenter un éventuel courant de probabilité comme **champ de courant**, jamais comme trajectoires individuelles ;
- [ ] ajouter plusieurs niveaux d’isosurface uniquement après validation visuelle sur cas connus ;
- [ ] étendre les familles d’orbitales réelles si l’UI en a réellement besoin ;
- [ ] toute future expérience de spectroscopie ou de transitions photon–hydrogène doit repartir d’un modèle neuf, documenté et testé.

---

# 6. Release `v5.0.0`

- [ ] créer les tags de jalon encore nécessaires ;
- [ ] publier `v5.0.0` uniquement lorsque validation globale, accessibilité, Firefox/WebKit, régressions visuelles, performance, documentation scientifique, licence et CI critique sont qualifiées.
