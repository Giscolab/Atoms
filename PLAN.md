# Atoms — Plan restant

Ce document contient uniquement le travail encore ouvert. Toute tâche implémentée et validée doit en disparaître.

**Projet :** `Giscolab/Atoms`
**Cible :** Showcase OpenAI puis évolutions non bloquantes.

## 1. Showcase

- [ ] soumettre Atoms au Showcase OpenAI.

## 2. Import / export scientifique — évolution non bloquante

- [ ] définir un schéma JSON versionné avec unités, base, état quantique, paramètres de rendu et seed explicites ;
- [ ] valider strictement les imports et refuser les unités ou états ambigus ;
- [ ] raccorder l’import à l’interface uniquement lorsque le schéma est stabilisé ;
- [ ] permettre l’export d’un état JSON reproductible incluant les paramètres et la seed ;
- [ ] permettre une capture PNG depuis l’application ;
- [ ] ajouter un export scientifique CSV/JSON seulement s’il apporte une utilité réelle.

## 3. Évolutions facultatives

Ces pistes ne bloquent pas le Showcase et peuvent être traitées après `v5.0.0`.

- [ ] représenter un éventuel courant de probabilité comme **champ de courant**, jamais comme trajectoires individuelles ;
- [ ] ajouter plusieurs niveaux d’isosurface uniquement après validation visuelle sur des cas connus ;
- [ ] étendre les familles d’orbitales réelles seulement si l’interface en a réellement besoin ;
- [ ] envisager le code-splitting du renderer/Three.js uniquement si une mesure de chargement réel montre un bénéfice justifiant la complexité supplémentaire ;
- [ ] toute future expérience de spectroscopie ou de transitions photon–hydrogène doit repartir d’un modèle scientifique autonome, documenté et testé.
