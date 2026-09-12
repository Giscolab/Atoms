# Atoms — Plan restant

Ce document contient uniquement le travail encore ouvert. Toute tâche implémentée et validée doit en disparaître.

**Projet :** `Giscolab/Atoms`
**Cible :** Showcase puis release `v5.0.0`

## 1. Showcase

- [ ] finaliser le dossier de soumission : couverture, titre, tagline, description, cas d’usage, stack, processus de construction avec Codex/GPT et URLs publiques ;
- [ ] soumettre Atoms au Showcase OpenAI.

## 2. Statut de licence

- [ ] décider explicitement du statut de licence du dépôt : publier une licence choisie par le propriétaire ou conserver volontairement l’absence de licence.

Cette décision appartient au propriétaire du projet et ne doit pas être prise automatiquement par un agent.

## 3. Import / export scientifique — évolution non bloquante

- [ ] définir un schéma JSON versionné avec unités, base, état quantique, paramètres de rendu et seed explicites ;
- [ ] valider strictement les imports et refuser les unités ou états ambigus ;
- [ ] raccorder l’import à l’interface uniquement lorsque le schéma est stabilisé ;
- [ ] permettre l’export d’un état JSON reproductible incluant les paramètres et la seed ;
- [ ] permettre une capture PNG depuis l’application ;
- [ ] ajouter un export scientifique CSV/JSON seulement s’il apporte une utilité réelle.

## 4. Évolutions facultatives

Ces pistes ne bloquent ni le Showcase ni la release si elles ne sont pas retenues.

- [ ] représenter un éventuel courant de probabilité comme **champ de courant**, jamais comme trajectoires individuelles ;
- [ ] ajouter plusieurs niveaux d’isosurface uniquement après validation visuelle sur des cas connus ;
- [ ] étendre les familles d’orbitales réelles seulement si l’interface en a réellement besoin ;
- [ ] toute future expérience de spectroscopie ou de transitions photon–hydrogène doit repartir d’un modèle scientifique autonome, documenté et testé.

## 5. Release `v5.0.0`

- [ ] créer le tag `v5.0.0` et publier la release lorsque le statut de licence est décidé et que le dernier commit destiné à la release est vert en CI et sur GitHub Pages.
