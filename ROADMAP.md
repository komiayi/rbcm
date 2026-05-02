# 🗺️ ROADMAP — rbcm

Ce document décrit l'état actuel du projet et les étapes de développement prévues.

---

## ✅ Réalisé (v0.1 — actuelle)

### Fondations de l'application
- [x] Architecture Shiny modulaire (UI / serveur / utilitaires)
- [x] Interface CSS personnalisée
- [x] Système de navigation entre modules

### Importation et préparation des données
- [x] Importation de fichiers CSV
- [x] Importation de fichiers Excel (xlsx, xls)
- [x] Visualisation tabulaire des données importées (DT)
- [x] Sélection interactive des variables d'intérêt (exposition, médiateurs, issue)

### Vérification des hypothèses
- [x] Tests de normalité (Shapiro-Wilk, Kolmogorov-Smirnov)
- [x] Visualisation des distributions (histogrammes, QQ-plots)
- [x] Calcul et affichage des matrices de corrélation
- [x] Détection visuelle des structures de corrélation

### Visualisation
- [x] Graphiques dynamiques avec ggplot2 / plotly
- [x] Modules interactifs de visualisation des corrélations
- [x] Export des graphiques (PNG, PDF)

---

## 🚧 En cours (v0.2)

### Optimisation des estimateurs CNC et CC

**Objectif :** rendre les méthodes utilisables sur des matrices à haute dimension (>100 000 paramètres) dans un contexte Shiny.

**Défi identifié :** la fonction `optim()` de base R consomme trop de mémoire dans l'environnement Shiny lors de l'estimation des paramètres pour de grandes matrices.

**Pistes en cours d'évaluation :**

- [ ] **Benchmarking des optimiseurs alternatifs**
  - [ ] Comparaison `optim()` vs `nlminb()` vs `optimx::optimx()`
  - [ ] Évaluation de `nloptr` (algorithmes non-linéaires modernes)
  - [ ] Mesure du compromis temps / mémoire / précision

- [ ] **Réécriture en C++ via Rcpp**
  - [ ] Profiling pour identifier les goulots d'étranglement
  - [ ] Réécriture de la fonction objectif en C++
  - [ ] Tests de régression pour garantir l'équivalence numérique

- [ ] **Architecture asynchrone**
  - [ ] Intégration de `future` + `promises`
  - [ ] Calculs longs déplacés hors du serveur Shiny principal
  - [ ] Indicateurs de progression utilisateur

- [ ] **Stratégie de chunking**
  - [ ] Traitement par blocs pour matrices très larges
  - [ ] Gestion de la mémoire par lots

---

## 📋 Prévu (v0.3)

### Méthodes statistiques complètes

- [ ] Implémentation finale de la méthode **CC (Constant Correlation)**
- [ ] Implémentation finale de la méthode **CNC (Non-Constant Correlation)**
- [ ] Calcul des intervalles de confiance (bootstrap)
- [ ] Tests d'hypothèses sur les effets de médiation
- [ ] Décomposition des effets directs et indirects

### Modules de visualisation des résultats

- [ ] Diagrammes de médiation (path diagrams)
- [ ] Visualisation comparative CC vs CNC
- [ ] Graphiques d'évolution des corrélations entre médiateurs
- [ ] Heatmaps interactives

### Documentation utilisateur

- [ ] Tutoriel interactif intégré à l'application
- [ ] Vidéo de démonstration
- [ ] Cas d'étude reproductible (données simulées)

---

## 🎯 Objectif v1.0

### Application stable et publiable

- [ ] Tests automatisés (`testthat`)
- [ ] Intégration continue (GitHub Actions)
- [ ] Hébergement sur shinyapps.io ou Posit Connect
- [ ] Article méthodologique soumis

### Package R complémentaire

- [ ] Création d'un package R `rbcm` avec les fonctions principales
- [ ] Soumission au CRAN
- [ ] Vignette et documentation Roxygen2

---

## 💡 Ambitions futures (v2.0+)

- [ ] Extension à plus de 2 médiateurs corrélés
- [ ] Support des médiateurs binaires et catégoriels
- [ ] Intégration avec d'autres formats de données (-omiques)
- [ ] API REST pour usage hors Shiny
- [ ] Version multilingue (FR / EN)

---

## 📊 Indicateurs de progression

| Version | Statut | Date cible |
|---------|--------|------------|
| v0.1 | ✅ Disponible | Avril 2026 |
| v0.2 | 🚧 En développement | Été 2026 |
| v0.3 | 📋 Planifiée | Automne 2026 |
| v1.0 | 🎯 Cible | Fin 2026 |

---

## 🤝 Comment contribuer

Les retours sur l'optimisation sont particulièrement bienvenus à ce stade. N'hésitez pas à ouvrir une [issue](../../issues) pour :

- Suggérer une approche d'optimisation
- Signaler un bug
- Proposer une fonctionnalité
- Partager un cas d'usage

---

*Roadmap mise à jour : avril 2026*
