# rbcm — Regression-Based methods for Correlated Mediators

> **Application R Shiny pour l'analyse de médiation causale en présence de médiateurs corrélés.**

[![Status](https://img.shields.io/badge/status-active%20development-orange)]()
[![R](https://img.shields.io/badge/R-≥4.0-blue)]()
[![Shiny](https://img.shields.io/badge/Shiny-app-brightgreen)]()

---
![Aperçu de l'application rbcm](docs/screenshot.png)
## 📖 Description

`rbcm` est une application R Shiny implémentant deux méthodes paramétriques originales pour l'identification des effets de médiation causale lorsque les médiateurs sont corrélés :

- **CC (Constant Correlation)** — méthode adaptée lorsque la structure de corrélation entre médiateurs est stable
- **CNC (Non-Constant Correlation)** — méthode généralisée pour les structures de corrélation hétérogènes

Ces méthodes ont été développées dans le cadre de mon mémoire de maîtrise en statistique à l'UQAM (2026), avec une application aux données de méthylation de l'ADN.

L'objectif de cette application est de **rendre ces méthodes accessibles aux chercheurs non-statisticiens** via une interface interactive, sans nécessiter de programmation R avancée.

---
## 🔗 Travaux liés

Les méthodes CC et CNC implémentées dans cette application ont été développées 
et validées dans le cadre du projet de recherche [`dna_mediation`](https://github.com/komiayi/dna_mediation), 
qui présente l'application des méthodes à un jeu de données réelles 
(traumatismes infantiles, méthylation de l'ADN, réactivité au cortisol).

## ✨ Fonctionnalités actuelles

* 🎯 **Estimation des effets de médiation par les méthodes CC et CNC** — 
  uniquement implémentées dans cette application, à ma connaissance
* 🔍 **Vérification automatisée des hypothèses** : normalité, indépendance, 
  structure de corrélation des résidus
* 📊 **Visualisations interactives** : diagrammes de corrélation, 
  distributions, intervalles de confiance
* 📥 **Importation flexible** : CSV, Excel, fichiers délimités
* 📑 **Documentation contextuelle** intégrée dans chaque module
---

## 🛠️ Technologies utilisées

| Composant | Technologie |
|-----------|-------------|
| Langage principal | R (≥ 4.0) |
| Framework web | Shiny |
| Interface | HTML, CSS personnalisé |
| Visualisation | ggplot2, plotly |
| Calcul statistique | optim (en cours de migration vers Rcpp) |
| Versionnage | Git / GitHub |

---

## 🚧 En développement actif

L'optimisation de l'estimateur des paramètres CNC/CC fait actuellement l'objet d'un travail technique :

- **Défi identifié** : la fonction d'optimisation `optim()` de R consomme une mémoire significative dans le contexte Shiny pour des matrices de haute dimension (typiquement > 100 000 paramètres).
- **Approches évaluées** :
  - Migration vers `nlminb()` ou `optimx::optimx()` pour une meilleure efficacité mémoire
  - Réécriture de la fonction objectif en C++ via `Rcpp` pour des gains de performance significatifs
  - Architecture asynchrone avec `future` + `promises` pour décharger le serveur Shiny
  - Traitement par blocs (chunking) pour les matrices de très grande dimension

Voir [ROADMAP.md](ROADMAP.md) pour le détail des prochaines étapes.

---

## 🚀 Installation

**1. Cloner le dépôt** (dans un terminal) :

```bash
git clone https://github.com/komiayi/rbcm.git
cd rbcm
```

**2. Installer les dépendances** (dans R / RStudio) :

```r
install.packages(c("shiny", "ggplot2", "plotly", "readxl", "DT"))
```

**3. Lancer l'application** (dans R / RStudio) :

```r
shiny::runApp()
```

---

## 📚 Références scientifiques

Les méthodes implémentées sont issues de mon mémoire de maîtrise :

> Ayi, K. R. (2025). *Analyse de médiation causale pour des médiateurs non 
> causalement liés* [Mémoire de maîtrise, Université du Québec à Montréal]. 
> Archipel UQAM. https://archipel.uqam.ca/19950

📄 [Télécharger le PDF complet](http://archipel.uqam.ca/19950/1/M19270.pdf)

**Directeur de recherche :** Pr Karim Oualkacha, Département de mathématiques, UQAM.

### Présentations associées
- Colloque ESPUM — Méthodes Quantitatives en Santé (2025)
- Mediation Research Days, UQAM (2024)
- Congrès SSC, Université Carleton, Ottawa (2023)

---

## 🤝 Contributions

Le projet est en développement actif. Les retours, suggestions et issues sont les bienvenus via l'onglet [Issues](../../issues).

---

## Licence

Distribué sous licence MIT. Voir [`LICENSE`](LICENSE) pour les détails complets.

---

## 👤 Auteur

**Komi Roger Ayi**  
Biostatisticien — Analyste de données en santé  
Montréal, Québec  
[LinkedIn](https://www.linkedin.com/in/komi-ayi) • [GitHub](https://github.com/komiayi) • [Portfolio](https://komiayi.github.io)

---
![Last commit](https://img.shields.io/github/last-commit/komiayi/rbcm)
