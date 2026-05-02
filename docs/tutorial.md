# 📘 Tutoriel — Utiliser rbcm

Ce guide explique comment utiliser l'application R Shiny `rbcm` pour analyser un cas de médiation causale avec deux médiateurs corrélés.

---

## 🚀 Lancement de l'application

### Option 1 — Depuis RStudio

```r
# Cloner le dépôt et ouvrir le projet dans RStudio
# Puis exécuter dans la console :
shiny::runApp()
```

### Option 2 — Depuis R en ligne de commande

```bash
cd rbcm
Rscript -e "shiny::runApp()"
```

L'application s'ouvre automatiquement dans votre navigateur à l'adresse `http://127.0.0.1:XXXX`.

---

## 📂 Préparation des données

Votre fichier de données doit être au format **CSV** ou **RData** et contenir au minimum les variables suivantes :

| Variable | Type | Description |
|----------|------|-------------|
| Exposition (X) | Numérique (idéalement binaire 0/1) | Variable causale d'intérêt |
| Médiateurs (M1, M2) | Numérique continu | Variables intermédiaires sur le chemin causal |
| Réponse (Y) | Numérique continu | Variable d'issue |
| Covariables | Numérique ou catégoriel | Variables de contrôle (optionnel) |

**Un fichier exemple est disponible** dans `inst/extdata/exemple_methylation.csv` : il simule un scénario d'épigénétique avec 500 sujets, deux médiateurs (intensités de méthylation `m1` et `m2`), une exposition binaire `x`, et trois covariables (`age`, `sex`, `bmi`).

---

## 🎯 Étape par étape

### Étape 1 — Importation et configuration

Dans l'onglet **Analyse Descriptive → Vue d'ensemble** :

1. **Téléchargez votre fichier** via le bouton « Fichier de données »
2. **Renseignez les noms exacts** de vos variables dans les champs prévus :
   - Exposition (X)
   - Réponse (Y)
   - Médiateurs (séparés par des virgules : `m1,m2`)
   - Médiateur primaire (ex : `m1`)
   - Covariables des médiateurs et de la réponse
3. **Choisissez les options** :
   - Nombre de réplications bootstrap (par défaut 20, recommandé : 200+ pour une vraie analyse)
   - Inclure des interactions ? (TRUE/FALSE)
4. **Cliquez sur « Charger les données »**

Avec le fichier exemple, utilisez ces paramètres :
- Exposition : `x`
- Réponse : `y`
- Médiateurs : `m1,m2`
- Médiateur primaire : `m1`
- Covariables médiateur primaire : `age`
- Covariables médiateur secondaire : `sex`
- Covariables réponse : `bmi`

### Étape 2 — Vérification des données

Onglet **Analyse Descriptive → Résumé statistique** :
- Vérifiez l'aperçu des données
- Examinez les statistiques descriptives quantitatives et qualitatives
- Inspectez les résidus des modèles ajustés

### Étape 3 — Vérification des hypothèses

Onglet **Analyse Descriptive → Graphiques** :
- **Boxplots** : repérez les valeurs aberrantes
- **Histogrammes** : vérifiez la forme des distributions
- **QQ-plots** : évaluez visuellement la normalité
- **Tests de normalité** :
  - Univariés (Shapiro-Wilk, Kolmogorov-Smirnov, Jarque-Bera selon la taille)
  - Multivarié (Henze-Zirkler)

⚠️ **Important** : les méthodes CC et CNC supposent une normalité multivariée des résidus des médiateurs. Si l'hypothèse est rejetée, envisagez une transformation des données.

### Étape 4 — Estimation des effets de médiation

Onglet **Médiation** :

Le tableau affiche les estimations pour les méthodes :
- **CC** : corrélation constante entre les médiateurs
- **CNCm** : corrélation non-constante (moyenne des deux groupes)
- **CNCr** : corrélation non-constante (estimée par les racines)

Pour chaque effet :
- **\(\zeta\)** : effet direct
- **\(\delta\)** : effet indirect

Les écarts-types et intervalles de confiance à 95 % sont obtenus par bootstrap.

---

## 📊 Interprétation des résultats

| Effet | Notation | Interprétation |
|-------|----------|----------------|
| Effet direct (DE) | \(\zeta\) | Effet de X sur Y *non* médié par M1 et M2 |
| Effet indirect (IE) | \(\delta\) | Effet de X sur Y *à travers* M1 et M2 |
| Effet total | DE + IE | Effet causal total de X sur Y |

**Avec les données d'exemple**, les valeurs vraies (utilisées pour la simulation) sont :
- Effet direct ≈ 0.40
- Effet indirect total ≈ 0.45 (= 0.5 × 0.6 + 0.3 × 0.5)

Vos estimations devraient s'approcher de ces valeurs, avec des intervalles de confiance qui les contiennent.

---

## ⚠️ Limitations actuelles

- Le calcul des effets via `optim()` peut être lent ou échouer pour de très grandes matrices (> 100 000 paramètres). Voir [ROADMAP.md](../ROADMAP.md) pour le plan d'optimisation.
- L'application gère actuellement deux médiateurs uniquement.
- Les médiateurs doivent être continus.

---

## 🆘 Dépannage

| Problème | Solution |
|----------|----------|
| « Aucun fichier chargé » | Vérifiez l'extension (.csv ou .RData) |
| Variables non trouvées | Assurez-vous que les noms saisis correspondent **exactement** aux noms de colonnes |
| Erreur lors de l'estimation | Vérifiez que les médiateurs sont continus et que la taille d'échantillon est suffisante (n ≥ 50) |
| Application qui se bloque | Réduisez le nombre de réplications bootstrap |

---

## 📚 Pour aller plus loin

- Consultez l'onglet **Méthodes** dans l'application pour les détails mathématiques
- Lisez le mémoire de référence (à venir)
- Ouvrez une [issue GitHub](https://github.com/komiayi/rbcm/issues) pour toute question

---

*Dernière mise à jour : avril 2026*
