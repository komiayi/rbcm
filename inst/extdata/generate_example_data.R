# =============================================================================
# Génération d'un jeu de données simulé pour rbcm
#
# Simule un scénario de médiation causale avec deux médiateurs corrélés,
# inspiré d'un contexte d'épigénétique (méthylation de l'ADN simplifiée).
#
# Variables générées :
#   - id          : identifiant du sujet
#   - x           : exposition binaire (0 = non exposé, 1 = exposé)
#   - age         : covariable continue (années)
#   - sex         : covariable binaire (0 = M, 1 = F)
#   - bmi         : covariable continue (indice de masse corporelle)
#   - m1          : médiateur primaire continu (intensité méthylation locus 1)
#   - m2          : médiateur secondaire continu (intensité méthylation locus 2)
#   - y           : réponse continue (biomarqueur de santé)
#
# Auteur : Komi Roger Ayi
# =============================================================================

set.seed(2026)

n <- 500

# ---- Covariables ----
age <- round(rnorm(n, mean = 50, sd = 10), 1)
sex <- rbinom(n, 1, 0.5)
bmi <- round(rnorm(n, mean = 25, sd = 4), 1)

# ---- Exposition ----
# Probabilité d'exposition dépendant légèrement de l'âge
prob_exp <- plogis(-1 + 0.02 * (age - 50))
x <- rbinom(n, 1, prob_exp)

# ---- Médiateurs corrélés ----
# Effet causal de x sur m1 et m2, avec corrélation résiduelle
library(mvtnorm)

# Matrice de variance-covariance des erreurs des médiateurs
rho_med <- 0.4  # corrélation entre les médiateurs
Sigma_med <- matrix(c(1, rho_med, rho_med, 1), 2, 2)

# Effets de l'exposition sur les médiateurs
alpha1 <- 0.5   # x -> m1
alpha2 <- 0.3   # x -> m2

# Génération
errors_med <- rmvnorm(n, mean = c(0, 0), sigma = Sigma_med)
m1 <- 0.2 + alpha1 * x + 0.01 * (age - 50) + errors_med[, 1]
m2 <- 0.1 + alpha2 * x - 0.05 * sex + errors_med[, 2]

# ---- Réponse Y ----
# Effets directs et indirects
beta_x <- 0.4   # effet direct
beta_m1 <- 0.6  # effet de m1 sur y
beta_m2 <- 0.5  # effet de m2 sur y

y <- 1 + beta_x * x + beta_m1 * m1 + beta_m2 * m2 +
     0.02 * age + 0.05 * bmi + rnorm(n, 0, 1)

# ---- Construction du data frame ----
exemple_methylation <- data.frame(
  id = 1:n,
  x = x,
  age = age,
  sex = sex,
  bmi = bmi,
  m1 = round(m1, 4),
  m2 = round(m2, 4),
  y = round(y, 4)
)

# ---- Sauvegarde ----
write.csv(exemple_methylation, "inst/extdata/exemple_methylation.csv", row.names = FALSE)
save(exemple_methylation, file = "inst/extdata/exemple_methylation.RData")

cat("Jeu de données généré :", n, "observations\n")
cat("Sauvegardé dans inst/extdata/exemple_methylation.csv et .RData\n")
cat("\nValeurs vraies des paramètres pour validation :\n")
cat("  alpha1 (x -> m1) :", alpha1, "\n")
cat("  alpha2 (x -> m2) :", alpha2, "\n")
cat("  beta_x (effet direct) :", beta_x, "\n")
cat("  beta_m1 (m1 -> y) :", beta_m1, "\n")
cat("  beta_m2 (m2 -> y) :", beta_m2, "\n")
cat("  rho (corr. médiateurs) :", rho_med, "\n")
cat("  Effet indirect total attendu :", alpha1 * beta_m1 + alpha2 * beta_m2, "\n")
