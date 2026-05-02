
# Normality test

test_normalite_clean <- function(vos_donnees) {
  n <- length(vos_donnees)
  
  if (n <= 50) {
    res <- shapiro.test(vos_donnees)
    method <- "Shapiro-Wilk"
  } else if (n <= 2000) {
    res <- ks.test(vos_donnees, "pnorm", mean = mean(vos_donnees), sd = sd(vos_donnees))
    method <- "Kolmogorov-Smirnov"
  } else {
    res <- tseries::jarque.bera.test(vos_donnees)
    method <- "Jarque-Bera"
  }
  
  # On crée un petit tableau propre
  data.frame(
    "Statistical test" = method,
    "Statistic" = round(res$statistic, 4),
    "p-value" = format.pval(res$p.value, digits = 4),
    "Interpretation" = ifelse(res$p.value > 0.05, "Normal Distribution", "Non-normal Distribution"),
    check.names = FALSE
  )
}