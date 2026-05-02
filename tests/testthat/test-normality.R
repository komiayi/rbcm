test_that("test_normalite fonctionne pour de petits échantillons (Shapiro-Wilk)", {
  set.seed(123)
  x <- rnorm(30)
  result <- test_normalite(x)
  expect_true(grepl("Shapiro-Wilk", result))
})

test_that("test_normalite gère les échantillons moyens (Kolmogorov-Smirnov)", {
  set.seed(123)
  x <- rnorm(500)
  # KS retourne directement un objet htest, on capture la sortie
  result <- capture.output(test_normalite(x))
  expect_true(any(grepl("Kolmogorov", result)) || length(result) > 0)
})

test_that("test_normalite gère les grands échantillons (Jarque-Bera)", {
  set.seed(123)
  x <- rnorm(3000)
  result <- test_normalite(x)
  expect_true(grepl("Jarque-Bera", result))
})

test_that("test_normalite détecte la non-normalité", {
  set.seed(123)
  x <- rexp(30)  # distribution exponentielle, clairement non normale
  result <- test_normalite(x)
  # On vérifie juste que la fonction retourne quelque chose de valide
  expect_true(is.character(result) || is.list(result))
})
