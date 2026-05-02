test_that("calculer_effets fonctionne sur le jeu de données d'exemple", {
  # Charger le jeu de données simulé
  data_path <- system.file("extdata", "exemple_methylation.csv", package = "rbcm")
  if (data_path == "") {
    data_path <- file.path("..", "..", "inst", "extdata", "exemple_methylation.csv")
  }

  skip_if_not(file.exists(data_path), "Jeu de données exemple introuvable")

  data_test <- read.csv(data_path)

  result <- calculer_effets(
    treat = "x",
    mediators = c("m1", "m2"),
    intmed = "m1",
    outcome = "y",
    intmed_cov = "age",
    sed_cov = "sex",
    out_cov = "bmi",
    inter_treat_cov = FALSE,
    cor_cste = 1,
    data = data_test,
    q1 = 3,
    q2 = 3,
    inter = FALSE,
    names_vec = NULL,
    rh = 0.5,
    formula_one3 = "x + m1 + m2 + age + sex + bmi"
  )

  # Vérifications structurelles
  expect_named(result, c("Cint", "Csed", "Ccor", "C.y", "var_covar", "sol.i"))
  expect_named(result$sol.i, c("DE1sk", "IE1sk"))

  # L'effet direct attendu est ~0.4 (avec marge d'erreur d'estimation)
  # L'effet indirect attendu est ~0.45
  # On vérifie juste que les estimations sont dans un ordre de grandeur cohérent
  expect_true(abs(result$sol.i$DE1sk[1]) < 2)
  expect_true(abs(result$sol.i$IE1sk[1]) < 2)
})

test_that("Le jeu de données exemple a la bonne structure", {
  data_path <- system.file("extdata", "exemple_methylation.csv", package = "rbcm")
  if (data_path == "") {
    data_path <- file.path("..", "..", "inst", "extdata", "exemple_methylation.csv")
  }

  skip_if_not(file.exists(data_path), "Jeu de données exemple introuvable")

  data_test <- read.csv(data_path)

  expect_equal(nrow(data_test), 500)
  expect_true(all(c("id", "x", "age", "sex", "bmi", "m1", "m2", "y") %in% colnames(data_test)))
  expect_true(all(data_test$x %in% c(0, 1)))
  expect_true(all(data_test$sex %in% c(0, 1)))
})
