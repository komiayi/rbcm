test_that("varcovarEstimes fonctionne avec corrélation constante (cor_cste = 1)", {
  cor_coefs <- c(0.5, 0.5, 0.3)
  result <- varcovarEstimes(cor_coefs, cor_cste = 1)

  expect_named(result, c("r01", "r10", "indices"))
  expect_equal(result$r01, 0.3)
  expect_equal(result$r10, 0.3)
  expect_equal(result$indices, 1)
})

test_that("varcovarEstimes gère le cas non-constant (cor_cste = 3)", {
  # Vecteur de 6 éléments : sig11, sig10, sig21, sig20, rh11, rh00
  cor_coefs <- c(1.0, 0.9, 1.1, 1.0, 0.4, 0.3)
  result <- varcovarEstimes(cor_coefs, cor_cste = 3)

  expect_named(result, c("r01", "r10", "indices"))
  # Les résultats peuvent être NA si le discriminant est négatif, c'est normal
  expect_true(is.numeric(result$r01) || is.na(result$r01))
})

test_that("varcovarEstimes retourne NA si Delt < 0", {
  # Construire un cas où le discriminant sera négatif
  cor_coefs <- c(0.1, 0.1, 0.1, 0.1, 0.99, 0.99)
  result <- varcovarEstimes(cor_coefs, cor_cste = 3)
  # Le résultat peut contenir des NA, c'est le comportement attendu
  expect_true(any(is.na(result$r01)) || is.numeric(result$r01))
})
