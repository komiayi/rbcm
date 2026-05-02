library(testthat)

# Charger toutes les fonctions du dossier R/
for (f in list.files("../../R", pattern = "\\.R$", full.names = TRUE)) {
  source(f, encoding = "UTF-8")
}

test_check("rbcm")
