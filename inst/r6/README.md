# R6 Refactoring — `rbcm` engine

This branch (`r6-refactor`) contains a **complete refactoring** of the CC and CNC method engine into an **object-oriented architecture** based on the R6 class system, with critical functions optimized in **C++ via Rcpp and RcppArmadillo**.

> ⚠️ **Status**: refactor completed locally, currently being integrated into the main branch and the Shiny application.

---

## 🎯 Goals of the refactoring

The CC/CNC causal mediation analysis from the Master's thesis relies on computationally expensive steps: bound-constrained log-likelihood optimization, computation of direct and indirect effects on extended design matrices, and Hessian-based variance-covariance estimation. The initial procedural implementation had two main limitations:

1. **Performance**: the log-likelihood functions, called hundreds of times by `optim()`, became a bottleneck as data dimensions grew.
2. **Maintainability**: the CC and CNC logic, design matrices, and effect computations were spread across several scripts that were difficult to evolve.

This refactoring addresses both limitations and facilitates integration with the `rbcm` Shiny application.

---

## 🏗️ Architecture

### `RbCorelatedMediators` class (R6)

```r
RbCorelatedMediators <- R6::R6Class(
  "Rbcm",
  public  = list(...),    # User-facing API
  private = list(...)     # Internal helpers (validation, design matrices, etc.)
)
```

#### Public members

- `params`: model parameters (mediators, exposure, response, covariates)
- `estimates_params`: estimated coefficients after optimization
- `data`: input data matrix
- `mediators_model_fit`, `outcome_model_fit`: returned `optim()` objects
- `method`: `"CC"` or `"CNC"`
- `direct`, `indirect`: estimated effect vectors

#### Main public methods

| Method | Role |
| --- | --- |
| `initialize(data, mediators, response, exposition, covariates)` | Validation and data preparation |
| `parmetors_estimators(method, interaction, outc_interactions)` | CC or CNC parameter estimation |
| `effests_Estimate()` | Computation of direct and indirect effects |
| `varcovarH()` | Variance-covariance matrix computation |

#### Private methods (helpers)

Column validation, design matrix construction for mediators and outcome, parameter initialization, `rho01` computation, optimization bounds, coefficient extraction, and numerical utilities (`safe_sqrt`).

### Rcpp/RcppArmadillo optimization

Performance-critical functions are implemented in C++ using **RcppArmadillo** (high-performance C++ linear algebra):

- **`log_lik_med_arma`**: log-likelihood of correlated mediators (CC and CNC), called intensively by `optim()`.
- The C++ code is compiled via `Rcpp::sourceCpp("logmed_rccp.cpp")` and called directly from R.

---

## 📁 Folder structure

```
inst/r6/
├── README.md                # This file
├── RbcmR6.R                 # R6 class + helpers (1150 lines)
└── logmed_rccp.cpp          # C++ log-likelihood function (RcppArmadillo)
```

---

## 🚀 Usage example

```r
library(R6)
library(Rcpp)
library(RcppArmadillo)

source("inst/r6/RbcmR6.R")

# Initialization
model <- RbCorelatedMediators$new(
  data       = my_data,
  mediators  = list(primary = "M1", secondary = "M2"),
  response   = list(outcome = "Y"),
  exposition = list(exposure = "T"),
  covariates = list(cov_primary   = c("X1", "X2"),
                    cov_secondary = c("X3"),
                    cov_outcome   = c("X1", "X2", "X3"))
)

# Parameter estimation (CC or CNC method)
model$parmetors_estimators(method = "CNC")

# Display detailed results
model$estimates_results()

# Compute direct and indirect effects
model$effests_Estimate()

# Variance-covariance matrix
model$varcovarH()
```

---

## 🔧 Work in progress

- [ ] Rename methods to consistent English naming (`parmetors_estimators` → `estimate_parameters`, `effests_Estimate` → `compute_effects`).
- [ ] Migrate error messages to English throughout the class.
- [ ] Add unit tests (`testthat`) for each public method.
- [ ] `roxygen2` documentation for public methods.
- [ ] Integration into the `rbcm` Shiny application (progressive replacement of the procedural backend).
- [ ] Comparative profiling: pure R vs Rcpp implementation (to be published in the main README).

---

## 📚 Methodological references

The CC and CNC methods implemented here come from my Master's thesis:

> Ayi, K. R. (2025). *Analyse de médiation causale pour des médiateurs non causalement liés* [Master's thesis, Université du Québec à Montréal]. Archipel UQAM. https://archipel.uqam.ca/19950

**Supervision**: Prof. Karim Oualkacha (Department of Mathematics, UQAM) and Prof. Geneviève Lefebvre (Biostatistics, UQAM).

---

## 🛠️ Technical stack

| Component | Technology |
| --- | --- |
| Base language | R 4.x |
| OO system | R6 |
| Numerical optimization | `optim()` (BFGS, L-BFGS-B with bounds) |
| High-performance linear algebra | C++ via Rcpp and RcppArmadillo |
| Error handling | `tryCatch`, defensive validations |
