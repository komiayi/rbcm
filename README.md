# rbcm — Regression-Based methods for Correlated Mediators

> An R Shiny application for causal mediation analysis when mediators are **correlated but not causally linked**.

[![Status](https://img.shields.io/badge/status-active%20development-orange)](https://github.com/komiayi/rbcm)
[![R](https://img.shields.io/badge/R-%E2%89%A5%204.0-blue?logo=r&logoColor=white)](https://www.r-project.org/)
[![Shiny](https://img.shields.io/badge/Shiny-app-brightgreen?logo=rstudioide&logoColor=white)](https://shiny.posit.co/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
![Last commit](https://img.shields.io/github/last-commit/komiayi/rbcm)

<!--
  Once a screenshot of the application interface is available,
  uncomment the line below and place the file at docs/screenshot.png
-->
<!-- ![rbcm — application interface preview](docs/screenshot.png) -->

---

##  Overview

`rbcm` is an R Shiny application implementing two parametric methods for the identification of natural direct and indirect effects in causal mediation analysis when mediators are **correlated but not causally linked** — that is, when their statistical association arises from an unmeasured common cause rather than from a direct causal pathway between them.

The application provides two original estimation strategies:

- **CC (Constant Correlation)** — a method that assumes a stable correlation structure between mediator residuals, simplifying effect identification.
- **CNC (Non-Constant Correlation)** — a generalized parametric method that models residual correlation as a function of the exposure level, relaxing the often unrealistic constant-correlation assumption.

Both methods were developed as part of my Master's thesis in Statistics at the Université du Québec à Montréal (2025), with an applied case study on DNA methylation and cortisol stress reactivity following childhood trauma.

The objective of this application is to **make these methods accessible to applied researchers without requiring advanced R programming skills**, through an interactive interface that integrates data import, assumption checking, estimation, and visualization within a single workflow.

---

##  Related work

The methodological foundations and a complete empirical application are documented in a companion repository:

> **[`komiayi/dna_mediation`](https://github.com/komiayi/dna_mediation)** — research project applying the CC and CNC methods to a real-world dataset linking childhood trauma, DNA methylation (KITLG and JAZF1 loci), and cortisol stress reactivity.

---

## Current features

-  **CC and CNC estimation of mediation effects** — a novel implementation made publicly available.
-  **Automated assumption checks** — normality, independence, and residual correlation structure.
-  **Interactive visualizations** — correlation diagrams, distribution plots, and confidence intervals.
-  **Flexible data import** — CSV, Excel, and delimited file formats.
-  **Contextual documentation** integrated within each module of the interface.

---

## ️ Technology stack

| Component         | Technology                                          |
| ----------------- | --------------------------------------------------- |
| Core language     | R (≥ 4.0)                                           |
| Web framework     | Shiny                                               |
| Front-end         | HTML, custom CSS                                    |
| Visualization     | ggplot2, plotly                                     |
| Optimization      | `optim` (migration to Rcpp in progress)             |
| Version control   | Git / GitHub                                        |

---

##  Active development

Optimization of the CC and CNC estimators is currently the focus of ongoing technical work:

- **Identified bottleneck.** R's `optim()` routine consumes substantial memory in a Shiny runtime context when applied to high-dimensional matrices (typically beyond 100 000 parameters).
- **Approaches under evaluation:**
  - Migrating the optimization layer to `nlminb()` or `optimx::optimx()` for improved memory efficiency.
  - Rewriting the objective function in C++ via `Rcpp` for substantial performance gains.
  - Adopting an asynchronous architecture with `future` and `promises` to offload heavy computations from the Shiny server.
  - Implementing block-wise (chunked) processing for very large matrices.

See [`ROADMAP.md`](ROADMAP.md) for the full technical roadmap.

---

## 🚀 Installation

**1. Clone the repository** (from a terminal):

```bash
git clone https://github.com/komiayi/rbcm.git
cd rbcm
```

**2. Install the R dependencies** (within R or RStudio):

```r
install.packages(c("shiny", "ggplot2", "plotly", "readxl", "DT"))
```

**3. Launch the application** (within R or RStudio):

```r
shiny::runApp()
```

---

##  References

The methods implemented in this application are derived from my Master's thesis:

> Ayi, K. R. (2025). *Analyse de médiation causale pour des médiateurs non causalement liés* [Causal mediation analysis for non-causally-linked mediators] [Master's thesis, Université du Québec à Montréal]. Archipel UQAM. https://archipel.uqam.ca/19950

📄 [Download the full PDF](http://archipel.uqam.ca/19950/1/M19270.pdf)

** Supervision:** Prof. Karim Oualkacha (Department of Mathematics, UQAM) and Prof. Geneviève Lefebvre (Department of Mathematics, UQAM).

### Associated presentations

- ESPUM Symposium — Quantitative Methods in Health Research (2025)
- Mediation Research Days, UQAM (2024)
- SSC Annual Meeting, Carleton University, Ottawa (2023)

---

##  How to cite

If you use this software in academic or applied work, please cite both the software and the underlying thesis. A [`CITATION.cff`](CITATION.cff) file is provided at the root of the repository, enabling GitHub's *Cite this repository* feature (visible in the **About** panel) to generate citations automatically in APA, BibTeX, and other standard formats.

---

##  Contributing

This project is under active development. Feedback, methodological suggestions, and bug reports are welcome through the [Issues](https://github.com/komiayi/rbcm/issues) tab. For substantial contributions, please open an issue first to discuss the proposed changes.

---

##  License

Distributed under the MIT License. See [`LICENSE`](LICENSE) for the complete terms.

---

##  Author

**Komi Roger Ayi**
Biostatistician — Health Data Analyst
Montréal, Québec, Canada

[Portfolio](https://komiayi.github.io) · [LinkedIn](https://www.linkedin.com/in/komi-ayi) · [GitHub](https://github.com/komiayi)
