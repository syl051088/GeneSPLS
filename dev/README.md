# GeneSPLS: Sparse Partial Least Squares for Genomic Data Analysis

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

## Description

<div align="justify">

**GeneSPLS** is a comprehensive toolkit designed for performing Sparse Partial Least Squares (SPLS) regression analyses on genomic data. The project primarily focuses on eQTL (expression quantitative trait loci) mapping where the goal is to identify associations between genetic variants (SNPs) and gene expression levels.

In this project, we utilize source data from the **Geuvadis consortium**, which provides an extensive collection of RNA-seq gene expression data and corresponding high-density genotype information from European populations. The dataset encompasses expression levels for thousands of genes across multiple individuals, allowing for the exploration of genetic regulation mechanisms. Our data preprocessing scripts harmonize and integrate these datasets to create both univariate and multivariate outcome matrices for analysis.

This repository includes:

- **Simulated Data Analyses:** Two R Markdown files (`test_simulated_univariate.Rmd` and `test_simulated_multivariate.Rmd`) that demonstrate SPLS regression using simulated datasets. These scripts showcase both univariate and multivariate outcome scenarios.
- **Real Data Analyses:** Two R Markdown files (`test_real_univariate.Rmd` and `test_real_multivariate.Rmd`) that apply SPLS regression to real genomic data. For the real data, gene expression (including AKT3 expression) and genotype data (from the Geuvadis consortium) are pre-processed, tuned, and analyzed.
- **Efficient Implementation:** A collection of C++ source files (via Rcpp) that implement the core algorithms, including cross-validation (`cv_spls_cpp`) and the main SPLS function (`spls_cpp`). The algorithms are inspired by the kernel-based PLS methodology (see Rännär et al., 1994) and further enhanced to incorporate sparsity following ideas from Chun and Keles (2010).

### Primary Purpose

The central research question addressed by this project is: **Which genomic variants (SNPs) are significantly associated with gene expression differences, and how can we effectively identify these associations using a sparse dimension reduction approach?** By employing SPLS regression, we aim to reduce the high dimensionality inherent in genomic data while simultaneously selecting the most relevant variables, thereby enhancing the interpretability and predictive performance in eQTL analyses. Special attention is given to key genes such as AKT3, providing insights into their regulatory mechanisms.

### Mathematical Formulation

The underlying PLS regression model is defined as:
  
$$
{\huge
Y = X B_{PLS} + F,
}
$$

where:  
- $X$ is the $n \times p$ predictor matrix (SNP genotypes),  
- $Y$ is the $n \times q$ response matrix (gene expression),  
- $B_{PLS}$ is the matrix of regression coefficients,  
- $F$ represents the residual errors.

For the sparse formulation, we solve an optimization problem that imposes an \(L_1\) penalty to encourage sparsity in the direction vectors. A generalized formulation is:

$$
{\huge
\min_{w,c} \left\lbrace -\kappa\, w^T M w + \frac{1-\kappa}{2}\|c-w\|^2 + \lambda_1 \|c\|_1 + \lambda_2 \|c\|_2^2 \right\rbrace \quad \text{subject to } \|w\|_2 = 1,
}
$$

where $M = X^T Y Y^T X$, and $\\kappa$, $\\lambda_1$, and $\\lambda_2$ are tuning parameters. Cross-validation (via `cv_spls_cpp`) is used to determine the optimal values for the sparsity (typically through the candidate eta values) and the number of latent components $K$.

This formulation combines dimension reduction with variable selection, allowing the model to focus on the most relevant genetic variants driving gene expression differences.

## Dependencies & Prerequisites

- **R**: Tested on R version 4.3 (or a recent stable release).
- **Rcpp & RcppArmadillo**: For C++ integration and efficient matrix operations.
- **R Packages**:
  - `microbenchmark` – for performance benchmarking.
  - `data.table` (optional) – for fast data reading when dealing with large files.
  - Additional dependencies may be required as specified in the individual scripts.

Ensure that you have the GNU General Public License v3.0 for licensing.

## Walkthrough & Usage

### Simulated Data Analysis
- **test_simulated_univariate.Rmd**: Demonstrates SPLS regression on simulated data with a single continuous outcome. It covers data simulation, parameter tuning using cross-validation, model fitting using both the R and Rcpp implementations, and benchmarking.
- **test_simulated_multivariate.Rmd**: Extends the analysis to simulated gene expression data with multivariate outcomes. It follows similar steps as the univariate analysis but emphasizes the multivariate model's consistency and efficiency.

### Real Data Analysis
- **test_real_univariate.Rmd**: Applies SPLS regression to real eQTL data, specifically focusing on the AKT3 gene expression levels. The script includes data loading, preprocessing (filtering SNPs by standard deviation), parameter tuning, model fitting, and result comparison.
- **test_real_multivariate.Rmd**: Analyzes real gene expression data with multivariate outcomes. It demonstrates how to select a subset of genes as outcomes and compares SPLS implementations for performance.

### Data Preprocessing
A separate script (`Process_GeneExpression_Genotype.R`) is provided to harmonize and process the genotype and gene expression datasets (sourced from the Geuvadis consortium). This script:
- Reads and preprocesses large-scale RNA-seq data and genotype information.
- Extracts specific gene expression levels (e.g., AKT3).
- Generates a covariate matrix from genotype data.
- Produces output files for downstream analysis (e.g., `akt3_expression.csv`, `genotype_matrix.rds`, and `gene_expression_matrix.rds`).

## Project Structure

```
GeneSPLS
├── DESCRIPTION
├── dev
│   ├── test_simulated_univariate.Rmd
│   ├── widekernelpls_cv_fit.cpp
│   ├── X_genotype_matrix.rda
│   └── Y_AKT3_expression.rda
├── GeneSPLS.Rproj
├── LICENSE
├── man
│   ├── correctp.Rd
│   ├── spls_cpp.Rd
│   ├── spls_dv.Rd
│   ├── ust.Rd
│   └── widekernelpls_fit.Rd
├── NAMESPACE
├── R
│   └── RcppExports.R
├── README.md
├── src
│   ├── correctp.cpp
│   ├── correctp.h
│   ├── cv.spls.cpp
│   ├── cv.spls.h
│   ├── GeneSPLS.dll
│   ├── GeneSPLS.so
│   ├── Makevars
│   ├── Makevars.win
│   ├── pls.cpp
│   ├── pls.h
│   ├── RcppExports.cpp
│   ├── spls.cpp
│   ├── spls_dv.cpp
│   ├── spls_dv.h
│   ├── spls.h
│   └── spls.o
├── tests
│   ├── testthat
│   │   ├── test.Rmd
│   │   └── test-widekernelpls.R
│   └── testthat.R
└── vignettes
    └── Tutorial.Rmd
```

## Contributing

Contributions are welcome! Please follow these guidelines when contributing:

- Fork the repository and create a feature branch.
- Ensure your code follows the established style and includes tests.
- Submit a pull request with a clear description of your changes.
- For major changes, please open an issue first to discuss the proposed modifications.

## License

This project is licensed under the [GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html).

## Contact

For questions, bug reports, or suggestions, please reach out to the project maintainer:

- **Developer/Maintainer:** StatGPT
- **GitHub:** [syl051088/GeneSPLS](https://github.com/syl051088/GeneSPLS)

## References

1. Rännär, S., Lindgren, F., Geladi, P., & Wold, S. (1994). *A PLS Kernel Algorithm for Data Sets with Many Variables and Fewer Objects: Part 1 – Theory and Algorithm*. Journal of Chemometrics.
2. Chun, H., & Keleş, S. (2010). *Sparse Partial Least Squares Regression for Simultaneous Dimension Reduction and Variable Selection*. Journal of the Royal Statistical Society: Series B.
3. Additional references as cited in the code documentation and in the accompanying literature.
</div>
---

*Current Year: 2025  
Current Date: 2025-04-03*
