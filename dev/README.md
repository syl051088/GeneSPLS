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

## Geuvadis Consortium: Background and Dataset Overview

The **Geuvadis Consortium** (short for *Genetic EUropean VAriation in DISease*) was an EU-funded research collaboration launched in 2010, aimed at exploring how genetic variation impacts gene expression and human disease ([Sharing capacity across Europe in high-throughput sequencing technology to explore genetic variation in health and disease | GEUVADIS | Project | News & Multimedia | FP7 | CORDIS | European Commission](https://cordis.europa.eu/project/id/261123/reporting/de#:~:text=between%20a%20wide%20range%20of,the%20Consortium%20has%20obtained%20high)).

### Geuvadis Dataset Characteristics

One of the major outcomes of Geuvadis is a **large-scale reference dataset** integrating genomic and transcriptomic information. This dataset is widely used in human genomics as it provides a rich resource for studying **expression quantitative trait loci (eQTLs)** and gene regulation. Key features of the Geuvadis dataset include:

- **Cohort Size & Populations:** The project profiled **≈465 individuals** (after quality control, 462 remained) drawn from **five distinct populations** ([
            Transcriptome and genome sequencing uncovers functional variation in humans - PMC
        ](https://pmc.ncbi.nlm.nih.gov/articles/PMC3918453/#:~:text=We%20performed%20mRNA%20and%20small,%C3%97%2010%E2%88%9216%20for%20mRNA%2C%20p)). These populations were the CEPH Utah residents (**CEU**), Finns (**FIN**), British (**GBR**), Toscani Italians (**TSI**), and Yoruba Nigerians (**YRI**), with roughly 89–95 individuals per group ([
            Transcriptome and genome sequencing uncovers functional variation in humans - PMC
        ](https://pmc.ncbi.nlm.nih.gov/articles/PMC3918453/#:~:text=We%20performed%20mRNA%20and%20small,%C3%97%2010%E2%88%9216%20for%20mRNA%2C%20p)). Notably, these were a subset of the 1000 Genomes Project samples, meaning each participant had well-characterized genetic backgrounds ([
            Transcriptome and genome sequencing uncovers functional variation in humans - PMC
        ](https://pmc.ncbi.nlm.nih.gov/articles/PMC3918453/#:~:text=In%20this%20study%2C%20we%20characterize,properties%20of%20causal%20functional%20variants)) ([
            Editors’ pick: transcriptomes of 1000 genomes - PMC
        ](https://pmc.ncbi.nlm.nih.gov/articles/PMC3831822/#:~:text=time%2C%20used%20the%20full%20capacity,Italians%2C%20and%20Yoruba%20from%20Nigeria)).

- **Tissue Source:** All gene expression measurements were performed on **lymphoblastoid cell lines (LCLs)**. These LCLs are immortalized cell lines derived from participants’ blood, providing a consistent *in vitro* tissue across all individuals ([
            Transcriptome and genome sequencing uncovers functional variation in humans - PMC
        ](https://pmc.ncbi.nlm.nih.gov/articles/PMC3918453/#:~:text=We%20performed%20mRNA%20and%20small,%C3%97%2010%E2%88%9216%20for%20mRNA%2C%20p)). Focusing on a single cell type ensured that differences observed were due to genetic and regulatory variation rather than different tissues.

- **Data Types:** Geuvadis generated **whole transcriptome RNA sequencing (RNA-seq)** data for each individual, including both messenger RNA (mRNA) and microRNA (**miRNA**) profiles ([
            Transcriptome and genome sequencing uncovers functional variation in humans - PMC
        ](https://pmc.ncbi.nlm.nih.gov/articles/PMC3918453/#:~:text=understanding%20the%20genetic%20basis%20of,to%20infer%20putative%20causal%20variants)) ([
            Transcriptome and genome sequencing uncovers functional variation in humans - PMC
        ](https://pmc.ncbi.nlm.nih.gov/articles/PMC3918453/#:~:text=We%20performed%20mRNA%20and%20small,%C3%97%2010%E2%88%9216%20for%20mRNA%2C%20p)). In total, high-throughput sequencing was done for both **mRNA** (capturing gene expression levels and splicing) and **small RNA** (capturing miRNA expression) from each sample. Complementing the RNA-seq, **genotype data** for each individual were obtained from the 1000 Genomes project or high-density SNP arrays ([
            Transcriptome and genome sequencing uncovers functional variation in humans - PMC
        ](https://pmc.ncbi.nlm.nih.gov/articles/PMC3918453/#:~:text=%28TSI%29%20and%20Yoruba%20%28YRI%29,coding%20and)). In fact, 421 of the 462 individuals had full genome sequences available (Phase 1 of 1000 Genomes), and the remainder’s genotypes were imputed from SNP chips ([
            Transcriptome and genome sequencing uncovers functional variation in humans - PMC
        ](https://pmc.ncbi.nlm.nih.gov/articles/PMC3918453/#:~:text=%28TSI%29%20and%20Yoruba%20%28YRI%29,coding%20and)). This provided a complete genomic variation map alongside the expression data. The combination of **paired genetic and expression data** in multiple populations made the Geuvadis dataset extremely valuable for eQTL mapping and integrative analyses ([
            Editors’ pick: transcriptomes of 1000 genomes - PMC
        ](https://pmc.ncbi.nlm.nih.gov/articles/PMC3831822/#:~:text=Altogether%2C%20Lappalainen%20et%20al,The%20authors%20concluded%20that%20in)).

### Relation to the GeneSPLS Project

The **GeneSPLS project** utilizes data from Geuvadis, effectively leveraging a *subset* of this larger consortium dataset for its analyses. In practice, this means that GeneSPLS did not collect new experimental data but rather uses the existing Geuvadis resource (for example, the gene expression levels and corresponding genotypes for a selection of **373 overlapping individuals**) as input to its models. By drawing on a portion of Geuvadis data, GeneSPLS can demonstrate and validate its methodology on real human genomic data while benefiting from the rigorous quality control and population diversity of the Geuvadis project. It should be emphasized that any results or figures in GeneSPLS are grounded in the Geuvadis dataset – a well-established reference. 

## Mathematical Formulation

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
