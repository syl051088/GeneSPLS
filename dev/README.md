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

The **Geuvadis Consortium** (short for *Genetic EUropean VAriation in DISease*) was an EU-funded research collaboration launched in 2010, aimed at exploring how genetic variation impacts gene expression and human disease ([CORDIS: GEUVADIS Project Report](https://cordis.europa.eu/project/id/261123/reporting/de)).

### Geuvadis Dataset Characteristics

One of the major outcomes of Geuvadis is a **large-scale reference dataset** integrating genomic and transcriptomic information. This dataset is widely used in human genomics as it provides a rich resource for studying **expression quantitative trait loci (eQTLs)** and gene regulation. Key features of the Geuvadis dataset include:

- **Cohort Size & Populations:** The project profiled **≈465 individuals** (after quality control, 462 remained) drawn from **five distinct populations**. These populations were the CEPH Utah residents (**CEU**), Finns (**FIN**), British (**GBR**), Toscani Italians (**TSI**), and Yoruba Nigerians (**YRI**), with roughly 89–95 individuals per group. Notably, these were a subset of the 1000 Genomes Project samples, meaning each participant had well-characterized genetic backgrounds ([PMC: Transcriptomes of 1000 Genomes](https://pmc.ncbi.nlm.nih.gov/articles/PMC3831822/)).

- **Tissue Source:** All gene expression measurements were performed on **lymphoblastoid cell lines (LCLs)**. These LCLs are immortalized cell lines derived from participants’ blood, providing a consistent *in vitro* tissue across all individuals. Focusing on a single cell type ensured that differences observed were due to genetic and regulatory variation rather than different tissues.

- **Data Types:** Geuvadis generated **whole transcriptome RNA sequencing (RNA-seq)** data for each individual, including both messenger RNA (mRNA) and microRNA (**miRNA**) profiles. In total, high-throughput sequencing was done for both **mRNA** (capturing gene expression levels and splicing) and **small RNA** (capturing miRNA expression) from each sample. Complementing the RNA-seq, **genotype data** for each individual were obtained from the 1000 Genomes project or high-density SNP arrays ([PMC: Genotype Data](https://pmc.ncbi.nlm.nih.gov/articles/PMC3918453/)). In fact, 421 of the 462 individuals had full genome sequences available (Phase 1 of 1000 Genomes), and the remainder’s genotypes were imputed from SNP chips. This provided a complete genomic variation map alongside the expression data. The combination of **paired genetic and expression data** in multiple populations made the Geuvadis dataset extremely valuable for eQTL mapping and integrative analyses.

### Relation to the GeneSPLS Project

The **GeneSPLS project** utilizes data from Geuvadis, effectively leveraging a *subset* of this larger consortium dataset for its analyses. In practice, this means that GeneSPLS did not collect new experimental data but rather uses the existing Geuvadis resource (for example, the gene expression levels and corresponding genotypes for a selection of **373 overlapping individuals**) as input to its models. By drawing on a portion of Geuvadis data, GeneSPLS can demonstrate and validate its methodology on real human genomic data while benefiting from the rigorous quality control and population diversity of the Geuvadis project. It should be emphasized that any results or figures in GeneSPLS are grounded in the Geuvadis dataset – a well-established reference. 

## Mathematical Formulation

---

### 1. Theoretical Overview of eQTL Analysis

In expression quantitative trait loci (eQTL) analysis, the objective is to identify associations between genetic variants (SNPs) and gene expression levels. This relationship is typically modeled using a standard linear regression framework:

$$
{\huge
Y = X B + F,
}
$$

where:
- $X \in \mathbb{R}^{n \times p}$ is the genotype matrix (with $n$ individuals and $p$ SNPs),
- $Y \in \mathbb{R}^{n \times q}$ is the gene expression matrix (with $q$ genes, where $q$ can be 1 or more),
- $B \in \mathbb{R}^{p \times q}$ is the matrix of regression coefficients,
- $F \in \mathbb{R}^{n \times q}$ represents the error terms.

The challenge in eQTL analysis arises because $p$ is typically much larger than $n$, and the SNPs (columns of $X$) are highly correlated. Dimensionality reduction methods such as Partial Least Squares (PLS) offer an effective solution in this context.

---

### 2. Partial Least Squares (PLS) Regression with the Wide Kernel Algorithm

PLS regression seeks to find latent components that capture the covariance between $X$ and $Y$. In our context, this helps model the relationship between SNPs and gene expression. The model assumes that both $X$ and $Y$ can be decomposed as:

$$
{\huge
X = T P^T + E, \quad Y = T Q^T + F,
}
$$

where:
- $T \in \mathbb{R}^{n \times K}$ is the matrix of latent scores,
- $P \in \mathbb{R}^{p \times K}$ and $Q \in \mathbb{R}^{q \times K}$ are the loading matrices,
- $E$ and $F$ are residual matrices,
- $K$ is the number of latent components.

#### Wide Kernel PLS Derivation

For datasets where $p \gg n$, the wide kernel PLS algorithm avoids direct computations in the high-dimensional space by utilizing kernel matrices. The key steps are:

1. **Kernel Formation:**

Compute the cross-product matrices:
   
$$
{\huge
XX^T \quad \text{and } \quad YY^T.
}
$$

Form the composite kernel matrix:

$$
{\huge
M_{\text{kernel}} = (XX^T)(YY^T).
}
$$

2. **Latent Score Computation:**

Solve the eigenvalue problem:

$$
{\huge
M_{\text{kernel}} t = \lambda t,
}
$$

where $t$ is a score vector (a column of $T$) and $\lambda$ is the associated eigenvalue. The iterative process, often implemented using the power method, extracts these score vectors efficiently.

3. **Weight Vector and Loadings:**  

Once a score vector $t$ is obtained, compute the corresponding weight vector:

$$
{\huge
w = \frac{X^T t}{\|X^T t\|}.
}
$$

The loading matrices are then derived as:
$$
{\huge
P = X^T T (T^T T)^{-1}, \quad Q = Y^T T (T^T T)^{-1}.
}
$$

4. **Regression Coefficients:**

Finally, the PLS regression coefficients are given by:

$$
{\huge
B_{\text{PLS}} = W (P^T W)^{-1} Q^T,
}
$$

where $W$ is the matrix whose columns are the weight vectors $w$.

#### Connection to `pls.cpp`

In our source code (`pls.cpp`), the function `widekernelpls_fit` implements the wide kernel PLS algorithm:
- **Input:**  
  - $X$ (predictor matrix) and $Y$ (response matrix).
  - $ncomp$ (number of latent components, $K$).
  - Options to center $X$ and $Y$ (outputs stored as `Xmeans` and `Ymeans`).
- **Processing:**  
  - Centers $X$ and $Y$ and computes $XX^T$ and $YY^T$.
  - Iteratively computes latent score vectors $T$ using the kernel $XX^T YY^T$ with normalization.
  - Derives the weight matrix $W = X^T U$ and computes the projection matrix $R$ that maps $X$ to the latent space.
  - Computes regression coefficients for each number of components, stored in `coefficients`.
- **Output:**  
  - `coefficients`: Represents $B_{\text{PLS}}$.
  - `projection`: The projection matrix $R$.
  - `Xmeans` and `Ymeans`: The centering constants used.

---

### 3. Sparse Partial Least Squares (SPLS) Regression

To further refine the model for eQTL analysis, SPLS introduces sparsity in the weight vectors, thereby selecting only the most relevant SNPs. The SPLS model is formulated as:

$$
{\huge
\min_{w,c} \left\lbrace -\kappa\, w^T M w + \frac{1-\kappa}{2}\|c-w\|^2 + \lambda_1 \|c\|_1 + \lambda_2 \|c\|_2^2 \right\rbrace \quad \text{subject to } \|w\|_2 = 1,
}
$$

where:
- $M = X^T Y Y^T X$,
- $w$ is the direction vector,
- $c$ is an auxiliary vector that is encouraged to be sparse,
- $\kappa$ controls the trade-off between maximizing covariance and matching $w$ to $c$,
- $\lambda_1$ and $\lambda_2$ are the regularization parameters for the $L_1$ and $L_2$ penalties, respectively.

Cross-validation (via `cv_spls_cpp`) is used to determine the optimal values for the sparsity (typically through the candidate eta values) and the number of latent components $K$.

#### Derivation Details for SPLS

The derivation proceeds by alternating between:
1. **Optimizing $w$ with $c$ fixed:**

Solve:

$$
{\huge
\min_{w} \left( -\kappa\, w^T M w + \frac{1-\kappa}{2}\|c - w\|^2 \right) \quad \text{subject to } \|w\|_2 = 1
}
$$

This can be tackled via Lagrange multipliers, leading to an update that resembles an eigen-decomposition of an adjusted matrix.

2. **Optimizing $c$ with $w$ fixed:**  

Solve:

$$
{\huge
\min_{c} \left( \frac{1-\kappa}{2} \|c - w\|^2 + \lambda_1 \|c\|_1 + \lambda_2 \|c\|_2^2 \right)
}
$$
   
This step is analogous to solving an elastic net problem and can be solved by soft-thresholding, which yields many zero entries in $c$.

The iterative process converges to a sparse estimate of the direction vector, which in turn leads to a sparse regression coefficient matrix $B_{\text{SPLS}}$.

#### Connection to `spls.cpp`

In our `spls.cpp` implementation, the function `spls_cpp` performs the SPLS algorithm:
- **Input:**  
  - $x$ and $y$ (predictor and response matrices).
  - $K$ (number of latent components).
  - Sparsity parameter `eta` (analogous to $\lambda_1$), along with $\kappa$, and additional parameters controlling convergence.
- **Processing:**  
  - The function centers and scales the data.
  - It then alternates between updating the weight vector $w$ and the sparse surrogate $c$ using the derivations outlined above.
  - An active set $A$ is determined based on the non-zero elements in $c$ (or the estimated coefficients).
  - For the subset of variables in $A$, a restricted PLS (via `widekernelpls_fit`) is performed to obtain the final regression coefficients.
- **Output:**  
  - `betahat`: The sparse regression coefficients $B_{\text{SPLS}}$.
  - `A`: The active set of predictors (indices corresponding to non-zero loadings).
  - `projection`: The projection matrix computed on the reduced set.

In GeneSPLS, cross-validation (via `cv_spls_cpp`) is used to determine the optimal values for the sparsity (typically through the candidate eta values) and the number of latent components $K$.

---

### 4. Integration in the eQTL Context

Both PLS and SPLS serve to overcome the challenges in eQTL analysis by reducing dimensionality and, in the case of SPLS, selecting the most relevant SNPs. The wide kernel approach in PLS efficiently computes the latent structure even when $p \gg n$, while SPLS refines the model by enforcing sparsity. In our GeneSPLS package, the implementations in `pls.cpp` and `spls.cpp` follow these theoretical derivations closely, ensuring that the outputs (regression coefficients, projection matrices, centering means, and active sets) are directly interpretable in terms of the underlying mathematical models.

---

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

---

## Results from test_real_univariate.Rmd

In this analysis, the GeneSPLS model was applied to real eQTL data focusing on the expression of the AKT3 gene. The primary aim was to identify SNPs that exhibit non-zero effects on AKT3 expression. After tuning the SPLS models over a range of sparsity levels (optimal `η` from 0.95 to 0.99 based on cross-validation), a total of 7 SNPs were selected as important features. Notably, **4 SNPs** (listed in **bold** in the table below) were consistently selected across all SPLS model fits, indicating robust associations. These SNPs map to two distinct genomic regions: three lie on chromosome 1p36 and the remainder on chromosome 9p22, as summarized in **Table 1**.

**Table 1.** Genomic location, gene annotation, and distance to AKT3 of the 7 selected SNPs (GRCh38). SNPs in **bold** were consistently selected in all SPLS models.

| **SNP (rsID)**   | **Chromosome: Position** | **Gene (Annotation)**       | **Distance to AKT3**                          |
|------------------|--------------------------|-----------------------------|-----------------------------------------------|
| **rs6691809**    | 1:7008964                | CAMTA1 (intronic variant)   | ~236.5 Mb upstream of AKT3                    |
| **rs6577402**    | 1:7014668                | CAMTA1 (intronic variant)   | ~236.5 Mb upstream of AKT3                    |
| rs2412208       | 1:7032722                | CAMTA1 (intronic variant)   | ~236.5 Mb upstream of AKT3                    |
| **rs693196**     | 9:15202419               | TTC39B (intronic variant)   | Not applicable (different chromosome)         |
| **rs7876026**    | 9:15202817               | TTC39B (intronic variant)   | Not applicable (different chromosome)         |
| rs508314        | 9:15189274               | TTC39B (intronic variant)   | Not applicable (different chromosome)         |
| rs1407977       | 9:15188108               | TTC39B (missense variant)   | Not applicable (different chromosome)         |

Most of the chromosome 1 SNPs lie within the _CAMTA1_ gene (Calmodulin-binding transcription activator 1). _CAMTA1_ is a known anti-tumor gene that plays a role in cell cycle regulation by inhibiting AKT phosphorylation ([CAMTA1 calmodulin binding transcription activator 1 [Homo sapiens (human)] - Gene - NCBI](https://www.ncbi.nlm.nih.gov/gene/23261#:~:text=Title%3A%20CAMTA1%2C%20a%20novel%20antitumor,glioma%20by%20inhibiting%20AKT%20phosphorylation)). This is biologically relevant because phosphorylation of AKT is a key activation step in the AKT signaling pathway, which promotes cell growth and survival. Higher levels of phospho-AKT are associated with more aggressive cell proliferation ([
            AKT3 promotes prostate cancer proliferation cells through regulation of Akt, B-Raf & TSC1/TSC2 - PMC
        ](https://pmc.ncbi.nlm.nih.gov/articles/PMC4694976/#:~:text=match%20at%20L354%20,of%20prostate%20cancer%20cells%20both)), and AKT3 (one of the three AKT isoforms) is expressed in many tissues and contributes to cell proliferation – silencing AKT3 can suppress cell growth in cancer models ([
            AKT3 promotes prostate cancer proliferation cells through regulation of Akt, B-Raf & TSC1/TSC2 - PMC
        ](https://pmc.ncbi.nlm.nih.gov/articles/PMC4694976/#:~:text=match%20at%20L354%20,of%20prostate%20cancer%20cells%20both)). The _CAMTA1_ locus is located on chromosome 1p36 (around 7.0 Mb ([GWAS Catalog - EMBL-EBI](https://www.ebi.ac.uk/gwas/search?query=CAMTA1#:~:text=GWAS%20Catalog%20,V%20rs11588880%20%28rs147096589%29))), approximately 236.5 Mb upstream of the _AKT3_ gene at 1q43-q44 (around 243.5 Mb ([Gene: AKT3 (ENSG00000117020) - Summary - Homo_sapiens - Ensembl genome browser 113](http://www.ensembl.org/id/ENSG00000117020#:~:text=Location))). While _CAMTA1_ and _AKT3_ are not neighboring genes, the fact that _CAMTA1_ can modulate AKT activity provides a plausible functional link between the chromosome 1p36 variants and AKT-related cell cycle control. In other words, the SPLS-selected SNPs in _CAMTA1_ point to a mechanism whereby altered CAMTA1 activity might influence AKT phosphorylation status and thus downstream AKT3 signaling.

By contrast, the remaining SNPs reside on chromosome 9 in the _TTC39B_ gene (tetratricopeptide repeat domain 39B) at 9p22.3. Two of these (**rs693196** and **rs7876026**) were consistently selected across all models (Table 1), and all four chromosome 9 variants fall within intronic regions of _TTC39B_ (except for rs1407977, which is a missense variant in an exon). _TTC39B_ encodes a protein involved in metabolic processes – for instance, it has been implicated in regulation of cholesterol homeostasis ([TTC39B tetratricopeptide repeat domain 39B [Homo sapiens (human)] - Gene - NCBI](https://www.ncbi.nlm.nih.gov/gene/158219#:~:text=Summary%20Predicted%20to%20be%20involved,of%20Genome%20Resources%2C%20Apr%202025)) – and is not known to be directly connected to the AKT pathway. The enrichment of SPLS-selected SNPs in _TTC39B_ therefore suggests a second, independent locus that may influence the outcome through a different biological mechanism. In summary, our SPLS analysis identified two genomic regions of interest: a chromosome 1 locus in _CAMTA1_ linked to cell cycle regulation via the AKT signaling axis, and a chromosome 9 locus in _TTC39B_ that may represent a novel association deserving further investigation. 

---

## Conclusion

The present study aimed to leverage the GeneSPLS methodology for eQTL analysis, specifically targeting the AKT3 gene. However, the genomic annotations of the SNPs selected by the model reveal a discrepancy: the non-zero coefficient SNPs identified by GeneSPLS are not in the cis-regulatory region of AKT3. Instead, the four chromosome 1 SNPs are located in the CAMTA1 gene region, which is approximately 236 Mb upstream of AKT3, and the two SNPs on chromosome 9 belong to a different genomic context (within TTC39B).

This outcome suggests several potential interpretations:

1. **Trans-eQTL Signals:** Although cis-eQTLs are generally expected to exhibit the strongest associations with gene expression, the possibility of trans-regulatory effects cannot be entirely dismissed. However, trans-eQTLs are usually subtler and require robust validation in larger datasets. The strong signals observed here may instead reflect model artifacts rather than genuine trans effects.

2. **Model and Data Considerations:** The divergence from the anticipated cis-regulatory pattern may indicate issues with data preprocessing, such as overly stringent filtering thresholds or batch effects, or limitations in the parameter tuning process. It is critical to ensure that the filtering steps have not inadvertently removed true cis-variants for AKT3. The tuning parameters for sparsity and the number of latent components might also need adjustment to capture the expected regulatory variants.

3. **Biological Implications:** If validated, the detected associations could hint at a more complex regulatory architecture for AKT3, potentially involving long-range chromosomal interactions or indirect regulatory networks. Nevertheless, given the substantial genomic distance between the selected SNPs and AKT3, the likelihood of such mechanisms should be carefully scrutinized.

In conclusion, while the GeneSPLS framework demonstrates robust computational performance and an innovative approach to dimensionality reduction in high-dimensional genomic data, the current output for AKT3 does not conform to conventional expectations of cis-eQTL mapping. Future work should include:
- A thorough review of data preprocessing and parameter tuning procedures.
- Independent validation of the detected SNPs using larger cohorts or alternative statistical methods.
- Exploration of potential trans-regulatory mechanisms if trans-eQTL signals are indeed present.

Overall, these findings emphasize the importance of integrating biological context with statistical modeling and suggest that further refinement and validation of the GeneSPLS approach are warranted before drawing definitive conclusions about the genetic regulation of AKT3.

---

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
