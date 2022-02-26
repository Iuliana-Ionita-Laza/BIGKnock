[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/Iuliana-Ionita-Laza/BIGKnock/workflows/R-CMD-check/badge.svg)](https://github.com/Iuliana-Ionita-Laza/BIGKnock/actions)

# BIGKnock 
This is an R package of performing biobank-scale gene-based association test via knockoffs.

## Description
The package contain functions for knockoff generation of gene and enhancers under biobank-scale data and conduct gene-based association test for related samples under fitted null GLMM.

## Prerequisites
R (recommended version >= 3.6.0)

## Dependencies
BIGKnock depends on R packages SKAT, Matrix, MASS, SPAtest, CompQuadForm and irlba. Make sure to install those packages before installing BIGKnock.
    
[SAIGE/SAIGE-Gene package](https://github.com/weizhouUMICH/SAIGE) is used for creating sparse GRM and fitting the null GLMMs.

## Installation
library(devtools) 

devtools::install_github("Iuliana-Ionita-Laza/BIGKnock")

The current version is 0.1 (February 21, 2022).

## Usage
Please see the BIGKnock <a href="https://github.com/Iuliana-Ionita-Laza/BIGKnock/blob/master/docs/BIGKnock_0.1.pdf"> **user manual** </a> for detailed usage of BIGKnock package. 


We implement [SAIGE/SAIGE-Gene](https://github.com/weizhouUMICH/SAIGE) package for Step 1: Construct the sparse GRM and Step 2: Fit the null GLMMs. Please see the <a href="https://htmlpreview.github.io/?https://github.com/Iuliana-Ionita-Laza/BIGKnock/blob/master/docs/BIGKnock_vignette.html">**tutorial**</a> for more details.

## Contact
If you have any questions about BIGKnock please contact

- <sm4857@cumc.columbia.edu>

If you want to submit a issue concerning the software please do so using the **BIGKnock** [Github repository](https://github.com/Iuliana-Ionita-Laza/BIGKnock/issues).


## Citation
* Ma, S., Wang, C., Khan, A., Liu, L., Dalgleish, J., Kiryluk, K., He, Z. and Ionita-Laza, I. (2022). "Causal gene prioritization via knockoff analysis of biobank-scale data with applications to UK Biobank". 
* He, Z., Guen, Y. L., Liu, L., Lee, J., Ma, S., Yang, A. C.,  Liu. X., Rutledge, J., Losada, P. M., Song, B., Belloy, M. E., Butler, R. R., Longo, F. M., Tang, H., Mormino, E. C., Wyss-Coray, T., Greicius, M. D. and Ionita-Laza, I. (2021) ["Genome-wide analysis of common and rare variants via multiple knockoffs at biobank scale, with an application to Alzheimer disease genetics".](https://doi.org/10.1016/j.ajhg.2021.10.009) _American Journal of Human Genetics_, **108**, 2336-2353.


## License
This software is licensed under GPL-3.
