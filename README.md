# GS3DKViz
Vignette to visualize the promoter-enhancer links corresponding to the 3D window with minimum p-value using the Gviz and GenomicInteractions packages in R.

The code to generate promoter-enhancer interaction plots is available as a vignette.

## Description
GS3DKViz is vignette for visualizing the promoter-enhancer interaction plots from GeneScan3DKnock analysis.

## Prerequisites
R (recommended version 4.0.3)

## Dependencies
Gviz,GenomicInteractions,magrittr, S4Vectors,GenomicRanges,IRanges,rlang,dplyr,tibble,tidyr,InteractionSet,janitor

## Installation
setRepositories(ind=c(1,2)) #this allows for Bioconductor dependencies to be installed correctly

remotes::install_github("Iuliana-Ionita-Laza/GS3DKViz",dependencies=T)

## Version
The current version is '1.0.2' (January 13, 2021).

## Citation

Ma, S., Dalgleish, J. L ., Lee, J., Wang, C., Liu, L., Gill, R., Buxbaum, J. D., Chung, W., Aschard, H., Silverman, E. K., Cho, M. H., He, Z. and Ionita-Laza, I. "Improved gene-based testing by integrating long-range chromatin interactions and knockoff statistics". medRxiv 2021.07.14.21260405 (2021）doi: 10.1101/2021.07.14.21260405


## License
This software is licensed under GPLv3.
