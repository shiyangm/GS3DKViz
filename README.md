# GS3DKViz
Vignette to visualize the promoter-enhancer links corresponding to the 3D window with minimum p-value using the Gviz and GenomicInteractions packages in R.

The code to generate promoter-enhancer interaction plots is available as a vignette.

## Description
GS3DKViz is vignette for visualizing the promoter-enhancer interaction plots from GeneScan3DKnock analysis.

## promoter-enhancer links
![Workflow](https://user-images.githubusercontent.com/57265092/99107266-8c690a80-25b3-11eb-8fe1-ceb388bffa38.jpg)

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

Shiyang Ma, James Dalgleish, Justin Lee, Chen Wang, Linxi Liu, Richard Gill, Wendy Chung, Hugues Aschard, Edwin K. Silverman, Michael H. Cho, Zihuai He, Iuliana Ionita-Laza, "A unified knockoff framework for gene-based testing with joint analysis of coding and regulatory variation", 2020

## License
This software is licensed under GPLv3.
