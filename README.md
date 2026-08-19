# Spatial O-GalNAc glycoproteomics reveals region-specific regulation across the mouse brain

## Contents

- [Overview](#overview)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [GlycoDB](#GlycoDB)
- [License](./LICENSE)
- [Issues](https://github.com/mkumf/glyco_mousebrain_p21/issues)
- [Citation](#citation)

# Overview

Repository for the project *Spatial O-GalNAc glycoproteomics reveals region-specific regulation across the mouse brain*, including scripts used for data pre-/processing and for generating the central figures.

# Repo Contents

- [preprocessing](./preprocessing): `R`/`Delphi` scripts for data preprocessing. 
- [processing_plots](./processing_plots): `R` scripts for glycopeptide abundance data processing and for the figures presented in the paper.
- [manual_glycodb](./manual_glycodb): manual for the website associated with the paper, [GlycoDB.org](https://glycodb.org/).

# System Requirements

## Hardware Requirements

The analysis was done on a 2019 MacBook Pro  
RAM: 16 GB  
CPU: 2,4 GHz Quad-Core Intel Core i5

## Software Requirements

The `R` scripts were run using  

MacOS: Sequoia 15.7.9  
RStudio: 2026.05.0+218  
R: 4.5.3 (2026-03-11)  
Bioconductor: 3.22

All packages and dependencies were installed from [CRAN](https://cran.r-project.org/web/packages/index.html) respectively [Bioconductor](https://www.bioconductor.org/packages/release/BiocViews.html#___Software) (if available) or directly from GitHub using the respective instructions.

# GlycoDB

The results can be accessed online at [GlycoDB.org](https://glycodb.org/).

# Citation

If you find this work useful please cite the paper.
