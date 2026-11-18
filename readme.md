# AviCRED: Avian Checklist Reliability Evaluation Data-quality Framework

This repository contains the R implementation of the AviCRED framework, designed to assess the data-quality and reliability of citizen-science checklists from the Indian EBD dataset (2010-25). 

*Note: This repository is currently structured for peer-review evaluation accompanying a manuscript submission.*

---

## Computational Environment & Requirements
This analysis was developed and executed using R (v4.4.0). 

### Core Package Versions
The primary package versions used to generate the results in the manuscript are:

* `terra_1.8-60`
* `patchwork_1.3.2`
* `sf_1.0-16`
* `ggalluvial_0.12.6`
* `ComplexUpset_1.3.3`
* `extrafont_0.20`
* `scales_1.4.0`
* `tidyr_1.3.1`
* `data.table_1.17.8`
* `caret_7.0-1`
* `ggridges_0.5.7`
* `lattice_0.22-6`
* `dplyr_1.2.1`
* `DBI_1.3.0`
* `dbplyr_2.6.0`
* `DescTools_0.99.60`
* `duckdb_1.5.6`
* `ggplot2_4.0.3`
* `ggspatial_1.1.11`
* `glmmTMB_1.1.15.2`
* `grid_4.6.1`
* `Hmisc_5.3-0`
* `isotree_0.6.1-5`
* `mgcv_1.9-4`
* `pROC_1.19.1`
* `purrr_1.2.2`
* `readr_2.2.0`
* `showtext_0.9-8`
* `stringr_1.6.0`
* `sysfonts_0.8.9`
* `tidyr_1.3.2`
* `tidyverse_2.0.0`
* `treemapify_2.6.1`
* `viridis_0.6.5`
* `yardstick_1.4.0`
* `ComplexUpset_1.3.3`
* `ggridges_0.5.7`
* `lattice_0.22-6`

---

## Getting Started & Installation

### 1. Clone or Download the Repository
Download this folder to your local machine and open the directory in RStudio.

### 2. Install Dependencies
Before running the scripts, execute the automated dependency script to ensure all required libraries are mapped and installed on your system:
```R

source("install_dependencies.R")