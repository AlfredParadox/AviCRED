# install_dependencies.R

# Run the script to automatically install all packages required for the AviCRED analysis:

required_packages <- c(
  "caret", "ComplexUpset", "data.table", "DBI", "dbplyr", "DescTools",   
  "dplyr", "duckdb", "extrafont", "ggplot2", "ggspatial", "glmmTMB",     
  "grid", "Hmisc", "isotree", "mgcv", "patchwork", "pROC",        
  "purrr", "readr", "scales", "sf", "showtext", "stringr",    
  "sysfonts", "terra", "tidyr", "tidyverse", "treemapify", "viridis",    
  "yardstick"   
)

# Checks if the user already has the package otherwise installs them:
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

print("All package dependencies have been verified and installed successfully")