# install_dependencies.R

# Run the script to automatically install all packages required for the AviCRED analysis:

required_packages <- c(
  "duckdb", "DBI", "dplyr", "tidyr", "data.table", 
  "isotree", "sf", "terra", "hexbin", "lattice","yardstick",
  "ggplot2", "treemapify", "ggalluvial", "ComplexUpset", 
  "scales", "extrafont", "patchwork", "ggridges", "caret"
)

# Checks if the user already has the package otherwise installs them:
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

print("All package dependencies have been verified and installed successfully")