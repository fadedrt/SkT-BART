# ============================================================================
# Data Loading Script for Skew-t BART Examples
# ============================================================================

# ---------------------------------------------------------------------------
# Dataset 1: Biomass
# ---------------------------------------------------------------------------
# Source: modeldata package (CRAN)
# Description:
#   Let xi = (xi1, ..., xi5)T denote the elemental composition of the ith biomass fuel sample.
#   The dataset comprises n = 536 observations. The five elemental components are:
#     - Carbon (C)
#     - Hydrogen (H)
#     - Nitrogen (N)
#     - Oxygen (O)
#     - Sulfur (S)
#   The response variable yi represents the Higher Heating Value (HHV), which is an
#   indicator of the energy content of biomass fuels.
#   This dataset can be used to model HHV as a function of the elemental components
#   using the BART framework.

library(modeldata)
data(biomass)

# Extract predictors and response
x_biomass <- as.matrix(biomass[c(3, 4, 5, 6, 7)])
y_biomass <- biomass[[8]]

# ---------------------------------------------------------------------------
# Dataset 2: Forest Fires
# ---------------------------------------------------------------------------
# Source: UCI Machine Learning Repository
# URL: https://archive.ics.uci.edu/ml/datasets/Forest+Fires
# Description:
#   Let xi = (xi1, ..., xi7)T denote the vector of environmental and spatial predictors
#   for the ith observation in the Forest Fires dataset. The dataset has 517 records
#   collected from Montesinho Park. The covariates include:
#     - Spatial coordinates (X, Y)
#     - Initial Spread Index (ISI)
#     - Temperature (Temp)
#     - Relative humidity (RH)
#     - Wind speed (Wind)
#     - Rainfall (Rain)
#   The response variable yi represents the Fine Fuel Moisture Code (FFMC), a key
#   indicator for assessing forest fire danger levels.

# Assuming the CSV has been downloaded locally as "forestfires.csv"
forestfires <- read.csv("forestfires.csv")

# Extract predictors and response
x_fire <- as.matrix(forestfires[c("X", "Y", "ISI", "temp", "RH", "wind", "rain")])
y_fire <- forestfires$FFMC

# ---------------------------------------------------------------------------
# End of data.R
# ---------------------------------------------------------------------------

