# Visceral Leishmaniasis Transmission Dynamics in Kenya 🌍🦟

**Institution:** International Centre of Physiology and Ecology (icipe) | Eneza Training Program  
**Scope:** Subcounty-level spatial epidemiology and predictive modeling (2020-2025)

## 📌 Project Overview
This repository contains the data science and spatial modeling pipeline analyzing the environmental and socioeconomic drivers of Visceral Leishmaniasis (VL) in Kenya. The project transitions from a rigorous Exploratory Data Analysis (EDA) of spatial clustering into a machine learning predictive framework to identify undocumented disease hotspots.

## 📊 Dataset & Scope
* **Confirmed Cases:** 6,977 cases across 67 sub-counties.
* **Spatial Clustering:** Significant spatial autocorrelation confirmed via Moran's I (0.2246, p < 0.001), isolating the epicenter in West Pokot (Sigor).
* **Predictors Used (11):** * **Environmental:** Temperature, Rainfall, RH, Greenness, Wetness (Sourced via Earth Engine, CHIRPS, ERA5).
  * **Socioeconomic:** Population Density, NightLights, Wealth Index, Healthcare Access (km), Water Distance (km) (Sourced via WorldPop, HDX, EOG).

## 🚀 Pipeline Architecture
1. **Exploratory Data Analysis (EDA):** Temporal seasonality profiling, KDE hotspot mapping, and socioeconomic vulnerability intersections.
2. **Feature Engineering:** Multicollinearity resolution (VIF validation) and variable selection.
3. **Predictive Modeling (In Progress):** * Population-weighted pseudo-absence generation.
   * 3-month lagged environmental data extraction to capture vector breeding lifecycles.
   * Baseline Logistic Regression (GLM) for interpretable odds ratios.
   * Tree-based algorithms (Random Forest, XGBoost) for non-linear threshold mapping.

## 📁 Repository Structure
* `/outputs`: Contains the 11 finalized data visualizations and spatial maps generated during the EDA phase.
* `VL_Kenya_Model_Dataset_v6_COMPLETE.csv`: The cleaned, 16-column modeling dataset.
*(Note: Massive raster covariates and GADM shapefiles are intentionally excluded via `.gitignore` to preserve repository stability).*
