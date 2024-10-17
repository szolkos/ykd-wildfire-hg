# Yukon-Kuskokwim Delta (YKD) wildfire mercury (Hg) emissions
## Introduction
Source code for the manuscript: Substantial mercury releases and local deposition from permafrost peatland wildfires in southwestern Alaska. Supporting data are available from [Zolkos et al. 2024](https://doi.org/10.7910/DVN/ZA5JBZ).

## Authors
- [Scott Zolkos](https://www.researchgate.net/profile/Scott-Zolkos)
- Benjamin M. Geyman
- Stefano Potter
- Michael Moubarak
- Brendan M. Rogers
- Natalie Baillargeon
- Sharmila Dey
- Sarah M. Ludwig
- Sierra Melton
- Edauri Navarro
- Ann McElvein
- Prentiss H. Balcom
- Susan M. Natali
- Seeta Sistla
- Elsie M. Sunderland

## Data
### The following data associated with this manuscript is available in this repository:
#### Data used in generation of results
- *polaris_data.csv* (called in script '1_compile.R'): Soils data from The Polaris Project (henceforth, 'Polaris') used in this study.  

- *thg_ykd_al.csv* (called in script '1_compile.R'): Hg data measured on Polaris soils for this study.  

- *ykd_veg_biomass.csv* (called in script '1_compile.R'): Polaris data on vegetation biomass in the YK-Delta (collected during 2017 & 2018 field expeditions).  

- *ykd_veg_hgt.csv* (called in script '1_compile.R'): Hg data measured on Polaris vegetation for this study.  

- *ykd_hg_ub_site_locs_validated_2023_10_06.csv* (called in script '1_compile.R'): Data on soil characteristics and SOC & Hg content, for use in calculating averages for sites within the same 30 m Landsat pixel.   

- *moubarak_ykd_2015_burn_depth_2023_01_24.csv* (called in script '2_statistics.R'): Burn depth data from Moubarak et al. 2023 Biogeosciences.  

- *burn_depth_preds_2023_01_20.csv* (called in script '5_burn_depth_rf_mdl.R'): Environmental data derived from satellite remote sensing for predicting burn depth; data extracted using Google Earth Engine script '3_burn_depth_predictors.js'.  

- *burn_depth_preds_gee_2023_01_22.csv* (called in script '5_burn_depth_rf_mdl.R'): Large (80 MB) dataset of per-pixel (n > 400,000) burn depth predictions and values of satellite remote sensing indices used to burn depth. Request from lead author.  

- *preds_for_soc_modeling_2023_10_06.csv* (called in script '6_soc_rf_mdl.R'): Environmental data derived from satellite remote sensing for predicting SOC; data extracted using Google Earth Engine script '4_soc_predictors.js'.  

- *soc_preds_from_gee_2023_10_06.csv* (called in script '6_soc_rf_mdl.R'): Large (137 MB) dataset of per-pixel (n > 400,000) SOC predictions and values of satellite remote sensing indices used to predict SOC. Request from lead author.  

## Scripts
#### *Note: These scripts process the data summarized above, which was compiled as detailed in the Methods of the manuscript. Users of these scripts will need to update the working directory pathway within scripts.*  
- *1_compile.R*: Load and prepare data on soils (from unburned sites) and vegetation for subsequent analyses.  

- *2_statistics.R*: Compute statistics reported in the manuscript.  

- *3_burn_depth_predictors.js*: Google Earth Engine script for deriving environmental conditions from satellite remote sensing data and extracting the environmental data for use in development of predictive models of burn depth.  

- *4_soc_predictors.js*: Google Earth Engine script for deriving environmental conditions from satellite remote sensing data and extracting the environmental data for use in development of predictive models of soil organic carbon stores.

- *5_burn_depth_rf_mdl.R*: Algorithm to predict YK Delta wildfire burn depth from field measurements of burn depth (Moubarak et al. 2022) and various remote sensing metrics of fire severity and terrain qualities.  

- *6_soc_rf_mdl.R*: Algorithm to predict YK Delta wildfire carbon stores from field measurements of carbon and various remote sensing metrics of vegetation and terrain characteristics.  

- *7_hg_emissions_uncertainty.R*: Estimate uncertainty in wildfire Hg emissions.  

- *8_peat_burned_area.R*: Estimate northern peatland wildfire Hg emissions from 2002-2022.
