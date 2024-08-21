# Yukon-Kuskokwim Delta (YKD) wildfire mercury (Hg) emissions
## Introduction
Source code for the manuscript: Substantial mercury releases and local deposition from permafrost peatland wildfires in southwestern Alaska. 

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
- *.csv*: .  

- *.csv*: .  

- *.csv*: .  

- *.csv*: .  

- *.csv*: .  

- *.csv*: .  

- *.csv*: .  

- *.csv*: .  
 
- *.csv*: .  

- *.csv*: .  

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
