#================================================================================================================================#
# Filename: 7_hg_emissions_uncertainty.R ####
# Author: Scott Zolkos
# Created: 2023-10-03
# Revised: 2024-08-20
# Background: Estimate uncertainty in wildfire Hg emissions.
# References
## None
#================================================================================================================================#

#-----------------------------------#
# 1. Prepare working environment ####
#-----------------------------------#

# Install R package management tool (pacman), as needed
if("pacman" %in% installed.packages() == FALSE){install.packages("pacman")}

# Load necessary packages
pacman::p_load(tidyverse, # functionality
               tidyr,
               dplyr,
               purrr,
               lubridate,
               parallel,
               foreign # functionality | e.g., read .dbf files
)

#---------------------------#
# 2. Load & prepare data ####
#---------------------------#

# Read in burn depth and SOC data, which were exported from ArcGIS
bd_soc_df <- read.dbf(paste0(root_path, "/", "burn_depth_soc_2023_10_07.dbf"))
options(digits = 10)
bd_soc <- bd_soc_df %>% select(num, lat_dd, lon_dd, burn_depth, soc_pred)
names(bd_soc) <- c("num", "lat_dd", "lon_dd", "bd", "soc")

# Read in veg Hg and RHgC data
veg_hg_smry2 = read.csv(paste0(root_path, "/", "veg_hg_smry2.csv"), header=T)
zea = read.csv(paste0(root_path, "/", "zea.csv"), header=T)

# Add columns for summary stats
bd_soc$mean_hg_em <- NA
bd_soc$median_hg_em <- NA
bd_soc$sd_hg_em <- NA
bd_soc$ci_lwr_hg_em <- NA
bd_soc$ci_upr_hg_em <- NA
bd_soc$mean_soil_hg <- NA
bd_soc$median_soil_hg <- NA
bd_soc$sd_soil_hg <- NA
bd_soc$ci_lwr_soil_hg <- NA
bd_soc$ci_upr_soil_hg <- NA
    
#------------------------------------#
# 3. Perform uncertainty analysis ####
#------------------------------------#

# Subset dataset for testing
#bd_soc <- bd_soc[1:475,] # 0.1% of database

# PFT Hg summary data
extract_veg_hg_values <- function(plant_type) {
  hg_mean <- veg_hg_smry2 %>% filter(pft == plant_type) %>% select(hg_stores_mean) %>% as.numeric()
  hg_sd <- veg_hg_smry2 %>% filter(pft == plant_type) %>% select(hg_stores_sd) %>% as.numeric()
  return(list(mean = hg_mean, sd = hg_sd))
}

lichen_vals <- extract_veg_hg_values("lichen")
moss_vals <- extract_veg_hg_values("moss")
ev_shrub_vals <- extract_veg_hg_values("evergreen shrub")
dec_shrub_vals <- extract_veg_hg_values("deciduous shrub")
forb_vals <- extract_veg_hg_values("forb")
gram_vals <- extract_veg_hg_values("graminoid")
  
# Function to perform a single simulation of Hg release
mc_fire_hg_release <- function(burn_depth, SOC){

  # Soil Hg release
    # SOC stores and uncertainty (kgC/m2) for pixel
      # Generate random values so that no negative value is obtained
        # Initialize an empty data frame to store results
          soc_results <- data.frame(value = numeric(0))
        # Initialize a flag to check for negative values
          neg_flag <- TRUE
        # Generate random values until no negative values are obtained
          while(neg_flag){
            # Generate a single random value from the normal distribution
              soc <- rnorm(1, mean=SOC, sd=2.428046) # sd(obs_preds_val$Pred) # uncertainty in modeled SOC stores using SD of SOC predicted from LOOCV (2.428046 kgC/m2)
            # Append the random value to the results data frame
              soc <- bind_rows(soc_results, data.frame(value = soc))
            # Check if the generated value is negative
              neg_flag <- soc < 0
          }
      
    # RHgC- samples from measured values
      rhgc <- sample(zea$rhgc_ugg, size=1)
      
    # Calculate soil Hg stores (mgHg/m2/30cm) = RHgC (µgHg/gC) * SOC (kgC/m2) ####
      soil_hg <- rhgc*soc
      
    # Burn depth (cm)
      # Generate random values so that no negative value is obtained
        # Initialize an empty data frame to store results
          bd_results <- data.frame(value = numeric(0))
        # Initialize a flag to check for negative values
          neg_flag <- TRUE
        # Generate random values until no negative values are obtained
          while(neg_flag){
            # Generate a single random value from the normal distribution
              #fin_mdl # RMSE = 3.131946, 2.5% CI = quantile(obs_preds_val$Pred, 0.025), 95% CI = quantile(obs_preds_val$Pred, 0.975)
              #sd_bz <- sd(obs_preds_val$Pred)/30 # uncertainty in modeled burn depth using SD of burn depth predicted from LOOCV (2.643267 cm)
              bd <- rnorm(1, mean=burn_depth, sd=2.643267)
              bd <- bd/30 # ÷ 30 b/c Hg stores are for top 30 cm
            # Append the random value to the results data frame
              bd <- bind_rows(bd_results, data.frame(value = bd))
            # Check if the generated value is negative
              neg_flag <- bd < 0
          }
      
    # Soil Hg release
      soil_hg_release <- soil_hg * bd
      
  # Vegetation Hg release
    # Vegetation Hg stores (µg Hg/m2) observed in data from this study
      # Lichen
        li <- rnorm(1, lichen_vals$mean, lichen_vals$sd)
        li <- li[li>0] # ensure no Hg stores < 0
      # Moss
        mo <- rnorm(1, moss_vals$mean, moss_vals$sd)
        mo <- mo[mo>0]
      # Evergreen shrubs
        ev <- rnorm(1, ev_shrub_vals$mean, ev_shrub_vals$sd)
        ev <- ev[ev>0]
      # Deciduous shrubs
        de <- rnorm(1, dec_shrub_vals$mean, dec_shrub_vals$sd)
        de <- de[de>0]
      # Forbs
        fo <- rnorm(1, forb_vals$mean, forb_vals$sd)
        fo <- fo[fo>0]
      # Graminoids
        gr <- rnorm(1, gram_vals$mean, gram_vals$sd)
        gr <- gr[gr>0]
    
    # Fraction veg combusted (from Frost et al. 2020 ERL Table 2)
      # Lichen
        li_fvc <- rnorm(1, 1-(1.3/81.1), 0.031)
        li_fvc <- li_fvc[li_fvc <= 1] # set maximum possible value for proportion to 1
      # Moss
        mo_fvc <- rnorm(1, 1-(10/12.8), 0.068) # Sphagnum
        mo_fvc <- mo_fvc[mo_fvc > 0]
      # Evergreen shrubs
        ev_fvc <- rnorm(1, 1-(8.8/25.0), 0.060)
        ev_fvc <- ev_fvc[ev_fvc <= 1]
      # Deciduous shrubs
        de_fvc <- rnorm(1, 1-(8.8/12.2), 0.041)
        de_fvc <- de_fvc[de_fvc > 0]
      # Forbs
        fo_fvc <- 0
      # Graminoids
        gr_fvc <- 0

    # Veg Hg release (mg/m2); ÷ 1000 converts µg to mg
      veg_hg_release <- sum(
        ((li/1000)*li_fvc),
        ((mo/1000)*mo_fvc), 
        ((ev/1000)*ev_fvc), 
        ((de/1000)*de_fvc), 
        ((fo/1000)*fo_fvc), 
        ((gr/1000)*gr_fvc))
  
  # Total wildfire Hg release (mg/m2)
    hg_release = soil_hg_release + veg_hg_release
  
  # Return results
    return(list(hg_release, soil_hg))
}

# Set seed for reproducibility
set.seed(99)
      
# Set the number of Monte Carlo simulations
n_simulations <- 100

# Progress bar for uncertainty analysis
n_iter <- nrow(bd_soc)
pb <- txtProgressBar(min = 0,
                     max = n_iter,
                     style = 3,
                     width = n_iter/n_iter*20, # Needed to avoid multiple printings
                     char = "=")
init <- numeric(n_iter)
end <- numeric(n_iter)

# Run uncertainty analysis of wildfire Hg release for every pixel
j=1
for(i in 1:nrow(bd_soc)){
# progress bar
  init[j] <- Sys.time()
  
# Run Monte Carlo for n simulations
  #mc_results <- replicate(n_simulations, mc_fire_hg_release(burn_depth=bd_soc[i,4], SOC=bd_soc[i,5]))
  mc_results <- replicate(n_simulations, as.numeric(unlist(mc_fire_hg_release(burn_depth=bd_soc[i,4], SOC=bd_soc[i,5]))))
  hg_em <- mc_results[1,]
  soil_hg <- mc_results[2,]
  
# Convert results to tidy dataframe
  md_df <- as.data.frame(matrix(nrow=n_simulations, ncol=3))
  names(md_df) <- c("simulation","hg_release_mgm2","soil_hg")
  md_df[1] <- 1:n_simulations
  md_df[2] <- hg_em
  md_df[3] <- soil_hg
      
# Summarize results
  stats <- md_df %>%
    summarise(
      mean_hg_em = mean(hg_release_mgm2),
      median_hg_em = median(hg_release_mgm2),
      sd_hg_em = sd(hg_release_mgm2),
      ci_lwr_hg_em = quantile(hg_release_mgm2, 0.05),
      ci_upr_hg_em = quantile(hg_release_mgm2, 0.95),
      mean_soil_hg = mean(soil_hg),
      median_soil_hg = median(soil_hg),
      sd_soil_hg = sd(soil_hg),
      ci_lwr_soil_hg = quantile(soil_hg, 0.05),
      ci_upr_soil_hg = quantile(soil_hg, 0.95)
    )
  
# Add summary stats to df
  bd_soc[j,6] <- stats$mean_hg_em
  bd_soc[j,7] <- stats$median_hg_em
  bd_soc[j,8] <- stats$sd_hg_em
  bd_soc[j,9] <- stats$ci_lwr_hg_em
  bd_soc[j,10] <- stats$ci_upr_hg_em
  bd_soc[j,11] <- stats$mean_soil_hg
  bd_soc[j,12] <- stats$median_soil_hg
  bd_soc[j,13] <- stats$sd_soil_hg
  bd_soc[j,14] <- stats$ci_lwr_soil_hg
  bd_soc[j,15] <- stats$ci_upr_soil_hg
  
  # close progress bar, print % complete
    end[j] <- Sys.time()
    setTxtProgressBar(pb, j)
    time <- round(seconds_to_period(sum(end - init)), 0)
    est <- n_iter * (mean(end[end != 0] - init[init != 0])) - time
    remaining <- round(seconds_to_period(est), 0)
  # Estimated remaining time based on the mean time that took to run the previous iterations
    if(round((j/n_iter*100),2) %in% seq(1,100,1)){
      cat(paste(" // Execution time:", time, " // Estimated time remaining:", remaining), "")
    }
  j <- j+1
}

#----------------------------------#
# 4. Summarize & export results ####
#----------------------------------#

# Summarize results
summary(bd_soc$mean_hg_em)
hist(bd_soc$mean_hg_em, breaks=20, xlim=c(0, max(bd_soc$mean_hg_em)))
summary(bd_soc$mean_soil_hg)
hist(bd_soc$mean_soil_hg, breaks=20, xlim=c(0, max(bd_soc$mean_soil_hg)))

# Total Hg emissions (kg) & uncertainty from YK Delta 2015 wildfire
pixel_res <- 30 # grid cell resolution
cf_mg_to_kg <- 1/1000000 # mg to kg
sum(bd_soc$mean_hg_em)*(pixel_res^2)*cf_mg_to_kg
sum(bd_soc$ci_lwr_hg_em)*(pixel_res^2)*cf_mg_to_kg
sum(bd_soc$ci_upr_hg_em)*(pixel_res^2)*cf_mg_to_kg

# Export results for mapping in ArcGIS
#write.csv(bd_soc, paste0(root_path, "/", "Hg_emissions_uncertainty_", Sys.Date(), ".csv"), row.names=F)
