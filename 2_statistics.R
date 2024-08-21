#================================================================================================================================#
# Filename: 2_statistics.R ####
# Author: Scott Zolkos (sgzolkos@gmail.com)
# Created: 2023-08-28
# Revised: 2024-08-20
# Background: Compute statistics reported in the manuscript.
# References
## None
#================================================================================================================================#

#-----------------------------------------------------------#
# 1. Summary stats of soils data for Sec. 3.1 & Table S5 ####
#-----------------------------------------------------------#
# Fibric layer depth (mean ± SD)
length(na.omit(df$depth_cm))
mean(na.omit(df$depth_cm)); sd(na.omit(df$depth_cm))
# Bulk density (mean ± SD)
length(na.omit(df$bd_gcm3))
mean(na.omit(df$bd_gcm3)); sd(na.omit(df$bd_gcm3))
# OC % (mean ± SD)
length(na.omit(df$c_prop))
mean(na.omit(df$c_prop))*100; sd(na.omit(df$c_prop))*100
# OC stores (mean ± SD)
length(na.omit(df$c_kg_m2))
mean(na.omit(df$c_kg_m2)); sd(na.omit(df$c_kg_m2))
# Hg (mean ± SD)
length(na.omit(df$hg_ngg))
mean(na.omit(df$hg_ngg)); sd(na.omit(df$hg_ngg))
# Hg stores (mean ± SD)
length(na.omit(df$hg_mg_m2))
mean(na.omit(df$hg_mg_m2))*1000; sd(na.omit(df$hg_mg_m2))*1000

#------------------------------------------------------#
# 2-5. Summary stats of veg data for Tables S3 & S5 ####
#--------------------------------------------------------------------------------------#
# 2. Calculate mean aboveground biomass (AGB, g/m2) per plant functional type (PFT) ####
#--------------------------------------------------------------------------------------#
# Summed biomass by PFT at each site
veg_sum <- veg %>% 
  group_by(site, pft) %>% 
  dplyr::summarise(tot_agb_gm2 = sum(spp_agb_gm2)) %>% 
  as.data.frame()
veg_sum$pft <- reorder(veg_sum$pft, -veg_sum$tot_agb_gm2, FUN=median)

# Total AGB by site
site_agb <- veg_sum %>% group_by(site) %>% dplyr::summarise(site_agb_gm2 = sum(tot_agb_gm2)) %>% as.data.frame()
mean(site_agb$site_agb_gm2); sd(site_agb$site_agb_gm2)

# Median and SD of summed total AGB by PFT, for all sites
veg_agb_smry <- veg_sum %>% 
  group_by(pft) %>% 
  dplyr::summarise(agb_mean = mean(tot_agb_gm2), agb_sd = sd(tot_agb_gm2)) %>% 
  arrange(desc(agb_mean)) %>% as.data.frame()

#--------------------------------------#
# 3. Calculate mean Hg conc per PFT ####
#--------------------------------------#
# Median and SD of summed total AGB by PFT, for all sites ####
veg_hg_smry <- veg_hg %>% 
  group_by(pft) %>% 
  dplyr::summarise(hg_mean = mean(hg_ppb), hg_sd = sd(hg_ppb)) %>% 
  arrange(desc(hg_mean)) %>% 
  as.data.frame()

#------------------------------------------------#
# 4. Calculate mean Hg stores (µg/m2) per PFT ####
#------------------------------------------------#
# Mean biomass Hg by PFT at each site
veg_hg_mean <- veg_hg %>% 
  group_by(site, pft) %>% 
  dplyr::summarise(hg_mean = mean(hg_ppb)) %>% 
  as.data.frame()
veg_hg_mean$pft <- reorder(veg_hg_mean$pft, -veg_hg_mean$hg_ppb, FUN=median)

# (1) Calculate mean hg_ppb by site & PFT
veg_hg_site_mean <- veg_hg_mean %>% 
  group_by(site, pft) %>% 
  dplyr::summarise(hg_mean = mean(hg_ppb)) %>% 
  as.data.frame()

# (2) Merge (1) with total veg AGB by site & PFT (NAs result where no [Hg] data for a given site & PFT)
veg_hg_site_agb <- left_join(veg_sum, veg_hg_site_mean, by=c("site","pft"))

# (3) Calculate Hg stores for (2)
veg_hg_site_agb$veg_hg_ugm2 <- veg_hg_site_agb$tot_agb_gm2 * veg_hg_site_agb$hg_mean / 1000

# (4) Calculate mean & SD for Hg stores by site & PFT
veg_hg_smry2 <- veg_hg_site_agb %>% 
  select(site, pft, veg_hg_ugm2) %>% 
  na.omit() %>%
  group_by(pft) %>% 
  dplyr::summarise(hg_stores_mean = mean(veg_hg_ugm2), hg_stores_sd = sd(veg_hg_ugm2)) %>% 
  arrange(desc(hg_stores_mean)) %>% 
  as.data.frame()

#------------------------------------------#
# 5. Vegetation data for Tables S3 & S5 ####
#------------------------------------------#
tables_s3_s5_veg <- left_join(veg_agb_smry, veg_hg_smry, by="pft") %>%
  left_join(., veg_hg_smry2, by="pft")

print(tables_s3_s5_veg)
#------------------------------------------#

#---------------------------------------------------------#
# 6. Hg stores in soils vs. veg (reported in Sec. 3.1) ####
#---------------------------------------------------------#
# Data
soil_hg_stores <- as.numeric(na.omit(df$hg_mg_m2*1000)) # n = 50 (length(na.omit(df$hg_mg_m2)))
veg_hg_stores <- veg_hg_site_agb %>% 
  select(site, pft, veg_hg_ugm2) %>% 
  na.omit() %>%
  group_by(site) %>% 
  dplyr::summarise(hg_stores_mean = sum(veg_hg_ugm2), hg_stores_sd = sd(veg_hg_ugm2)) %>% 
  arrange(desc(hg_stores_mean)) %>% 
  as.data.frame() %>%
  select(hg_stores_mean)
# Run a Bartlett's Test to determine if variances are homogenous (p < 0.05 = not homogenous)
bartlett.test(list(soil_hg_stores, veg_hg_stores))
# Run Welch's t-test with unequal variances
t.test(soil_hg_stores, veg_hg_stores, var.equal=F)

#----------------------------------------#
# 7. Stats reported in Sec. 3.2 & 3.3 ####
#----------------------------------------#
# Burn depth ####
bd <- na.omit(read.csv(paste0(root_path, "/", "moubarak_ykd_2015_burn_depth_2023_01_24.csv"), header=T)) %>%
  mutate(site = factor(site))
length(na.omit(bd$burn_depth))
mean(bd$burn_depth); sd(bd$burn_depth)

# RHgC ####
length(na.omit(df$rhgc_ugg))
mean(na.omit(df$rhgc_ugg)); sd(na.omit(df$rhgc_ugg))

# Hg emissions ####
hgem_df <- read.dbf(paste0(root_path, "/", "Hg_emissions_uncertainty_2023_10_30_no_wtr", ".dbf"))

# Rename columns
names(hgem_df) <- c("num", "lat_dd", "lon_dd", "bd", "soc", 
                    "mean_hg_em", "median_hg_em", "sd_hg_em", "ci_lwr_hg_em", "ci_upr_hg_em", 
                    "mean_soil_hg", "median_soil_hg", "sd_soil_hg", "ci_lwr_soil_hg", "ci_upr_soil_hg")

# 5th and 95th percentiles for entire study area
# Burn depth (cm)
hist(hgem_df$bd, breaks=40)
mean(hgem_df$bd); quantile(hgem_df$bd, 0.05); quantile(hgem_df$bd, 0.95)
# SOC stores (kgC/m2)
hist(hgem_df$soc, breaks=40)
mean(hgem_df$soc); quantile(hgem_df$soc, 0.05); quantile(hgem_df$soc, 0.95)
# Hg stores (mgHg/m2)
hist(hgem_df$mean_soil_hg, breaks=40)
mean(hgem_df$mean_soil_hg); quantile(hgem_df$mean_soil_hg, 0.05); quantile(hgem_df$mean_soil_hg, 0.95)
# Hg emissions (mg/m2)
hist(hgem_df$mean_hg_em, breaks=40, xlab=expression("Hg emissions (mg m"^"-2"*")"))
mean(hgem_df$mean_hg_em); quantile(hgem_df$mean_hg_em, 0.05); quantile(hgem_df$mean_hg_em, 0.95)
# Total Hg emissions (kg) & uncertainty from YK Delta 2015 wildfire
pixel_res <- 30 # grid cell resolution
cf_mg_to_kg <- 1/1000000 # mg to kg
sum(hgem_df$mean_hg_em)*(pixel_res^2)*cf
sum(hgem_df$ci_lwr_hg_em)*(pixel_res^2)*cf
sum(hgem_df$ci_upr_hg_em)*(pixel_res^2)*cf



#------------------------------------#
# Other stats from raw field data ####
#------------------------------------#

# Mean core depth ± SD
length(na.omit(ykd_hg$depth_cm))
mean(na.omit(ykd_hg$depth_cm)); sd(na.omit(ykd_hg$depth_cm))

# Bulk density
length(na.omit(df$bd_gcm3))
mean(na.omit(df$bd_gcm3)); sd(na.omit(df$bd_gcm3))

# Mean C content
length(na.omit(ykd_hg$c_prop))
mean(na.omit(ykd_hg$c_prop)); sd(na.omit(ykd_hg$c_prop))