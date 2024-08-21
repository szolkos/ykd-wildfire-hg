#================================================================================================================================#
# Filename: 1_compile.R ####
# Author: Scott Zolkos (sgzolkos@gmail.com)
# Created: 2022-10-21
# Revised: 2024-08-20
# Background: Load and prepare data on soils (from unburned sites) and vegetation for subsequent analyses.
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
               here, # functionality | file referencing
               dplyr, # functionality
               reshape2, # functionality | data manipulation
               zoo, # functionality | time series structuring
               broom, # functionality | 'tidy' up prediction outputs
               janitor, # functionality | cleaning & examining data
               readxl, # functionality | read files in .xls or .xlsx format
               lubridate, # functionality | date manipulation
               foreign, # functionality | read .dbf files
               gdata, # functionality | bind columns with different number of rows
               Hmisc, # functionality
               onewaytests, # statistics | Welch's ANOVA for unequal variances
               ggplot2, # graphics | 'declaratively' creating graphics
               gridExtra, # graphics
               RColorBrewer, # graphics | plotting with R color brewer
               plotrix, # graphics | plotting, labeling, axis and color scaling functions, etc.
               ggpubr, # graphics | multi-panel ggplot (e.g. 'ggarrange' function,
               corrplot, # graphics | graphically display correlation matrices
               ecodist # graphics | ellipses in RDA
)

# Set location of root filepath directory
root_path <- here("data")

#---------------------------#
# 2. Load & prepare data ####
#---------------------------#

# Read in soil data from The Polaris Project (henceforth, 'Polaris')
polaris <- read.csv(paste0(root_path, "/", "polaris_data.csv"), header=T)
names(polaris)[names(polaris) == "vial_num"] <- "sample"
polaris <- droplevels(subset(polaris, select=c("last_name","site_id","sample","substrate","lat_dd","lon_dd","om_prop","c_prop","gwc_prop","bd_gcm3","depth_cm","c_kg_m2","burn_hx","burn_year","landcover","landcover_details")))

# Read in Hg data measured on Polaris soils for this study
hg_al <- read.csv(paste0(root_path, "/","thg_ykd_al.csv"), header=T)
hg_al <- droplevels(subset(hg_al, select=c("name","hg_ng","hg_ngg"), hg_al$name != "blank" & hg_al$name != "mess4"))
names(hg_al) <- c("sample","hg_ng","hg_ngg")

# Read in Polaris data on vegetation biomass in the YK-Delta (collected during 2017 & 2018 field expeditions)
## Note: 'spp_agb_gm2' = Summed dry weight of aboveground species or functional group, divided by harvest plot area (0.09 or 0.25 m2)
## Note: 'tot_agb_gm2' = Summed dry weight of aboveground vascular and nonvascular plants and lichen, divided by harvest plot area (0.09 or 0.25 m2)
veg <- read.csv(paste0(root_path, "/","ykd_veg_biomass.csv"), header=T) %>% 
  select(yr, site, sample_id, pft, spp_agb_gm2, tot_agb_gm2) %>%
  mutate(
    site = factor(site),
    sample_id = factor(sample_id),
    pft = factor(pft)
  )

# Read in Hg data measured on Polaris vegetation for this study
veg_hg <- read.csv(paste0(root_path, "/","ykd_veg_hgt.csv"), header=T) %>% 
  select(yr, site, plot, sample_id, pft, hg_ppb) %>%
  mutate(
    site = factor(site),
    plot = factor(plot),
    sample_id = factor(sample_id),
    pft = factor(pft)
  )

# Bind
ykd_hg <- merge(polaris, hg_al, by="sample", all.x=T)
ykd_hg$landcover <- as.factor(ykd_hg$landcover)

# Subset data then recombine, in order to omit n=9 samples from 1970s burn collected by Jardine
u_df <- droplevels(subset(ykd_hg, ykd_hg$burn_hx=="u"))
b_df <- droplevels(subset(ykd_hg, ykd_hg$burn_hx=="b" & ykd_hg$burn_year == "2015"))
ykd_hg <- rbind(u_df, b_df)
  
# Calculate Hg stores
# Hg (mg/m2) = SOC (kgC/m2) * RHgC (µgHg/gC) * cf1 (1,000 g/kg) ÷ cf2 (1,000 µg/mg)
## simplifies to: Hg (mg/m2) = SOC (kgC/m2) * RHgC (µgHg/gC)
depth_cm <- 30 # depth of soil core
ykd_hg$hg_mg_m2 <- ykd_hg$bd_gcm3 * ykd_hg$hg_ngg * depth_cm * 10000/1000000

# Omit lake sediments
ykd_hg <- droplevels(subset(ykd_hg, ykd_hg$substrate != "sediment"))
  
# Subset and export all data (not averaged) from unburned sites
ykd_hg_ub <- droplevels(subset(ykd_hg, ykd_hg$burn_hx == "u"))
#write.csv(ykd_hg_ub, paste0(getwd(), "/", "ykd_hg_ub_no_site_avg.csv"), row.names=F)
dim(droplevels(na.omit(subset(ykd_hg_ub, select=c("hg_mg_m2"))))) # for n unburned mmnts (n = 126)

# NEXT STEP DONE OUTSIDE OF R: In e.g. ArcGIS, determine which sites are within same 30 m landsat pixel (use WGS 1984 UTM Zone 3N)

# Create dataframe with average Hg & SOC values for samples collected within same 30 m Landsat pixel
df <- read.csv(paste0(root_path, "/", "ykd_hg_ub_site_locs_validated_2023_10_06.csv"), header=T) %>%
  group_by(siteid_sample) %>%
  summarise(lat_dd = mean(lat_dd),
            lon_dd = mean(lon_dd),
            depth_cm = mean(depth_cm),
            bd_gcm3 = mean(bd_gcm3),
            gwc_prop = mean(gwc_prop),
            c_prop = mean(c_prop),
            c_kg_m2 = mean(c_kg_m2),
            hg_ngg = mean(hg_ngg),
            hg_mg_m2 = mean(hg_mg_m2),
            ) %>%
  as.data.frame() %>%
  mutate(siteid_sample = factor(siteid_sample))
  
# RHgC calculations
# Calculate g C per sample: C (%) x bulk density (g/cm3) x core volume (cm3)
df$C_g <- df$c_prop * df$bd_gcm3 * df$depth_cm * (3.1415926*3^2)
# Calculate µg Hg per sample: Hg (ng/g) x bulk density (g/cm3) x core volume (cm3) ÷ 1000 (cf, ng=>µg)
df$Hg_ug <- df$hg_ngg * df$bd_gcm3 * df$depth_cm * (3.1415926*3^2) / 1000
# Calculate RHgC
df$rhgc_ugg <- df$Hg_ug/df$C_g

# Subset Lat & Lon & other relevant data, extract remote sensing data at those points using Google Earth Engine
soc_ub <- df %>% select(siteid_sample, lat_dd, lon_dd, c_kg_m2) %>% na.omit()
hg_ub <- df %>% select(siteid_sample, lat_dd, lon_dd, hg_mg_m2) %>% na.omit()
#write.csv(soc_ub, paste0(getwd(), "/", "soc_ub.csv"), row.names=F)
#write.csv(hg_ub, paste0(getwd(), "/", "hg_ub.csv"), row.names=F)