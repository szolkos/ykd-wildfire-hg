#===========================================================================================================#
# Filename: 6_peat_burned_area.R ####
# Author: Scott Zolkos (sgzolkos@gmail.com)
# Created: 2024-06-13
# Revised: 2024-08-20
# Background: Estimate northern peatland wildfire Hg emissions from 2002-2022.
# References
## Loboda et al. 2024: ABBA v2 (https://daac.ornl.gov/cgi-bin/dsviewer.pl?ds_id=2328)
#===========================================================================================================#

#-----------------------------------#
# 1. Prepare working environment ####
#-----------------------------------#

# Install packages, as needed
if("pacman" %in% installed.packages() == FALSE){install.packages("pacman")}

# Load necessary packages
pacman::p_load(tidyverse,
               here, # functionality | file referencing
               reshape2, # functionality | data manipulation
               zoo, # functionality | time series structuring
               broom, # functionality | 'tidy' up prediction outputs
               janitor, # functionality | cleaning & examining data
               readxl, # functionality | read files in .xls or .xlsx format
               lubridate, # functionality | date manipulation
               foreign, # functionality | read .dbf file
               purrr, # functionality | read multiple .dbf files
               MARSS, # statistics | Multivariate Auto-Regressive State Space Models
               ggplot2 # plotting
)

# Set location of root filepath directory
root_path <- here("data")

# Define conversion factors
cf_m2_to_ha <- 0.0001
cf_m2_to_km2 <- 0.000001

# Initialize an empty list to store results
results_list <- list()

years <- c(2002:2022)

#-----------------------------------------------------#
# 2. Compute annual peatland widlfire Hg emissions ####
#-----------------------------------------------------#

# Loop through each year from 2002 to 2022
for(year in years){
  # Construct file path for the current year to access burned area from ABBA v2 dataset (Loboda et al. 2024)
  file_path <- file.path(root_path, "burned_area", paste0("abba_", year, "_rpj_rtp_ba.dbf"))
  
  # Read data for the current year
  current_year_data <- read.dbf(file_path)
  
  # Calculate total annual burned area (km2) and Hg emissions (Mg)
  current_year_summary <- current_year_data %>%
    filter(peat >= 0) %>%
    mutate(
      ba_peat = (463.313^2) * cf_m2_to_km2 * (peat / 100),
      ba_nonpeat = (463.313^2) * cf_m2_to_km2 * (1 - (peat / 100)),
      ba_tot = ba_peat + ba_nonpeat,
      hg_em_med_peat = ba_peat * 1.36 / 1000, # divide by 1000 to convert kg to Mg
      hg_em_lo_peat = ba_peat * 0.80 / 1000,
      hg_em_hi_peat = ba_peat * 1.80 / 1000,
      hg_em_nonpeat = ba_nonpeat * 0.15 / 1000,
      hg_em_med_tot = hg_em_med_peat + hg_em_nonpeat,
      hg_em_lo_tot = hg_em_lo_peat + hg_em_nonpeat,
      hg_em_hi_tot = hg_em_hi_peat + hg_em_nonpeat
    ) %>%
    select(-grid_code) %>%
    summarise(across(
      c(
        ba_peat, ba_nonpeat, ba_tot,
        hg_em_lo_peat, hg_em_med_peat, hg_em_hi_peat, hg_em_nonpeat,
        hg_em_lo_tot, hg_em_med_tot, hg_em_hi_tot
      ), sum, na.rm = TRUE
    )) %>%
    mutate(year = year, dataset = "ABBA") %>%
    select(year, dataset, everything())
  
  # Append the current year's summary to the results list
  results_list[[year - 2001]] <- current_year_summary
}

# Combine annual summaries into a single dataframe
all_yrs <- bind_rows(results_list)
all_yrs$hg_em_peat_fxn <- all_yrs$hg_em_med_peat / all_yrs$hg_em_med_tot

# Print the final dataframe
print(all_yrs)

scale_factor <- max(all_yrs$hg_em_peat_fxn)/max(all_yrs$hg_em_med_peat)

# Plot
ggplot(all_yrs, aes(x = year)) +
  # Scales
  #scale_y_continuous(limits=c(0,140), breaks=seq(0,140,20)) +
  scale_y_continuous(limits=c(0,72), breaks=seq(0,70,10),
                     sec.axis=sec_axis(~.*scale_factor, breaks=seq(0,1,0.2), name="Peatland fraction of total Hg emissions")) +
  scale_x_continuous(limits=c(2000,2022), breaks=seq(2000,2020,5)) +
  # Themes
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "in"),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        panel.background=element_rect(fill='white'),
        panel.border=element_rect(colour="black", fill=NA, size=1),
        plot.background=element_rect(colour="white", size=1),
        #text=element_text(size=15),
        axis.title.y.left=element_text(size=15, margin=ggplot2::margin(t=0, r=10, b=0, l=0)),
        axis.title.y.right=element_text(size=15, margin=ggplot2::margin(t=0, r=0, b=0, l=10), colour="blue4"),
        axis.text.y.left=element_text(size=13.5, angle=0, hjust=0.5, colour="black"),
        axis.text.y.right=element_text(size=13.5, angle=0, hjust=0.5, colour="blue4"),
        #axis.text.y.right=element_text(color="blue4"),
        #axis.title.x=element_text(size=15, margin=ggplot2::margin(t=0, r=0, b=0, l=0)),
        axis.text.x=element_text(size=13.5, angle=0, hjust=0.5, colour="black"),
        legend.position="none") +
  # Total Hg emissions
  #geom_line(aes(y = hg_em_med_tot), color = "gray20") +
  #geom_ribbon(aes(ymin = hg_em_lo_tot, ymax = hg_em_hi_tot), fill = "gray50", alpha = 0.2) +
  # Peatland Hg emissions
  geom_line(aes(y = hg_em_med_peat), color = "black", lwd=0.7) +
  geom_ribbon(aes(ymin = hg_em_lo_peat, ymax = hg_em_hi_peat), fill = "gray50", alpha = 0.3) +
  # Peatland Hg fraction
  geom_line(aes(y = hg_em_peat_fxn/scale_factor), color="blue4", lty=2, lwd=0.7) +
  # Axis labels
  labs(x = NULL, y=expression(Peatland~Hg~emissions~(Mg~y^-1)))

#--------------------------#
# 3. Summary statistics ####
#--------------------------#
# Mean peatland annual northern Hg emissions
all_yrs %>%
  summarise(
    mean_hg_em_lo_peat = mean(hg_em_lo_peat, na.rm = TRUE),
    mean_hg_em_med_peat = mean(hg_em_med_peat, na.rm = TRUE),
    mean_hg_em_hi_peat = mean(hg_em_hi_peat, na.rm = TRUE)
  ) %>%
  print()

# Mean total annual northern Hg emissions
all_yrs %>%
  summarise(
    mean_hg_em_lo_tot = mean(hg_em_lo_tot, na.rm = TRUE),
    mean_hg_em_med_tot = mean(hg_em_med_tot, na.rm = TRUE),
    mean_hg_em_hi_tot = mean(hg_em_hi_tot, na.rm = TRUE)
  ) %>%
  print()

# Peatland fraction
all_yrs %>%
  summarise(mean_hg_em_peat_fxn = mean(hg_em_peat_fxn, na.rm = TRUE)) %>%
  print()