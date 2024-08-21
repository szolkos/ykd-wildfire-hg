#================================================================================================================================#
# Filename: 3_burn_depth_rf_mdl.R ####
# Author: Scott Zolkos (sgzolkos@gmail.com) adapted code written by Stefano Potter (Potter et al. 2022)
# Created: 2023-01-06
# Revised: 2024-08-20
# Background: Algorithm to predict YK Delta wildfire burn depth from field measurements of burn depth (Moubarak et al. 2022) and 
## various remote sensing metrics of fire severity and terrain qualities.
# References
## Moubarak et al. 2022 (https://bg.copernicus.org/preprints/bg-2022-144/)
## Potter et al. 2022 (https://doi.org/10.5194/egusphere-2022-364)
#================================================================================================================================#

#-----------------------------------#
# 1. Prepare working environment ####
#-----------------------------------#

# Install R package management tool (pacman), as needed
if("pacman" %in% installed.packages() == FALSE){install.packages("pacman")}

# Load necessary packages
pacman::p_load(tidyverse, # functionality
               caret,
               doSNOW,
               ranger,
               parallel, # for corrplot
               psych, # for corrplot
               pdp, # for partial dependency plots
               gridExtra,
               ggpubr
)

# Set location of root filepath directory
root_path <- here("data")

# Create the outpath in case it doesn't exist
run_no <- 1
met <- "LOOCV"
out_path <- paste0(root_path, "/Machine learning/Burn depth/", Sys.Date(), " ", met, "_", run_no, "/")
dir.create(out_path, recursive = F)

# Set up cluster - requires doSNOW, and parallel packages, allows to use multiple cpus
cores <- detectCores() - 2 #if you want to use your computer for things other than the model, subtract 1 or 2.
cl <-makeCluster(cores)
registerDoSNOW(cl)

#---------------------------#
# 2. Load & prepare data ####
#---------------------------#

# Set a random seed, for reproducibility 
seed <- 57
set.seed(seed)

# Prepare data
preds_bd <- read.csv(paste0(root_path, "/", "burn_depth_preds_2023_01_20.csv"), header=T)
preds_bd <- merge(x=bd, y=preds_bd, by="site", all.x=T) %>% 
  select(-c(bd_mean, lat_dd, lon_dd)) %>% # "bd_mean" = mean bulk density (not burn depth!)
  dplyr::rename(
    bd = burn_depth, # rename "burn_depth" column to "bd"
    id = site
    ) %>%
  select(id, bd, everything()) %>% # Move response variable into first column of df
  droplevels() %>%
  na.omit()

#---------------------------------------------------#
# 3. Inspect data (e.g., covariate correlations) ####
#---------------------------------------------------#

# Notes on transforming data for tree-based models
## Normalizing values to similar scale can improve model performance and training stability) | https://developers.google.com/machine-learning/data-prep/transform/normalization
## However, while normalizing is not critical for tree-based models, it is important for distance-based (e.g. k-nearest-neighbors & k-means are Euclidean distance) | https://www.kdnuggets.com/2022/07/random-forest-algorithm-need-normalization.html
## Don't center & scale, either | https://stackoverflow.com/questions/8961586/do-i-need-to-normalize-or-scale-data-for-randomforest-r-package
# Note on covariate correlations for machine learning methods
## Especially for boosting methods, should reduce data to relatively uncorrelated variables where correlation within [-0.7 < r < 0.7] | https://stats.stackexchange.com/questions/141619/wont-highly-correlated-variables-in-random-forest-distort-accuracy-and-feature

# Investigate response variable, transform as needed
hist(preds_bd$bd, breaks=10); hist(log(preds_bd$bd), breaks=10)

# Evaluate correlations between variables
cor_df <- preds_bd %>% select(-c(id))
pval <- corr.test(cor_df, adjust="none", method="pearson")$p # "pearson" assesses linear relationship, uses raw data; "spearman" can also assess monotonic (non-linear), using rank-ordered variables
pairs.panels(cor_df,
             smooth = T, # If TRUE, draws loess smooths
             scale = F, # If TRUE, scales the correlation text font
             density = F, # If TRUE, adds density plots and histograms
             ellipses = F, # If TRUE, draws ellipses
             method = "pearson", # Correlation method (also "spearman" or "kendall")
             pch = 21,
             cex=0.4,
             lm = F, # If TRUE, plots linear fit rather than the LOESS (smoothed) fit
             cor = T, # If TRUE, reports correlations
             jiggle = F, # If TRUE, data points are jittered
             factor = 2, # Jittering factor
             hist.col = 4, # Histograms color
             stars = T, # If TRUE, adds significance level with stars
             ci = T,
             cex.cor=3)

# Subset non-correlated predictors
preds_bd <- preds_bd %>% 
  select(c(id, bd, dnbr, tc_wet, slope, rugosity)) # NDWI & TC wetness are functionally similar, so former omitted
  
# Shuffle predictors
rf_bd <- preds_bd[sample(1:nrow(preds_bd)), ] 

# Get unique ids for each row
ids <- rf_bd$id

# Store for modeling as model dataframe 'rf_bd' (Random Forests burn depth)
rf_bd <- rf_bd %>% 
  select(-id) %>% # Don't model with id
  dplyr::rename(y = bd) # Rename the first column (y) to y

#------------------------------------#
# 4. Develop Random Forests model ####
#------------------------------------#

# Look up parameters which can be tuned for ranger Random Forests model
modelLookup('ranger')

# Set method, either "LOOCV" or "repeatedcv"; LOOCV generally better for small datasets (n ~ ≤ 100)
## More info on LOOCV: https://topepo.github.io/caret/recursive-feature-elimination.html; https://machinelearningmastery.com/loocv-for-evaluating-machine-learning-algorithms/
## More info on k-fold cross validation: https://machinelearningmastery.com/k-fold-cross-validation/
met <- "LOOCV"

#______________________#
# Feature selection ####
#______________________#

# Set parameters for recursive feature elimination (RFE)
set.seed(seed)
control <- rfeControl(functions = rfFuncs,
                      method = met,
                      number = 100, # Set if 'met' = 'LOOCV', aim for 1000, depending on compute limitations
                      verbose = FALSE,
                      allowParallel = T)

# Name of target (i.e. 'y' variable you are predicting)
outcomeName <- names(rf_bd)[1]

# Names of predictors
predictors <- names(rf_bd)[!names(rf_bd) %in% outcomeName]

# RFE. 'sizes' is how many variables at once you will check (# of features that will be retained)
## More info: http://topepo.github.io/caret/recursive-feature-elimination.html, Sec. 20.5 = Helper Functions
pred_importance <- rfe(as.matrix(rf_bd[,predictors]),
                       as.matrix(rf_bd[,outcomeName]),
                       metric = "RMSE", # for regression, use "RMSE" or "Rsquared" as metric for optimal model selection
                       maximize = FALSE, # FALSE = model w/ minimum value of defined metric selected as optimal model (use for RMSE; for R2, use TRUE)
                       rfeControl = control,
                       sizes = c(1:length(rf_bd[,predictors])))
  
# Save selected variables to a csv
write_csv(as_tibble(pred_importance$results), file.path(out_path, "rfe.csv"))
    
# Plot predictor importance
# Store results
pi <- pred_importance$results
# set graphical parameters
miny <- min(pi$RMSE)*0.4
maxy <- max(pi$RMSE)*1.2
scale_factor <- max(pi$Rsquared)/max(pi$RMSE)
y_lims <- round(c(miny, maxy),1)
y_breaks <- unique(round(seq(0,maxy,maxy*0.1),1))
sec_y_breaks <- seq(min(pi$Rsquared-0.1), max(pi$Rsquared+0.1), 0.1)
ylab <- "RMSE (LOOCV)"
y_lims <- c(0,4)
y_breaks <- seq(0,4,1)
      
# Plot
fs_plot <- ggplot(data=pi, aes(y=RMSE, x=Variables)) +
  # Scales
    scale_y_continuous(limits=y_lims, breaks=y_breaks, sec.axis=sec_axis(~.*scale_factor, name=expression(italic(R)^2))) +
    scale_x_continuous(position="bottom", limits=c(0.8,max(pi$Variables)+max(pi$Variables)*0.05), breaks=seq(0,max(pi$Variables),1)) +
  # Themes
    theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "in"),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          panel.background=element_rect(fill='white'),
          panel.border=element_rect(colour="black", fill=NA, linewidth=1),
          plot.background=element_rect(colour="white", linewidth=1),
          text=element_text(size=18),
          axis.text.y=element_text(angle=0, hjust=0.5, colour="black"), axis.text.y.right=element_text(color="blue"),
          axis.title.y.left=element_text(margin=unit(c(t=0, r=10, b=0, l=0), "pt")),
          axis.title.y.right=element_text(margin=unit(c(t=0, r=0, b=0, l=10), "pt"), color="blue"),
          axis.text.x=element_text(angle=0, hjust=0.5, colour="black"),
          axis.title.x=element_text(margin=unit(c(t=10, r=0, b=0, l=0), "pt")),
          legend.position="none") +
  # Data
    # R2
      geom_line(data=pi, aes(x=Variables, y=Rsquared/scale_factor), color="blue", size=1) +
        #geom_errorbar(data=pi, aes(ymin=(Rsquared-(RsquaredSD/2))/scale_factor, ymax=(Rsquared+(RsquaredSD/2))/scale_factor), color="blue", position=position_dodge()) +
      geom_point(data=pi, aes(y=Rsquared/scale_factor, x=Variables), pch=21, fill="blue", col="blue", size=3, stroke=1, alpha=1) +
    # RMSE
      #geom_path(color="black", size=1) +
      geom_line(data=pi, aes(x=Variables, y=RMSE), color="black", size=1) +
      #geom_errorbar(data=pi, aes(ymin=RMSE-(RMSESD/2), ymax=RMSE+(RMSESD/2)), color="black", position=position_dodge()) + # 
      geom_point(data=pi, aes(y=RMSE, x=Variables), pch=21, fill="black", col="black", size=3, stroke=1, alpha=1) +
  # Axis labels
    labs(y=ylab, x=expression(italic(n)~variables))

# Print
print(fs_plot)

# Export
ggsave(fs_plot, filename = file.path(out_path, paste0("feature selection.png")), device = 'png', dpi = 150, width = 7, height = 5)

#________________________#
# Variable importance ####
#________________________#

# Compile all variables
pred_imp <- pred_importance$optVariables %>%
  as.data.frame() %>%
  dplyr::rename(Var = ".") %>%
  rownames_to_column(., "VarOrder") %>%
  mutate(., VarOrder = as.numeric(VarOrder), Var = as.character(Var))

# Subset the optimal variables, using minimum RMSE
opt <- pred_importance$optVariables
opt[[length(opt) + 1]] <- "y" #add in the y
df_opt <- rf_bd %>% dplyr::select(opt) #select the good variables
df_opt <- df_opt %>% dplyr::select(y, everything()) #put y as first column again
        
# Use random forests to get feature importance 
control <- trainControl(method = met, number = 100)
mod <- train(y~., data = df_opt, method = 'rf', trControl = control, importance = T)

# Store & plot variable importance
importance <- varImp(mod, scale = FALSE)
plot(importance)

# Export plot
gsave(imp_plot, filename = file.path(out_path, paste0("variable importance.png")), device = 'png', dpi = 150, width = 6, height = 4)

# Make a table of importance
importance <- as.matrix(varImp(mod)$importance) %>%
  as.data.frame() %>%
  rownames_to_column(., "Variables")

# Export as csv
write_csv(importance, file.path(out_path, 'rfe_importance.csv'))

#__________________#
# Develop model ####
#__________________#

# First, using optimal variables (determined above), search for optimal model parameters
# Here, RMSE is used to select the the best parameters for each model
# See https://topepo.github.io/caret/model-training-and-tuning.html#basic-parameter-tuning
set.seed(seed)
fitControl <- trainControl(method = met,
                           number = 100,
                           savePredictions = 'final',
                           index = createResample(df_opt$y, 25),
                           allowParallel = TRUE)
 
# Next, get the best model parameters with a random grid search of length 10
# i.e., 10 random numbers selected for any given model parameter and tested for each parameter combination
ranger_bd <- train(y ~., data = df_opt, method = 'ranger', tuneLength = 10, metric = 'RMSE', trControl = fitControl)

# Lastly, get final model fits to estimate model performance
set.seed(seed)
fitControl <- trainControl(method = met,
                           number = 100,
                           savePredictions = 'final',  
                           allowParallel = T)

# Use the best parameters for each model to see how well it is performing (assess model below)
mdl_eval <- train(y ~., 
                            data = df_opt, 
                            method = 'ranger', 
                            tuneGrid=data.frame(.mtry = ranger_bd$bestTune$mtry, # the # of randomly selected predictors (Kuhn et al. 2022)
                                                .splitrule = ranger_bd$bestTune$splitrule, # splitting rule, default = 'variance' (Wright et al. 2022)
                                                .min.node.size = ranger_bd$bestTune$min.node.size # min n observations in each node. implicitly controls tree depth, b/c larger n obs in each node => smaller trees | https://stats.stackexchange.com/questions/158583/what-does-node-size-refer-to-in-the-random-forest (Wright et al. 2022)
                                                ), trControl = fitControl)

#_________________#
# Assess model ####
#_________________#

# Inspect variables stores in model
ls(mdl_eval)

# Get R^2
mdl_eval$results$Rsquared

# Observed vs predicted for calibration data
compare <- as_tibble(mdl_eval$pred)
resids_obs <- as.numeric(tapply(compare$obs, compare$rowIndex, FUN=mean, na.rm=T))
resids_preds <- as.numeric(tapply(compare$pred, compare$rowIndex, FUN=mean, na.rm=T))
resids_obs_preds <- tibble(obs = resids_obs, pred = resids_preds)
compare <- resids_obs_preds %>% mutate(Resids = obs - pred)

# Plot calibration data residuals
resids <- ggplot(compare, aes(x=pred, y=Resids))  +
  theme_bw() +
  theme(text=element_text(size=18)) +
  geom_point(pch=21, size=6, col="white", fill="black", stroke=1) +
  labs(x = "Predicted", y = "Residual")
print(resids)

# Export
final_model <- "ranger"
ggsave(resids, filename = file.path(out_path,  paste0(final_model, '_residuals.png')), device = 'png', dpi = 150, width = 10, height = 10)
      
# Store final model for assessment
fin_mdl <- train(y ~., 
                 data = df_opt, 
                 method = 'ranger', 
                 tuneGrid=data.frame(.mtry = ranger_bd$bestTune$mtry, 
                                     .splitrule = ranger_bd$bestTune$splitrule, 
                                     .min.node.size = ranger_bd$bestTune$min.node.size
                 ))
      
# If desired, export data in native R format (Rds is used to save a single R object)
saveRDS(fin_mdl, file=file.path(out_path, paste0("full_model_", final_model ,".rds")))

# Get full predictions and Rsq for desired model
full_pred <- predict(fin_mdl, newdata=df_opt %>% select(-y))
full_rsq = rsq(df_opt$y, full_pred)

# Obs & preds (calibration data)
obs_preds <- tibble(Obs = df_opt$y, Pred = full_pred) #, id = ids)

# Export observations & predictions, if desired (e.g., for mapping in GEE)
#write_csv(obs_preds, file.path(out_path, paste0(final_model, '_full_model_ob_pred.csv')))

    
# Predicted vs. observed (full dataset) ####

# Summary stats
obs_preds_lm <- summary(lm(obs_preds$Pred ~ obs_preds$Obs))
#ls(obs_preds_lm)
slope_val <- round(coef(obs_preds_lm)[2],2)
int_val <- round(coef(obs_preds_lm)[1],2)
tval <- round(summary(lm(obs_preds$Pred ~ obs_preds$Obs))$coefficients[6],2)
fstat_val <- round(summary(lm(obs_preds$Pred ~ obs_preds$Obs))$fstatistic[1],1)
dfval <- round(summary(lm(obs_preds$Pred ~ obs_preds$Obs))$df[2],0)
pval <- round(summary(lm(obs_preds$Pred ~ obs_preds$Obs))$coefficients[8],3)
if(pval < 0.01){pval <- "< 0.01"}
rsq_val <- round(summary(lm(obs_preds$Pred ~ obs_preds$Obs))$r.squared,2)
      
# Plot
obs_preds_plot <- ggplot(obs_preds, aes(x = Obs, y =  Pred)) + 
  # Scales
    scale_x_continuous(limits=c(0,20)) +
    scale_y_continuous(limits=c(0,20)) +
  # Themes
    theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "in"),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          panel.background=element_rect(fill='white'),
          panel.border=element_rect(colour="black", fill=NA, size=1),
          plot.background=element_rect(colour="white", size=1),
          plot.title=element_text(size=10), #
          text=element_text(size=18),
          axis.title.y.left=element_text(margin=unit(c(t=0, r=10, b=0, l=0), "pt")),
          axis.text.y=element_text(angle=0, hjust=0.5, colour="black"),
          axis.title.x=element_text(margin=unit(c(t=10, r=0, b=0, l=0), "pt")),
          axis.text.x=element_text(angle=0, hjust=0.5, colour="black"),
          legend.position="none") +
  # Reference lines
    geom_abline(intercept=0, slope=1, color='blue', linetype='solid', size=0.6) +
    geom_abline(intercept=int_val, slope=slope_val, color='black', linetype='dashed', size=0.6) +
  # Data
    geom_point(pch=21, size=6, col="white", fill="black", stroke=1) +
  # Labels
    ggtitle(paste0("Model = ", final_model, "\n", paste0("F = ", fstat_val, ", df = ", dfval, ", p = ", pval, ", R2 = ", rsq_val), "\n", paste0("pred = ", slope_val, "*obs + ", int_val))) +
    labs(x = "Observed", y = "Predicted")

# Obs & preds (LOOCV)
loocv_obs <- as.numeric(tapply(mdl_eval[[5]]$obs, mdl_eval[[5]]$rowIndex, FUN=mean, na.rm=T))
loocv_preds <- as.numeric(tapply(mdl_eval[[5]]$pred, mdl_eval[[5]]$rowIndex, FUN=mean, na.rm=T))
obs_preds_val <- tibble(Obs = loocv_obs, Pred = loocv_preds)


# Predicted vs. observed (LOOCV) ####

# Summary stats
obs_preds_val_lm <- summary(lm(obs_preds_val$Pred ~ obs_preds_val$Obs))
#ls(obs_preds_val_lm)
slope_val <- round(coef(obs_preds_val_lm)[2],2)
int_val <- round(coef(obs_preds_val_lm)[1],2)
tval <- round(summary(lm(obs_preds_val$Pred ~ obs_preds_val$Obs))$coefficients[6],2)
fstat_val <- round(summary(lm(obs_preds_val$Pred ~ obs_preds_val$Obs))$fstatistic[1],1)
dfval <- round(summary(lm(obs_preds_val$Pred ~ obs_preds_val$Obs))$df[2],0)
pval <- round(summary(lm(obs_preds_val$Pred ~ obs_preds_val$Obs))$coefficients[8],3)
if(pval < 0.01){pval <- "< 0.01"}
rsq_val <- round(summary(lm(obs_preds_val$Pred ~ obs_preds_val$Obs))$r.squared,2)
          
# Plot
obs_preds_val_plot <- ggplot(obs_preds_val, aes(x = Obs, y =  Pred)) + 
  # Scales
    #scale_x_continuous(limits=c(floor(min(obs_preds_val$Obs)), ceiling(max(obs_preds_val$Obs)))) +
    #scale_y_continuous(limits=c(floor(min(obs_preds_val$Obs)), ceiling(max(obs_preds_val$Obs)))) +
    scale_x_continuous(limits=c(0,20)) +
    scale_y_continuous(limits=c(0,20)) +
  # Themes
    theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "in"),
          panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          panel.background=element_rect(fill='white'),
          panel.border=element_rect(colour="black", fill=NA, size=1),
          plot.background=element_rect(colour="white", size=1),
          plot.title=element_text(size=10), #
          text=element_text(size=18),
          axis.title.y.left=element_text(margin=unit(c(t=0, r=10, b=0, l=0), "pt")),
          axis.text.y=element_text(angle=0, hjust=0.5, colour="black"),
          axis.title.x=element_text(margin=unit(c(t=10, r=0, b=0, l=0), "pt")),
          axis.text.x=element_text(angle=0, hjust=0.5, colour="black"),
          legend.position="none") +
  # Reference lines
    geom_abline(intercept=0, slope=1, color='blue', linetype='solid', size=0.6) +
    geom_abline(intercept=int_val, slope=slope_val, color='black', linetype='dashed', size=0.6) +
  # Data
    geom_point(pch=21, size=6, col="white", fill="black", stroke=1) +
  # Labels
    ggtitle(paste0("Model = ", final_model, "\n", paste0("F = ",fstat_val, ", df = ", dfval, ", p = ", pval, ", R2 = ", rsq_val), "\n", paste0("pred = ", slope_val, "*obs + ", int_val))) +
    labs(x = "Observed", y = "Predicted")

cal_val_plots <- ggarrange(obs_preds_plot, obs_preds_val_plot)
print(cal_val_plots)

# Export full model and LOOCV cal-val plots
ggsave(cal_val_plots, filename = file.path(out_path, paste0(" cal-val plots ", final_model, ".png")), device = 'png', dpi = 150, width = 10, height = 5)

#---------------------------#
# 5. Predict to new data ####
#---------------------------#
  
#________________________#
# Read & inspect data ####
#________________________#
preds_gee_bd <- read.csv(paste0(root_path, "/", "burn_depth_preds_gee_2023_01_22.csv"), header=T)
options(digits = 20) # Modify global R options for # sig digits
head(preds_gee_bd)

# Clean data
preds_gee_bd <- preds_gee_bd %>%
  # Rename columns
  dplyr::rename(
    pixel = system.index,
    dnbr = dNBR,
    rugosity = elevation_stdDev,
    slope = slope,
    tc_wet = wetness,
    geo = .geo
  ) %>%
  # Extract lat and lon
  separate(., geo, c("text", "lon", "lat"), sep=",") %>%
  separate(., lon, c("text2", "lon"), sep="-") %>%
  separate(., lat, c("lat", "text3"), sep="]}") %>%
  mutate(.,
    lon_dd = as.numeric(lon)*-1,
    lat_dd = as.numeric(lat)
  ) %>% 
  select(-c(text, text2, text3)) %>%
  mutate(num = as.numeric(row_number())) %>%
  select(num, everything())

#________________________#
# Predict to new data ####
#________________________#
# Omit NAs
preds_gee_bd2 <- preds_gee_bd %>% na.omit()

# Predict to new data
preds_gee_bd2$burn_depth_pred <- predict(fin_mdl, newdata=preds_gee_bd2 %>% select(-c(num,pixel, lat_dd, lon_dd)))

# Merge predicted values to original df with all pixel ids
preds_gee_bd3 <- merge(x=preds_gee_bd, y=preds_gee_bd2 %>% select(num, burn_depth_pred), by="num", all.x=T) %>%
  arrange(., num) %>%
  select(., c(num, lat_dd, lon_dd, dnbr, tc_wet, slope, rugosity, burn_depth_pred)) %>%
  na.omit()

# Inspect data
head(preds_gee_bd3)
hist(preds_gee_bd3$burn_depth_pred, breaks=20)
mean(preds_gee_bd3$burn_depth_pred)
sd(preds_gee_bd3$burn_depth_pred)
summary(preds_gee_bd3)

# Export
#write.csv(preds_gee_bd3, paste0(root_path, "/", "burn_depth_for_gee_", Sys.Date(), ".csv"), row.names=F)

# Reset global R options for # sig digits to default (5)
options(digits = 5)
        
# Stop the cluster
stopCluster(cl)

# Clear data from active memory
rm(preds_gee_bd); rm(preds_gee_bd2)

#--------------------------------------------------------#
# 6. Variable importance and partial dependence plots ####
#--------------------------------------------------------#
# Load best model
fin_mdl_burn_depth <- readRDS(file.path(paste0(root_path, "/Machine learning/Burn depth/2023-01-24 LOOCV_6/full_model_ranger.rds")))
summary(fin_mdl_burn_depth)
fin_mdl_burn_depth$finalModel
fin_mdl_burn_depth$bestTune
fin_mdl_burn_depth$results
fin_mdl_burn_depth$pred
fin_mdl_burn_depth$trainingData

#________________________#
# Variable importance ####
#________________________#
# Extract variable importance
imp_burn_depth <- varImp(fin_mdl_burn_depth, scale = FALSE)

# Convert the importance to a data frame
imp_burn_depth_df <- as.data.frame(imp_burn_depth$importance)

# Sum the importance values
imp_burn_depth_tot <- sum(imp_burn_depth_df$Overall)

# Normalize the importance values
imp_burn_depth_df$Normalized <- (imp_burn_depth_df$Overall / imp_burn_depth_tot) * 100

#_____________________________#
# Partial dependence plots ####
#_____________________________#
# Extract data from the trellis object
pdp_burn_depth_rugosity <- partial(fin_mdl_burn_depth, pred.var="rugosity", plot=TRUE, chull=TRUE, train=fin_mdl_burn_depth$trainingData)$panel.args[[1]] %>% 
  as.data.frame()  
pdp_burn_depth_dnbr <- partial(fin_mdl_burn_depth, pred.var="dnbr", plot=TRUE, chull=TRUE, train=fin_mdl_burn_depth$trainingData)$panel.args[[1]] %>% 
  as.data.frame()  
pdp_burn_depth_tc_wet <- partial(fin_mdl_burn_depth, pred.var="tc_wet", plot=TRUE, chull=TRUE, train=fin_mdl_burn_depth$trainingData)$panel.args[[1]] %>% 
  as.data.frame()  
pdp_burn_depth_slope <- partial(fin_mdl_burn_depth, pred.var="slope", plot=TRUE, chull=TRUE, train=fin_mdl_burn_depth$trainingData)$panel.args[[1]] %>% 
  as.data.frame()

pdp_burn_depth_rugosity_plot <- ggplot(pdp_burn_depth_rugosity, aes(y=y, x=x)) +
  # Scales
  scale_y_continuous(limits=c(9,14), breaks=seq(9,14,1)) +
  scale_x_continuous(limits=c(1,5), breaks=seq(0,5,1)) +
  # Themes
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "in"),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.background=element_rect(colour="white", size=1),
        panel.border=element_rect(colour="black", fill=NA, size=1),
        text=element_text(size=12),
        axis.text.y=element_text(angle=0, hjust=0.5, colour="black"),
        axis.title.y.left=element_text(margin=margin(t=0, r=10, b=0, l=0)),
        axis.title.y.right=element_text(margin=margin(t=0, r=0, b=0, l=10)),
        axis.text.x=element_text(angle=0, hjust=0.5, colour="black"),
        axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0)),
        panel.background=element_rect(fill='white'),
        legend.position="none")+
        #strip.text = element_blank()) +
  # Data
  geom_line(data=pdp_burn_depth_rugosity, aes(y=y, x=x), size=1) +
  # Axis labels
  labs(y="Mean predicted burn depth", x="Rugosity")

pdp_burn_depth_dnbr_plot <- ggplot(pdp_burn_depth_dnbr, aes(y=y, x=x)) +
  # Scales
  scale_y_continuous(limits=c(9,12), breaks=seq(9,12,1)) +
  #scale_x_continuous(limits=c(1,5), breaks=seq(0,5,1)) +
  # Themes
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "in"),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.background=element_rect(colour="white", size=1),
        panel.border=element_rect(colour="black", fill=NA, size=1),
        text=element_text(size=12),
        axis.text.y=element_text(angle=0, hjust=0.5, colour="black"),
        axis.title.y.left=element_text(margin=margin(t=0, r=10, b=0, l=0)),
        axis.title.y.right=element_text(margin=margin(t=0, r=0, b=0, l=10)),
        axis.text.x=element_text(angle=0, hjust=0.5, colour="black"),
        axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0)),
        panel.background=element_rect(fill='white'),
        legend.position="none")+
  #strip.text = element_blank()) +
  # Data
  geom_line(data=pdp_burn_depth_dnbr, aes(y=y, x=x), size=1) +
  # Axis labels
  labs(y="Mean predicted burn depth", x="dNBR")

pdp_burn_depth_tcw_plot <- ggplot(pdp_burn_depth_tc_wet, aes(y=y, x=x)) +
  # Scales
  scale_y_continuous(limits=c(9.9,11), breaks=seq(10,11,0.2)) +
  scale_x_continuous(limits=c(-0.15,-0.05), breaks=seq(-0.14, -0.04, 0.02)) +
  # Themes
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "in"),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.background=element_rect(colour="white", size=1),
        panel.border=element_rect(colour="black", fill=NA, size=1),
        text=element_text(size=12),
        axis.text.y=element_text(angle=0, hjust=0.5, colour="black"),
        axis.title.y.left=element_text(margin=margin(t=0, r=10, b=0, l=0)),
        axis.title.y.right=element_text(margin=margin(t=0, r=0, b=0, l=10)),
        axis.text.x=element_text(angle=0, hjust=0.5, colour="black"),
        axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0)),
        panel.background=element_rect(fill='white'),
        legend.position="none")+
  #strip.text = element_blank()) +
  # Data
  geom_line(data=pdp_burn_depth_tc_wet, aes(y=y, x=x), size=1) +
  # Axis labels
  labs(y="Mean predicted burn depth", x="Tasseled cap wetness")

pdp_burn_depth_slope_plot <- ggplot(pdp_burn_depth_slope, aes(y=y, x=x)) +
  # Scales
  scale_y_continuous(limits=c(9.5,10.6), breaks=seq(9.5,10.5,0.25)) +
  scale_x_continuous(limits=c(0,6), breaks=seq(0,6,1)) +
  # Themes
  theme(plot.margin=unit(c(0.1,0.1,0.1,0.1), "in"),
        panel.grid.minor=element_blank(),
        panel.grid.major=element_blank(),
        plot.background=element_rect(colour="white", size=1),
        panel.border=element_rect(colour="black", fill=NA, size=1),
        text=element_text(size=12),
        axis.text.y=element_text(angle=0, hjust=0.5, colour="black"),
        axis.title.y.left=element_text(margin=margin(t=0, r=10, b=0, l=0)),
        axis.title.y.right=element_text(margin=margin(t=0, r=0, b=0, l=10)),
        axis.text.x=element_text(angle=0, hjust=0.5, colour="black"),
        axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0)),
        panel.background=element_rect(fill='white'),
        legend.position="none")+
  #strip.text = element_blank()) +
  # Data
  geom_line(data=pdp_burn_depth_slope, aes(y=y, x=x), size=1) +
  # Axis labels
  labs(y="Mean predicted burn depth", x="Slope")

# Arrange the plots in a 2x2 grid
grid.arrange(pdp_burn_depth_rugosity_plot, 
             pdp_burn_depth_dnbr_plot, 
             pdp_burn_depth_tcw_plot, 
             pdp_burn_depth_slope_plot, 
             nrow=2)