#-------------------------------------------------------------------------------
# Final Project
# Script:    script_main.R
# Title:     Estimating above-ground biomass (AGB)
# Author:    Radhika Dhuri
# Date:      2026-03-05
#
# Conclusion:
# - Develop allometry model for AGB estimation using in-situ measurements
# - Estimate AGB from lidar height map.
#-------------------------------------------------------------------------------


# Loading packages
library(randomForest)
library(caret)
library(sf)
library(terra)
library(corrplot)
install.packages('corrplot')


########################## Step 1: Pre-processing #################################
# Read raster data
# c band sentinal 
c_vv <- rast("S1_gamma_VV.nc")
c_vv
hist(c_vv[[1]]) # if the value between -25 to 5 then db; 0 to 1 then linear
c_vh <- rast("S1_gamma_VH.nc")
c_vh
# l and p band fsar
fsar_l <- rast("FSAR_L_band_gamma.tif")
fsar_l
hist(fsar_l)
fsar_p <- rast("FSAR_P_band_gamma.tif")
fsar_p
# read and resample gedi (given in shp format with other supporting file like cpg, prj etc. ) 
# to same projection as the rest

gedi <- vect("gedi/gedi_agb.shp")

gedi_utm <- project(gedi, crs(c_vv))
gedi_utm
# aggregate temporal steps in sentinal c band to temporal average at linear scale
c_vv_agg <- app(c_vv, median, na.rm=TRUE) #median is better than mean in radar 
c_vh_agg <- app(c_vh, median, na.rm=TRUE)
# change spatial resolution from 2m to 25m for fsar
f_l <- aggregate(fsar_l, fact=12.5, fun=mean, na.rm=TRUE)

f_p <- aggregate(fsar_p, fact=12.5, fun=mean, na.rm=TRUE)
# divide fsar bands in vv and vh
l_hh <- f_l[[1]]
l_hh
l_vh <- f_l[[2]]
l_vh
l_vv <- f_l[[3]]
l_vv

p_hh <- f_p[[1]]
p_hh
p_vh <- f_p[[2]]
p_vh
p_vv <- f_p[[3]]
p_vv
# Calculate backscatter cross ratio and radar vegetation index for fsar l and p
# compute the cross ratio (still need to check the scale)
# CR = VH - VV (on logarithmic (dB) scale)
# CR = VH / VV (on linear scale)
cr_l <- l_vh / l_vv
plot(cr_l, range=c(0, 1), col=rev(topo.colors(15)), main = "Cross Ratio FSAR L-Band")
cr_p <- p_vh / p_vv
plot(cr_p, range=c(0, 1), col=rev(topo.colors(15)), main = "Cross Ratio FSAR P-Band")
cr_c <- c_vh_agg / c_vv_agg
plot(cr_c, range=c(0, 1), col=rev(topo.colors(15)), main = "Cross Ratio Sentinel C-Band")
# radar vegetation index (still check the formula)
rvi_l <- 8*l_vh / (l_hh + l_vv + 2*l_vh)  
plot(rvi_l, range=c(0, 2), main = "RVI FSAR L-Band", col=rev(topo.colors(50)))
rvi_p <- 8*l_vh / (l_hh + l_vv + 2*l_vh) 
hist(rvi_p)
plot(rvi_p, range=c(0, 2), main = "RVI FSAR L-Band", col=rev(topo.colors(50)))


########################## Step 2: Correlation analysis  #################################
#Combine the raster
#Sentinel-1
c_combi <- c(c_vh_agg, c_vv_agg, cr_c)
names(c_combi) <- c("VH", "VV", "CR")

#FSAR-L
l_combi <- c(l_vh, l_hh, l_vv, cr_l, rvi_l)
names(l_combi) <- c("LHV", "LHH", "LVV", "LCR", "LRVI")

#FSAR-L
p_combi <- c(p_vh, p_hh, p_vv, cr_p, rvi_p)
names(p_combi) <- c("PHV", "PHH", "PVV", "PCR", "PRVI")

extract_at_gedi <- function(raster_composite, gedi_points) {
  # Extract values
  extracted <- extract(raster_composite, gedi_points, bind = TRUE)
  # Convert to data frame
  df <- as.data.frame(extracted)
  df <- na.omit(df)
  return(df)
}

# Extract for each band
S1_gedi_df <- extract_at_gedi(c_combi, gedi_utm)
L_gedi_df <- extract_at_gedi(l_combi, gedi_utm)
P_gedi_df <- extract_at_gedi(p_combi, gedi_utm)

# Check column names - identify which column is AGB
names(S1_gedi_df)
S1_gedi_df

# Set up plotting area
par(mfrow = c(2, 3))
names(P_gedi_df)
# Histograms for C-band
hist(S1_gedi_df$VH, main = "C-band VH", xlab = "Backscatter", col = "lightblue")
hist(S1_gedi_df$VV, main = "C-band VV", xlab = "Backscatter", col = "lightblue")
hist(S1_gedi_df$agb, main = "GEDI AGB", xlab = "Biomass (Mg/ha)", col = "lightgreen")
hist(L_gedi_df$LHV, main = "p-band VH", xlab = "Backscatter", col = "lightblue")
hist(L_gedi_df$LVV, main = "p-band VV", xlab = "Backscatter", col = "lightblue")
hist(L_gedi_df$agb, main = "GEDI AGB", xlab = "Biomass (Mg/ha)", col = "lightgreen")
hist(P_gedi_df$PHV, main = "p-band VH", xlab = "Backscatter", col = "lightblue")
hist(P_gedi_df$PVV, main = "p-band VV", xlab = "Backscatter", col = "lightblue")
hist(P_gedi_df$agb, main = "GEDI AGB", xlab = "Biomass (Mg/ha)", col = "lightgreen")

# Correlation matrix (linear)
S1_cor <- cor(S1_gedi_df[, 4:ncol(S1_gedi_df)], use = "complete.obs", method = "spearman")
L_cor <- cor(L_gedi_df[, 4:ncol(L_gedi_df)], use = "complete.obs", method = "spearman")
P_cor <- cor(P_gedi_df[, 4:ncol(P_gedi_df)], use = "complete.obs", method = "spearman")
L_cor
# Correlation matrix (decible)
S1_gedi_df_log <- 10*log10(S1_gedi_df[, 4:ncol(S1_gedi_df)])
L_gedi_df_log <- 10*log10(L_gedi_df[, 4:ncol(L_gedi_df)])
P_gedi_df_log <- 10*log10(P_gedi_df[, 4:ncol(P_gedi_df)])
cor_matrix_S1_log <- cor(S1_gedi_df_log, use = "complete.obs", method = "spearman")
cor_matrix_L_log <- cor(L_gedi_df_log, use = "complete.obs", method = "spearman")
cor_matrix_P_log <- cor(P_gedi_df_log, use = "complete.obs", method = "spearman")
cor_matrix_L_log

corrplot(
  cor_matrix_L_log,
  method = "color",
  type = "upper",
  order = "hclust",
  addCoef.col = "black",
  tl.col = "black",
  tl.srt = 45,
  diag = FALSE
)
########################## Step 3: Linear regression with cross validation #################################
#linear model
lm_P_log <- lm(log(agb) ~ PHV, data = P_gedi_df)
summary(lm_P_log)
#cross validation
# Set up cross-validation
ctrl <- trainControl(method = "cv",     # cross-validation
                     number = 5,         # number of folds
                     verboseIter = TRUE) # show progress

# Train model of HV with CV
set.seed(123)
model_phv_cv <- train(log(agb) ~ PHV, 
                      data = P_gedi_df,
                      method = "lm",        # linear regression
                      trControl = ctrl)

# View results
print(model_phv_cv)
print(model_phv_cv$results)  # CV performance metrics

# Train model of HH with CV
set.seed(123)
model_phh_cv <- train(log(agb) ~ PHH, 
                      data = P_gedi_df,
                      method = "lm",        # linear regression
                      trControl = ctrl)

# View results
print(model_phh_cv)
print(model_phh_cv$results)  # CV performance metrics

# Train model of VV with CV
set.seed(123)
model_pvv_cv <- train(log(agb) ~ PVV, 
                      data = P_gedi_df,
                      method = "lm",        # linear regression
                      trControl = ctrl)

# View results
print(model_pvv_cv)
print(model_pvv_cv$results)  # CV performance metrics

########################## Step 3.5: AGB estimation #################################
#combine l and p band values
all_fsar <- cbind(
  L_gedi_df[, c("LHV", "LHH", "LVV", "LCR", "LRVI","agb")],
  P_gedi_df[, c("PHV", "PHH", "PVV", "PCR", "PRVI")]
)
names(all_fsar)

#random forest
ctrl <- trainControl(
  method = "cv",           # Use cross-validation
  number = 5,             # Number of folds (k=10)
  summaryFunction = defaultSummary,  # Use regression metrics (RMSE, R²)
  savePredictions = "final"
)

rf_model <- train(
  agb ~ LHV + LHH + LVV + LCR + LRVI + PHV + PHH + PVV + PCR + PRVI,  # All P-band variables
  data = all_fsar,
  method = "rf",
  trControl = ctrl,
)
rf_model

########################## Step 4: Random Forest with cross validation #################################

########################## Step 4.5: AGB estimation #################################


########################## Step 5: Validation #################################


