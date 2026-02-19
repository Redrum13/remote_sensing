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
library(ggplot2)
library(gridExtra)
library(randomForest)
library(caret)
library(sf)
library(terra)
library(corrplot)


########################## Step 1: Pre-processing #################################
# Read raster data
lvis <- rast("lvis_agb_mean.tif")
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
names(c_combi) <- c("SVH", "SVV", "SCR")

#FSAR-L
l_combi <- c(l_vh, l_hh, l_vv, cr_l, rvi_l)
names(l_combi) <- c("LHV", "LHH", "LVV", "LCR", "LRVI")
head (l_combi)
#FSAR-L
p_combi <- c(p_vh, p_hh, p_vv, cr_p, rvi_p)
names(p_combi) <- c("PHV", "PHH", "PVV", "PCR", "PRVI")

extract_at_gedi <- function(raster_composite, gedi_points) {
  # Extract values
  extracted <- extract(raster_composite, gedi_points, bind = TRUE)
  # Convert to data frame
  df <- as.data.frame(extracted)
  return(df)
}

head(as.data.frame(gedi_utm))

# Extract for each band
S1_gedi_df <- extract_at_gedi(c_combi, gedi_utm)
L_gedi_df <- extract_at_gedi(l_combi, gedi_utm)
P_gedi_df <- extract_at_gedi(p_combi, gedi_utm)

# Check column names - identify which column is AGB
names(L_gedi_df)
S1_gedi_df

# Combine dataframes 
alldata <- cbind(S1_gedi_df, L_gedi_df)
alldata <- cbind(alldata, P_gedi_df)
alldata <- na.omit(alldata)
alldata <- alldata[, !duplicated(names(alldata))]
str(alldata)
log_cols <- c("LHV", "PHV", "LHH", "PHH", "LVV", "PVV", "SVH", "SVV")
for(col in log_cols) {
  alldata[paste0(col, "_log")] <- 10*log10(alldata[[col]])
}

head(alldata)

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
# P-band cross validation

names(P_gedi_df)
P_gedi_df <-na.omit(P_gedi_df)
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
                      data = fsar_P_gedi_df,
                      method = "lm",        # linear regression
                      trControl = ctrl)
# View results
print(model_phh_cv)
print(model_phh_cv$results)  # CV performance metrics

# Train model of VV with CV
set.seed(123)
model_pvv_cv <- train(log(agb) ~ PVV, 
                      data = fsar_P_gedi_df,
                      method = "lm",        # linear regression
                      trControl = ctrl)
# View results
print(model_pvv_cv)
print(model_pvv_cv$results)  # CV performance metrics



# L-band cross validation

# Train model of HV with CV
set.seed(123)
model_lhv_cv <- train(log(agb) ~ LHV, 
                      data = fsar_L_gedi_df,
                      method = "lm",        # linear regression
                      trControl = ctrl)

# View results
print(model_lhv_cv)
print(model_lhv_cv$results)  # CV performance metrics

# Train model of HH with CV
set.seed(123)
model_lhh_cv <- train(log(agb) ~ LHH, 
                      data = fsar_L_gedi_df,
                      method = "lm",        # linear regression
                      trControl = ctrl)

# View model_lhh_cv
print(model_phh_cv)
print(model_lhh_cv$results)  # CV performance metrics

# Train model of VV with CV
set.seed(123)
model_lvv_cv <- train(log(agb) ~ LVV, 
                      data = fsar_L_gedi_df,
                      method = "lm",        # linear regression
                      trControl = ctrl)

# View results
print(model_lvv_cv)
print(model_lvv_cv$results)  # CV performance metrics

#combine to one table
models <- list(
  lvv = model_lvv_cv,
  lhh = model_lhh_cv,
  lhv = model_lhv_cv,
  pvv = model_pvv_cv,
  phh = model_phh_cv,
  phv = model_phv_cv
)

# Combine their results into one table with a model column
library(dplyr)

linear_reg_results_overview <- bind_rows(
  lapply(names(models), function(nm) {
    df <- models[[nm]]$results
    df$model <- nm  # add a column for model name
    df
  })
)

########################## Step 3.5: AGB estimation #################################


########################## Step 4: Random Forest with cross validation #################################
#combine l and p band values
all_fsar <- cbind(
  L_gedi_df[, c("LHV", "LHH", "LVV", "LCR", "LRVI","agb")],
  P_gedi_df[, c("PHV", "PHH", "PVV", "PCR", "PRVI")]
)
names(all_fsar)
all_fsar <- all_fsar[, !duplicated(names(all_fsar))]

#random forest
ctrl <- trainControl(
  method = "cv",           # Use cross-validation
  number = 5,             # Number of folds (k=5)
  summaryFunction = defaultSummary,  # Use regression metrics (RMSE, R²)
  savePredictions = "final"
)

rf_model <- train(
  10*log10(agb) ~ LHV + LHH + LVV + LCR + LRVI + PHV + PHH + PVV + PCR + PRVI,  
  data = alldata,
  method = "rf",
  trControl = ctrl
)
rf_model
names(alldata)
rf_model_log <- train(
  10*log10(agb) ~ LHV_log + LHH_log + LVV_log + PHV_log + PHH_log + PVV_log,  
  data = alldata,
  method = "rf",
  trControl = ctrl
)
rf_model_log

########################## Step 4.5: AGB estimation #################################
alldata_combi <- c(l_combi, p_combi)
alldata_combi_df <- as.data.frame(alldata_combi, xy=TRUE)
alldata_combi_df <- na.omit(alldata_combi_df)

names(alldata_combi_df) <- c("lon", "lat", names(alldata_combi))

head(as.data.frame(alldata_combi))
# estimation using rf

prediction <- predict(rf_model, newdata=alldata_combi_df)

alldata_combi_df$agb_pred <- prediction

head(alldata_combi_df)

plot(alldata_combi_df)

r_pred <- rast(alldata_combi_df, type="xyz", crs = crs(alldata_combi))

plot(10^(r_pred$agb_pred/10), main = "AGB estimation RF linear")

r_pred_50 <- resample(r_pred, lvis, method='bilinear')

plot(10^(r_pred_50$agb_pred/10), main = "AGB estimation RF linear", range=c(0,800))

plot(lvis, main = "AGB LVIS", range=c(0,800))

lvis
########################## Step 5: Validation #################################


validation_df <- as.data.frame(c(r_pred_50$agb_pred, lvis))
validation_df <- na.omit(validation_df)
compareGeom(lvis, r_pred_50)
head(validation_df)


# Calculate validation metrics
validation_stats <- data.frame(
  n_pairs = nrow(validation_df),
  rmse = sqrt(mean((validation_df$lvis_agb_mean - validation_df$agb_pred)^2)),
  r2 = cor(validation_df$lvis_agb_mean, validation_df$agb_pred)^2,
  mae = mean(abs(validation_df$lvis_agb_mean - validation_df$agb_pred)),
  bias = mean(validation_df$agb_pred - validation_df$lvis_agb_mean),
  rmse_percent = 100 * sqrt(mean((validation_df$lvis_agb_mean - validation_df$agb_pred)^2)) / mean(validation_df$lvis_agb_mean)
)

print(validation_stats)

# Create beautiful comparison plots
p1 <- ggplot(comparison_df, aes(x = lvis_agb, y = rf_predicted)) +
  geom_point(alpha = 0.5, color = "darkgreen") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", alpha = 0.3) +
  annotate("text", x = max(comparison_df$lvis_agb)*0.7, 
           y = min(comparison_df$rf_predicted)*1.1,
           label = paste0("R² = ", round(validation_stats$r2, 3),
                          "\nRMSE = ", round(validation_stats$rmse, 2),
                          "\nn = ", validation_stats$n_pairs),
           hjust = 0, size = 5) +
  labs(title = "RF Predictions vs LVIS AGB",
       x = "LVIS AGB (Mg/ha)", y = "RF Predicted AGB (Mg/ha)") +
  theme_minimal() +
  coord_fixed()

# Residual plot
comparison_df$residuals <- comparison_df$rf_predicted - comparison_df$lvis_agb

p2 <- ggplot(comparison_df, aes(x = lvis_agb, y = residuals)) +
  geom_point(alpha = 0.5, color = "darkgreen") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 1) +
  geom_smooth(method = "loess", color = "blue", se = TRUE) +
  labs(title = "Residuals vs LVIS AGB",
       x = "LVIS AGB (Mg/ha)", y = "Residuals (RF - LVIS)") +
  theme_minimal()

# Density plot of errors
p3 <- ggplot(comparison_df, aes(x = residuals)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "forestgreen", alpha = 0.7) +
  geom_density(color = "darkgreen", size = 1) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Distribution of Prediction Errors",
       x = "Prediction Error (RF - LVIS)", y = "Density") +
  theme_minimal()

# Arrange plots
grid.arrange(p1, p2, p3, ncol = 2, nrow = 2)

