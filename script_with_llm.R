#-------------------------------------------------------------------------------
# Final Project
# Script:    finalProject_v1.R
# Title:     Estimating Above-Ground Biomass (AGB)
# Author:    Radhika Dhuri
# Date:      2026-03-05
#
# Conclusion:
# - Develop allometry model for AGB estimation using in-situ measurements
# - Estimate AGB from lidar height map
#-------------------------------------------------------------------------------

library(ggplot2)
library(gridExtra)
library(randomForest)
library(caret)
library(terra)
library(corrplot)
library(dplyr) # MOVED: was loaded mid-script in Step 3; load all packages at the top

# ==============================================================================
# Step 1: Pre-processing
# ==============================================================================

# --- Load rasters ---
lvis   <- rast("lvis_agb_mean.tif")
c_vv   <- rast("S1_gamma_VV.nc")
c_vh   <- rast("S1_gamma_VH.nc")
fsar_l <- rast("FSAR_L_band_gamma.tif")
fsar_p <- rast("FSAR_P_band_gamma.tif")

# Check scale of C-band: values -25 to 5 = dB; 0 to 1 = linear
hist(c_vv[[1]])

# BUG: hist(fsar_l) plots all 3 bands without labelling which is which.
# Specify a single layer for clarity.
hist(fsar_l[[1]])

# --- Load and reproject GEDI AGB point data ---
# GEDI is provided as shapefile with supporting files (.cpg, .prj, etc.)
gedi     <- vect("gedi/gedi_agb.shp")
gedi_utm <- project(gedi, "EPSG:32732")

# --- Aggregate Sentinel-1 temporal steps to a single median image ---
# Median is preferred over mean for radar data (more robust to outliers).
# BUG (RISK): If c_vv/c_vh are in dB scale, aggregating with median is physically
# incorrect. Convert to linear first: app(c_vv, function(x) 10^(x/10)), take the
# median, then convert back if needed. Verify the scale using hist() above.
c_vv_agg <- app(c_vv, median, na.rm = TRUE)
c_vh_agg <- app(c_vh, median, na.rm = TRUE)

# --- Resample F-SAR from 2 m to 25 m spatial resolution ---
# BUG: fact=12.5 is not a valid integer — terra requires a whole number.
# fact=12.5 may error or silently round, giving the wrong resolution.
# Use fact=13 or use resample() with an explicit 25 m template raster.
template <- rast(ext(fsar_l), resolution = 25, crs = crs(fsar_l))
f_l <- resample(fsar_l, template, method="bilinear")
f_p <- resample(fsar_p, template, method="bilinear")

# --- Split F-SAR bands into HH, VH, VV polarisations ---
l_hh <- f_l[[1]]
l_vh <- f_l[[2]]
l_vv <- f_l[[3]]

p_hh <- f_p[[1]]
p_vh <- f_p[[2]]
p_vv <- f_p[[3]]

# --- Compute Cross Ratio (CR) ---
# CR = VH / VV on linear scale  |  CR = VH - VV on dB scale
cr_l <- l_vh / l_vv
cr_p <- p_vh / p_vv
cr_c <- c_vh_agg / c_vv_agg

plot(cr_l, range = c(0, 1), col = rev(topo.colors(15)), main = "Cross Ratio FSAR L-Band")
plot(cr_p, range = c(0, 1), col = rev(topo.colors(15)), main = "Cross Ratio FSAR P-Band")
plot(cr_c, range = c(0, 1), col = rev(topo.colors(15)), main = "Cross Ratio Sentinel C-Band")

# --- Compute Radar Vegetation Index (RVI) ---
# Full-pol RVI: 8*VH / (HH + VV + 2*VH)
rvi_l <- 8 * l_vh / (l_hh + l_vv + 2 * l_vh)
rvi_p <- 8 * p_vh / (p_hh + p_vv + 2 * p_vh)

plot(rvi_l, range = c(0, 2), col = rev(topo.colors(50)), main = "RVI FSAR L-Band")
hist(rvi_p)
# BUG: main label said "RVI FSAR L-Band" — this is the P-band plot.
plot(rvi_p, range = c(0, 2), col = rev(topo.colors(50)), main = "RVI FSAR P-Band")

# ==============================================================================
# Step 2: Correlation Analysis
# ==============================================================================

# --- Combine rasters per sensor ---
c_combi <- c(c_vh_agg, c_vv_agg, cr_c)
names(c_combi) <- c("SVH", "SVV", "SCR")

l_combi <- c(l_vh, l_hh, l_vv, cr_l, rvi_l)
names(l_combi) <- c("LHV", "LHH", "LVV", "LCR", "LRVI")

p_combi <- c(p_vh, p_hh, p_vv, cr_p, rvi_p)
names(p_combi) <- c("PHV", "PHH", "PVV", "PCR", "PRVI")

# --- Helper function: extract raster values at GEDI points ---
extract_at_gedi <- function(raster_composite, gedi_points) {
  extracted <- extract(raster_composite, gedi_points, bind = TRUE)
  return(as.data.frame(extracted))
}

# --- Extract values at GEDI locations ---
S1_gedi_df <- extract_at_gedi(c_combi, gedi_utm)
L_gedi_df  <- extract_at_gedi(l_combi, gedi_utm)
P_gedi_df  <- extract_at_gedi(p_combi, gedi_utm)

# --- Combine all into one dataframe ---
alldata <- cbind(S1_gedi_df, L_gedi_df, P_gedi_df)
alldata <- na.omit(alldata)
alldata <- alldata[, !duplicated(names(alldata))]

# --- Convert backscatter columns to dB scale ---
log_cols <- c("LHV", "PHV", "LHH", "PHH", "LVV", "PVV", "SVH", "SVV")
for (col in log_cols) {
  alldata[paste0(col, "_log")] <- 10 * log10(alldata[[col]])
}

# --- Histograms for visual inspection ---
par(mfrow = c(2, 3))

# BUG: Original code referenced S1_gedi_df$VH and S1_gedi_df$VV, which don't exist.
# The correct column names (from c_combi) are SVH and SVV.
hist(S1_gedi_df$SVH, main = "C-band VH",  xlab = "Backscatter",     col = "lightblue")
hist(S1_gedi_df$SVV, main = "C-band VV",  xlab = "Backscatter",     col = "lightblue")
hist(S1_gedi_df$agb, main = "GEDI AGB",   xlab = "Biomass (Mg/ha)", col = "lightgreen")

# BUG: Original histogram titles said "p-band" for L_gedi_df data. Fixed to "L-band".
hist(L_gedi_df$LHV,  main = "L-band VH",  xlab = "Backscatter",     col = "lightblue")
hist(L_gedi_df$LVV,  main = "L-band VV",  xlab = "Backscatter",     col = "lightblue")
hist(L_gedi_df$agb,  main = "GEDI AGB",   xlab = "Biomass (Mg/ha)", col = "lightgreen")

hist(P_gedi_df$PHV,  main = "P-band VH",  xlab = "Backscatter",     col = "lightblue")
hist(P_gedi_df$PVV,  main = "P-band VV",  xlab = "Backscatter",     col = "lightblue")
hist(P_gedi_df$agb,  main = "GEDI AGB",   xlab = "Biomass (Mg/ha)", col = "lightgreen")

par(mfrow = c(1, 1))

# --- Correlation matrices (linear scale) ---
S1_cor <- cor(S1_gedi_df[, 4:ncol(S1_gedi_df)], use = "complete.obs", method = "spearman")
L_cor  <- cor(L_gedi_df[,  4:ncol(L_gedi_df)],  use = "complete.obs", method = "spearman")
P_cor  <- cor(P_gedi_df[,  4:ncol(P_gedi_df)],  use = "complete.obs", method = "spearman")

# --- Correlation matrices (dB scale) ---
# BUG: Original code applied 10*log10() to ALL columns including 'agb',
# which converts AGB to dB — physically meaningless. Only backscatter columns
# should be log-transformed. Fixed below by excluding 'agb' before transforming.
S1_log_df <- S1_gedi_df[, 4:ncol(S1_gedi_df)]
S1_log_df[setdiff(names(S1_log_df), "agb")] <- 10 * log10(S1_log_df[setdiff(names(S1_log_df), "agb")])

L_log_df  <- L_gedi_df[, 4:ncol(L_gedi_df)]
L_log_df[setdiff(names(L_log_df), "agb")]  <- 10 * log10(L_log_df[setdiff(names(L_log_df), "agb")])

P_log_df  <- P_gedi_df[, 4:ncol(P_gedi_df)]
P_log_df[setdiff(names(P_log_df), "agb")]  <- 10 * log10(P_log_df[setdiff(names(P_log_df), "agb")])

cor_matrix_S1_log <- cor(S1_log_df, use = "complete.obs", method = "spearman")
cor_matrix_L_log  <- cor(L_log_df,  use = "complete.obs", method = "spearman")
cor_matrix_P_log  <- cor(P_log_df,  use = "complete.obs", method = "spearman")

# --- Corrplot ---
corrplot(
  cor_matrix_P_log,
  method      = "color",
  type        = "upper",
  order       = "hclust",
  addCoef.col = "black",
  tl.col      = "black",
  tl.srt      = 45,
  diag        = FALSE
)

# ==============================================================================
# Step 3: Linear Regression with Cross-Validation
# ==============================================================================

P_gedi_df_omit <- na.omit(P_gedi_df)
L_gedi_df_omit <- na.omit(L_gedi_df)

ctrl <- trainControl(
  method      = "cv",
  number      = 5,
  verboseIter = TRUE
)

# --- P-band models ---
set.seed(123)
model_phv_cv <- train(log(agb) ~ PHV, data = P_gedi_df_omit, method = "lm", trControl = ctrl)
print(model_phv_cv)
print(model_phv_cv$results)

set.seed(123)
model_phh_cv <- train(log(agb) ~ PHH, data = P_gedi_df_omit, method = "lm", trControl = ctrl)
print(model_phh_cv)
print(model_phh_cv$results)

set.seed(123)
model_pvv_cv <- train(log(agb) ~ PVV, data = P_gedi_df_omit, method = "lm", trControl = ctrl)
print(model_pvv_cv)
print(model_pvv_cv$results)

# --- L-band models ---
set.seed(123)
model_lhv_cv <- train(log(agb) ~ LHV, data = L_gedi_df_omit, method = "lm", trControl = ctrl)
print(model_lhv_cv)
print(model_lhv_cv$results)

set.seed(123)
model_lhh_cv <- train(log(agb) ~ LHH, data = L_gedi_df_omit, method = "lm", trControl = ctrl)
# BUG: Original code printed model_phh_cv here instead of model_lhh_cv.
# This showed P-band HH results in place of L-band HH results.
print(model_lhh_cv)
print(model_lhh_cv$results)

set.seed(123)
model_lvv_cv <- train(log(agb) ~ LVV, data = L_gedi_df_omit, method = "lm", trControl = ctrl)
print(model_lvv_cv)
print(model_lvv_cv$results)

# --- Combine all model results into one summary table ---
models <- list(
  lvv = model_lvv_cv,
  lhh = model_lhh_cv,
  lhv = model_lhv_cv,
  pvv = model_pvv_cv,
  phh = model_phh_cv,
  phv = model_phv_cv
)

linear_reg_results_overview <- bind_rows(
  lapply(names(models), function(nm) {
    df       <- models[[nm]]$results
    df$model <- nm
    df
  })
) %>% select(model, everything())

print(linear_reg_results_overview)

# ==============================================================================
# Step 3.5: AGB Estimation (Linear Regression)
# ==============================================================================



# ==============================================================================
# Step 4: Random Forest with Cross-Validation
# ==============================================================================

# --- Combine L-band and P-band features ---
all_fsar <- cbind(
  L_gedi_df[, c("LHV", "LHH", "LVV", "LCR", "LRVI","agb")],
  P_gedi_df[, c("PHV", "PHH", "PVV", "PCR", "PRVI")]
)
all_fsar <- all_fsar[, !duplicated(names(all_fsar))]

ctrl_rf <- trainControl(
  method          = "cv",
  number          = 5,
  summaryFunction = defaultSummary,
  savePredictions = "final"
)

# RF model using linear-scale backscatter as predictors.
# BUG: Original formula used 10*log10(agb) as the response variable — this
# converts AGB to dB, which has no physical justification. Use log(agb) for a
# log-linear model, or agb directly. Fixed to log(agb) to match other models.
rf_model <- train(
  log(agb) ~ LHV + LHH + LVV + LCR + LRVI + PHV + PHH + PVV + PCR + PRVI,
  data      = alldata,
  method    = "rf",
  trControl = ctrl_rf
)
rf_model

# RF model using dB-scale backscatter predictors
rf_model_log <- train(
  log(agb) ~ LHV_log + LHH_log + LVV_log + PHV_log + PHH_log + PVV_log,
  data      = alldata,
  method    = "rf",
  trControl = ctrl_rf
)
print(rf_model_log)

# ==============================================================================
# Step 4.5: AGB Estimation (Random Forest)
# ==============================================================================

alldata_combi    <- c(l_combi, p_combi)
alldata_combi_df <- as.data.frame(alldata_combi, xy = TRUE)
alldata_combi_df <- na.omit(alldata_combi_df)
names(alldata_combi_df) <- c("lon", "lat", names(alldata_combi))

# Predict AGB across the full raster extent
prediction                <- predict(rf_model, newdata = alldata_combi_df)
alldata_combi_df$agb_pred <- prediction

# Back-transform: model trained on log(agb), so back-transform with exp()
# BUG: Original code used 10^(agb_pred/10) which reverses a 10*log10() transform,
# not a log() transform. Since the model uses log(agb), the correct inverse is exp().
r_pred <- rast(alldata_combi_df, type = "xyz", crs = crs(alldata_combi))

plot(exp(r_pred$agb_pred), main = "AGB Estimation RF (25 m)", range = c(0, 800))

# Resample to LVIS resolution for comparison
r_pred_50 <- resample(r_pred, lvis, method = "bilinear")

plot(exp(r_pred_50$agb_pred), main = "AGB Estimation RF (LVIS resolution)", range = c(0, 800))
plot(lvis, main = "AGB LVIS", range = c(0, 800))

# ==============================================================================
# Step 5: Validation
# ==============================================================================

validation_df <- as.data.frame(c(r_pred_50$agb_pred, lvis))
validation_df <- na.omit(validation_df)

# Back-transform predicted values from log scale
validation_df$agb_pred_bt <- exp(validation_df$agb_pred)

validation_stats <- data.frame(
  n_pairs      = nrow(validation_df),
  rmse         = sqrt(mean((validation_df$lvis_agb_mean - validation_df$agb_pred_bt)^2)),
  r2           = cor(validation_df$lvis_agb_mean, validation_df$agb_pred_bt)^2,
  mae          = mean(abs(validation_df$lvis_agb_mean - validation_df$agb_pred_bt)),
  bias         = mean(validation_df$agb_pred_bt - validation_df$lvis_agb_mean),
  rmse_percent = 100 * sqrt(mean((validation_df$lvis_agb_mean - validation_df$agb_pred_bt)^2)) /
    mean(validation_df$lvis_agb_mean)
)
print(validation_stats)

# BUG: ggplot blocks referenced 'comparison_df' which was never created anywhere.
# Built from validation_df here with the expected column names.
comparison_df              <- validation_df
comparison_df$lvis_agb     <- comparison_df$lvis_agb_mean
comparison_df$rf_predicted <- comparison_df$agb_pred_bt

# --- Scatter plot: RF predictions vs LVIS AGB ---
p1 <- ggplot(comparison_df, aes(x = lvis_agb, y = rf_predicted)) +
  geom_point(alpha = 0.5, color = "darkgreen") +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", alpha = 0.3) +
  annotate("text",
           x     = max(comparison_df$lvis_agb) * 0.7,
           y     = min(comparison_df$rf_predicted) * 1.1,
           label = paste0("R² = ", round(validation_stats$r2, 3),
                          "\nRMSE = ", round(validation_stats$rmse, 2),
                          "\nn = ", validation_stats$n_pairs),
           hjust = 0, size = 5) +
  labs(title = "RF Predictions vs LVIS AGB",
       x = "LVIS AGB (Mg/ha)", y = "RF Predicted AGB (Mg/ha)") +
  theme_minimal() +
  coord_fixed()

# --- Residual plot ---
comparison_df$residuals <- comparison_df$rf_predicted - comparison_df$lvis_agb
p2 <- ggplot(comparison_df, aes(x = lvis_agb, y = residuals)) +
  geom_point(alpha = 0.5, color = "darkgreen") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 1) +
  geom_smooth(method = "loess", color = "blue", se = TRUE) +
  labs(title = "Residuals vs LVIS AGB",
       x = "LVIS AGB (Mg/ha)", y = "Residuals (RF - LVIS)") +
  theme_minimal()

# --- Error distribution plot ---
# NOTE: ..density.. is deprecated in newer ggplot2; replaced with after_stat(density).
p3 <- ggplot(comparison_df, aes(x = residuals)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "forestgreen", alpha = 0.7) +
  geom_density(color = "darkgreen", linewidth = 1) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
  labs(title = "Distribution of Prediction Errors",
       x = "Prediction Error (RF - LVIS)", y = "Density") +
  theme_minimal()

grid.arrange(p1, p2, p3, ncol = 2, nrow = 2)
