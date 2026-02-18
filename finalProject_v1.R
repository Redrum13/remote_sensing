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
library(sf)
library(terra)

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

########################## Step 3: Linear regression with cross validation #################################

########################## Step 3.5: AGB estimation #################################


########################## Step 4: Random Forest with cross validation #################################

########################## Step 4.5: AGB estimation #################################


########################## Step 5: Validation #################################


