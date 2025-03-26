#########################################
##' Used-habitat calibration (UHC) plots
##' to validate turaco habitat selection
##' models
##' Nicholas J. Russo and Antoine Tekam
##' nicholasrusso@g.harvard.edu
#########################################

## Clear environment
rm(list=ls())
graphics.off()

##' Load in packages
## Plotting
library(ggplot2)
library(ggpubr)
library(viridis)

## Movement analysis and data management
library(lubridate)
library(dplyr)
library(amt)
library(move)
library(adehabitatHR)

## Spatial products
library(terra)
library(sf)

## Load in turaco data, use Movebank format
turacos <- read.csv("~/Turaco Data/Turaco_GPS_16Dec2024.csv")
turacos <- turacos[!is.na(turacos$location.long),] # get rid of NA coordinates
turacos$tag.local.identifier <- as.character(turacos$tag.local.identifier)
turacos$timestamps <- as.POSIXct(strptime(turacos$study.local.timestamp, "%Y-%m-%d %H:%M:%S", tz="GMT"))

# Turacos
df93 <- turacos[which(turacos$tag.local.identifier=="8893"),]
df8895 <- turacos[which(turacos$tag.local.identifier=="8895"),]
df51 <- turacos[which(turacos$tag.local.identifier=="11851"),]
df53 <- turacos[which(turacos$tag.local.identifier=="11853"),]

turacos <- rbind(df93, df8895, df51, df53)

#' Load in environmental layers
ch <- terra::rast("~/Turaco Data/ch.tif")
vc <- terra::rast("~/vc.tif")
d50 <- terra::rast("~/Turaco Data/d50.tif")
d500 <- terra::rast("~/Turaco Data/d500.tif")

## Plant Volume Density (PVD)
pvd15_20 <- terra::rast("~/Turaco Data/pvd_15_20.tif")
pvd15_20[pvd15_20 > 40] <- NA # remove outliers

#' Land Cover
swamp <- terra::rast("~/Turaco Data/swamp.tif")
names(swamp) <- "Swamp"

# Prep issf function
prep_issf <- function(bat){
  # Format as steps
  steps <- bat |>
    make_track(utm.easting, utm.northing, timestamps, crs = 32612) |>
    track_resample(minutes(120), tolerance = minutes(10)) |>
    steps()
  # Split into train (80%) and test (20%)
  set.seed(1)
  steps$train <- rbinom(n = nrow(steps),
                        size = 1, prob = 0.8)
  train <- steps[steps$train == 1, ]
  test <- steps[steps$train == 0, ]
  # Generate available steps
  train_dat <- train |>
    random_steps(n_control = 100) |>
    # Attach covariates
    amt::extract_covariates(ch) |>
    amt::extract_covariates(vc) |>
    amt::extract_covariates(d50) |>
    amt::extract_covariates(d500) |>
    amt::extract_covariates(pvd15_20) |>
    amt::extract_covariates(swamp) |>
    mutate(log_sl_ = log(sl_ + 1),
           cos_ta_ = cos(ta_)) |>
    # Drop 'train' column
    dplyr::select(-train) |> 
    # Get rid of any NAs (sometimes available steps fall outside of raster)
    na.omit()
  
  test_dat <- test |>
    random_steps(n_control = 100) |>
    # Attach covariates
    amt::extract_covariates(ch) |>
    amt::extract_covariates(vc) |>
    amt::extract_covariates(d50) |>
    amt::extract_covariates(d500) |>
    amt::extract_covariates(pvd15_20) |>
    amt::extract_covariates(swamp) |>
    mutate(log_sl_ = log(sl_ + 1),
           cos_ta_ = cos(ta_)) |>
    # Drop 'train' column
    dplyr::select(-train) |>
    # Get rid of any NAs (sometimes available steps fall outside of raster)
    na.omit()
  
  # Scale and center variables (training data)
  train_dat$CanopyHeight10m = scale(train_dat$CanopyHeight10m)
  train_dat$VCI = scale(train_dat$VCI)
  train_dat$dist2gap50 = scale(train_dat$dist2gap50)
  train_dat$dist2gap500 = scale(train_dat$dist2gap500)
  train_dat$pvd15_20 = scale(train_dat$pvd15_20)
  # Scale and center variables (test data)
  test_dat$CanopyHeight10m = scale(test_dat$CanopyHeight10m)
  test_dat$VCI = scale(test_dat$VCI)
  test_dat$dist2gap50 = scale(test_dat$dist2gap50)
  test_dat$dist2gap500 = scale(test_dat$dist2gap500)
  test_dat$pvd15_20 = scale(test_dat$pvd15_20)
  # Swamp variable as factor
  train_dat$Swamp <- as.factor(train_dat$Swamp)
  test_dat$Swamp <- as.factor(test_dat$Swamp)
  # Prep model
  issf <- fit_issf(train_dat, case_ ~ CanopyHeight10m + I(CanopyHeight10m^2) + VCI +
                     dist2gap50 + dist2gap500 + pvd15_20 + Swamp +
                     log_sl_ + cos_ta_ + strata(step_id_), model = TRUE)
  # Create uhc object
  uhc_issf <- prep_uhc(object = issf, test_dat = test_dat,
                       n_samp = 1000, verbose = TRUE)
  # Return the prepped uhc object
  return(uhc_issf)
}

##################################
#' Function to prep issf without 
#' swamp variable
##################################
prep_issf_noSwamp <- function(bat){
  # Format as steps
  steps <- bat |>
    make_track(utm.easting, utm.northing, timestamps, crs = 32612) |>
    track_resample(minutes(120), tolerance = minutes(10)) |>
    steps()
  # Split into train (80%) and test (20%)
  set.seed(1)
  steps$train <- rbinom(n = nrow(steps),
                        size = 1, prob = 0.8)
  train <- steps[steps$train == 1, ]
  test <- steps[steps$train == 0, ]
  # Generate available steps
  train_dat <- train |>
    random_steps(n_control = 100) |>
    # Attach covariates
    amt::extract_covariates(ch) |>
    amt::extract_covariates(vc) |>
    amt::extract_covariates(d50) |>
    amt::extract_covariates(d500) |>
    amt::extract_covariates(pvd15_20) |>
    mutate(log_sl_ = log(sl_ + 1),
           cos_ta_ = cos(ta_)) |>
    # Drop 'train' column
    dplyr::select(-train) |> 
    # Get rid of any NAs (sometimes available steps fall outside of raster)
    na.omit()
  
  test_dat <- test |>
    random_steps(n_control = 100) |>
    # Attach covariates
    amt::extract_covariates(ch) |>
    amt::extract_covariates(vc) |>
    amt::extract_covariates(d50) |>
    amt::extract_covariates(d500) |>
    amt::extract_covariates(pvd15_20) |>
    mutate(log_sl_ = log(sl_ + 1),
           cos_ta_ = cos(ta_)) |>
    # Drop 'train' column
    dplyr::select(-train) |>
    # Get rid of any NAs (sometimes available steps fall outside of raster)
    na.omit()
  
  # Scale and center variables (training data)
  train_dat$CanopyHeight10m = scale(train_dat$CanopyHeight10m)
  train_dat$VCI = scale(train_dat$VCI)
  train_dat$dist2gap50 = scale(train_dat$dist2gap50)
  train_dat$dist2gap500 = scale(train_dat$dist2gap500)
  train_dat$pvd15_20 = scale(train_dat$pvd15_20)
  # Scale and center variables (test data)
  test_dat$CanopyHeight10m = scale(test_dat$CanopyHeight10m)
  test_dat$VCI = scale(test_dat$VCI)
  test_dat$dist2gap50 = scale(test_dat$dist2gap50)
  test_dat$dist2gap500 = scale(test_dat$dist2gap500)
  test_dat$pvd15_20 = scale(test_dat$pvd15_20)
  # Prep model
  issf <- fit_issf(train_dat, case_ ~ CanopyHeight10m + I(CanopyHeight10m^2) + VCI +
                     dist2gap50 + dist2gap500 + pvd15_20 +
                     log_sl_ + cos_ta_ + strata(step_id_), model = TRUE)
  # Create uhc object
  uhc_issf <- prep_uhc(object = issf, test_dat = test_dat,
                       n_samp = 1000, verbose = TRUE)
  # Return the prepped uhc object
  return(uhc_issf)
}



# Prep issfs
issf_8893 <- prep_issf(df93)
issf_8895 <- prep_issf_noSwamp(df8895)
issf_11851 <- prep_issf(df51)
issf_11853 <- prep_issf_noSwamp(df53)

# Coerce to data.frame
uhc_8893 <- as.data.frame(issf_8893)
uhc_8895 <- as.data.frame(issf_8895)
uhc_11851 <- as.data.frame(issf_11851)
uhc_11853 <- as.data.frame(issf_11853)

# Simplify sampled lines to confidence envelopes
conf_8893 <- conf_envelope(uhc_8893)
conf_8895 <- conf_envelope(uhc_8895)
conf_11851 <- conf_envelope(uhc_11851)
conf_11853 <- conf_envelope(uhc_11853)

###########################################
# Make Used-habitat calibration (UHC) plots
###########################################

# 8893
windows()
par(mfrow=c(3, 2))
plot(conf_8893)

# 8895
windows()
par(mfrow=c(3, 2))
plot(conf_8895)

# 11851
windows()
par(mfrow=c(3, 2))
plot(conf_11851)

# 11853
windows()
par(mfrow=c(3, 2))
plot(conf_11853)
