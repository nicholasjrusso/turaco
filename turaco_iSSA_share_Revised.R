#########################################
##' Turaco habitat selection
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

## Spatial statistics
library(spatstat)



## Load in turaco data, use Movebank format
turacos <- read.csv("~/Turaco Data/Turaco_GPS_Sep2024.csv")
turacos <- turacos[!is.na(turacos$location.long),] # get rid of NA coordinates
turacos$tag.local.identifier <- as.character(turacos$tag.local.identifier)
turacos$timestamps <- as.POSIXct(strptime(turacos$study.local.timestamp, "%Y-%m-%d %H:%M:%S", tz="GMT"))

# Turacos
df93 <- turacos[which(turacos$tag.local.identifier=="8893"),]
df8895 <- turacos[which(turacos$tag.local.identifier=="8895"),]
df51 <- turacos[which(turacos$tag.local.identifier=="11851"),]
df53 <- turacos[which(turacos$tag.local.identifier=="11853"),]

turacos <- rbind(df93, df8895, df51, df53)

# Create tracks from data frames
trackIt <- function(bat, step_duration) {
  # Make the track
  bat.trk <- make_track(bat, .x=utm.easting, .y=utm.northing, .t=timestamps, id=tag.local.identifier,
                        crs="EPSG:32633")
  
  #' Now we resample to 30 minute intervals, with a tolerance of 10 minutes
  bat.trk <- track_resample(bat.trk, minutes(step_duration), tolerance = minutes(10))
  
  #' If there are bursts, we may want to filter bursts with very few locations. For example, to calculate a turning angle, 
  #' we need at least three locations. So we often will want to filter out bursts with at least 3 observations:
  bat.trk <- filter_min_n_burst(bat.trk, 3)
  
  # Return the track
  return(bat.trk)
}

# Resample each track to 2-hour intervals
t93.trk <- trackIt(df93, 120)
t95.trk <- trackIt(df8895, 120) # Not a lot
t51.trk <- trackIt(df51, 120)
t53.trk <- trackIt(df53, 120) # Not a lot

#' Load in environmental layers
ch <- terra::rast("~/Turaco Data/ch.tif")
vc <- terra::rast("~/Turaco Data/vc.tif")
d50 <- terra::rast("~/Turaco Data/d50.tif")
d500 <- terra::rast("~/Turaco Data/d500.tif")

## Plant Volume Density (PVD)
pvd15_20 <- terra::rast("~/Turaco Data/pvd_15_20.tif")
pvd15_20[pvd15_20 > 40] <- NA # remove outliers

#' Land Cover
swamp <- terra::rast("~/Turaco Data/swamp.tif")
names(swamp) <- "Swamp"

#' Create one raster stack for 1) structural variables, and 
#' 2) structure variables + swamp
veg.stack <- c(ch, vc, d50, d500, pvd15_20)
veg.stack.swamp <- c(ch, vc, d50, d500, pvd15_20, swamp)

# Sort steps into bursts
ssf_93 <- steps_by_burst(t93.trk)
ssf_93 <- filter(ssf_93, !is.na(ta_))

ssf_95 <- steps_by_burst(t95.trk)
ssf_95 <- filter(ssf_95, !is.na(ta_))

ssf_51 <- steps_by_burst(t51.trk)
ssf_51 <- filter(ssf_51, !is.na(ta_))

ssf_53 <- steps_by_burst(t53.trk)
ssf_53 <- filter(ssf_53, !is.na(ta_))

#' Create random steps. We typically get a warning that "Step-lengths or turning angles contained NA, which were removed", 
#' because of the missing turning angles at the start of a burst.
set.seed(2)

# Generate random steps
ssf_93 <- random_steps(ssf_93, n_control = 100)
ssf_95 <- random_steps(ssf_95, n_control = 100)
ssf_51 <- random_steps(ssf_51, n_control = 100)
ssf_53 <- random_steps(ssf_53, n_control = 100)

# Function for extracting covariates
extractIt <- function(dat, covar1, covar2, covar3, covar4, covar5, covar6) {
  # Extract covariates
  dat <- amt::extract_covariates(dat, covar1, where = "both")
  dat <- amt::extract_covariates(dat, covar2, where = "both")
  dat <- amt::extract_covariates(dat, covar3, where = "both")
  dat <- amt::extract_covariates(dat, covar4, where = "both")
  dat <- amt::extract_covariates(dat, covar5, where = "both")
  dat <- amt::extract_covariates(dat, covar6, where = "both")
  # Add variable for hour
  dat <- mutate(dat, "hour" = hour(t1_) + minute(t1_) / 60)
  #' Remove NA's
  dat <- dat[complete.cases(dat),]
  # Return data frame
  return(dat)
}

# Turaco SSFs
ssf_93 <- extractIt(ssf_93, ch, vc, pvd15_20, d50, d500, swamp)
ssf_95 <- extractIt(ssf_95, ch, vc, pvd15_20, d50, d500, swamp)
ssf_51 <- extractIt(ssf_51, ch, vc, pvd15_20, d50, d500, swamp)
ssf_53 <- extractIt(ssf_53, ch, vc, pvd15_20, d50, d500, swamp)

# Add ID column to each data frame
ssf_93$ID <- rep("8893", length(ssf_93$x1_))
ssf_95$ID <- rep("8895", length(ssf_95$x1_))
ssf_51$ID <- rep("11851", length(ssf_51$x1_))
ssf_53$ID <- rep("11853", length(ssf_53$x1_))

# Scale and center variables
scaleIt <- function(d) {
  # Scale and center variables
  d %>% mutate(CanopyHeight10m_end = scale(CanopyHeight10m_end),
               VCI_end = scale(VCI_end),
               dist2gap50_end = scale(dist2gap50_end),
               dist2gap500_end = scale(dist2gap500_end),
               pvd15_20_end = scale(pvd15_20_end),
               case = as.numeric(case_))
}

# Center and scale
m93 <- scaleIt(ssf_93)
m95 <- scaleIt(ssf_95)
m51 <- scaleIt(ssf_51)
m53 <- scaleIt(ssf_53)

#######
# 8893
#######
fit8893 <- fit_issf(case_ ~ CanopyHeight10m_end + I(CanopyHeight10m_end^2) + VCI_end + 
                      dist2gap50_end + dist2gap500_end + pvd15_20_end + Swamp_end +
                      log(sl_ + 1) + cos(ta_) + strata(step_id_),
                    data = m93, model = TRUE) # global model 3


summary(fit8893)

#######
# 8895
#######
# No convergence when swamp is included
fit8895 <- fit_issf(case_ ~ CanopyHeight10m_end + I(CanopyHeight10m_end^2) + VCI_end + 
                      dist2gap50_end + dist2gap500_end + pvd15_20_end + 
                      log(sl_ + 1) + cos(ta_) + strata(step_id_),
                    data = m95, model = TRUE) # global model 3

summary(fit8895)

#######
# 11851
#######
fit11851 <- fit_issf(case_ ~ CanopyHeight10m_end + I(CanopyHeight10m_end^2) + VCI_end + 
                       dist2gap50_end + dist2gap500_end + pvd15_20_end + Swamp_end +
                       log(sl_ + 1) + cos(ta_) + strata(step_id_),
                     data = m51, model = TRUE) # global model 3

summary(fit11851)

#######
# 11853
#######
# No convergence when swamp is included
fit11853 <- fit_issf(case_ ~ CanopyHeight10m_end + I(CanopyHeight10m_end^2) + VCI_end + 
                       dist2gap50_end + dist2gap500_end + pvd15_20_end +
                       log(sl_ + 1) + cos(ta_) + strata(step_id_),
                     data = m53, model = TRUE) # global model 3

summary(fit11853)

#################################################
## Figures
#################################################
#########################################
## Confidence intervals
#########################################
getConfints <- function(mod, ID){
  # Pull CI from model and create data frame
  conf <- as.data.frame(mod$conf.int)
  # Add bird ID, species, sex, and season
  conf$ID <- rep(ID, length(conf$`exp(coef)`))
  # Return data frame
  return(conf)
}

# Confidence intervals
confints93 <- getConfints(summary(fit8893), "8893")
confints95 <- getConfints(summary(fit8895), "8895")
confints51 <- getConfints(summary(fit11851), "11851")
confints53 <- getConfints(summary(fit11853), "11853")


get_Coefs <- function(num, name, numBirds){
  # Rbind the confints
  Coefs<- rbind(confints93[num,], confints95[num,], 
                confints51[num,], confints53[num,]
  )
  # Add the covariate name
  Coefs$covar <- rep(name, numBirds)
  # log-transform covariate estimates
  Coefs[c("coef", "-coef", "lower.95", "upper.95")] <- 
    lapply(Coefs[c("exp(coef)", "exp(-coef)", "lower .95", "upper .95")], function(x) log(x))
  # Return data frame
  return(Coefs)
}

get_CoefsSwamp <- function(num, name, numBirds){
  # Rbind the confints
  Coefs<- rbind(confints93[num,],
                confints51[num,]
  )
  # Add the covariate name
  Coefs$covar <- rep(name, numBirds)
  # log-transform covariate estimates
  Coefs[c("coef", "-coef", "lower.95", "upper.95")] <- 
    lapply(Coefs[c("exp(coef)", "exp(-coef)", "lower .95", "upper .95")], function(x) log(x))
  # Return data frame
  return(Coefs)
}

### Get all coefs in a data frame
# Canopy Height
ssfCoefs.canopyHeight <- get_Coefs(1, "Canopy Height", 4)
if (!("lower" %in% names(ssfCoefs.canopyHeight))) {
  ssfCoefs.canopyHeight <- ssfCoefs.canopyHeight %>%
    mutate(
      lower = coalesce(!!!rlang::syms(c("lower.95", "lower .95"))),
      upper = coalesce(!!!rlang::syms(c("upper.95", "upper .95")))
    )
}

# Create effect category: Positive (both CI bounds > 0),
# Negative (both < 0), otherwise No effect.
ssfCoefs.canopyHeight <- ssfCoefs.canopyHeight %>%
  mutate(
    Sign = case_when(
      lower > 0 & upper > 0 ~ "Positive",
      lower < 0 & upper < 0 ~ "Negative",
      TRUE                  ~ "No effect"
    ),
    Sign = factor(Sign, levels = c("Positive","Negative","No effect"))
  )

# VCI
ssfCoefs.VCI <- get_Coefs(3, "VCI", 4)
if (!("lower" %in% names(ssfCoefs.VCI))) {
  ssfCoefs.VCI <- ssfCoefs.VCI %>%
    mutate(
      lower = coalesce(!!!rlang::syms(c("lower.95", "lower .95"))),
      upper = coalesce(!!!rlang::syms(c("upper.95", "upper .95")))
    )
}

# Create effect category: Positive (both CI bounds > 0),
# Negative (both < 0), otherwise No effect.
ssfCoefs.VCI <- ssfCoefs.VCI %>%
  mutate(
    Sign = case_when(
      lower > 0 & upper > 0 ~ "Positive",
      lower < 0 & upper < 0 ~ "Negative",
      TRUE                  ~ "No effect"
    ),
    Sign = factor(Sign, levels = c("Positive","Negative","No effect"))
  )

# Dist to gap 50 m2
ssfCoefs.d50 <- get_Coefs(4, "Dist. to Gap 50 m", 4)
if (!("lower" %in% names(ssfCoefs.d50))) {
  ssfCoefs.d50 <- ssfCoefs.d50 %>%
    mutate(
      lower = coalesce(!!!rlang::syms(c("lower.95", "lower .95"))),
      upper = coalesce(!!!rlang::syms(c("upper.95", "upper .95")))
    )
}

# Create effect category: Positive (both CI bounds > 0),
# Negative (both < 0), otherwise No effect.
ssfCoefs.d50 <- ssfCoefs.d50 %>%
  mutate(
    Sign = case_when(
      lower > 0 & upper > 0 ~ "Positive",
      lower < 0 & upper < 0 ~ "Negative",
      TRUE                  ~ "No effect"
    ),
    Sign = factor(Sign, levels = c("Positive","Negative","No effect"))
  )

# Dist. to gap 500 m2
ssfCoefs.d500 <- get_Coefs(5, "Dist. to Gap 500 m", 4)
if (!("lower" %in% names(ssfCoefs.d500))) {
  ssfCoefs.d500 <- ssfCoefs.d500 %>%
    mutate(
      lower = coalesce(!!!rlang::syms(c("lower.95", "lower .95"))),
      upper = coalesce(!!!rlang::syms(c("upper.95", "upper .95")))
    )
}

# Create effect category: Positive (both CI bounds > 0),
# Negative (both < 0), otherwise No effect.
ssfCoefs.d500 <- ssfCoefs.d500 %>%
  mutate(
    Sign = case_when(
      lower > 0 & upper > 0 ~ "Positive",
      lower < 0 & upper < 0 ~ "Negative",
      TRUE                  ~ "No effect"
    ),
    Sign = factor(Sign, levels = c("Positive","Negative","No effect"))
  )

# PVD
ssfCoefs.pvd <- get_Coefs(6, "PVD 15-20", 4)
if (!("lower" %in% names(ssfCoefs.pvd))) {
  ssfCoefs.pvd <- ssfCoefs.pvd %>%
    mutate(
      lower = coalesce(!!!rlang::syms(c("lower.95", "lower .95"))),
      upper = coalesce(!!!rlang::syms(c("upper.95", "upper .95")))
    )
}

# Create effect category: Positive (both CI bounds > 0),
# Negative (both < 0), otherwise No effect.
ssfCoefs.pvd <- ssfCoefs.pvd %>%
  mutate(
    Sign = case_when(
      lower > 0 & upper > 0 ~ "Positive",
      lower < 0 & upper < 0 ~ "Negative",
      TRUE                  ~ "No effect"
    ),
    Sign = factor(Sign, levels = c("Positive","Negative","No effect"))
  )

# Swamp
ssfCoefs.swamp <- get_CoefsSwamp(7, "Swamp", 2)
if (!("lower" %in% names(ssfCoefs.swamp))) {
  ssfCoefs.swamp <- ssfCoefs.swamp %>%
    mutate(
      lower = coalesce(!!!rlang::syms(c("lower.95", "lower .95"))),
      upper = coalesce(!!!rlang::syms(c("upper.95", "upper .95")))
    )
}

# Create effect category: Positive (both CI bounds > 0),
# Negative (both < 0), otherwise No effect.
ssfCoefs.swamp <- ssfCoefs.swamp %>%
  mutate(
    Sign = case_when(
      lower > 0 & upper > 0 ~ "Positive",
      lower < 0 & upper < 0 ~ "Negative",
      TRUE                  ~ "No effect"
    ),
    Sign = factor(Sign, levels = c("Positive","Negative","No effect"))
  )



### Make SSF plots for 2 hr intervals
## Plot with NO x and y axis labels
ssfPlots <- function(df){
  # Arrange data frame by coefficients
  level_order <- df %>% arrange(coef)
  level_order <- level_order$ID
  
  # Make the plot
  plot <- ggplot(df, aes(x = factor(ID, level = level_order), y = coef, color = Sign)) +
    geom_pointrange(aes(ymin = lower.95, ymax = upper.95),
                    position = position_dodge(width = 0.7), size = 0.8) +
    geom_hline(yintercept = 0, lty = 4) +
    scale_color_manual(
      values = c("Positive" = "forestgreen",
                 "Negative" = "firebrick2",
                 "No effect" = "black")
    ) +
    labs(x = "", y="") +
    coord_flip() +
    theme_light(base_size = 20)
}

# Make the plots
ch.plot <- ssfPlots(ssfCoefs.canopyHeight)
vc.plot <- ssfPlots(ssfCoefs.VCI)
d50.plot <- ssfPlots(ssfCoefs.d50)
d500.plot <- ssfPlots(ssfCoefs.d500)
pvd15_20.plot <- ssfPlots(ssfCoefs.pvd)
swamp.plot <- ssfPlots(ssfCoefs.swamp)

#Define consistent legend across panels
effect_levels <- c("Positive","Negative","No effect")
effect_colors <- c("Positive"="forestgreen",
                   "Negative"="firebrick2",
                   "No effect"="black")

effect_scale <- scale_color_manual(
  values = effect_colors,
  limits = effect_levels,   # enforce all 3 levels in legend
  breaks = effect_levels,
  drop   = FALSE,
  name   = NULL
)

# Make legend keys show points (not errorbar lines), even if a level is absent:
effect_guides <- guides(
  color = guide_legend(
    override.aes = list(shape = 16, size = 1.2, linetype = 0)
  )
)
# acc effect scalse
ch.plot      <- ch.plot      + effect_scale + effect_guides + theme(legend.position = "bottom")
vc.plot      <- vc.plot      + effect_scale + effect_guides + theme(legend.position = "none")
d50.plot     <- d50.plot     + effect_scale + effect_guides + theme(legend.position = "none")
d500.plot    <- d500.plot    + effect_scale + effect_guides + theme(legend.position = "none")
pvd15_20.plot<- pvd15_20.plot+ effect_scale + effect_guides + theme(legend.position = "none")
swamp.plot   <- swamp.plot   + effect_scale + effect_guides + theme(legend.position = "none")

# Plots all together
# All on one plot
windows()
p_grid <- ggarrange(ch.plot + ggtitle(bquote(~bold('A)')~'Canopy Height')) + theme(legend.position="none", plot.title = element_text(size = 18)),
          vc.plot + ggtitle(bquote(~bold('B)')~'Vertical Complexity Index')) + theme(legend.position = "none", plot.title = element_text(size = 18)), 
          d50.plot + ggtitle(bquote(~bold('C)')~'Dist. to Small Gaps (50'~m^2~")")) + theme(legend.position = "none", plot.title = element_text(size = 18)), 
          d500.plot + ggtitle(bquote(~bold('D)')~'Dist. to Large Gaps (500'~m^2~")")) + theme(legend.position = "none", plot.title = element_text(size = 18)),
          pvd15_20.plot + ggtitle(bquote(~bold('E)')~'Plant Vol. Dens. (15-20 m)')) + theme(legend.position = "none", plot.title = element_text(size = 18)),
          swamp.plot + ggtitle(bquote(~bold('E)')~'Swamp Habitat')) + theme(legend.position = "none", plot.title = element_text(size = 18)),
          common.legend = TRUE, legend = "bottom",
          nrow=3,
          ncol=2)

p_annot<- annotate_figure(
  p_grid,
  left   = text_grob("Turaco ID", rot = 90, size = 24), 
  bottom = text_grob(expression(beta ~ "coefficient"), size = 24)
)

# 3) Extract a legend from one plot (the one with legend visible)
leg <- cowplot::get_legend(ch.plot)  # uses the bottom legend we set

# 4) Stack the annotated grid OVER the legend (legend ends up below x-label)
final <- cowplot::plot_grid(
  p_annot,
  cowplot::ggdraw(leg),
  ncol = 1,
  rel_heights = c(1, 0.10)  # tweak legend height as needed
)

final

##############################################
## Canopy Height Relative selection strength
##############################################
s1 <- data.frame(
  CanopyHeight10m_end = seq(from = -2, to = 2, length.out = 100),
  VCI_end = 0,
  dist2gap50_end = 0,
  dist2gap500_end = 0,
  pvd15_20_end = 0,
  Swamp_end = 0,
  sl_ = log(100),
  ta_ = 1)

s2 <- data.frame(
  CanopyHeight10m_end = 0,
  VCI_end = 0,
  dist2gap50_end = 0,
  dist2gap500_end = 0,
  pvd15_20_end = 0,
  Swamp_end = 0,
  sl_ = log(100),
  ta_ = 1)


# Calculate log-RSS for all birds
lr_93 <- log_rss(fit8893, s1, s2)
lr_95 <- log_rss(fit8895, s1, s2)
lr_51 <- log_rss(fit11851, s1, s2)
lr_53 <- log_rss(fit11853, s1, s2)

# Create data frame with RSS predictions for all birds
rss <- cbind(lr_93$df$CanopyHeight10m_end_x1, 
             lr_93$df$log_rss,
             lr_95$df$log_rss,
             lr_51$df$log_rss,
             lr_53$df$log_rss)
rss <- as.data.frame(rss)
colnames(rss) <- c('CanopyHeight', 'T8893', 'T8895', 'T11851', 'T11853')
head(rss)

# Colors for lines
colors <- c("T8893" = "#440154FF", 
            "T8895" = "#3B528BFF",
            "T11851" = "#5DC863FF",
            "T11853" = "#FDE725FF")

# Plot the RSS figure
windows()
RSS <- ggplot(rss, aes(x = CanopyHeight)) +
  geom_line(size = 2, aes(y=T8893, color = "T8893")) +
  geom_line(size = 2, aes(y=T8895, color = "T8895")) +
  geom_line(size = 2, aes(y=T11851, color = "T11851")) +
  geom_line(size = 2, aes(y=T11853, color = "T11853")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  labs(x = "Canopy Height (Scaled)",
       y = "log-RSS vs Mean Canopy Height",
       color = "Legend") +
  scale_color_manual(name = "Turaco ID", values = colors) +
  theme_light(base_size = 20) +
  theme(legend.position = "bottom")
RSS

##############################
## R session info
##############################
sessionInfo()
