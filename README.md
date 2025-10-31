# Habitat selection of great blue turacos and resulting patterns of seed dispersal

Nicholas J. Russo and Antoine S.A. Tekam

nicholasrusso@g.harvard.edu

This repository contains the data and code necessary to reproduce the results of Tekam et al. (DOI: )

### Scripts
- `turaco_iSSA_share_GitHub.R`: Script that produces all analyses and figures included in the manuscript, including:
	- Selection for all structural covariates based on Integrated Step Selection Analyses 		(iSSAs; Fig. 2)
	- Relative Selection Strength for Canopy Height (Fig. 3)
	- Seed shadow simulations (Fig. 4)
 - `turaco_UHC_share_GitHub.R`: Script that produces used-habitat calibration plots used to validate iSSA models and produce Figures S2-S5 in the manuscript.
- `turaco_iSSA_share_Revised.R`: Script that performs all the functions of `turaco_iSSA_share_GitHub.R` but with an updated Fig. 2 that colors points according to the sign of the relationship

### Data folder
## Turaco Data
- `Turaco_GPS_16Dec2024.csv`: Locations of four great blue turacos in Movebank format. 
Relevant columns:                       
  	- "timestamp": time in GMT
	- "location.long": longitude
	- "location.lat": latitude
	- "tag.local.identifier": e-obs tag number
	- "individual.local.identifier": same as e-obs tag number (one tag per individual)
	- "utm.easting": easting
	- "utm.northing": northing
	- "study.local.timestamp": time in West Africa Time
- `ch.tif`: Canopy Height, 10 m resolution (site level)
- `vc.tif`: Vertical Complexity Index, 10 m resolution (site level)
- `d50.tif`: Distance to gap of size >= 50 m, 10 m resolution (site level)
- `d500.tif`: Distance to gap of size >=500 m, 10 m resolution (site level)
- `pvd_15_20.tif`: Plant Volume Density 15-20 m from the ground, 10 m resolution (site level)
- `swamp_forest.tif`: Landscape classified as either swamp (1) or other habitat (0) 


################
Session Info:

R version 4.4.2 (2024-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 22631)

Matrix products: default


locale:
[1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
[3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.utf8    

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] viridis_0.6.5          viridisLite_0.4.2      ggpubr_0.6.0          
 [4] spatstat_3.2-1         spatstat.linnet_3.2-2  spatstat.model_3.3-2  
 [7] rpart_4.1.23           spatstat.explore_3.3-3 nlme_3.1-166          
[10] spatstat.random_3.3-2  spatstat.geom_3.3-4    spatstat.univar_3.1-1 
[13] spatstat.data_3.1-4    sf_1.0-19              terra_1.7-83          
[16] ctmm_1.2.0             adehabitatHR_0.4.22    adehabitatLT_0.3.28   
[19] CircStats_0.2-6        boot_1.3-31            MASS_7.3-61           
[22] adehabitatMA_0.3.17    ade4_1.7-22            move_4.2.5            
[25] raster_3.6-30          sp_2.1-4               geosphere_1.5-20      
[28] amt_0.2.2.0            dplyr_1.1.4            lubridate_1.9.3       
[31] ggplot2_3.5.1   
