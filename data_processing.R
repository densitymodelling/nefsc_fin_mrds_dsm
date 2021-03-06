# get the data into a dsm-friendly format
# data provided by Doug Sigourney, NOAA NEFSC.
# based on data from the paper at https://peerj.com/articles/8226/


# Load data files
# Sightings data formatted for mrds
mrds_Sightings <- read.csv("data/Fin_Sightings_mrds_for_dsm.csv")
# segment data
Fin_Whale_Data <- read.csv("data/Final_Fin_Whale_Data_for_dsm.csv")


# rename columns, data munging
# rescale distances
mrds_Sightings$distance <- mrds_Sightings$distance/1000
# remove -9999 observer
mrds_Sightings <- subset(mrds_Sightings, observer != -9999)

# get column names
mrds_Sightings$Sample.Label <- mrds_Sightings$segment
mrds_Sightings$segment <- NULL
Fin_Whale_Data$Sample.Label <- Fin_Whale_Data$segment
Fin_Whale_Data$segment <- NULL

Fin_Whale_Data$Sample.Label <- paste0(Fin_Whale_Data$Survey, "-", Fin_Whale_Data$Sample.Label)
mrds_Sightings$Sample.Label <- paste0(mrds_Sightings$survey, "-", mrds_Sightings$Sample.Label)

# identifier for which detectioin function this effort is for
Fin_Whale_Data$ddfobj <- Fin_Whale_Data$Survey
Fin_Whale_Data$Survey <- NULL

# same detection function covariate names in both data sets
Fin_Whale_Data$SUBJ_WAVG <- Fin_Whale_Data$Subjective
Fin_Whale_Data$Subjective <- NULL
Fin_Whale_Data$beaufort <- Fin_Whale_Data$Beaufort
Fin_Whale_Data$Beaufort <- NULL


# Species_Hab_Num is species identifier:
#   1   fin whale
#   2   fin or sei
#   3   sei whale
# we only want to fit the spatial model to the fin data only
fin_segs <- subset(Fin_Whale_Data, Species_Hab_Num==1)
fin_segs <- fin_segs[, c("Sample.Label", "Effort", "DIST125", "DEPTH",
                         "DIST2SHORE", "SST", "ddfobj", "SUBJ_WAVG",
                         "beaufort")]

# truncation distances
w_ship <- 6
w_plane <- 0.9


# don't need all columns in detection data
mrds_Sightings <- mrds_Sightings[, c("object", "observer", "detected",
                                     "distance", "survey", "beaufort",
                                     "SUBJ_WAVG", "size", "Sample.Label",
                                     "Species_Hab_Num")]

# shipboard survey data
ship_detections <- subset(mrds_Sightings,
                          survey == 1 &          # ship detections
                          distance <= w_ship)  # within truncation

# aerial survey data
# all sightings from the front team
plane_detections <- subset(mrds_Sightings,
                           survey == 2 &          # plane detections
                           distance <= w_plane &  # within truncation
                           observer == 1 &        # front team only
                           detected == 1)         # and detected!



## we're going to fit detection functions with all detections (fin, fin/sei, sei)
# but only build the spatial model with fin detections
# build observation data.frame that reflects this
fin_obs <- mrds_Sightings

## filter only fin detections
fin_obs <- fin_obs[fin_obs$Species_Hab_Num == 1, ]
fin_obs$Species_Hab_Num <- NULL
# only want one of each observation in the obs frame
fin_obs <- fin_obs[!duplicated(fin_obs$object), ]
fin_obs$ddfobj <- fin_obs$survey
## only get observations within truncation
# note this is different from what is in the Sigourney et al
# supplementary materials (but this is being corrected!)
fin_obs <- subset(fin_obs, (ddfobj == 1 & distance < w_ship) |
                           (ddfobj == 2 & distance < w_plane))

# select only the columns we want
fin_obs <- fin_obs[, c("object", "Sample.Label", "observer", "ddfobj", "size", "distance")]



# Import data set with values for all grid cells
All_Covs <- read.csv("data/Mean_Summer_Covariates_By_Grid_Cell_Final.csv")
All_Covs <- subset(All_Covs, DIST125!='NA' & DEPTH!='NA'
                   & DIST2SHORE!='NA' & SST!='NA')
predgrid <- All_Covs[, c("Area", "DIST125", "DEPTH",
                         "DIST2SHORE", "SST", "LAT", "LON")]
predgrid$off.set <- predgrid$Area
predgrid$Area <- NULL


# project lat/long for plotting
library(sf)
proj <- "+proj=omerc +lat_0=35 +lonc=-75 +alpha=40 +k=0.9996 +x_0=0 +y_0=0 +gamma=40 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
coords <- st_sfc(st_multipoint(as.matrix(predgrid[,c("LON", "LAT")])))
st_crs(coords) <- 4326
coords <- st_coordinates(st_transform(coords,
                                      crs=proj))
predgrid$x <- coords[,1]
predgrid$y <- coords[,2]

# save all data
save(ship_detections, plane_detections, fin_obs, fin_segs, predgrid,
     file="findata.RData")
