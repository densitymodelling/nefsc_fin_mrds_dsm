# Using double observer methods with density surface models
# based on code from Doug Sigourney https://peerj.com/articles/8226/
# data provided by Doug Sigourney, NOAA NEFSC.

library(mrds)
# need at least dsm version 2.3.1.9007 from github!
library(dsm)
library(Distance)

# Load data files
# Sightings data formatted for mrds
mrds_Sightings <- read.csv("Fin_Sightings_mrds_for_dsm.csv")
Fin_Whale_Data <- read.csv("Final_Fin_Whale_Data_for_dsm.csv")

# data formatting for DSM
mrds_Sightings$distance <- mrds_Sightings$distance/1000
mrds_Sightings$Sample.Label <- mrds_Sightings$segment
Fin_Whale_Data$Sample.Label <- Fin_Whale_Data$segment
Fin_Whale_Data$ddfobj <- Fin_Whale_Data$Survey

# is this right?
Fin_Whale_Data$SUBJ_WAVG <- Fin_Whale_Data$Subjective
Fin_Whale_Data$beaufort <- Fin_Whale_Data$Beaufort

# species identifier: 1=Fin, 2=fin/sei, 3=sei
Final_Cov_Mat_Fin <- subset(Fin_Whale_Data, Species_Hab_Num==1)


# Truncation Distances (1=Shipboard, 2=Aerial...here and throughout)
W_1 <- 6
W_2 <- 0.9

# shipboard survey data
All_Dat_Ship <- subset(mrds_Sightings, survey==1 & distance>-1 & distance < W_1)

# aerial survey data
All_Dat_Plane <- subset(mrds_Sightings, survey==2 & distance>-1 & distance < W_2)
# all sightings from the front team
All_Dat_Plane_1 <- subset(All_Dat_Plane, observer==1 & detected==1)

# build observation data.frame
obs <- rbind(All_Dat_Ship, All_Dat_Plane_1)

## model fitting

# ship surveys - MRDS
Ship.mrds <- ddf(dsmodel = ~mcds(key = "hr", formula = ~beaufort + SUBJ_WAVG),
                 mrmodel = ~glm(~distance),
                 data = All_Dat_Ship, method = "io",
                 meta.data = list(width = W_1))
summary(Ship.mrds)

# aerial surveys - MCDS
Plane.ds <- ddf(dsmodel = ~mcds(key = "hr", formula = ~beaufort),
                data = All_Dat_Plane_1, meta.data = list(width = W_2))
summary(Plane.ds)

# might need remove sei, fin/sei observations

# Fit DSM
# this is as in Sigourney et al, setting group sizes to be 1
# then inflating with average group size later
Best.Model <- dsm(count ~ s(DIST125,bs="ts", k=5) +
                          s(DEPTH,bs="ts", k=5) +
                          s(DIST2SHORE,bs="ts", k=5) +
                          s(SST,bs="ts", k=5),
                  ddf.obj=list(Ship.mrds, Plane.ds),
                  family=tw(link="log"), method="REML",
                  segment.data=Final_Cov_Mat_Fin, group=TRUE,
                  observation.data=obs)

obs_exp(Best.Model, "beaufort", c(0, 1, 2, 3, 4, 5))



# Make predictions

# Import data set with values for all grid cells
All_Covs <- read.csv("Mean_Summer_Covariates_By_Grid_Cell_Final.csv")
All_Covs <- subset(All_Covs, DIST125!='NA' & DEPTH!='NA'
                   & DIST2SHORE!='NA' & SST!='NA')
pred_mat <- All_Covs[, c("grid", "Area", "DIST125", "DEPTH",
                         "DIST2SHORE", "SST")]
pred_mat$off.set <- pred_mat$Area


#Make predictions using predict function
predict.fun <- predict(Best.Model, newdata=pred_mat, off.set=pred_mat$Area)

# without correction
sum(predict.fun)

# mutiply each prediction by estimate of average group size
sum(predict.fun*1.38)


# same model but using observed rather than estimated group sizes
group.model <- dsm(count ~ s(DIST125,bs="ts", k=5) +
                           s(DEPTH,bs="ts", k=5) +
                           s(DIST2SHORE,bs="ts", k=5) +
                           s(SST,bs="ts", k=5),
                   ddf.obj=list(Ship.mrds, Plane.ds),
                   family=tw(link="log"), method="REML",
                   segment.data=Final_Cov_Mat_Fin,
                   observation.data=obs)

obs_exp(group.model, "beaufort", c(0, 1, 2, 3, 4, 5))
# better?
sum(predict(group.model, newdata=pred_mat, off.set=pred_mat$Area))

# propagate uncertainty
vp <- dsm_varprop(group.model, pred_mat)
vp

