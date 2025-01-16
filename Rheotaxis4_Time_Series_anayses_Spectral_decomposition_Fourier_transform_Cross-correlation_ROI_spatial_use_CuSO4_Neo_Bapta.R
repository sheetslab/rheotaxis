#START HERE FOR TIME SERIES ANALYSES

library(tidyverse)
library(dplyr)
library(ggplot2)
#library(broom)
#library(lubridate)
#library(forcats)
#library(tsibble)

library(mgcv)
library(MASS)

library(splines)
library(zoo)
library(tsModel)
library(gganimate)
library(spectral)
library(circular)


#all these files should have been run through rheotaxis_file_prep_concatenate.R to unify format
# master.all1 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/CuNeo/master.all.CuNeo.csv") #CuSO4 v Neomycin - this one is fine?
master.all2 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/Old_shake/master.all.shake.csv") #Shake - no recovery! - this one is fine?
master.all3 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/60shake/master.data.shake.rec.csv") #SHAKE - RECOVERY 
# master.all4 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/60bapta_test//master.data.bapta.noKS.csv") #BAPTA
# master.all5 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/sixty_flow_test/no_KS/master.data.60s.CuNeo.csv") #60s flow test - CuSO4 v Neomycin
# master.all6 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/sixty_flow_test/no_KS/master.data.60s.shake.csv") #60s flow test - Shake

master.all <- master.all3

#make sure that this data is a factor and not simply a string of characters
master.all$Rec_Time <- as.factor(master.all$Rec_Time) #shake-recovery only
master.all$Space_bin <- as.factor(master.all$Space_bin) 
master.all$Individual <- as.factor(master.all$Individual)
master.all$Treatment <- as.factor(master.all$Treatment)
master.all$Stimulus <- as.factor(master.all$Stimulus)
str(master.all) # check and fix data categorization 

#master.all <- master.all <- master.all[, c(2:19)] # if necessary, remove extra first column added when imported file - WITHOUT SPACE BIN
master.all <- master.all[, c(2:21)] # if necessary, remove extra first column added when imported file - WITH SPACE BIN


#these are used in TIME SERIES 
#no Space_bin for omnibus analyses (others?)
rheo_agg1 <- aggregate(list(master.all[ , c(1,3,5,6,9:13,17,18)]), by = list(master.all$Rheotaxis, master.all$Treatment, master.all$Time), mean, na.rm=TRUE) #Is this an issue with NAs? does it skip the calculation and return NA for entire aggregate if one NA is present???
names (rheo_agg1)[names(rheo_agg1) == "Group.1"] <- "Rheotaxis"
names (rheo_agg1)[names(rheo_agg1) == "Group.2"] <- "Treatment"
names (rheo_agg1)[names(rheo_agg1) == "Group.3"] <- "Time"
rheo_agg1 <- rheo_agg1[ , c(3,2,6:13,1)] #delete extra columns and rearrange
rheo_agg1$Fish_mean_angle_rad <- atan2(rheo_agg1$Fish_angle_sin,rheo_agg1$Fish_angle_cos) #calculate mean body angle (radians) from mean sine and cosine of angle
rheo_agg1$Fish_mean_angle_deg <- deg(rheo_agg1$Fish_mean_angle_rad) #convert radian to degrees - good sanity check
rheo_agg1$Fish_mean_angle_deg <- ifelse(rheo_agg1$Fish_mean_angle_deg<0, rheo_agg1$Fish_mean_angle_deg +360, rheo_agg1$Fish_mean_angle_deg) # NEED TO FIX +/- 180 to 0-360
rheo_agg1$Mean_angle_dispersion <- sqrt((rheo_agg1$Fish_angle_sin^2)+(rheo_agg1$Fish_angle_cos^2))
rheo_agg1a <- rheo_agg1 #change name of data frame to save multiple version for different data sets


# rheo_agg2 <- aggregate(list(master.all3[ , c(1,3,5,6,9:13,17,18,20)]), by = list(master.all3$Rheotaxis, master.all3$Treatment, master.all3$Time, master.all3$Rec_Time), mean, na.rm=TRUE) #Is this an issue with NAs? does it skip the calculation and return NA for entire aggregate if one NA is present???
# names (rheo_agg2)[names(rheo_agg2) == "Group.1"] <- "Rheotaxis"
# names (rheo_agg2)[names(rheo_agg2) == "Group.2"] <- "Treatment"
# names (rheo_agg2)[names(rheo_agg2) == "Group.3"] <- "Time"
# names (rheo_agg2)[names(rheo_agg2) == "Group.4"] <- "Rec_Time"
# rheo_agg2 <- rheo_agg2[ , c(3,2,4,7:14,1)] #delete extra columns and rearrange
# rheo_agg2$Fish_mean_angle_rad <- atan2(rheo_agg2$Fish_angle_sin,rheo_agg2$Fish_angle_cos) #calculate mean body angle (radians) from mean sine and cosine of angle
# rheo_agg2$Fish_mean_angle_deg <- deg(rheo_agg2$Fish_mean_angle_rad) #convert radian to degrees - good sanity check
# rheo_agg2$Fish_mean_angle_deg <- ifelse(rheo_agg2$Fish_mean_angle_deg<0, rheo_agg2$Fish_mean_angle_deg +360, rheo_agg2$Fish_mean_angle_deg) # NEED TO FIX +/- 180 to 0-360
# rheo_agg2$Mean_angle_dispersion <- sqrt((rheo_agg2$Fish_angle_sin^2)+(rheo_agg2$Fish_angle_cos^2))
# rheo_agg2a <- rheo_agg2 #change name of data frame to save multiple version for different data sets


# DECOMPOSE TIME SERIES ---####
#intermediate plots generated are helpful in understanding decomposition process

#Control
#decompose
rheo_agg1_CTL <- filter(rheo_agg1, rheo_agg1$Treatment=="Control")
rheo_agg1_CTL_rheo <- filter(rheo_agg1_CTL, rheo_agg1_CTL$Rheotaxis==1)
#rheo_agg1_CTL_rheo  <- filter(rheo_agg1_CTL_rheo, rheo_agg1_CTL_rheo$Time<=30000)
ts_CTL_rheo_agg1_SBmov <- zoo(rheo_agg1_CTL_rheo$SB_mov, rheo_agg1_CTL_rheo$Time)
ts_CTL_rheo_agg1_SBvel <- zoo(rheo_agg1_CTL_rheo$SB_vel, rheo_agg1_CTL_rheo$Time)
ts_CTL_rheo_agg1_SBacc <- zoo(rheo_agg1_CTL_rheo$SB_acc, rheo_agg1_CTL_rheo$Time)
ts_CTL_rheo_agg1_angle_rad <- zoo(rheo_agg1_CTL_rheo$Fish_mean_angle_rad, rheo_agg1_CTL_rheo$Time)
ts_CTL_rheo_agg1_angle_dispersion <- zoo(rheo_agg1_CTL_rheo$Mean_angle_dispersion, rheo_agg1_CTL_rheo$Time)

ts_CTL_rheo_agg1_SBmov1 <- ts(ts_CTL_rheo_agg1_SBmov, frequency = 10)
plot(ts_CTL_rheo_agg1_SBmov1)
ts_CTL_rheo_agg1_SBvel1 <- ts(ts_CTL_rheo_agg1_SBvel, frequency = 10)
ts_CTL_rheo_agg1_SBvel1 <- na.locf(ts_CTL_rheo_agg1_SBvel1) #must account for NAs!!!
plot(ts_CTL_rheo_agg1_SBvel1)
ts_CTL_rheo_agg1_SBacc1 <- ts(ts_CTL_rheo_agg1_SBacc, frequency = 10)
ts_CTL_rheo_agg1_SBacc1 <- na.locf(ts_CTL_rheo_agg1_SBacc1) #must account for NAs!!!
plot(ts_CTL_rheo_agg1_SBacc1)
ts_CTL_rheo_agg1_angle_rad1 <- ts(ts_CTL_rheo_agg1_angle_rad, frequency = 10)
plot(ts_CTL_rheo_agg1_angle_rad1)
ts_CTL_rheo_agg1_angle_dispersion1 <- ts(ts_CTL_rheo_agg1_angle_dispersion, frequency = 10)
plot(ts_CTL_rheo_agg1_angle_dispersion1)

plot(cbind(ts_CTL_rheo_agg1_SBmov1, ts_CTL_rheo_agg1_angle_rad1))
plot(cbind(ts_CTL_rheo_agg1_SBmov1, ts_CTL_rheo_agg1_angle_dispersion1))
plot(cbind(ts_CTL_rheo_agg1_SBvel1, ts_CTL_rheo_agg1_angle_rad1))
plot(cbind(ts_CTL_rheo_agg1_SBvel1, ts_CTL_rheo_agg1_angle_dispersion1))
plot(cbind(ts_CTL_rheo_agg1_SBacc1, ts_CTL_rheo_agg1_angle_rad1))
plot(cbind(ts_CTL_rheo_agg1_SBacc1, ts_CTL_rheo_agg1_angle_dispersion1))
plot(cbind(ts_CTL_rheo_agg1_angle_rad1, ts_CTL_rheo_agg1_angle_dispersion1))

#cross-correlation
ctl1 <- ccf(ts_CTL_rheo_agg1_SBmov1, ts_CTL_rheo_agg1_angle_rad1, na.action = na.pass) #must account for one NA value
ctl2 <- ccf(ts_CTL_rheo_agg1_SBmov1, ts_CTL_rheo_agg1_angle_dispersion1, na.action = na.pass)
ctl3 <- ccf(ts_CTL_rheo_agg1_SBvel1, ts_CTL_rheo_agg1_angle_rad1, na.action = na.pass)
ctl4 <- ccf(ts_CTL_rheo_agg1_SBvel1, ts_CTL_rheo_agg1_angle_dispersion1, na.action = na.pass)
ctl5 <- ccf(ts_CTL_rheo_agg1_SBacc1, ts_CTL_rheo_agg1_angle_rad1, na.action = na.pass)
ctl6 <- ccf(ts_CTL_rheo_agg1_SBacc1, ts_CTL_rheo_agg1_angle_dispersion1, na.action = na.pass)
# 
# 
# #CuSO4
# #decompose
# rheo_agg1_Cu <- filter(rheo_agg1, rheo_agg1$Treatment=="CuSO4")
# rheo_agg1_Cu_rheo <- filter(rheo_agg1_Cu, rheo_agg1_Cu$Rheotaxis==1)
# #rheo_agg1_Cu_rheo  <- filter(rheo_agg1_Cu_rheo, rheo_agg1_Cu_rheo$Time<=30000)
# ts_Cu_rheo_agg1_SBmov <- zoo(rheo_agg1_Cu_rheo$SB_mov, rheo_agg1_Cu_rheo$Time)
# ts_Cu_rheo_agg1_SBvel <- zoo(rheo_agg1_Cu_rheo$SB_vel, rheo_agg1_Cu_rheo$Time)
# ts_Cu_rheo_agg1_SBacc <- zoo(rheo_agg1_Cu_rheo$SB_acc, rheo_agg1_Cu_rheo$Time)
# ts_Cu_rheo_agg1_angle_rad <- zoo(rheo_agg1_Cu_rheo$Fish_mean_angle_rad, rheo_agg1_Cu_rheo$Time)
# ts_Cu_rheo_agg1_angle_dispersion <- zoo(rheo_agg1_Cu_rheo$Mean_angle_dispersion, rheo_agg1_Cu_rheo$Time)
# 
# ts_Cu_rheo_agg1_SBmov1 <- ts(ts_Cu_rheo_agg1_SBmov, frequency = 10)
# plot(ts_Cu_rheo_agg1_SBmov1)
# ts_Cu_rheo_agg1_SBvel1 <- ts(ts_Cu_rheo_agg1_SBvel, frequency = 10)
# plot(ts_Cu_rheo_agg1_SBvel1)
# ts_Cu_rheo_agg1_SBacc1 <- ts(ts_Cu_rheo_agg1_SBacc, frequency = 10)
# plot(ts_Cu_rheo_agg1_SBacc1)
# ts_Cu_rheo_agg1_angle_rad1 <- ts(ts_Cu_rheo_agg1_angle_rad, frequency = 10)
# plot(ts_Cu_rheo_agg1_angle_rad1)
# ts_Cu_rheo_agg1_angle_dispersion1 <- ts(ts_Cu_rheo_agg1_angle_dispersion, frequency = 10)
# plot(ts_Cu_rheo_agg1_angle_dispersion1)
# 
# plot(cbind(ts_Cu_rheo_agg1_SBmov1, ts_Cu_rheo_agg1_angle_rad1))
# plot(cbind(ts_Cu_rheo_agg1_SBmov1, ts_Cu_rheo_agg1_angle_dispersion1))
# plot(cbind(ts_Cu_rheo_agg1_SBvel1, ts_Cu_rheo_agg1_angle_rad1))
# plot(cbind(ts_Cu_rheo_agg1_SBvel1, ts_Cu_rheo_agg1_angle_dispersion1))
# plot(cbind(ts_Cu_rheo_agg1_SBacc1, ts_Cu_rheo_agg1_angle_rad1))
# plot(cbind(ts_Cu_rheo_agg1_SBacc1, ts_Cu_rheo_agg1_angle_dispersion1))
# plot(cbind(ts_Cu_rheo_agg1_angle_rad1, ts_Cu_rheo_agg1_angle_dispersion1))
# 
# #cross-correlationm
# cu1 <- ccf(ts_Cu_rheo_agg1_SBmov1, ts_Cu_rheo_agg1_angle_rad1, na.action = na.pass) #must account for one NA value
# cu2 <- ccf(ts_Cu_rheo_agg1_SBmov1, ts_Cu_rheo_agg1_angle_dispersion1, na.action = na.pass)
# cu3 <- ccf(ts_Cu_rheo_agg1_SBvel1, ts_Cu_rheo_agg1_angle_rad1, na.action = na.pass)
# cu4 <- ccf(ts_Cu_rheo_agg1_SBvel1, ts_Cu_rheo_agg1_angle_dispersion1, na.action = na.pass)
# cu5 <- ccf(ts_Cu_rheo_agg1_SBacc1, ts_Cu_rheo_agg1_angle_rad1, na.action = na.pass)
# cu6 <- ccf(ts_Cu_rheo_agg1_SBacc1, ts_Cu_rheo_agg1_angle_dispersion1, na.action = na.pass)
# 
# 
# #Neomycin
# #decompose
# rheo_agg1_Neo <- filter(rheo_agg1, rheo_agg1$Treatment=="Neo")
# rheo_agg1_Neo_rheo <- filter(rheo_agg1_Neo, rheo_agg1_Neo$Rheotaxis==1)
# #rheo_agg1_Neo_rheo  <- filter(rheo_agg1_Neo_rheo, rheo_agg1_Neo_rheo$Time<=30000)
# ts_Neo_rheo_agg1_SBmov <- zoo(rheo_agg1_Neo_rheo$SB_mov, rheo_agg1_Neo_rheo$Time)
# ts_Neo_rheo_agg1_SBvel <- zoo(rheo_agg1_Neo_rheo$SB_vel, rheo_agg1_Neo_rheo$Time)
# ts_Neo_rheo_agg1_SBacc <- zoo(rheo_agg1_Neo_rheo$SB_acc, rheo_agg1_Neo_rheo$Time)
# ts_Neo_rheo_agg1_angle_rad <- zoo(rheo_agg1_Neo_rheo$Fish_mean_angle_rad, rheo_agg1_Neo_rheo$Time)
# ts_Neo_rheo_agg1_angle_dispersion <- zoo(rheo_agg1_Neo_rheo$Mean_angle_dispersion, rheo_agg1_Neo_rheo$Time)
# 
# ts_Neo_rheo_agg1_SBmov1 <- ts(ts_Neo_rheo_agg1_SBmov, frequency = 10)
# plot(ts_Neo_rheo_agg1_SBmov1)
# ts_Neo_rheo_agg1_SBvel1 <- ts(ts_Neo_rheo_agg1_SBvel, frequency = 10)
# plot(ts_Neo_rheo_agg1_SBvel1)
# ts_Neo_rheo_agg1_SBacc1 <- ts(ts_Neo_rheo_agg1_SBacc, frequency = 10)
# plot(ts_Neo_rheo_agg1_SBacc1)
# ts_Neo_rheo_agg1_angle_rad1 <- ts(ts_Neo_rheo_agg1_angle_rad, frequency = 10)
# plot(ts_Neo_rheo_agg1_angle_rad1)
# ts_Neo_rheo_agg1_angle_dispersion1 <- ts(ts_Neo_rheo_agg1_angle_dispersion, frequency = 10)
# plot(ts_Neo_rheo_agg1_angle_dispersion1)
# 
# plot(cbind(ts_Neo_rheo_agg1_SBmov1, ts_Neo_rheo_agg1_angle_rad1))
# plot(cbind(ts_Neo_rheo_agg1_SBmov1, ts_Neo_rheo_agg1_angle_dispersion1))
# plot(cbind(ts_Neo_rheo_agg1_SBvel1, ts_Neo_rheo_agg1_angle_rad1))
# plot(cbind(ts_Neo_rheo_agg1_SBvel1, ts_Neo_rheo_agg1_angle_dispersion1))
# plot(cbind(ts_Neo_rheo_agg1_SBacc1, ts_Neo_rheo_agg1_angle_rad1))
# plot(cbind(ts_Neo_rheo_agg1_SBacc1, ts_Neo_rheo_agg1_angle_dispersion1))
# plot(cbind(ts_Neo_rheo_agg1_angle_rad1, ts_Neo_rheo_agg1_angle_dispersion1))
# 
# #Cross-correlation
# neo1 <- ccf(ts_Neo_rheo_agg1_SBmov1, ts_Neo_rheo_agg1_angle_rad1, na.action = na.pass) #must account for one NA value
# neo2 <- ccf(ts_Neo_rheo_agg1_SBmov1, ts_Neo_rheo_agg1_angle_dispersion1, na.action = na.pass) 
# neo3 <- ccf(ts_Neo_rheo_agg1_SBvel1, ts_Neo_rheo_agg1_angle_rad1, na.action = na.pass) 
# neo4 <- ccf(ts_Neo_rheo_agg1_SBvel1, ts_Neo_rheo_agg1_angle_dispersion1, na.action = na.pass)
# neo5 <- ccf(ts_Neo_rheo_agg1_SBacc1, ts_Neo_rheo_agg1_angle_rad1, na.action = na.pass) 
# neo6 <- ccf(ts_Neo_rheo_agg1_SBacc1, ts_Neo_rheo_agg1_angle_dispersion1, na.action = na.pass) 


#BAPTA
#decompose
# rheo_agg1_bap <- filter(rheo_agg1, rheo_agg1$Treatment=="Bapta")
# rheo_agg1_bap_rheo <- filter(rheo_agg1_bap, rheo_agg1_bap$Rheotaxis==1)
# #rheo_agg1_bap_rheo  <- filter(rheo_agg1_bap_rheo, rheo_agg1_bap_rheo$Time<=30000)
# ts_bap_rheo_agg1_SBmov <- zoo(rheo_agg1_bap_rheo$SB_mov, rheo_agg1_bap_rheo$Time)
# ts_bap_rheo_agg1_SBvel <- zoo(rheo_agg1_bap_rheo$SB_vel, rheo_agg1_bap_rheo$Time)
# ts_bap_rheo_agg1_SBacc <- zoo(rheo_agg1_bap_rheo$SB_acc, rheo_agg1_bap_rheo$Time)
# ts_bap_rheo_agg1_angle_rad <- zoo(rheo_agg1_bap_rheo$Fish_mean_angle_rad, rheo_agg1_bap_rheo$Time)
# ts_bap_rheo_agg1_angle_dispersion <- zoo(rheo_agg1_bap_rheo$Mean_angle_dispersion, rheo_agg1_bap_rheo$Time)
# 
# ts_bap_rheo_agg1_SBmov1 <- ts(ts_bap_rheo_agg1_SBmov, frequency = 10)
# plot(ts_bap_rheo_agg1_SBmov1)
# ts_bap_rheo_agg1_SBvel1 <- ts(ts_bap_rheo_agg1_SBvel, frequency = 10)
# plot(ts_bap_rheo_agg1_SBvel1)
# ts_bap_rheo_agg1_SBacc1 <- ts(ts_bap_rheo_agg1_SBacc, frequency = 10)
# plot(ts_bap_rheo_agg1_SBacc1)
# ts_bap_rheo_agg1_angle_rad1 <- ts(ts_bap_rheo_agg1_angle_rad, frequency = 10)
# plot(ts_bap_rheo_agg1_angle_rad1)
# ts_bap_rheo_agg1_angle_dispersion1 <- ts(ts_bap_rheo_agg1_angle_dispersion, frequency = 10)
# plot(ts_bap_rheo_agg1_angle_dispersion1)
# 
# plot(cbind(ts_bap_rheo_agg1_SBmov1, ts_bap_rheo_agg1_angle_rad1))
# plot(cbind(ts_bap_rheo_agg1_SBmov1, ts_bap_rheo_agg1_angle_dispersion1))
# plot(cbind(ts_bap_rheo_agg1_SBvel1, ts_bap_rheo_agg1_angle_rad1))
# plot(cbind(ts_bap_rheo_agg1_SBvel1, ts_bap_rheo_agg1_angle_dispersion1))
# plot(cbind(ts_bap_rheo_agg1_SBacc1, ts_bap_rheo_agg1_angle_rad1))
# plot(cbind(ts_bap_rheo_agg1_SBacc1, ts_bap_rheo_agg1_angle_dispersion1))
# plot(cbind(ts_bap_rheo_agg1_angle_rad1, ts_bap_rheo_agg1_angle_dispersion1))
# 
# #cross-correlation
# bap1 <- ccf(ts_bap_rheo_agg1_SBmov1, ts_bap_rheo_agg1_angle_rad1, na.action = na.pass) #must account for one NA value
# bap2 <- ccf(ts_bap_rheo_agg1_SBmov1, ts_bap_rheo_agg1_angle_dispersion1, na.action = na.pass)
# bap3 <- ccf(ts_bap_rheo_agg1_SBvel1, ts_bap_rheo_agg1_angle_rad1, na.action = na.pass)
# bap4 <- ccf(ts_bap_rheo_agg1_SBvel1, ts_bap_rheo_agg1_angle_dispersion1, na.action = na.pass)
# bap5 <- ccf(ts_bap_rheo_agg1_SBacc1, ts_bap_rheo_agg1_angle_rad1, na.action = na.pass)
# bap6 <- ccf(ts_bap_rheo_agg1_SBacc1, ts_bap_rheo_agg1_angle_dispersion1, na.action = na.pass)



#SHAKE
#decompose
rheo_agg1_shk <- filter(rheo_agg1, rheo_agg1$Treatment=="Shake")
rheo_agg1_shk_rheo <- filter(rheo_agg1_shk, rheo_agg1_shk$Rheotaxis==1)
#rheo_agg1_shk_rheo  <- filter(rheo_agg1_shk_rheo, rheo_agg1_shk_rheo$Time<=30000)
ts_shk_rheo_agg1_SBmov <- zoo(rheo_agg1_shk_rheo$SB_mov, rheo_agg1_shk_rheo$Time)
ts_shk_rheo_agg1_SBvel <- zoo(rheo_agg1_shk_rheo$SB_vel, rheo_agg1_shk_rheo$Time)
ts_shk_rheo_agg1_SBacc <- zoo(rheo_agg1_shk_rheo$SB_acc, rheo_agg1_shk_rheo$Time)
ts_shk_rheo_agg1_angle_rad <- zoo(rheo_agg1_shk_rheo$Fish_mean_angle_rad, rheo_agg1_shk_rheo$Time)
ts_shk_rheo_agg1_angle_dispersion <- zoo(rheo_agg1_shk_rheo$Mean_angle_dispersion, rheo_agg1_shk_rheo$Time)

ts_shk_rheo_agg1_SBmov1 <- ts(ts_shk_rheo_agg1_SBmov, frequency = 10)
plot(ts_shk_rheo_agg1_SBmov1)
ts_shk_rheo_agg1_SBvel1 <- ts(ts_shk_rheo_agg1_SBvel, frequency = 10)
ts_shk_rheo_agg1_SBvel1 <- na.locf(ts_shk_rheo_agg1_SBvel1) #must account for NAs!!!
plot(ts_shk_rheo_agg1_SBvel1)
ts_shk_rheo_agg1_SBacc1 <- ts(ts_shk_rheo_agg1_SBacc, frequency = 10)
ts_shk_rheo_agg1_SBacc1 <- na.locf(ts_shk_rheo_agg1_SBacc1) #must account for NAs!!!
plot(ts_shk_rheo_agg1_SBacc1)
ts_shk_rheo_agg1_angle_rad1 <- ts(ts_shk_rheo_agg1_angle_rad, frequency = 10)
plot(ts_shk_rheo_agg1_angle_rad1)
ts_shk_rheo_agg1_angle_dispersion1 <- ts(ts_shk_rheo_agg1_angle_dispersion, frequency = 10)
plot(ts_shk_rheo_agg1_angle_dispersion1)

plot(cbind(ts_shk_rheo_agg1_SBmov1, ts_shk_rheo_agg1_angle_rad1))
plot(cbind(ts_shk_rheo_agg1_SBmov1, ts_shk_rheo_agg1_angle_dispersion1))
plot(cbind(ts_shk_rheo_agg1_SBvel1, ts_shk_rheo_agg1_angle_rad1))
plot(cbind(ts_shk_rheo_agg1_SBvel1, ts_shk_rheo_agg1_angle_dispersion1))
plot(cbind(ts_shk_rheo_agg1_SBacc1, ts_shk_rheo_agg1_angle_rad1))
plot(cbind(ts_shk_rheo_agg1_SBacc1, ts_shk_rheo_agg1_angle_dispersion1))
plot(cbind(ts_shk_rheo_agg1_angle_rad1, ts_shk_rheo_agg1_angle_dispersion1))

#cross-correlation
shk1 <- ccf(ts_shk_rheo_agg1_SBmov1, ts_shk_rheo_agg1_angle_rad1, na.action = na.pass) #must account for one NA value
shk2 <- ccf(ts_shk_rheo_agg1_SBmov1, ts_shk_rheo_agg1_angle_dispersion1, na.action = na.pass)
shk3 <- ccf(ts_shk_rheo_agg1_SBvel1, ts_shk_rheo_agg1_angle_rad1, na.action = na.pass)
shk4 <- ccf(ts_shk_rheo_agg1_SBvel1, ts_shk_rheo_agg1_angle_dispersion1, na.action = na.pass)
shk5 <- ccf(ts_shk_rheo_agg1_SBacc1, ts_shk_rheo_agg1_angle_rad1, na.action = na.pass)
shk6 <- ccf(ts_shk_rheo_agg1_SBacc1, ts_shk_rheo_agg1_angle_dispersion1, na.action = na.pass)



#SHAKE-RECOVERY

#Rec_time subsets 
rheo_agg2_0hr <- filter(rheo_agg2, rheo_agg2$Rec_Time=="0hr")
rheo_agg2_2hr <- filter(rheo_agg2, rheo_agg2$Rec_Time=="2hr")
rheo_agg2_4hr <- filter(rheo_agg2, rheo_agg2$Rec_Time=="4hr")
rheo_agg2_8hr <- filter(rheo_agg2, rheo_agg2$Rec_Time=="8hr")
rheo_agg2_48hr <- filter(rheo_agg2, rheo_agg2$Rec_Time=="48hr")


#Control
#decompose
rheo_agg2_CTL8 <- filter(rheo_agg2_8hr, rheo_agg2_8hr$Treatment=="Control") #adjust rec-time 
rheo_agg2_CTL8_rheo <- filter(rheo_agg2_CTL8, rheo_agg2_CTL8$Rheotaxis==1) #adjust rec-time 

#rheo_agg2_CTL8_rheo  <- filter(rheo_agg2_CTL8_rheo, rheo_agg2_CTL8_rheo$Time<=30000)
ts_CTL8_rheo_agg2_SBmov <- zoo(rheo_agg2_CTL8_rheo$SB_mov, rheo_agg2_CTL8_rheo$Time)
ts_CTL8_rheo_agg2_SBvel <- zoo(rheo_agg2_CTL8_rheo$SB_vel, rheo_agg2_CTL8_rheo$Time)
ts_CTL8_rheo_agg2_SBacc <- zoo(rheo_agg2_CTL8_rheo$SB_acc, rheo_agg2_CTL8_rheo$Time)
ts_CTL8_rheo_agg2_angle_rad <- zoo(rheo_agg2_CTL8_rheo$Fish_mean_angle_rad, rheo_agg2_CTL8_rheo$Time)
ts_CTL8_rheo_agg2_angle_dispersion <- zoo(rheo_agg2_CTL8_rheo$Mean_angle_dispersion, rheo_agg2_CTL8_rheo$Time)

ts_CTL8_rheo_agg2_SBmov1 <- ts(ts_CTL8_rheo_agg2_SBmov, frequency = 10)
ts_CTL8_rheo_agg2_SBmov1 <- na.locf(ts_CTL8_rheo_agg2_SBmov1) #must account for NAs!!!
plot(ts_CTL8_rheo_agg2_SBmov1)
ts_CTL8_rheo_agg2_SBvel1 <- ts(ts_CTL8_rheo_agg2_SBvel, frequency = 10)
ts_CTL8_rheo_agg2_SBvel1 <- na.locf(ts_CTL8_rheo_agg2_SBvel1) #must account for NAs!!!
plot(ts_CTL8_rheo_agg2_SBvel1)
ts_CTL8_rheo_agg2_SBacc1 <- ts(ts_CTL8_rheo_agg2_SBacc, frequency = 10)
ts_CTL8_rheo_agg2_SBacc1 <- na.locf(ts_CTL8_rheo_agg2_SBacc1) #must account for NAs!!!
plot(ts_CTL8_rheo_agg2_SBacc1)
ts_CTL8_rheo_agg2_angle_rad1 <- ts(ts_CTL8_rheo_agg2_angle_rad, frequency = 10)
ts_CTL8_rheo_agg2_angle_rad1 <- na.locf(ts_CTL8_rheo_agg2_angle_rad1) #must account for NAs!!!
plot(ts_CTL8_rheo_agg2_angle_rad1)
ts_CTL8_rheo_agg2_angle_dispersion1 <- ts(ts_CTL8_rheo_agg2_angle_dispersion, frequency = 10)
ts_CTL8_rheo_agg2_angle_dispersion1 <- na.locf(ts_CTL8_rheo_agg2_angle_dispersion1) #must account for NAs!!!
plot(ts_CTL8_rheo_agg2_angle_dispersion1)

plot(cbind(ts_CTL8_rheo_agg2_SBmov1, ts_CTL8_rheo_agg2_angle_rad1))
plot(cbind(ts_CTL8_rheo_agg2_SBmov1, ts_CTL8_rheo_agg2_angle_dispersion1))
plot(cbind(ts_CTL8_rheo_agg2_SBvel1, ts_CTL8_rheo_agg2_angle_rad1))
plot(cbind(ts_CTL8_rheo_agg2_SBvel1, ts_CTL8_rheo_agg2_angle_dispersion1))
plot(cbind(ts_CTL8_rheo_agg2_SBacc1, ts_CTL8_rheo_agg2_angle_rad1))
plot(cbind(ts_CTL8_rheo_agg2_SBacc1, ts_CTL8_rheo_agg2_angle_dispersion1))
plot(cbind(ts_CTL8_rheo_agg2_angle_rad1, ts_CTL8_rheo_agg2_angle_dispersion1))

#cross-correlation
ctl8_1 <- ccf(ts_CTL8_rheo_agg2_SBmov1, ts_CTL8_rheo_agg2_angle_rad1, na.action = na.pass) #must account for one NA value
ctl8_2 <- ccf(ts_CTL8_rheo_agg2_SBmov1, ts_CTL8_rheo_agg2_angle_dispersion1, na.action = na.pass)
ctl8_3 <- ccf(ts_CTL8_rheo_agg2_SBvel1, ts_CTL8_rheo_agg2_angle_rad1, na.action = na.pass)
ctl8_4 <- ccf(ts_CTL8_rheo_agg2_SBvel1, ts_CTL8_rheo_agg2_angle_dispersion1, na.action = na.pass)
ctl8_5 <- ccf(ts_CTL8_rheo_agg2_SBacc1, ts_CTL8_rheo_agg2_angle_rad1, na.action = na.pass)
ctl8_6 <- ccf(ts_CTL8_rheo_agg2_SBacc1, ts_CTL8_rheo_agg2_angle_dispersion1, na.action = na.pass)



#shake
#decompose
rheo_agg2_shk8 <- filter(rheo_agg2_8hr, rheo_agg2_8hr$Treatment=="Shake") #adjust rec-time 
rheo_agg2_shk8_rheo <- filter(rheo_agg2_shk8, rheo_agg2_shk8$Rheotaxis==1) #adjust rec-time 
#rheo_agg2_shk8_rheo  <- filter(rheo_agg2_shk8_rheo, rheo_agg2_shk8_rheo$Time<=30000)
ts_shk8_rheo_agg2_SBmov <- zoo(rheo_agg2_shk8_rheo$SB_mov, rheo_agg2_shk8_rheo$Time)
ts_shk8_rheo_agg2_SBvel <- zoo(rheo_agg2_shk8_rheo$SB_vel, rheo_agg2_shk8_rheo$Time)
ts_shk8_rheo_agg2_SBacc <- zoo(rheo_agg2_shk8_rheo$SB_acc, rheo_agg2_shk8_rheo$Time)
ts_shk8_rheo_agg2_angle_rad <- zoo(rheo_agg2_shk8_rheo$Fish_mean_angle_rad, rheo_agg2_shk8_rheo$Time)
ts_shk8_rheo_agg2_angle_dispersion <- zoo(rheo_agg2_shk8_rheo$Mean_angle_dispersion, rheo_agg2_shk8_rheo$Time)

ts_shk8_rheo_agg2_SBmov1 <- ts(ts_shk8_rheo_agg2_SBmov, frequency = 10)
ts_shk8_rheo_agg2_SBmov1 <- na.locf(ts_shk8_rheo_agg2_SBmov1)  #must account for NAs!!!
plot(ts_shk8_rheo_agg2_SBmov1)
ts_shk8_rheo_agg2_SBvel1 <- ts(ts_shk8_rheo_agg2_SBvel, frequency = 10)
ts_shk8_rheo_agg2_SBvel1 <- na.locf(ts_shk8_rheo_agg2_SBvel1)  #must account for NAs!!!
plot(ts_shk8_rheo_agg2_SBvel1)
ts_shk8_rheo_agg2_SBacc1 <- ts(ts_shk8_rheo_agg2_SBacc, frequency = 10)
ts_shk8_rheo_agg2_SBacc1 <- na.locf(ts_shk8_rheo_agg2_SBacc1)  #must account for NAs!!!
plot(ts_shk8_rheo_agg2_SBacc1)
ts_shk8_rheo_agg2_angle_rad1 <- ts(ts_shk8_rheo_agg2_angle_rad, frequency = 10)
ts_shk8_rheo_agg2_angle_rad1 <- na.locf(ts_shk8_rheo_agg2_angle_rad1)  #must account for NAs!!!
plot(ts_shk8_rheo_agg2_angle_rad1)
ts_shk8_rheo_agg2_angle_dispersion1 <- ts(ts_shk8_rheo_agg2_angle_dispersion, frequency = 10)
ts_shk8_rheo_agg2_angle_dispersion1 <- na.locf(ts_shk8_rheo_agg2_angle_dispersion1)  #must account for NAs!!!
plot(ts_shk8_rheo_agg2_angle_dispersion1)

plot(cbind(ts_shk8_rheo_agg2_SBmov1, ts_shk8_rheo_agg2_angle_rad1))
plot(cbind(ts_shk8_rheo_agg2_SBmov1, ts_shk8_rheo_agg2_angle_dispersion1))
plot(cbind(ts_shk8_rheo_agg2_SBvel1, ts_shk8_rheo_agg2_angle_rad1))
plot(cbind(ts_shk8_rheo_agg2_SBvel1, ts_shk8_rheo_agg2_angle_dispersion1))
plot(cbind(ts_shk8_rheo_agg2_SBacc1, ts_shk8_rheo_agg2_angle_rad1))
plot(cbind(ts_shk8_rheo_agg2_SBacc1, ts_shk8_rheo_agg2_angle_dispersion1))
plot(cbind(ts_shk8_rheo_agg2_angle_rad1, ts_shk8_rheo_agg2_angle_dispersion1))

#cross-correlation
shk8_1 <- ccf(ts_shk8_rheo_agg2_SBmov1, ts_shk8_rheo_agg2_angle_rad1, na.action = na.pass) #must account for one NA value
shk8_2 <- ccf(ts_shk8_rheo_agg2_SBmov1, ts_shk8_rheo_agg2_angle_dispersion1, na.action = na.pass)
shk8_3 <- ccf(ts_shk8_rheo_agg2_SBvel1, ts_shk8_rheo_agg2_angle_rad1, na.action = na.pass)
shk8_4 <- ccf(ts_shk8_rheo_agg2_SBvel1, ts_shk8_rheo_agg2_angle_dispersion1, na.action = na.pass)
shk8_5 <- ccf(ts_shk8_rheo_agg2_SBacc1, ts_shk8_rheo_agg2_angle_rad1, na.action = na.pass)
shk8_6 <- ccf(ts_shk8_rheo_agg2_SBacc1, ts_shk8_rheo_agg2_angle_dispersion1, na.action = na.pass)




####---FIGURES - DECOMPOSITION OF TIME SERIES into: observed, TREND, SEASONALITY/PERIODICITY, noise ---####
plot(decompose(ts_CTL_rheo_agg1_SBmov1))
plot(decompose(ts_CTL_rheo_agg1_SBvel1))
plot(decompose(ts_CTL_rheo_agg1_SBacc1))
plot(decompose(ts_CTL_rheo_agg1_angle_rad1))
plot(decompose(ts_CTL_rheo_agg1_angle_dispersion1))

# ts_Cu_rheo_agg1_SBmov2 <- na.locf(ts_Cu_rheo_agg1_SBmov1) #must account for NAs???
# ts_Cu_rheo_agg1_SBvel2 <- na.locf(ts_Cu_rheo_agg1_SBvel1)
# ts_Cu_rheo_agg1_SBacc2 <- na.locf(ts_Cu_rheo_agg1_SBacc1)
# ts_Cu_rheo_agg1_angle_rad2 <- na.locf(ts_Cu_rheo_agg1_angle_rad1)
# ts_Cu_rheo_agg1_angle_dispersion2 <- na.locf(ts_Cu_rheo_agg1_angle_dispersion1)

# plot(decompose(ts_Cu_rheo_agg1_SBmov1)) #nope, you can run it
# plot(decompose(ts_Cu_rheo_agg1_SBvel1))
# plot(decompose(ts_Cu_rheo_agg1_SBacc1))
# plot(decompose(ts_Cu_rheo_agg1_angle_rad1))
# plot(decompose(ts_Cu_rheo_agg1_angle_dispersion1))
# 
# plot(decompose(ts_Neo_rheo_agg1_SBmov1))
# plot(decompose(ts_Neo_rheo_agg1_SBvel1))
# plot(decompose(ts_Neo_rheo_agg1_SBacc1))
# plot(decompose(ts_Neo_rheo_agg1_angle_rad1))
# plot(decompose(ts_Neo_rheo_agg1_angle_dispersion1))

# plot(decompose(ts_bap_rheo_agg1_SBmov1)) 
# plot(decompose(ts_bap_rheo_agg1_SBvel1))
# plot(decompose(ts_bap_rheo_agg1_SBacc1))
# plot(decompose(ts_bap_rheo_agg1_angle_rad1))
# plot(decompose(ts_bap_rheo_agg1_angle_dispersion1))

plot(decompose(ts_shk_rheo_agg1_SBmov1))
plot(decompose(ts_shk_rheo_agg1_SBvel1))
plot(decompose(ts_shk_rheo_agg1_SBacc1))
plot(decompose(ts_shk_rheo_agg1_angle_rad1))
plot(decompose(ts_shk_rheo_agg1_angle_dispersion1))

# plot(decompose(ts_CTL8_rheo_agg2_SBmov1))
# plot(decompose(ts_CTL8_rheo_agg2_SBvel1))
# plot(decompose(ts_CTL8_rheo_agg2_SBacc1))
# plot(decompose(ts_CTL8_rheo_agg2_angle_rad1))
# plot(decompose(ts_CTL8_rheo_agg2_angle_dispersion1))
# 
# plot(decompose(ts_shk8_rheo_agg2_SBmov1)) 
# plot(decompose(ts_shk8_rheo_agg2_SBvel1))
# plot(decompose(ts_shk8_rheo_agg2_SBacc1))
# plot(decompose(ts_shk8_rheo_agg2_angle_rad1))
# plot(decompose(ts_shk8_rheo_agg2_angle_dispersion1))


#SB movement
# plot(decompose(ts_CTL_rheo_agg1_SBmov1))
# 
deco_ts_CTL_rheo_agg1_SBmov1 <- decompose(ts_CTL_rheo_agg1_SBmov1)
obs_CTL_mov <- deco_ts_CTL_rheo_agg1_SBmov1$x
trd_CTL_mov  <- deco_ts_CTL_rheo_agg1_SBmov1$trend
seas_CTL_mov  <- deco_ts_CTL_rheo_agg1_SBmov1$seasonal
ran_CTL_mov  <- deco_ts_CTL_rheo_agg1_SBmov1$random

# deco_ts_Cu_rheo_agg1_SBmov2 <- decompose(ts_Cu_rheo_agg1_SBmov2)
# obs_Cu_mov <- deco_ts_Cu_rheo_agg1_SBmov2$x
# trd_Cu_mov  <- deco_ts_Cu_rheo_agg1_SBmov2$trend
# seas_Cu_mov  <- deco_ts_Cu_rheo_agg1_SBmov2$seasonal
# ran_Cu_mov  <- deco_ts_Cu_rheo_agg1_SBmov2$random
# 
# deco_ts_Neo_rheo_agg1_SBmov1 <- decompose(ts_Neo_rheo_agg1_SBmov1)
# obs_Neo_mov <- deco_ts_Neo_rheo_agg1_SBmov1$x
# trd_Neo_mov  <- deco_ts_Neo_rheo_agg1_SBmov1$trend
# seas_Neo_mov  <- deco_ts_Neo_rheo_agg1_SBmov1$seasonal
# ran_Neo_mov  <- deco_ts_Neo_rheo_agg1_SBmov1$random

# deco_ts_bap_rheo_agg1_SBmov1 <- decompose(ts_bap_rheo_agg1_SBmov1)
# obs_bap_mov <- deco_ts_bap_rheo_agg1_SBmov1$x
# trd_bap_mov  <- deco_ts_bap_rheo_agg1_SBmov1$trend
# seas_bap_mov  <- deco_ts_bap_rheo_agg1_SBmov1$seasonal
# ran_bap_mov  <- deco_ts_bap_rheo_agg1_SBmov1$random

deco_ts_shk_rheo_agg1_SBmov1 <- decompose(ts_shk_rheo_agg1_SBmov1)
obs_shk_mov <- deco_ts_shk_rheo_agg1_SBmov1$x
trd_shk_mov  <- deco_ts_shk_rheo_agg1_SBmov1$trend
seas_shk_mov  <- deco_ts_shk_rheo_agg1_SBmov1$seasonal
ran_shk_mov  <- deco_ts_shk_rheo_agg1_SBmov1$random

# deco_ts_CTL8_rheo_agg2_SBmov1 <- decompose(ts_CTL8_rheo_agg2_SBmov1)
# obs_CTL8_mov <- deco_ts_CTL8_rheo_agg2_SBmov1$x
# trd_CTL8_mov  <- deco_ts_CTL8_rheo_agg2_SBmov1$trend
# seas_CTL8_mov  <- deco_ts_CTL8_rheo_agg2_SBmov1$seasonal
# ran_CTL8_mov  <- deco_ts_CTL8_rheo_agg2_SBmov1$random
# 
# deco_ts_shk8_rheo_agg2_SBmov1 <- decompose(ts_shk8_rheo_agg2_SBmov1)
# obs_shk8_mov <- deco_ts_shk8_rheo_agg2_SBmov1$x
# trd_shk8_mov  <- deco_ts_shk8_rheo_agg2_SBmov1$trend
# seas_shk8_mov  <- deco_ts_shk8_rheo_agg2_SBmov1$seasonal
# ran_shk8_mov  <- deco_ts_shk8_rheo_agg2_SBmov1$random


# all_obs_mov <- cbind(obs_Cu_mov, obs_Neo_mov, obs_CTL_mov) #this order helps with graph layer ordering bottom -> top :Cu-blue, Neo-gr, Ctl-gray
# all_trd_mov <- cbind(trd_Cu_mov, trd_Neo_mov, trd_CTL_mov)
# all_seas_mov <- cbind(seas_Cu_mov, seas_Neo_mov, seas_CTL_mov)
# all_ran_mov <- cbind(ran_Cu_mov, ran_Neo_mov, ran_CTL_mov)
# 
# all_obs_mov <- cbind(obs_bap_mov, obs_CTL_mov) #this order helps with graph layer ordering bottom -> top :bap-blue, Neo-gr, Ctl-gray
# all_trd_mov <- cbind(trd_bap_mov, trd_CTL_mov)
# all_seas_mov <- cbind(seas_bap_mov, seas_CTL_mov)
# all_ran_mov <- cbind(ran_bap_mov, ran_CTL_mov)

all_obs_mov <- cbind(obs_shk_mov, obs_CTL_mov) #this order helps with graph layer ordering bottom -> top :shk-blue, Neo-gr, Ctl-gray
all_trd_mov <- cbind(trd_shk_mov, trd_CTL_mov)
all_seas_mov <- cbind(seas_shk_mov, seas_CTL_mov)
all_ran_mov <- cbind(ran_shk_mov, ran_CTL_mov)

# all_obs_mov <- cbind(obs_shk8_mov, obs_CTL8_mov) #this order helps with graph layer ordering bottom -> top :shk-blue, Neo-gr, Ctl-gray
# all_trd_mov <- cbind(trd_shk8_mov, trd_CTL8_mov)
# all_seas_mov <- cbind(seas_shk8_mov, seas_CTL8_mov)
# all_ran_mov <- cbind(ran_shk8_mov, ran_CTL8_mov)

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/SpecDeco_SBmov_obs_shake.pdf", width = 5, height = 4)
plot.ts(all_obs_mov,
        plot.type = "single", 
        # col=c("blue", "green", "gray"),
        col=c("red", "gray"),
        ylab="SB distance moved (mm)",
        main="Observed",
        xlim=c(0,40),
        ylim=c(0,5)
)
abline(v=10, col="gray", lty=5)
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/SpecDeco_SBmov_trd_shake.pdf", width = 5, height = 4)
plot.ts(all_trd_mov,
        plot.type = "single", 
        # col=c("blue", "green", "gray"),
        col=c("red", "gray"),
        ylab="SB distance moved (mm)",
        main="Overall Trend",
        xlim=c(0,40),
        ylim=c(0,2.5)
        )
abline(v=10, col="gray", lty=5)
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/SpecDeco_SBmov_seas_shake.pdf", width = 5, height = 4)
plot.ts(all_seas_mov,
        plot.type = "single", 
        # col=c("blue", "green", "gray"),
        col=c("red", "gray"),
        ylab="SB distance moved (mm)",
        main="Seasonal Fluctuation",
        xlim=c(0,5),
        ylim=c(-0.1,0.15)
)
abline(v=10, col="gray", lty=5)
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/SpecDeco_SBmov_ran_shake.pdf", width = 5, height = 4)
plot.ts(all_ran_mov,
        plot.type = "single", 
        # col=c("blue", "green", "gray"),
        col=c("red", "gray"),
        ylab="SB distance moved (mm)",
        main="Random Noise",
        xlim=c(0,40),
        ylim=c(-2,3)
)
abline(v=10, col="gray", lty=5)
dev.off()


#SB velocity
# plot(decompose(ts_CTL_rheo_agg1_SBvel1))
# 
deco_ts_CTL_rheo_agg1_SBvel1 <- decompose(ts_CTL_rheo_agg1_SBvel1)
obs_CTL_vel <- deco_ts_CTL_rheo_agg1_SBvel1$x
trd_CTL_vel  <- deco_ts_CTL_rheo_agg1_SBvel1$trend
seas_CTL_vel  <- deco_ts_CTL_rheo_agg1_SBvel1$seasonal
ran_CTL_vel  <- deco_ts_CTL_rheo_agg1_SBvel1$random

# deco_ts_Cu_rheo_agg1_SBvel2 <- decompose(ts_Cu_rheo_agg1_SBvel2)
# obs_Cu_vel <- deco_ts_Cu_rheo_agg1_SBvel2$x
# trd_Cu_vel  <- deco_ts_Cu_rheo_agg1_SBvel2$trend
# seas_Cu_vel  <- deco_ts_Cu_rheo_agg1_SBvel2$seasonal
# ran_Cu_vel  <- deco_ts_Cu_rheo_agg1_SBvel2$random
# 
# deco_ts_Neo_rheo_agg1_SBvel1 <- decompose(ts_Neo_rheo_agg1_SBvel1)
# obs_Neo_vel <- deco_ts_Neo_rheo_agg1_SBvel1$x
# trd_Neo_vel  <- deco_ts_Neo_rheo_agg1_SBvel1$trend
# seas_Neo_vel  <- deco_ts_Neo_rheo_agg1_SBvel1$seasonal
# ran_Neo_vel  <- deco_ts_Neo_rheo_agg1_SBvel1$random

# deco_ts_bap_rheo_agg1_SBvel1 <- decompose(ts_bap_rheo_agg1_SBvel1)
# obs_bap_vel <- deco_ts_bap_rheo_agg1_SBvel1$x
# trd_bap_vel  <- deco_ts_bap_rheo_agg1_SBvel1$trend
# seas_bap_vel  <- deco_ts_bap_rheo_agg1_SBvel1$seasonal
# ran_bap_vel  <- deco_ts_bap_rheo_agg1_SBvel1$random

deco_ts_shk_rheo_agg1_SBvel1 <- decompose(ts_shk_rheo_agg1_SBvel1)
obs_shk_vel <- deco_ts_shk_rheo_agg1_SBvel1$x
trd_shk_vel  <- deco_ts_shk_rheo_agg1_SBvel1$trend
seas_shk_vel  <- deco_ts_shk_rheo_agg1_SBvel1$seasonal
ran_shk_vel  <- deco_ts_shk_rheo_agg1_SBvel1$random

# deco_ts_CTL8_rheo_agg2_SBvel1 <- decompose(ts_CTL8_rheo_agg2_SBvel1)
# obs_CTL8_vel <- deco_ts_CTL8_rheo_agg2_SBvel1$x
# trd_CTL8_vel  <- deco_ts_CTL8_rheo_agg2_SBvel1$trend
# seas_CTL8_vel  <- deco_ts_CTL8_rheo_agg2_SBvel1$seasonal
# ran_CTL8_vel  <- deco_ts_CTL8_rheo_agg2_SBvel1$random
# 
# deco_ts_shk8_rheo_agg2_SBvel1 <- decompose(ts_shk8_rheo_agg2_SBvel1)
# obs_shk8_vel <- deco_ts_shk8_rheo_agg2_SBvel1$x
# trd_shk8_vel  <- deco_ts_shk8_rheo_agg2_SBvel1$trend
# seas_shk8_vel  <- deco_ts_shk8_rheo_agg2_SBvel1$seasonal
# ran_shk8_vel  <- deco_ts_shk8_rheo_agg2_SBvel1$random

# all_obs_vel <- cbind(obs_Cu_vel, obs_Neo_vel, obs_CTL_vel) #this order helps with graph layer ordering bottom -> top :Cu-blue, Neo-gr, Ctl-gray
# all_trd_vel <- cbind(trd_Cu_vel, trd_Neo_vel, trd_CTL_vel)
# all_seas_vel <- cbind(seas_Cu_vel, seas_Neo_vel, seas_CTL_vel)
# all_ran_vel <- cbind(ran_Cu_vel, ran_Neo_vel, ran_CTL_vel)

# all_obs_vel <- cbind(obs_bap_vel, obs_CTL_vel) #this order helps with graph layer ordering bottom -> top :bap-blue, Neo-gr, Ctl-gray
# all_trd_vel <- cbind(trd_bap_vel, trd_CTL_vel)
# all_seas_vel <- cbind(seas_bap_vel, seas_CTL_vel)
# all_ran_vel <- cbind(ran_bap_vel, ran_CTL_vel)

all_obs_vel <- cbind(obs_shk_vel, obs_CTL_vel) #this order helps with graph layer ordering bottom -> top :shk-blue, Neo-gr, Ctl-gray
all_trd_vel <- cbind(trd_shk_vel, trd_CTL_vel)
all_seas_vel <- cbind(seas_shk_vel, seas_CTL_vel)
all_ran_vel <- cbind(ran_shk_vel, ran_CTL_vel)

# all_obs_vel <- cbind(obs_shk8_vel, obs_CTL8_vel) #this order helps with graph layer ordering bottom -> top :shk-blue, Neo-gr, Ctl-gray
# all_trd_vel <- cbind(trd_shk8_vel, trd_CTL8_vel)
# all_seas_vel <- cbind(seas_shk8_vel, seas_CTL8_vel)
# all_ran_vel <- cbind(ran_shk8_vel, ran_CTL8_vel)

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/SpecDeco_SBvel_obs_shake.pdf", width = 5, height = 4)
plot.ts(all_obs_vel,
        plot.type = "single", 
        # col=c("blue", "green", "gray"),
        col=c("red", "gray"),
        ylab="SB velocity (mm/sec)",
        main="Observed",
        xlim=c(0,40),
        ylim=c(-35,30)
)
abline(v=10, col="gray", lty=5)
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/SpecDeco_SBvel_trd_shake.pdf", width = 5, height = 4)
plot.ts(all_trd_vel,
        plot.type = "single", 
        # col=c("blue", "green", "gray"),
        col=c("red", "gray"),
        ylab="SB velocity (mm/sec)",
        main="Overall Trend",
        xlim=c(0,40),
        ylim=c(-7,6)
)
abline(v=10, col="gray", lty=5)
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/SpecDeco_SBvel_seas_shake.pdf", width = 5, height = 4)
plot.ts(all_seas_vel,
        plot.type = "single", 
        # col=c("blue", "green", "gray"),
        col=c("red", "gray"),
        ylab="SB velocity (mm/sec)",
        main="Seasonal Fluctuation",
        xlim=c(0,5),
        ylim=c(-1.5,2)
)
abline(v=10, col="gray", lty=5)
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/SpecDeco_SBvel_ran_shake.pdf", width = 5, height = 4)
plot.ts(all_ran_vel,
        plot.type = "single", 
        # col=c("blue", "green", "gray"),
        col=c("red", "gray"),
        ylab="SB velocity (mm/sec)",
        main="Random Noise",
        xlim=c(0,40),
        ylim=c(-40,30)
)
abline(v=10, col="gray", lty=5)
dev.off()


#SB acceleration
# plot(decompose(ts_CTL_rheo_agg1_SBacc1))
# 
deco_ts_CTL_rheo_agg1_SBacc1 <- decompose(ts_CTL_rheo_agg1_SBacc1)
obs_CTL_acc <- deco_ts_CTL_rheo_agg1_SBacc1$x
trd_CTL_acc  <- deco_ts_CTL_rheo_agg1_SBacc1$trend
seas_CTL_acc  <- deco_ts_CTL_rheo_agg1_SBacc1$seasonal
ran_CTL_acc  <- deco_ts_CTL_rheo_agg1_SBacc1$random

# deco_ts_Cu_rheo_agg1_SBacc2 <- decompose(ts_Cu_rheo_agg1_SBacc2)
# obs_Cu_acc <- deco_ts_Cu_rheo_agg1_SBacc2$x
# trd_Cu_acc  <- deco_ts_Cu_rheo_agg1_SBacc2$trend
# seas_Cu_acc  <- deco_ts_Cu_rheo_agg1_SBacc2$seasonal
# ran_Cu_acc  <- deco_ts_Cu_rheo_agg1_SBacc2$random
# 
# deco_ts_Neo_rheo_agg1_SBacc1 <- decompose(ts_Neo_rheo_agg1_SBacc1)
# obs_Neo_acc <- deco_ts_Neo_rheo_agg1_SBacc1$x
# trd_Neo_acc  <- deco_ts_Neo_rheo_agg1_SBacc1$trend
# seas_Neo_acc  <- deco_ts_Neo_rheo_agg1_SBacc1$seasonal
# ran_Neo_acc  <- deco_ts_Neo_rheo_agg1_SBacc1$random

# deco_ts_bap_rheo_agg1_SBacc1 <- decompose(ts_bap_rheo_agg1_SBacc1)
# obs_bap_acc <- deco_ts_bap_rheo_agg1_SBacc1$x
# trd_bap_acc  <- deco_ts_bap_rheo_agg1_SBacc1$trend
# seas_bap_acc  <- deco_ts_bap_rheo_agg1_SBacc1$seasonal
# ran_bap_acc  <- deco_ts_bap_rheo_agg1_SBacc1$random
# 
deco_ts_shk_rheo_agg1_SBacc1 <- decompose(ts_shk_rheo_agg1_SBacc1)
obs_shk_acc <- deco_ts_shk_rheo_agg1_SBacc1$x
trd_shk_acc  <- deco_ts_shk_rheo_agg1_SBacc1$trend
seas_shk_acc  <- deco_ts_shk_rheo_agg1_SBacc1$seasonal
ran_shk_acc  <- deco_ts_shk_rheo_agg1_SBacc1$random

# deco_ts_CTL8_rheo_agg2_SBacc1 <- decompose(ts_CTL8_rheo_agg2_SBacc1)
# obs_CTL8_acc <- deco_ts_CTL8_rheo_agg2_SBacc1$x
# trd_CTL8_acc  <- deco_ts_CTL8_rheo_agg2_SBacc1$trend
# seas_CTL8_acc  <- deco_ts_CTL8_rheo_agg2_SBacc1$seasonal
# ran_CTL8_acc  <- deco_ts_CTL8_rheo_agg2_SBacc1$random
# 
# deco_ts_shk8_rheo_agg2_SBacc1 <- decompose(ts_shk8_rheo_agg2_SBacc1)
# obs_shk8_acc <- deco_ts_shk8_rheo_agg2_SBacc1$x
# trd_shk8_acc  <- deco_ts_shk8_rheo_agg2_SBacc1$trend
# seas_shk8_acc  <- deco_ts_shk8_rheo_agg2_SBacc1$seasonal
# ran_shk8_acc  <- deco_ts_shk8_rheo_agg2_SBacc1$random


# all_obs_acc <- cbind(obs_Cu_acc, obs_Neo_acc, obs_CTL_acc) #this order helps with graph layer ordering bottom -> top :Cu-blue, Neo-gr, Ctl-gray
# all_trd_acc <- cbind(trd_Cu_acc, trd_Neo_acc, trd_CTL_acc)
# all_seas_acc <- cbind(seas_Cu_acc, seas_Neo_acc, seas_CTL_acc)
# all_ran_acc <- cbind(ran_Cu_acc, ran_Neo_acc, ran_CTL_acc)

# all_obs_acc <- cbind(obs_bap_acc, obs_CTL_acc) #this order helps with graph layer ordering bottom -> top :bap-blue, Neo-gr, Ctl-gray
# all_trd_acc <- cbind(trd_bap_acc, trd_CTL_acc)
# all_seas_acc <- cbind(seas_bap_acc, seas_CTL_acc)
# all_ran_acc <- cbind(ran_bap_acc, ran_CTL_acc)
# 
all_obs_acc <- cbind(obs_shk_acc, obs_CTL_acc) #this order helps with graph layer ordering bottom -> top :shk-blue, Neo-gr, Ctl-gray
all_trd_acc <- cbind(trd_shk_acc, trd_CTL_acc)
all_seas_acc <- cbind(seas_shk_acc, seas_CTL_acc)
all_ran_acc <- cbind(ran_shk_acc, ran_CTL_acc)

# all_obs_acc <- cbind(obs_shk48_acc, obs_CTL48_acc) #this order helps with graph layer ordering bottom -> top :shk-blue, Neo-gr, Ctl-gray
# all_trd_acc <- cbind(trd_shk48_acc, trd_CTL48_acc)
# all_seas_acc <- cbind(seas_shk8_acc, seas_CTL8_acc)
# all_ran_acc <- cbind(ran_shk8_acc, ran_CTL8_acc)


pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/SpecDeco_SBacc_obs_shake.pdf", width = 5, height = 4)
plot.ts(all_obs_acc,
        plot.type = "single", 
        # col=c("blue", "green", "gray"),
        col=c("red", "gray"),
        ylab="SB acceleration (mm/sec^2)",
        main="Observed",
        xlim=c(0,40),
        ylim=c(-500,400)
)
abline(v=10, col="gray", lty=5)
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/SpecDeco_SBacc_trd_shake.pdf", width = 5, height = 4)
plot.ts(all_trd_acc,
        plot.type = "single", 
        # col=c("blue", "green", "gray"),
        col=c("red", "gray"),
        ylab="SB acceleration (mm/sec^2)",
        main="Overall Trend",
        xlim=c(0,40),
        ylim=c(-80,80)
)
abline(v=10, col="gray", lty=5)
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/SpecDeco_SBacc_seas_shake.pdf", width = 5, height = 4)
plot.ts(all_seas_acc,
        plot.type = "single", 
        # col=c("blue", "green", "gray"),
        col=c("red", "gray"),
        ylab="SB acceleration (mm/sec^2)",
        main="Seasonal Fluctuation",
        xlim=c(0,5),
        ylim=c(-20,20)
)
abline(v=10, col="gray", lty=5)
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/SpecDeco_SBacc_ran_shake.pdf", width = 5, height = 4)
plot.ts(all_ran_acc,
        plot.type = "single", 
        # col=c("blue", "green", "gray"),
        col=c("red", "gray"),
        ylab="SB acceleration (mm/sec^2)",
        main="Random Noise",
        xlim=c(0,40),
        ylim=c(-550,300)
)
abline(v=10, col="gray", lty=5)
dev.off()



#mean body angle
# plot(decompose(ts_CTL_rheo_agg1_angle_rad1))

deco_ts_CTL_rheo_agg1_angle_rad1 <- decompose(ts_CTL_rheo_agg1_angle_rad1)
obs_CTL_angle_rad <- deco_ts_CTL_rheo_agg1_angle_rad1$x
trd_CTL_angle_rad  <- deco_ts_CTL_rheo_agg1_angle_rad1$trend
seas_CTL_angle_rad  <- deco_ts_CTL_rheo_agg1_angle_rad1$seasonal
ran_CTL_angle_rad  <- deco_ts_CTL_rheo_agg1_angle_rad1$random

# deco_ts_Cu_rheo_agg1_angle_rad2 <- decompose(ts_Cu_rheo_agg1_angle_rad2)
# obs_Cu_angle_rad <- deco_ts_Cu_rheo_agg1_angle_rad2$x
# trd_Cu_angle_rad  <- deco_ts_Cu_rheo_agg1_angle_rad2$trend
# seas_Cu_angle_rad  <- deco_ts_Cu_rheo_agg1_angle_rad2$seasonal
# ran_Cu_angle_rad  <- deco_ts_Cu_rheo_agg1_angle_rad2$random
# 
# deco_ts_Neo_rheo_agg1_angle_rad1 <- decompose(ts_Neo_rheo_agg1_angle_rad1)
# obs_Neo_angle_rad <- deco_ts_Neo_rheo_agg1_angle_rad1$x
# trd_Neo_angle_rad  <- deco_ts_Neo_rheo_agg1_angle_rad1$trend
# seas_Neo_angle_rad  <- deco_ts_Neo_rheo_agg1_angle_rad1$seasonal
# ran_Neo_angle_rad  <- deco_ts_Neo_rheo_agg1_angle_rad1$random
# 
# deco_ts_bap_rheo_agg1_angle_rad1 <- decompose(ts_bap_rheo_agg1_angle_rad1)
# obs_bap_angle_rad <- deco_ts_bap_rheo_agg1_angle_rad1$x
# trd_bap_angle_rad  <- deco_ts_bap_rheo_agg1_angle_rad1$trend
# seas_bap_angle_rad  <- deco_ts_bap_rheo_agg1_angle_rad1$seasonal
# ran_bap_angle_rad  <- deco_ts_bap_rheo_agg1_angle_rad1$random

deco_ts_shk_rheo_agg1_angle_rad1 <- decompose(ts_shk_rheo_agg1_angle_rad1)
obs_shk_angle_rad <- deco_ts_shk_rheo_agg1_angle_rad1$x
trd_shk_angle_rad  <- deco_ts_shk_rheo_agg1_angle_rad1$trend
seas_shk_angle_rad  <- deco_ts_shk_rheo_agg1_angle_rad1$seasonal
ran_shk_angle_rad  <- deco_ts_shk_rheo_agg1_angle_rad1$random

# deco_ts_CTL8_rheo_agg2_angle_rad1 <- decompose(ts_CTL8_rheo_agg2_angle_rad1)
# obs_CTL8_angle_rad <- deco_ts_CTL8_rheo_agg2_angle_rad1$x
# trd_CTL8_angle_rad  <- deco_ts_CTL8_rheo_agg2_angle_rad1$trend
# seas_CTL8_angle_rad  <- deco_ts_CTL8_rheo_agg2_angle_rad1$seasonal
# ran_CTL8_angle_rad  <- deco_ts_CTL8_rheo_agg2_angle_rad1$random
# 
# deco_ts_shk8_rheo_agg2_angle_rad1 <- decompose(ts_shk8_rheo_agg2_angle_rad1)
# obs_shk8_angle_rad <- deco_ts_shk8_rheo_agg2_angle_rad1$x
# trd_shk8_angle_rad  <- deco_ts_shk8_rheo_agg2_angle_rad1$trend
# seas_shk8_angle_rad  <- deco_ts_shk8_rheo_agg2_angle_rad1$seasonal
# ran_shk8_angle_rad  <- deco_ts_shk8_rheo_agg2_angle_rad1$random


# all_obs_angle_rad <- cbind(obs_Cu_angle_rad, obs_Neo_angle_rad, obs_CTL_angle_rad) #this order helps with graph layer ordering bottom -> top :Cu-blue, Neo-gr, Ctl-gray
# all_trd_angle_rad <- cbind(trd_Cu_angle_rad, trd_Neo_angle_rad, trd_CTL_angle_rad)
# all_seas_angle_rad <- cbind(seas_Cu_angle_rad, seas_Neo_angle_rad, seas_CTL_angle_rad)
# all_ran_angle_rad <- cbind(ran_Cu_angle_rad, ran_Neo_angle_rad, ran_CTL_angle_rad)

# all_obs_angle_rad <- cbind(obs_bap_angle_rad, obs_CTL_angle_rad) #this order helps with graph layer ordering bottom -> top :bap-blue, Neo-gr, Ctl-gray
# all_trd_angle_rad <- cbind(trd_bap_angle_rad, trd_CTL_angle_rad)
# all_seas_angle_rad <- cbind(seas_bap_angle_rad, seas_CTL_angle_rad)
# all_ran_angle_rad <- cbind(ran_bap_angle_rad, ran_CTL_angle_rad)

all_obs_angle_rad <- cbind(obs_shk_angle_rad, obs_CTL_angle_rad) #this order helps with graph layer ordering bottom -> top :shk-blue, Neo-gr, Ctl-gray
all_trd_angle_rad <- cbind(trd_shk_angle_rad, trd_CTL_angle_rad)
all_seas_angle_rad <- cbind(seas_shk_angle_rad, seas_CTL_angle_rad)
all_ran_angle_rad <- cbind(ran_shk_angle_rad, ran_CTL_angle_rad)

# all_obs_angle_rad <- cbind(obs_shk48_angle_rad, obs_CTL48_angle_rad) #this order helps with graph layer ordering bottom -> top :shk-blue, Neo-gr, Ctl-gray
# all_trd_angle_rad <- cbind(trd_shk0_angle_rad, trd_CTL0_angle_rad)
# all_seas_angle_rad <- cbind(seas_shk8_angle_rad, seas_CTL8_angle_rad)
# all_ran_angle_rad <- cbind(ran_shk8_angle_rad, ran_CTL8_angle_rad)

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/SpecDeco_SBang_obs_shake.pdf", width = 5, height = 4)
plot.ts(all_obs_angle_rad,
        plot.type = "single", 
        # col=c("blue", "green", "gray"),
        col=c("red", "gray"),
        ylab="Mean Body Angle (radians)",
        main="Observed",
        xlim=c(0,40),
        ylim=c(-0.7,0.7)
)
abline(v=10, col="gray", lty=5)
abline(h=0, col="gray", lty=3)
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/SpecDeco_SBang_trd_shake.pdf", width = 5, height = 4)
plot.ts(all_trd_angle_rad,
        plot.type = "single", 
        # col=c("blue", "green", "gray"),
        col=c("red", "gray"),
        ylab="Mean Body Angle (radians)",
        main="Overall Trend",
        xlim=c(0,40),
        ylim=c(-0.45,0.45)
)
abline(v=10, col="gray", lty=5)
abline(h=0, col="gray", lty=3)
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/SpecDeco_SBang_seas_shake.pdf", width = 5, height = 4)
plot.ts(all_seas_angle_rad,
        plot.type = "single", 
        # col=c("blue", "green", "gray"),
        col=c("red", "gray"),
        ylab="Mean Body Angle (radians)",
        main="Seasonal Fluctuation",
        xlim=c(0,5),
        ylim=c(-0.03,0.03)
)
abline(v=10, col="gray", lty=5)
abline(h=0, col="gray", lty=3)
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/SpecDeco_SBang_ran_shake.pdf", width = 5, height = 4)
plot.ts(all_ran_angle_rad,
        plot.type = "single", 
        # col=c("blue", "green", "gray"),
        col=c("red", "gray"),
        ylab="Mean Body Angle (radians)",
        main="Random Noise",
        xlim=c(0,40),
        ylim=c(-0.6,0.6)
)
abline(v=10, col="gray", lty=5)
abline(h=0, col="gray", lty=3)
dev.off()


#mean anglular dispersion (angular variance)
# plot(decompose(ts_CTL_rheo_agg1_angle_dispersion1))
# 
deco_ts_CTL_rheo_agg1_angle_dispersion1 <- decompose(ts_CTL_rheo_agg1_angle_dispersion1)
obs_CTL_angle_dispersion <- deco_ts_CTL_rheo_agg1_angle_dispersion1$x
trd_CTL_angle_dispersion  <- deco_ts_CTL_rheo_agg1_angle_dispersion1$trend
seas_CTL_angle_dispersion  <- deco_ts_CTL_rheo_agg1_angle_dispersion1$seasonal
ran_CTL_angle_dispersion  <- deco_ts_CTL_rheo_agg1_angle_dispersion1$random

# deco_ts_Cu_rheo_agg1_angle_dispersion2 <- decompose(ts_Cu_rheo_agg1_angle_dispersion2)
# obs_Cu_angle_dispersion <- deco_ts_Cu_rheo_agg1_angle_dispersion2$x
# trd_Cu_angle_dispersion  <- deco_ts_Cu_rheo_agg1_angle_dispersion2$trend
# seas_Cu_angle_dispersion  <- deco_ts_Cu_rheo_agg1_angle_dispersion2$seasonal
# ran_Cu_angle_dispersion  <- deco_ts_Cu_rheo_agg1_angle_dispersion2$random
# 
# deco_ts_Neo_rheo_agg1_angle_dispersion1 <- decompose(ts_Neo_rheo_agg1_angle_dispersion1)
# obs_Neo_angle_dispersion <- deco_ts_Neo_rheo_agg1_angle_dispersion1$x
# trd_Neo_angle_dispersion  <- deco_ts_Neo_rheo_agg1_angle_dispersion1$trend
# seas_Neo_angle_dispersion  <- deco_ts_Neo_rheo_agg1_angle_dispersion1$seasonal
# ran_Neo_angle_dispersion  <- deco_ts_Neo_rheo_agg1_angle_dispersion1$random

# deco_ts_bap_rheo_agg1_angle_dispersion1 <- decompose(ts_bap_rheo_agg1_angle_dispersion1)
# obs_bap_angle_dispersion <- deco_ts_bap_rheo_agg1_angle_dispersion1$x
# trd_bap_angle_dispersion  <- deco_ts_bap_rheo_agg1_angle_dispersion1$trend
# seas_bap_angle_dispersion  <- deco_ts_bap_rheo_agg1_angle_dispersion1$seasonal
# ran_bap_angle_dispersion  <- deco_ts_bap_rheo_agg1_angle_dispersion1$random

deco_ts_shk_rheo_agg1_angle_dispersion1 <- decompose(ts_shk_rheo_agg1_angle_dispersion1)
obs_shk_angle_dispersion <- deco_ts_shk_rheo_agg1_angle_dispersion1$x
trd_shk_angle_dispersion  <- deco_ts_shk_rheo_agg1_angle_dispersion1$trend
seas_shk_angle_dispersion  <- deco_ts_shk_rheo_agg1_angle_dispersion1$seasonal
ran_shk_angle_dispersion  <- deco_ts_shk_rheo_agg1_angle_dispersion1$random

# deco_ts_CTL8_rheo_agg2_angle_dispersion1 <- decompose(ts_CTL8_rheo_agg2_angle_dispersion1)
# obs_CTL8_angle_dispersion <- deco_ts_CTL8_rheo_agg2_angle_dispersion1$x
# trd_CTL8_angle_dispersion  <- deco_ts_CTL8_rheo_agg2_angle_dispersion1$trend
# seas_CTL8_angle_dispersion  <- deco_ts_CTL8_rheo_agg2_angle_dispersion1$seasonal
# ran_CTL8_angle_dispersion  <- deco_ts_CTL8_rheo_agg2_angle_dispersion1$random
# 
# deco_ts_shk8_rheo_agg2_angle_dispersion1 <- decompose(ts_shk8_rheo_agg2_angle_dispersion1)
# obs_shk8_angle_dispersion <- deco_ts_shk8_rheo_agg2_angle_dispersion1$x
# trd_shk8_angle_dispersion  <- deco_ts_shk8_rheo_agg2_angle_dispersion1$trend
# seas_shk8_angle_dispersion  <- deco_ts_shk8_rheo_agg2_angle_dispersion1$seasonal
# ran_shk8_angle_dispersion  <- deco_ts_shk8_rheo_agg2_angle_dispersion1$random


# all_obs_angle_dispersion <- cbind(obs_Cu_angle_dispersion, obs_Neo_angle_dispersion, obs_CTL_angle_dispersion) #this order helps with graph layer ordering bottom -> top :Cu-blue, Neo-gr, Ctl-gray
# all_trd_angle_dispersion <- cbind(trd_Cu_angle_dispersion, trd_Neo_angle_dispersion, trd_CTL_angle_dispersion)
# all_seas_angle_dispersion <- cbind(seas_Cu_angle_dispersion, seas_Neo_angle_dispersion, seas_CTL_angle_dispersion)
# all_ran_angle_dispersion <- cbind(ran_Cu_angle_dispersion, ran_Neo_angle_dispersion, ran_CTL_angle_dispersion)

# all_obs_angle_dispersion <- cbind(obs_bap_angle_dispersion, obs_CTL_angle_dispersion) #this order helps with graph layer ordering bottom -> top :bap-blue, Neo-gr, Ctl-gray
# all_trd_angle_dispersion <- cbind(trd_bap_angle_dispersion, trd_CTL_angle_dispersion)
# all_seas_angle_dispersion <- cbind(seas_bap_angle_dispersion, seas_CTL_angle_dispersion)
# all_ran_angle_dispersion <- cbind(ran_bap_angle_dispersion, ran_CTL_angle_dispersion)

all_obs_angle_dispersion <- cbind(obs_shk_angle_dispersion, obs_CTL_angle_dispersion) #this order helps with graph layer ordering bottom -> top :shk-blue, Neo-gr, Ctl-gray
all_trd_angle_dispersion <- cbind(trd_shk_angle_dispersion, trd_CTL_angle_dispersion)
all_seas_angle_dispersion <- cbind(seas_shk_angle_dispersion, seas_CTL_angle_dispersion)
all_ran_angle_dispersion <- cbind(ran_shk_angle_dispersion, ran_CTL_angle_dispersion)

# all_obs_angle_dispersion <- cbind(obs_shk8_angle_dispersion, obs_CTL8_angle_dispersion) #this order helps with graph layer ordering bottom -> top :shk-blue, Neo-gr, Ctl-gray
# all_trd_angle_dispersion <- cbind(trd_shk0_angle_dispersion, trd_CTL0_angle_dispersion)
# all_seas_angle_dispersion <- cbind(seas_shk48_angle_dispersion, seas_CTL48_angle_dispersion)
# all_ran_angle_dispersion <- cbind(ran_shk8_angle_dispersion, ran_CTL8_angle_dispersion)


pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/SpecDeco_SBresL_obs_shake.pdf", width = 5, height = 4)
plot.ts(all_obs_angle_dispersion,
        plot.type = "single", 
        # col=c("blue", "green", "gray"),
        col=c("red", "gray"),
        ylab="Mean Angular Dispersion (radians)",
        main="Observed",
        xlim=c(0,40),
        ylim=c(0.5,1)
)
abline(v=10, col="gray", lty=5)
abline(h=1, col="gray", lty=3)
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/SpecDeco_SBresL_trd_shake.pdf", width = 5, height = 4)
plot.ts(all_trd_angle_dispersion,
        plot.type = "single", 
        # col=c("blue", "green", "gray"),
        col=c("red", "gray"),
        ylab="Mean Angular Dispersion (radians)",
        main="Overall Trend",
        xlim=c(0,40),
        ylim=c(0.85,1)
)
abline(v=10, col="gray", lty=5)
abline(h=1, col="gray", lty=3)
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/SpecDeco_SBresL_seas_shake.pdf", width = 5, height = 4)
plot.ts(all_seas_angle_dispersion,
        plot.type = "single", 
        # col=c("blue", "green", "gray"),
        col=c("red", "gray"),
        ylab="Mean Angular Dispersion (radians)",
        main="Seasonal Fluctuation",
        xlim=c(0,5),
        ylim=c(-0.02,0.01)
)
abline(v=10, col="gray", lty=5)
abline(h=0, col="gray", lty=3)
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/SpecDeco_SBresL_ran_shake.pdf", width = 5, height = 4)
plot.ts(all_ran_angle_dispersion,
        plot.type = "single", 
        # col=c("blue", "green", "gray"),
        col=c("red", "gray"),
        ylab="Mean Angular Dispersion (radians)",
        main="Random Noise",
        xlim=c(0,40),
        ylim=c(-0.4,0.1)
)
abline(v=10, col="gray", lty=5)
abline(h=0, col="gray", lty=3)
dev.off()




####---FOURIER ANALYSIS - TEMPORAL COMPONENTS OF PERIODICITY ---####
#Controls
spec_ts_ctl_mov <- spectrum(ts_CTL_rheo_agg1_SBmov1, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
delta <- 1/.1 
spec_ctl_mov_X <- spec_ts_ctl_mov$freq/delta
spec_ctl_mov_Y <- 2*spec_ts_ctl_mov$spec #multiply spectrum by 2 to make equal to variance
plot.mov <- plot(spec_ctl_mov_X, spec_ctl_mov_Y, xlab="Period (sec)", ylab="Spectral Density", main= "SB_movement_CTL", type="l")

# spec_ts_Cu_mov <- spectrum(ts_Cu_rheo_agg1_SBmov2, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
# delta <- 1/.1 
# spec_Cu_mov_X <- spec_ts_Cu_mov$freq/delta
# spec_Cu_mov_Y <- 2*spec_ts_Cu_mov$spec #multiply spectrum by 2 to make equal to variance
# plot.mov <- plot(spec_Cu_mov_X, spec_Cu_mov_Y, xlab="Period (sec)", ylab="Spectral Density", main= "SB_movement_CuSO4", type="l")
# 
# spec_ts_Neo_mov <- spectrum(ts_Neo_rheo_agg1_SBmov1, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
# delta <- 1/.1 
# spec_Neo_mov_X <- spec_ts_Neo_mov$freq/delta
# spec_Neo_mov_Y <- 2*spec_ts_Neo_mov$spec #multiply spectrum by 2 to make equal to variance
# plot_mov <- plot(spec_Neo_mov_X, spec_Neo_mov_Y, xlab="Period (sec)", ylab="Spectral Density", main= "SB_movement_Neo", type="l")

spec_ts_ctl_vel <- spectrum(ts_CTL_rheo_agg1_SBvel1, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
delta <- 1/.1 
spec_ctl_vel_X <- spec_ts_ctl_vel$freq/delta
spec_ctl_vel_Y <- 2*spec_ts_ctl_vel$spec #multiply spectrum by 2 to make equal to variance
plot.vel <- plot(spec_ctl_vel_X, spec_ctl_vel_Y, xlab="Period (sec)", ylab="Spectral Density", main= "SB_velocity_CTL", type="l")

spec_ts_ctl_acc <- spectrum(ts_CTL_rheo_agg1_SBacc1, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
delta <- 1/.1 
spec_ctl_acc_X <- spec_ts_ctl_acc$freq/delta
spec_ctl_acc_Y <- 2*spec_ts_ctl_acc$spec #multiply spectrum by 2 to make equal to variance
plot.acc <- plot(spec_ctl_acc_X, spec_ctl_acc_Y, xlab="Period (sec)", ylab="Spectral Density", main= "SB_acceleration_CTL", type="l")

spec_ts_ctl_ang <- spectrum(ts_CTL_rheo_agg1_angle_rad1, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
delta <- 1/.1 
spec_ctl_ang_X <- spec_ts_ctl_ang$freq/delta
spec_ctl_ang_Y <- 2*spec_ts_ctl_ang$spec #multiply spectrum by 2 to make equal to variance
plot.ang <- plot(spec_ctl_ang_X, spec_ctl_ang_Y, xlab="Period (sec)", ylab="Spectral Density", main= "Body_angle_CTL", type="l")

spec_ts_ctl_dis <- spectrum(ts_CTL_rheo_agg1_angle_dispersion1, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
delta <- 1/.1 
spec_ctl_dis_X <- spec_ts_ctl_dis$freq/delta
spec_ctl_dis_Y <- 2*spec_ts_ctl_dis$spec #multiply spectrum by 2 to make equal to variance
plot.dis <- plot(spec_ctl_dis_X, spec_ctl_dis_Y, xlab="Period (sec)", ylab="Spectral Density", main= "Body_angle_dispersion_CTL", type="l")


#CuSO4
# spec_ts_Cu_mov <- spectrum(ts_Cu_rheo_agg1_SBmov2, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
# delta <- 1/.1 
# spec_Cu_mov_X <- spec_ts_Cu_mov$freq/delta
# spec_Cu_mov_Y <- 2*spec_ts_Cu_mov$spec #multiply spectrum by 2 to make equal to variance
# plot.mov <- plot(spec_Cu_mov_X, spec_Cu_mov_Y, xlab="Period (sec)", ylab="Spectral Density", main= "SB_movement_CuSO4", type="l")
# 
# spec_ts_Cu_vel <- spectrum(ts_Cu_rheo_agg1_SBvel2, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
# delta <- 1/.1 
# spec_Cu_vel_X <- spec_ts_Cu_vel$freq/delta
# spec_Cu_vel_Y <- 2*spec_ts_Cu_vel$spec #multiply spectrum by 2 to make equal to variance
# plot.vel <- plot(spec_Cu_vel_X, spec_Cu_vel_Y, xlab="Period (sec)", ylab="Spectral Density", main= "SB_velocity_CuSO4", type="l")
# 
# spec_ts_Cu_acc <- spectrum(ts_Cu_rheo_agg1_SBacc2, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
# delta <- 1/.1 
# spec_Cu_acc_X <- spec_ts_Cu_acc$freq/delta
# spec_Cu_acc_Y <- 2*spec_ts_Cu_acc$spec #multiply spectrum by 2 to make equal to variance
# plot.acc <- plot(spec_Cu_acc_X, spec_Cu_acc_Y, xlab="Period (sec)", ylab="Spectral Density", main= "SB_acceleration_CuSO4", type="l")
# 
# spec_ts_Cu_ang <- spectrum(ts_Cu_rheo_agg1_angle_rad1, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
# delta <- 1/.1 
# spec_Cu_ang_X <- spec_ts_Cu_ang$freq/delta
# spec_Cu_ang_Y <- 2*spec_ts_Cu_ang$spec #multiply spectrum by 2 to make equal to variance
# plot.ang <- plot(spec_Cu_ang_X, spec_Cu_ang_Y, xlab="Period (sec)", ylab="Spectral Density", main= "Body_angle_CuSO4", type="l")
# 
# spec_ts_Cu_dis <- spectrum(ts_Cu_rheo_agg1_angle_dispersion1, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
# delta <- 1/.1 
# spec_Cu_dis_X <- spec_ts_Cu_dis$freq/delta
# spec_Cu_dis_Y <- 2*spec_ts_Cu_dis$spec #multiply spectrum by 2 to make equal to variance
# plot.dis <- plot(spec_Cu_dis_X, spec_Cu_dis_Y, xlab="Period (sec)", ylab="Spectral Density", main= "Body_angle_dispersion_CuSO4", type="l")
# 
# 
# #Neomycin
# spec_ts_Neo_mov <- spectrum(ts_Neo_rheo_agg1_SBmov1, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
# delta <- 1/.1 
# spec_Neo_mov_X <- spec_ts_Neo_mov$freq/delta
# spec_Neo_mov_Y <- 2*spec_ts_Neo_mov$spec #multiply spectrum by 2 to make equal to variance
# plot_mov <- plot(spec_Neo_mov_X, spec_Neo_mov_Y, xlab="Frequency (Hz)", ylab="Spectral Density", main= "SB_movement_Neo", type="l")
# 
# spec_ts_Neo_vel <- spectrum(ts_Neo_rheo_agg1_SBvel1, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
# delta <- 1/.1 
# spec_Neo_vel_X <- spec_ts_Neo_vel$freq/delta
# spec_Neo_vel_Y <- 2*spec_ts_Neo_vel$spec #multiply spectrum by 2 to make equal to variance
# plot_vel <- plot(spec_Neo_vel_X, spec_Neo_vel_Y, xlab="Period (sec)", ylab="Spectral Density", main= "SB_velocity_Neo", type="l")
# 
# spec_ts_Neo_acc <- spectrum(ts_Neo_rheo_agg1_SBacc1, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
# delta <- 1/.1 
# spec_Neo_acc_X <- spec_ts_Neo_acc$freq/delta
# spec_Neo_acc_Y <- 2*spec_ts_Neo_acc$spec #multiply spectrum by 2 to make equal to variance
# plot_acc <- plot(spec_Neo_acc_X, spec_Neo_acc_Y, xlab="Period (sec)", ylab="Spectral Density", main= "SB_acceleration_Neo", type="l")
# 
# spec_ts_Neo_ang <- spectrum(ts_Neo_rheo_agg1_angle_rad1, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
# delta <- 1/.1 
# spec_Neo_ang_X <- spec_ts_Neo_ang$freq/delta
# spec_Neo_ang_Y <- 2*spec_ts_Neo_ang$spec #multiply spectrum by 2 to make equal to variance
# plot_ang <- plot(spec_Neo_ang_X, spec_Neo_ang_Y, xlab="Period (sec)", ylab="Spectral Density", main= "Body_angle_Neo", type="l")
# 
# spec_ts_Neo_dis <- spectrum(ts_Neo_rheo_agg1_angle_dispersion1, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
# delta <- 1/.1 
# spec_Neo_dis_X <- spec_ts_Neo_dis$freq/delta
# spec_Neo_dis_Y <- 2*spec_ts_Neo_dis$spec #multiply spectrum by 2 to make equal to variance
# plot_dis <- plot(spec_Neo_dis_X, spec_Neo_dis_Y, xlab="Period (sec)", ylab="Spectral Density", main= "Body_angle_dispersion_Neo", type="l")
# 
# 
# #BAPTA
# spec_ts_bap_mov <- spectrum(ts_bap_rheo_agg1_SBmov1, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
# delta <- 1/.1 
# spec_bap_mov_X <- spec_ts_bap_mov$freq/delta
# spec_bap_mov_Y <- 2*spec_ts_bap_mov$spec #multiply spectrum by 2 to make equal to variance
# plot.mov <- plot(spec_bap_mov_X, spec_bap_mov_Y, xlab="Period (sec)", ylab="Spectral Density", main= "SB_movement_BAPTA", type="l")
# 
# ts_bap_rheo_agg1_SBvel1a <- na.omit(ts_bap_rheo_agg1_SBvel1) #remove NA values from time series
# spec_ts_bap_vel <- spectrum(ts_bap_rheo_agg1_SBvel1a, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
# delta <- 1/.1 
# spec_bap_vel_X <- spec_ts_bap_vel$freq/delta
# spec_bap_vel_Y <- 2*spec_ts_bap_vel$spec #multiply spectrum by 2 to make equal to variance
# plot.vel <- plot(spec_bap_vel_X, spec_bap_vel_Y, xlab="Period (sec)", ylab="Spectral Density", main= "SB_velocity_BAPTA", type="l")
# 
# ts_bap_rheo_agg1_SBacc1a <- na.omit(ts_bap_rheo_agg1_SBacc1) #remove NA values from time series
# spec_ts_bap_acc <- spectrum(ts_bap_rheo_agg1_SBacc1a, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
# delta <- 1/.1 
# spec_bap_acc_X <- spec_ts_bap_acc$freq/delta
# spec_bap_acc_Y <- 2*spec_ts_bap_acc$spec #multiply spectrum by 2 to make equal to variance
# plot.acc <- plot(spec_bap_acc_X, spec_bap_acc_Y, xlab="Period (sec)", ylab="Spectral Density", main= "SB_acceleration_BAPTA", type="l")
# 
# spec_ts_bap_ang <- spectrum(ts_bap_rheo_agg1_angle_rad1, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
# delta <- 1/.1 
# spec_bap_ang_X <- spec_ts_bap_ang$freq/delta
# spec_bap_ang_Y <- 2*spec_ts_bap_ang$spec #multiply spectrum by 2 to make equal to variance
# plot.ang <- plot(spec_bap_ang_X, spec_bap_ang_Y, xlab="Period (sec)", ylab="Spectral Density", main= "Body_angle_BAPTA", type="l")
# 
# spec_ts_bap_dis <- spectrum(ts_bap_rheo_agg1_angle_dispersion1, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
# delta <- 1/.1 
# spec_bap_dis_X <- spec_ts_bap_dis$freq/delta
# spec_bap_dis_Y <- 2*spec_ts_bap_dis$spec #multiply spectrum by 2 to make equal to variance
# plot.dis <- plot(spec_bap_dis_X, spec_bap_dis_Y, xlab="Period (sec)", ylab="Spectral Density", main= "Body_angle_dispersion_BAPTA", type="l")



#SHAKE
spec_ts_shk_mov <- spectrum(ts_shk_rheo_agg1_SBmov1, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
delta <- 1/.1 
spec_shk_mov_X <- spec_ts_shk_mov$freq/delta
spec_shk_mov_Y <- 2*spec_ts_shk_mov$spec #multiply spectrum by 2 to make equal to variance
plot.mov <- plot(spec_shk_mov_X, spec_shk_mov_Y, xlab="Period (sec)", ylab="Spectral Density", main= "SB_movement_shake", type="l")

ts_shk_rheo_agg1_SBvel1a <- na.omit(ts_shk_rheo_agg1_SBvel1) #remove NA values from time series
spec_ts_shk_vel <- spectrum(ts_shk_rheo_agg1_SBvel1a, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
delta <- 1/.1 
spec_shk_vel_X <- spec_ts_shk_vel$freq/delta
spec_shk_vel_Y <- 2*spec_ts_shk_vel$spec #multiply spectrum by 2 to make equal to variance
plot.vel <- plot(spec_shk_vel_X, spec_shk_vel_Y, xlab="Period (sec)", ylab="Spectral Density", main= "SB_velocity_shake", type="l")

ts_shk_rheo_agg1_SBacc1a <- na.omit(ts_shk_rheo_agg1_SBacc1) #remove NA values from time series
spec_ts_shk_acc <- spectrum(ts_shk_rheo_agg1_SBacc1a, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
delta <- 1/.1 
spec_shk_acc_X <- spec_ts_shk_acc$freq/delta
spec_shk_acc_Y <- 2*spec_ts_shk_acc$spec #multiply spectrum by 2 to make equal to variance
plot.acc <- plot(spec_shk_acc_X, spec_shk_acc_Y, xlab="Period (sec)", ylab="Spectral Density", main= "SB_acceleration_shake", type="l")

spec_ts_shk_ang <- spectrum(ts_shk_rheo_agg1_angle_rad1, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
delta <- 1/.1 
spec_shk_ang_X <- spec_ts_shk_ang$freq/delta
spec_shk_ang_Y <- 2*spec_ts_shk_ang$spec #multiply spectrum by 2 to make equal to variance
plot.ang <- plot(spec_shk_ang_X, spec_shk_ang_Y, xlab="Period (sec)", ylab="Spectral Density", main= "Body_angle_shake", type="l")

spec_ts_shk_dis <- spectrum(ts_shk_rheo_agg1_angle_dispersion1, log="no", spans=c(7, 7), plot=FALSE) #convert time series to frequency by dividing my sampling interval (100ms, 1s, 10s, etc)
delta <- 1/.1 
spec_shk_dis_X <- spec_ts_shk_dis$freq/delta
spec_shk_dis_Y <- 2*spec_ts_shk_dis$spec #multiply spectrum by 2 to make equal to variance
plot.dis <- plot(spec_shk_dis_X, spec_shk_dis_Y, xlab="Period (sec)", ylab="Spectral Density", main= "Body_angle_dispersion_shake", type="l")



#FINAL FIGURES - overlays --FREQUENCY COMPONENTS (TEMPORAL CHANGES) OF SEASONALITY/PERIODICITY DATA!
#SBmov
pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/FreqComp_SBmov_shake.pdf", width = 5, height = 4)
# plot(x=spec_Neo_mov_X, y=spec_Neo_mov_Y, col="green", xlab="Frequency (Hz)", ylab="Spectral Density", xaxs="i", yaxs="i", xlim=c(0, 0.5), ylim=c(0, 0.05), main= "SB_movement", type="l")
plot(x=spec_shk_mov_X, y=spec_shk_mov_Y, col="red", xlab="Frequency (Hz)", ylab="Spectral Density", xaxs="i", yaxs="i", xlim=c(0, 0.5), ylim=c(0, 0.11), main= "SB_movement", type="l")
# lines(x=spec_Cu_mov_X, y=spec_Cu_mov_Y, col="blue")
lines(x=spec_ctl_mov_X, y=spec_ctl_mov_Y, col="gray")
dev.off()

#SBvel
pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/FreqComp_SBvel_shake.pdf", width = 5, height = 4)
# plot(x=spec_Neo_vel_X, y=spec_Neo_vel_Y, col="green", xlab="Frequency (Hz)", ylab="Spectral Density", xaxs="i", yaxs="i", xlim=c(0, 0.5), ylim=c(0, 8), main= "SB_velocity", type="l")
plot(x=spec_shk_vel_X, y=spec_shk_vel_Y, col="red", xlab="Frequency (Hz)", ylab="Spectral Density", xaxs="i", yaxs="i", xlim=c(0, 0.5), ylim=c(0, 8), main= "SB_velocity", type="l")
# lines(x=spec_Cu_vel_X, y=spec_Cu_vel_Y, col="blue")
lines(x=spec_ctl_vel_X, y=spec_ctl_vel_Y, col="gray")
dev.off()


#SBacc
pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/FreqComp_SBacc_shake.pdf", width = 5, height = 4)
# plot(x=spec_Neo_acc_X, y=spec_Neo_acc_Y, col="green", xlab="Frequency (Hz)", ylab="Spectral Density", xaxs="i", yaxs="i", xlim=c(0, 0.5), ylim=c(0, 1000), main= "SB_acceleration", type="l")
plot(x=spec_shk_acc_X, y=spec_shk_acc_Y, col="red", xlab="Frequency (Hz)", ylab="Spectral Density", xaxs="i", yaxs="i", xlim=c(0, 0.5), ylim=c(0, 1800), main= "SB_acceleration", type="l")
# lines(x=spec_Cu_acc_X, y=spec_Cu_acc_Y, col="blue")
lines(x=spec_ctl_acc_X, y=spec_ctl_acc_Y, col="gray")
dev.off()


#Mean body angle
pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/FreqComp_SBang_shake.pdf", width = 5, height = 4)
# plot(x=spec_Neo_ang_X, y=spec_Neo_ang_Y, col="green", xlab="Frequency (Hz)", ylab="Spectral Density", xaxs="i", yaxs="i", xlim=c(0, 0.5), ylim=c(0, 0.006), main= "Mean Body Angle", type="l")
plot(x=spec_shk_ang_X, y=spec_shk_ang_Y, col="red", xlab="Frequency (Hz)", ylab="Spectral Density", xaxs="i", yaxs="i", xlim=c(0, 0.5), ylim=c(0, 0.045), main= "Mean Body Angle", type="l")
# lines(x=spec_Cu_ang_X, y=spec_Cu_ang_Y, col="blue")
lines(x=spec_ctl_ang_X, y=spec_ctl_ang_Y, col="gray")
dev.off()


#body angle dispersion (variance)
pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/FreqComp_SBresL_shake.pdf", width = 5, height = 4)
# plot(x=spec_Neo_dis_X, y=spec_Neo_dis_Y, col="green", xlab="Frequency (Hz)", ylab="Spectral Density", xaxs="i", yaxs="i", xlim=c(0, 0.5), ylim=c(0, 0.0010), main= "Angular Dispersion", type="l")
plot(x=spec_shk_dis_X, y=spec_shk_dis_Y, col="red", xlab="Frequency (Hz)", ylab="Spectral Density", xaxs="i", yaxs="i", xlim=c(0, 0.5), ylim=c(0, 0.0015), main= "Mean Resultant Length", type="l")
# lines(x=spec_Cu_dis_X, y=spec_Cu_dis_Y, col="blue")
lines(x=spec_ctl_dis_X, y=spec_ctl_dis_Y, col="gray")
dev.off()





####---FIGURES - CROSS-CORRELATION OVERLAYS ---####
#note that arrangement is changed to overly lines and make treatments stands out from shakes
#CCF
pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/CCF_Mov-Ang_shake.pdf", width = 5, height = 4)
#plot(x=cu1$lag, y=cu1$acf, xlab="Lag", ylab="ACF", xlim=c(-2,2), ylim=c(-0.5, 0.5), main= "CCF: movement~angle", type="b", col = rgb(0, 0, 1, 1))
plot(x=shk1$lag, y=shk1$acf, xlab="Lag", ylab="ACF", xlim=c(-2,2), ylim=c(-0.35, 0.55), main= "CCF: movement~angle", type="b", col = rgb(1, 0, 0, 1))
#lines(x=neo1$lag, y=neo1$acf, type="b", col = rgb(0, 1, 0, 1))
lines(x=ctl1$lag, y=ctl1$acf, type="b", col = rgb(0.5, 0.5, 0.5, 1))
abline(h=0, col = rgb(0, 0, 0, 0.25), lty=1)
abline(h=c(-0.12, 0.12), col = rgb(1, 0, 0, .5), lty=3)
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/CCF_Mov-ResL_shake.pdf", width = 5, height = 4)
# plot(x=shk2$lag, y=shk2$acf, xlab="Lag", ylab="ACF", xlim=c(-2,2), ylim=c(0, 0.5), main= "CCF: movement~dispersion", type="b", col = rgb(0, 0, 1, 1))
plot(x=shk2$lag, y=shk2$acf, xlab="Lag", ylab="ACF", xlim=c(-2,2), ylim=c(-0.2, 0.4), main= "CCF: movement~dispersion", type="b", col = rgb(1, 0, 0, 1))
#lines(x=neo2$lag, y=neo2$acf, type="b", col = rgb(0, 1, 0, 1))
lines(x=ctl2$lag, y=ctl2$acf, type="b", col = rgb(0.5, 0.5, 0.5, 1))
abline(h=0, col = rgb(0, 0, 0, 0.25), lty=1)
abline(h=c(-0.12, 0.12), col = rgb(1, 0, 0, .5), lty=3)
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/CCF_Vel-Ang_shake.pdf", width = 5, height = 4)
# plot(x=shk3$lag, y=shk3$acf, xlab="Lag", ylab="ACF", xlim=c(-2,2), ylim=c(-0.4, 0.4), main= "CCF: velocity~angle", type="b", col = rgb(0, 0, 1, 1))
plot(x=shk3$lag, y=shk3$acf, xlab="Lag", ylab="ACF", xlim=c(-2,2), ylim=c(-0.2, 0.2), main= "CCF: velocity~angle", type="b", col = rgb(1, 0, 0, 1))
#lines(x=neo3$lag, y=neo3$acf, type="b", col = rgb(0, 1, 0, 1))
lines(x=ctl3$lag, y=ctl3$acf, type="b", col = rgb(0.5, 0.5, 0.5, 1))
abline(h=0, col = rgb(0, 0, 0, 0.25), lty=1)
abline(h=c(-0.12, 0.12), col = rgb(1, 0, 0, .5), lty=3)
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/CCF_Vel-ResL_shake.pdf", width = 5, height = 4)
# plot(x=shk8$lag, y=shk8$acf, xlab="Lag", ylab="ACF", xlim=c(-2,2), ylim=c(-0.2, 0.3), main= "CCF: velocity~dispersion", type="b", col = rgb(0, 0, 1, 1))
plot(x=shk4$lag, y=shk4$acf, xlab="Lag", ylab="ACF", xlim=c(-2,2), ylim=c(-0.2, 0.25), main= "CCF: velocity~dispersion", type="b", col = rgb(1, 0, 0, 1))
#lines(x=neo4$lag, y=neo4$acf, type="b", col = rgb(0, 1, 0, 1))
lines(x=ctl4$lag, y=ctl4$acf, type="b", col = rgb(0.5, 0.5, 0.5, 1))
abline(h=0, col = rgb(0, 0, 0, 0.25), lty=1)
abline(h=c(-0.12, 0.12), col = rgb(1, 0, 0, .5), lty=3)
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/CCF_Acc-Ang_shake.pdf", width = 5, height = 4)
# plot(x=shk5$lag, y=shk5$acf, xlab="Lag", ylab="ACF", xlim=c(-2,2), ylim=c(-0.3, 0.4), main= "CCF: acceleration~angle", type="b", col = rgb(0, 0, 1, 1))
plot(x=shk5$lag, y=shk5$acf, xlab="Lag", ylab="ACF", xlim=c(-2,2), ylim=c(-0.2, 0.2), main= "CCF: acceleration~angle", type="b", col = rgb(1, 0, 0, 1))
#lines(x=neo5$lag, y=neo5$acf, type="b", col = rgb(0, 1, 0, 1))
lines(x=ctl5$lag, y=ctl5$acf, type="b", col = rgb(0.5, 0.5, 0.5, 1))
abline(h=0, col = rgb(0, 0, 0, 0.25), lty=1)
abline(h=c(-0.12, 0.12), col = rgb(1, 0, 0, .5), lty=3)
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/CCF_Acc-ResL_shake.pdf", width = 5, height = 4)
# plot(x=shk6$lag, y=shk6$acf, xlab="Lag", ylab="ACF", xlim=c(-2,2), ylim=c(-0.2, 0.2), main= "CCF:acceleration~dispersion", type="b", col = rgb(0, 0, 1, 1))
plot(x=shk6$lag, y=shk6$acf, xlab="Lag", ylab="ACF", xlim=c(-2,2), ylim=c(-0.25, 0.3), main= "CCF:acceleration~dispersion", type="b", col = rgb(1, 0, 0, 1))
#lines(x=neo6$lag, y=neo6$acf, type="b", col = rgb(0, 1, 0, 1))
lines(x=ctl6$lag, y=ctl6$acf, type="b", col = rgb(0.5, 0.5, 0.5, 1))
abline(h=0, col = rgb(0, 0, 0, 0.25), lty=1)
abline(h=c(-0.12, 0.12), col = rgb(1, 0, 0, .5), lty=3)
dev.off()










####---Just for funzies! TOTAL SPATIAL USE (ROI) DURING TIME SERIES ---####
library(plyr)
counts5 <- ddply(master.all, .(master.all$Time, master.all$Treatment, master.all$Space_bin), nrow) #ALL data!! - RHEOTAXIS AND NOT RHEOTAXIS

names (counts5)[names(counts5) == "master.all$Treatment"] <- "Treatment"
names (counts5)[names(counts5) == "master.all$Time"] <- "Time"
names (counts5)[names(counts5) == "master.all$Space_bin"] <- "Space"
names (counts5)[names(counts5) == "V1"] <- "Counts"
counts5a <- counts5[!is.na(counts5$Space), ] #remove NA from Space_bin

# OR try dplyr
#counts5z <- master.all %>% (Time, Treatment, Space_bin)

#first must change counts to PROPORTION in each bin/all bins
time.vec2 <- unique(counts5a$Time) #extract unique elements and return a vector without duplicates
treat.vec2 <- unique(counts5a$Treatment)
counts5b <- data.frame()

for(i in 1:length(treat.vec2)){
  data.temp <- counts5a[counts5a$Treatment==treat.vec2[i], ]
  for(j in 1:length(time.vec2)){
    data.temp2 <- data.temp[data.temp$Time==time.vec2[j], ]
    sum.counts <- sum(data.temp2$Counts)
    data.temp2$Prop <- data.temp2$Counts/sum.counts
    data.temp2$Sum_time <- data.temp2$Prop*sum.counts
    counts5b <- rbind(counts5b, data.temp2)
    
  }
  
}

#remove post-stimulus data if desired (based on 100 ms time bins)
counts5b <- filter(counts5b, counts5b$Time<=30000) 


#lots of setup
loess_Ctl_f <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Front", ]) 
loess_Ctl_b <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Back", ])
loess_Ctl_c <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Center", ]) 
loess_Ctl_l <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Left", ])
loess_Ctl_r <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Right", ]) 

# loess_Cu_f <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Front", ]) 
# loess_Cu_b <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Back", ])
# loess_Cu_c <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Center", ]) 
# loess_Cu_l <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Left", ])
# loess_Cu_r <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Right", ]) 
# loess_Neo_f <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Front", ]) 
# loess_Neo_b <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Back", ])
# loess_Neo_c <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Center", ]) 
# loess_Neo_l <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Left", ])
# loess_Neo_r <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Right", ]) 
# 
# loess_bap_f <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="Bapta" & counts5b$Space=="Front", ]) 
# loess_bap_b <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="Bapta" & counts5b$Space=="Back", ])
# loess_bap_c <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="Bapta" & counts5b$Space=="Center", ]) 
# loess_bap_l <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="Bapta" & counts5b$Space=="Left", ])
# loess_bap_r <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="Bapta" & counts5b$Space=="Right", ]) 

loess_shk_f <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Front", ]) 
loess_shk_b <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Back", ])
loess_shk_c <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Center", ]) 
loess_shk_l <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Left", ])
loess_shk_r <- loess(Prop ~ Time, data=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Right", ]) 


smoothed_Ctl_f <- predict(loess_Ctl_f)
smoothed_Ctl_b <- predict(loess_Ctl_b)
smoothed_Ctl_c <- predict(loess_Ctl_c)
smoothed_Ctl_l <- predict(loess_Ctl_l)
smoothed_Ctl_r <- predict(loess_Ctl_r)

# smoothed_Cu_f <- predict(loess_Cu_f)
# smoothed_Cu_b <- predict(loess_Cu_b)
# smoothed_Cu_c <- predict(loess_Cu_c)
# smoothed_Cu_l <- predict(loess_Cu_l)
# smoothed_Cu_r <- predict(loess_Cu_r)
# smoothed_Neo_f <- predict(loess_Neo_f)
# smoothed_Neo_b <- predict(loess_Neo_b)
# smoothed_Neo_c <- predict(loess_Neo_c)
# smoothed_Neo_l <- predict(loess_Neo_l)
# smoothed_Neo_r <- predict(loess_Neo_r)

# smoothed_bap_f <- predict(loess_bap_f)
# smoothed_bap_b <- predict(loess_bap_b)
# smoothed_bap_c <- predict(loess_bap_c)
# smoothed_bap_l <- predict(loess_bap_l)
# smoothed_bap_r <- predict(loess_bap_r)

smoothed_shk_f <- predict(loess_shk_f)
smoothed_shk_b <- predict(loess_shk_b)
smoothed_shk_c <- predict(loess_shk_c)
smoothed_shk_l <- predict(loess_shk_l)
smoothed_shk_r <- predict(loess_shk_r)


pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/Shake_graphs/ROI_prop_rheo+non_Shake.pdf", width = 5, height = 4)
plot(1, type="n", xlab="Time", ylab="Proportion of Time in Space Bin", xlim=c(0, 45000), ylim=c(0, .7), main = "Space_bin_TS_alldata")
lines(smoothed_Ctl_f, x=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Front", ]$Time, col="green", lty=1)
lines(smoothed_Ctl_b, x=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Back", ]$Time, col="black", lty=1)
lines(smoothed_Ctl_c, x=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Center", ]$Time, col="red", lty=1)
lines(smoothed_Ctl_l, x=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Left", ]$Time, col="blue", lty=1)
lines(smoothed_Ctl_r, x=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Right", ]$Time, col="light blue", lty=1)
# 
# lines(smoothed_Cu_f, x=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Front", ]$Time, col="green", lty=2)
# lines(smoothed_Cu_b, x=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Back", ]$Time, col="black", lty=2)
# lines(smoothed_Cu_c, x=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Center", ]$Time, col="red", lty=2)
# lines(smoothed_Cu_l, x=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Left", ]$Time, col="blue", lty=2)
# lines(smoothed_Cu_r, x=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Right", ]$Time, col="light blue", lty=2)
# lines(smoothed_Neo_f, x=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Front", ]$Time, col="green", lty=3)
# lines(smoothed_Neo_b, x=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Back", ]$Time, col="black", lty=3)
# lines(smoothed_Neo_c, x=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Center", ]$Time, col="red", lty=3)
# lines(smoothed_Neo_l, x=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Left", ]$Time, col="blue", lty=3)
# lines(smoothed_Neo_r, x=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Right", ]$Time, col="light blue", lty=3)

# lines(smoothed_shk_f, x=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Front", ]$Time, col="green", lty=2)
# lines(smoothed_shk_b, x=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Back", ]$Time, col="black", lty=2)
# lines(smoothed_shk_c, x=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Center", ]$Time, col="red", lty=2)
# lines(smoothed_shk_l, x=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Left", ]$Time, col="blue", lty=2)
# lines(smoothed_shk_r, x=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Right", ]$Time, col="light blue", lty=2)

lines(smoothed_shk_f, x=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Front", ]$Time, col="green", lty=2)
lines(smoothed_shk_b, x=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Back", ]$Time, col="black", lty=2)
lines(smoothed_shk_c, x=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Center", ]$Time, col="red", lty=2)
lines(smoothed_shk_l, x=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Left", ]$Time, col="blue", lty=2)
lines(smoothed_shk_r, x=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Right", ]$Time, col="light blue", lty=2)

# legend(31000, 0.45, legend=c("Center", "Front","Left", "Right", "Back"),
#        col=c("red", "green", "blue", "light blue", "black"), lty=1, cex=0.8)
# legend(31000, 0.2, legend=c("Control", "CuSO4","Neo"),
#        lty=c(1,2,3), cex=0.8)
# legend(41000, 0.6, legend=c("Center", "Front","Left", "Right", "Back"),
#        col=c("red", "green", "blue", "light blue", "black"), lty=1, cex=0.8)
legend(41000, 0.35, legend=c("Control", "Shake"),
       lty=c(1,2,3), cex=0.8)
abline(v=10000, col="gray", lty=4)
dev.off()

# #same dataframe as proportion but looking at counts in a time bin
# #lots of setup
# loess_Ctl_f <- loess(Counts ~ Time, data=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Front", ]) 
# loess_Ctl_b <- loess(Counts ~ Time, data=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Back", ])
# loess_Ctl_c <- loess(Counts ~ Time, data=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Center", ]) 
# loess_Ctl_l <- loess(Counts ~ Time, data=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Left", ])
# loess_Ctl_r <- loess(Counts ~ Time, data=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Right", ]) 
# 
# loess_Cu_f <- loess(Counts ~ Time, data=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Front", ]) 
# loess_Cu_b <- loess(Counts ~ Time, data=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Back", ])
# loess_Cu_c <- loess(Counts ~ Time, data=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Center", ]) 
# loess_Cu_l <- loess(Counts ~ Time, data=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Left", ])
# loess_Cu_r <- loess(Counts ~ Time, data=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Right", ]) 
# loess_Neo_f <- loess(Counts ~ Time, data=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Front", ]) 
# loess_Neo_b <- loess(Counts ~ Time, data=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Back", ])
# loess_Neo_c <- loess(Counts ~ Time, data=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Center", ]) 
# loess_Neo_l <- loess(Counts ~ Time, data=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Left", ])
# loess_Neo_r <- loess(Counts ~ Time, data=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Right", ]) 
# 
# loess_bap_f <- loess(Counts ~ Time, data=counts5b[counts5b$Treatment=="Bapta" & counts5b$Space=="Front", ]) 
# loess_bap_b <- loess(Counts ~ Time, data=counts5b[counts5b$Treatment=="Bapta" & counts5b$Space=="Back", ])
# loess_bap_c <- loess(Counts ~ Time, data=counts5b[counts5b$Treatment=="Bapta" & counts5b$Space=="Center", ]) 
# loess_bap_l <- loess(Counts ~ Time, data=counts5b[counts5b$Treatment=="Bapta" & counts5b$Space=="Left", ])
# loess_bap_r <- loess(Counts ~ Time, data=counts5b[counts5b$Treatment=="Bapta" & counts5b$Space=="Right", ]) 
# 
# 
# smoothed_Ctl_f <- predict(loess_Ctl_f)
# smoothed_Ctl_b <- predict(loess_Ctl_b)
# smoothed_Ctl_c <- predict(loess_Ctl_c)
# smoothed_Ctl_l <- predict(loess_Ctl_l)
# smoothed_Ctl_r <- predict(loess_Ctl_r)
# 
# smoothed_Cu_f <- predict(loess_Cu_f)
# smoothed_Cu_b <- predict(loess_Cu_b)
# smoothed_Cu_c <- predict(loess_Cu_c)
# smoothed_Cu_l <- predict(loess_Cu_l)
# smoothed_Cu_r <- predict(loess_Cu_r)
# smoothed_Neo_f <- predict(loess_Neo_f)
# smoothed_Neo_b <- predict(loess_Neo_b)
# smoothed_Neo_c <- predict(loess_Neo_c)
# smoothed_Neo_l <- predict(loess_Neo_l)
# smoothed_Neo_r <- predict(loess_Neo_r)
# 
# smoothed_bap_f <- predict(loess_bap_f)
# smoothed_bap_b <- predict(loess_bap_b)
# smoothed_bap_c <- predict(loess_bap_c)
# smoothed_bap_l <- predict(loess_bap_l)
# smoothed_bap_r <- predict(loess_bap_r)
# 
# #final figure of TS in ROI
# plot(1, type="n", xlab="Time", ylab="Number of Occurances in Space Bin", xlim=c(0, 35000), ylim=c(0, 40), main = "Space_bin_TS_alldata")
# lines(smoothed_Ctl_f, x=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Front", ]$Time, col="green", lty=1)
# lines(smoothed_Ctl_b, x=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Back", ]$Time, col="black", lty=1)
# lines(smoothed_Ctl_c, x=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Center", ]$Time, col="red", lty=1)
# lines(smoothed_Ctl_l, x=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Left", ]$Time, col="blue", lty=1)
# lines(smoothed_Ctl_r, x=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Right", ]$Time, col="light blue", lty=1)
# 
# lines(smoothed_Cu_f, x=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Front", ]$Time, col="green", lty=2)
# lines(smoothed_Cu_b, x=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Back", ]$Time, col="black", lty=2)
# lines(smoothed_Cu_c, x=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Center", ]$Time, col="red", lty=2)
# lines(smoothed_Cu_l, x=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Left", ]$Time, col="blue", lty=2)
# lines(smoothed_Cu_r, x=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Right", ]$Time, col="light blue", lty=2)
# lines(smoothed_Neo_f, x=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Front", ]$Time, col="green", lty=3)
# lines(smoothed_Neo_b, x=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Back", ]$Time, col="black", lty=3)
# lines(smoothed_Neo_c, x=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Center", ]$Time, col="red", lty=3)
# lines(smoothed_Neo_l, x=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Left", ]$Time, col="blue", lty=3)
# lines(smoothed_Neo_r, x=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Right", ]$Time, col="light blue", lty=3)
# 
# lines(smoothed_bap_f, x=counts5b[counts5b$Treatment=="Bapta" & counts5b$Space=="Front", ]$Time, col="green", lty=3)
# lines(smoothed_bap_b, x=counts5b[counts5b$Treatment=="Bapta" & counts5b$Space=="Back", ]$Time, col="black", lty=3)
# lines(smoothed_bap_c, x=counts5b[counts5b$Treatment=="Bapta" & counts5b$Space=="Center", ]$Time, col="red", lty=3)
# lines(smoothed_bap_l, x=counts5b[counts5b$Treatment=="Bapta" & counts5b$Space=="Left", ]$Time, col="blue", lty=3)
# lines(smoothed_bap_r, x=counts5b[counts5b$Treatment=="Bapta" & counts5b$Space=="Right", ]$Time, col="light blue", lty=3)
# 
# 
# # legend(31000, 0.45, legend=c("Center", "Front","Left", "Right", "Back"),
# #        col=c("red", "green", "blue", "light blue", "black"), lty=1, cex=0.8)
# # legend(31000, 0.2, legend=c("Control", "CuSO4","Neo"),
# #        lty=c(1,2,3), cex=0.8)
# legend(31000, 30, legend=c("Center", "Front","Left", "Right", "Back"),
#        col=c("red", "green", "blue", "light blue", "black"), lty=1, cex=0.8)
# legend(31000, 20, legend=c("Control", "Bapta"),
#        lty=c(1,2,3), cex=0.8)
# abline(v=10000, col="gray", lty=4)



#same dataframe as proportion but looking at total time of rheotaxis in a time bin -

loess_Ctl_f <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Front", ]) 
loess_Ctl_b <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Back", ])
loess_Ctl_c <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Center", ]) 
loess_Ctl_l <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Left", ])
loess_Ctl_r <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Right", ]) 

# loess_Cu_f <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Front", ]) 
# loess_Cu_b <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Back", ])
# loess_Cu_c <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Center", ]) 
# loess_Cu_l <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Left", ])
# loess_Cu_r <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Right", ]) 
# loess_Neo_f <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Front", ]) 
# loess_Neo_b <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Back", ])
# loess_Neo_c <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Center", ]) 
# loess_Neo_l <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Left", ])
# loess_Neo_r <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Right", ]) 

# loess_bap_f <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="Bapta" & counts5b$Space=="Front", ]) 
# loess_bap_b <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="Bapta" & counts5b$Space=="Back", ])
# loess_bap_c <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="Bapta" & counts5b$Space=="Center", ]) 
# loess_bap_l <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="Bapta" & counts5b$Space=="Left", ])
# loess_bap_r <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="Bapta" & counts5b$Space=="Right", ]) 

loess_shk_f <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Front", ]) 
loess_shk_b <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Back", ])
loess_shk_c <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Center", ]) 
loess_shk_l <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Left", ])
loess_shk_r <- loess(Sum_time ~ Time, data=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Right", ]) 


smoothed_Ctl_f <- predict(loess_Ctl_f)
smoothed_Ctl_b <- predict(loess_Ctl_b)
smoothed_Ctl_c <- predict(loess_Ctl_c)
smoothed_Ctl_l <- predict(loess_Ctl_l)
smoothed_Ctl_r <- predict(loess_Ctl_r)

# smoothed_Cu_f <- predict(loess_Cu_f)
# smoothed_Cu_b <- predict(loess_Cu_b)
# smoothed_Cu_c <- predict(loess_Cu_c)
# smoothed_Cu_l <- predict(loess_Cu_l)
# smoothed_Cu_r <- predict(loess_Cu_r)
# smoothed_Neo_f <- predict(loess_Neo_f)
# smoothed_Neo_b <- predict(loess_Neo_b)
# smoothed_Neo_c <- predict(loess_Neo_c)
# smoothed_Neo_l <- predict(loess_Neo_l)
# smoothed_Neo_r <- predict(loess_Neo_r)

# smoothed_bap_f <- predict(loess_bap_f)
# smoothed_bap_b <- predict(loess_bap_b)
# smoothed_bap_c <- predict(loess_bap_c)
# smoothed_bap_l <- predict(loess_bap_l)
# smoothed_bap_r <- predict(loess_bap_r)

smoothed_shk_f <- predict(loess_shk_f)
smoothed_shk_b <- predict(loess_shk_b)
smoothed_shk_c <- predict(loess_shk_c)
smoothed_shk_l <- predict(loess_shk_l)
smoothed_shk_r <- predict(loess_shk_r)

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/ROI_sum_time_rheo+non_shake.pdf", width = 5, height = 4)
plot(1, type="n", xlab="Time", ylab="Total Time in Space Bin", xlim=c(0, 45000), ylim=c(0, 35), main = "Space_bin_TS_alldata")
lines(smoothed_Ctl_f, x=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Front", ]$Time, col="green", lty=1)
lines(smoothed_Ctl_b, x=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Back", ]$Time, col="black", lty=1)
lines(smoothed_Ctl_c, x=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Center", ]$Time, col="red", lty=1)
lines(smoothed_Ctl_l, x=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Left", ]$Time, col="blue", lty=1)
lines(smoothed_Ctl_r, x=counts5b[counts5b$Treatment=="Control" & counts5b$Space=="Right", ]$Time, col="light blue", lty=1)

# lines(smoothed_Cu_f, x=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Front", ]$Time, col="green", lty=2)
# lines(smoothed_Cu_b, x=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Back", ]$Time, col="black", lty=2)
# lines(smoothed_Cu_c, x=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Center", ]$Time, col="red", lty=2)
# lines(smoothed_Cu_l, x=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Left", ]$Time, col="blue", lty=2)
# lines(smoothed_Cu_r, x=counts5b[counts5b$Treatment=="CuSO4" & counts5b$Space=="Right", ]$Time, col="light blue", lty=2)
# lines(smoothed_Neo_f, x=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Front", ]$Time, col="green", lty=3)
# lines(smoothed_Neo_b, x=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Back", ]$Time, col="black", lty=3)
# lines(smoothed_Neo_c, x=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Center", ]$Time, col="red", lty=3)
# lines(smoothed_Neo_l, x=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Left", ]$Time, col="blue", lty=3)
# lines(smoothed_Neo_r, x=counts5b[counts5b$Treatment=="Neo" & counts5b$Space=="Right", ]$Time, col="light blue", lty=3)

lines(smoothed_shk_f, x=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Front", ]$Time, col="green", lty=2)
lines(smoothed_shk_b, x=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Back", ]$Time, col="black", lty=2)
lines(smoothed_shk_c, x=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Center", ]$Time, col="red", lty=2)
lines(smoothed_shk_l, x=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Left", ]$Time, col="blue", lty=2)
lines(smoothed_shk_r, x=counts5b[counts5b$Treatment=="Shake" & counts5b$Space=="Right", ]$Time, col="light blue", lty=2)

# legend(31000, 120, legend=c("Center", "Front","Left", "Right", "Back"),
#        col=c("red", "green", "blue", "light blue", "black"), lty=1, cex=0.8)
# legend(31000, 80, legend=c("Control", "CuSO4","Neo"),
#        lty=c(1,2,3), cex=0.8)
# legend(41000, 40, legend=c("Center", "Front","Left", "Right", "Back"),
#        col=c("red", "green", "blue", "light blue", "black"), lty=1, cex=0.8)
legend(41000, 20, legend=c("Control", "Shake"),
       lty=c(1,2,3), cex=0.8)
abline(v=10000, col="gray", lty=4)
dev.off()






####--- SPATIAL USE (ROI) DURING RHEOTAXIS OVER TIME SERIES ---####

counts2 <- ddply(master.all, .(master.all$Time, master.all$Treatment, master.all$Rheotaxis, master.all$Space_bin), nrow) #RHEOTAXIS only
names (counts2)[names(counts2) == "master.all$Rheotaxis"] <- "Rheotaxis"
names (counts2)[names(counts2) == "master.all$Treatment"] <- "Treatment"
names (counts2)[names(counts2) == "master.all$Time"] <- "Time"
names (counts2)[names(counts2) == "master.all$Space_bin"] <- "Space"
names (counts2)[names(counts2) == "V1"] <- "Counts"
counts2$Rheotaxis2 <- counts2$Rheotaxis
counts2$Rheotaxis2[is.na(counts2$Rheotaxis2)] <- 0 #replace NA with 0 but keep original NA data

counts6a <- counts2[counts2$Rheotaxis2==1, ] #rheotaxis only
# counts6a <- counts6a[counts6a$Time<=30000, ] #trim post stimulus data based on 100 ms time bins
counts6a <- counts6a[!is.na(counts6a$Space), ] #remove NA from Space_bin


#repeat for RHEOTAXIS data only
#first must change counts to proportion in each bin/all bins
time.vec2 <- unique(counts6a$Time) #extract unique elements and return a vector without duplicates
treat.vec2 <- unique(counts6a$Treatment)
counts6b <- data.frame()

for(i in 1:length(treat.vec2)){
  data.temp <- counts6a[counts6a$Treatment==treat.vec2[i], ]
  for(j in 1:length(time.vec2)){
    data.temp2 <- data.temp[data.temp$Time==time.vec2[j], ]
    sum.counts <- sum(data.temp2$Counts)
    data.temp2$Prop <- data.temp2$Counts/sum.counts
    data.temp2$Sum_time <- data.temp2$Prop*sum.counts
    counts6b <- rbind(counts6b, data.temp2)
    
  }
  
}




#lots of setup
loess_Ctl_f <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Front", ]) 
loess_Ctl_b <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Back", ])
loess_Ctl_c <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Center", ]) 
loess_Ctl_l <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Left", ])
loess_Ctl_r <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Right", ]) 
# 
# loess_Cu_f <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Front", ]) 
# loess_Cu_b <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Back", ])
# loess_Cu_c <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Center", ]) 
# loess_Cu_l <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Left", ])
# loess_Cu_r <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Right", ]) 
# loess_Neo_f <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Front", ]) 
# loess_Neo_b <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Back", ])
# loess_Neo_c <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Center", ]) 
# loess_Neo_l <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Left", ])
# loess_Neo_r <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Right", ]) 

# loess_bap_f <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Front", ]) 
# loess_bap_b <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Back", ])
# loess_bap_c <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Center", ]) 
# loess_bap_l <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Left", ])
# loess_bap_r <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Right", ]) 

loess_shk_f <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="Shake" & counts6b$Space=="Front", ]) 
loess_shk_b <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="Shake" & counts6b$Space=="Back", ])
loess_shk_c <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="Shake" & counts6b$Space=="Center", ]) 
loess_shk_l <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="Shake" & counts6b$Space=="Left", ])
loess_shk_r <- loess(Prop ~ Time, data=counts6b[counts6b$Treatment=="Shake" & counts6b$Space=="Right", ]) 


smoothed_Ctl_f <- predict(loess_Ctl_f)
smoothed_Ctl_b <- predict(loess_Ctl_b)
smoothed_Ctl_c <- predict(loess_Ctl_c)
smoothed_Ctl_l <- predict(loess_Ctl_l)
smoothed_Ctl_r <- predict(loess_Ctl_r)

# smoothed_Cu_f <- predict(loess_Cu_f)
# smoothed_Cu_b <- predict(loess_Cu_b)
# smoothed_Cu_c <- predict(loess_Cu_c)
# smoothed_Cu_l <- predict(loess_Cu_l)
# smoothed_Cu_r <- predict(loess_Cu_r)
# smoothed_Neo_f <- predict(loess_Neo_f)
# smoothed_Neo_b <- predict(loess_Neo_b)
# smoothed_Neo_c <- predict(loess_Neo_c)
# smoothed_Neo_l <- predict(loess_Neo_l)
# smoothed_Neo_r <- predict(loess_Neo_r)

# smoothed_shk_f <- predict(loess_shk_f)
# smoothed_shk_b <- predict(loess_shk_b)
# smoothed_shk_c <- predict(loess_shk_c)
# smoothed_shk_l <- predict(loess_shk_l)
# smoothed_shk_r <- predict(loess_shk_r)

smoothed_shk_f <- predict(loess_shk_f)
smoothed_shk_b <- predict(loess_shk_b)
smoothed_shk_c <- predict(loess_shk_c)
smoothed_shk_l <- predict(loess_shk_l)
smoothed_shk_r <- predict(loess_shk_r)


pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/ROI_prop_rheo_shake.pdf", width = 5, height = 4)
plot(1, type="n", xlab="Time", ylab="Proportion of Time in Space Bin", xlim=c(0, 45000), ylim=c(0, 1), main = "Space_bin_TS_rheo_data")
lines(smoothed_Ctl_f, x=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Front", ]$Time, col="green", lty=1)
lines(smoothed_Ctl_b, x=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Back", ]$Time, col="black", lty=1)
lines(smoothed_Ctl_c, x=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Center", ]$Time, col="red", lty=1)
lines(smoothed_Ctl_l, x=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Left", ]$Time, col="blue", lty=1)
lines(smoothed_Ctl_r, x=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Right", ]$Time, col="light blue", lty=1)

# lines(smoothed_Cu_f, x=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Front", ]$Time, col="green", lty=2)
# lines(smoothed_Cu_b, x=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Back", ]$Time, col="black", lty=2)
# lines(smoothed_Cu_c, x=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Center", ]$Time, col="red", lty=2)
# lines(smoothed_Cu_l, x=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Left", ]$Time, col="blue", lty=2)
# lines(smoothed_Cu_r, x=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Right", ]$Time, col="light blue", lty=2)
# lines(smoothed_Neo_f, x=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Front", ]$Time, col="green", lty=3)
# lines(smoothed_Neo_b, x=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Back", ]$Time, col="black", lty=3)
# lines(smoothed_Neo_c, x=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Center", ]$Time, col="red", lty=3)
# lines(smoothed_Neo_l, x=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Left", ]$Time, col="blue", lty=3)
# lines(smoothed_Neo_r, x=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Right", ]$Time, col="light blue", lty=3)

# lines(smoothed_bap_f, x=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Front", ]$Time, col="green", lty=2)
# lines(smoothed_bap_b, x=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Back", ]$Time, col="black", lty=2)
# lines(smoothed_bap_c, x=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Center", ]$Time, col="red", lty=2)
# lines(smoothed_bap_l, x=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Left", ]$Time, col="blue", lty=2)
# lines(smoothed_bap_r, x=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Right", ]$Time, col="light blue", lty=2)

lines(smoothed_shk_f, x=counts6b[counts6b$Treatment=="Shake" & counts6b$Space=="Front", ]$Time, col="green", lty=2)
lines(smoothed_shk_b, x=counts6b[counts6b$Treatment=="Shake" & counts6b$Space=="Back", ]$Time, col="black", lty=2)
lines(smoothed_shk_c, x=counts6b[counts6b$Treatment=="Shake" & counts6b$Space=="Center", ]$Time, col="red", lty=2)
lines(smoothed_shk_l, x=counts6b[counts6b$Treatment=="Shake" & counts6b$Space=="Left", ]$Time, col="blue", lty=2)
lines(smoothed_shk_r, x=counts6b[counts6b$Treatment=="Shake" & counts6b$Space=="Right", ]$Time, col="light blue", lty=2)

# legend(31000, 0.45, legend=c("Center", "Front","Left", "Right", "Back"),
#        col=c("red", "green", "blue", "light blue", "black"), lty=1, cex=0.8)
# legend(31000, 0.2, legend=c("Control", "CuSO4","Neo"),
#        lty=c(1,2,3), cex=0.8)
# legend(41000, 0.95, legend=c("Center", "Front","Left", "Right", "Back"),
#        col=c("red", "green", "blue", "light blue", "black"), lty=1, cex=0.8)
legend(41000, 0.5, legend=c("Control", "Shake"),
       lty=c(1,2,3), cex=0.8)
abline(v=10000, col="gray", lty=4)
dev.off()

# #same dataframe as proportion but looking at counts in a time bin
# #lots of setup
# loess_Ctl_f <- loess(Counts ~ Time, data=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Front", ]) 
# loess_Ctl_b <- loess(Counts ~ Time, data=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Back", ])
# loess_Ctl_c <- loess(Counts ~ Time, data=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Center", ]) 
# loess_Ctl_l <- loess(Counts ~ Time, data=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Left", ])
# loess_Ctl_r <- loess(Counts ~ Time, data=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Right", ]) 
# 
# loess_Cu_f <- loess(Counts ~ Time, data=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Front", ]) 
# loess_Cu_b <- loess(Counts ~ Time, data=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Back", ])
# loess_Cu_c <- loess(Counts ~ Time, data=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Center", ]) 
# loess_Cu_l <- loess(Counts ~ Time, data=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Left", ])
# loess_Cu_r <- loess(Counts ~ Time, data=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Right", ]) 
# loess_Neo_f <- loess(Counts ~ Time, data=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Front", ]) 
# loess_Neo_b <- loess(Counts ~ Time, data=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Back", ])
# loess_Neo_c <- loess(Counts ~ Time, data=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Center", ]) 
# loess_Neo_l <- loess(Counts ~ Time, data=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Left", ])
# loess_Neo_r <- loess(Counts ~ Time, data=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Right", ]) 
# 
# loess_bap_f <- loess(Counts ~ Time, data=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Front", ]) 
# loess_bap_b <- loess(Counts ~ Time, data=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Back", ])
# loess_bap_c <- loess(Counts ~ Time, data=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Center", ]) 
# loess_bap_l <- loess(Counts ~ Time, data=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Left", ])
# loess_bap_r <- loess(Counts ~ Time, data=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Right", ]) 
# 
# 
# smoothed_Ctl_f <- predict(loess_Ctl_f)
# smoothed_Ctl_b <- predict(loess_Ctl_b)
# smoothed_Ctl_c <- predict(loess_Ctl_c)
# smoothed_Ctl_l <- predict(loess_Ctl_l)
# smoothed_Ctl_r <- predict(loess_Ctl_r)
# 
# smoothed_Cu_f <- predict(loess_Cu_f)
# smoothed_Cu_b <- predict(loess_Cu_b)
# smoothed_Cu_c <- predict(loess_Cu_c)
# smoothed_Cu_l <- predict(loess_Cu_l)
# smoothed_Cu_r <- predict(loess_Cu_r)
# smoothed_Neo_f <- predict(loess_Neo_f)
# smoothed_Neo_b <- predict(loess_Neo_b)
# smoothed_Neo_c <- predict(loess_Neo_c)
# smoothed_Neo_l <- predict(loess_Neo_l)
# smoothed_Neo_r <- predict(loess_Neo_r)
# 
# smoothed_bap_f <- predict(loess_bap_f)
# smoothed_bap_b <- predict(loess_bap_b)
# smoothed_bap_c <- predict(loess_bap_c)
# smoothed_bap_l <- predict(loess_bap_l)
# smoothed_bap_r <- predict(loess_bap_r)
# 
# #final figure of TS in ROI
# plot(1, type="n", xlab="Time", ylab="Number of Occurances in Space Bin", xlim=c(0, 35000), ylim=c(0, 20), main = "Space_bin_TS_rheo_data")
# lines(smoothed_Ctl_f, x=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Front", ]$Time, col="green", lty=1)
# lines(smoothed_Ctl_b, x=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Back", ]$Time, col="black", lty=1)
# lines(smoothed_Ctl_c, x=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Center", ]$Time, col="red", lty=1)
# lines(smoothed_Ctl_l, x=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Left", ]$Time, col="blue", lty=1)
# lines(smoothed_Ctl_r, x=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Right", ]$Time, col="light blue", lty=1)
# 
# lines(smoothed_Cu_f, x=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Front", ]$Time, col="green", lty=2)
# lines(smoothed_Cu_b, x=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Back", ]$Time, col="black", lty=2)
# lines(smoothed_Cu_c, x=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Center", ]$Time, col="red", lty=2)
# lines(smoothed_Cu_l, x=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Left", ]$Time, col="blue", lty=2)
# lines(smoothed_Cu_r, x=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Right", ]$Time, col="light blue", lty=2)
# lines(smoothed_Neo_f, x=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Front", ]$Time, col="green", lty=3)
# lines(smoothed_Neo_b, x=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Back", ]$Time, col="black", lty=3)
# lines(smoothed_Neo_c, x=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Center", ]$Time, col="red", lty=3)
# lines(smoothed_Neo_l, x=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Left", ]$Time, col="blue", lty=3)
# lines(smoothed_Neo_r, x=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Right", ]$Time, col="light blue", lty=3)
# 
# lines(smoothed_bap_f, x=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Front", ]$Time, col="green", lty=3)
# lines(smoothed_bap_b, x=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Back", ]$Time, col="black", lty=3)
# lines(smoothed_bap_c, x=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Center", ]$Time, col="red", lty=3)
# lines(smoothed_bap_l, x=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Left", ]$Time, col="blue", lty=3)
# lines(smoothed_bap_r, x=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Right", ]$Time, col="light blue", lty=3)
# 
# 
# # legend(31000, 0.45, legend=c("Center", "Front","Left", "Right", "Back"),
# #        col=c("red", "green", "blue", "light blue", "black"), lty=1, cex=0.8)
# # legend(31000, 0.2, legend=c("Control", "CuSO4","Neo"),
# #        lty=c(1,2,3), cex=0.8)
# legend(31000, 30, legend=c("Center", "Front","Left", "Right", "Back"),
#        col=c("red", "green", "blue", "light blue", "black"), lty=1, cex=0.8)
# legend(31000, 20, legend=c("Control", "Bapta"),
#        lty=c(1,2,3), cex=0.8)
# abline(v=10000, col="gray", lty=4)



#same dataframe as proportion but looking at total time of rheotaxis in a time bin -

loess_Ctl_f <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Front", ]) 
loess_Ctl_b <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Back", ])
loess_Ctl_c <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Center", ]) 
loess_Ctl_l <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Left", ])
loess_Ctl_r <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Right", ]) 

# loess_Cu_f <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Front", ]) 
# loess_Cu_b <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Back", ])
# loess_Cu_c <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Center", ]) 
# loess_Cu_l <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Left", ])
# loess_Cu_r <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Right", ]) 
# loess_Neo_f <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Front", ]) 
# loess_Neo_b <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Back", ])
# loess_Neo_c <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Center", ]) 
# loess_Neo_l <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Left", ])
# loess_Neo_r <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Right", ]) 

# loess_bap_f <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Front", ]) 
# loess_bap_b <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Back", ])
# loess_bap_c <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Center", ]) 
# loess_bap_l <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Left", ])
# loess_bap_r <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Right", ]) 

loess_shk_f <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="Shake" & counts6b$Space=="Front", ]) 
loess_shk_b <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="Shake" & counts6b$Space=="Back", ])
loess_shk_c <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="Shake" & counts6b$Space=="Center", ]) 
loess_shk_l <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="Shake" & counts6b$Space=="Left", ])
loess_shk_r <- loess(Sum_time ~ Time, data=counts6b[counts6b$Treatment=="Shake" & counts6b$Space=="Right", ]) 


smoothed_Ctl_f <- predict(loess_Ctl_f)
smoothed_Ctl_b <- predict(loess_Ctl_b)
smoothed_Ctl_c <- predict(loess_Ctl_c)
smoothed_Ctl_l <- predict(loess_Ctl_l)
smoothed_Ctl_r <- predict(loess_Ctl_r)
# 
# smoothed_Cu_f <- predict(loess_Cu_f)
# smoothed_Cu_b <- predict(loess_Cu_b)
# smoothed_Cu_c <- predict(loess_Cu_c)
# smoothed_Cu_l <- predict(loess_Cu_l)
# smoothed_Cu_r <- predict(loess_Cu_r)
# smoothed_Neo_f <- predict(loess_Neo_f)
# smoothed_Neo_b <- predict(loess_Neo_b)
# smoothed_Neo_c <- predict(loess_Neo_c)
# smoothed_Neo_l <- predict(loess_Neo_l)
# smoothed_Neo_r <- predict(loess_Neo_r)

# smoothed_bap_f <- predict(loess_bap_f)
# smoothed_bap_b <- predict(loess_bap_b)
# smoothed_bap_c <- predict(loess_bap_c)
# smoothed_bap_l <- predict(loess_bap_l)
# smoothed_bap_r <- predict(loess_bap_r)

smoothed_shk_f <- predict(loess_shk_f)
smoothed_shk_b <- predict(loess_shk_b)
smoothed_shk_c <- predict(loess_shk_c)
smoothed_shk_l <- predict(loess_shk_l)
smoothed_shk_r <- predict(loess_shk_r)

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/ROI_sum_time_rheo_shake.pdf", width = 5, height = 4)
plot(1, type="n", xlab="Time", ylab="Total Time in Space Bin", xlim=c(0, 45000), ylim=c(0, 1), main = "Space_bin_TS_rheo_data")
lines(smoothed_Ctl_f, x=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Front", ]$Time, col="green", lty=1)
lines(smoothed_Ctl_b, x=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Back", ]$Time, col="black", lty=1)
lines(smoothed_Ctl_c, x=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Center", ]$Time, col="red", lty=1)
lines(smoothed_Ctl_l, x=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Left", ]$Time, col="blue", lty=1)
lines(smoothed_Ctl_r, x=counts6b[counts6b$Treatment=="Control" & counts6b$Space=="Right", ]$Time, col="light blue", lty=1)

# lines(smoothed_Cu_f, x=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Front", ]$Time, col="green", lty=2)
# lines(smoothed_Cu_b, x=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Back", ]$Time, col="black", lty=2)
# lines(smoothed_Cu_c, x=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Center", ]$Time, col="red", lty=2)
# lines(smoothed_Cu_l, x=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Left", ]$Time, col="blue", lty=2)
# lines(smoothed_Cu_r, x=counts6b[counts6b$Treatment=="CuSO4" & counts6b$Space=="Right", ]$Time, col="light blue", lty=2)
# lines(smoothed_Neo_f, x=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Front", ]$Time, col="green", lty=3)
# lines(smoothed_Neo_b, x=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Back", ]$Time, col="black", lty=3)
# lines(smoothed_Neo_c, x=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Center", ]$Time, col="red", lty=3)
# lines(smoothed_Neo_l, x=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Left", ]$Time, col="blue", lty=3)
# lines(smoothed_Neo_r, x=counts6b[counts6b$Treatment=="Neo" & counts6b$Space=="Right", ]$Time, col="light blue", lty=3)

# lines(smoothed_bap_f, x=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Front", ]$Time, col="green", lty=2)
# lines(smoothed_bap_b, x=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Back", ]$Time, col="black", lty=2)
# lines(smoothed_bap_c, x=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Center", ]$Time, col="red", lty=2)
# lines(smoothed_bap_l, x=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Left", ]$Time, col="blue", lty=2)
# lines(smoothed_bap_r, x=counts6b[counts6b$Treatment=="Bapta" & counts6b$Space=="Right", ]$Time, col="light blue", lty=2)

lines(smoothed_shk_f, x=counts6b[counts6b$Treatment=="Shake" & counts6b$Space=="Front", ]$Time, col="green", lty=2)
lines(smoothed_shk_b, x=counts6b[counts6b$Treatment=="Shake" & counts6b$Space=="Back", ]$Time, col="black", lty=2)
lines(smoothed_shk_c, x=counts6b[counts6b$Treatment=="Shake" & counts6b$Space=="Center", ]$Time, col="red", lty=2)
lines(smoothed_shk_l, x=counts6b[counts6b$Treatment=="Shake" & counts6b$Space=="Left", ]$Time, col="blue", lty=2)
lines(smoothed_shk_r, x=counts6b[counts6b$Treatment=="Shake" & counts6b$Space=="Right", ]$Time, col="light blue", lty=2)

# legend(31000, 120, legend=c("Center", "Front","Left", "Right", "Back"),
#        col=c("red", "green", "blue", "light blue", "black"), lty=1, cex=0.8)
# legend(31000, 80, legend=c("Control", "CuSO4","Neo"),
#        lty=c(1,2,3), cex=0.8)
# legend(41000, 1, legend=c("Center", "Front","Left", "Right", "Back"),
#        col=c("red", "green", "blue", "light blue", "black"), lty=1, cex=0.8)
legend(41000, .5, legend=c("Control", "Shake"),
       lty=c(1,2,3), cex=0.8)
abline(v=10000, col="gray", lty=4)
dev.off()







  
  
