#install.packages
library(tidyverse)
library(circular)
library(bpnreg)
library(viridis)


##################################################
# START HERE FOR ANALYSES --- just use master.all.csv and load into master.all data frame 
##################################################
setwd("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results")

#all these files should have been run through rheotaxis_file_prep_concatenate.R to unify format
#master.all <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/60shake/master.data.shake.rec.csv") #SHAKE - RECOVERY 
# master.all <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/60bapta_test//master.data.bapta.noKS.csv") #BAPTA
master.all <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/sixty_flow_test/no_KS/master.data.60s.shake.csv") #60s flow test - Shake
# master.all <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/sixty_flow_test/no_KS/master.data.60s.CuNeo.csv") #60s flow test - CuSO4 v Neomycin
# master.all <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/CuNeo/master.all.CuNeo.csv") #CuSO4 v Neomycin 
# master.all <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/Old_shake/master.all.shake.csv") #Shake - no recovery! 


# master.all$File <- as.character(master.all$File)
master.all$Time <- as.integer(master.all$Time)
master.all$Individual <- as.factor(master.all$Individual)
master.all$Treatment <- as.factor(master.all$Treatment)
master.all$Stimulus <- as.factor(master.all$Stimulus)
master.all$Space_bin <- as.factor(master.all$Space_bin)
# master.all$Rec_Time <- as.factor(master.all$Rec_Time) #only for shake - recovery

View(master.all)
str(master.all)


# 11/23/21 must run this to divide Stimulus column of master.all.csv into two 10s bins =  Pre, 10s, 20s, Post bins are equal!  
master.all$Stimulus2 <- ifelse(master.all$Time<10001, "Pre_Stim",
                              ifelse(master.all$Time<20001, "Stim_10s",
                                     ifelse(master.all$Time<30001, "Stim_20s", "Post_stim")))
master.all$Stimulus2 <- as.factor(master.all$Stimulus2)

# 11/29/21 must run this to divide Stimulus column of master.all.csv into six 10s bins =  Pre, 10s, 20s, 30s, 40s, 50s, 60s Post bins are equal!  
# master.all$Stimulus3 <- ifelse(master.all$Time<10001, "Pre_Stim",
#                                ifelse(master.all$Time<20001, "Stim_10s",
#                                       ifelse(master.all$Time<30001, "Stim_20s", 
#                                              ifelse(master.all$Time<40001, "Stim_30s",
#                                                     ifelse(master.all$Time<50001, "Stim_40s",
#                                                            ifelse(master.all$Time<60001, "Stim_50s",
#                                                                   ifelse(master.all$Time<70001, "Stim_60s", "Post_Stim")))))))
# master.all$Stimulus3 <- as.factor(master.all$Stimulus3)

master.all <-  master.all[, c(2:21)] 

#counts3_agg2 <- aggregate(list(master.all[ , c(1,3,5,6,9:13,17,19)]), by = list(master.all$Rheotaxis, master.all$Treatment, master.all$Time), mean, na.rm=TRUE) 
# master.all.indiv_mean <- aggregate(list(master.all[, c(5,6,16:19)]), by = list(master.all$Individual,master.all$Treatment, master.all$Stimulus, master.all$Rec_Time), mean, na.rm=TRUE) #Shake - recovery
master.all.indiv.mean <- aggregate(list(master.all[, c(2:4,12:13)]), by = list(master.all$Individual,master.all$Treatment, master.all$Stimulus2), mean, na.rm=TRUE) #might have to adjust column numbers to select correctly
# master.all.indiv.mean <- aggregate(list(master.all[, c(3,4,21,13:14)]), by = list(master.all$Individual,master.all$Treatment, master.all$Stimulus3), mean, na.rm=TRUE) #on ly for 60s flow stimulus

# master.all_indiv_mean <- aggregate(list(master.all[, c(2:6)]), by = list(master.all$Individual,master.all$Treatment, master.all$Stimulus, master.all$Rec_Time), mean, na.rm=TRUE)
names (master.all.indiv.mean)[names(master.all.indiv.mean) == "Group.1"] <- "Individual"
names (master.all.indiv.mean)[names(master.all.indiv.mean) == "Group.2"] <- "Treatment"
names (master.all.indiv.mean)[names(master.all.indiv.mean) == "Group.3"] <- "Stimulus"
names (master.all.indiv.mean)[names(master.all.indiv.mean) == "Fish_angle_sin"] <- "Angle_sin"
names (master.all.indiv.mean)[names(master.all.indiv.mean) == "Fish_angle_cos"] <- "Angle_cos"
master.all.indiv.mean$Angle_rad2 <- atan2(master.all.indiv.mean$Angle_sin, master.all.indiv.mean$Angle_cos)

initial.angle <- master.all.indiv.mean$Angle_rad2 #convert +/- Radians (+/-180 degress) to 2pi or 6.28 radians 
final.angle <- initial.angle + 2*pi
master.all.indiv.mean$Angle_rad <- ifelse(initial.angle < 0, final.angle, initial.angle)
master.all.indiv.mean1 <- master.all.indiv.mean[ , c(1:3,7:10)] #delete extra columns

rheotaxis_data_all <- master.all.indiv.mean1


####---CuSO4-NEOMYCIN SUBSETS---####

#CONTROL: Stimulus dataframes
rheo_CTL <- filter(rheotaxis_data_all, rheotaxis_data_all$Treatment == "Control")
rheo_CTL <- na.omit(rheo_CTL) #remove NA 
rheo_CTL_ps <- filter(rheo_CTL, rheo_CTL$Stimulus == "Pre_Stim")
rheo_CTL_10s <- filter(rheo_CTL, rheo_CTL$Stimulus == "Stim_10s")
rheo_CTL_20s <- filter(rheo_CTL, rheo_CTL$Stimulus == "Stim_20s")
rheo_CTL_30s <- filter(rheo_CTL, rheo_CTL$Stimulus == "Stim_30s")
rheo_CTL_40s <- filter(rheo_CTL, rheo_CTL$Stimulus == "Stim_40s")
rheo_CTL_50s <- filter(rheo_CTL, rheo_CTL$Stimulus == "Stim_50s")
rheo_CTL_60s <- filter(rheo_CTL, rheo_CTL$Stimulus == "Stim_60s")
rheo_CTL_post <- filter(rheo_CTL, rheo_CTL$Stimulus == "Post_Stim")


rheo_CTL_ps_rad <- select(rheo_CTL_ps, Angle_rad)
rheo_CTL_10s_rad <- select(rheo_CTL_10s, Angle_rad)
rheo_CTL_20s_rad <- select(rheo_CTL_20s, Angle_rad)
rheo_CTL_30s_rad <- select(rheo_CTL_30s, Angle_rad)
rheo_CTL_40s_rad <- select(rheo_CTL_40s, Angle_rad)
rheo_CTL_50s_rad <- select(rheo_CTL_50s, Angle_rad)
rheo_CTL_60s_rad <- select(rheo_CTL_60s, Angle_rad)
rheo_CTL_post_rad <- select(rheo_CTL_post, Angle_rad)


#CUSO4: Stimulus dataframes
rheo_CuSO4 <- filter(rheotaxis_data_all, rheotaxis_data_all$Treatment == "CuSO4")
rheo_CuSO4 <- na.omit(rheo_CuSO4) #remove NA 
rheo_CuSO4_ps <- filter(rheo_CuSO4, rheo_CuSO4$Stimulus == "Pre_Stim")
rheo_CuSO4_10s <- filter(rheo_CuSO4, rheo_CuSO4$Stimulus == "Stim_10s")
rheo_CuSO4_20s <- filter(rheo_CuSO4, rheo_CuSO4$Stimulus == "Stim_20s")
rheo_CuSO4_30s <- filter(rheo_CuSO4, rheo_CuSO4$Stimulus == "Stim_30s")
rheo_CuSO4_40s <- filter(rheo_CuSO4, rheo_CuSO4$Stimulus == "Stim_40s")
rheo_CuSO4_50s <- filter(rheo_CuSO4, rheo_CuSO4$Stimulus == "Stim_50s")
rheo_CuSO4_60s <- filter(rheo_CuSO4, rheo_CuSO4$Stimulus == "Stim_60s")
rheo_CuSO4_post <- filter(rheo_CuSO4, rheo_CuSO4$Stimulus == "Post_Stim")


rheo_CuSO4_ps_rad <- select(rheo_CuSO4_ps, Angle_rad)
rheo_CuSO4_10s_rad <- select(rheo_CuSO4_10s, Angle_rad)
rheo_CuSO4_20s_rad <- select(rheo_CuSO4_20s, Angle_rad)
rheo_CuSO4_30s_rad <- select(rheo_CuSO4_30s, Angle_rad)
rheo_CuSO4_40s_rad <- select(rheo_CuSO4_40s, Angle_rad)
rheo_CuSO4_50s_rad <- select(rheo_CuSO4_50s, Angle_rad)
rheo_CuSO4_60s_rad <- select(rheo_CuSO4_60s, Angle_rad)
rheo_CuSO4_post_rad <- select(rheo_CuSO4_post, Angle_rad)


#NEOMYCIN: Stimulus dataframes
rheo_Neo <- filter(rheotaxis_data_all, rheotaxis_data_all$Treatment == "Neo")
rheo_Neo <- na.omit(rheo_Neo) #remove NA 
rheo_Neo_ps <- filter(rheo_Neo, rheo_Neo$Stimulus == "Pre_Stim")
rheo_Neo_10s <- filter(rheo_Neo, rheo_Neo$Stimulus == "Stim_10s")
rheo_Neo_20s <- filter(rheo_Neo, rheo_Neo$Stimulus == "Stim_20s")
rheo_Neo_30s <- filter(rheo_Neo, rheo_Neo$Stimulus == "Stim_30s")
rheo_Neo_40s <- filter(rheo_Neo, rheo_Neo$Stimulus == "Stim_40s")
rheo_Neo_50s <- filter(rheo_Neo, rheo_Neo$Stimulus == "Stim_50s")
rheo_Neo_60s <- filter(rheo_Neo, rheo_Neo$Stimulus == "Stim_60s")
rheo_Neo_post <- filter(rheo_Neo, rheo_Neo$Stimulus == "Post_Stim")

rheo_Neo_ps_rad <- select(rheo_Neo_ps, Angle_rad)
rheo_Neo_10s_rad <- select(rheo_Neo_10s, Angle_rad)
rheo_Neo_20s_rad <- select(rheo_Neo_20s, Angle_rad)
rheo_Neo_30s_rad <- select(rheo_Neo_30s, Angle_rad)
rheo_Neo_40s_rad <- select(rheo_Neo_40s, Angle_rad)
rheo_Neo_50s_rad <- select(rheo_Neo_50s, Angle_rad)
rheo_Neo_60s_rad <- select(rheo_Neo_60s, Angle_rad)
rheo_Neo_post_rad <- select(rheo_Neo_post, Angle_rad)



####---BAPTA SUBSETS---####

#CONTROL: Stimulus dataframes
rheo_CTL <- filter(master.all.indiv.mean1, master.all.indiv.mean1$Treatment == "Control")
rheo_CTL_ps <- filter(rheo_CTL, rheo_CTL$Stimulus == "Pre_Stim")
rheo_CTL_10s <- filter(rheo_CTL, rheo_CTL$Stimulus == "Stim_10s")
rheo_CTL_20s <- filter(rheo_CTL, rheo_CTL$Stimulus == "Stim_20s")
rheo_CTL_post <- filter(rheo_CTL, rheo_CTL$Stimulus == "Post_Stim")

rheo_CTL_ps_rad <- select(rheo_CTL_ps, Angle_rad)
rheo_CTL_10s_rad <- select(rheo_CTL_10s, Angle_rad)
rheo_CTL_20s_rad <- select(rheo_CTL_20s, Angle_rad)
rheo_CTL_post_rad <- select(rheo_CTL_post, Angle_rad)

#BAPTA subsets
rheo_bap <- filter(master.all.indiv.mean1, master.all.indiv.mean1$Treatment ==  "Bapta")
rheo_bap_ps <- filter(rheo_bap, rheo_bap$Stimulus == "Pre_Stim")
rheo_bap_10s <- filter(rheo_bap, rheo_bap$Stimulus == "Stim_10s")
rheo_bap_20s <- filter(rheo_bap, rheo_bap$Stimulus == "Stim_20s")
rheo_bap_post <- filter(rheo_bap, rheo_bap$Stimulus == "Post_Stim")

rheo_bap_ps_rad <- select(rheo_bap_ps, Angle_rad)
rheo_bap_10s_rad <- select(rheo_bap_10s, Angle_rad)
rheo_bap_20s_rad <- select(rheo_bap_20s, Angle_rad)
rheo_bap_post_rad <- select(rheo_bap_post, Angle_rad)



####---SHAKE-RECOVERY SUBSETS---####
#CONTROLS

rheo_CTL <- filter(rheotaxis_data_all, rheotaxis_data_all$Treatment == "Control")
rheo_CTL <- na.omit(rheo_CTL) #remove NA 
rheo_CTL_ps <- filter(rheo_CTL, rheo_CTL$Stimulus == "Pre_Stim")
rheo_CTL_10s <- filter(rheo_CTL, rheo_CTL$Stimulus == "Stim_10s")
rheo_CTL_20s <- filter(rheo_CTL, rheo_CTL$Stimulus == "Stim_20s")
# rheo_CTL_30s <- filter(rheo_CTL, rheo_CTL$Stimulus == "Stim_30s")
# rheo_CTL_40s <- filter(rheo_CTL, rheo_CTL$Stimulus == "Stim_40s")
# rheo_CTL_50s <- filter(rheo_CTL, rheo_CTL$Stimulus == "Stim_50s")
# rheo_CTL_60s <- filter(rheo_CTL, rheo_CTL$Stimulus == "Stim_60s")
rheo_CTL_post <- filter(rheo_CTL, rheo_CTL$Stimulus == "Post_Stim")
rheo_CTL_post <- filter(rheo_CTL, rheo_CTL$Stimulus == "Post_stim")

rheo_CTL_ps_rad <- select(rheo_CTL_ps, Angle_rad)
rheo_CTL_10s_rad <- select(rheo_CTL_10s, Angle_rad)
rheo_CTL_20s_rad <- select(rheo_CTL_20s, Angle_rad)
# rheo_CTL_30s_rad <- select(rheo_CTL_30s, Angle_rad)
# rheo_CTL_40s_rad <- select(rheo_CTL_40s, Angle_rad)
# rheo_CTL_50s_rad <- select(rheo_CTL_50s, Angle_rad)
# rheo_CTL_60s_rad <- select(rheo_CTL_60s, Angle_rad)
rheo_CTL_post_rad <- select(rheo_CTL_post, Angle_rad)


# rheo_CTL_ps0 <- filter(rheo_CTL_ps, rheo_CTL_ps$Rec_Time == "0hr")
# rheo_CTL_ps2 <- filter(rheo_CTL_ps, rheo_CTL_ps$Rec_Time == "2hr")
# rheo_CTL_ps4 <- filter(rheo_CTL_ps, rheo_CTL_ps$Rec_Time == "4hr")
# rheo_CTL_ps8 <- filter(rheo_CTL_ps, rheo_CTL_ps$Rec_Time == "8hr")
# rheo_CTL_ps48 <- filter(rheo_CTL_ps, rheo_CTL_ps$Rec_Time == "48hr")
# 
# rheo_CTL_ps0_rad <- select(rheo_CTL_ps0, Angle_rad)
# rheo_CTL_ps2_rad <- select(rheo_CTL_ps2, Angle_rad)
# rheo_CTL_ps4_rad <- select(rheo_CTL_ps4, Angle_rad)
# rheo_CTL_ps8_rad <- select(rheo_CTL_ps8, Angle_rad)
# rheo_CTL_ps48_rad <- select(rheo_CTL_ps48, Angle_rad)
# 
# rheo_CTL_10s0 <- filter(rheo_CTL_10s, rheo_CTL_10s$Rec_Time == "0hr")
# rheo_CTL_10s2 <- filter(rheo_CTL_10s, rheo_CTL_10s$Rec_Time == "2hr")
# rheo_CTL_10s4 <- filter(rheo_CTL_10s, rheo_CTL_10s$Rec_Time == "4hr")
# rheo_CTL_10s8 <- filter(rheo_CTL_10s, rheo_CTL_10s$Rec_Time == "8hr")
# rheo_CTL_10s48 <- filter(rheo_CTL_10s, rheo_CTL_10s$Rec_Time == "48hr")
# 
# rheo_CTL_10s0_rad <- select(rheo_CTL_10s0, Angle_rad)
# rheo_CTL_10s2_rad <- select(rheo_CTL_10s2, Angle_rad)
# rheo_CTL_10s4_rad <- select(rheo_CTL_10s4, Angle_rad)
# rheo_CTL_10s8_rad <- select(rheo_CTL_10s8, Angle_rad)
# rheo_CTL_10s48_rad <- select(rheo_CTL_10s48, Angle_rad)
# 
# rheo_CTL_20s0 <- filter(rheo_CTL_20s, rheo_CTL_20s$Rec_Time == "0hr")
# rheo_CTL_20s2 <- filter(rheo_CTL_20s, rheo_CTL_20s$Rec_Time == "2hr")
# rheo_CTL_20s4 <- filter(rheo_CTL_20s, rheo_CTL_20s$Rec_Time == "4hr")
# rheo_CTL_20s8 <- filter(rheo_CTL_20s, rheo_CTL_20s$Rec_Time == "8hr")
# rheo_CTL_20s48 <- filter(rheo_CTL_20s, rheo_CTL_20s$Rec_Time == "48hr")
# 
# rheo_CTL_20s0_rad <- select(rheo_CTL_20s0, Angle_rad)
# rheo_CTL_20s2_rad <- select(rheo_CTL_20s2, Angle_rad)
# rheo_CTL_20s4_rad <- select(rheo_CTL_20s4, Angle_rad)
# rheo_CTL_20s8_rad <- select(rheo_CTL_20s8, Angle_rad)
# rheo_CTL_20s48_rad <- select(rheo_CTL_20s48, Angle_rad)
# 
# rheo_CTL_post0 <- filter(rheo_CTL_post, rheo_CTL_post$Rec_Time == "0hr")
# rheo_CTL_post2 <- filter(rheo_CTL_post, rheo_CTL_post$Rec_Time == "2hr")
# rheo_CTL_post4 <- filter(rheo_CTL_post, rheo_CTL_post$Rec_Time == "4hr")
# rheo_CTL_post8 <- filter(rheo_CTL_post, rheo_CTL_post$Rec_Time == "8hr")
# rheo_CTL_post48 <- filter(rheo_CTL_post, rheo_CTL_post$Rec_Time == "48hr")
# 
# rheo_CTL_post0_rad <- select(rheo_CTL_post0, Angle_rad)
# rheo_CTL_post2_rad <- select(rheo_CTL_post2, Angle_rad)
# rheo_CTL_post4_rad <- select(rheo_CTL_post4, Angle_rad)
# rheo_CTL_post8_rad <- select(rheo_CTL_post8, Angle_rad)
# rheo_CTL_post48_rad <- select(rheo_CTL_post48, Angle_rad)
# 
# 
# #Shake-Recovery STIMULUS subsets

rheo_shake <- filter(rheotaxis_data_all, rheotaxis_data_all$Treatment == "Shake")
rheo_shake <- na.omit(rheo_shake) #remove NA 
rheo_shake_ps <- filter(rheo_shake, rheo_shake$Stimulus == "Pre_Stim")
rheo_shake_10s <- filter(rheo_shake, rheo_shake$Stimulus == "Stim_10s")
rheo_shake_20s <- filter(rheo_shake, rheo_shake$Stimulus == "Stim_20s")
# rheo_shake_30s <- filter(rheo_shake, rheo_shake$Stimulus == "Stim_30s")
# rheo_shake_40s <- filter(rheo_shake, rheo_shake$Stimulus == "Stim_40s")
# rheo_shake_50s <- filter(rheo_shake, rheo_shake$Stimulus == "Stim_50s")
# rheo_shake_60s <- filter(rheo_shake, rheo_shake$Stimulus == "Stim_60s")
rheo_shake_post <- filter(rheo_shake, rheo_shake$Stimulus == "Post_Stim")
rheo_shake_post <- filter(rheo_shake, rheo_shake$Stimulus == "Post_stim")

rheo_shake_ps_rad <- select(rheo_shake_ps, Angle_rad)
rheo_shake_10s_rad <- select(rheo_shake_10s, Angle_rad)
rheo_shake_20s_rad <- select(rheo_shake_20s, Angle_rad)
# rheo_shake_30s_rad <- select(rheo_shake_30s, Angle_rad)
# rheo_shake_40s_rad <- select(rheo_shake_40s, Angle_rad)
# rheo_shake_50s_rad <- select(rheo_shake_50s, Angle_rad)
# rheo_shake_60s_rad <- select(rheo_shake_60s, Angle_rad)
rheo_shake_post_rad <- select(rheo_shake_post, Angle_rad)
# 
# 
# rheo_shake_ps0 <- filter(rheo_shake_ps, rheo_shake_ps$Rec_Time == "0hr")
# rheo_shake_ps2 <- filter(rheo_shake_ps, rheo_shake_ps$Rec_Time == "2hr")
# rheo_shake_ps4 <- filter(rheo_shake_ps, rheo_shake_ps$Rec_Time == "4hr")
# rheo_shake_ps8 <- filter(rheo_shake_ps, rheo_shake_ps$Rec_Time == "8hr")
# rheo_shake_ps48 <- filter(rheo_shake_ps, rheo_shake_ps$Rec_Time == "48hr")
# 
# rheo_shake_ps0_rad <- select(rheo_shake_ps0, Angle_rad)
# rheo_shake_ps2_rad <- select(rheo_shake_ps2, Angle_rad)
# rheo_shake_ps4_rad <- select(rheo_shake_ps4, Angle_rad)
# rheo_shake_ps8_rad <- select(rheo_shake_ps8, Angle_rad)
# rheo_shake_ps48_rad <- select(rheo_shake_ps48, Angle_rad)
# 
# rheo_shake_10s0 <- filter(rheo_shake_10s, rheo_shake_10s$Rec_Time == "0hr")
# rheo_shake_10s2 <- filter(rheo_shake_10s, rheo_shake_10s$Rec_Time == "2hr")
# rheo_shake_10s4 <- filter(rheo_shake_10s, rheo_shake_10s$Rec_Time == "4hr")
# rheo_shake_10s8 <- filter(rheo_shake_10s, rheo_shake_10s$Rec_Time == "8hr")
# rheo_shake_10s48 <- filter(rheo_shake_10s, rheo_shake_10s$Rec_Time == "48hr")
# 
# rheo_shake_10s0_rad <- select(rheo_shake_10s0, Angle_rad)
# rheo_shake_10s2_rad <- select(rheo_shake_10s2, Angle_rad)
# rheo_shake_10s4_rad <- select(rheo_shake_10s4, Angle_rad)
# rheo_shake_10s8_rad <- select(rheo_shake_10s8, Angle_rad)
# rheo_shake_10s48_rad <- select(rheo_shake_10s48, Angle_rad)
# 
# rheo_shake_20s0 <- filter(rheo_shake_20s, rheo_shake_20s$Rec_Time == "0hr")
# rheo_shake_20s2 <- filter(rheo_shake_20s, rheo_shake_20s$Rec_Time == "2hr")
# rheo_shake_20s4 <- filter(rheo_shake_20s, rheo_shake_20s$Rec_Time == "4hr")
# rheo_shake_20s8 <- filter(rheo_shake_20s, rheo_shake_20s$Rec_Time == "8hr")
# rheo_shake_20s48 <- filter(rheo_shake_20s, rheo_shake_20s$Rec_Time == "48hr")
# 
# rheo_shake_20s0_rad <- select(rheo_shake_20s0, Angle_rad)
# rheo_shake_20s2_rad <- select(rheo_shake_20s2, Angle_rad)
# rheo_shake_20s4_rad <- select(rheo_shake_20s4, Angle_rad)
# rheo_shake_20s8_rad <- select(rheo_shake_20s8, Angle_rad)
# rheo_shake_20s48_rad <- select(rheo_shake_20s48, Angle_rad)
# 
# rheo_shake_post0 <- filter(rheo_shake_post, rheo_shake_post$Rec_Time == "0hr")
# rheo_shake_post2 <- filter(rheo_shake_post, rheo_shake_post$Rec_Time == "2hr")
# rheo_shake_post4 <- filter(rheo_shake_post, rheo_shake_post$Rec_Time == "4hr")
# rheo_shake_post8 <- filter(rheo_shake_post, rheo_shake_post$Rec_Time == "8hr")
# rheo_shake_post48 <- filter(rheo_shake_post, rheo_shake_post$Rec_Time == "48hr")
# 
# rheo_shake_post0_rad <- select(rheo_shake_post0, Angle_rad)
# rheo_shake_post2_rad <- select(rheo_shake_post2, Angle_rad)
# rheo_shake_post4_rad <- select(rheo_shake_post4, Angle_rad)
# rheo_shake_post8_rad <- select(rheo_shake_post8, Angle_rad)
# rheo_shake_post48_rad <- select(rheo_shake_post48, Angle_rad)



library(circular)


####---CUSO4 -NEOMYCIN GRAPHS ---####

setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/Control_Pre_rheo.eps")
mean.circular(rheo_CTL_ps_rad)*180/pi
mean.rheo_CTL_ps_deg <- mean.circular(rheo_CTL_ps_rad)*180/pi
mean.circular(rheo_CTL_ps_rad)
mean.rheo_CTL_ps_rad <- mean.circular(rheo_CTL_ps_rad)
rho.circular(rheo_CTL_ps_rad)
rho.rheo_CTL_ps_rad <- rho.circular(rheo_CTL_ps_rad)
sd.circular(rheo_CTL_ps_rad)*180/pi
sd.rheo_CTL_ps_rad <- sd.circular(rheo_CTL_ps_rad)*180/pi
plot_CTL_ps <- plot.circular(rheo_CTL_ps_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: Pre-Stimulus", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_ps_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_ps_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/Control_10s_rheo.eps")  
mean.circular(rheo_CTL_10s_rad)*180/pi
mean.rheo_CTL_10s_deg <- mean.circular(rheo_CTL_10s_rad)*180/pi
mean.circular(rheo_CTL_10s_rad)
mean.rheo_CTL_10s_rad <- mean.circular(rheo_CTL_10s_rad)
rho.circular(rheo_CTL_10s_rad)
rho.rheo_CTL_10s_rad <- rho.circular(rheo_CTL_10s_rad)
sd.circular(rheo_CTL_10s_rad)*180/pi
sd.rheo_CTL_10s_rad <- sd.circular(rheo_CTL_10s_rad)*180/pi
plot_CTL_10s <- plot.circular(rheo_CTL_10s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 10s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_10s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_10s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/Control_20s_rheo.eps")  
mean.circular(rheo_CTL_20s_rad)*180/pi
mean.rheo_CTL_20s_deg <- mean.circular(rheo_CTL_20s_rad)*180/pi
mean.circular(rheo_CTL_20s_rad)
mean.rheo_CTL_20s_rad <- mean.circular(rheo_CTL_20s_rad)
rho.circular(rheo_CTL_20s_rad)
rho.rheo_CTL_20s_rad <- rho.circular(rheo_CTL_20s_rad)
sd.circular(rheo_CTL_20s_rad)*180/pi
sd.rheo_CTL_20s_rad <- sd.circular(rheo_CTL_20s_rad)*180/pi
plot_CTL_20s <- plot.circular(rheo_CTL_20s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 20s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_20s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_20s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/Control_30s_rheo.eps")  
mean.circular(rheo_CTL_30s_rad)*180/pi
mean.rheo_CTL_30s_deg <- mean.circular(rheo_CTL_30s_rad)*180/pi
mean.circular(rheo_CTL_30s_rad)
mean.rheo_CTL_30s_rad <- mean.circular(rheo_CTL_30s_rad)
rho.circular(rheo_CTL_30s_rad)
rho.rheo_CTL_30s_rad <- rho.circular(rheo_CTL_30s_rad)
sd.circular(rheo_CTL_30s_rad)*180/pi
sd.rheo_CTL_30s_rad <- sd.circular(rheo_CTL_30s_rad)*180/pi
plot_CTL_30s <- plot.circular(rheo_CTL_30s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 30s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_30s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_30s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/Control_40s_rheo.eps")  
mean.circular(rheo_CTL_40s_rad)*180/pi
mean.rheo_CTL_40s_deg <- mean.circular(rheo_CTL_40s_rad)*180/pi
mean.circular(rheo_CTL_40s_rad)
mean.rheo_CTL_40s_rad <- mean.circular(rheo_CTL_40s_rad)
rho.circular(rheo_CTL_40s_rad)
rho.rheo_CTL_40s_rad <- rho.circular(rheo_CTL_40s_rad)
sd.circular(rheo_CTL_40s_rad)*180/pi
sd.rheo_CTL_40s_rad <- sd.circular(rheo_CTL_40s_rad)*180/pi
plot_CTL_40s <- plot.circular(rheo_CTL_40s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 40s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_40s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_40s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/Control_50s_rheo.eps")  
mean.circular(rheo_CTL_50s_rad)*180/pi
mean.rheo_CTL_50s_deg <- mean.circular(rheo_CTL_50s_rad)*180/pi
mean.circular(rheo_CTL_50s_rad)
mean.rheo_CTL_50s_rad <- mean.circular(rheo_CTL_50s_rad)
rho.circular(rheo_CTL_50s_rad)
rho.rheo_CTL_50s_rad <- rho.circular(rheo_CTL_50s_rad)
sd.circular(rheo_CTL_50s_rad)*180/pi
sd.rheo_CTL_50s_rad <- sd.circular(rheo_CTL_50s_rad)*180/pi
plot_CTL_50s <- plot.circular(rheo_CTL_50s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 50s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_50s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_50s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/Control_60s_rheo.eps")  
mean.circular(rheo_CTL_60s_rad)*180/pi
mean.rheo_CTL_60s_deg <- mean.circular(rheo_CTL_60s_rad)*180/pi
mean.circular(rheo_CTL_60s_rad)
mean.rheo_CTL_60s_rad <- mean.circular(rheo_CTL_60s_rad)
rho.circular(rheo_CTL_60s_rad)
rho.rheo_CTL_60s_rad <- rho.circular(rheo_CTL_60s_rad)
sd.circular(rheo_CTL_60s_rad)*180/pi
sd.rheo_CTL_60s_rad <- sd.circular(rheo_CTL_60s_rad)*180/pi
plot_CTL_60s <- plot.circular(rheo_CTL_60s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 60s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_60s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_60s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/Control_Post_rheo.eps")  
mean.circular(rheo_CTL_post_rad)*180/pi
mean.rheo_CTL_post_deg <- mean.circular(rheo_CTL_post_rad)*180/pi
mean.circular(rheo_CTL_post_rad)
mean.rheo_CTL_post_rad <- mean.circular(rheo_CTL_post_rad)
rho.circular(rheo_CTL_post_rad)
rho.rheo_CTL_post_rad <- rho.circular(rheo_CTL_post_rad)
sd.circular(rheo_CTL_post_rad)*180/pi
sd.rheo_CTL_post_rad <- sd.circular(rheo_CTL_post_rad)*180/pi
plot_CTL_post <- plot.circular(rheo_CTL_post_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                               start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                               tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                               xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                               template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: Post-Stimulus", 
                               sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_post_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_post_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd



setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/CuSO4_Pre_rheo.eps")
mean.circular(rheo_CuSO4_ps_rad)*180/pi
mean.rheo_CuSO4_ps_deg <- mean.circular(rheo_CuSO4_ps_rad)*180/pi
mean.circular(rheo_CuSO4_ps_rad)
mean.rheo_CuSO4_ps_rad <- mean.circular(rheo_CuSO4_ps_rad)
rho.circular(rheo_CuSO4_ps_rad)
rho.rheo_CuSO4_ps_rad <- rho.circular(rheo_CuSO4_ps_rad)
sd.circular(rheo_CuSO4_ps_rad)*180/pi
sd.rheo_CuSO4_ps_rad <- sd.circular(rheo_CuSO4_ps_rad)*180/pi
plot_CuSO4_ps <- plot.circular(rheo_CuSO4_ps_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "CuSO4 Fish: Pre-Stimulus", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CuSO4_ps_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CuSO4_ps_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/CuSO4_10s_rheo.eps")  
mean.circular(rheo_CuSO4_10s_rad)*180/pi
mean.rheo_CuSO4_10s_deg <- mean.circular(rheo_CuSO4_10s_rad)*180/pi
mean.circular(rheo_CuSO4_10s_rad)
mean.rheo_CuSO4_10s_rad <- mean.circular(rheo_CuSO4_10s_rad)
rho.circular(rheo_CuSO4_10s_rad)
rho.rheo_CuSO4_10s_rad <- rho.circular(rheo_CuSO4_10s_rad)
sd.circular(rheo_CuSO4_10s_rad)*180/pi
sd.rheo_CuSO4_10s_rad <- sd.circular(rheo_CuSO4_10s_rad)*180/pi
plot_CuSO4_10s <- plot.circular(rheo_CuSO4_10s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "CuSO4 Fish: 10s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CuSO4_10s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CuSO4_10s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/CuSO4_20s_rheo.eps")  
mean.circular(rheo_CuSO4_20s_rad)*180/pi
mean.rheo_CuSO4_20s_deg <- mean.circular(rheo_CuSO4_20s_rad)*180/pi
mean.circular(rheo_CuSO4_20s_rad)
mean.rheo_CuSO4_20s_rad <- mean.circular(rheo_CuSO4_20s_rad)
rho.circular(rheo_CuSO4_20s_rad)
rho.rheo_CuSO4_20s_rad <- rho.circular(rheo_CuSO4_20s_rad)
sd.circular(rheo_CuSO4_20s_rad)*180/pi
sd.rheo_CuSO4_20s_rad <- sd.circular(rheo_CuSO4_20s_rad)*180/pi
plot_CuSO4_20s <- plot.circular(rheo_CuSO4_20s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "CuSO4 Fish: 20s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CuSO4_20s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CuSO4_20s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/CuSO4_30s_rheo.eps")  
mean.circular(rheo_CuSO4_30s_rad)*180/pi
mean.rheo_CuSO4_30s_deg <- mean.circular(rheo_CuSO4_30s_rad)*180/pi
mean.circular(rheo_CuSO4_30s_rad)
mean.rheo_CuSO4_30s_rad <- mean.circular(rheo_CuSO4_30s_rad)
rho.circular(rheo_CuSO4_30s_rad)
rho.rheo_CuSO4_30s_rad <- rho.circular(rheo_CuSO4_30s_rad)
sd.circular(rheo_CuSO4_30s_rad)*180/pi
sd.rheo_CuSO4_30s_rad <- sd.circular(rheo_CuSO4_30s_rad)*180/pi
plot_CuSO4_30s <- plot.circular(rheo_CuSO4_30s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                                start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                                tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                                xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                                template = NULL, zero = pi/2, rotation = "clock", main = "CuSO4 Fish: 30s-Stimulus", 
                                sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CuSO4_30s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CuSO4_30s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/CuSO4_40s_rheo.eps")  
mean.circular(rheo_CuSO4_40s_rad)*180/pi
mean.rheo_CuSO4_40s_deg <- mean.circular(rheo_CuSO4_40s_rad)*180/pi
mean.circular(rheo_CuSO4_40s_rad)
mean.rheo_CuSO4_40s_rad <- mean.circular(rheo_CuSO4_40s_rad)
rho.circular(rheo_CuSO4_40s_rad)
rho.rheo_CuSO4_40s_rad <- rho.circular(rheo_CuSO4_40s_rad)
sd.circular(rheo_CuSO4_40s_rad)*180/pi
sd.rheo_CuSO4_40s_rad <- sd.circular(rheo_CuSO4_40s_rad)*180/pi
plot_CuSO4_40s <- plot.circular(rheo_CuSO4_40s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                                start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                                tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                                xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                                template = NULL, zero = pi/2, rotation = "clock", main = "CuSO4 Fish: 40s-Stimulus", 
                                sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CuSO4_40s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CuSO4_40s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/CuSO4_50s_rheo.eps")  
mean.circular(rheo_CuSO4_50s_rad)*180/pi
mean.rheo_CuSO4_50s_deg <- mean.circular(rheo_CuSO4_50s_rad)*180/pi
mean.circular(rheo_CuSO4_50s_rad)
mean.rheo_CuSO4_50s_rad <- mean.circular(rheo_CuSO4_50s_rad)
rho.circular(rheo_CuSO4_50s_rad)
rho.rheo_CuSO4_50s_rad <- rho.circular(rheo_CuSO4_50s_rad)
sd.circular(rheo_CuSO4_50s_rad)*180/pi
sd.rheo_CuSO4_50s_rad <- sd.circular(rheo_CuSO4_50s_rad)*180/pi
plot_CuSO4_50s <- plot.circular(rheo_CuSO4_50s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                                start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                                tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                                xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                                template = NULL, zero = pi/2, rotation = "clock", main = "CuSO4 Fish: 50s-Stimulus", 
                                sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CuSO4_50s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CuSO4_50s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/CuSO4_60s_rheo.eps")  
mean.circular(rheo_CuSO4_60s_rad)*180/pi
mean.rheo_CuSO4_60s_deg <- mean.circular(rheo_CuSO4_60s_rad)*180/pi
mean.circular(rheo_CuSO4_60s_rad)
mean.rheo_CuSO4_60s_rad <- mean.circular(rheo_CuSO4_60s_rad)
rho.circular(rheo_CuSO4_60s_rad)
rho.rheo_CuSO4_60s_rad <- rho.circular(rheo_CuSO4_60s_rad)
sd.circular(rheo_CuSO4_60s_rad)*180/pi
sd.rheo_CuSO4_60s_rad <- sd.circular(rheo_CuSO4_60s_rad)*180/pi
plot_CuSO4_60s <- plot.circular(rheo_CuSO4_60s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                                start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                                tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                                xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                                template = NULL, zero = pi/2, rotation = "clock", main = "CuSO4 Fish: 60s-Stimulus", 
                                sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CuSO4_60s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CuSO4_60s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/CuSO4_Post_rheo.eps")  
mean.circular(rheo_CuSO4_post_rad)*180/pi
mean.rheo_CuSO4_post_deg <- mean.circular(rheo_CuSO4_post_rad)*180/pi
mean.circular(rheo_CuSO4_post_rad)
mean.rheo_CuSO4_post_rad <- mean.circular(rheo_CuSO4_post_rad)
rho.circular(rheo_CuSO4_post_rad)
rho.rheo_CuSO4_post_rad <- rho.circular(rheo_CuSO4_post_rad)
sd.circular(rheo_CuSO4_post_rad)*180/pi
sd.rheo_CuSO4_post_rad <- sd.circular(rheo_CuSO4_post_rad)*180/pi
plot_CuSO4_post <- plot.circular(rheo_CuSO4_post_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                               start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                               tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                               xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                               template = NULL, zero = pi/2, rotation = "clock", main = "CuSO4 Fish: Post-Stimulus", 
                               sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CuSO4_post_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CuSO4_post_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd




setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/Neo_Pre_rheo.eps")
mean.circular(rheo_Neo_ps_rad)*180/pi
mean.rheo_Neo_ps_deg <- mean.circular(rheo_Neo_ps_rad)*180/pi
mean.circular(rheo_Neo_ps_rad)
mean.rheo_Neo_ps_rad <- mean.circular(rheo_Neo_ps_rad)
rho.circular(rheo_Neo_ps_rad)
rho.rheo_Neo_ps_rad <- rho.circular(rheo_Neo_ps_rad)
sd.circular(rheo_Neo_ps_rad)*180/pi
sd.rheo_Neo_ps_rad <- sd.circular(rheo_Neo_ps_rad)*180/pi
plot_Neo_ps <- plot.circular(rheo_Neo_ps_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Neo Fish: Pre-Stimulus", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_Neo_ps_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_Neo_ps_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/Neo_10s_rheo.eps")  
mean.circular(rheo_Neo_10s_rad)*180/pi
mean.rheo_Neo_10s_deg <- mean.circular(rheo_Neo_10s_rad)*180/pi
mean.circular(rheo_Neo_10s_rad)
mean.rheo_Neo_10s_rad <- mean.circular(rheo_Neo_10s_rad)
rho.circular(rheo_Neo_10s_rad)
rho.rheo_Neo_10s_rad <- rho.circular(rheo_Neo_10s_rad)
sd.circular(rheo_Neo_10s_rad)*180/pi
sd.rheo_Neo_10s_rad <- sd.circular(rheo_Neo_10s_rad)*180/pi
plot_Neo_10s <- plot.circular(rheo_Neo_10s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Neo Fish: 10s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_Neo_10s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_Neo_10s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/Neo_20s_rheo.eps")  
mean.circular(rheo_Neo_20s_rad)*180/pi
mean.rheo_Neo_20s_deg <- mean.circular(rheo_Neo_20s_rad)*180/pi
mean.circular(rheo_Neo_20s_rad)
mean.rheo_Neo_20s_rad <- mean.circular(rheo_Neo_20s_rad)
rho.circular(rheo_Neo_20s_rad)
rho.rheo_Neo_20s_rad <- rho.circular(rheo_Neo_20s_rad)
sd.circular(rheo_Neo_20s_rad)*180/pi
sd.rheo_Neo_20s_rad <- sd.circular(rheo_Neo_20s_rad)*180/pi
plot_Neo_20s <- plot.circular(rheo_Neo_20s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Neo Fish: 20s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_Neo_20s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_Neo_20s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/Neo_30s_rheo.eps")  
mean.circular(rheo_Neo_30s_rad)*180/pi
mean.rheo_Neo_30s_deg <- mean.circular(rheo_Neo_30s_rad)*180/pi
mean.circular(rheo_Neo_30s_rad)
mean.rheo_Neo_30s_rad <- mean.circular(rheo_Neo_30s_rad)
rho.circular(rheo_Neo_30s_rad)
rho.rheo_Neo_30s_rad <- rho.circular(rheo_Neo_30s_rad)
sd.circular(rheo_Neo_30s_rad)*180/pi
sd.rheo_Neo_30s_rad <- sd.circular(rheo_Neo_30s_rad)*180/pi
plot_Neo_30s <- plot.circular(rheo_Neo_30s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Neo Fish: 30s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_Neo_30s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_Neo_30s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/Neo_40s_rheo.eps")  
mean.circular(rheo_Neo_40s_rad)*180/pi
mean.rheo_Neo_40s_deg <- mean.circular(rheo_Neo_40s_rad)*180/pi
mean.circular(rheo_Neo_40s_rad)
mean.rheo_Neo_40s_rad <- mean.circular(rheo_Neo_40s_rad)
rho.circular(rheo_Neo_40s_rad)
rho.rheo_Neo_40s_rad <- rho.circular(rheo_Neo_40s_rad)
sd.circular(rheo_Neo_40s_rad)*180/pi
sd.rheo_Neo_40s_rad <- sd.circular(rheo_Neo_40s_rad)*180/pi
plot_Neo_40s <- plot.circular(rheo_Neo_40s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Neo Fish: 40s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_Neo_40s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_Neo_40s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/Neo_50s_rheo.eps")  
mean.circular(rheo_Neo_50s_rad)*180/pi
mean.rheo_Neo_50s_deg <- mean.circular(rheo_Neo_50s_rad)*180/pi
mean.circular(rheo_Neo_50s_rad)
mean.rheo_Neo_50s_rad <- mean.circular(rheo_Neo_50s_rad)
rho.circular(rheo_Neo_50s_rad)
rho.rheo_Neo_50s_rad <- rho.circular(rheo_Neo_50s_rad)
sd.circular(rheo_Neo_50s_rad)*180/pi
sd.rheo_Neo_50s_rad <- sd.circular(rheo_Neo_50s_rad)*180/pi
plot_Neo_50s <- plot.circular(rheo_Neo_50s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Neo Fish: 50s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_Neo_50s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_Neo_50s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/Neo_60s_rheo.eps")  
mean.circular(rheo_Neo_60s_rad)*180/pi
mean.rheo_Neo_60s_deg <- mean.circular(rheo_Neo_60s_rad)*180/pi
mean.circular(rheo_Neo_60s_rad)
mean.rheo_Neo_60s_rad <- mean.circular(rheo_Neo_60s_rad)
rho.circular(rheo_Neo_60s_rad)
rho.rheo_Neo_60s_rad <- rho.circular(rheo_Neo_60s_rad)
sd.circular(rheo_Neo_60s_rad)*180/pi
sd.rheo_Neo_60s_rad <- sd.circular(rheo_Neo_60s_rad)*180/pi
plot_Neo_60s <- plot.circular(rheo_Neo_60s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Neo Fish: 60s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_Neo_60s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_Neo_60s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/Neo_Post_rheo.eps")  
mean.circular(rheo_Neo_post_rad)*180/pi
mean.rheo_Neo_post_deg <- mean.circular(rheo_Neo_post_rad)*180/pi
mean.circular(rheo_Neo_post_rad)
mean.rheo_Neo_post_rad <- mean.circular(rheo_Neo_post_rad)
rho.circular(rheo_Neo_post_rad)
rho.rheo_Neo_post_rad <- rho.circular(rheo_Neo_post_rad)
sd.circular(rheo_Neo_post_rad)*180/pi
sd.rheo_Neo_post_rad <- sd.circular(rheo_Neo_post_rad)*180/pi
plot_Neo_post <- plot.circular(rheo_Neo_post_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                               start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                               tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                               xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                               template = NULL, zero = pi/2, rotation = "clock", main = "Neo Fish: Post-Stimulus", 
                               sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_Neo_post_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_Neo_post_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 
# ####--BAPTA- GRAPHS---####
# 
# setEPS()                                             # Set postscript arguments
# postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/Control_Pre_rheo.eps")
# mean.circular(rheo_CTL_ps_rad)*180/pi
# mean.rheo_CTL_ps_deg <- mean.circular(rheo_CTL_ps_rad)*180/pi
# mean.circular(rheo_CTL_ps_rad)
# mean.rheo_CTL_ps_rad <- mean.circular(rheo_CTL_ps_rad)
# rho.circular(rheo_CTL_ps_rad)
# rho.rheo_CTL_ps_rad <- rho.circular(rheo_CTL_ps_rad)
# sd.circular(rheo_CTL_ps_rad)*180/pi
# sd.rheo_CTL_ps_rad <- sd.circular(rheo_CTL_ps_rad)*180/pi
# plot_CTL_ps <- plot.circular(rheo_CTL_ps_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
#                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
#                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
#                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
#                              template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: Pre-Stimulus", 
#                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
# arrows.circular(mean.rheo_CTL_ps_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
#                 shrink = rho.rheo_CTL_ps_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
# dev.off() #expot eps file to setwd
# 
# 
# setEPS()                                             # Set postscript arguments
# postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/Control_10s_rheo.eps")  
# mean.circular(rheo_CTL_10s_rad)*180/pi
# mean.rheo_CTL_10s_deg <- mean.circular(rheo_CTL_10s_rad)*180/pi
# mean.circular(rheo_CTL_10s_rad)
# mean.rheo_CTL_10s_rad <- mean.circular(rheo_CTL_10s_rad)
# rho.circular(rheo_CTL_10s_rad)
# rho.rheo_CTL_10s_rad <- rho.circular(rheo_CTL_10s_rad)
# sd.circular(rheo_CTL_10s_rad)*180/pi
# sd.rheo_CTL_10s_rad <- sd.circular(rheo_CTL_10s_rad)*180/pi
# plot_CTL_10s <- plot.circular(rheo_CTL_10s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
#                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
#                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
#                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
#                              template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 10s-Stimulus", 
#                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
# arrows.circular(mean.rheo_CTL_10s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
#                 shrink = rho.rheo_CTL_10s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
# dev.off() #expot eps file to setwd
# 
# 
# setEPS()                                             # Set postscript arguments
# postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/Control_20s_rheo.eps")  
# mean.circular(rheo_CTL_20s_rad)*180/pi
# mean.rheo_CTL_20s_deg <- mean.circular(rheo_CTL_20s_rad)*180/pi
# mean.circular(rheo_CTL_20s_rad)
# mean.rheo_CTL_20s_rad <- mean.circular(rheo_CTL_20s_rad)
# rho.circular(rheo_CTL_20s_rad)
# rho.rheo_CTL_20s_rad <- rho.circular(rheo_CTL_20s_rad)
# sd.circular(rheo_CTL_20s_rad)*180/pi
# sd.rheo_CTL_20s_rad <- sd.circular(rheo_CTL_20s_rad)*180/pi
# plot_CTL_20s <- plot.circular(rheo_CTL_20s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
#                               start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
#                               tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
#                               xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
#                               template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 20s-Stimulus", 
#                               sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
# arrows.circular(mean.rheo_CTL_20s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
#                 shrink = rho.rheo_CTL_20s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
# dev.off() #expot eps file to setwd
# 
# 
# setEPS()                                             # Set postscript arguments
# postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/Control_Post_rheo.eps")  
# mean.circular(rheo_CTL_post_rad)*180/pi
# mean.rheo_CTL_post_deg <- mean.circular(rheo_CTL_post_rad)*180/pi
# mean.circular(rheo_CTL_post_rad)
# mean.rheo_CTL_post_rad <- mean.circular(rheo_CTL_post_rad)
# rho.circular(rheo_CTL_post_rad)
# rho.rheo_CTL_post_rad <- rho.circular(rheo_CTL_post_rad)
# sd.circular(rheo_CTL_post_rad)*180/pi
# sd.rheo_CTL_post_rad <- sd.circular(rheo_CTL_post_rad)*180/pi
# plot_CTL_post <- plot.circular(rheo_CTL_post_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
#                               start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
#                               tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
#                               xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
#                               template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: Post-Stimulus", 
#                               sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
# arrows.circular(mean.rheo_CTL_post_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
#                 shrink = rho.rheo_CTL_post_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
# dev.off() #expot eps file to setwd
# 
# 
# 
# setEPS()                                             # Set postscript arguments
# postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/Bapta_Pre_rheo.eps")
# mean.circular(rheo_bap_ps_rad)*180/pi
# mean.rheo_bap_ps_deg <- mean.circular(rheo_bap_ps_rad)*180/pi
# mean.circular(rheo_bap_ps_rad)
# mean.rheo_bap_ps_rad <- mean.circular(rheo_bap_ps_rad)
# rho.circular(rheo_bap_ps_rad)
# rho.rheo_bap_ps_rad <- rho.circular(rheo_bap_ps_rad)
# sd.circular(rheo_bap_ps_rad)*180/pi
# sd.rheo_bap_ps_rad <- sd.circular(rheo_bap_ps_rad)*180/pi
# plot_bap_ps <- plot.circular(rheo_bap_ps_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
#                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
#                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
#                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
#                              template = NULL, zero = pi/2, rotation = "clock", main = "Bapta Fish: Pre-Stimulus", 
#                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
# arrows.circular(mean.rheo_bap_ps_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
#                 shrink = rho.rheo_bap_ps_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
# dev.off() #expot eps file to setwd
# 
# 
# setEPS()                                             # Set postscript arguments
# postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/Bapta_10s_rheo.eps")  
# mean.circular(rheo_bap_10s_rad)*180/pi
# mean.rheo_bap_10s_deg <- mean.circular(rheo_bap_10s_rad)*180/pi
# mean.circular(rheo_bap_10s_rad)
# mean.rheo_bap_10s_rad <- mean.circular(rheo_bap_10s_rad)
# rho.circular(rheo_bap_10s_rad)
# rho.rheo_bap_10s_rad <- rho.circular(rheo_bap_10s_rad)
# sd.circular(rheo_bap_10s_rad)*180/pi
# sd.rheo_bap_10s_rad <- sd.circular(rheo_bap_10s_rad)*180/pi
# plot_bap_10s <- plot.circular(rheo_bap_10s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
#                               start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
#                               tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
#                               xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
#                               template = NULL, zero = pi/2, rotation = "clock", main = "Bapta Fish: 10s-Stimulus", 
#                               sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
# arrows.circular(mean.rheo_bap_10s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
#                 shrink = rho.rheo_bap_10s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
# dev.off() #expot eps file to setwd
# 
# 
# setEPS()                                             # Set postscript arguments
# postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/Bapta_20s_rheo.eps")  
# mean.circular(rheo_bap_20s_rad)*180/pi
# mean.rheo_bap_20s_deg <- mean.circular(rheo_bap_20s_rad)*180/pi
# mean.circular(rheo_bap_20s_rad)
# mean.rheo_bap_20s_rad <- mean.circular(rheo_bap_20s_rad)
# rho.circular(rheo_bap_20s_rad)
# rho.rheo_bap_20s_rad <- rho.circular(rheo_bap_20s_rad)
# sd.circular(rheo_bap_20s_rad)*180/pi
# sd.rheo_bap_20s_rad <- sd.circular(rheo_bap_20s_rad)*180/pi
# plot_bap_20s <- plot.circular(rheo_bap_20s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
#                               start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
#                               tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
#                               xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
#                               template = NULL, zero = pi/2, rotation = "clock", main = "Bapta Fish: 20s-Stimulus", 
#                               sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
# arrows.circular(mean.rheo_bap_20s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
#                 shrink = rho.rheo_bap_20s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
# dev.off() #expot eps file to setwd
# 
# 
# setEPS()                                             # Set postscript arguments
# postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/Bapta_Post_rheo.eps")  
# mean.circular(rheo_bap_post_rad)*180/pi
# mean.rheo_bap_post_deg <- mean.circular(rheo_bap_post_rad)*180/pi
# mean.circular(rheo_bap_post_rad)
# mean.rheo_bap_post_rad <- mean.circular(rheo_bap_post_rad)
# rho.circular(rheo_bap_post_rad)
# rho.rheo_bap_post_rad <- rho.circular(rheo_bap_post_rad)
# sd.circular(rheo_bap_post_rad)*180/pi
# sd.rheo_bap_post_rad <- sd.circular(rheo_bap_post_rad)*180/pi
# plot_bap_post <- plot.circular(rheo_bap_post_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
#                                start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
#                                tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
#                                xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
#                                template = NULL, zero = pi/2, rotation = "clock", main = "Bapta Fish: Post-Stimulus", 
#                                sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
# arrows.circular(mean.rheo_bap_post_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
#                 shrink = rho.rheo_bap_post_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
# dev.off() #expot eps file to setwd
# 


####---SHAKE -GRAPHS---####


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/Control_Pre_rheo.eps")
mean.circular(rheo_CTL_ps_rad)*180/pi
mean.rheo_CTL_ps_deg <- mean.circular(rheo_CTL_ps_rad)*180/pi
mean.circular(rheo_CTL_ps_rad)
mean.rheo_CTL_ps_rad <- mean.circular(rheo_CTL_ps_rad)
rho.circular(rheo_CTL_ps_rad)
rho.rheo_CTL_ps_rad <- rho.circular(rheo_CTL_ps_rad)
sd.circular(rheo_CTL_ps_rad)*180/pi
sd.rheo_CTL_ps_rad <- sd.circular(rheo_CTL_ps_rad)*180/pi
plot_CTL_ps <- plot.circular(rheo_CTL_ps_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: Pre-Stimulus", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_ps_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_ps_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/Control_10s_rheo.eps")  
mean.circular(rheo_CTL_10s_rad)*180/pi
mean.rheo_CTL_10s_deg <- mean.circular(rheo_CTL_10s_rad)*180/pi
mean.circular(rheo_CTL_10s_rad)
mean.rheo_CTL_10s_rad <- mean.circular(rheo_CTL_10s_rad)
rho.circular(rheo_CTL_10s_rad)
rho.rheo_CTL_10s_rad <- rho.circular(rheo_CTL_10s_rad)
sd.circular(rheo_CTL_10s_rad)*180/pi
sd.rheo_CTL_10s_rad <- sd.circular(rheo_CTL_10s_rad)*180/pi
plot_CTL_10s <- plot.circular(rheo_CTL_10s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 10s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_10s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_10s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/Control_20s_rheo.eps")  
mean.circular(rheo_CTL_20s_rad)*180/pi
mean.rheo_CTL_20s_deg <- mean.circular(rheo_CTL_20s_rad)*180/pi
mean.circular(rheo_CTL_20s_rad)
mean.rheo_CTL_20s_rad <- mean.circular(rheo_CTL_20s_rad)
rho.circular(rheo_CTL_20s_rad)
rho.rheo_CTL_20s_rad <- rho.circular(rheo_CTL_20s_rad)
sd.circular(rheo_CTL_20s_rad)*180/pi
sd.rheo_CTL_20s_rad <- sd.circular(rheo_CTL_20s_rad)*180/pi
plot_CTL_20s <- plot.circular(rheo_CTL_20s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 20s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_20s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_20s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/Control_30s_rheo.eps")  
mean.circular(rheo_CTL_30s_rad)*180/pi
mean.rheo_CTL_30s_deg <- mean.circular(rheo_CTL_30s_rad)*180/pi
mean.circular(rheo_CTL_30s_rad)
mean.rheo_CTL_30s_rad <- mean.circular(rheo_CTL_30s_rad)
rho.circular(rheo_CTL_30s_rad)
rho.rheo_CTL_30s_rad <- rho.circular(rheo_CTL_30s_rad)
sd.circular(rheo_CTL_30s_rad)*180/pi
sd.rheo_CTL_30s_rad <- sd.circular(rheo_CTL_30s_rad)*180/pi
plot_CTL_30s <- plot.circular(rheo_CTL_30s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 30s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_30s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_30s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/Control_40s_rheo.eps")  
mean.circular(rheo_CTL_40s_rad)*180/pi
mean.rheo_CTL_40s_deg <- mean.circular(rheo_CTL_40s_rad)*180/pi
mean.circular(rheo_CTL_40s_rad)
mean.rheo_CTL_40s_rad <- mean.circular(rheo_CTL_40s_rad)
rho.circular(rheo_CTL_40s_rad)
rho.rheo_CTL_40s_rad <- rho.circular(rheo_CTL_40s_rad)
sd.circular(rheo_CTL_40s_rad)*180/pi
sd.rheo_CTL_40s_rad <- sd.circular(rheo_CTL_40s_rad)*180/pi
plot_CTL_40s <- plot.circular(rheo_CTL_40s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 40s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_40s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_40s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/Control_50s_rheo.eps")  
mean.circular(rheo_CTL_50s_rad)*180/pi
mean.rheo_CTL_50s_deg <- mean.circular(rheo_CTL_50s_rad)*180/pi
mean.circular(rheo_CTL_50s_rad)
mean.rheo_CTL_50s_rad <- mean.circular(rheo_CTL_50s_rad)
rho.circular(rheo_CTL_50s_rad)
rho.rheo_CTL_50s_rad <- rho.circular(rheo_CTL_50s_rad)
sd.circular(rheo_CTL_50s_rad)*180/pi
sd.rheo_CTL_50s_rad <- sd.circular(rheo_CTL_50s_rad)*180/pi
plot_CTL_50s <- plot.circular(rheo_CTL_50s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 50s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_50s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_50s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/Control_60s_rheo.eps")  
mean.circular(rheo_CTL_60s_rad)*180/pi
mean.rheo_CTL_60s_deg <- mean.circular(rheo_CTL_60s_rad)*180/pi
mean.circular(rheo_CTL_60s_rad)
mean.rheo_CTL_60s_rad <- mean.circular(rheo_CTL_60s_rad)
rho.circular(rheo_CTL_60s_rad)
rho.rheo_CTL_60s_rad <- rho.circular(rheo_CTL_60s_rad)
sd.circular(rheo_CTL_60s_rad)*180/pi
sd.rheo_CTL_60s_rad <- sd.circular(rheo_CTL_60s_rad)*180/pi
plot_CTL_60s <- plot.circular(rheo_CTL_60s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 60s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_60s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_60s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/Control_Post_rheo.eps")  
mean.circular(rheo_CTL_post_rad)*180/pi
mean.rheo_CTL_post_deg <- mean.circular(rheo_CTL_post_rad)*180/pi
mean.circular(rheo_CTL_post_rad)
mean.rheo_CTL_post_rad <- mean.circular(rheo_CTL_post_rad)
rho.circular(rheo_CTL_post_rad)
rho.rheo_CTL_post_rad <- rho.circular(rheo_CTL_post_rad)
sd.circular(rheo_CTL_post_rad)*180/pi
sd.rheo_CTL_post_rad <- sd.circular(rheo_CTL_post_rad)*180/pi
plot_CTL_post <- plot.circular(rheo_CTL_post_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                               start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                               tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                               xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                               template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: Post-Stimulus", 
                               sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_post_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_post_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd





setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/shake_Pre_rheo.eps")
mean.circular(rheo_shake_ps_rad)*180/pi
mean.rheo_shake_ps_deg <- mean.circular(rheo_shake_ps_rad)*180/pi
mean.circular(rheo_shake_ps_rad)
mean.rheo_shake_ps_rad <- mean.circular(rheo_shake_ps_rad)
rho.circular(rheo_shake_ps_rad)
rho.rheo_shake_ps_rad <- rho.circular(rheo_shake_ps_rad)
sd.circular(rheo_shake_ps_rad)*180/pi
sd.rheo_shake_ps_rad <- sd.circular(rheo_shake_ps_rad)*180/pi
plot_shake_ps <- plot.circular(rheo_shake_ps_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "shake Fish: Pre-Stimulus", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_ps_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_ps_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/shake_10s_rheo.eps")  
mean.circular(rheo_shake_10s_rad)*180/pi
mean.rheo_shake_10s_deg <- mean.circular(rheo_shake_10s_rad)*180/pi
mean.circular(rheo_shake_10s_rad)
mean.rheo_shake_10s_rad <- mean.circular(rheo_shake_10s_rad)
rho.circular(rheo_shake_10s_rad)
rho.rheo_shake_10s_rad <- rho.circular(rheo_shake_10s_rad)
sd.circular(rheo_shake_10s_rad)*180/pi
sd.rheo_shake_10s_rad <- sd.circular(rheo_shake_10s_rad)*180/pi
plot_shake_10s <- plot.circular(rheo_shake_10s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "shake Fish: 10s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_10s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_10s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/shake_20s_rheo.eps")  
mean.circular(rheo_shake_20s_rad)*180/pi
mean.rheo_shake_20s_deg <- mean.circular(rheo_shake_20s_rad)*180/pi
mean.circular(rheo_shake_20s_rad)
mean.rheo_shake_20s_rad <- mean.circular(rheo_shake_20s_rad)
rho.circular(rheo_shake_20s_rad)
rho.rheo_shake_20s_rad <- rho.circular(rheo_shake_20s_rad)
sd.circular(rheo_shake_20s_rad)*180/pi
sd.rheo_shake_20s_rad <- sd.circular(rheo_shake_20s_rad)*180/pi
plot_shake_20s <- plot.circular(rheo_shake_20s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "shake Fish: 20s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_20s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_20s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/shake_30s_rheo.eps")  
mean.circular(rheo_shake_30s_rad)*180/pi
mean.rheo_shake_30s_deg <- mean.circular(rheo_shake_30s_rad)*180/pi
mean.circular(rheo_shake_30s_rad)
mean.rheo_shake_30s_rad <- mean.circular(rheo_shake_30s_rad)
rho.circular(rheo_shake_30s_rad)
rho.rheo_shake_30s_rad <- rho.circular(rheo_shake_30s_rad)
sd.circular(rheo_shake_30s_rad)*180/pi
sd.rheo_shake_30s_rad <- sd.circular(rheo_shake_30s_rad)*180/pi
plot_shake_30s <- plot.circular(rheo_shake_30s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "shake Fish: 30s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_30s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_30s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/shake_40s_rheo.eps")  
mean.circular(rheo_shake_40s_rad)*180/pi
mean.rheo_shake_40s_deg <- mean.circular(rheo_shake_40s_rad)*180/pi
mean.circular(rheo_shake_40s_rad)
mean.rheo_shake_40s_rad <- mean.circular(rheo_shake_40s_rad)
rho.circular(rheo_shake_40s_rad)
rho.rheo_shake_40s_rad <- rho.circular(rheo_shake_40s_rad)
sd.circular(rheo_shake_40s_rad)*180/pi
sd.rheo_shake_40s_rad <- sd.circular(rheo_shake_40s_rad)*180/pi
plot_shake_40s <- plot.circular(rheo_shake_40s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "shake Fish: 40s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_40s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_40s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/shake_50s_rheo.eps")  
mean.circular(rheo_shake_50s_rad)*180/pi
mean.rheo_shake_50s_deg <- mean.circular(rheo_shake_50s_rad)*180/pi
mean.circular(rheo_shake_50s_rad)
mean.rheo_shake_50s_rad <- mean.circular(rheo_shake_50s_rad)
rho.circular(rheo_shake_50s_rad)
rho.rheo_shake_50s_rad <- rho.circular(rheo_shake_50s_rad)
sd.circular(rheo_shake_50s_rad)*180/pi
sd.rheo_shake_50s_rad <- sd.circular(rheo_shake_50s_rad)*180/pi
plot_shake_50s <- plot.circular(rheo_shake_50s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "shake Fish: 50s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_50s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_50s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/shake_60s_rheo.eps")  
mean.circular(rheo_shake_60s_rad)*180/pi
mean.rheo_shake_60s_deg <- mean.circular(rheo_shake_60s_rad)*180/pi
mean.circular(rheo_shake_60s_rad)
mean.rheo_shake_60s_rad <- mean.circular(rheo_shake_60s_rad)
rho.circular(rheo_shake_60s_rad)
rho.rheo_shake_60s_rad <- rho.circular(rheo_shake_60s_rad)
sd.circular(rheo_shake_60s_rad)*180/pi
sd.rheo_shake_60s_rad <- sd.circular(rheo_shake_60s_rad)*180/pi
plot_shake_60s <- plot.circular(rheo_shake_60s_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "shake Fish: 60s-Stimulus", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_60s_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_60s_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


setEPS()                                             # Set postscript arguments
postscript("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/shake_Post_rheo.eps")  
mean.circular(rheo_shake_post_rad)*180/pi
mean.rheo_shake_post_deg <- mean.circular(rheo_shake_post_rad)*180/pi
mean.circular(rheo_shake_post_rad)
mean.rheo_shake_post_rad <- mean.circular(rheo_shake_post_rad)
rho.circular(rheo_shake_post_rad)
rho.rheo_shake_post_rad <- rho.circular(rheo_shake_post_rad)
sd.circular(rheo_shake_post_rad)*180/pi
sd.rheo_shake_post_rad <- sd.circular(rheo_shake_post_rad)*180/pi
plot_shake_post <- plot.circular(rheo_shake_post_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                               start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                               tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                               xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                               template = NULL, zero = pi/2, rotation = "clock", main = "shake Fish: Post-Stimulus", 
                               sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_post_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_post_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd











####---SHAKE-RECOVERY-GRAPHS---####
#Control fish angular stats: Stimulus dataframes
#Pre-stimulus -NO -REC



#Control fish angular stats: Stimulus dataframes
#Pre-stimulus
# 0hr-rec
setEPS()                                             # Set postscript arguments
postscript("CTL_pre_0hr.eps")  
mean.circular(rheo_CTL_ps0_rad)*180/pi
mean.rheo_CTL_ps0_deg <- mean.circular(rheo_CTL_ps0_rad)*180/pi
mean.circular(rheo_CTL_ps0_rad)
mean.rheo_CTL_ps0_rad <- mean.circular(rheo_CTL_ps0_rad)
rho.circular(rheo_CTL_ps0_rad)
rho.rheo_CTL_ps0_rad <- rho.circular(rheo_CTL_ps0_rad)
sd.circular(rheo_CTL_ps0_rad)*180/pi
sd.rheo_CTL_ps0_rad <- sd.circular(rheo_CTL_ps0_rad)*180/pi
plot_CTL_ps <- plot.circular(rheo_CTL_ps0_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: Pre-Stimulus, 0hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_ps0_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_ps0_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd



# 2hr-rec
setEPS()                                             # Set postscript arguments
postscript("CTL_pre_2hr.eps")  
mean.circular(rheo_CTL_ps2_rad)*180/pi
mean.rheo_CTL_ps2_deg <- mean.circular(rheo_CTL_ps2_rad)*180/pi
mean.circular(rheo_CTL_ps2_rad)
mean.rheo_CTL_ps2_rad <- mean.circular(rheo_CTL_ps2_rad)
rho.circular(rheo_CTL_ps2_rad)
rho.rheo_CTL_ps2_rad <- rho.circular(rheo_CTL_ps2_rad)
sd.circular(rheo_CTL_ps2_rad)*180/pi
sd.rheo_CTL_ps2_rad <- sd.circular(rheo_CTL_ps2_rad)*180/pi
plot_CTL_ps <- plot.circular(rheo_CTL_ps2_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: Pre-Stimulus, 2hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_ps2_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_ps2_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 4hr-rec
setEPS()                                             # Set postscript arguments
postscript("CTL_pre_4hr.eps")  
mean.circular(rheo_CTL_ps4_rad)*180/pi
mean.rheo_CTL_ps4_deg <- mean.circular(rheo_CTL_ps4_rad)*180/pi
mean.circular(rheo_CTL_ps4_rad)
mean.rheo_CTL_ps4_rad <- mean.circular(rheo_CTL_ps4_rad)
rho.circular(rheo_CTL_ps4_rad)
rho.rheo_CTL_ps4_rad <- rho.circular(rheo_CTL_ps4_rad)
sd.circular(rheo_CTL_ps4_rad)*180/pi
sd.rheo_CTL_ps4_rad <- sd.circular(rheo_CTL_ps4_rad)*180/pi
plot_CTL_ps <- plot.circular(rheo_CTL_ps4_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: Pre-Stimulus, 4hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_ps4_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_ps4_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 8hr-rec
setEPS()                                             # Set postscript arguments
postscript("CTL_pre_8hr.eps")  
mean.circular(rheo_CTL_ps8_rad)*180/pi
mean.rheo_CTL_ps8_deg <- mean.circular(rheo_CTL_ps8_rad)*180/pi
mean.circular(rheo_CTL_ps8_rad)
mean.rheo_CTL_ps8_rad <- mean.circular(rheo_CTL_ps8_rad)
rho.circular(rheo_CTL_ps8_rad)
rho.rheo_CTL_ps8_rad <- rho.circular(rheo_CTL_ps8_rad)
sd.circular(rheo_CTL_ps8_rad)*180/pi
sd.rheo_CTL_ps8_rad <- sd.circular(rheo_CTL_ps8_rad)*180/pi
plot_CTL_ps <- plot.circular(rheo_CTL_ps8_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: Pre-Stimulus, 8hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_ps8_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_ps8_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 48hr-rec
setEPS()                                             # Set postscript arguments
postscript("CTL_pre_48hr.eps") 
mean.circular(rheo_CTL_ps48_rad)*180/pi
mean.rheo_CTL_ps48_deg <- mean.circular(rheo_CTL_ps48_rad)*180/pi
mean.circular(rheo_CTL_ps48_rad)
mean.rheo_CTL_ps48_rad <- mean.circular(rheo_CTL_ps48_rad)
rho.circular(rheo_CTL_ps48_rad)
rho.rheo_CTL_ps48_rad <- rho.circular(rheo_CTL_ps48_rad)
sd.circular(rheo_CTL_ps48_rad)*180/pi
sd.rheo_CTL_ps48_rad <- sd.circular(rheo_CTL_ps48_rad)*180/pi
plot_CTL_ps <- plot.circular(rheo_CTL_ps48_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: Pre-Stimulus, 48hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_ps48_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_ps48_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd



##Stimulus-10s
# 0hr-rec
setEPS()                                             # Set postscript arguments
postscript("CTL_10s_0hr.eps")  
mean.circular(rheo_CTL_10s0_rad)*180/pi
mean.rheo_CTL_10s0_deg <- mean.circular(rheo_CTL_10s0_rad)*180/pi
mean.circular(rheo_CTL_10s0_rad)
mean.rheo_CTL_10s0_rad <- mean.circular(rheo_CTL_10s0_rad)
rho.circular(rheo_CTL_10s0_rad)
rho.rheo_CTL_10s0_rad <- rho.circular(rheo_CTL_10s0_rad)
sd.circular(rheo_CTL_10s0_rad)*180/pi
sd.rheo_CTL_10s0_rad <- sd.circular(rheo_CTL_10s0_rad)*180/pi
plot_CTL_10s <- plot.circular(rheo_CTL_10s0_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 10s-Stimulus, 0hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_10s0_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_10s0_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd



# 2hr-rec
setEPS()                                             # Set postscript arguments
postscript("CTL_10s_2hr.eps")  
mean.circular(rheo_CTL_10s2_rad)*180/pi
mean.rheo_CTL_10s2_deg <- mean.circular(rheo_CTL_10s2_rad)*180/pi
mean.circular(rheo_CTL_10s2_rad)
mean.rheo_CTL_10s2_rad <- mean.circular(rheo_CTL_10s2_rad)
rho.circular(rheo_CTL_10s2_rad)
rho.rheo_CTL_10s2_rad <- rho.circular(rheo_CTL_10s2_rad)
sd.circular(rheo_CTL_10s2_rad)*180/pi
sd.rheo_CTL_10s2_rad <- sd.circular(rheo_CTL_10s2_rad)*180/pi
plot_CTL_10s <- plot.circular(rheo_CTL_10s2_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 10s-Stimulus, 2hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_10s2_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_10s2_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 4hr-rec
setEPS()                                             # Set postscript arguments
postscript("CTL_10s_4hr.eps")  
mean.circular(rheo_CTL_10s4_rad)*180/pi
mean.rheo_CTL_10s4_deg <- mean.circular(rheo_CTL_10s4_rad)*180/pi
mean.circular(rheo_CTL_10s4_rad)
mean.rheo_CTL_10s4_rad <- mean.circular(rheo_CTL_10s4_rad)
rho.circular(rheo_CTL_10s4_rad)
rho.rheo_CTL_10s4_rad <- rho.circular(rheo_CTL_10s4_rad)
sd.circular(rheo_CTL_10s4_rad)*180/pi
sd.rheo_CTL_10s4_rad <- sd.circular(rheo_CTL_10s4_rad)*180/pi
plot_CTL_10s <- plot.circular(rheo_CTL_10s4_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 10s-Stimulus, 4hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_10s4_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_10s4_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 8hr-rec
setEPS()                                             # Set postscript arguments
postscript("CTL_10s_8hr.eps")  
mean.circular(rheo_CTL_10s8_rad)*180/pi
mean.rheo_CTL_10s8_deg <- mean.circular(rheo_CTL_10s8_rad)*180/pi
mean.circular(rheo_CTL_10s8_rad)
mean.rheo_CTL_10s8_rad <- mean.circular(rheo_CTL_10s8_rad)
rho.circular(rheo_CTL_10s8_rad)
rho.rheo_CTL_10s8_rad <- rho.circular(rheo_CTL_10s8_rad)
sd.circular(rheo_CTL_10s8_rad)*180/pi
sd.rheo_CTL_10s8_rad <- sd.circular(rheo_CTL_10s8_rad)*180/pi
plot_CTL_10s <- plot.circular(rheo_CTL_10s8_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 10s-Stimulus, 8hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_10s8_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_10s8_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 48hr-rec
setEPS()                                             # Set postscript arguments
postscript("CTL_10s_48hr.eps") 
mean.circular(rheo_CTL_10s48_rad)*180/pi
mean.rheo_CTL_10s48_deg <- mean.circular(rheo_CTL_10s48_rad)*180/pi
mean.circular(rheo_CTL_10s48_rad)
mean.rheo_CTL_10s48_rad <- mean.circular(rheo_CTL_10s48_rad)
rho.circular(rheo_CTL_10s48_rad)
rho.rheo_CTL_10s48_rad <- rho.circular(rheo_CTL_10s48_rad)
sd.circular(rheo_CTL_10s48_rad)*180/pi
sd.rheo_CTL_10s48_rad <- sd.circular(rheo_CTL_10s48_rad)*180/pi
plot_CTL_10s <- plot.circular(rheo_CTL_10s48_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 10s-Stimulus, 48hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_10s48_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_10s48_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd



#Stimulus-20s
# 0hr-rec
setEPS()                                             # Set postscript arguments
postscript("CTL_20s_0hr.eps")  
mean.circular(rheo_CTL_20s0_rad)*180/pi
mean.rheo_CTL_20s0_deg <- mean.circular(rheo_CTL_20s0_rad)*180/pi
mean.circular(rheo_CTL_20s0_rad)
mean.rheo_CTL_20s0_rad <- mean.circular(rheo_CTL_20s0_rad)
rho.circular(rheo_CTL_20s0_rad)
rho.rheo_CTL_20s0_rad <- rho.circular(rheo_CTL_20s0_rad)
sd.circular(rheo_CTL_20s0_rad)*180/pi
sd.rheo_CTL_20s0_rad <- sd.circular(rheo_CTL_20s0_rad)*180/pi
plot_CTL_20s <- plot.circular(rheo_CTL_20s0_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 20s-Stimulus, 0hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_20s0_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_20s0_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd



# 2hr-rec
setEPS()                                             # Set postscript arguments
postscript("CTL_20s_2hr.eps")  
mean.circular(rheo_CTL_20s2_rad, na.rm=TRUE)*180/pi
mean.rheo_CTL_20s2_deg <- mean.circular(rheo_CTL_20s2_rad, na.rm=TRUE)*180/pi
mean.circular(rheo_CTL_20s2_rad, na.rm=TRUE)
mean.rheo_CTL_20s2_rad <- mean.circular(rheo_CTL_20s2_rad, na.rm=TRUE)
rho.circular(rheo_CTL_20s2_rad, na.rm=TRUE)
rho.rheo_CTL_20s2_rad <- rho.circular(rheo_CTL_20s2_rad, na.rm=TRUE)
sd.circular(rheo_CTL_20s2_rad, na.rm=TRUE)*180/pi
sd.rheo_CTL_20s2_rad <- sd.circular(rheo_CTL_20s2_rad, na.rm=TRUE)*180/pi
plot_CTL_20s <- plot.circular(rheo_CTL_20s2_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 20s-Stimulus, 2hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_20s2_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_20s2_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 4hr-rec
setEPS()                                             # Set postscript arguments
postscript("CTL_20s_4hr.eps")  
mean.circular(rheo_CTL_20s4_rad)*180/pi
mean.rheo_CTL_20s4_deg <- mean.circular(rheo_CTL_20s4_rad)*180/pi
mean.circular(rheo_CTL_20s4_rad)
mean.rheo_CTL_20s4_rad <- mean.circular(rheo_CTL_20s4_rad)
rho.circular(rheo_CTL_20s4_rad)
rho.rheo_CTL_20s4_rad <- rho.circular(rheo_CTL_20s4_rad)
sd.circular(rheo_CTL_20s4_rad)*180/pi
sd.rheo_CTL_20s4_rad <- sd.circular(rheo_CTL_20s4_rad)*180/pi
plot_CTL_20s <- plot.circular(rheo_CTL_20s4_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 20s-Stimulus, 4hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_20s4_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_20s4_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 8hr-rec
setEPS()                                             # Set postscript arguments
postscript("CTL_20s_8hr.eps")  
mean.circular(rheo_CTL_20s8_rad)*180/pi
mean.rheo_CTL_20s8_deg <- mean.circular(rheo_CTL_20s8_rad)*180/pi
mean.circular(rheo_CTL_20s8_rad)
mean.rheo_CTL_20s8_rad <- mean.circular(rheo_CTL_20s8_rad)
rho.circular(rheo_CTL_20s8_rad)
rho.rheo_CTL_20s8_rad <- rho.circular(rheo_CTL_20s8_rad)
sd.circular(rheo_CTL_20s8_rad)*180/pi
sd.rheo_CTL_20s8_rad <- sd.circular(rheo_CTL_20s8_rad)*180/pi
plot_CTL_20s <- plot.circular(rheo_CTL_20s8_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 20s-Stimulus, 8hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_20s8_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_20s8_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 48hr-rec
setEPS()                                             # Set postscript arguments
postscript("CTL_20s_48hr.eps") 
mean.circular(rheo_CTL_20s48_rad)*180/pi
mean.rheo_CTL_20s48_deg <- mean.circular(rheo_CTL_20s48_rad)*180/pi
mean.circular(rheo_CTL_20s48_rad)
mean.rheo_CTL_20s48_rad <- mean.circular(rheo_CTL_20s48_rad)
rho.circular(rheo_CTL_20s48_rad)
rho.rheo_CTL_20s48_rad <- rho.circular(rheo_CTL_20s48_rad)
sd.circular(rheo_CTL_20s48_rad)*180/pi
sd.rheo_CTL_20s48_rad <- sd.circular(rheo_CTL_20s48_rad)*180/pi
plot_CTL_20s <- plot.circular(rheo_CTL_20s48_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: 20s-Stimulus, 48hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_20s48_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_20s48_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd



#Post-stimulus
# 0hr-rec
setEPS()                                             # Set postscript arguments
postscript("CTL_post_0hr.eps")  
mean.circular(rheo_CTL_post0_rad)*180/pi
mean.rheo_CTL_post0_deg <- mean.circular(rheo_CTL_post0_rad)*180/pi
mean.circular(rheo_CTL_post0_rad)
mean.rheo_CTL_post0_rad <- mean.circular(rheo_CTL_post0_rad)
rho.circular(rheo_CTL_post0_rad)
rho.rheo_CTL_post0_rad <- rho.circular(rheo_CTL_post0_rad)
sd.circular(rheo_CTL_post0_rad)*180/pi
sd.rheo_CTL_post0_rad <- sd.circular(rheo_CTL_post0_rad)*180/pi
plot_CTL_post <- plot.circular(rheo_CTL_post0_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: post-Stimulus, 0hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_post0_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_post0_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd



# 2hr-rec
setEPS()                                             # Set postscript arguments
postscript("CTL_post_2hr.eps")  
mean.circular(rheo_CTL_post2_rad)*180/pi
mean.rheo_CTL_post2_deg <- mean.circular(rheo_CTL_post2_rad)*180/pi
mean.circular(rheo_CTL_post2_rad)
mean.rheo_CTL_post2_rad <- mean.circular(rheo_CTL_post2_rad)
rho.circular(rheo_CTL_post2_rad)
rho.rheo_CTL_post2_rad <- rho.circular(rheo_CTL_post2_rad)
sd.circular(rheo_CTL_post2_rad)*180/pi
sd.rheo_CTL_post2_rad <- sd.circular(rheo_CTL_post2_rad)*180/pi
plot_CTL_post <- plot.circular(rheo_CTL_post2_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: post-Stimulus, 2hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_post2_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_post2_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 4hr-rec
setEPS()                                             # Set postscript arguments
postscript("CTL_post_4hr.eps")  
mean.circular(rheo_CTL_post4_rad)*180/pi
mean.rheo_CTL_post4_deg <- mean.circular(rheo_CTL_post4_rad)*180/pi
mean.circular(rheo_CTL_post4_rad)
mean.rheo_CTL_post4_rad <- mean.circular(rheo_CTL_post4_rad)
rho.circular(rheo_CTL_post4_rad)
rho.rheo_CTL_post4_rad <- rho.circular(rheo_CTL_post4_rad)
sd.circular(rheo_CTL_post4_rad)*180/pi
sd.rheo_CTL_post4_rad <- sd.circular(rheo_CTL_post4_rad)*180/pi
plot_CTL_post <- plot.circular(rheo_CTL_post4_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: post-Stimulus, 4hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_post4_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_post4_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 8hr-rec
setEPS()                                             # Set postscript arguments
postscript("CTL_post_8hr.eps")  
mean.circular(rheo_CTL_post8_rad)*180/pi
mean.rheo_CTL_post8_deg <- mean.circular(rheo_CTL_post8_rad)*180/pi
mean.circular(rheo_CTL_post8_rad)
mean.rheo_CTL_post8_rad <- mean.circular(rheo_CTL_post8_rad)
rho.circular(rheo_CTL_post8_rad)
rho.rheo_CTL_post8_rad <- rho.circular(rheo_CTL_post8_rad)
sd.circular(rheo_CTL_post8_rad)*180/pi
sd.rheo_CTL_post8_rad <- sd.circular(rheo_CTL_post8_rad)*180/pi
plot_CTL_post <- plot.circular(rheo_CTL_post8_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: post-Stimulus, 8hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_post8_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_post8_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 48hr-rec
setEPS()                                             # Set postscript arguments
postscript("CTL_post_48hr.eps") 
mean.circular(rheo_CTL_post48_rad)*180/pi
mean.rheo_CTL_post48_deg <- mean.circular(rheo_CTL_post48_rad)*180/pi
mean.circular(rheo_CTL_post48_rad)
mean.rheo_CTL_post48_rad <- mean.circular(rheo_CTL_post48_rad)
rho.circular(rheo_CTL_post48_rad)
rho.rheo_CTL_post48_rad <- rho.circular(rheo_CTL_post48_rad)
sd.circular(rheo_CTL_post48_rad)*180/pi
sd.rheo_CTL_post48_rad <- sd.circular(rheo_CTL_post48_rad)*180/pi
plot_CTL_post <- plot.circular(rheo_CTL_post48_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Control Fish: post-Stimulus, 48hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_CTL_post48_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_CTL_post48_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


###SHAKE###
#Shake fish angular stats: Stimulus dataframes
#Pre-stimulus
# 0hr-rec
setEPS()                                             # Set postscript arguments
postscript("shake_pre_0hr.eps")  
mean.circular(rheo_shake_ps0_rad)*180/pi
mean.rheo_shake_ps0_deg <- mean.circular(rheo_shake_ps0_rad)*180/pi
mean.circular(rheo_shake_ps0_rad)
mean.rheo_shake_ps0_rad <- mean.circular(rheo_shake_ps0_rad)
rho.circular(rheo_shake_ps0_rad)
rho.rheo_shake_ps0_rad <- rho.circular(rheo_shake_ps0_rad)
sd.circular(rheo_shake_ps0_rad)*180/pi
sd.rheo_shake_ps0_rad <- sd.circular(rheo_shake_ps0_rad)*180/pi
plot_shake_ps <- plot.circular(rheo_shake_ps0_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Shake Fish: Pre-Stimulus, 0hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_ps0_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_ps0_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd



# 2hr-rec
setEPS()                                             # Set postscript arguments
postscript("shake_pre_2hr.eps")  
mean.circular(rheo_shake_ps2_rad)*180/pi
mean.rheo_shake_ps2_deg <- mean.circular(rheo_shake_ps2_rad)*180/pi
mean.circular(rheo_shake_ps2_rad)
mean.rheo_shake_ps2_rad <- mean.circular(rheo_shake_ps2_rad)
rho.circular(rheo_shake_ps2_rad)
rho.rheo_shake_ps2_rad <- rho.circular(rheo_shake_ps2_rad)
sd.circular(rheo_shake_ps2_rad)*180/pi
sd.rheo_shake_ps2_rad <- sd.circular(rheo_shake_ps2_rad)*180/pi
plot_shake_ps <- plot.circular(rheo_shake_ps2_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Shake Fish: Pre-Stimulus, 2hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_ps2_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_ps2_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 4hr-rec
setEPS()                                             # Set postscript arguments
postscript("shake_pre_4hr.eps")  
mean.circular(rheo_shake_ps4_rad)*180/pi
mean.rheo_shake_ps4_deg <- mean.circular(rheo_shake_ps4_rad)*180/pi
mean.circular(rheo_shake_ps4_rad)
mean.rheo_shake_ps4_rad <- mean.circular(rheo_shake_ps4_rad)
rho.circular(rheo_shake_ps4_rad)
rho.rheo_shake_ps4_rad <- rho.circular(rheo_shake_ps4_rad)
sd.circular(rheo_shake_ps4_rad)*180/pi
sd.rheo_shake_ps4_rad <- sd.circular(rheo_shake_ps4_rad)*180/pi
plot_shake_ps <- plot.circular(rheo_shake_ps4_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Shake Fish: Pre-Stimulus, 4hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_ps4_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_ps4_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 8hr-rec
setEPS()                                             # Set postscript arguments
postscript("shake_pre_8hr.eps")  
mean.circular(rheo_shake_ps8_rad)*180/pi
mean.rheo_shake_ps8_deg <- mean.circular(rheo_shake_ps8_rad)*180/pi
mean.circular(rheo_shake_ps8_rad)
mean.rheo_shake_ps8_rad <- mean.circular(rheo_shake_ps8_rad)
rho.circular(rheo_shake_ps8_rad)
rho.rheo_shake_ps8_rad <- rho.circular(rheo_shake_ps8_rad)
sd.circular(rheo_shake_ps8_rad)*180/pi
sd.rheo_shake_ps8_rad <- sd.circular(rheo_shake_ps8_rad)*180/pi
plot_shake_ps <- plot.circular(rheo_shake_ps8_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Shake Fish: Pre-Stimulus, 8hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_ps8_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_ps8_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 48hr-rec
setEPS()                                             # Set postscript arguments
postscript("shake_pre_48hr.eps") 
mean.circular(rheo_shake_ps48_rad)*180/pi
mean.rheo_shake_ps48_deg <- mean.circular(rheo_shake_ps48_rad)*180/pi
mean.circular(rheo_shake_ps48_rad)
mean.rheo_shake_ps48_rad <- mean.circular(rheo_shake_ps48_rad)
rho.circular(rheo_shake_ps48_rad)
rho.rheo_shake_ps48_rad <- rho.circular(rheo_shake_ps48_rad)
sd.circular(rheo_shake_ps48_rad)*180/pi
sd.rheo_shake_ps48_rad <- sd.circular(rheo_shake_ps48_rad)*180/pi
plot_shake_ps <- plot.circular(rheo_shake_ps48_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                             start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                             tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                             xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                             template = NULL, zero = pi/2, rotation = "clock", main = "Shake Fish: Pre-Stimulus, 48hr", 
                             sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_ps48_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_ps48_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd



##Stimulus-10s
# 0hr-rec
setEPS()                                             # Set postscript arguments
postscript("shake_10s_0hr.eps")  
mean.circular(rheo_shake_10s0_rad)*180/pi
mean.rheo_shake_10s0_deg <- mean.circular(rheo_shake_10s0_rad)*180/pi
mean.circular(rheo_shake_10s0_rad)
mean.rheo_shake_10s0_rad <- mean.circular(rheo_shake_10s0_rad)
rho.circular(rheo_shake_10s0_rad)
rho.rheo_shake_10s0_rad <- rho.circular(rheo_shake_10s0_rad)
sd.circular(rheo_shake_10s0_rad)*180/pi
sd.rheo_shake_10s0_rad <- sd.circular(rheo_shake_10s0_rad)*180/pi
plot_shake_10s <- plot.circular(rheo_shake_10s0_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Shake Fish: 10s-Stimulus, 0hr", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_10s0_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_10s0_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd



# 2hr-rec
setEPS()                                             # Set postscript arguments
postscript("shake_10s_2hr.eps")  
mean.circular(rheo_shake_10s2_rad)*180/pi
mean.rheo_shake_10s2_deg <- mean.circular(rheo_shake_10s2_rad)*180/pi
mean.circular(rheo_shake_10s2_rad)
mean.rheo_shake_10s2_rad <- mean.circular(rheo_shake_10s2_rad)
rho.circular(rheo_shake_10s2_rad)
rho.rheo_shake_10s2_rad <- rho.circular(rheo_shake_10s2_rad)
sd.circular(rheo_shake_10s2_rad)*180/pi
sd.rheo_shake_10s2_rad <- sd.circular(rheo_shake_10s2_rad)*180/pi
plot_shake_10s <- plot.circular(rheo_shake_10s2_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Shake Fish: 10s-Stimulus, 2hr", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_10s2_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_10s2_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 4hr-rec
setEPS()                                             # Set postscript arguments
postscript("shake_10s_4hr.eps")  
mean.circular(rheo_shake_10s4_rad)*180/pi
mean.rheo_shake_10s4_deg <- mean.circular(rheo_shake_10s4_rad)*180/pi
mean.circular(rheo_shake_10s4_rad)
mean.rheo_shake_10s4_rad <- mean.circular(rheo_shake_10s4_rad)
rho.circular(rheo_shake_10s4_rad)
rho.rheo_shake_10s4_rad <- rho.circular(rheo_shake_10s4_rad)
sd.circular(rheo_shake_10s4_rad)*180/pi
sd.rheo_shake_10s4_rad <- sd.circular(rheo_shake_10s4_rad)*180/pi
plot_shake_10s <- plot.circular(rheo_shake_10s4_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Shake Fish: 10s-Stimulus, 4hr", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_10s4_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_10s4_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 8hr-rec
setEPS()                                             # Set postscript arguments
postscript("shake_10s_8hr.eps")  
mean.circular(rheo_shake_10s8_rad)*180/pi
mean.rheo_shake_10s8_deg <- mean.circular(rheo_shake_10s8_rad)*180/pi
mean.circular(rheo_shake_10s8_rad)
mean.rheo_shake_10s8_rad <- mean.circular(rheo_shake_10s8_rad)
rho.circular(rheo_shake_10s8_rad)
rho.rheo_shake_10s8_rad <- rho.circular(rheo_shake_10s8_rad)
sd.circular(rheo_shake_10s8_rad)*180/pi
sd.rheo_shake_10s8_rad <- sd.circular(rheo_shake_10s8_rad)*180/pi
plot_shake_10s <- plot.circular(rheo_shake_10s8_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Shake Fish: 10s-Stimulus, 8hr", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_10s8_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_10s8_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 48hr-rec
setEPS()                                             # Set postscript arguments
postscript("shake_10s_48hr.eps") 
mean.circular(rheo_shake_10s48_rad)*180/pi
mean.rheo_shake_10s48_deg <- mean.circular(rheo_shake_10s48_rad)*180/pi
mean.circular(rheo_shake_10s48_rad)
mean.rheo_shake_10s48_rad <- mean.circular(rheo_shake_10s48_rad)
rho.circular(rheo_shake_10s48_rad)
rho.rheo_shake_10s48_rad <- rho.circular(rheo_shake_10s48_rad)
sd.circular(rheo_shake_10s48_rad)*180/pi
sd.rheo_shake_10s48_rad <- sd.circular(rheo_shake_10s48_rad)*180/pi
plot_shake_10s <- plot.circular(rheo_shake_10s48_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Shake Fish: 10s-Stimulus, 48hr", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_10s48_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_10s48_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd



#Stimulus-20s
# 0hr-rec
setEPS()                                             # Set postscript arguments
postscript("shake_20s_0hr.eps")  
mean.circular(rheo_shake_20s0_rad)*180/pi
mean.rheo_shake_20s0_deg <- mean.circular(rheo_shake_20s0_rad)*180/pi
mean.circular(rheo_shake_20s0_rad)
mean.rheo_shake_20s0_rad <- mean.circular(rheo_shake_20s0_rad)
rho.circular(rheo_shake_20s0_rad)
rho.rheo_shake_20s0_rad <- rho.circular(rheo_shake_20s0_rad)
sd.circular(rheo_shake_20s0_rad)*180/pi
sd.rheo_shake_20s0_rad <- sd.circular(rheo_shake_20s0_rad)*180/pi
plot_shake_20s <- plot.circular(rheo_shake_20s0_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Shake Fish: 20s-Stimulus, 0hr", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_20s0_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_20s0_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd



# 2hr-rec
setEPS()                                             # Set postscript arguments
postscript("shake_20s_2hr.eps")  
mean.circular(rheo_shake_20s2_rad)*180/pi
mean.rheo_shake_20s2_deg <- mean.circular(rheo_shake_20s2_rad)*180/pi
mean.circular(rheo_shake_20s2_rad)
mean.rheo_shake_20s2_rad <- mean.circular(rheo_shake_20s2_rad)
rho.circular(rheo_shake_20s2_rad)
rho.rheo_shake_20s2_rad <- rho.circular(rheo_shake_20s2_rad)
sd.circular(rheo_shake_20s2_rad)*180/pi
sd.rheo_shake_20s2_rad <- sd.circular(rheo_shake_20s2_rad, na.rm=TRUE)*180/pi
plot_shake_20s <- plot.circular(rheo_shake_20s2_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Shake Fish: 20s-Stimulus, 2hr", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_20s2_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_20s2_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 4hr-rec
setEPS()                                             # Set postscript arguments
postscript("shake_20s_4hr.eps")  
mean.circular(rheo_shake_20s4_rad, na.rm=TRUE)*180/pi
mean.rheo_shake_20s4_deg <- mean.circular(rheo_shake_20s4_rad, na.rm=TRUE)*180/pi
mean.circular(rheo_shake_20s4_rad, na.rm=TRUE)
mean.rheo_shake_20s4_rad <- mean.circular(rheo_shake_20s4_rad, na.rm=TRUE)
rho.circular(rheo_shake_20s4_rad, na.rm=TRUE)
rho.rheo_shake_20s4_rad <- rho.circular(rheo_shake_20s4_rad, na.rm=TRUE)
sd.circular(rheo_shake_20s4_rad, na.rm=TRUE)*180/pi
sd.rheo_shake_20s4_rad <- sd.circular(rheo_shake_20s4_rad, na.rm=TRUE)*180/pi
plot_shake_20s <- plot.circular(rheo_shake_20s4_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Shake Fish: 20s-Stimulus, 4hr", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_20s4_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_20s4_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 8hr-rec
setEPS()                                             # Set postscript arguments
postscript("shake_20s_8hr.eps")  
mean.circular(rheo_shake_20s8_rad)*180/pi
mean.rheo_shake_20s8_deg <- mean.circular(rheo_shake_20s8_rad)*180/pi
mean.circular(rheo_shake_20s8_rad)
mean.rheo_shake_20s8_rad <- mean.circular(rheo_shake_20s8_rad)
rho.circular(rheo_shake_20s8_rad)
rho.rheo_shake_20s8_rad <- rho.circular(rheo_shake_20s8_rad)
sd.circular(rheo_shake_20s8_rad)*180/pi
sd.rheo_shake_20s8_rad <- sd.circular(rheo_shake_20s8_rad)*180/pi
plot_shake_20s <- plot.circular(rheo_shake_20s8_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Shake Fish: 20s-Stimulus, 8hr", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_20s8_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_20s8_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 48hr-rec
setEPS()                                             # Set postscript arguments
postscript("shake_20s_48hr.eps") 
mean.circular(rheo_shake_20s48_rad, na.rm=TRUE)*180/pi
mean.rheo_shake_20s48_deg <- mean.circular(rheo_shake_20s48_rad, na.rm=TRUE)*180/pi
mean.circular(rheo_shake_20s48_rad, na.rm=TRUE)
mean.rheo_shake_20s48_rad <- mean.circular(rheo_shake_20s48_rad, na.rm=TRUE)
rho.circular(rheo_shake_20s48_rad, na.rm=TRUE)
rho.rheo_shake_20s48_rad <- rho.circular(rheo_shake_20s48_rad, na.rm=TRUE)
sd.circular(rheo_shake_20s48_rad, na.rm=TRUE)*180/pi
sd.rheo_shake_20s48_rad <- sd.circular(rheo_shake_20s48_rad, na.rm=TRUE)*180/pi
plot_shake_20s <- plot.circular(rheo_shake_20s48_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                              start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                              tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                              xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                              template = NULL, zero = pi/2, rotation = "clock", main = "Shake Fish: 20s-Stimulus, 48hr", 
                              sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_20s48_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_20s48_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd



#Post-stimulus
# 0hr-rec
setEPS()                                             # Set postscript arguments
postscript("shake_post_0hr.eps")  
mean.circular(rheo_shake_post0_rad)*180/pi
mean.rheo_shake_post0_deg <- mean.circular(rheo_shake_post0_rad)*180/pi
mean.circular(rheo_shake_post0_rad)
mean.rheo_shake_post0_rad <- mean.circular(rheo_shake_post0_rad)
rho.circular(rheo_shake_post0_rad)
rho.rheo_shake_post0_rad <- rho.circular(rheo_shake_post0_rad)
sd.circular(rheo_shake_post0_rad)*180/pi
sd.rheo_shake_post0_rad <- sd.circular(rheo_shake_post0_rad)*180/pi
plot_shake_post <- plot.circular(rheo_shake_post0_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                               start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                               tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                               xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                               template = NULL, zero = pi/2, rotation = "clock", main = "Shake Fish: post-Stimulus, 0hr", 
                               sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_post0_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_post0_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd



# 2hr-rec
setEPS()                                             # Set postscript arguments
postscript("shake_post_2hr.eps")  
mean.circular(rheo_shake_post2_rad, na.rm=TRUE)*180/pi
mean.rheo_shake_post2_deg <- mean.circular(rheo_shake_post2_rad, na.rm=TRUE)*180/pi
mean.circular(rheo_shake_post2_rad, na.rm=TRUE)
mean.rheo_shake_post2_rad <- mean.circular(rheo_shake_post2_rad, na.rm=TRUE)
rho.circular(rheo_shake_post2_rad, na.rm=TRUE)
rho.rheo_shake_post2_rad <- rho.circular(rheo_shake_post2_rad, na.rm=TRUE)
sd.circular(rheo_shake_post2_rad, na.rm=TRUE)*180/pi
sd.rheo_shake_post2_rad <- sd.circular(rheo_shake_post2_rad, na.rm=TRUE)*180/pi
plot_shake_post <- plot.circular(rheo_shake_post2_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                               start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                               tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                               xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                               template = NULL, zero = pi/2, rotation = "clock", main = "Shake Fish: post-Stimulus, 2hr", 
                               sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_post2_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_post2_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 4hr-rec
setEPS()                                             # Set postscript arguments
postscript("shake_post_4hr.eps")  
mean.circular(rheo_shake_post4_rad)*180/pi
mean.rheo_shake_post4_deg <- mean.circular(rheo_shake_post4_rad)*180/pi
mean.circular(rheo_shake_post4_rad)
mean.rheo_shake_post4_rad <- mean.circular(rheo_shake_post4_rad)
rho.circular(rheo_shake_post4_rad)
rho.rheo_shake_post4_rad <- rho.circular(rheo_shake_post4_rad)
sd.circular(rheo_shake_post4_rad)*180/pi
sd.rheo_shake_post4_rad <- sd.circular(rheo_shake_post4_rad)*180/pi
plot_shake_post <- plot.circular(rheo_shake_post4_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                               start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                               tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                               xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                               template = NULL, zero = pi/2, rotation = "clock", main = "Shake Fish: post-Stimulus, 4hr", 
                               sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_post4_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_post4_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 8hr-rec
setEPS()                                             # Set postscript arguments
postscript("shake_post_8hr.eps")  
mean.circular(rheo_shake_post8_rad)*180/pi
mean.rheo_shake_post8_deg <- mean.circular(rheo_shake_post8_rad)*180/pi
mean.circular(rheo_shake_post8_rad)
mean.rheo_shake_post8_rad <- mean.circular(rheo_shake_post8_rad)
rho.circular(rheo_shake_post8_rad)
rho.rheo_shake_post8_rad <- rho.circular(rheo_shake_post8_rad)
sd.circular(rheo_shake_post8_rad)*180/pi
sd.rheo_shake_post8_rad <- sd.circular(rheo_shake_post8_rad)*180/pi
plot_shake_post <- plot.circular(rheo_shake_post8_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                               start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                               tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                               xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                               template = NULL, zero = pi/2, rotation = "clock", main = "Shake Fish: post-Stimulus, 8hr", 
                               sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_post8_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_post8_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd


# 48hr-rec
setEPS()                                             # Set postscript arguments
postscript("shake_post_48hr.eps") 
mean.circular(rheo_shake_post48_rad)*180/pi
mean.rheo_shake_post48_deg <- mean.circular(rheo_shake_post48_rad)*180/pi
mean.circular(rheo_shake_post48_rad)
mean.rheo_shake_post48_rad <- mean.circular(rheo_shake_post48_rad)
rho.circular(rheo_shake_post48_rad)
rho.rheo_shake_post48_rad <- rho.circular(rheo_shake_post48_rad)
sd.circular(rheo_shake_post48_rad)*180/pi
sd.rheo_shake_post48_rad <- sd.circular(rheo_shake_post48_rad)*180/pi
plot_shake_post <- plot.circular(rheo_shake_post48_rad, pch = 20, cex = 0.5, stack = TRUE, axes = TRUE, 
                               start.sep = 0.04, sep = 0.1, shrink = 1, bins = 180, ticks = TRUE, 
                               tcl = 0.1, tcl.text = 0.15, col = 1, tol = 1, uin = NULL, 
                               xlim = c(-0.75, 0.75), ylim = c(-0.75, 0.75), digits = 2, units = "degrees", 
                               template = NULL, zero = pi/2, rotation = "clock", main = "Shake Fish: post-Stimulus, 48hr", 
                               sub = NULL, xlab = "", ylab = "", control.circle = circle.control(lwd = 2))
arrows.circular(mean.rheo_shake_post48_rad, y = NULL, x0 = 0, y0 = 0, na.rm = FALSE, 
                shrink = rho.rheo_shake_post48_rad, plot.info = NULL, zero = pi/2, rotation = "clock", length = .15, angle = 30, lwd = 1.5)
dev.off() #expot eps file to setwd




####---CIRCULAR STATS CUSO4 - NEOMYCIN ---####
# use radian data

# library(circular)
# 
# #full data set, not curated as suggested above
# #V-test specifies angle = 0°
# rayleigh.test(rheo_CTL$Angle_rad)
# rayleigh.test(rheo_CTL$Angle_rad, mu = circular(0))
# rayleigh.test(rheo_CuSO4$Angle_rad)
# rayleigh.test(rheo_CuSO4$Angle_rad, mu = circular(0))
# rayleigh.test(rheo_Neo$Angle_rad)
# rayleigh.test(rheo_Neo$Angle_rad, mu = circular(0))
# 
# #pre-stim, post-stim are uniform and not distributed around 0°, unlike threatments
# rayleigh.test(rheo_CTL_ps_rad)
# rayleigh.test(rheo_CTL_ps_rad, mu = circular(0))
# rayleigh.test(rheo_CTL_10s_rad)
# rayleigh.test(rheo_CTL_10s_rad, mu = circular(0))
# rayleigh.test(rheo_CTL_20s_rad)
# rayleigh.test(rheo_CTL_20s_rad, mu = circular(0))
# rayleigh.test(rheo_CTL_post_rad)
# rayleigh.test(rheo_CTL_post_rad, mu = circular(0))
# 
# #pre-stim, are uniform and not distributed around 0°, unlike threatments. post-stim are unimodal but not clustered aroun 0°
# rayleigh.test(rheo_CuSO4_ps_rad)
# rayleigh.test(rheo_CuSO4_ps_rad, mu = circular(0))
# rayleigh.test(rheo_CuSO4_10s_rad)
# rayleigh.test(rheo_CuSO4_10s_rad, mu = circular(0))
# rayleigh.test(rheo_CuSO4_20s_rad)
# rayleigh.test(rheo_CuSO4_20s_rad, mu = circular(0))
# rayleigh.test(rheo_CuSO4_post_rad)
# rayleigh.test(rheo_CuSO4_post_rad, mu = circular(0))

# #pre-stim, post-stim are uniform and not distributed around 0°, unlike threatments
# rayleigh.test(rheo_Neo_ps_rad)
# rayleigh.test(rheo_Neo_ps_rad, mu = circular(0))
# rayleigh.test(rheo_Neo_10s_rad)
# rayleigh.test(rheo_Neo_10s_rad, mu = circular(0))
# rayleigh.test(rheo_Neo_20s_rad)
# rayleigh.test(rheo_Neo_20s_rad, mu = circular(0))
# rayleigh.test(rheo_Neo_post_rad)
# rayleigh.test(rheo_Neo_post_rad, mu = circular(0))
# 
# 
# #Fitak & Jonsen 2017 - compares 10 models (DeltaAIC) and gives best fit: NOTE these tests are done on full data set, not curated. 
# # Look at AIC wieghts for probabilities in model fit. Values fluctuate each time it is run (convergence = 0), so is there a burn-in time, run for 1000 iterations? No change with 1Billion iterations
# 
# #uniform M1 (35.2%)
# circ_mle(rheo_CTL_ps_rad)
# 
# #modified unimodal M2C (99%)
# circ_mle(rheo_CTL_10s_rad)
# 
# #axial bimodal M4B (100%) 
# circ_mle(rheo_CTL_20s_rad)
# 
# #uniform M5B (27.7%)
# circ_mle(rheo_CTL_post_rad)
# 
# 
# 
# #homogeneous bimodal M5A (23.1%)
# circ_mle(rheo_CuSO4_ps_rad)
# 
# #modified unimodal M2C (54.1)
# circ_mle(rheo_CuSO4_10s_rad)
# 
# #modified unimodal M2C (99.2%)
# circ_mle(rheo_CuSO4_20s_rad)
# 
# #bimodal M2A (35.1%)
# circ_mle(rheo_CuSO4_post_rad)
# 
# 
# 
# #uniform M1 (28.2%)
# circ_mle(rheo_Neo_ps_rad)
# 
# #modified unimodal M2C (98.4%)
# circ_mle(rheo_Neo_10s_rad)
# 
# #modified unimodal M2C (80%) 
# circ_mle(rheo_Neo_20s_rad)
# 
# #symmetric uninmodal M2B (58.8%)
# circ_mle(rheo_Neo_post_rad)




# 
# watson.wheeler.test(list(rheo_CTL_10s$Angle_rad,rheo_CTL_20s$Angle_rad))
# # data:  1 and 2
# # W = 1.49, df = 2, p-value = 0.4747
# 
# watson.wheeler.test(list(rheo_CTL_10s$Angle_rad,rheo_CuSO4_10s$Angle_rad))
# # data:  1 and 2
# # W = 21.128, df = 2, p-value = 2.583e-05
# 
# watson.wheeler.test(list(rheo_CTL_10s$Angle_rad,rheo_CuSO4_20s$Angle_rad))
# # data:  1 and 2
# # W = 4.2573, df = 2, p-value = 0.119
# 
# watson.wheeler.test(list(rheo_CTL_10s$Angle_rad,rheo_Neo_10s$Angle_rad))
# # data:  1 and 2
# # W = 6.6655, df = 2, p-value = 0.0357
# 
# watson.wheeler.test(list(rheo_CTL_10s$Angle_rad,rheo_Neo_20s$Angle_rad))
# # Test Statistic: 0.1251
# # P-value > 0.10
# 
# 
# 
# watson.wheeler.test(list(rheo_CTL_20s$Angle_rad,rheo_CuSO4_10s$Angle_rad))
# # data:  1 and 2
# # W = 19.122, df = 2, p-value = 7.041e-05
# 
# watson.wheeler.test(list(rheo_CTL_20s$Angle_rad,rheo_CuSO4_20s$Angle_rad))
# # data:  1 and 2
# # W = 5.9597, df = 2, p-value = 0.0508 
# 
# watson.wheeler.test(list(rheo_CTL_20s$Angle_rad,rheo_Neo_10s$Angle_rad))
# # data:  1 and 2
# # W = 5.2022, df = 2, p-value = 0.07419
# 
# watson.wheeler.test(list(rheo_CTL_20s$Angle_rad,rheo_Neo_20s$Angle_rad))
# # data:  1 and 2
# # W = 1.5961, df = 2, p-value = 0.4502
# 
# 
# 
# 
# watson.wheeler.test(list(rheo_CuSO4_10s$Angle_rad,rheo_CuSO4_20s$Angle_rad))
# # data:  1 and 2
# # W = 6.0014, df = 2, p-value = 0.04975
# 
# watson.wheeler.test(list(rheo_CuSO4_10s$Angle_rad,rheo_Neo_10s$Angle_rad))
# # data:  1 and 2
# # W = 5.6737, df = 2, p-value = 0.05861 
# 
# watson.wheeler.test(list(rheo_CuSO4_10s$Angle_rad,rheo_Neo_20s$Angle_rad))
# # data:  1 and 2
# # W = 9.4268, df = 2, p-value = 0.008974
# 
# 
# 
# watson.wheeler.test(list(rheo_CuSO4_20s$Angle_rad,rheo_Neo_10s$Angle_rad))
# # data:  1 and 2
# # W = 2.3357, df = 2, p-value = 0.311 
# 
# watson.wheeler.test(list(rheo_CuSO4_20s$Angle_rad,rheo_Neo_20s$Angle_rad))
# # data:  1 and 2
# # W = 3.3, df = 2, p-value = 0.1921
# 
# 
# 
# 
# watson.wheeler.test(list(rheo_Neo_10s$Angle_rad,rheo_Neo_20s$Angle_rad))
# # data:  1 and 2
# # W = 0.61196, df = 2, p-value = 0.7364


####---CIRCULAR STATS BAPTA---####
# use radian data

library(circular)


#full data set, not curated as suggested above
#V-test specifies angle = 0°
rayleigh.test(rheo_CTL$Angle_rad)
rayleigh.test(rheo_CTL$Angle_rad, mu = circular(0))
rayleigh.test(rheo_shake$Angle_rad)
rayleigh.test(rheo_shake$Angle_rad, mu = circular(0))


#pre-stim, post-stim are uniform and not distributed around 0°, unlike 10s & 20s stimulus threatments
rayleigh.test(rheo_CTL_ps_rad)
rayleigh.test(rheo_CTL_ps_rad, mu = circular(0))
rayleigh.test(rheo_CTL_10s_rad)
rayleigh.test(rheo_CTL_10s_rad, mu = circular(0))
rayleigh.test(rheo_CTL_20s_rad)
rayleigh.test(rheo_CTL_20s_rad, mu = circular(0))
rayleigh.test(rheo_CTL_post_rad)
rayleigh.test(rheo_CTL_post_rad, mu = circular(0))

#pre-stim, are uniform and not distributed around 0°, unlike threatments. post-stim are unimodal but not clustered aroun 0°
rayleigh.test(rheo_shake_ps_rad)
rayleigh.test(rheo_shake_ps_rad, mu = circular(0))
rayleigh.test(rheo_shake_10s_rad)
rayleigh.test(rheo_shake_10s_rad, mu = circular(0))
rayleigh.test(rheo_shake_20s_rad)
rayleigh.test(rheo_shake_20s_rad, mu = circular(0))
rayleigh.test(rheo_shake_post_rad)
rayleigh.test(rheo_shake_post_rad, mu = circular(0))



#Fitak & Jonsen 2017 - compares 10 models (DeltaAIC) and gives best fit: NOTE these tests are done on full data set, not curated. 
# Look at AIC weights for probabilities in model fit. Values fluctuate each time it is run (convergence = 0), so is there a burn-in time, run for 1000 iterations? No change with 1Billion iterations
#refer to table 1 of Fitak & Johnsen 2017 for model code and distribution type for best fit

library(CircMLE) 

# # THE FOLLOWING DESCRIBES THE DISTRIBUTION OF EACH DATA SUBSET AND AIC WEIGHT
# # control
# #uniform M1 (37%)
# circ_mle(rheo_CTL_ps_rad)
# 
# # symmetric unimodal M2B (57%)
# circ_mle(rheo_CTL_10s_rad)
# 
# # symmetric unimodal M2B  (84%) 
# circ_mle(rheo_CTL_20s_rad)
# 
# # symmetric unimodal M2B  (26%)
# circ_mle(rheo_CTL_post_rad)
# 
# 
# # bapta
# #  uniform M1 (30%)
# circ_mle(rheo_bap_ps_rad)
# 
# # homogeneous bimodal M5A (86%)
# circ_mle(rheo_bap_10s_rad)
# 
# # symmetric unimodal M2B (100%)
# circ_mle(rheo_bap_20s_rad)
# 
# # symmetric unimodal M5A (50%)
# circ_mle(rheo_bap_post_rad)



# Watson tests for fit to von Mises or circular uniform distribution; p>0.05 = fail = homogeneous
watson.test(rheo_CTL_ps_rad)


watson.test(rheo_CTL_10s_rad)


watson.test(rheo_CTL_20s_rad)


watson.test(rheo_CTL_post_rad)



watson.test(rheo_shake_ps_rad)


watson.test(rheo_shake_10s_rad)


watson.test(rheo_shake_20s_rad)


watson.test(rheo_shake_post_rad)
 





watson.wheeler.test(list(rheo_CTL_10s$Angle_rad,rheo_CTL_20s$Angle_rad))



watson.wheeler.test(list(rheo_CTL_10s$Angle_rad,rheo_shake_10s$Angle_rad))



watson.wheeler.test(list(rheo_CTL_10s$Angle_rad,rheo_shake_20s$Angle_rad))



watson.wheeler.test(list(rheo_CTL_20s$Angle_rad,rheo_shake_10s$Angle_rad))



watson.wheeler.test(list(rheo_CTL_20s$Angle_rad,rheo_shake_20s$Angle_rad))



watson.wheeler.test(list(rheo_shake_10s$Angle_rad,rheo_shake_20s$Angle_rad))





####---CIRCULAR STATS SHAKE-RECOVERY---####
# use radian data

library(circular)

#full data set, not curated as suggested above
#V-test specifies angle = 0°
rayleigh.test(rheo_CTL$Angle_rad)
rayleigh.test(rheo_CTL$Angle_rad, mu = circular(0))
rayleigh.test(rheo_shake$Angle_rad)
rayleigh.test(rheo_shake$Angle_rad, mu = circular(0))

#0hr Recovery
#pre-stim, post-stim are uniform and not distributed around 0°, unlike 10s & 20s stimulus threatments
rayleigh.test(rheo_CTL_ps0_rad)
rayleigh.test(rheo_CTL_ps0_rad, mu = circular(0))
rayleigh.test(rheo_CTL_10s0_rad)
rayleigh.test(rheo_CTL_10s0_rad, mu = circular(0))
rayleigh.test(rheo_CTL_20s0_rad)
rayleigh.test(rheo_CTL_20s0_rad, mu = circular(0))
rayleigh.test(rheo_CTL_post0_rad)
rayleigh.test(rheo_CTL_post0_rad, mu = circular(0))

#pre-stim, are uniform and not distributed around 0°, unlike threatments. post-stim are unimodal but not clustered aroun 0°
rayleigh.test(rheo_shake_ps0_rad)
rayleigh.test(rheo_shake_ps0_rad, mu = circular(0))
rayleigh.test(rheo_shake_10s0_rad)
rayleigh.test(rheo_shake_10s0_rad, mu = circular(0))
rayleigh.test(rheo_shake_20s0_rad)
rayleigh.test(rheo_shake_20s0_rad, mu = circular(0))
rayleigh.test(rheo_shake_post0_rad)
rayleigh.test(rheo_shake_post0_rad, mu = circular(0))


#2hr Recovery
#pre-stim, post-stim are uniform and not distributed around 0°, unlike 10s & 20s stimulus threatments
rayleigh.test(rheo_CTL_ps2_rad)
rayleigh.test(rheo_CTL_ps2_rad, mu = circular(0))
rayleigh.test(rheo_CTL_10s2_rad)
rayleigh.test(rheo_CTL_10s2_rad, mu = circular(0))
rayleigh.test(rheo_CTL_20s2_rad)
rayleigh.test(rheo_CTL_20s2_rad, mu = circular(0))
rayleigh.test(rheo_CTL_post2_rad)
rayleigh.test(rheo_CTL_post2_rad, mu = circular(0))

#pre-stim, are uniform and not distributed around 0°, unlike threatments. post-stim are unimodal but not clustered aroun 0°
rayleigh.test(rheo_shake_ps2_rad)
rayleigh.test(rheo_shake_ps2_rad, mu = circular(0))
rayleigh.test(rheo_shake_10s2_rad)
rayleigh.test(rheo_shake_10s2_rad, mu = circular(0))
rayleigh.test(rheo_shake_20s2_rad)
rayleigh.test(rheo_shake_20s2_rad, mu = circular(0))
rayleigh.test(rheo_shake_post2_rad)
rayleigh.test(rheo_shake_post2_rad, mu = circular(0))


#4hr Recovery
#pre-stim, post-stim are uniform and not distributed around 0°, unlike 10s & 20s stimulus threatments
rayleigh.test(rheo_CTL_ps4_rad)
rayleigh.test(rheo_CTL_ps4_rad, mu = circular(0))
rayleigh.test(rheo_CTL_10s4_rad)
rayleigh.test(rheo_CTL_10s4_rad, mu = circular(0))
rayleigh.test(rheo_CTL_20s4_rad)
rayleigh.test(rheo_CTL_20s4_rad, mu = circular(0))
rayleigh.test(rheo_CTL_post4_rad)
rayleigh.test(rheo_CTL_post4_rad, mu = circular(0))

#pre-stim, are uniform and not distributed around 0°, unlike threatments. post-stim are unimodal but not clustered aroun 0°
rayleigh.test(rheo_shake_ps4_rad)
rayleigh.test(rheo_shake_ps4_rad, mu = circular(0))
rayleigh.test(rheo_shake_10s4_rad)
rayleigh.test(rheo_shake_10s4_rad, mu = circular(0))
rayleigh.test(rheo_shake_20s4_rad)
rayleigh.test(rheo_shake_20s4_rad, mu = circular(0))
rayleigh.test(rheo_shake_post4_rad)
rayleigh.test(rheo_shake_post4_rad, mu = circular(0))


#8hr Recovery
#pre-stim, post-stim are uniform and not distributed around 0°, unlike 10s & 20s stimulus threatments
rayleigh.test(rheo_CTL_ps8_rad)
rayleigh.test(rheo_CTL_ps8_rad, mu = circular(0))
rayleigh.test(rheo_CTL_10s8_rad)
rayleigh.test(rheo_CTL_10s8_rad, mu = circular(0))
rayleigh.test(rheo_CTL_20s8_rad)
rayleigh.test(rheo_CTL_20s8_rad, mu = circular(0))
rayleigh.test(rheo_CTL_post8_rad)
rayleigh.test(rheo_CTL_post8_rad, mu = circular(0))

#pre-stim, are uniform and not distributed around 0°, unlike threatments. post-stim are unimodal but not clustered aroun 0°
rayleigh.test(rheo_shake_ps8_rad)
rayleigh.test(rheo_shake_ps8_rad, mu = circular(0))
rayleigh.test(rheo_shake_10s8_rad)
rayleigh.test(rheo_shake_10s8_rad, mu = circular(0))
rayleigh.test(rheo_shake_20s8_rad)
rayleigh.test(rheo_shake_20s8_rad, mu = circular(0))
rayleigh.test(rheo_shake_post8_rad)
rayleigh.test(rheo_shake_post8_rad, mu = circular(0))


#48hr Recovery
#pre-stim, post-stim are uniform and not distributed around 0°, unlike 10s & 20s stimulus threatments
rayleigh.test(rheo_CTL_ps48_rad)
rayleigh.test(rheo_CTL_ps48_rad, mu = circular(0))
rayleigh.test(rheo_CTL_10s48_rad)
rayleigh.test(rheo_CTL_10s48_rad, mu = circular(0))
rayleigh.test(rheo_CTL_20s48_rad)
rayleigh.test(rheo_CTL_20s48_rad, mu = circular(0))
rayleigh.test(rheo_CTL_post48_rad)
rayleigh.test(rheo_CTL_post48_rad, mu = circular(0))

#pre-stim, are uniform and not distributed around 0°, unlike threatments. post-stim are unimodal but not clustered aroun 0°
rayleigh.test(rheo_shake_ps48_rad)
rayleigh.test(rheo_shake_ps48_rad, mu = circular(0))
rayleigh.test(rheo_shake_10s48_rad)
rayleigh.test(rheo_shake_10s48_rad, mu = circular(0))
rayleigh.test(rheo_shake_20s48_rad)
rayleigh.test(rheo_shake_20s48_rad, mu = circular(0))
rayleigh.test(rheo_shake_post48_rad)
rayleigh.test(rheo_shake_post48_rad, mu = circular(0))


#Fitak & Jonsen 2017 - compares 10 models (DeltaAIC) and gives best fit: NOTE these tests are done on full data set, not curated. 
# Look at AIC weights for probabilities in model fit. Values fluctuate each time it is run (convergence = 0), so is there a burn-in time, run for 1000 iterations? No change with 1Billion iterations
#refer to table 1 of Fitak & Johnsen 2017 for model code and distribution type for best fit

library(CircMLE) 

# # THE FOLLOWING DESCRIBES THE DISTRIBUTION OF EACH DATA SUBSET AND AIC WEIGHT
# # control
# #uniform M1 (37%)
# circ_mle(rheo_CTL_ps_rad)
# 
# # symmetric unimodal M2B (57%)
# circ_mle(rheo_CTL_10s_rad)
# 
# # symmetric unimodal M2B  (84%) 
# circ_mle(rheo_CTL_20s_rad)
# 
# # symmetric unimodal M2B  (26%)
# circ_mle(rheo_CTL_post_rad)
# 
# 
# # bapta
# #  uniform M1 (30%)
# circ_mle(rheo_bap_ps_rad)
# 
# # homogeneous bimodal M5A (86%)
# circ_mle(rheo_bap_10s_rad)
# 
# # symmetric unimodal M2B (100%)
# circ_mle(rheo_bap_20s_rad)
# 
# # symmetric unimodal M5A (50%)
# circ_mle(rheo_bap_post_rad)



# Watson tests for fit to von Mises or circular uniform distribution; p>0.05 = fail = homogeneous

#0hr
watson.test(rheo_CTL_ps0_rad)
# Test Statistic: 0.0829 
# P-value > 0.10 

watson.test(rheo_CTL_10s0_rad)
# Test Statistic: 2.8254 
# P-value < 0.01 

watson.test(rheo_CTL_20s0_rad)
# Test Statistic: 3.183 
# P-value < 0.01 

watson.test(rheo_CTL_post0_rad)
# Test Statistic: 0.1321 
# P-value > 0.10 


watson.test(rheo_shake_ps0_rad)
# Test Statistic: 0.0829 
# P-value > 0.10 

watson.test(rheo_shake_10s0_rad)
# Test Statistic: 2.8254 
# P-value < 0.01 

watson.test(rheo_shake_20s0_rad)
# Test Statistic: 3.183 
# P-value < 0.01 

watson.test(rheo_shake_post0_rad)
# Test Statistic: 0.1321 
# P-value > 0.10 


#2hr
watson.test(rheo_CTL_ps2_rad)
# Test Statistic: 0.0829 
# P-value > 0.10 

watson.test(rheo_CTL_10s2_rad)
# Test Statistic: 2.8254 
# P-value < 0.01 

watson.test(rheo_CTL_20s2_rad)
# Test Statistic: 3.183 
# P-value < 0.01 

watson.test(rheo_CTL_post2_rad)
# Test Statistic: 0.1321 
# P-value > 0.10 


watson.test(rheo_shake_ps2_rad)
# Test Statistic: 0.0829 
# P-value > 0.10 

watson.test(rheo_shake_10s2_rad)
# Test Statistic: 2.8254 
# P-value < 0.01 

watson.test(rheo_shake_20s2_rad)
# Test Statistic: 3.183 
# P-value < 0.01 

watson.test(rheo_shake_post2_rad)
# Test Statistic: 0.1321 
# P-value > 0.10 


#4hr
watson.test(rheo_CTL_ps4_rad)
# Test Statistic: 0.0829 
# P-value > 0.10 

watson.test(rheo_CTL_10s4_rad)
# Test Statistic: 2.8254 
# P-value < 0.01 

watson.test(rheo_CTL_20s4_rad)
# Test Statistic: 3.183 
# P-value < 0.01 

watson.test(rheo_CTL_post4_rad)
# Test Statistic: 0.1321 
# P-value > 0.10 


watson.test(rheo_shake_ps4_rad)
# Test Statistic: 0.0829 
# P-value > 0.10 

watson.test(rheo_shake_10s4_rad)
# Test Statistic: 2.8254 
# P-value < 0.01 

watson.test(rheo_shake_20s4_rad)
# Test Statistic: 3.183 
# P-value < 0.01 

watson.test(rheo_shake_post4_rad)
# Test Statistic: 0.1321 
# P-value > 0.10 


#8hr
watson.test(rheo_CTL_ps8_rad)
# Test Statistic: 0.0829 
# P-value > 0.10 

watson.test(rheo_CTL_10s8_rad)
# Test Statistic: 2.8254 
# P-value < 0.01 

watson.test(rheo_CTL_20s8_rad)
# Test Statistic: 3.183 
# P-value < 0.01 

watson.test(rheo_CTL_post8_rad)
# Test Statistic: 0.1321 
# P-value > 0.10 


watson.test(rheo_shake_ps8_rad)
# Test Statistic: 0.0829 
# P-value > 0.10 

watson.test(rheo_shake_10s8_rad)
# Test Statistic: 2.8254 
# P-value < 0.01 

watson.test(rheo_shake_20s8_rad)
# Test Statistic: 3.183 
# P-value < 0.01 

watson.test(rheo_shake_post8_rad)
# Test Statistic: 0.1321 
# P-value > 0.10 


#48hr
watson.test(rheo_CTL_ps48_rad)
# Test Statistic: 0.0829 
# P-value > 0.10 

watson.test(rheo_CTL_10s48_rad)
# Test Statistic: 2.8254 
# P-value < 0.01 

watson.test(rheo_CTL_20s48_rad)
# Test Statistic: 3.183 
# P-value < 0.01 

watson.test(rheo_CTL_post48_rad)
# Test Statistic: 0.1321 
# P-value > 0.10 


watson.test(rheo_shake_ps48_rad)
# Test Statistic: 0.0829 
# P-value > 0.10 

watson.test(rheo_shake_10s48_rad)
# Test Statistic: 2.8254 
# P-value < 0.01 

watson.test(rheo_shake_20s48_rad)
# Test Statistic: 3.183 
# P-value < 0.01 

watson.test(rheo_shake_post48_rad)
# Test Statistic: 0.1321 
# P-value > 0.10 





#0hr
watson.wheeler.test(list(rheo_CTL_10s0$Angle_rad,rheo_CTL_20s0$Angle_rad))
# data:  1 and 2
# W = 7.5531, df = 2, p-value = 0.0229


watson.wheeler.test(list(rheo_CTL_10s0$Angle_rad,rheo_shake_10s0$Angle_rad))
# data:  1 and 2
# W = 1.6772, df = 2, p-value = 0.4323


watson.wheeler.test(list(rheo_CTL_10s0$Angle_rad,rheo_shake_20s0$Angle_rad))
# data:  1 and 2
# W = 2.1076, df = 2, p-value = 0.3486


watson.wheeler.test(list(rheo_CTL_20s0$Angle_rad,rheo_shake_10s0$Angle_rad))
# data:  1 and 2
# W = 7.8733, df = 2, p-value = 0.01951


watson.wheeler.test(list(rheo_CTL_20s0$Angle_rad,rheo_shake_20s0$Angle_rad))
# data:  1 and 2
# W = 4.7242, df = 2, p-value = 0.09422


watson.wheeler.test(list(rheo_shake_10s0$Angle_rad,rheo_shake_20s0$Angle_rad))
# data:  1 and 2
# W = 6.6214, df = 2, p-value = 0.03649



#2hr
watson.wheeler.test(list(rheo_CTL_10s2$Angle_rad,rheo_CTL_20s2$Angle_rad))
# data:  1 and 2
# W = 7.5531, df = 2, p-value = 0.0229


watson.wheeler.test(list(rheo_CTL_10s2$Angle_rad,rheo_shake_10s2$Angle_rad))
# data:  1 and 2
# W = 1.6772, df = 2, p-value = 0.4323


watson.wheeler.test(list(rheo_CTL_10s2$Angle_rad,rheo_shake_20s2$Angle_rad))
# data:  1 and 2
# W = 2.1076, df = 2, p-value = 0.3486


watson.wheeler.test(list(rheo_CTL_20s2$Angle_rad,rheo_shake_10s2$Angle_rad))
# data:  1 and 2
# W = 7.8733, df = 2, p-value = 0.01951


watson.wheeler.test(list(rheo_CTL_20s2$Angle_rad,rheo_shake_20s2$Angle_rad))
# data:  1 and 2
# W = 4.7242, df = 2, p-value = 0.09422


watson.wheeler.test(list(rheo_shake_10s2$Angle_rad,rheo_shake_20s2$Angle_rad))
# data:  1 and 2
# W = 6.6214, df = 2, p-value = 0.03649



#4hr
watson.wheeler.test(list(rheo_CTL_10s4$Angle_rad,rheo_CTL_20s4$Angle_rad))
# data:  1 and 2
# W = 7.5531, df = 2, p-value = 0.0229


watson.wheeler.test(list(rheo_CTL_10s4$Angle_rad,rheo_shake_10s4$Angle_rad))
# data:  1 and 2
# W = 1.6772, df = 2, p-value = 0.4323


watson.wheeler.test(list(rheo_CTL_10s4$Angle_rad,rheo_shake_20s4$Angle_rad))
# data:  1 and 2
# W = 2.1076, df = 2, p-value = 0.3486


watson.wheeler.test(list(rheo_CTL_20s4$Angle_rad,rheo_shake_10s4$Angle_rad))
# data:  1 and 2
# W = 7.8733, df = 2, p-value = 0.01951


watson.wheeler.test(list(rheo_CTL_20s4$Angle_rad,rheo_shake_20s4$Angle_rad))
# data:  1 and 2
# W = 4.7242, df = 2, p-value = 0.09422


watson.wheeler.test(list(rheo_shake_10s4$Angle_rad,rheo_shake_20s4$Angle_rad))
# data:  1 and 2
# W = 6.6214, df = 2, p-value = 0.03649



#8hr
watson.wheeler.test(list(rheo_CTL_10s8$Angle_rad,rheo_CTL_20s8$Angle_rad))
# data:  1 and 2
# W = 7.5531, df = 2, p-value = 0.0229


watson.wheeler.test(list(rheo_CTL_10s8$Angle_rad,rheo_shake_10s8$Angle_rad))
# data:  1 and 2
# W = 1.6772, df = 2, p-value = 0.4323


watson.wheeler.test(list(rheo_CTL_10s8$Angle_rad,rheo_shake_20s8$Angle_rad))
# data:  1 and 2
# W = 2.1076, df = 2, p-value = 0.3486


watson.wheeler.test(list(rheo_CTL_20s8$Angle_rad,rheo_shake_10s8$Angle_rad))
# data:  1 and 2
# W = 7.8733, df = 2, p-value = 0.01951


watson.wheeler.test(list(rheo_CTL_20s8$Angle_rad,rheo_shake_20s8$Angle_rad))
# data:  1 and 2
# W = 4.7242, df = 2, p-value = 0.09422


watson.wheeler.test(list(rheo_shake_10s8$Angle_rad,rheo_shake_20s8$Angle_rad))
# data:  1 and 2
# W = 6.6214, df = 2, p-value = 0.03649



#48hr
watson.wheeler.test(list(rheo_CTL_10s48$Angle_rad,rheo_CTL_20s48$Angle_rad))
# data:  1 and 2
# W = 7.5531, df = 2, p-value = 0.0229


watson.wheeler.test(list(rheo_CTL_10s48$Angle_rad,rheo_shake_10s48$Angle_rad))
# data:  1 and 2
# W = 1.6772, df = 2, p-value = 0.4323


watson.wheeler.test(list(rheo_CTL_10s48$Angle_rad,rheo_shake_20s48$Angle_rad))
# data:  1 and 2
# W = 2.1076, df = 2, p-value = 0.3486


watson.wheeler.test(list(rheo_CTL_20s48$Angle_rad,rheo_shake_10s48$Angle_rad))
# data:  1 and 2
# W = 7.8733, df = 2, p-value = 0.01951


watson.wheeler.test(list(rheo_CTL_20s48$Angle_rad,rheo_shake_20s48$Angle_rad))
# data:  1 and 2
# W = 4.7242, df = 2, p-value = 0.09422


watson.wheeler.test(list(rheo_shake_10s48$Angle_rad,rheo_shake_20s48$Angle_rad))
# data:  1 and 2
# W = 6.6214, df = 2, p-value = 0.03649





####---EXTRA STATS TESTS ---####

# #rao test for homogeneity: p>0.05 = fail = homogeneous (chi square)
# rao.spacing.test(rheo_CTL_ps_rad)
# # Test Statistic = 119.486 
# # P-value > 0.10 
# 
# rao.spacing.test(rheo_CTL_10s_rad)
# # Test Statistic = 247.5589 
# # P-value < 0.001 
# 
# rao.spacing.test(rheo_CTL_20s_rad)
# # Test Statistic = 255.6702 
# # P-value < 0.001 
# 
# rao.spacing.test(rheo_CTL_post_rad)
# # Test Statistic = 152.2823 
# # P-value > 0.10 


# #akin to Kolmogorov-Smirnoff test for uniformity.
# kuiper.test(rheo_CTL_ps_rad)
# # Test Statistic:  1.2 
# # P-value > 0.15 
# 
# kuiper.test(rheo_CTL_10s_rad)
# # Test Statistic:  5.4837 
# # P-value < 0.01 
# 
# kuiper.test(rheo_CTL_20s_rad)
# # Test Statistic:  5.7291 
# # P-value < 0.01 
# 
# kuiper.test(rheo_CTL_post_rad)
# # Test Statistic:  1.3613 
# # P-value > 0.15 



# #watson.williams.test (homogeneity of means) - akin to an ANOVA for circular data. If use data subsets, it gives error from unequal concentration (variance) around mean angle (radians)
# watson.williams.test(list(
#   rheo_CTL_ps$Angle_rad, 
#   rheo_CTL_10s$Angle_rad, 
#   rheo_CTL_20s$Angle_rad, 
#   rheo_CTL_post$Angle_rad, 
#   rheo_bap_ps$Angle_rad,
#   rheo_bap_10s$Angle_rad,
#   rheo_bap_20s$Angle_rad,
#   rheo_bap_post$Angle_rad))


# #data:  1 and 2 and 3 and 4 and 5 and 6 and 7 and 8
# F = 6.7006, df1 = 7, df2 = 472, p-value = 1.397e-07
# sample estimates:
#   Circular Data: 
#   Type = angles 
# Units = radians 
# Template = none 
# Modulo = asis 
# Zero = 0 
# Rotation = counter 
# mean of 1   mean of 2   mean of 3   mean of 4   mean of 5   mean of 6   mean of 7   mean of 8 
# 1.03645921  0.01262942 -0.13930426  0.78818385  2.12222519  0.08371670  0.01891986  0.34061405 
# 
# Warning message:
#   In watson.williams.test.default(x, group) :
#   Concentration parameters (0.294, 3.795, 4.412, 0.391, 0.184, 3.048, 4.016, 1.018) not equal between groups. The test might not be applicable
# > 



# #watson-wheeler, non-parametric (chi square), test diffs in mean angle or variance around mean angle (homogeniety of angles)
# watson.wheeler.test(list(
#   rheo_CTL_ps$Angle_rad, 
#   rheo_CTL_10s$Angle_rad, 
#   rheo_CTL_20s$Angle_rad, 
#   rheo_CTL_post$Angle_rad, 
#   rheo_bap_ps$Angle_rad,
#   rheo_bap_10s$Angle_rad,
#   rheo_bap_20s$Angle_rad,
#   rheo_bap_post$Angle_rad))
# 
# # Watson-Wheeler test for homogeneity of angles
# # 
# # data:  1 and 2 and 3 and 4 and 5 and 6 and 7 and 8
# # W = 222.85, df = 14, p-value < 2.2e-16
# 
# 
# 
# # watson test for circular uniformity
# watson.test(rheo_bap_ps$Angle_rad)
# # Test Statistic: 0.0627 
# # P-value > 0.10 
# 
# watson.test(rheo_bap_10s$Angle_rad)
# # Test Statistic: 2.4427 
# # P-value < 0.01 
# 
# watson.test(rheo_bap_20s$Angle_rad)
# # Test Statistic: 3.1558 
# # P-value < 0.01 
# 
# watson.test(rheo_bap_post$Angle_rad)
# # Test Statistic: 0.7261 
# # P-value < 0.01 





# ########### THE FOLLOWING PAIRWISE COMPARISONS NEED A BONFERRONI CORRECTION??? OF ALPHA/6 = 0.00833) ############# p>0.05 = fail = homogeneous
# # Pre and post-stimulus means or variance are not different
# watson.two.test(rheo_CTL_ps$Angle_rad,rheo_CTL_post$Angle_rad)
# # Test Statistic: 0.0219 
# # P-value > 0.10 
# 
# watson.two.test(rheo_CTL_ps$Angle_rad,rheo_bap_ps$Angle_rad)
# # Test Statistic: 0.0503 
# # P-value > 0.10 
# 
# watson.two.test(rheo_bap_ps$Angle_rad,rheo_bap_post$Angle_rad)
# # Test Statistic: 0.3927 
# # P-value < 0.001 
# 

# 
# #RESULT IS THAT ONLY shake-10S IS DIFF THAN shake 20s CTL-10S AND CTL-20S
# watson.two.test(rheo_CTL_10s$Angle_rad,rheo_CTL_20s$Angle_rad)
# # Test Statistic: 0.2142 
# # 0.01 < P-value < 0.05 
# 
# watson.two.test(rheo_CTL_10s$Angle_rad,rheo_bap_10s$Angle_rad)
# # Test Statistic: 0.0627 
# # P-value > 0.10 
# 
# watson.two.test(rheo_CTL_10s$Angle_rad,rheo_bap_20s$Angle_rad)
# # Test Statistic: 0.0757 
# # P-value > 0.10 
# 
# 
# watson.two.test(rheo_CTL_20s$Angle_rad,rheo_bap_10s$Angle_rad)
# # Test Statistic: 0.2296 
# # 0.01 < P-value < 0.05 
# 
# watson.two.test(rheo_CTL_20s$Angle_rad,rheo_bap_20s$Angle_rad)
# # Test Statistic: 0.1711 
# # 0.05 < P-value < 0.10 
# 
# 
# watson.two.test(rheo_bap_10s$Angle_rad,rheo_bap_20s$Angle_rad)
# # Test Statistic: 0.1904 
# # 0.01 < P-value < 0.05 










