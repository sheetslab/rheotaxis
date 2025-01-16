
#library packages
library(tidyverse)
library(dplyr)
library(stats)
library(plotly)

library(viridis)

library(lmerTest) #maybe not?
library(lme4)
library(ggplot2)
#library(patchwork)

#do the following steps in Excel if easier
#NOTE: separate SimBA logs (Time_bins_movement_results.csv, & Time_bins_ML_results.csv) by experiment type (bapta, shake, Cu v Neo, etc) then combine into one master event summary file
#NOTE2: summary event file should have the data ordered like this: movement and velocity stuff (Time_bins_movement_results.csv) and event parameters (Time_bins_ML_results.csv)
#NOTE3: to the right of the file name coumn, insert new columns for Treatment and Stimulus (and Recovery time if necessary) and fill in accordingly. These are the "_movt_vel_event.csv" files used below.

#REARRANGE 10S BIN DATA INTO LONG FORMAT
# SimBA_all1 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/SimBA_logs_summary_data/Bapta_movt_vel_event.csv")
# SimBA_all1 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/SimBA_logs_summary_data/Shake_noRec_movt_vel_event.csv")
# SimBA_all1 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/SimBA_logs_summary_data/Shake_REC_movt_vel_event.csv")
names (SimBA_all1)[names(SimBA_all1) == "Animal.1..total.movement.bin.1..cm."] <- "Pre_Stim" 
names (SimBA_all1)[names(SimBA_all1) == "Animal.1..total.movement.bin.2..cm."] <- "Stim_10s" 
names (SimBA_all1)[names(SimBA_all1) == "Animal.1..total.movement.bin.3..cm."] <- "Stim_20s" 
names (SimBA_all1)[names(SimBA_all1) == "Animal.1..total.movement.bin.4..cm."] <- "Post_Stim" 
SimBA_all1a <- SimBA_all1[, c(1:6)] 
# SimBA_all1a <- SimBA_all1[, c(1:7)] # use only with Rec_Time
SimBA_all1b <- gather(SimBA_all1a, 'Pre_Stim', 'Stim_10s', 'Stim_20s', 'Post_Stim', key="Stimulus", value="movement")


# SimBA_all2 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/SimBA_logs_summary_data/Bapta_movt_vel_event.csv")
# SimBA_all2 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/SimBA_logs_summary_data/Shake_noRec_movt_vel_event.csv")
# SimBA_all2 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/SimBA_logs_summary_data/Shake_REC_movt_vel_event.csv")
names (SimBA_all2)[names(SimBA_all2) == "Animal.1..mean.velocity.1..cm."] <- "Pre_Stim" 
names (SimBA_all2)[names(SimBA_all2) == "Animal.1..mean.velocity.2..cm."] <- "Stim_10s" 
names (SimBA_all2)[names(SimBA_all2) == "Animal.1..mean.velocity.3..cm."] <- "Stim_20s" 
names (SimBA_all2)[names(SimBA_all2) == "Animal.1..mean.velocity.4..cm."] <- "Post_Stim" 
SimBA_all2a <- SimBA_all2[, c(1:2,7:10)] 
# SimBA_all2a <- SimBA_all2[, c(1:3,8:11)] # use only with Rec_Time
SimBA_all2b <- gather(SimBA_all2a, 'Pre_Stim', 'Stim_10s', 'Stim_20s', 'Post_Stim', key="Stimulus", value="velocity")
velocity <- SimBA_all2b[, c(4)] 
# velocity <- SimBA_all2b[, c(5)] # use only with Rec_Time

# SimBA_all3 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/SimBA_logs_summary_data/Bapta_movt_vel_event.csv")
# SimBA_all3 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/SimBA_logs_summary_data/Shake_noRec_movt_vel_event.csv")
# SimBA_all3 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/SimBA_logs_summary_data/Shake_REC_movt_vel_event.csv")
SimBA_all3 <- SimBA_all3[, c(1, 2, 13:40)] 
# SimBA_all3 <- SimBA_all3[, c(1:3, 14:41)] # use only with Rec_Time
SimBA_all3a <- SimBA_all3[, c(1:3,10,17,24)] 
# SimBA_all3a <- SimBA_all3[, c(1:4,11,18,25)] # use only with Rec_Time
names (SimBA_all3a)[names(SimBA_all3a) == "Rheotaxis_bin_no_1_number_of_events"] <- "Pre_Stim" 
names (SimBA_all3a)[names(SimBA_all3a) == "Rheotaxis_bin_no_2_number_of_events"] <- "Stim_10s" 
names (SimBA_all3a)[names(SimBA_all3a) == "Rheotaxis_bin_no_3_number_of_events"] <- "Stim_20s" 
names (SimBA_all3a)[names(SimBA_all3a) == "Rheotaxis_bin_no_4_number_of_events"] <- "Post_Stim" 
SimBA_all3aa <- gather(SimBA_all3a, 'Pre_Stim', 'Stim_10s', 'Stim_20s', 'Post_Stim', key="Stimulus", value="number_events")
number_events <- SimBA_all3aa[, c(4)] 
# number_events <- SimBA_all3aa[, c(5)] # use only with Rec_Time

SimBA_all3b <- SimBA_all3[, c(1,2,4,11,18,25)] 
# SimBA_all3b <- SimBA_all3[, c(1:3,5,12,19,26)] # use only with Rec_Time
names (SimBA_all3b)[names(SimBA_all3b) == "Rheotaxis_bin_no_1_total_event_duration"] <- "Pre_Stim" 
names (SimBA_all3b)[names(SimBA_all3b) == "Rheotaxis_bin_no_2_total_event_duration"] <- "Stim_10s" 
names (SimBA_all3b)[names(SimBA_all3b) == "Rheotaxis_bin_no_3_total_event_duration"] <- "Stim_20s" 
names (SimBA_all3b)[names(SimBA_all3b) == "Rheotaxis_bin_no_4_total_event_duration"] <- "Post_Stim" 
SimBA_all3bb <- gather(SimBA_all3b, 'Pre_Stim', 'Stim_10s', 'Stim_20s', 'Post_Stim', key="Stimulus", value="total_event_duration")
total_event_duration <- SimBA_all3bb[, c(4)] 
# total_event_duration <- SimBA_all3bb[, c(5)] # use only with Rec_Time

SimBA_all3c <- SimBA_all3[, c(1,2,5,12,19,26)] 
# SimBA_all3c <- SimBA_all3[, c(1:3,6,13,20,27)] # use only with Rec_Time
names (SimBA_all3c)[names(SimBA_all3c) == "Rheotaxis_bin_no_1_mean_event_duration"] <- "Pre_Stim" 
names (SimBA_all3c)[names(SimBA_all3c) == "Rheotaxis_bin_no_2_mean_event_duration"] <- "Stim_10s" 
names (SimBA_all3c)[names(SimBA_all3c) == "Rheotaxis_bin_no_3_mean_event_duration"] <- "Stim_20s" 
names (SimBA_all3c)[names(SimBA_all3c) == "Rheotaxis_bin_no_4_mean_event_duration"] <- "Post_Stim" 
SimBA_all3cc <- gather(SimBA_all3c, 'Pre_Stim', 'Stim_10s', 'Stim_20s', 'Post_Stim', key="Stimulus", value="mean_event_duration")
# mean_event_duration <- SimBA_all3cc[, c(4)] 
# mean_event_duration <- SimBA_all3cc[, c(5)] # use only with Rec_Time

SimBA_all3d <- SimBA_all3[, c(1,2,7,14,21,28)] 
# SimBA_all3d <- SimBA_all3[, c(1:3,8,15,22,29)] # use only with Rec_Time
names (SimBA_all3d)[names(SimBA_all3d) == "Rheotaxis_bin_no_1_time_first_occurrence"] <- "Pre_Stim" 
names (SimBA_all3d)[names(SimBA_all3d) == "Rheotaxis_bin_no_2_time_first_occurrence"] <- "Stim_10s" 
names (SimBA_all3d)[names(SimBA_all3d) == "Rheotaxis_bin_no_3_time_first_occurrence"] <- "Stim_20s" 
names (SimBA_all3d)[names(SimBA_all3d) == "Rheotaxis_bin_no_4_time_first_occurrence"] <- "Post_Stim" 
SimBA_all3dd <- gather(SimBA_all3d, 'Pre_Stim', 'Stim_10s', 'Stim_20s', 'Post_Stim', key="Stimulus", value="time_first_occurrence")
time_first_occurrence <- SimBA_all3dd[, c(4)] 
# time_first_occurrence <- SimBA_all3dd[, c(5)] # use only with Rec_Time

SimBA_all <- cbind(SimBA_all1b, velocity, number_events, total_event_duration, mean_event_duration, time_first_occurrence)
names (SimBA_all)[names(SimBA_all) == "Video_name"] <- "Individual" 
SimBA_all$Individual <- as.factor(SimBA_all$Individual)
SimBA_all$Treatment <- as.factor(SimBA_all$Treatment)
SimBA_all$Stimulus <- as.factor(SimBA_all$Stimulus)
# SimBA_all$Rec_Time <- as.factor(SimBA_all$Rec_Time) # use only with Rec_Time
SimBA_all$movement <- SimBA_all$movement*10 #convert cm to mm units
SimBA_all$velocity <- SimBA_all$velocity*10 #convert cm to mm units


# write.csv(SimBA_all, file=paste("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/SimBA_logs_summary_data/Bapta_all.csv"))
# write.csv(SimBA_all, file=paste("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/SimBA_logs_summary_data/Shake_noRec_all.csv"))
# write.csv(SimBA_all, file=paste("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/SimBA_logs_summary_data/Shake_REC_all.csv"))


# SimBA_all4 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/SimBA_logs_summary_data/Bapta_all.csv")
# SimBA_all5 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/SimBA_logs_summary_data/Shake_noRec_all.csv")
SimBA_all6 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/SimBA_logs_summary_data/Shake_REC_all.csv")
SimBA_all7 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/SimBA_logs_summary_data/cisplatin_movt_vel_event.csv")


SimBA_all <- SimBA_all7[, c(2:10)] 
SimBA_all <- SimBA_all6[, c(2:11)] # use only with Rec_Time

#fix data classifications
SimBA_all$Individual <- as.factor(SimBA_all$Individual)
SimBA_all$Treatment <- as.factor(SimBA_all$Treatment)
SimBA_all$Stimulus <- as.factor(SimBA_all$Stimulus)
SimBA_all$Rec_Time <- as.factor(SimBA_all$Rec_Time)
str(SimBA_all)



####--- old CUSO4 - NEO subsets---####
# C_N_rheo <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/CuNeo/CuSO4_Neo_rheo_events.csv")
# C_N_rheo <- select(C_N_rheo, File, Treatment, Stimulus, rheo_number_events, rheo_mean_duration, rheo_total_duration, rheo_latency) #select only relevant data
# names (C_N_rheo)[names(C_N_rheo) == "File"] <- "Individual" 
# 
# C_N_movt <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/CuNeo/CuSO4_Neo_movt_vel.csv")
# C_N_movt <- select(C_N_movt, File, Treatment, Stimulus, total_movement, mean_velocity) #select only relevant data
# names (C_N_movt)[names(C_N_movt) == "File"] <- "Individual" 
# C_N_movt1 <- select(C_N_movt, total_movement, mean_velocity) #select only relevant data
# 
# C_N_all <- bind_cols(C_N_rheo, C_N_movt1)
# str(C_N_all)

#fix data classifications
# C_N_all$Individual <- as.factor(C_N_rheo$Individual)
# C_N_all$Treatment <- as.factor(C_N_rheo$Treatment)
# C_N_all$Stimulus <- as.factor(C_N_rheo$Stimulus)

#subsets - CUSO4 - NEOMYCIN rheotaxis events
# C_N_all_nopost <- filter(C_N_all, C_N_all$Stimulus != "post_stim")
# C_N_all_stim_10s <- filter(C_N_all, C_N_all$Stimulus == "stim_10s")
# 
# C_N_all_CTL <- filter(C_N_all, C_N_all$Treatment == "Control")
# C_N_all_CTL_ps <- filter(C_N_all_CTL, C_N_all_CTL$Stimulus == "pre_stim")
# C_N_all_CTL_10s <- filter(C_N_all_CTL, C_N_all_CTL$Stimulus == "stim_10s")
# C_N_all_CTL_20s <- filter(C_N_all_CTL, C_N_all_CTL$Stimulus == "stim_20s")
# C_N_all_CTL_post <- filter(C_N_all_CTL, C_N_all_CTL$Stimulus == "post_stim")
# 
# C_N_all_CuSO4 <- filter(C_N_all, C_N_all$Treatment == "CuSO4")
# C_N_all_CuSO4_ps <- filter(C_N_all_CuSO4, C_N_all_CuSO4$Stimulus == "pre_stim")
# C_N_all_CuSO4_10s <- filter(C_N_all_CuSO4, C_N_all_CuSO4$Stimulus == "stim_10s")
# C_N_all_CuSO4_20s <- filter(C_N_all_CuSO4, C_N_all_CuSO4$Stimulus == "stim_20s")
# C_N_all_CuSO4_post <- filter(C_N_all_CuSO4, C_N_all_CuSO4$Stimulus == "post_stim")
# 
# C_N_all_Neo <- filter(C_N_all, C_N_all$Treatment == "Neo")
# C_N_all_Neo_ps <- filter(C_N_all_Neo, C_N_all_Neo$Stimulus == "pre_stim")
# C_N_all_Neo_10s <- filter(C_N_all_Neo, C_N_all_Neo$Stimulus == "stim_10s")
# C_N_all_Neo_20s <- filter(C_N_all_Neo, C_N_all_Neo$Stimulus == "stim_20s")
# C_N_all_Neo_post <- filter(C_N_all_Neo, C_N_all_Neo$Stimulus == "post_stim")
# 
# #subset for anova of pre-stimulus / no flow data
# C_N_all_ps <- filter(C_N_all, C_N_all$Stimulus == "pre_stim")
# C_N_all_10s <- filter(C_N_all, C_N_all$Stimulus == "stim_10s")



#BAPTA
# SimBA_all_nopost <- filter(SimBA_all, SimBA_all$Stimulus != "Post_Stim")
# SimBA_all_stim_10s <- filter(SimBA_all, SimBA_all$Stimulus == "Stim_10s")
# 
# SimBA_all_CTL <- filter(SimBA_all, SimBA_all$Treatment == "Control")
# SimBA_all_CTL_ps <- filter(SimBA_all_CTL, SimBA_all_CTL$Stimulus == "Pre_Stim")
# SimBA_all_CTL_10s <- filter(SimBA_all_CTL, SimBA_all_CTL$Stimulus == "Stim_10s")
# SimBA_all_CTL_20s <- filter(SimBA_all_CTL, SimBA_all_CTL$Stimulus == "Stim_20s")
# SimBA_all_CTL_post <- filter(SimBA_all_CTL, SimBA_all_CTL$Stimulus == "Post_Stim")
# 
# SimBA_all_Bapta <- filter(SimBA_all, SimBA_all$Treatment == "Bapta")
# SimBA_all_Bapta_ps <- filter(SimBA_all_Bapta, SimBA_all_Bapta$Stimulus == "Pre_Stim")
# SimBA_all_Bapta_10s <- filter(SimBA_all_Bapta, SimBA_all_Bapta$Stimulus == "Stim_10s")
# SimBA_all_Bapta_20s <- filter(SimBA_all_Bapta, SimBA_all_Bapta$Stimulus == "Stim_20s")
# SimBA_all_Bapta_post <- filter(SimBA_all_Bapta, SimBA_all_Bapta$Stimulus == "Post_Stim")
# 
# #subset for anova of pre-stimulus / no flow data
# SimBA_all_ps <- filter(SimBA_all, SimBA_all$Stimulus == "Pre_Stim")
# SimBA_all_10s <- filter(SimBA_all, SimBA_all$Stimulus == "Stim_10s")
# SimBA_all_20s <- filter(SimBA_all, SimBA_all$Stimulus == "Stim_20s")
# SimBA_all_post <- filter(SimBA_all, SimBA_all$Stimulus == "Post_Stim")



#SHAKE - no recovery
SimBA_all_nopost <- filter(SimBA_all, SimBA_all$Stimulus != "Post_Stim")
SimBA_all_stim_10s <- filter(SimBA_all, SimBA_all$Stimulus == "Stim_10s")

SimBA_all_CTL <- filter(SimBA_all, SimBA_all$Treatment == "Control")
SimBA_all_CTL_ps <- filter(SimBA_all_CTL, SimBA_all_CTL$Stimulus == "Pre_Stim")
SimBA_all_CTL_10s <- filter(SimBA_all_CTL, SimBA_all_CTL$Stimulus == "Stim_10s")
SimBA_all_CTL_20s <- filter(SimBA_all_CTL, SimBA_all_CTL$Stimulus == "Stim_20s")
SimBA_all_CTL_post <- filter(SimBA_all_CTL, SimBA_all_CTL$Stimulus == "Post_Stim")

SimBA_all_Shake <- filter(SimBA_all, SimBA_all$Treatment == "Shake")
SimBA_all_Shake_ps <- filter(SimBA_all_Shake, SimBA_all_Shake$Stimulus == "Pre_Stim")
SimBA_all_Shake_10s <- filter(SimBA_all_Shake, SimBA_all_Shake$Stimulus == "Stim_10s")
SimBA_all_Shake_20s <- filter(SimBA_all_Shake, SimBA_all_Shake$Stimulus == "Stim_20s")
SimBA_all_Shake_post <- filter(SimBA_all_Shake, SimBA_all_Shake$Stimulus == "Post_Stim")

#subset for anova of pre-stimulus / no flow data
SimBA_all_ps <- filter(SimBA_all, SimBA_all$Stimulus == "Pre_Stim")
SimBA_all_10s <- filter(SimBA_all, SimBA_all$Stimulus == "Stim_10s")
SimBA_all_20s <- filter(SimBA_all, SimBA_all$Stimulus == "Stim_20s")
SimBA_all_post <- filter(SimBA_all, SimBA_all$Stimulus == "Post_Stim")



#SHAKE - RECOVERY

#Rec_Time
SimBA_all_0hr <- filter(SimBA_all, SimBA_all$Rec_Time == "0hr")
SimBA_all_2hr <- filter(SimBA_all, SimBA_all$Rec_Time == "2hr")
SimBA_all_4hr <- filter(SimBA_all, SimBA_all$Rec_Time == "4hr")
SimBA_all_8hr <- filter(SimBA_all, SimBA_all$Rec_Time == "8hr")
SimBA_all_48hr <- filter(SimBA_all, SimBA_all$Rec_Time == "48hr")


SimBA_all_nopost <- filter(SimBA_all, SimBA_all$Stimulus != "Post_Stim")
SimBA_all_stim_10s <- filter(SimBA_all, SimBA_all$Stimulus == "Stim_10s")

SimBA_all_CTL <- filter(SimBA_all, SimBA_all$Treatment == "Control")
SimBA_all_CTL_ps <- filter(SimBA_all_CTL, SimBA_all_CTL$Stimulus == "Pre_Stim")
SimBA_all_CTL_10s <- filter(SimBA_all_CTL, SimBA_all_CTL$Stimulus == "Stim_10s")
SimBA_all_CTL_20s <- filter(SimBA_all_CTL, SimBA_all_CTL$Stimulus == "Stim_20s")
SimBA_all_CTL_post <- filter(SimBA_all_CTL, SimBA_all_CTL$Stimulus == "Post_Stim")

SimBA_all_Shake <- filter(SimBA_all, SimBA_all$Treatment == "Shake")
SimBA_all_Shake_ps <- filter(SimBA_all_Shake, SimBA_all_Shake$Stimulus == "Pre_Stim")
SimBA_all_Shake_10s <- filter(SimBA_all_Shake, SimBA_all_Shake$Stimulus == "Stim_10s")
SimBA_all_Shake_20s <- filter(SimBA_all_Shake, SimBA_all_Shake$Stimulus == "Stim_20s")
SimBA_all_Shake_post <- filter(SimBA_all_Shake, SimBA_all_Shake$Stimulus == "Post_Stim")

#subset for anova of pre-stimulus / no flow data
SimBA_all_ps <- filter(SimBA_all, SimBA_all$Stimulus == "Pre_Stim")
SimBA_all_10s <- filter(SimBA_all, SimBA_all$Stimulus == "Stim_10s")
SimBA_all_20s <- filter(SimBA_all, SimBA_all$Stimulus == "Stim_20s")
SimBA_all_post <- filter(SimBA_all, SimBA_all$Stimulus == "Post_Stim")




# #0hr -rec
# SimBA_all_CTL_ps_0hr <- filter(SimBA_all_CTL_ps, SimBA_all_CTL_ps$Rec_Time == "0hr")
# SimBA_all_CTL_10s_0hr <- filter(SimBA_all_CTL_10s, SimBA_all_CTL_10s$Rec_Time == "0hr")
# SimBA_all_CTL_20s_0hr <- filter(SimBA_all_CTL_20s, SimBA_all_CTL_20s$Rec_Time == "0hr")
# SimBA_all_CTL_post_0hr <- filter(SimBA_all_CTL_post, SimBA_all_CTL_post$Rec_Time == "0hr")
# 
# SimBA_all_Shake_ps_0hr <- filter(SimBA_all_Shake_ps, SimBA_all_Shake_ps$Rec_Time == "0hr")
# SimBA_all_Shake_10s_0hr <- filter(SimBA_all_Shake_10s, SimBA_all_Shake_10s$Rec_Time == "0hr")
# SimBA_all_Shake_20s_0hr <- filter(SimBA_all_Shake_20s, SimBA_all_Shake_20s$Rec_Time == "0hr")
# SimBA_all_Shake_post_0hr <- filter(SimBA_all_Shake_post, SimBA_all_Shake_post$Rec_Time == "0hr")
# # 
# # #2hr -rec
# SimBA_all_CTL_ps_2hr <- filter(SimBA_all_CTL_ps, SimBA_all_CTL_ps$Rec_Time == "2hr")
# SimBA_all_CTL_10s_2hr <- filter(SimBA_all_CTL_10s, SimBA_all_CTL_10s$Rec_Time == "2hr")
# SimBA_all_CTL_20s_2hr <- filter(SimBA_all_CTL_20s, SimBA_all_CTL_20s$Rec_Time == "2hr")
# SimBA_all_CTL_post_2hr <- filter(SimBA_all_CTL_post, SimBA_all_CTL_post$Rec_Time == "2hr")
# 
# SimBA_all_Shake_ps_2hr <- filter(SimBA_all_Shake_ps, SimBA_all_Shake_ps$Rec_Time == "2hr")
# SimBA_all_Shake_10s_2hr <- filter(SimBA_all_Shake_10s, SimBA_all_Shake_10s$Rec_Time == "2hr")
# SimBA_all_Shake_20s_2hr <- filter(SimBA_all_Shake_20s, SimBA_all_Shake_20s$Rec_Time == "2hr")
# SimBA_all_Shake_post_2hr <- filter(SimBA_all_Shake_post, SimBA_all_Shake_post$Rec_Time == "2hr")
# # 
# # #4hr -rec
# SimBA_all_CTL_ps_4hr <- filter(SimBA_all_CTL_ps, SimBA_all_CTL_ps$Rec_Time == "4hr")
# SimBA_all_CTL_10s_4hr <- filter(SimBA_all_CTL_10s, SimBA_all_CTL_10s$Rec_Time == "4hr")
# SimBA_all_CTL_20s_4hr <- filter(SimBA_all_CTL_20s, SimBA_all_CTL_20s$Rec_Time == "4hr")
# SimBA_all_CTL_post_4hr <- filter(SimBA_all_CTL_post, SimBA_all_CTL_post$Rec_Time == "4hr")
# 
# SimBA_all_Shake_ps_4hr <- filter(SimBA_all_Shake_ps, SimBA_all_Shake_ps$Rec_Time == "4hr")
# SimBA_all_Shake_10s_4hr <- filter(SimBA_all_Shake_10s, SimBA_all_Shake_10s$Rec_Time == "4hr")
# SimBA_all_Shake_20s_4hr <- filter(SimBA_all_Shake_20s, SimBA_all_Shake_20s$Rec_Time == "4hr")
# SimBA_all_Shake_post_4hr <- filter(SimBA_all_Shake_post, SimBA_all_Shake_post$Rec_Time == "4hr")
# # 
# # #8hr -rec
# SimBA_all_CTL_ps_8hr <- filter(SimBA_all_CTL_ps, SimBA_all_CTL_ps$Rec_Time == "8hr")
# SimBA_all_CTL_10s_8hr <- filter(SimBA_all_CTL_10s, SimBA_all_CTL_10s$Rec_Time == "8hr")
# SimBA_all_CTL_20s_8hr <- filter(SimBA_all_CTL_20s, SimBA_all_CTL_20s$Rec_Time == "8hr")
# SimBA_all_CTL_post_8hr <- filter(SimBA_all_CTL_post, SimBA_all_CTL_post$Rec_Time == "8hr")
# 
# SimBA_all_Shake_ps_8hr <- filter(SimBA_all_Shake_ps, SimBA_all_Shake_ps$Rec_Time == "8hr")
# SimBA_all_Shake_10s_8hr <- filter(SimBA_all_Shake_10s, SimBA_all_Shake_10s$Rec_Time == "8hr")
# SimBA_all_Shake_20s_8hr <- filter(SimBA_all_Shake_20s, SimBA_all_Shake_20s$Rec_Time == "8hr")
# SimBA_all_Shake_post_8hr <- filter(SimBA_all_Shake_post, SimBA_all_Shake_post$Rec_Time == "8hr")
# # 
# # #48hr -rec
# SimBA_all_CTL_ps_48hr <- filter(SimBA_all_CTL_ps, SimBA_all_CTL_ps$Rec_Time == "48hr")
# SimBA_all_CTL_10s_48hr <- filter(SimBA_all_CTL_10s, SimBA_all_CTL_10s$Rec_Time == "48hr")
# SimBA_all_CTL_20s_48hr <- filter(SimBA_all_CTL_20s, SimBA_all_CTL_20s$Rec_Time == "48hr")
# SimBA_all_CTL_post_48hr <- filter(SimBA_all_CTL_post, SimBA_all_CTL_post$Rec_Time == "48hr")
# 
# SimBA_all_Shake_ps_48hr <- filter(SimBA_all_Shake_ps, SimBA_all_Shake_ps$Rec_Time == "48hr")
# SimBA_all_Shake_10s_48hr <- filter(SimBA_all_Shake_10s, SimBA_all_Shake_10s$Rec_Time == "48hr")
# SimBA_all_Shake_20s_48hr <- filter(SimBA_all_Shake_20s, SimBA_all_Shake_20s$Rec_Time == "48hr")
# SimBA_all_Shake_post_48hr <- filter(SimBA_all_Shake_post, SimBA_all_Shake_post$Rec_Time == "48hr")
# 

###maybe use [] from base R as it is a more flexible and global syntax for subsetting !!!!!!

### this notation means that the linear model formula = mean duration is a function of number of events, or y(dependent) is a function of x(independent)

# t.test(formula = rheo_mean_duration ~ Treatment,
#        data = C_N_all,
#        alternative = 'two.sided')
# 
# cor.test(formula = ~ rheo_number_events + rheo_mean_duration,
#          data = C_N_all)


#linear model
# rheo.mean.dur.stim.lm <- lm(formula = rheo_mean_duration ~ Stimulus,
#                       data = C_N_all)
# rheo.mean.dur.stim.lm.summary <- summary(rheo.mean.dur.stim.lm)
# summary(rheo.mean.dur.stim.lm)
# 
# anova(rheo.mean.dur.stim.lm)
# rheo.mean.dur.stim.anova <- anova(rheo.mean.dur.stim.lm)
# 
# view(rheo.mean.dur.stim.anova)



# rheo.mean.dur.stim.treat.lm <- lm(formula = rheo_mean_duration ~ Stimulus + Treatment,
#                       data = C_N_all)
# summary(rheo.mean.dur.stim.treat.lm)
# 
# 
# 
# #anova
# rheo.mean.dur.stim.aov <- aov(formula = rheo_mean_duration ~ Stimulus,
#                       data = C_N_all)
# rheo.mean.dur.stim.aov.summary <- summary(rheo.mean.dur.stim.aov)
# rheo.mean.dur.stim.aov.tukey <- TukeyHSD(rheo.mean.dur.stim.aov)
# 
# rheo.mean.dur.treat.aov <- aov(formula = rheo_mean_duration ~ Treatment,
#                       data = C_N_all)
# 
# summary(rheo.mean.dur.treat.aov)
# TukeyHSD(rheo.mean.dur.treat.aov)



# #two-way anova
# rheo.mean.dur.stim.treat.aov2 <- aov(formula = rheo_mean_duration ~ Stimulus + Treatment,
#                       data = C_N_all)
# summary(rheo.mean.dur.stim.treat.aov2)
# TukeyHSD(rheo.mean.dur.stim.treat.aov2)
# 
# 
# rheo.mean.dur.stim.treat.lm2 <- lm(formula = rheo_mean_duration ~ Stimulus + Treatment,
#                        data = C_N_all)
# rheo.mean.dur.stim.treat.lm2.summary <- summary(rheo.mean.dur.stim.treat.lm2)
# summary(rheo.mean.dur.stim.treat.lm2)
# 
# 
# 
# #two-way anova with interactions (IF not significant then go with two-way, but if it is then go with LM)
# rheo.mean.dur.stim.treat.int.aov <- aov(formula = rheo_mean_duration ~ Stimulus * Treatment,
#                      data = C_N_all)
# rheo.mean.dur.stim.treat.int.aov.summary <- summary(rheo.event.stim.treat.int.aov)
# summary(rheo.mean.dur.stim.treat.int.aov)
# 
# rheo.mean.dur.stim.treat.int.lm <- lm(formula = rheo_mean_duration ~ Stimulus * Treatment,
#                           data = C_N_all)
# rheo.mean.dur.stim.treat.int.lm.summary <- summary(rheo.mean.dur.stim.treat.int.lm)
# summary(rheo.mean.dur.stim.treat.int.lm)
# 
# 
# #unbalanced design  = check balance of N
# with(C_N_all,
#      table(Stimulus, Treatment))


# #unbalanced design must use type 2 or 3 anova from car package. type 1 is typical of aov. must create a regression (lm) object. does it really matter??? Calculates SumSquares diff
# event.lm <- lm(formula = rheo_number_events ~ Stimulus + Treatment,
#                data = C_N_all)
# 
# 
# # Type I ANOVA - aov()
# event.I.aov <- aov(event.lm)
# summary(event.I.aov)
# TukeyHSD(event.I.aov)
# 
# # Type II ANOVA - Anova(type = 2)
# event.II.aov <- car::Anova(event.lm, type = 2)
# summary(event.II.aov)
# TukeyHSD(event.II.aov)
# 
# # Type III ANOVA - Anova(type = 3)
# event.III.aov <- car::Anova(event.lm, type = 3)
# summary(event.III.aov)
# TukeyHSD(event.III.aov)
# 
# 
# #add interesting calculated values back to datframe by assigning with $ operator
# names(rheo.event.int.aov)
# head(rheo.event.int.aov)
# View(rheo.event.int.aov)



############### 11/13/20 - chat with DOVI ########

# #repeated measures anova - use lme4 package. this should compare for diffs between Treatments given that each fish was measured repeatedly under diff stimulus conditions, right???
library(lme4)

rheo.mean.dur.stim.treat.int.random.lmer <- lmer(formula = number_events ~ Stimulus + Treatment +  (1|Individual), data = SimBA_all)
summary(rheo.mean.dur.stim.treat.int.random.lmer)

anova(rheo.mean.dur.stim.treat.int.random.lmer)


###this accounts for potential individual random effects using a GLMM (if get errors then use LMER)
rheo.mean.dur.stim.treat.random.glmer <- glmer(formula = mean_event_duration ~ Stimulus + Treatment + Stimulus*Treatment + (1|Individual), data = SimBA_all)
# summary(rheo.mean.dur.stim.treat.int.random.glmer)
# 
# anova(rheo.mean.dur.stim.treat.int.random.glmer)



###this aacounts for potential individual random effects using a LMM
# rheo.mean.dur.lmer.2 <- lmer(formula = rheo_mean_duration ~ Stimulus + Treatment + (1|Individual),
#                             data = C_N_all)
# summary(rheo.mean.dur.lmer.2)

###this is the package name and the call, so it means run this call from this package to make sure it's not overwritten or masked by another package. In this case it is baseR.
#base::anova(rheo_mean_duration.2)


##okay, so how do do I do a two-way mixed model ANOVA, between Treatment groups (independent) and within Stimulus groups (repeated measures)??? Also need post-hoc: t-test, etc.
## running the same call but with p-values (diff in lines 224 and 227)
##note that it ran a typeIII anova for unbalanced desgin - it figured it out automatically!!!!

anova(rheo.mean.dur.stim.treat.random.lmer)

###this post hoc test blows, try using an interaction term??
# TukeyHSD(rheo.mean.dur.lmer)


##RHEOTAXIS MEAN EVENT DURATION, FULL MODEL = FIXED (BETWEEN + WITHIN) + INTERACTION + RANDOM (INDIV)
rheo.mean.dur.stim.treat.int.random.lmer <- lmer(formula = rheo_mean_duration ~ Stimulus + Treatment + Stimulus * Treatment + (1|Individual),
                             data = C_N_all)
summary(rheo.mean.dur.stim.treat.int.random.lmer)

#to get p-values (LME4 wont do it -jerks)
anova(rheo.mean.dur.stim.treat.int.random.lmer)
#view(rheo_mean_duration)

rheo.mean.dur.glm <- glm(formula = rheo_mean_duration ~ Stimulus + Treatment + Stimulus * Treatment,
                                                 data = C_N_all)
rheo.mean.dur.glm.post <- posthoc(rheo.mean.dur.glm)
summary(rheo.mean.dur.glm.post)




# #one-way anova of pre-stim / no flow data only
# # #unbalanced design must use type 2 or 3 anova from car package. type 1 is typical of aov. must create a regression (lm) object. does it really matter??? Calculates SumSquares diff
event.lm <- lm(formula = number_events ~ Treatment,
               data = SimBA_all)

# # Type I ANOVA - aov()
event.I.aov <- aov(event.lm)
summary(event.I.aov)
TukeyHSD(event.I.aov)




####---SUMMARY FIGURES---####
#setwd("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/")


#CuSO4 v NEO
# ggplot(C_N_all_nopost, aes(x=Stimulus, y=total_movement, fill=Treatment)) + 
#   geom_dotplot(binaxis='y', stackdir='center', stackratio=.15, position=position_dodge(.85), binwidth = .25, dotsize=.85, stroke=.1, alpha=1) +
#   scale_fill_manual(values=c("gray", "blue", "green")) +
#   stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="pointrange", position=position_dodge(.85), shape=5, size=.5, color="red", stroke=.5, alpha=1) +
#   theme_classic() +
#   theme(plot.title = element_text(size=12)) +
#   ggtitle("") +
#   xlab("Stimulus") +
#   ylab("Total SB Movement (mm)")
# ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/CuNeo_total_SB_movt.eps", width = 6, height = 4)
# 
# ggplot(C_N_all_nopost, aes(x=Stimulus, y=rheo_number_events, fill=Treatment)) + 
#   geom_dotplot(binaxis='y', stackdir='center', stackratio=.15, position=position_dodge(.75), binwidth = .25, dotsize=.9, stroke=.1, alpha=1) +
#   scale_fill_manual(values=c("gray", "blue", "green")) +
#   stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="pointrange", position=position_dodge(.75), shape=5, size=.5, color="red", stroke=.5, alpha=1) +
#   theme_classic() +
#   theme(plot.title = element_text(size=12)) +
#   ggtitle("") +
#   xlab("Stimulus") +
#   ylab("Mean Number of Events")
# ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/CuNeo_mean_number_events.eps", width = 6, height = 4)
# 
# ggplot(C_N_all_nopost, aes(x=Stimulus, y=rheo_mean_duration, fill=Treatment)) + 
#   geom_dotplot(binaxis='y', stackdir='center', stackratio=.1, position=position_dodge(.85), binwidth = .25, dotsize=.5, stroke=.1, alpha=1) +
#   scale_fill_manual(values=c("gray", "blue", "green")) +
#   stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="pointrange", position=position_dodge(.85), shape=5, size=.5, color="red", stroke=.5, alpha=1) +
#   theme_classic() +
#   theme(plot.title = element_text(size=12)) +
#   ggtitle("") +
#   xlab("Stimulus") +
#   ylab("Mean Duration of Events (s)")
# ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/CuNeo_mean_event_duration.eps", width = 6, height = 4)
# 
# ggplot(C_N_all_nopost, aes(x=Stimulus, y=rheo_total_duration, fill=Treatment)) + 
#   geom_dotplot(binaxis='y', stackdir='center', stackratio=.1, position=position_dodge(.85), binwidth = .25, dotsize=.5, stroke=.1, alpha=1) +
#   scale_fill_manual(values=c("gray", "blue", "green")) +
#   stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="pointrange", position=position_dodge(.85), shape=5, size=.5, color="red", stroke=.5, alpha=1) +
#   theme_classic() +
#   theme(plot.title = element_text(size=12)) +
#   ggtitle("") +
#   xlab("Stimulus") +
#   ylab("Total Duration of Events (s)")
# ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/CuNeo_total_event_duration.eps", width = 6, height = 4)
# 
# 
# #think about how to best subset these data from first time bin
# C_N_all_stim_10s <- filter(C_N_all, C_N_all$Stimulus == "stim_10s")
# C_N_all_stim_10s_1 <- filter(C_N_all_stim_10s, C_N_all_stim_10s$rheo_latency>=10) # nope, just adjusts Y-axis scale
# C_N_all_stim_10s_1$rheo_lat_fix <- C_N_all_stim_10s_1$rheo_latency-10
# 
# 
# ggplot(C_N_all_stim_10s_1, aes(x=Treatment, y=rheo_lat_fix, fill=Treatment)) + 
#   geom_dotplot(binaxis='y', stackdir='center', stackratio=.4, binwidth = .3, dotsize=.5, stroke=.1, alpha=1) +
#   scale_fill_manual(values=c("gray", "blue", "green")) +
#   stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="pointrange", shape=5, size=.5, color="red", stroke=.5, alpha=1) +
#   theme_classic() +
#   theme(plot.title = element_text(size=12)) +
#   ggtitle("") +
#   xlab("Stimulus") +
#   ylab("Mean Latency to First Event")
# ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/CuNeo_mean_latency_1st_event.eps", width = 6, height = 4)



#BAPTA - note that Bapta comes before Control - alphabetical order so adjust order of scale_fill_manual values to match (3rd line)
# ggplot(SimBA_all, aes(x=Stimulus, y=movement, fill=Treatment)) + 
#   geom_dotplot(binaxis='y', stackdir='center', stackratio=1, position=position_dodge(.5), binwidth = .25, dotsize=20, stroke=.1, alpha=.25) +
#   scale_fill_manual(values=c("red", "gray")) +
#   stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="pointrange", position=position_dodge(1), shape=5, size=.5, color="blue", stroke=1, alpha=1) +
#   theme_classic() +
#   theme(plot.title = element_text(size=12)) +
#   ggtitle("") +
#   xlab("Stimulus") +
#   ylab("Total SB Movement (mm)")
# ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/bapta_total_SB_movement.pdf", width = 6, height = 4)
# 
# ggplot(SimBA_all, aes(x=Stimulus, y=velocity, fill=Treatment)) + 
#   geom_dotplot(binaxis='y', stackdir='center', stackratio=.3, position=position_dodge(.75), binwidth = .25, dotsize=2, stroke=.1, alpha=.25) +
#   scale_fill_manual(values=c("red", "gray")) +
#   stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="pointrange", position=position_dodge(1.35), shape=5, size=.5, color="blue", stroke=1, alpha=1) +
#   theme_classic() +
#   theme(plot.title = element_text(size=12)) +
#   ggtitle("") +
#   xlab("Stimulus") +
#   ylab("Total SB Velocity (mm s-1)")
# ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/bapta_total_SB_velocity.pdf", width = 6, height = 4)
# 
# ggplot(SimBA_all, aes(x=Stimulus, y=number_events, fill=Treatment)) + 
#   geom_dotplot(binaxis='y', stackdir='center', stackratio=.1, position=position_dodge(.75), binwidth = .25, dotsize=.9, stroke=.1, alpha=.1) +
#   scale_fill_manual(values=c("red", "gray")) +
#   stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="pointrange", position=position_dodge(1.25), shape=5, size=.5, color="blue", stroke=1, alpha=1) +
#   theme_classic() +
#   theme(plot.title = element_text(size=12)) +
#   ggtitle("") +
#   xlab("Stimulus") +
#   ylab("Mean Number of Events")
# ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/bapta_mean_number_events.pdf", width = 6, height = 4)
# 
# ggplot(SimBA_all, aes(x=Stimulus, y=mean_event_duration, fill=Treatment)) + 
#   geom_dotplot(binaxis='y', stackdir='center', stackratio=.1, position=position_dodge(.85), binwidth = .25, dotsize=1, stroke=.1, alpha=.1) +
#   scale_fill_manual(values=c("red", "gray")) +
#   stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="pointrange", position=position_dodge(1.25), shape=5, size=.5, color="blue", stroke=1, alpha=1) +
#   theme_classic() +
#   theme(plot.title = element_text(size=12)) +
#   ggtitle("") +
#   xlab("Stimulus") +
#   ylab("Mean Duration of Events (s)")
# ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/bapta_mean_event_duration.pdf", width = 6, height = 4)
# 
# ggplot(SimBA_all, aes(x=Stimulus, y=total_event_duration, fill=Treatment)) + 
#   geom_dotplot(binaxis='y', stackdir='center', stackratio=.1, position=position_dodge(.85), binwidth = .25, dotsize=1, stroke=.1, alpha=.1) +
#   scale_fill_manual(values=c("red", "gray")) +
#   stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="pointrange", position=position_dodge(1.3), shape=5, size=.5, color="blue", stroke=1, alpha=1) +
#   theme_classic() +
#   theme(plot.title = element_text(size=12)) +
#   ggtitle("") +
#   xlab("Stimulus") +
#   ylab("Total Duration of Events (s)")
# ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/bapta_total_event_duration.pdf", width = 6, height = 4)
# 
# 
# #think about how to best subset these data from first time bin
# # SimBA_all_stim_10s <- filter(SimBA_all, SimBA_all$Stimulus == "Stim_10s")
# SimBA_all_stim_10s_1 <- filter(SimBA_all_stim_10s, SimBA_all_stim_10s$time_first_occurrence>=10) # nope, just adjusts Y-axis scale
# SimBA_all_stim_10s_1$rheo_lat_fix <- SimBA_all_stim_10s_1$time_first_occurrence-10
# 
# 
# ggplot(SimBA_all_stim_10s_1, aes(x=Treatment, y=rheo_lat_fix, fill=Treatment)) + 
#   geom_dotplot(binaxis='y', stackdir='center', stackratio=.3,  binwidth = .3, dotsize=1, stroke=.1, alpha=.2) +
#   scale_fill_manual(values=c("red", "gray")) +
#   stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="pointrange", shape=5, size=.5, color="blue", stroke=1, alpha=1) +
#   theme_classic() +
#   theme(plot.title = element_text(size=12)) +
#   ggtitle("") +
#   xlab("Stimulus") +
#   ylab("Mean Latency to First Event")
# ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/bapta_mean_latency_1st_event.pdf", width = 6, height = 4)



# #SHAKE - no recovery
# ggplot(SimBA_all, aes(x=Stimulus, y=movement, fill=Treatment)) + 
#   geom_dotplot(binaxis='y', stackdir='center', stackratio=.9, position=position_dodge(.85), binwidth = .25, dotsize=15, stroke=.1, alpha=.1) +
#   scale_fill_manual(values=c("gray", "red")) +
#   stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="pointrange", position=position_dodge(1.25), shape=5, size=.5, color="blue", stroke=1, alpha=1) +
#   theme_classic() +
#   theme(plot.title = element_text(size=12)) +
#   ggtitle("") +
#   xlab("Stimulus") +
#   ylab("Total SB Movement (mm)")
# ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/shake_noRec_total_SB_movement.pdf", width = 6, height = 4)
# 
# #removed velocity graph cuz its the same as total movement but 1/10 scale on Y-axis
# 
# ggplot(SimBA_all, aes(x=Stimulus, y=number_events, fill=Treatment)) + 
#   geom_dotplot(binaxis='y', stackdir='center', stackratio=.1, position=position_dodge(.75), binwidth = .25, dotsize=.9, stroke=.1, alpha=.1) +
#   scale_fill_manual(values=c("gray", "red")) +
#   stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="pointrange", position=position_dodge(1.25), shape=5, size=.5, color="blue", stroke=1, alpha=1) +
#   theme_classic() +
#   theme(plot.title = element_text(size=12)) +
#   ggtitle("") +
#   xlab("Stimulus") +
#   ylab("Mean Number of Events")
# ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/shake_noRec_mean_number_events.pdf", width = 6, height = 4)
# 
# ggplot(SimBA_all, aes(x=Stimulus, y=mean_event_duration, fill=Treatment)) + 
#   geom_dotplot(binaxis='y', stackdir='center', stackratio=.1, position=position_dodge(.85), binwidth = .25, dotsize=1, stroke=.1, alpha=.1) +
#   scale_fill_manual(values=c("gray", "red")) +
#   stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="pointrange", position=position_dodge(1.25), shape=5, size=.5, color="blue", stroke=1, alpha=1) +
#   theme_classic() +
#   theme(plot.title = element_text(size=12)) +
#   ggtitle("") +
#   xlab("Stimulus") +
#   ylab("Mean Duration of Events (s)")
# ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/shake_noRec_mean_event_duration.pdf", width = 6, height = 4)
# 
# ggplot(SimBA_all, aes(x=Stimulus, y=total_event_duration, fill=Treatment)) + 
#   geom_dotplot(binaxis='y', stackdir='center', stackratio=.1, position=position_dodge(.85), binwidth = .25, dotsize=1, stroke=.1, alpha=.1) +
#   scale_fill_manual(values=c("gray", "red")) +
#   stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="pointrange", position=position_dodge(1.3), shape=5, size=.5, color="blue", stroke=1, alpha=1) +
#   theme_classic() +
#   theme(plot.title = element_text(size=12)) +
#   ggtitle("") +
#   xlab("Stimulus") +
#   ylab("Total Duration of Events (s)")
# ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/shake_noRec_total_event_duration.pdf", width = 6, height = 4)
# 
# 
# #think about how to best subset these data from first time bin
# # SimBA_all_stim_10s <- filter(SimBA_all, SimBA_all$Stimulus == "Stim_10s")
# SimBA_all_stim_10s_1 <- filter(SimBA_all_stim_10s, SimBA_all_stim_10s$time_first_occurrence>=10) # nope, just adjusts Y-axis scale
# SimBA_all_stim_10s_1$rheo_lat_fix <- SimBA_all_stim_10s_1$time_first_occurrence-10
# 
# 
# ggplot(SimBA_all_stim_10s_1, aes(x=Treatment, y=rheo_lat_fix, fill=Treatment)) + 
#   geom_dotplot(binaxis='y', stackdir='center', stackratio=.3,  binwidth = .3, dotsize=1, stroke=.1, alpha=.2) +
#   scale_fill_manual(values=c("gray", "red")) +
#   stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="pointrange", shape=5, size=.5, color="blue", stroke=1, alpha=1) +
#   theme_classic() +
#   theme(plot.title = element_text(size=12)) +
#   ggtitle("") +
#   xlab("Stimulus") +
#   ylab("Mean Latency to First Event")
# ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/shake_0hr_Rec_mean_latency_1st_event.pdf", width = 6, height = 4)



####---FIGURES- - -SHAKE - RECOVERY---####
ggplot(SimBA_all_0hr, aes(x=Stimulus, y=movement, fill=Treatment)) + 
  geom_dotplot(binaxis='y', stackdir='center', stackratio=.9, position=position_dodge(.85), binwidth = .25, dotsize=20, stroke=.1, alpha=.1) +
  scale_fill_manual(values=c("gray", "red")) +
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="pointrange", position=position_dodge(1.25), shape=5, size=.5, color="blue", stroke=1, alpha=1) +
  theme_classic() +
  theme(plot.title = element_text(size=12)) +
  ggtitle("0hr REC") +
  xlab("Stimulus") +
  ylab("Total SB Movement (mm)")
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/shake_0hr_Rec_total_SB_movement.pdf", width = 6, height = 4)

#removed velocity graph cuz its the same as total movement but 1/10 scale on Y-axis

ggplot(SimBA_all_0hr, aes(x=Stimulus, y=number_events, fill=Treatment)) + 
  geom_dotplot(binaxis='y', stackdir='center', stackratio=.1, position=position_dodge(.75), binwidth = .25, dotsize=2, stroke=.1, alpha=.1) +
  scale_fill_manual(values=c("gray", "red")) +
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="pointrange", position=position_dodge(1.25), shape=5, size=.5, color="blue", stroke=1, alpha=1) +
  theme_classic() +
  theme(plot.title = element_text(size=12)) +
  ggtitle("0hr REC") +
  xlab("Stimulus") +
  ylab("Mean Number of Events")
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/shake_0hr_Rec_mean_number_events.pdf", width = 6, height = 4)

ggplot(SimBA_all_0hr, aes(x=Stimulus, y=mean_event_duration, fill=Treatment)) + 
  geom_dotplot(binaxis='y', stackdir='center', stackratio=.1, position=position_dodge(.85), binwidth = .25, dotsize=1, stroke=.1, alpha=.1) +
  scale_fill_manual(values=c("gray", "red")) +
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="pointrange", position=position_dodge(1.25), shape=5, size=.5, color="blue", stroke=1, alpha=1) +
  theme_classic() +
  theme(plot.title = element_text(size=12)) +
  ggtitle("0hr REC") +
  xlab("Stimulus") +
  ylab("Mean Duration of Events (s)")
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/shake_0hr_Rec_mean_event_duration.pdf", width = 6, height = 4)

ggplot(SimBA_all_0hr, aes(x=Stimulus, y=total_event_duration, fill=Treatment)) + 
  geom_dotplot(binaxis='y', stackdir='center', stackratio=.1, position=position_dodge(.85), binwidth = .25, dotsize=1, stroke=.1, alpha=.1) +
  scale_fill_manual(values=c("gray", "red")) +
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="pointrange", position=position_dodge(1.3), shape=5, size=.5, color="blue", stroke=1, alpha=1) +
  theme_classic() +
  theme(plot.title = element_text(size=12)) +
  ggtitle("0hr REC") +
  xlab("Stimulus") +
  ylab("Total Duration of Events (s)")
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/shake_0hr_Rec_total_event_duration.pdf", width = 6, height = 4)


#think about how to best subset these data from first time bin
SimBA_all_0hr_stim_10s <- filter(SimBA_all_0hr, SimBA_all_0hr$Stimulus == "Stim_10s")
SimBA_all_0hr_stim_10s_1 <- filter(SimBA_all_0hr_stim_10s, SimBA_all_0hr_stim_10s$time_first_occurrence>=10) # nope, just adjusts Y-axis scale
SimBA_all_0hr_stim_10s_1$rheo_lat_fix <- SimBA_all_0hr_stim_10s_1$time_first_occurrence-10


ggplot(SimBA_all_0hr_stim_10s_1, aes(x=Treatment, y=rheo_lat_fix, fill=Treatment)) + 
  geom_dotplot(binaxis='y', stackdir='center', stackratio=.3,  binwidth = .3, dotsize=1, stroke=.1, alpha=.2) +
  scale_fill_manual(values=c("gray", "red")) +
  stat_summary(fun.data="mean_se", fun.args = list(mult=1), geom="pointrange", shape=5, size=.5, color="blue", stroke=1, alpha=1) +
  theme_classic() +
  theme(plot.title = element_text(size=12)) +
  ggtitle("0hr REC") +
  xlab("Stimulus") +
  ylab("Mean Latency to First Event")
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/shake_0hr_Rec_mean_latency_1st_event.pdf", width = 6, height = 4)




###interaction is significant so could do individual t-tests but too conservative??? if not significant then no need to do pairwise comparisons!!! No-need for bonferroni.
#t.test on subsets - the significant diffs are what drives the effects shown earlier - mostly the two stimulus categories 10v20s, shake v CTL.
#intercept (TreatmentControl), StimPre_stim, Stim

t.test(formula = rheo_mean_duration ~ Treatment,
       data = C_N_all,
       alternative = 'two.sided')






### next time, how do we deal with lack of independence between 10 and 20s period b/c 10s response could influence 20s response (think of energetics/ O2 consumption, cortisol, etc)


####

##RHEOTAXIS NUMBER OF EVENTS, FULL MODEL = FIXED (BETWEEN + WITHIN) + INTERACTION + RANDOM (INDIV)
rheo.num.event.stim.treat.int.random.lmer <- lmer(formula = rheo_number_events ~ Stimulus + Treatment + Stimulus * Treatment + (1|Individual),
                                                 data = C_N_all)
summary(rheo.num.event.stim.treat.int.random.lmer)


#to get p-values (LME4 wont do it -jerks)
anova(rheo.num.event.stim.treat.int.random.lmer)


# #one-way anova of pre-stim / no flow data only
# # #unbalanced design must use type 2 or 3 anova from car package. type 1 is typical of aov. must create a regression (lm) object. does it really matter??? Calculates SumSquares diff
event.lm <- lm(formula = rheo_number_events ~ Treatment,
               data = C_N_all_ps)

# # Type I ANOVA - aov()
event.I.aov <- aov(event.lm)
summary(event.I.aov)
TukeyHSD(event.I.aov)




##RHEOTAXIS LATENCY TO FIRST EVENT, FULL MODEL = FIXED (BETWEEN + WITHIN) + INTERACTION + RANDOM (INDIV)

# NOTE - rerun with time corrected data (subtract 10s, 20s from stimulus times)??? if so, then dont see full range???
# rheo.latency.stim.treat.int.random.lmer <- lmer(formula = rheo_latency ~ Stimulus + Treatment + Stimulus * Treatment + (1|Individual),
#                                                   data = C_N_all)
# summary(rheo.latency.stim.treat.int.random.lmer)

rheo.latency.stim.treat.int.random.lmer <- lmer(formula = rheo_latency ~ Stimulus + Treatment + Stimulus * Treatment + (1|Individual),
                                                data = C_N_all)
summary(rheo.latency.stim.treat.int.random.lmer)


str(C_N_all_10s)


#to get p-values (LME4 wont do it -jerks)
# anova(rheo.latency.stim.treat.int.random.lmer)
anova(rheo.latency.stim.treat.int.random.lmer)


# #one-way anova of pre-stim / no flow data only
# # #unbalanced design must use type 2 or 3 anova from car package. type 1 is typical of aov. must create a regression (lm) object. does it really matter??? Calculates SumSquares diff
event.lm <- lm(formula = rheo_latency ~ Treatment + (1|Individual),
               data = C_N_all_10s)

# # Type I ANOVA - aov()
event.I.aov <- aov(event.lm)
summary(event.I.aov)
TukeyHSD(event.I.aov)


event.lm <- lm(formula = rheo_latency ~ Treatment,
               data = C_N_all_10s)

# # Type III ANOVA - Anova(type = 3)
event.III.aov <- aov(event.lm, type = 3)
summary(event.III.aov)
TukeyHSD(event.III.aov)




##RHEOTAXIS MEAN INTERVAL, FULL MODEL = FIXED (BETWEEN + WITHIN) + INTERACTION + RANDOM (INDIV)

# NOTE - rerun with time corrected data (subtract 10s, 20s from stimulus times)??? if so, then dont see full range???
# rheo.latency.stim.treat.int.random.lmer <- lmer(formula = rheo_latency ~ Stimulus + Treatment + Stimulus * Treatment + (1|Individual),
#                                                   data = C_N_all)
# summary(rheo.latency.stim.treat.int.random.lmer)

rheo.interval.stim.treat.random.lmer <- lmer(formula = rheo_mean_interval ~ Stimulus + Treatment + (1|Individual),
                                            data = C_N_all)
summary(rheo.interval.stim.treat.random.lmer)

#to get p-values (LME4 wont do it -jerks)
# anova(rheo.latency.stim.treat.int.random.lmer)
anova(rheo.interval.stim.treat.random.lmer)


# #one-way anova of pre-stim / no flow data only
# # #unbalanced design must use type 2 or 3 anova from car package. type 1 is typical of aov. must create a regression (lm) object. does it really matter??? Calculates SumSquares diff
event.lm <- lm(formula = rheo_mean_interval ~ Treatment,
               data = C_N_all_ps)

# # Type I ANOVA - aov()
event.I.aov <- aov(event.lm)
summary(event.I.aov)
TukeyHSD(event.I.aov)




##RHEOTAXIS TOTAL MOVEMENT/MEAN VELOCITY, FULL MODEL = FIXED (BETWEEN + WITHIN) + INTERACTION + RANDOM (INDIV)
# # rheo.movt.stim.treat.int.random.lmer <- lmer(formula = total_movement ~ Stimulus + Treatment + Stimulus * Treatment + (1|Individual),
#                                                   data = C_N_all)
# summary(rheo.movt.stim.treat.int.random.lmer)

rheo.movt.stim.treat.random.lmer <- lmer(formula = total_movement ~ Stimulus + Treatment+ Stimulus * Treatment + (1|Individual),
                                             data = C_N_all)
summary(rheo.movt.stim.treat.random.lmer)

#to get p-values (LME4 wont do it -jerks)
# anova(rheo.movt.stim.treat.int.random.lmer)
anova(rheo.movt.stim.treat.random.lmer)

# #one-way anova of pre-stim / no flow data only
# # #unbalanced design must use type 2 or 3 anova from car package. type 1 is typical of aov. must create a regression (lm) object. does it really matter??? Calculates SumSquares diff
event.lm <- lm(formula = total_movement ~ Treatment,
               data = C_N_all_ps)

# # Type I ANOVA - aov()
event.I.aov <- aov(event.lm)
summary(event.I.aov)
TukeyHSD(event.I.aov)




##RHEOTAXIS TOTAL EVENT DURATION, FULL MODEL = FIXED (BETWEEN + WITHIN) + INTERACTION + RANDOM (INDIV)
rheo.total.dur.stim.treat.int.random.lmer <- lmer(formula = rheo_total_duration ~ Stimulus + Treatment + Stimulus * Treatment + (1|Individual),
                                                 data = C_N_all)
summary(rheo.total.dur.stim.treat.int.random.lmer)

#to get p-values (LME4 wont do it -jerks)
anova(rheo.total.dur.stim.treat.int.random.lmer)
#view(rheo_mean_duration)

# #one-way anova of pre-stim / no flow data only
# # #unbalanced design must use type 2 or 3 anova from car package. type 1 is typical of aov. must create a regression (lm) object. does it really matter??? Calculates SumSquares diff
event.lm <- lm(formula = rheo_total_duration ~ Treatment,
               data = C_N_all_ps)

# # Type I ANOVA - aov()
event.I.aov <- aov(event.lm)
summary(event.I.aov)
TukeyHSD(event.I.aov)




####---3/8/21 tests---####
#make sure to library LmerTest not Lme4 if you want p-values!!!!
##RHEOTAXIS proportion time in ROI, FULL MODEL = FIXED (BETWEEN + WITHIN) + INTERACTION + RANDOM (INDIV) 
# rheo.pre.roi.treat.int.random.lmer <- lmer(formula = Time ~ Treatment + Space_bin + Treatment * Space_bin + (1|Individual),
#                                            data = master.pre.rheo) 
# summary(rheo.pre.roi.treat.int.random.lmer)
rheo.pre.roi.treat.int.random.lmer <- lmer(formula = Time ~ Treatment + Space_bin + Treatment * Space_bin + (1|Individual),
                                           data = master.pre) 
summary(rheo.pre.roi.treat.int.random.lmer)

#to get p-values (LME4 wont do it -jerks --USE LmerTest not Lme4)
anova(rheo.pre.roi.treat.int.random.lmer)


#now for stimulus data
rheo.stim.roi.treat.int.random.lmer <- lmer(formula = Time ~ Treatment + Space_bin + Treatment * Space_bin + (1|Individual),
                                           data = master.stim.rheo) 
summary(rheo.stim.roi.treat.int.random.lmer)


anova(rheo.stim.roi.treat.int.random.lmer)




####---Stats ROI fix-11/12/21---#### 

# number of bins in ROI/number individuals in treatment group to normalize comparisons - - oops dont need to standardize
# Use master.all subsets created with rheotaxis3_spatial_use_1D_2D_.R files

#create vectors for a loop 
space.vec <- unique(master.stim.rheo$Space_bin) #how many different things are in the column Space bin?
indiv.vec <- unique(master.stim.rheo$Individual)
#num.bin.rheo.indiv <- table(master.stim.rheo$Individual) #number of rhoetaxis bins for each individual - use if want to compare one time in one ROI to all bins for each indiv fish
roi.data.comp <- data.frame() #create empty data frame to fill with subsequent loop

for(i in 1:length(indiv.vec)){
  for(j in 1:length(space.vec)){

indiv.1 <- indiv.vec[i]
space.1 <- space.vec[j]
treatment.1 <- master.stim.rheo[master.stim.rheo$Individual==indiv.1, 3][1]
bin.count.roi <- nrow(master.stim.rheo[master.stim.rheo$Space_bin==space.1 & master.stim.rheo$Individual==indiv.1, ])
#bin.count.roi <- bin.count.roi/num.bin.rheo.indiv[i] #uncomment this line and 594 to calculate proportion of rheo event in one ROI out of all ROI for individual
#bin.count.roi <- bin.count.roi/200 # normalize by number of time bins
row.vec.1 <- cbind.data.frame(treatment.1, indiv.1, space.1, bin.count.roi)
roi.data.comp <- rbind(roi.data.comp, row.vec.1)
}}

roi.comp.model <- glm(formula = bin.count.roi ~ treatment.1 + space.1 + treatment.1 * space.1, 
                      data = roi.data.comp) 
summary(roi.comp.model)

# summary(glht(roi.comp.model, mcp(space.1="Tukey", treatment.1="Tukey", space.1&)))




# roi.comp.model2 <- aov(formula = bin.count.roi ~ treatment.1 + space.1 + treatment.1 * space.1,
#                       data = roi.data.comp) 
# summary(roi.comp.model)

#Create a post-hoc model
# model_posthoc<-with(data, glm(cont_y ~ groups, family=gaussian))
roi.comp.model.posthoc <- with(roi.data.comp, glm(bin.count.roi ~ treatment.1 + space.1 + treatment.1 * space.1, family=gaussian))



#STILL WORKING ON THIS FOR SPECIFIC POST HOC PAIRWISE COMPARISONS
##Determine the post-hoc comparisons of interest
summary(glht(roi.comp.model.posthoc, 
             linfct = mcp(groups=
                            #Is the difference between these groups different from zero?
                            c("(Back.Control)-(Back.CuSO4)=0",
                              "(Back.Control)-(Back.Neo)=0",
                              "(Back.CuSO4)-(Back.Neo)=0",
                              "(Front.Control)-(Front.CuSO4)=0",
                              "(Front.Control)-(Front.Neo)=0",
                              "(Front.CuSO4)-(Front.Neo)=0",
                              "(Center.Control)-(Center.CuSO4)=0",
                              "(Center.Control)-(Center.Neo)=0",
                              "(Center.CuSO4)-(Center.Neo)=0",
                              "(Left.Control)-(Left.CuSO4)=0",
                              "(Left.Control)-(Left.Neo)=0",
                              "(Left.CuSO4)-(Left.Neo)=0",
                              "(Right.Control)-(Right.CuSO4)=0",
                              "(Right.Control)-(Right.Neo)=0",
                              "(Right.CuSO4)-(Right.Neo)=0"))),test = adjusted("holm"))


# MM <- glm(Y ~ Treatment + 0, data = DeIdentifiedExample)
GG <- posthoc(roi.comp.model)
summary(GG)

library(lsmeans)
emmeans(roi.comp.model, pairwise ~ treatment.1 | space.1)
emmeans(roi.comp.model, pairwise ~ space.1 | treatment.1)


