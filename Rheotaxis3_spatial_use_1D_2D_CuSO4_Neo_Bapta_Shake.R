setwd("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results")

library(tidyverse)
library(dplyr)
library(plotly)
library(ggplot2)
library(viridis)

# library(stats)
# library(lmerTest)
# library(lme4)



####---DATA-SUBSETTING---####
#preliminary analysis steps require you to temporarily subset data into manageable and useful bits insterad of using the whoe damn gigantic file = time reduction in processing
# if you use the "####--blah---####" it will allow you to jump to different sections using the nav menu at the bottom right corner of this editor window pane

#all these files should have been run through rheotaxis_file_prep_concatenate.R to unify format
master.all1 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/CuNeo/master.all.CuNeo.csv") #CuSO4 v Neomycin - this one is fine?
# 
# master.all2 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/Old_shake/master.all.shake.csv") #Shake - no recovery! - this one is fine?
# master.all3 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/60shake/master.data.shake.rec.csv") #SHAKE - RECOVERY 
# 
master.all4 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/60bapta_test//master.data.bapta.noKS.csv") #BAPTA

master.all5 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/sixty_flow_test/no_KS/master.data.60s.CuNeo.csv") #60s flow test - CuSO4 v Neomycin
# master.all6 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/sixty_flow_test/no_KS/master.data.60s.shake.csv") #60s flow test - Shake


master.all <- master.all4

#might need to add Space_bins if somehow left off saved file
# master.all$Space_bin <- ifelse(master.all$SwimBladder_y<175, "Front",
#                                ifelse(master.all$SwimBladder_y>725, "Back", 
#                                       ifelse(master.all$SwimBladder_x<125 | master.all$SwimBladder_y<175 | master.all$SwimBladder_y>725, "Left",
#                                              ifelse(master.all$SwimBladder_x>775 | master.all$SwimBladder_y<175 | master.all$SwimBladder_y>725, "Right","Center"))))

#make sure that this data is a factor and not simply a string of characters
#master.all3$Rec_Time <- as.factor(master.all3$Rec_Time) #shake-recovery only
master.all$Space_bin <- as.factor(master.all$Space_bin) 
master.all$Individual <- as.factor(master.all$Individual)
master.all$Treatment <- as.factor(master.all$Treatment)
master.all$Stimulus <- as.factor(master.all$Stimulus)
str(master.all) # check and fix data categorization 


#master.all <- master.all <- master.all[, c(2:19)] # if necessary, remove extra first column added when imported file - WITHOUT SPACE BIN
master.all <- master.all <- master.all[, c(2:20)] # if necessary, remove extra first column added when imported file - WITH SPACE BIN


master.rheo <- filter(master.all, master.all$Rheotaxis==1) #must library dplyr if get Error in filter(master.all, master.all$Rheotaxis == 1) : missing values in 'filter' (cuz its using base::filter)
#View(master.rheo) #check a data frame susbset to make sure it has what you want in it


####---CuSO4-NEOMYCIN SUBSETS---####
master.ctl <- filter(master.all, master.all$Treatment=="Control")
master.cu <- filter(master.all, master.all$Treatment=="CuSO4")
master.neo <- filter(master.all, master.all$Treatment=="Neo")

master.pre <- filter(master.all,master.all$Stimulus=="Pre_Stim")
master.stim <- filter(master.all,master.all$Stimulus=="Stimulus")

master.pre.rheo <- filter(master.pre, master.pre$Rheotaxis==1)
master.stim.rheo <- filter(master.stim, master.stim$Rheotaxis==1)

master.ctl.pre <- filter(master.ctl,master.ctl$Stimulus=="Pre_Stim")
master.ctl.stim <- filter(master.ctl,master.ctl$Stimulus=="Stimulus")
master.cu.pre <- filter(master.cu,master.cu$Stimulus=="Pre_Stim")
master.cu.stim <- filter(master.cu,master.cu$Stimulus=="Stimulus")
master.neo.pre <- filter(master.neo,master.neo$Stimulus=="Pre_Stim")
master.neo.stim <- filter(master.neo,master.neo$Stimulus=="Stimulus")

master.shk.pre.rheo <- filter(master.ctl.pre,master.ctl.pre$Rheotaxis==1)
master.ctl.stim.rheo <- filter(master.ctl.stim,master.ctl.stim$Rheotaxis==1)
master.cu.pre.rheo <- filter(master.cu.pre,master.cu.pre$Rheotaxis==1)
master.cu.stim.rheo <- filter(master.cu.stim,master.cu.stim$Rheotaxis==1)
master.neo.pre.rheo <- filter(master.neo.pre,master.neo.pre$Rheotaxis==1)
master.neo.stim.rheo <- filter(master.neo.stim,master.neo.stim$Rheotaxis==1)


####---BAPTA SUBSETS---####
master.ctl <- filter(master.all, master.all$Treatment=="Control")
master.bap <- filter(master.all, master.all$Treatment=="Bapta")

master.pre <- filter(master.all,master.all$Stimulus=="Pre_Stim")
master.stim <- filter(master.all,master.all$Stimulus=="Stimulus")

master.pre.rheo <- filter(master.pre, master.pre$Rheotaxis==1)
master.stim.rheo <- filter(master.stim, master.stim$Rheotaxis==1)

master.ctl.pre <- filter(master.ctl,master.ctl$Stimulus=="Pre_Stim")
master.ctl.stim <- filter(master.ctl,master.ctl$Stimulus=="Stimulus")
master.bap.pre <- filter(master.bap,master.bap$Stimulus=="Pre_Stim")
master.bap.stim <- filter(master.bap,master.bap$Stimulus=="Stimulus")

master.ctl.pre.rheo <- filter(master.ctl.pre,master.ctl.pre$Rheotaxis==1)
master.ctl.stim.rheo <- filter(master.ctl.stim,master.ctl.stim$Rheotaxis==1)
master.bap.pre.rheo <- filter(master.bap.pre,master.bap.pre$Rheotaxis==1)
master.bap.stim.rheo <- filter(master.bap.stim,master.bap.stim$Rheotaxis==1)



####--- SHAKE - SUBSETS---####
master.ctl <- filter(master.all, master.all$Treatment=="Control")
master.shk <- filter(master.all, master.all$Treatment=="Shake")

master.pre <- filter(master.all,master.all$Stimulus=="Pre_Stim")
master.stim <- filter(master.all,master.all$Stimulus=="Stimulus")

master.pre.rheo <- filter(master.pre, master.pre$Rheotaxis==1)
master.stim.rheo <- filter(master.stim, master.stim$Rheotaxis==1)

master.ctl.pre <- filter(master.ctl,master.ctl$Stimulus=="Pre_Stim")
master.ctl.stim <- filter(master.ctl,master.ctl$Stimulus=="Stimulus")
master.shk.pre <- filter(master.shk,master.shk$Stimulus=="Pre_Stim")
master.shk.stim <- filter(master.shk,master.shk$Stimulus=="Stimulus")

master.ctl.pre.rheo <- filter(master.ctl.pre,master.ctl.pre$Rheotaxis==1)
master.ctl.stim.rheo <- filter(master.ctl.stim,master.ctl.stim$Rheotaxis==1)
master.shk.pre.rheo <- filter(master.shk.pre,master.shk.pre$Rheotaxis==1)
master.shk.stim.rheo <- filter(master.shk.stim,master.shk.stim$Rheotaxis==1)


####--- SHAKE-RECOVERY - SUBSETS---####
master.ctl.pre0 <- filter(master.ctl.pre,master.ctl.pre$Rec_Time=="0hr")
master.ctl.stim0 <- filter(master.ctl.stim,master.ctl.stim$Rec_Time=="0hr")
master.shk.pre0 <- filter(master.shk.pre,master.shk.pre$Rec_Time=="0hr")
master.shk.stim0 <- filter(master.shk.stim,master.shk.stim$Rec_Time=="0hr")

master.ctl.pre0.rheo <- filter(master.ctl.pre0,master.ctl.pre0$Rheotaxis==1)
master.ctl.stim0.rheo <- filter(master.ctl.stim0,master.ctl.stim0$Rheotaxis==1)
master.shk.pre0.rheo <- filter(master.shk.pre0,master.shk.pre0$Rheotaxis==1)
master.shk.stim0.rheo <- filter(master.shk.stim0,master.shk.stim0$Rheotaxis==1)

master.ctl.pre2 <- filter(master.ctl.pre,master.ctl.pre$Rec_Time=="2hr")
master.ctl.stim2 <- filter(master.ctl.stim,master.ctl.stim$Rec_Time=="2hr")
master.shk.pre2 <- filter(master.shk.pre,master.shk.pre$Rec_Time=="2hr")
master.shk.stim2 <- filter(master.shk.stim,master.shk.stim$Rec_Time=="2hr")

master.ctl.pre2.rheo <- filter(master.ctl.pre2,master.ctl.pre2$Rheotaxis==1)
master.ctl.stim2.rheo <- filter(master.ctl.stim2,master.ctl.stim2$Rheotaxis==1)
master.shk.pre2.rheo <- filter(master.shk.pre2,master.shk.pre2$Rheotaxis==1)
master.shk.stim2.rheo <- filter(master.shk.stim2,master.shk.stim2$Rheotaxis==1)

master.ctl.pre4 <- filter(master.ctl.pre,master.ctl.pre$Rec_Time=="4hr")
master.ctl.stim4 <- filter(master.ctl.stim,master.ctl.stim$Rec_Time=="4hr")
master.shk.pre4 <- filter(master.shk.pre,master.shk.pre$Rec_Time=="4hr")
master.shk.stim4 <- filter(master.shk.stim,master.shk.stim$Rec_Time=="4hr")

master.ctl.pre4.rheo <- filter(master.ctl.pre4,master.ctl.pre4$Rheotaxis==1)
master.ctl.stim4.rheo <- filter(master.ctl.stim4,master.ctl.stim4$Rheotaxis==1)
master.shk.pre4.rheo <- filter(master.shk.pre4,master.shk.pre4$Rheotaxis==1)
master.shk.stim4.rheo <- filter(master.shk.stim4,master.shk.stim4$Rheotaxis==1)

master.ctl.pre8 <- filter(master.ctl.pre,master.ctl.pre$Rec_Time=="8hr")
master.ctl.stim8 <- filter(master.ctl.stim,master.ctl.stim$Rec_Time=="8hr")
master.shk.pre8 <- filter(master.shk.pre,master.shk.pre$Rec_Time=="8hr")
master.shk.stim8 <- filter(master.shk.stim,master.shk.stim$Rec_Time=="8hr")

master.ctl.pre8.rheo <- filter(master.ctl.pre8,master.ctl.pre8$Rheotaxis==1)
master.ctl.stim8.rheo <- filter(master.ctl.stim8,master.ctl.stim8$Rheotaxis==1)
master.shk.pre8.rheo <- filter(master.shk.pre8,master.shk.pre8$Rheotaxis==1)
master.shk.stim8.rheo <- filter(master.shk.stim8,master.shk.stim8$Rheotaxis==1)

master.ctl.pre48 <- filter(master.ctl.pre,master.ctl.pre$Rec_Time=="48hr")
master.ctl.stim48 <- filter(master.ctl.stim,master.ctl.stim$Rec_Time=="48hr")
master.shk.pre48 <- filter(master.shk.pre,master.shk.pre$Rec_Time=="48hr")
master.shk.stim48 <- filter(master.shk.stim,master.shk.stim$Rec_Time=="48hr")

master.ctl.pre48.rheo <- filter(master.ctl.pre48,master.ctl.pre48$Rheotaxis==1)
master.ctl.stim48.rheo <- filter(master.ctl.stim48,master.ctl.stim48$Rheotaxis==1)
master.shk.pre48.rheo <- filter(master.shk.pre48,master.shk.pre48$Rheotaxis==1)
master.shk.stim48.rheo <- filter(master.shk.stim48,master.shk.stim48$Rheotaxis==1)



####---SINGLE AXIS DENSITY PLOTS---####
#remove na.rm=TRUE??? added 3/11

#CuSO4 v NEO
# these are handy if you imagine the visualizing the data down its X- axis (or Y-axis)
pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/Density_mean_SB_X_CuNeo.pdf", width = 5, height = 4) #creates a PDF
plot(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 5], na.rm=TRUE), col="gray75", xlim=c(0,1000), ylim=c(0,0.015), main="Mean SBx")
lines(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 5], na.rm=TRUE), col="gray50")
lines(density(master.all[master.all$Treatment=="CuSO4" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 5], na.rm=TRUE), col="lightblue")
lines(density(master.all[master.all$Treatment=="CuSO4" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 5], na.rm=TRUE), col="blue")
lines(density(master.all[master.all$Treatment=="Neo" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 5], na.rm=TRUE), col="lightgreen")
lines(density(master.all[master.all$Treatment=="Neo" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 5], na.rm=TRUE), col="green")
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/Density_mean_SB_Y_CuNeo.pdf", width = 5, height = 4)
plot(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 6], na.rm=TRUE), col="gray75", xlim=c(0,1000), ylim=c(0,0.0035), main="Mean SBy")
lines(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 6], na.rm=TRUE), col="gray50")
lines(density(master.all[master.all$Treatment=="CuSO4" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 6], na.rm=TRUE), col="lightblue")
lines(density(master.all[master.all$Treatment=="CuSO4" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 6], na.rm=TRUE), col="blue")
lines(density(master.all[master.all$Treatment=="Neo" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 6], na.rm=TRUE), col="lightgreen")
lines(density(master.all[master.all$Treatment=="Neo" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 6], na.rm=TRUE), col="green")
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/Density_mean_velocity_CuNeo.pdf", width = 5, height = 4)
plot(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 10], na.rm=TRUE), col="gray75", xlim=c(-50,50), ylim=c(0,0.2), main="SB Velocity")
lines(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 10], na.rm=TRUE), col="gray50")
lines(density(master.all[master.all$Treatment=="CuSO4" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 10], na.rm=TRUE), col="lightblue")
lines(density(master.all[master.all$Treatment=="CuSO4" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 10], na.rm=TRUE), col="blue")
lines(density(master.all[master.all$Treatment=="Neo" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 10], na.rm=TRUE), col="lightgreen")
lines(density(master.all[master.all$Treatment=="Neo" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 10], na.rm=TRUE), col="green")
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/Density_mean_acceleration_CuNeo.pdf", width = 5, height = 4)
plot(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 11], na.rm=TRUE), col="gray75", xlim=c(-1000,1000), ylim=c(0, 0.01), main="SB Acceleration")
lines(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 11], na.rm=TRUE), col="gray50")
lines(density(master.all[master.all$Treatment=="CuSO4" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 11], na.rm=TRUE), col="lightblue")
lines(density(master.all[master.all$Treatment=="CuSO4" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 11], na.rm=TRUE), col="blue")
lines(density(master.all[master.all$Treatment=="Neo" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 11], na.rm=TRUE), col="lightgreen")
lines(density(master.all[master.all$Treatment=="Neo" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 11], na.rm=TRUE), col="green")
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/Density_mean_body_angle_CuNeo.pdf", width = 5, height = 4)
plot(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 14], na.rm=TRUE), col="gray75", xlim=c(-1,1), ylim=c(0,2.5), main="Body Angle")
lines(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 14], na.rm=TRUE), col="gray50")
lines(density(master.all[master.all$Treatment=="CuSO4" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 14], na.rm=TRUE), col="lightblue")
lines(density(master.all[master.all$Treatment=="CuSO4" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 14], na.rm=TRUE), col="blue")
lines(density(master.all[master.all$Treatment=="Neo" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 14], na.rm=TRUE), col="lightgreen")
lines(density(master.all[master.all$Treatment=="Neo" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 14], na.rm=TRUE), col="green")
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/Density_mean_resL_CuNeo.pdf", width = 5, height = 4)
plot(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 16], na.rm=TRUE), col="gray75", xlim=c(0.99,1), ylim=c(0,1000), main="Mean Resultant Length")
lines(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 16], na.rm=TRUE), col="gray50")
lines(density(master.all[master.all$Treatment=="CuSO4" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 16], na.rm=TRUE), col="lightblue")
lines(density(master.all[master.all$Treatment=="CuSO4" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 16], na.rm=TRUE), col="blue")
lines(density(master.all[master.all$Treatment=="Neo" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 16], na.rm=TRUE), col="lightgreen")
lines(density(master.all[master.all$Treatment=="Neo" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 16], na.rm=TRUE), col="green")
dev.off()


#BAPTA
# these are handy if you imagine the visualizing the data down its X- axis (or Y-axis)
pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/Density_mean_SB_X_bapta.pdf", width = 5, height = 4) #creates a PDF
plot(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1,  5], na.rm=TRUE), col="gray75", xlim=c(0,1000), ylim=c(0,0.003), main="Mean SBx")
lines(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 5], na.rm=TRUE), col="gray50")
lines(density(master.all[master.all$Treatment=="Bapta" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 5], na.rm=TRUE), col="red")
lines(density(master.all[master.all$Treatment=="Bapta" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 5], na.rm=TRUE), col="red4")
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/Density_mean_SB_Y_bapta.pdf", width = 5, height = 4)
plot(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 6], na.rm=TRUE), col="gray75", xlim=c(0,1000), ylim=c(0,0.004), main="Mean SBy")
lines(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 6], na.rm=TRUE), col="gray50")
lines(density(master.all[master.all$Treatment=="Bapta" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 6], na.rm=TRUE), col="red")
lines(density(master.all[master.all$Treatment=="Bapta" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 6], na.rm=TRUE), col="red4")
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/Density_mean_velocity_bapta.pdf", width = 5, height = 4)
plot(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 10], na.rm=TRUE), col="gray75", xlim=c(-50,50), ylim=c(0,0.35), main="SB Velocity")
lines(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 10], na.rm=TRUE), col="gray50")
lines(density(master.all[master.all$Treatment=="Bapta" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 10], na.rm=TRUE), col="red")
lines(density(master.all[master.all$Treatment=="Bapta" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 10], na.rm=TRUE), col="red4")
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/Density_mean_acceleration_bapta.pdf", width = 5, height = 4)
plot(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 11], na.rm=TRUE), col="gray75", xlim=c(-1000,1000), ylim=c(0, 0.018), main="SB Acceleration")
lines(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 11], na.rm=TRUE), col="gray50")
lines(density(master.all[master.all$Treatment=="Bapta" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 11], na.rm=TRUE), col="red")
lines(density(master.all[master.all$Treatment=="Bapta" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 11], na.rm=TRUE), col="red4")
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/Density_mean_body_angle_bapta.pdf", width = 5, height = 4)
plot(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 14], na.rm=TRUE), col="gray75", xlim=c(-1,1), ylim=c(0,4), main="Body Angle")
lines(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 14], na.rm=TRUE), col="gray50")
lines(density(master.all[master.all$Treatment=="Bapta" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 14], na.rm=TRUE), col="red")
lines(density(master.all[master.all$Treatment=="Bapta" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 14], na.rm=TRUE), col="red4")
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/Density_mean_resL_bapta.pdf", width = 5, height = 4)
plot(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 16], na.rm=TRUE), col="gray75", xlim=c(0.99,1), ylim=c(0,2500), main="Mean Resultant Length")
lines(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 16], na.rm=TRUE), col="gray50")
lines(density(master.all[master.all$Treatment=="Bapta" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 16], na.rm=TRUE), col="red")
lines(density(master.all[master.all$Treatment=="Bapta" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 16], na.rm=TRUE), col="red4")
dev.off()


#SHAKE-recovery (add; & master.all$SRec_Time=="0hr")
# these are handy if you imagine the visualizing the data down its X- axis (or Y-axis)
pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/Density_mean_SB_X_shake.pdf", width = 5, height = 4) #creates a PDF
plot(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1,  5], na.rm=TRUE), col="gray75", xlim=c(0,1000), ylim=c(0,0.003), main="Mean SBx")
lines(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 5], na.rm=TRUE), col="gray50")
lines(density(master.all[master.all$Treatment=="Shake" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 5], na.rm=TRUE), col="red")
lines(density(master.all[master.all$Treatment=="Shake" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 5], na.rm=TRUE), col="red4")
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/Density_mean_SB_Y_shake.pdf", width = 5, height = 4)
plot(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 6], na.rm=TRUE), col="gray75", xlim=c(0,1000), ylim=c(0,0.004), main="Mean SBy")
lines(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 6], na.rm=TRUE), col="gray50")
lines(density(master.all[master.all$Treatment=="Shake" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 6], na.rm=TRUE), col="red")
lines(density(master.all[master.all$Treatment=="Shake" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 6], na.rm=TRUE), col="red4")
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/Density_mean_velocity_shake.pdf", width = 5, height = 4)
plot(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 10], na.rm=TRUE), col="gray75", xlim=c(-50,50), ylim=c(0,0.35), main="SB Velocity")
lines(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 10], na.rm=TRUE), col="gray50")
lines(density(master.all[master.all$Treatment=="Shake" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 10], na.rm=TRUE), col="red")
lines(density(master.all[master.all$Treatment=="Shake" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 10], na.rm=TRUE), col="red4")
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/Density_mean_acceleration_shake.pdf", width = 5, height = 4)
plot(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 11], na.rm=TRUE), col="gray75", xlim=c(-1000,1000), ylim=c(0, 0.018), main="SB Acceleration")
lines(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 11], na.rm=TRUE), col="gray50")
lines(density(master.all[master.all$Treatment=="Shake" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 11], na.rm=TRUE), col="red")
lines(density(master.all[master.all$Treatment=="Shake" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 11], na.rm=TRUE), col="red4")
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/Density_mean_body_angle_shake.pdf", width = 5, height = 4)
plot(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 14], na.rm=TRUE), col="gray75", xlim=c(-1,1), ylim=c(0,4), main="Body Angle")
lines(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 14], na.rm=TRUE), col="gray50")
lines(density(master.all[master.all$Treatment=="Shake" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 14], na.rm=TRUE), col="red")
lines(density(master.all[master.all$Treatment=="Shake" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 14], na.rm=TRUE), col="red4")
dev.off()

pdf("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/Density_mean_resL_shake.pdf", width = 5, height = 4)
plot(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 16], na.rm=TRUE), col="gray75", xlim=c(0.99,1), ylim=c(0,2500), main="Mean Resultant Length")
lines(density(master.all[master.all$Treatment=="Control" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 16], na.rm=TRUE), col="gray50")
lines(density(master.all[master.all$Treatment=="Shake" & master.all$Stimulus=="Pre_Stim" & master.all$Rheotaxis==1, 16], na.rm=TRUE), col="red")
lines(density(master.all[master.all$Treatment=="Shake" & master.all$Stimulus=="Stimulus" & master.all$Rheotaxis==1, 16], na.rm=TRUE), col="red4")
dev.off()



####---X-Y 2D DENSITY Plots / HEATMAPS---####

# library(ggplot2)
# library(hrbrthemes)
# library(viridis)
# library(plotly)


#CuSO4 - NEOMYCIN
XYdensity.ctl.pre <- ggplot(data=master.ctl.pre, aes(x=SwimBladder_x, y=SwimBladder_y) ) +
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,7e-06)) +
  labs(title = "Control Pre-Stimulus") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.ctl.pre
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/XYdensity.ctl.pre.CuNeo.pdf", width = 5, height = 4) #creates a PDF

XYdensity.ctl.pre.rheo <- ggplot(data=master.ctl.pre.rheo, aes(x=SwimBladder_x, y=SwimBladder_y) ) +
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,7e-06)) +
  labs(title = "Control Pre-Stimulus 0°Orientation") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.ctl.pre.rheo
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/XYdensity.ctl.pre.rheo.CuNeo.pdf", width = 5, height = 4) #creates a PDF

XYdensity.ctl.stim <- ggplot(data=master.ctl.stim, aes(x=SwimBladder_x, y=SwimBladder_y) ) + 
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,7e-06)) +
  labs(title = "Control Stimulus") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.ctl.stim
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/XYdensity.ctl.stim.CuNeo.pdf", width = 5, height = 4) #creates a PDF

XYdensity.ctl.stim.rheo <- ggplot(data=master.ctl.stim.rheo, aes(x=SwimBladder_x, y=SwimBladder_y) ) + 
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,7e-06)) +
  labs(title = "Control Stimulus Rheotaxis") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.ctl.stim.rheo
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/XYdensity.ctl.stim.rheo.CuNeo.pdf", width = 5, height = 4) #creates a PDF


XYdensity.cu.pre <- ggplot(data=master.cu.pre, aes(x=SwimBladder_x, y=SwimBladder_y) ) + 
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,7e-06)) +
  labs(title = "CuSO4 Pre-Stimulus") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.cu.pre
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/XYdensity.cu.pre.CuNeo.pdf", width = 5, height = 4) #creates a PDF

XYdensity.cu.pre.rheo <- ggplot(data=master.cu.pre.rheo, aes(x=SwimBladder_x, y=SwimBladder_y) ) + 
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,7e-06)) +
  labs(title = "CuSO4 Pre-Stimulus 0° Orientation") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.cu.pre.rheo
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/XYdensity.cu.pre.rheo.CuNeo.pdf", width = 5, height = 4) #creates a PDF

XYdensity.cu.stim <- ggplot(data=master.cu.stim, aes(x=SwimBladder_x, y=SwimBladder_y) ) + 
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,7e-06)) +
  labs(title = "CuSO4 Stimulus") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.cu.stim
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/XYdensity.cu.stim.CuNeo.pdf", width = 5, height = 4) #creates a PDF

XYdensity.cu.stim.rheo <- ggplot(data=master.cu.stim.rheo, aes(x=SwimBladder_x, y=SwimBladder_y) ) + 
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,7e-06)) +
  labs(title = "CuSO4 Stimulus Rheotaxis") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.cu.stim.rheo
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/XYdensity.cu.stim.rheo.CuNeo.pdf", width = 5, height = 4) #creates a PDF


XYdensity.neo.pre <- ggplot(data=master.neo.pre, aes(x=SwimBladder_x, y=SwimBladder_y) ) + 
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,7e-06)) +
  labs(title = "Neomycin Pre-Stimulus") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.neo.pre
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/XYdensity.neo.pre.CuNeo.pdf", width = 5, height = 4) #creates a PDF

XYdensity.neo.pre.rheo <- ggplot(data=master.neo.pre.rheo, aes(x=SwimBladder_x, y=SwimBladder_y) ) + 
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,7e-06)) +
  labs(title = "Neomycin Pre-Stimulus 0° Orientation") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.neo.pre.rheo
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/XYdensity.neo.pre.rheo.CuNeo.pdf", width = 5, height = 4) #creates a PDF

XYdensity.neo.stim <- ggplot(data=master.neo.stim, aes(x=SwimBladder_x, y=SwimBladder_y) ) + 
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,7e-06)) +
  labs(title = "Neomycin Stimulus") +
  theme_classic() +
  theme(
    legend.position ='right'
  ) 
XYdensity.neo.stim
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/XYdensity.neo.stim.CuNeo.pdf", width = 5, height = 4) #creates a PDF

XYdensity.neo.stim.rheo <- ggplot(data=master.neo.stim.rheo, aes(x=SwimBladder_x, y=SwimBladder_y) ) + 
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,7e-06)) +
  labs(title = "Neomycin Stimulus Rheotaxis") +
  theme_classic() +
  theme(
    legend.position ='right'
  ) 
XYdensity.neo.stim.rheo
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/CuSO4vNEO_graphs/XYdensity.neo.stim.rheo.CuNeo.pdf", width = 5, height = 4) #creates a PDF



#BAPTA
XYdensity.ctl.pre <- ggplot(data=master.ctl.pre, aes(x=SwimBladder_x, y=SwimBladder_y) ) +
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,7e-06)) +
  labs(title = "Control Pre-Stimulus") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.ctl.pre
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/XYdensity.ctl.pre.pdf", width = 5, height = 4) #creates a PDF

XYdensity.ctl.pre.rheo <- ggplot(data=master.ctl.pre.rheo, aes(x=SwimBladder_x, y=SwimBladder_y) ) +
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,7e-06)) +
  labs(title = "Control Pre-Stimulus 0°Orientation") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.ctl.pre.rheo
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/XYdensity.ctl.pre.rheo.pdf", width = 5, height = 4)

XYdensity.ctl.stim <- ggplot(data=master.ctl.stim, aes(x=SwimBladder_x, y=SwimBladder_y) ) + 
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,7e-06)) +
  labs(title = "Control Stimulus") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.ctl.stim
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/XYdensity.ctl.stim.pdf", width = 5, height = 4)

XYdensity.ctl.stim.rheo <- ggplot(data=master.ctl.stim.rheo, aes(x=SwimBladder_x, y=SwimBladder_y) ) + 
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,7e-06)) +
  labs(title = "Control Stimulus Rheotaxis") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.ctl.stim.rheo
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/XYdensity.ctl.stim.rheo.pdf", width = 5, height = 4)

XYdensity.bap.pre <- ggplot(data=master.bap.pre, aes(x=SwimBladder_x, y=SwimBladder_y) ) + 
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,7e-06)) +
  labs(title = "Bapta Pre-Stimulus") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.bap.pre
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/XYdensity.bap.pre.pdf", width = 5, height = 4)

XYdensity.bap.pre.rheo <- ggplot(data=master.bap.pre.rheo, aes(x=SwimBladder_x, y=SwimBladder_y) ) + 
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,7e-06)) +
  labs(title = "Bapta Pre-Stimulus 0° Orientation") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.bap.pre.rheo
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/XYdensity.bap.pre.rheo.pdf", width = 5, height = 4)

XYdensity.bap.stim <- ggplot(data=master.bap.stim, aes(x=SwimBladder_x, y=SwimBladder_y) ) + 
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,7e-06)) +
  labs(title = "Bapta Stimulus") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.bap.stim
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/XYdensity.bap.stim.pdf", width = 5, height = 4)

XYdensity.bap.stim.rheo <- ggplot(data=master.bap.stim.rheo, aes(x=SwimBladder_x, y=SwimBladder_y) ) + 
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,7e-06)) +
  labs(title = "Bapta Stimulus Rheotaxis") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.bap.stim.rheo
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/bapta_graphs/XYdensity.bap.stim.rheo.pdf", width = 5, height = 4)




#SHAKE
XYdensity.ctl.pre <- ggplot(data=master.ctl.pre, aes(x=SwimBladder_x, y=SwimBladder_y) ) +
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,10e-06)) +
  labs(title = "Control Pre-Stimulus") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.ctl.pre
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/XYdensity.ctl.pre.pdf", width = 5, height = 4) #creates a PDF

XYdensity.ctl.pre.rheo <- ggplot(data=master.ctl.pre.rheo, aes(x=SwimBladder_x, y=SwimBladder_y) ) +
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,10e-06)) +
  labs(title = "Control Pre-Stimulus 0°Orientation") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.ctl.pre.rheo
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/XYdensity.ctl.pre.rheo.pdf", width = 5, height = 4)

XYdensity.ctl.stim <- ggplot(data=master.ctl.stim, aes(x=SwimBladder_x, y=SwimBladder_y) ) + 
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,10e-06)) +
  labs(title = "Control Stimulus") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.ctl.stim
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/XYdensity.ctl.stim.pdf", width = 5, height = 4)

XYdensity.ctl.stim.rheo <- ggplot(data=master.ctl.stim.rheo, aes(x=SwimBladder_x, y=SwimBladder_y) ) + 
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,10e-06)) +
  labs(title = "Control Stimulus Rheotaxis") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.ctl.stim.rheo
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/XYdensity.ctl.stim.rheo.pdf", width = 5, height = 4)

XYdensity.shk.pre <- ggplot(data=master.shk.pre, aes(x=SwimBladder_x, y=SwimBladder_y) ) + 
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,10e-06)) +
  labs(title = "Shake Pre-Stimulus") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.shk.pre
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/XYdensity.shk.pre.pdf", width = 5, height = 4)

XYdensity.shk.pre.rheo <- ggplot(data=master.shk.pre.rheo, aes(x=SwimBladder_x, y=SwimBladder_y) ) + 
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,10e-06)) +
  labs(title = "Shake Pre-Stimulus 0° Orientation") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.shk.pre.rheo
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/XYdensity.shk.pre.rheo.pdf", width = 5, height = 4)

XYdensity.shk.stim <- ggplot(data=master.shk.stim, aes(x=SwimBladder_x, y=SwimBladder_y) ) + 
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,10e-06)) +
  labs(title = "Shake Stimulus") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.shk.stim
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/XYdensity.shk.stim.pdf", width = 5, height = 4)

XYdensity.shk.stim.rheo <- ggplot(data=master.shk.stim.rheo, aes(x=SwimBladder_x, y=SwimBladder_y) ) + 
  xlim(0,1000) +
  ylim(0,1000) +
  stat_density_2d(aes(fill = ..density..), 
                  geom = "raster", 
                  contour = FALSE) +
  coord_cartesian(xlim =c(60, 900), ylim = c(900, 60)) +
  scale_fill_viridis(limits = c(0,10e-06)) +
  labs(title = "Shake Stimulus Rheotaxis") +
  theme_classic() +
  theme(
    legend.position='right'
  ) 
XYdensity.shk.stim.rheo
ggsave("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/graphs/shake_graphs/XYdensity.shk.stim.rheo.pdf", width = 5, height = 4)



