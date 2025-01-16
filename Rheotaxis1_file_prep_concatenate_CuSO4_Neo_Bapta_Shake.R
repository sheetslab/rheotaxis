setwd("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/")

library(tidyverse)
library(dplyr)

# library(plotly)
# library(ggplot2)
# library(viridis)
# 
# library(stats)
# library(lmerTest)
# library(lme4)


#this will allow you to read h5 data output files from DLC before you import them into SimBA
# library(rhdf5)
# h5ls('20200318_AB_7dpf_neo_0002DLC_resnet50_1LZF_DLCma_modelJul29shuffle1_100000_bx.h5') #the file you open should be in the working directory you set in line 1
# datah5 <- h5read('20200318_AB_7dpf_neo_0002DLC_resnet50_1LZF_DLCma_modelJul29shuffle1_100000_bx.h5', "/df_with_missing/_i_table/index/H5I_DATASET")

  
#this loop will concatenate all machine_learning csv files from SimBA into one master file  

####------------------STRAIGHT OUTTA SIMBA - loop--------------####
library(circular)

frame_rate.vec <- c(60) #this should combine files of different frame rates if necessary, make sure the files are stored in a folder with "60" in the name
#frame_rate.vec <- c(200,60) #this should combine files of different frame rates if necessary, make sure the files are stored in a folder with "60" and "200" in the name
master.all <- data.frame()

for(h in 1:length(frame_rate.vec)){
  setwd(paste("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/60shake/no_KS/",frame_rate.vec[h],"no_recovery",sep="")) #This looks for folders with "fps" in the name (200fps or 60fps) for source files - can change to a folder path "shake/no_KS/recovery/"
        frame_rate <- frame_rate.vec[h] # number of frames per second
        Time_bin <- 100 # number of msecs
        px_per_mm <- 27.67 #pixels per millimeter conversion from SimBA video parameters
        data.files <- list.files()
        master.data.file <- data.frame()
        for(j in 1:length(data.files)){
          dataKN <- read.csv(data.files[j])
          
          # Treatment <- ifelse(grepl("Neo", data.files[j]), "Neo",
          #                            ifelse(grepl("CuSO4", data.files[j]), "CuSO4", "Control"))
          
          
          # Treatment <- ifelse(grepl("Bapta", data.files[j]), "Bapta","Control")

        
          Treatment <- ifelse(grepl("shake", data.files[j]), "Shake",
                              ifelse(grepl("Shake", data.files[j]), "Shake", "Control"))

          # Rec_Time <- ifelse(grepl("0hr", data.files[j]), "0hr", #make sure to remove comment out on line 132 # dataKN2_agg$TRec_Time <- rep(Treatment, times = nrow(dataKN2_agg))
          #                     ifelse(grepl("2hr", data.files[j]), "2hr",
          #                         ifelse(grepl("4hr", data.files[j]), "4hr",
          #                             ifelse(grepl("8hr", data.files[j]), "8hr", "48hr"))))

          
          dataKN2 <- dataKN[ ,c(8,9,29,30,49,50,64:67,93:94)] #select out and concatenate specific columns in SimBA machine learning csv files according to column number
          
          dataKN2[which(dataKN2$Zebrafish_SwimBladder_x==0 & dataKN2$Zebrafish_SwimBladder_y==0), ] <- NA # this replaces NA for times when DLC does not track fish and you get 0,0 coordinates
          
          
          ms_frame <- (1/frame_rate)*1000 #milliseconds per frame
          numFr_bin <- Time_bin/ms_frame #number of frames per bin to aggregate within our time bins
          num_bin <- nrow(dataKN2)/numFr_bin #number of bins in the data set
          bin_vec <- 1:num_bin
          # t1 <- 1+(bin_num-1)*numFr_bin
          # t2 <- 5+(bin_num-1)*numFr_bin
          dataKN2$binnum <- rep(bin_vec, each=numFr_bin) #adding a column with the bin numbers for each frame
          
          
          #listKN <- as.list(colnames(dataKN2)[2:15])
          #listKN <- paste("dataKN$",listKN,sep="")
          # func2 <- function(x){abs(x[length(x)]-x[1])}
          #dataKN_agg <- aggregate(list(dataKN[ ,2:15]),by = list(dataKN$binnum), mean) #within Bin Average instead of between bin average
          # dataKN2_agg <- aggregate(list(dataKN2[ ,c(1:6,9:12)]),by = list(dataKN2$binnum), mean)
          dataKN2_agg <- aggregate(list(dataKN2[ ,c(1:2,9:11)]),by = list(dataKN2$binnum), mean, na.action=na.omit) #na.action=na.omit ignore NAs in calculation
          dataKN2_agg$body_angle_rad <- atan2(dataKN2_agg$Fish_angle_sin,dataKN2_agg$Fish_angle_cos) #calculate mean body angle 
          dataKN2_agg$body_angle_deg <- deg(dataKN2_agg$body_angle_rad)
          dataKN2_agg$body_angle_deg <- ifelse(dataKN2_agg$body_angle_deg<0, dataKN2_agg$body_angle_deg +360, dataKN2_agg$body_angle_deg) # NEED TO FIX +/- 180 to 0-360
          dataKN2_agg$body_angle_disp <- sqrt((dataKN2_agg$Fish_angle_sin^2)+(dataKN2_agg$Fish_angle_cos^2))
          median_disp <- median(dataKN2_agg$body_angle_disp, na.rm=T)
          avg_disp <- mean(dataKN2_agg$body_angle_disp, na.rm=T)
          
          
          Xmov.vec <- vector() #this section takes SB_x position to determine movement in X axis between bins
          Xmov.vec <- c(Xmov.vec, 0)
          for(r in 2:nrow(dataKN2_agg)){
            # rrr <- dataKN2_agg$Zebrafish_SwimBladder_x[i]-dataKN2_agg$Zebrafish_SwimBladder_x[i-1]
            rrr <- ifelse(is.na(dataKN2_agg$Zebrafish_SwimBladder_x[r]) | is.na(dataKN2_agg$Zebrafish_SwimBladder_x[r-1]), NA, dataKN2_agg$Zebrafish_SwimBladder_x[r]-dataKN2_agg$Zebrafish_SwimBladder_x[r-1])
            Xmov.vec <- c(Xmov.vec, rrr)
          }
          dataKN2_agg$SB_Xmov <- (Xmov.vec)/(px_per_mm)

          Ymov.vec <- vector() #SB_y position to determine movement in Y axis between bins
          Ymov.vec <- c(Ymov.vec, 0)
          for(s in 2:nrow(dataKN2_agg)){
            #sss <- dataKN2_agg$Zebrafish_SwimBladder_y[i]-dataKN2_agg$Zebrafish_SwimBladder_y[i-1]
            sss <- ifelse(is.na(dataKN2_agg$Zebrafish_SwimBladder_y[s]) | is.na(dataKN2_agg$Zebrafish_SwimBladder_y[s-1]), NA, dataKN2_agg$Zebrafish_SwimBladder_y[s]-dataKN2_agg$Zebrafish_SwimBladder_y[s-1])
            Ymov.vec <- c(Ymov.vec, sss)
          }
          dataKN2_agg$SB_Ymov <- (-1)*(Ymov.vec)/(px_per_mm) #need to invert Y- axis to account for DLC origin in upper left corner
    
          dataKN2_agg$SB_mov <- sqrt((dataKN2_agg$SB_Xmov^2)+(dataKN2_agg$SB_Ymov^2))

          SBvel.vec <- vector() #SB_mov (X,Y) position to determine velocity between bins
          SBvel.vec <- c(SBvel.vec, 0)
          for(t in 2:nrow(dataKN2_agg)){
            # ttt <- dataKN2_agg$SB_mov[i]-dataKN2_agg$SB_mov[i-1]
            ttt <- ifelse(is.na(dataKN2_agg$SB_mov[t]) | is.na(dataKN2_agg$SB_mov[t-1]), NA, dataKN2_agg$SB_mov[t]-dataKN2_agg$SB_mov[t-1])
            SBvel.vec <- c(SBvel.vec, ttt)
          }
          dataKN2_agg$SB_vel <- (SBvel.vec)/(Time_bin/1000)        
          
          SBacc.vec <- vector() #SB_vel (X,Y) velocity to determine acceleration between bins
          SBacc.vec <- c(SBacc.vec, 0)
          for(u in 2:nrow(dataKN2_agg)){
            # uuu <- dataKN2_agg$SB_vel[i]-dataKN2_agg$SB_vel[i-1]
            uuu <- ifelse(is.na(dataKN2_agg$SB_vel[u]) | is.na(dataKN2_agg$SB_vel[u-1]), NA, dataKN2_agg$SB_vel[u]-dataKN2_agg$SB_vel[u-1])
            SBacc.vec <- c(SBacc.vec, uuu)
          }
          dataKN2_agg$SB_acc <- (SBacc.vec)/(Time_bin/1000) 
          
          
          dataKN2_agg$Group.1 <- dataKN2_agg$Group.1 *100 #calculate and rename bins
          colnames(dataKN2_agg)[1] <- "Time" #changes name of bins?
          dataKN2_agg$Individual <- rep(data.files[j], times = nrow(dataKN2_agg))
          dataKN2_agg$Treatment <- rep(Treatment, times = nrow(dataKN2_agg))
          # dataKN2_agg$TRec_Time <- rep(Treatment, times = nrow(dataKN2_agg))
          print(j) 
          master.data.file <- rbind(master.data.file, dataKN2_agg) 

        }
        
        #run this line when loop is done to save file for each folder (200fps, 60fps, etc)
        write.csv(master.data.file, file=paste("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/", frame_rate, "master.data.csv", sep = "_"))
        master.all <- rbind(master.all, master.data.file) #combine previous "fps" files into new data frame
}



#this can be used to fix old concatenated files into a unified format
# master.all1 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/CuNeo/master.all.CuNeo.csv") #CuSO4 v Neomycin - this one is fine?
# master.all2 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/Old_shake/master.all.shake.csv") #Shake - no recovery! - this one is fine?
# master.all3 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/60shake/master.data.shake.rec.csv") #SHAKE - RECOVERY 
# master.all4 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/60bapta_test//master.data.bapta.noKS.csv") #BAPTA
# master.all5 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/sixty_flow_test/no_KS/master.data.60s.CuNeo.csv") #60s flow test - CuSO4 v Neomycin
# master.all6 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/sixty_flow_test/no_KS/master.data.60s.shake.csv") #60s flow test - Shake
master.all7 <- read.csv("/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/_60_master.data.csv")

master.all <- master.all7 


#create a new column of Stimulus types based on the time stamp (bin number)
master.all$Stimulus <- ifelse(master.all$Time<10001, "Pre_Stim",
                              ifelse(master.all$Time<30001, "Stimulus", "Post_Stim"))
                                  
master.all$Rheotaxis <- ifelse(master.all$Probability_Rheotaxis<0.5, 0, 1) #create a new column of Rheotaxis (yes=1, no=0) based on probablity caluculated by SimBA

#rename columns
names (master.all)[names(master.all) == "Zebrafish_SwimBladder_x"] <- "SwimBladder_x" #rename columns to somethign more manageable
names (master.all)[names(master.all) == "Zebrafish_SwimBladder_y"] <- "SwimBladder_y"
names (master.all)[names(master.all) == "body_angle_rad"] <- "Fish_mean_angle_rad"
names (master.all)[names(master.all) == "body_angle_deg"] <- "Fish_mean_angle_deg"
names (master.all)[names(master.all) == "body_angle_disp"] <- "Mean_angle_dispersion" #MAYBE CHANGE TO MEAN RESULTANT LENGTH????


#create spatial ROIs based on X,Y position of swim bladder .
master.all$Space_bin <- ifelse(master.all$SwimBladder_y<175, "Front",
                               ifelse(master.all$SwimBladder_y>725, "Back", 
                                      ifelse(master.all$SwimBladder_x<125 | master.all$SwimBladder_y<175 | master.all$SwimBladder_y>725, "Left",
                                             ifelse(master.all$SwimBladder_x>775 | master.all$SwimBladder_y<175 | master.all$SwimBladder_y>725, "Right","Center"))))

#format data type - but wont seem to save??
master.all$Space_bin <- as.factor(master.all$Space_bin) #make sure that these data are a factor and not simply a string of characters
master.all$Individual <- as.factor(master.all$Individual)
master.all$Treatment <- as.factor(master.all$Treatment)
master.all$Stimulus <- as.factor(master.all$Stimulus)
# master.all$Rec_Time <- as.factor(master.all$Rec_Time) #shake-recovery only
str(master.all) # check and fix data categorization 


# master.all <- master.all[, c(1,15:17,2,3,10:14,4,5,7:9,6,18)] #selecting out specific columns based on position number
# master.all <- master.all[, c(2,16,17,19,3,4,11:15,5,6,8:10,7,20,21,18)] #selecting out specific columns based on position number #shake-recovery time
master.all <- master.all[, c(2,16,17,18,3,4,11:15,5,6,8:10,7,19,20)] #selecting out specific columns based on position number
#master.all <- master.all[, c(3:21)]


#write.csv(disp.all, file="disp.all.csv")
# write.csv(master.all, file="/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/master.all.csv")
# write.csv(master.all, file="/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/60shake/no_KS//master.data.shake.rec.csv") #master.all1 #SHAKE - RECOVERY no rheotaxis & probability columns should fix
# write.csv(master.all, file="/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/60bapta_test//master.data.bapta.noKS.csv") #master.all2  #BAPTA
# write.csv(master.all, file="/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/sixty_flow_test/no_KS/master.data.60s.shake.csv") #master.all3  #60s flow test - Shake
# write.csv(master.all, file="/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/sixty_flow_test/no_KS/master.data.60s.CuNeo.csv") #master.all4  #60s flow test - CuSO4 v Neomycin
# write.csv(master.all, file="/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/CuNeo/master.all.CuNeo.csv") #master.all5  #CuSO4 v Neomycin - this one is fine?
write.csv(master.all, file="/Users/kylenewton/Desktop/RHEOTAXIS:SHAKE/Rheotaxis_data/machine_results/Old_shake/master.all.shake.fixed.csv")

# this should be it, the master.all.csv file is the master data file that you will subset and analyze 



####---DATA-SUBSETTING---####
#preliminary analysis steps require you to temporarily subset data into manageable and useful bits instead of using the whole damn gigantic file = time reduction in processing
# if you use the "####--blah---####" it will allow you to jump to different sections using the nav menu at the bottom right corner of this editor window pane

master.rheo <- filter(master.all, master.all$Rheotaxis==1) #must library dplyr if get Error in filter(master.all, master.all$Rheotaxis == 1) : missing values in 'filter' (cuz its using base::filter)
#View(master.rheo) #check a data frame susbset to make sure it has what you want in it


#CuSO4, NEOMYCIN  - USE CORRECT master.all.csv
master.ctl <- filter(master.all, master.all$Treatment=="Control")
master.cu <- filter(master.all, master.all$Treatment=="CuSO4")
master.neo <-filter(master.all, master.all$Treatment=="Neo")

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

master.ctl.pre.rheo <- filter(master.ctl.pre,master.ctl.pre$Rheotaxis==1)
master.ctl.stim.rheo <- filter(master.ctl.stim,master.ctl.stim$Rheotaxis==1)
master.cu.pre.rheo <- filter(master.cu.pre,master.cu.pre$Rheotaxis==1)
master.cu.stim.rheo <- filter(master.cu.stim,master.cu.stim$Rheotaxis==1)
master.neo.pre.rheo <- filter(master.neo.pre,master.neo.pre$Rheotaxis==1)
master.neo.stim.rheo <- filter(master.neo.stim,master.neo.stim$Rheotaxis==1)
str()



#BAPTA - USE CORRECT master.all.csv
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



#SHAKE - NO RECOVERY - USE CORRECT master.all.csv
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


#SHAKE-RECOVERY - USE CORRECT master.all.csv
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
master.shk.pre48.rheo <- filter(master.shk.pre48,master.shk.pre48$Rheotaxis==1) # only six data points?????
# master.shk.pre48.norheo <- filter(master.shk.pre48,master.shk.pre48$Rheotaxis==0) 
master.shk.stim48.rheo <- filter(master.shk.stim48,master.shk.stim48$Rheotaxis==1)












