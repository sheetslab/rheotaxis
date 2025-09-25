The "Rheotaxis" repository gathers ressources, codes, files and protocol guidelines corresponding to the the paper Lateral line ablation by ototoxic compounds results in distinct rheotaxis profiles in larval zebrafish (Newton & al. 2023) and the companion protocol (S.Cohen-Bodenes and al. 2025).

# 3D printing the microflume and the chamber 
The files for 3D printing the microflume and the chamber are available under 

# Microflume flow calibration 
The python code for calibrating the flow rate is available under "CalibrationVideos-checkpoint.ipynb".

#Custom Feature extraction script
Our custom python feature extraction script is available under Zebrafish.py

#Configuration file
You can use our configuration file to load our model with its parameters directly : 

#Example folder
The folder ""provides example videos and associated DLC and Simba files.

# rheotaxis behavioral data analysis 
Protocol and code files for the Rheotaxis Behavioral Set up (Newton & al. 2023). 

This repository contains the R code to replicate the data analysis from the paper Lateral line ablation by ototoxic compounds results in distinct rheotaxis profiles in larval zebrafish (Newton & al. 2023).

The R code generated for the analyses during this study have been initially made available in the Open Science Framework repository, https://osf.io/rvyfz/. Deep Lab Cut software package: https://github.com/DeepLabCut. SimBA: https://github.com/sgoldenlab/simba

1. Download and install the R software on your desktop.
https://mirror.las.iastate.edu/CRAN/

2. Description of the R files.

* File 1: "Rheotaxis1_file_prep_concatenate_CuSO4_Neo_Bapta_Shake.R"
This R code reads the output h5 data from DeepLabCat before importing them into Simba.
The loop will concatenate all machine_learning csv files from SimBA into one master file.
The output "csv" file from the file 1 is the master csv that will be used as intput data file for all the forthcoming analysis on the other R files.

* File 2: "Rheotaxis2_circular_data_CuSO4_Neo_Bapta_Shake.R"
This R code generates circular graphs corresponding to the Fig. 3 of the paper. 

* File 3: "Rheotaxis3_spatial_use_1D_2D_CuSO4_Neo_Bapta_Shake.R"
This R code generates spatial heatmaps corresponding to the Fig. 5 of the paper.

* File 4: "Rheotaxis4_Time_Series_anayses_Spectral_decomposition_Fourier_transform_Cross-correlation_ROI_spatial_use_CuSO4_Neo_Bapta.R" This R code generates spectral analysis graphs corresponding to the Fig. 7 of the paper."

* File 5: "Rheotaxis5_SimBA_summary_movt_GLMMtests1"
This R code performs the statistical analysis (Anova and GLMM) and outputs the corresponding p-values.

2. Run the code to make the data analysis.

I. Installing appropriate libraries and packages. 
   
1) Open RStudio
2) The first few lines of code are where you should set your working drive so that R will save files to a directory specified by youa.MacOS = setwd("/Users/kylenewton/Desktop/Rheotaxis_data") b.WinPC = setwd("C:/Users/kylenewton/Desktop/Rheotaxis_data") c.Pro Tip: from now on, NEVER use spaces in file or directory names! Use the underscore, dash, or period.
4) You must first install the appropriate package into R. Once it is installed you should not have to do it again unless you want to update the package.
   a. In the R Console (lower left pane of RStudio) type: install.packages(“package_name”)  then pressreturnb.
   Or in Source Editor (upper left pane): type install.packages(“package_name”) then press command+enter keys to run the line of code (the cursor must be on the line you want to execute or several lines of code can be selected)c.Or in the (bottom right pane) you can click Packages > Install then search for whatever you want.
    Notice that R autocompletes packages available through the R-CRAN repository.
 6) Once installed, the next few lines of code should load or “library” the packages relevant to your analysis. You must library the package every time you start a new session of R or open a new R data file with the library call:
   a.library(package_name)
   b.command+enter to run the line of code
7) Open your data file and import it into a data frame or a temporary object you can give a descriptive name.a.data_rheotaxis <- read.csv('original_data_file.csv')
   b.If the file is not in the working directory, then yhou need to provide the full path to the file (e.g. “Users/kylenewton/Desktop/Complicated_crap/'original_data_file.csv')
8) Now you can inspect, wrangle, and clean up your data, create data subsets, and do the preprocessing.
  
