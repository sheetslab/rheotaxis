The "Rheotaxis" repository gathers ressources, codes, files and protocol guidelines corresponding to the the paper Lateral line ablation by ototoxic compounds results in distinct rheotaxis profiles in larval zebrafish (Newton & al. 2023).

# 3D printing the microflume and the chamber 
The files for 3D printing the microflume and the chamber are available under 

# Microflume flow calibration 
The python code and associated formula in excel file are available under 

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
   
1) Open RStudio (don’t open R – its way harder to use) then open the R file appropriate for the type of analysis you want to do (file_name.R – not the RData extension)
2) Notice that the active code looks different from comments which have the hashtag and are italicized and grayed out (#install.packages). Comments are useful to make notes for what each line or block of code means, how to modify it, etc. You can toggle lines of code into comments and it will inactivate the code but keep it available as a reference for when you want to track modifications to your code
3) The first few lines of code are where you should set your working drive so that R will save files to a directory specified by youa.MacOS = setwd("/Users/kylenewton/Desktop/Rheotaxis_data") b.WinPC = setwd("C:/Users/kylenewton/Desktop/Rheotaxis_data") c.Pro Tip: from now on, NEVER use spaces in file or directory names! Use the underscore, dash, or period.
4) 4)You must first install the appropriate package into R. Once it is installed you should not have to do it againunless you want to update the package.
   a. In the R Console (lower left pane of RStudio) type: install.packages(“package_name”)  then pressreturnb.
   Or in Source Editor (upper left pane): type install.packages(“package_name”) then press command+enter keys to run the line of code (the cursor must be on the line you want to execute or several lines of code can be selected)c.Or in the (bottom right pane) you can click Packages > Install then search for whatever you want.
   Notice that R autocompletes packages available through the R-CRAN repository.
 5) Once installed, the next few lines of code should load or “library” the packages relevant to your analysis. You must library the package every time you start a new session of R or open a new R data file with the library call:
   a.library(package_name)
   b.command+enter to run the line of code
6) Open your data file and import it into a data frame or a temporary object you can give a descriptive name.a.data_rheotaxis <- read.csv('original_data_file.csv')
   b.If the file is not in the working directory, then yhou need to provide the full path to the file (e.g. “Users/kylenewton/Desktop/Complicated_crap/'original_data_file.csv')
8) Now you can inspect, wrangle, and clean up your data, create data subsets, and do a lot of prep work. Initially, this takes a LOT of time (hours to days for large datasets), but it makes life much easier on the back end. This is when you can find and replace missing values, do some exploratory basic stats, etc.
a.Refer to the cheat sheets for Base R and Tidyverse (Data Import, Data Visualization, Data Tidying) – there are many ways to do the same task but some will be more intuitive than others.
b. I would suggest spending a lot of time getting familiar with data wrangling because it will give you the confidence to execute commands and know that your results are what you intended. Start with a small and intuitive dataset that you care about so that you can check your work then move up to a larger dataset and work with confidence.
c.The key to performing analyses that you want is data subsetting! Well that and using the correct commands.


II. Running the analysis.

1) File 1- Circular graphs representing the mean resultant vector before and during water flow stimulation.

Edit the following lines of code to run the code with your own master csv file.

1) Adjust the data column names.
After running the master csv file
- edit the datatype according to your column names
eg. master.all&Individual <- as.factor(master.all$Individual) for a column named "Individual"

- reconstitute the selection of columns for the aggregate. The aggregate corresponds to the columns you want to select and stack together to perform your analysis. Each number (1,2,3...) is associated with the corresponding column. 
eg. c(1,3,5,6,9:13,17,18) is an aggregate that selects column 1, 3,5,6, 9 to 13, 17 and 18.

- Edit all the columns names that do not fit your own column names eg. "Rheotaxis", "Genotype", "Treatment". Delete for example the column "Treatment" of you do not have any specific treatment in your set up. 

2) optional : edit the name of the variables to fit your own experimental set up.
Eg. We used "Cu" as a variable refering to "Copper". 


3) Edit after each "$" to fit to your column names eg $Genotype

4) Change the names of the variable in "....". Eg. Genotype = "mutant". Change in the "" according to your own genotype as named on your csv file columns (e.g. Genotype="wildtype)



III. Plotting the graphs.

At the end of each code in each files lines of code enable to plot directly the graphs and save them as PDF. 
Edit the color line to have a graph in 3 colors instead of 2.
  
