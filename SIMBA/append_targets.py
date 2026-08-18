import os, glob
import pandas as pd

#PATH TO THE FOLDER CONTAINING YOUR NEW FEATURE FILES (THAT DO NOT HAVE THE TARGET BEHAVIOR ANNOTATIONS: I.E., MISSING YOUR RHEOTAXIS COLUMN)
featureFilePath = r'Z:\Classifiers\Anogenital_081120\project_folder\csv\features_extracted'

#PATH TO THE TARGET FOLDER CONTAINING FILES WITH ANNOTATATIONS (FILES THAT DO HAVE THE TARGET BEHAVIOR ANNOTATIONS, BUT THE OLD, WRONG FEATURES)
targetFilePath = r'Z:\DeepLabCut\DLC_extract\Troubleshooting\SimBA_projects_troubleshooting\Simon_031618\project_folder\csv\targets_inserted'

#PATH TO THE FOLDER WHERE YOUR NEW FILES (WITH YOUR NEW FEATURES, AND YOUR RHETAXIS ANNOTATIONS) MAKE SURE THIS FOLDER EXIST  BEFORE RUNNING
outputPath = r'Z:\Classifiers\Anogenital_081120\project_folder\csv\targets_inserted'

#THE NAME OF YOUR CLASSIFIER (IE, Rheotaxsis)
classifierName = 'anogenital_prediction'

featureFiles = glob.glob(featureFilePath + '/*.csv')
targetFiles = glob.glob(targetFilePath + '/*.csv')

for file in targetFiles:
    currTargetDf = pd.read_csv(file, usecols = [classifierName, 'scorer'])
    indexNos, classifications = list(currTargetDf.scorer), list(currTargetDf[classifierName])
    print(indexNos)
    fName = os.path.basename(file)
    try:
        featureFile = pd.read_csv(os.path.join(featureFilePath, fName), index_col=0)
    except Exception:
        print('File not Found')
        print(file)
        continue
    print(file)
    featureFile = featureFile.iloc[indexNos]
    featureFile[classifierName] = classifications
    featureFile.index.names = ['scorer']
    outPath = os.path.join(outputPath, fName)
    featureFile.to_csv(outPath)
    print('Saved ' + outPath)


