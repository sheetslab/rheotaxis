import numpy as np
import cv2
import os
import pandas as pd
from scipy import ndimage
from configparser import ConfigParser, MissingSectionHeaderError, NoSectionError
import glob
from simba.drop_bp_cords import getBpNames
from pylab import *
import random

configini = r"C:\Users\Fish_Behavior\Desktop\SIMBA\1LZF_DLCma_model\project_folder\project_config.ini"
frameSetting = 2     #### 0 = Create frames and video, 1 = Creates frames,  2 = Creates video

config = ConfigParser()
configFile = str(configini)
try:
    config.read(configFile)
except MissingSectionHeaderError:
    print('ERROR:  Not a valid project_config file. Please check the project_config.ini path.')
animalsNo = config.getint('General settings', 'animal_no')
projectPath = config.get('General settings', 'project_path')
frames_dir_out = config.get('Frame settings', 'frames_dir_out')
csv_dir_in = os.path.join(projectPath, 'csv', 'features_extracted')
frames_dir_out = os.path.join(frames_dir_out, 'sklearn_results')
poseConfSetting = config.get('create ensemble settings', 'pose_estimation_body_parts')
if not os.path.exists(frames_dir_out):
    os.makedirs(frames_dir_out)
vidInfPath = os.path.join(projectPath, 'logs', 'video_info.csv')
mulltiAnimalIDList= config.get('Multi animal IDs', 'id_list')
mulltiAnimalIDList = mulltiAnimalIDList.split(",")
if mulltiAnimalIDList[0] != '':
    mulltiAnimalStatus = True
if mulltiAnimalIDList[0] == '':
    mulltiAnimalStatus = False
vidinfDf = pd.read_csv(vidInfPath)
target_names, colorList_animal_1, colorList_animal_2, loopy = [], [], [], 0
Xcols, Ycols, Pcols = getBpNames(configini)
cmap = cm.get_cmap('hot', len(Xcols) + 1)
for i in range(cmap.N):
    rgb = list((cmap(i)[:3]))
    rgb = [i * 255 for i in rgb]
    rgb.reverse()
    colorList_animal_1.append(rgb)
filesFound = glob.glob(csv_dir_in + '/*.csv')
print('Processing ' + str(len(filesFound)) + ' videos ...')

for i in filesFound:
    currentVideo = i
    loopy += 1
    CurrentVideoName = os.path.basename(currentVideo)
    if frameSetting == 1:
        videoFrameDir = os.path.join(frames_dir_out, CurrentVideoName.replace('.csv', ''))
        if not os.path.exists(videoFrameDir):
            os.makedirs(videoFrameDir)
    CurrentVideoRow = vidinfDf.loc[vidinfDf['Video'] == str(CurrentVideoName.replace('.csv', ''))]
    try:
        fps = int(CurrentVideoRow['fps'])
    except TypeError:
        print('Error: make sure all the videos that are going to be analyzed are represented in the project_folder/logs/video_info.csv file')
    currentDf = pd.read_csv(currentVideo)
    currentDf = currentDf.fillna(0)
    currentDf = currentDf.loc[:, ~currentDf.columns.str.contains('^Unnamed')]
    currentDf = currentDf.reset_index()
    animalBpHeaderList, animalBpHeaderListY, animalBpHeaderListX = ([], [], [])
    animal1_BPsX, animal1_BPsY = (currentDf[Xcols], currentDf[Ycols])
    for i in range(len(animal1_BPsX.columns)):
        animalBpHeaderListX.append(animal1_BPsX.columns[i])
        animalBpHeaderListY.append(animal1_BPsY.columns[i])
        animalBpHeaderList.append(animal1_BPsX.columns[i])
        animalBpHeaderList.append(animal1_BPsY.columns[i])
    animalBpHeaderListX, animalBpHeaderListY, animalBpHeaderList = ([x for x in animalBpHeaderListX if "Tail_end" not in x], [x for x in animalBpHeaderListY if "Tail_end" not in x], [x for x in animalBpHeaderList if "Tail_end" not in x])
    if os.path.exists(os.path.join(projectPath,'videos', CurrentVideoName.replace('.csv', '.mp4'))):
        videoPathName = os.path.join(projectPath,'videos', CurrentVideoName.replace('.csv', '.mp4'))
    elif os.path.exists(os.path.join(projectPath,'videos', CurrentVideoName.replace('.csv', '.avi'))):
        videoPathName = os.path.join(projectPath,'videos', CurrentVideoName.replace('.csv', '.avi'))
    cap = cv2.VideoCapture(videoPathName)
    width, height, frames = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH)), int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT)), int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    outputFileName = os.path.join(frames_dir_out, CurrentVideoName)
    if height < width:
        videoHeight, videoWidth = width, height
    if height >= width:
        videoHeight, videoWidth = height, width
    recWidth = videoWidth + int(videoWidth / 2)
    writer = cv2.VideoWriter(outputFileName.replace('.csv', '.avi'), fourcc, fps, (recWidth, videoHeight))
    currRow = 0
    mySpaceScale, myRadius, myResolution, myFontScale = 60, 12, 1500, 1.2
    maxResDimension = max(width, height)
    circleScale, fontScale, spacingScale = int(myRadius / (myResolution / maxResDimension)), float(myFontScale / (myResolution / maxResDimension)), int(mySpaceScale / (myResolution / maxResDimension))
    while (cap.isOpened()):
        ret, frame = cap.read()
        IDlabelLoc, rotationFlag = [], False
        sideImg = np.zeros((videoHeight, int(videoWidth / 2), 3))
        if ret == True:
            if (animalsNo == 1) and (poseConfSetting != 'user_defined'):
                currAnimal1 = currentDf.loc[currentDf.index[currRow], animalBpHeaderList]
                currAnimal1 = np.array(currAnimal1).astype(int)
                currAnimal1 = np.reshape(currAnimal1, (-1, 2))
                M1polyglon_array_hull = cv2.convexHull((currAnimal1.astype(int)))
                cv2.drawContours(frame, [M1polyglon_array_hull.astype(int)], 0, (255, 255, 255), 2)
             for cords in range(len(animalBpHeaderListX)):
                 currXval = animal1_BPsX.loc[animal1_BPsX.index[currRow], animalBpHeaderListX[cords]]
                 currYval = animal1_BPsY.loc[animal1_BPsY.index[currRow], animalBpHeaderListY[cords]]
                 if animalBpHeaderListX[cords] in animalBpHeaderList:
                     color = colorList_animal_1[cords]
                 cv2.circle(frame, (int(currXval), int(currYval)), circleScale, color, -1, lineType=cv2.LINE_AA)

            directionsList = ['East', 'North', 'North-east', 'North-west', 'South', 'South-east', 'South-west', 'West']
            clockwiseRot = currentDf.loc[currentDf.index[currRow], 'Fish_clockwise_angle_degrees']
            directionList = [currentDf.loc[currentDf.index[currRow], 'Direction_E'], currentDf.loc[currentDf.index[currRow], 'Direction_N'], currentDf.loc[currentDf.index[currRow], 'Direction_NE'], currentDf.loc[currentDf.index[currRow], 'Direction_NW'], currentDf.loc[currentDf.index[currRow], 'Direction_S'], currentDf.loc[currentDf.index[currRow], 'Direction_SE'], currentDf.loc[currentDf.index[currRow], 'Direction_SW'], currentDf.loc[currentDf.index[currRow], 'Direction_SW']]
            for direction in range(len(directionList)):
                currval = directionList[direction]
                if currval == 1:
                    currDirection = directionsList[direction]
            movementXrelY = currentDf.loc[currentDf.index[currRow], 'All_bps_X_relative_2_y_sum_15']
            tailMovement = currentDf.loc[currentDf.index[currRow], 'Lower_tail_movement']
            currentAngleDegrees = round(currentDf.loc[currentDf.index[currRow], 'Fish_clockwise_angle_degrees'], 1)
            currentAngleRadians = round(currentDf.loc[currentDf.index[currRow], 'Fish_angle_radians'], 1)
            currCompN, currCompNE, currCompE, currCompSE, currCompS, currCompSW, currCompW, currCompNW  = currentDf.loc[currentDf.index[currRow], 'Direction_N'], currentDf.loc[currentDf.index[currRow], 'Direction_NE'], currentDf.loc[currentDf.index[currRow], 'Direction_E'], currentDf.loc[currentDf.index[currRow], 'Direction_SE'], currentDf.loc[currentDf.index[currRow], 'Direction_S'], currentDf.loc[currentDf.index[currRow], 'Direction_SW'], currentDf.loc[currentDf.index[currRow], 'Direction_W'], currentDf.loc[currentDf.index[currRow], 'Direction_NW']
            currentVelocity = round(currentDf.loc[currentDf.index[currRow], 'SwimBladder_velocity'], 1)
            dispersion = round(currentDf.loc[currentDf.index[currRow], 'Angular_dispersion'], 3)
            totalDistanceTravelled = round(currentDf.loc[currentDf.index[currRow], 'SwimBladder_distance travelled'], 1)
            cv2.putText(sideImg, 'Current angle (deg): ' + str(currentAngleDegrees), (10, 60), cv2.FONT_HERSHEY_SIMPLEX, fontScale, (0, 255, 0), 2)
            cv2.putText(sideImg,'Velocity (mm / s): ' + str(currentVelocity), (10, 90), cv2.FONT_HERSHEY_SIMPLEX, fontScale, (0, 255, 0), 2)
            cv2.putText(sideImg, 'Distance travelled (mm): ' + str(totalDistanceTravelled), (10, 120), cv2.FONT_HERSHEY_SIMPLEX, fontScale, (0, 255, 0), 2)
            if (currCompN == 1) or (currCompNE == 1) or (currCompNW == 1):
                rheotaxisStatus = 'Y'
            else:
                rheotaxisStatus = 'N'
            cv2.putText(sideImg, 'Rheotaxis: ' + str(rheotaxisStatus), (10, 150), cv2.FONT_HERSHEY_SIMPLEX, fontScale, (0, 255, 0), 2, lineType=cv2.LINE_AA)
            #cv2.putText(sideImg, 'Current angular dispersion: ' + str(dispersion), (10, 180), cv2.FONT_HERSHEY_SIMPLEX, fontScale, (0, 255, 0), 2)
            outImage = np.concatenate((sideImg, frame), axis=1)
            outImage = outImage.astype(np.uint8)


            #
            #
            # cv2.imshow('frame', outImage)
            # if cv2.waitKey(2) & 0xFF == ord('q'):
            #    break
        if (frameSetting == 0) or (frameSetting == 1):
            frameName = os.path.join(videoFrameDir, str(currRow) + '.png')
            cv2.imwrite(frameName, outImage)
        if (frameSetting == 0) or (frameSetting == 2):
            writer.write(outImage)
        currRow += 1
        print('Processed frame ' + str(currRow) + ' / ' + str(len(currentDf)))
        if (frame is None) and (frameSetting == 0 or 2):
            print('Video ' + str(os.path.basename(CurrentVideoName.replace('.csv', '.mp4'))) + ' saved.')
            cap.release()
            break
        if (frame is None) and (frameSetting == 1):
            print('All images saved')