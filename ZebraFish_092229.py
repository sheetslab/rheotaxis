from __future__ import division
import os
import pandas as pd
import numpy as np
from configparser import ConfigParser, NoOptionError
import glob
import math
from numba import jit
from simba.rw_dfs import *
from simba.drop_bp_cords import *


inifile = r"Z:\DeepLabCut\DLC_extract\Troubleshooting\Zebrafish\Zebrafish\project_folder\project_config.ini"


def extract_features_userdef(inifile):
    windows_angular_dispersion_seconds = [10]
    config = ConfigParser()
    configFile = str(inifile)
    config.read(configFile)
    projectPath = config.get('General settings', 'project_path')
    csv_dir_in, csv_dir_out = os.path.join(projectPath, 'csv', 'outlier_corrected_movement_location'), os.path.join(projectPath, 'csv', 'features_extracted')
    vidInfPath = os.path.join(projectPath, 'logs', 'video_info.csv')
    vidinfDf = pd.read_csv(vidInfPath)
    poseConfigPath = os.path.join(os.path.join(projectPath, 'logs'), 'measures', 'pose_configs', 'bp_names', 'project_bp_names.csv')
    poseConfigDf = pd.read_csv(poseConfigPath, header=None)
    poseConfigDf = list(poseConfigDf[0])
    try:
        wfileType = config.get('General settings', 'workflow_file_type')
    except NoOptionError:
        wfileType = 'csv'

    def angle2pt_degrees(ax, ay, bx, by):
        angle_degrees = math.degrees(math.atan2(ax - bx, by - ay))
        return angle_degrees + 360 if angle_degrees < 0 else angle_degrees

    def angle2pt_radians(degrees):
        angle_radians = degrees * math.pi / 180
        return angle_radians

    def angle2pt_sin(angle_radians):
        angle_sin = math.sin(angle_radians)
        return angle_sin

    def angle2pt_cos(angle_radians):
        angle_cos = math.cos(angle_radians)
        return angle_cos

    @jit(nopython=True)
    def EuclidianDistCald(bp1xVals, bp1yVals, bp2xVals, bp2yVals):
        series = (np.sqrt((bp1xVals - bp2xVals) ** 2 + (bp1yVals - bp2yVals) ** 2))
        return series

    fileCounter = 0
    roll_windows_values = [10, 4, 2, 1, 0.1]

    filesFound = glob.glob(csv_dir_in + '/*.csv')
    print('Extracting features from ' + str(len(filesFound)) + ' files...')

    ########### CREATE PD FOR RAW DATA AND PD FOR MOVEMENT BETWEEN FRAMES ###########
    for currentFile in filesFound:
        roll_windows, angular_dispersion_windows = [], []
        currVidName = os.path.basename(currentFile.replace('.' + wfileType, ''))
        currVideoSettings = vidinfDf.loc[vidinfDf['Video'] == currVidName]
        try:
            currPixPerMM = float(currVideoSettings['pixels/mm'])
        except TypeError:
            print('Error: make sure all the videos that are going to be analyzed are represented in the project_folder/logs/video_info.csv file')
        fps = float(currVideoSettings['fps'])
        print('Processing ' + '"' + str(currVidName) + '".' + ' Fps: ' + str(fps) + ". mm/ppx: " + str(currPixPerMM))
        for i in range(len(roll_windows_values)):
            roll_windows.append(int(fps / roll_windows_values[i]))
        for i in range(len(windows_angular_dispersion_seconds)):
            angular_dispersion_windows.append(int(fps * windows_angular_dispersion_seconds[i]))
        bodypartNames = list(poseConfigDf)

        columnHeaders, columnHeadersShifted = [], []
        x_cols, y_cols, p_cols,x_cols_shifted, y_cols_shifted = [], [], [], [], []
        for bodypart in bodypartNames:
            colHead1, colHead2, colHead3 = (bodypart + '_x', bodypart + '_y', bodypart + '_p')
            colHead4, colHead5, colHead6 = (bodypart + '_x_shifted', bodypart + '_y_shifted', bodypart + '_p_shifted')
            columnHeaders.extend((colHead1, colHead2, colHead3))
            columnHeadersShifted.extend((colHead4, colHead5, colHead6))
            p_cols.append(colHead3)
            x_cols.append(colHead1)
            y_cols.append(colHead2)
            x_cols_shifted.append(colHead4)
            y_cols_shifted.append(colHead5)

        csv_df = read_df(currentFile, wfileType)
        csv_df.columns = columnHeaders
        csv_df = csv_df.fillna(0)
        csv_df = csv_df.apply(pd.to_numeric)

        # ########### CREATE SHIFTED DATAFRAME FOR DISTANCE CALCULATIONS ###########################################
        csv_df_shifted = csv_df.shift(periods=1)
        csv_df_shifted.columns = columnHeadersShifted
        csv_df_combined = pd.concat([csv_df, csv_df_shifted], axis=1, join='inner')
        csv_df_combined = csv_df_combined.fillna(0)

        print('Calculating single frame movements...')
        headerList = []
        ################## CALUCLATE "X relative to Y" movements ##################
        for bp in range(len(x_cols)):
            currX, currShiftedX, currY, currShiftedY = x_cols[bp], x_cols_shifted[bp], y_cols[bp], y_cols_shifted[bp]
            currXcolMoveName, currYcolMoveName, currXrelative2Y = str('Movement_bp_' + str(bp) + '_axis_x'), str('Movement_bp_' + str(bp) + '_axis_y'), str('Movement_bp_' + str(bp) + '_X_relative_2_Y')
            headerList.append(currXrelative2Y)
            csv_df_combined[currXcolMoveName] = (csv_df_combined[currX] - csv_df_combined[currShiftedX])
            csv_df_combined[currYcolMoveName] = (csv_df_combined[currY] - csv_df_combined[currShiftedY])
            csv_df_combined[currXrelative2Y] = (csv_df_combined[currXcolMoveName] - csv_df_combined[currYcolMoveName])
        csv_df_combined['All_bps_X_relative_2_y'] = csv_df_combined[headerList].sum(axis=1)
        csv_df_combined.drop(headerList, axis=1)

        ################## Calculate tail movements ##################
        csv_df_combined['Tail_1_movement'] = EuclidianDistCald(csv_df_combined['Zebrafish_Tail1_x'].values, csv_df_combined['Zebrafish_Tail1_y'].values, csv_df_combined['Zebrafish_Tail1_x_shifted'].values, csv_df_combined['Zebrafish_Tail1_y_shifted'].values) / currPixPerMM
        csv_df_combined['Tail_2_movement'] = EuclidianDistCald(csv_df_combined['Zebrafish_Tail2_x'].values, csv_df_combined['Zebrafish_Tail2_y'].values, csv_df_combined['Zebrafish_Tail2_x_shifted'].values, csv_df_combined['Zebrafish_Tail2_y_shifted'].values) / currPixPerMM
        csv_df_combined['Tail_3_movement'] = EuclidianDistCald(csv_df_combined['Zebrafish_Tail3_x'].values, csv_df_combined['Zebrafish_Tail3_y'].values, csv_df_combined['Zebrafish_Tail3_x_shifted'].values, csv_df_combined['Zebrafish_Tail3_y_shifted'].values) / currPixPerMM
        csv_df_combined['Tail_4_movement'] = EuclidianDistCald(csv_df_combined['Zebrafish_Tail4_x'].values, csv_df_combined['Zebrafish_Tail4_y'].values, csv_df_combined['Zebrafish_Tail4_x_shifted'].values, csv_df_combined['Zebrafish_Tail4_y_shifted'].values) / currPixPerMM
        csv_df_combined['SwimBladder_movement'] = EuclidianDistCald(csv_df_combined['Zebrafish_SwimBladder_x'].values, csv_df_combined['Zebrafish_SwimBladder_y'].values, csv_df_combined['Zebrafish_SwimBladder_x_shifted'].values, csv_df_combined['Zebrafish_SwimBladder_y_shifted'].values) / currPixPerMM
        csv_df_combined.loc[0, 'SwimBladder_movement'] = 0
        csv_df_combined['SwimBladder_distance travelled'] = csv_df_combined['SwimBladder_movement'].cumsum()

        csv_df_combined['Total_tail_movement'] = csv_df_combined['Tail_1_movement'] + csv_df_combined['Tail_2_movement'] + csv_df_combined['Tail_3_movement'] + csv_df_combined['Tail_4_movement']
        csv_df_combined['Lower_tail_movement'] = csv_df_combined['Tail_3_movement'] + csv_df_combined['Tail_4_movement']

        ########### CALC ROLLING WINDOWS MEANS AND SUMS ###########################################
        print('Calculating rolling windows: medians and sums X, Y movements...')
        for i in range(len(roll_windows_values)):
            currentColName = 'All_bps_X_relative_2_y_mean_' + str(roll_windows_values[i])
            csv_df_combined[currentColName] = csv_df_combined['All_bps_X_relative_2_y'].rolling(roll_windows[i], min_periods=1).mean()
            currentColName = 'All_bps_X_relative_2_y_sum_' + str(roll_windows_values[i])
            csv_df_combined[currentColName] = csv_df_combined['All_bps_X_relative_2_y'].rolling(roll_windows[i], min_periods=1).sum()

        csv_df_combined['SwimBladder_velocity'] = csv_df_combined['SwimBladder_movement'].rolling(int(fps), min_periods=1).sum()
        csv_df_combined['Fish_clockwise_angle_degrees'] = csv_df_combined.apply(lambda x: angle2pt_degrees(x['Zebrafish_SwimBladder_x'], x['Zebrafish_SwimBladder_y'], x['Zebrafish_Tail3_x'],x['Zebrafish_Tail3_y']), axis=1)
        csv_df_combined['Fish_angle_radians'] = angle2pt_radians(csv_df_combined['Fish_clockwise_angle_degrees'])
        csv_df_combined['Fish_angle_sin'] = csv_df_combined.apply(lambda x: angle2pt_sin(x['Fish_angle_radians']), axis=1)
        csv_df_combined['Fish_angle_cos'] = csv_df_combined.apply(lambda x: angle2pt_cos(x['Fish_angle_radians']), axis=1)
        csv_df_combined['Fish_angle_cumsum_sin'] = csv_df_combined['Fish_angle_sin'].cumsum()
        csv_df_combined['Fish_angle_cumsum_cos'] = csv_df_combined['Fish_angle_cos'].cumsum()

        compass_brackets = ["N", "NE", "E", "SE", "S", "SW", "W", "NW", "N"]
        compass_brackets_digits = ["0", "1", "2", "3", "4", "5", "6", "7", "0"]
        compass_lookup = list(round(csv_df_combined['Fish_clockwise_angle_degrees'] / 45))
        compass_lookup = [int(i) for i in compass_lookup]
        compasFaceList_bracket, compasFaceList_digit = [], []
        for compasDirection in compass_lookup:
            compasFaceList_bracket.append(compass_brackets[compasDirection])
            compasFaceList_digit.append(compass_brackets_digits[compasDirection])
        csv_df_combined['Compass_Direction'] = compasFaceList_bracket
        csv_df_combined['Compass_Digit'] = compasFaceList_digit


        for i in range(len(roll_windows_values)):
            currentColName = 'Mean_angle_time_window_' + str(roll_windows_values[i])
            csv_df_combined[currentColName] = csv_df_combined['Fish_clockwise_angle_degrees'].rolling(roll_windows[i], min_periods=1).mean()


        ########### CALC ROLLING WINDOWS: DIRECTION CHANGES ###########################################
        groupDf = pd.DataFrame()
        v = (csv_df_combined['Compass_Digit'] != csv_df_combined['Compass_Digit'].shift()).cumsum()
        u = csv_df_combined.groupby(v)['Compass_Digit'].agg(['all', 'count'])
        m = u['all'] & u['count'].ge(1)
        groupDf['groups'] = csv_df_combined.groupby(v).apply(lambda x: (x.index[0], x.index[-1]))[m]
        currdirectionList, DirectionSwitchIndexList, currdirectionListValue = [], [], []
        for indexes, row in groupDf.iterrows():
            currdirectionList.append(csv_df_combined.loc[row['groups'][0]]['Compass_Direction'])
            DirectionSwitchIndexList.append(row['groups'][1])
            currdirectionListValue.append(csv_df_combined.loc[row['groups'][0]]['Compass_Digit'])
        groupDf['Direction_switch'] = currdirectionList
        groupDf['Direction_value'] = currdirectionListValue
        for i in DirectionSwitchIndexList:
            csv_df_combined.at[i, 'DirectionSwitch'] = 1
        csv_df_combined['Compass_Digit_shifted'] = csv_df_combined.Compass_Digit.shift(-1)
        csv_df_combined = csv_df_combined.dropna(axis=0, how='all')
        csv_df_combined = csv_df_combined.fillna(0)
        switch_directionList = []
        for index, row in csv_df_combined.iterrows():
            if (row['Compass_Digit_shifted'] == '0') and ((row['Compass_Digit'] == '7')):
                switchDirection = 1
            else:
                switchDirection = int(row['Compass_Digit_shifted']) - int(row['Compass_Digit'])
            switch_directionList.append(switchDirection)
        csv_df_combined['SwitchDirectionValue'] = switch_directionList

        compassHotEnd = pd.get_dummies(csv_df_combined.Compass_Direction, prefix='Direction')
        compass_brackets = ["Direction_N", "Direction_NE", "Direction_E", "Direction_SE", "Direction_S", "Direction_SW", "Direction_W", "Direction_NW"]
        compassHotEnd = compassHotEnd.T.reindex(compass_brackets).T.fillna(0)
        csv_df_combined = pd.concat([csv_df_combined,compassHotEnd], axis=1)

        ########### CALC ROLLING WINDOWS MEANS AND SUMS ###########################################
        print('Calculating rolling windows: sums direction switches...')
        for i in range(len(roll_windows_values)):
            currentColName = 'Number_of_direction_switches_' + str(roll_windows_values[i])
            csv_df_combined[currentColName] = csv_df_combined['DirectionSwitch'].rolling(roll_windows[i], min_periods=1).sum()
            currentColName = 'Directionality_of_switches_switches_' + str(roll_windows_values[i])
            csv_df_combined[currentColName] = csv_df_combined['SwitchDirectionValue'].rolling(roll_windows[i], min_periods=1).sum()

        ############# CALCULATING ANGULAR DISPERSION ###############
        angularDispertionList = []
        for index, row in csv_df_combined.iterrows():
            X, Y = row['Fish_angle_cumsum_cos'] / (index + 1),  row['Fish_angle_cumsum_sin'] / (index + 1)
            angularDispertionList.append(math.sqrt(X**2 + Y**2))
        csv_df_combined['Angular_dispersion'] = angularDispertionList
        for window in range(len(angular_dispersion_windows)):
            currentColName = 'Angular_dispersion_window_' + str(windows_angular_dispersion_seconds[window])
            csv_df_combined[currentColName] = csv_df_combined['Angular_dispersion'].rolling(angular_dispersion_windows[window], min_periods=1).mean()
        csv_df_combined = csv_df_combined.drop(columnHeadersShifted, axis=1)
        csv_df_combined = csv_df_combined.drop(['Compass_Digit_shifted', 'DirectionSwitch', 'SwitchDirectionValue','Compass_Digit', 'Compass_Direction', 'Fish_angle_cumsum_sin', 'Fish_angle_cumsum_cos'], axis=1)
        csv_df_combined.to_csv(os.path.join(csv_dir_out, currVidName + '.' + wfileType))
        print('Features extracted for ' + currVidName + '...')
    print('All features extracted for files in project.')