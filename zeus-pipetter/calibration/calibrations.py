import os
import re
import time

import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import json

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
robot_name = os.environ['PIPETTER_NAME']

def treat_console_txt( file ):

    # this will convert the console txt to the csv file needed for calib."

    folder = '\\'.join(file.split('\\')[:-1]) + '\\'
    csv_name = file.split('\\')[-1].split('.')[0] + '.csv'

    t = []
    # read txt file
    with open(file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if 'INFO:main:Result:' in line:
                t.append(line)

    output = {}
    for line in t:
        vol = int(line.split('ethanol_')[-1].split("ul': {'weight':")[0])
        vol_lst = line.split("'volume': ")[-1].split(", 'liq")[0]
        vol_lst = vol_lst.replace('[', '')
        vol_lst = vol_lst.replace(']', '')
        vol_lst = vol_lst.replace(' ', '')
        print(vol, ',', vol_lst)
        output[vol] = [round(float(i) / 1000,2) for i in vol_lst.split(',')]

    # convert the dictionary to dataframe
    df = pd.DataFrame(output, dtype=object)
    df = df.T
    # set the first column as index
    df.set_index(df.iloc[:, 0])

    print(f'df here is {df}')
    # get number of columns
    num_cols = len(df.columns)
    df.columns = [f'm{i}' for i in range(num_cols)]

    # save the df to csv
    df.to_csv(folder + csv_name)

    # calculate mean
    df['mean'] = round(df.iloc[:, 1:].mean(axis=1),2)

    # get time string
    time_here = time.strftime("%Y-%m-%d-%H-%M-%S", time.localtime())
    # save results_all.csv
    df.to_csv(folder + f'results_all_{time_here}.csv')

    return df

class Interpolation:
    def __init__(self, solvent_name):
        self.solvent_name = solvent_name
        print('##Calibration for ', solvent_name, ' is initialized.')
        if robot_name == 'Robowski#1':
            self.root_folder = 'C:\\Users\\Chemiluminescence\Dropbox\\robochem\\data\\pipetter_calib\\'
            self.calib_file_dict = {'DCE':(self.root_folder+'2024-03-05-calib-DCE\\', 'results_all.csv'),
                                    'Dioxane':(self.root_folder+'2024-03-02-calib-dioxane\\', 'results_all.csv'),
                                    'ethanol':(self.root_folder+'2024-03-05-calib-ethanol\\', 'results_all.csv'),
                                    'MeCN':(self.root_folder+'2024-03-07-calib-MeCN\\', 'results_all.csv'),
                                    'DMF':(self.root_folder+'2024-03-09-calib-DMF\\', 'results_all.csv'),
                                    'hbrhac1v1':()}
        elif robot_name == 'Robowski#2':
            self.root_folder = data_folder+ '\\pipetter_calib\\robowski2'
            self.calib_file_dict = {'DCE':(),
                                    'Dioxane':(),
                                    'ethanol':(self.root_folder+'\\ethanol\\', 'results_all.csv'),
                                    'MeCN':(),
                                    'DMF':(),
                                    'hbrhac1v1':()}
    def import_csv(self):
        # print(f'self.sovent_name is : {self.solvent_name}.')
        csv_file = self.calib_file_dict[self.solvent_name][0] + self.calib_file_dict[self.solvent_name][1]
        df_calib = pd.read_csv(csv_file)

        # set the first column as index
        df_calib.set_index(df_calib.columns[0], inplace=True)

        # calculate mean
        if 'mean' not in df_calib.columns:
            df_calib['mean'] = df_calib.mean(numeric_only=True, axis=1)

        # save this new df to csv by overwriting the old one
        df_calib.to_csv(csv_file, index=True)

        return df_calib

    def interp(self, target_volume, purpose = 'mixing', verbose = True, correction = True):

        df_calib = self.import_csv()

        df_calib_50 = df_calib.loc[(df_calib.index >= 5) & (df_calib.index <= 50), :]
        df_calib_300 = df_calib.loc[(df_calib.index > 50) & (df_calib.index <= 300), :]
        df_calib_1000 = df_calib.loc[(df_calib.index > 300) & (df_calib.index <= 1000), :]

        set_volumes_50 = df_calib_50.index.tolist()
        set_volumes_300 = df_calib_300.index.tolist()
        set_volumes_1000 = df_calib_1000.index.tolist()
        actual_volumes_50 = df_calib_50['mean'].tolist()
        actual_volumes_300 = df_calib_300['mean'].tolist()
        actual_volumes_1000 = df_calib_1000['mean'].tolist()

        # print(f'set_volumes_50: {set_volumes_50}')
        # print(f'actual_columes_50: {actual_volumes_50}')

        interpolation_func_50 = interp1d(actual_volumes_50, set_volumes_50, kind='linear', fill_value="extrapolate")
        setting_volume_50 = round(float(interpolation_func_50(target_volume)), 2)

        interpolation_func_300 = interp1d(actual_volumes_300, set_volumes_300, kind='linear', fill_value="extrapolate")
        setting_volume_300 = round(float(interpolation_func_300(target_volume)), 2)

        interpolation_func_1000 = interp1d(actual_volumes_1000, set_volumes_1000, kind='linear', fill_value="extrapolate")
        setting_volume_1000 = round(float(interpolation_func_1000(target_volume)), 2)

        # print(setting_volume_50, setting_volume_300, setting_volume_1000)

        # # plot the mean of df_calib_50, df_calib_300, df_calib_1000 against the index
        # import matplotlib.pyplot as plt
        # import matplotlib
        # matplotlib.use('TkAgg')
        # plt.plot(df_calib_50.index, df_calib_50['mean_0128'], label='50', marker='o')
        # plt.plot(df_calib_300.index, df_calib_300['mean_0128'], label='300', marker='o')
        # plt.plot(df_calib_1000.index, df_calib_1000['mean_0128'], label='1000', marker='o')
        # plt.legend()
        # plt.show()

        if target_volume <= 50:
            if setting_volume_50 <= 50:
                if verbose: print(
                    f'{target_volume} ul is calibrated by setting_volume_50 for @@{purpose}, set to {setting_volume_50}.')
                return setting_volume_50
            elif setting_volume_50 <= 300:
                if verbose: print(
                    f'{target_volume} ul is calibrated by setting_volume_300 for @@{purpose}, set to {setting_volume_300}')
                return setting_volume_300
        elif target_volume <= 300:
            if setting_volume_300 <= 50:
                if verbose: print(
                    f'{target_volume} ul is calibrated by setting_volume_50 for @@{purpose}, set to {setting_volume_50}')
                return setting_volume_50
            elif setting_volume_300 <= 300:
                if verbose: print(
                    f'{target_volume} ul is calibrated by setting_volume_300 for @@{purpose}, set to {setting_volume_300}.')
                return setting_volume_300
            elif setting_volume_300 < 1000:
                if verbose: print(
                    f'{target_volume} ul is calibrated by setting_volume_1000 for @@{purpose}, set to {setting_volume_1000}.')
                return setting_volume_1000
        elif target_volume <= 1000:
            if setting_volume_1000 <= 300:
                if verbose: print(
                    f'{target_volume} ul is calibrated by setting_volume_300 for @@{purpose}, set to {setting_volume_300}.')
                return setting_volume_300
            elif setting_volume_1000 <= 1000:
                if verbose: print(
                    f'{target_volume} ul is calibrated by setting_volume_1000 for @@{purpose}, set to {setting_volume_1000}.')
                return setting_volume_1000

            elif target_volume <= 1030:
                if verbose: print(f'{target_volume} ul is calibrated by setting_volume_1000 for @@{purpose}, set to 1000ul.')

                return 1000
            else:
                raise ValueError(f'@@{purpose} target_volume is out of range.')

    def treat_json(self):
        json_file = self.calib_file_dict[self.solvent_name][0] + 'results_after_calib.json'
        # read json file into a dictionary
        with open(json_file, 'r') as f:
            data = json.load(f)

        data_dict = {}

        for key in data.keys():
           for dict in data[key]:
               for key, value in dict.items():
                   # get the digit from the key using regex
                   num = re.findall(r'\d+', key)
                   num = int(''.join([str(i) for i in num]))
                   data_dict[num]=value['volume']

        print(data_dict)

        for key, value in data_dict.items():
            list = data_dict[key]
            num_of_tests = len(list)
            break

        dict_df = {}

        for i in range(num_of_tests):
            dict_df[f'm_{i}'] = [] if f'm_{i}' not in dict_df.keys() else dict_df[f'm_{i}']
            for key, value in data_dict.items():
                dict_df[f'm_{i}'].append(value[i])

        # save the dictionary into a dataframe
        df_after_calib = pd.DataFrame(dict_df, dtype=object)
        # set the keys of data_list as index
        df_after_calib.index= [key for key, values in data_dict.items()]

        # delete the columns ['mean', 'std', 'error_%'] if they exist
        if 'mean' in df_after_calib.columns:
            del df_after_calib['mean']
        if 'std' in df_after_calib.columns:
            del df_after_calib['std']
        if 'error_%' in df_after_calib.columns:
            del df_after_calib['error_%']

        # calculate mean
        if 'mean' not in df_after_calib.columns:
            df_after_calib['mean'] = df_after_calib.mean(numeric_only=True, axis=1)
        # calculate std
        if 'std' not in df_after_calib.columns:
            df_after_calib['std'] = df_after_calib.std(numeric_only=True, axis=1)

        # calculate the error by dividing the offest by the mean
        if 'error_%' not in df_after_calib.columns:
            df_after_calib['error_%'] = np.abs((df_after_calib['mean'] - df_after_calib.index) / df_after_calib.index * 100)

        # save the df to csv
        with open(self.calib_file_dict[self.solvent_name][0] + 'results_from_json.csv', 'w') as f:
            df_after_calib.to_csv(f, index=True)


        return df_after_calib

    def cal_after_calib(self):
        file = self.calib_file_dict[self.solvent_name][0] + 'results_all_after_calib.csv'
        df = pd.read_csv(file)
        # set the first column as index
        df.set_index(df.columns[0], inplace=True)

        for column in df.columns:
            if column not in ['m0', 'm1', 'm2', 'm3', 'm4']:
                del df[column]
        # calculate mean
        if 'mean' not in df.columns:
            df['mean'] = df.mean(numeric_only=True, axis=1)
        # calculate std
        if 'std' not in df.columns:
            df['std'] = df.std(numeric_only=True, axis=1)
        # calculate the error by dividing the offest by the mean
        if 'error_%' not in df.columns:
            df['error_%'] = np.abs((df['mean'] - df.index) / df.index * 100)

        # save this new df to csv by overwriting the old one
        df.to_csv(file, index=True)

        return df

if __name__ == '__main__':

    file = 'C:\\Yankai\\Dropbox\\robochem\\data\\pipetter_calib\\robowski2\\ethanol\\console_text_20240426_after_calib.txt'
    df = treat_console_txt(file)

    calib = Interpolation('ethanol')
    calib.interp(20)
