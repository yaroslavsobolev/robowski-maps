"""Split the csv file into samller csv files, each containing the data for one well plate."""
import os
import PySimpleGUI as sg

# use pysimplegui to get the path of a csv file
csv_path = sg.popup_get_file('file path', no_window=True)
csv_folder_path = os.path.dirname(os.path.abspath(csv_path))

csvfile = open(csv_path, 'r').readlines()
header = csvfile[0]
csvfile = csvfile[1:]

plate_index = 0

for i in range(len(csvfile)):
    if (i+1) % 54 == 0:
        with open(csv_folder_path + f'\\outV_split_by_plate\\plate_{plate_index}_{i-53}-{i}.csv', 'w+') as f:
            f.writelines(header)
            f.writelines(csvfile[i-53 : i+1])
            print(f'plate_{plate_index}_{i-53}-{i}.csv')
        plate_index += 1