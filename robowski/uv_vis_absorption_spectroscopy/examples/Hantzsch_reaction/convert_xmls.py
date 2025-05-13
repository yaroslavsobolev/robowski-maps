import importlib
import glob
import os

xml_converter = importlib.import_module('zeus-pipetter.misc.nanodrop_xml_processing')
data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

# xml_folder = f'{data_folder}BPRF/calibrations/'
# xml_name = '2024_04_0_UV-Vis_bb021_F2.xml'
# print("Working on the following file: ")
# print(f'xml_folder: {xml_folder}')
# print(f'xml_name: {xml_name}')
# if "2023" not in xml_name:
#     df = xml_converter.treat_one_file(xml_folder=xml_folder,
#                         xml_name=xml_name)
# else:
#     print('This file has been processed already!')

xml_folder = f'{data_folder}Yaroslav/'
xml_name = '2024_04_08_UV-Vis__nilered_congored_1.xml'
print("Working on the following file: ")
print(f'xml_folder: {xml_folder}')
print(f'xml_name: {xml_name}')
if "2023" not in xml_name:
    df = xml_converter.treat_one_file(xml_folder=xml_folder,
                        xml_name=xml_name)
else:
    print('This file has been processed already!')