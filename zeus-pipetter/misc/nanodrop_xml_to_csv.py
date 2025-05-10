import sys

import pandas as pd, os, glob, numpy as np, time, datetime
from scipy.signal import savgol_filter
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from labellines import labelLine, labelLines


def remove_unwanted_columns(df):
    for column in df.columns:
        if 'remove' in column:
            print(column)
            try:
                df = df.drop(columns=[column])
            except:
                print(f"column {column} is not in df or it is already removed")
    return df


def treat_one_file(xml_folder, xml_name, rename = True):

    global df
    xml_path = xml_folder + xml_name
    tree = ET.parse(xml_path)
    workbook = tree.getroot()
    all_spectra, row_index, time_stamps = [], [], []
    # parse xml tree
    for worksheet in workbook[1:]:  ## root = workbook. workbook[0] is Styles
        table = list(worksheet)[0]  # worksheet has one child: table
        rows = list(table)  # table[0], table[1], table[2] are empty, table[3] is column data
        row_dict = {}
        for row in rows[5:]:
            # print(row.tag, row.attrib)
            row_data = list(row)
            wavelength = int(row_data[0][0].text)
            try:
                absorbance = float(row_data[1][0].text)
                row_dict[wavelength] = absorbance
            except ValueError:
                print("The format of the xml is incorrect! Please double check.")

            title = row_data[2][0].text
            if title:
                row_index.append(title)
            if row_data[3][0].text:
                time_here = row_data[3][0].text
                time_stamps.append(time_here)

        all_spectra.append(row_dict)

    assert len(all_spectra) == len(row_index), "wrong data structure of XML file!!"

    # save row_dict to dataframe
    df = pd.DataFrame(all_spectra, index=row_index, dtype=object)
    df = df.T

    # remove unwated columens
    df = remove_unwanted_columns(df=df)

    if xml_name == '2023_08_27_UV-Vis_reference.xml':
        aa = ['dioxane-1', 'dioxane-2',
              'pro3000x', 'pro2000x', 'pro1000x',
              'pro800x', 'pro400x', 'pro200x', 'pro100x',
              'sub3000x', 'sub2000x', 'sub1000x',
              'sub800x', 'sub400x', 'sub200x', 'sub100x']
        df = df.reindex(columns=aa)

    if '2023-09-06-run01' in xml_path:
        if 'plate_50_16c' in xml_path:
            pass
        if 'plate_51_16c' in xml_path:
            pass
        if 'plate_52_16c' in xml_path:
            df.columns.values[64] = '100'  # this will be later deleted
            df.drop(['100'], axis=1, inplace=True)  # it is done this way for workaround the duplicating column names
            df = df.loc[:, ~df.columns.duplicated()]

    if '2023-09-07-run01' in xml_path:
        if 'plate_74_36c' in xml_path:
            pass
        if 'plate_75_36c' in xml_path:
            pass
        if 'plate_76_36c' in xml_path:
            pass

    if 'reference_for_E1' in xml_path:
        if '54times' in xml_path:
            # renaming the columns from 0 through 53
            df.columns = [i for i in range(54)]

    # rename xml and csv file
    if rename:
        # change xml file name

        # convert time stamp to unix time
        unix_time = time.mktime(datetime.datetime.strptime(time_stamps[0], "%m/%d/%Y %I:%M:%S %p").timetuple())
        # convert unix time to human-readable time
        human_readable_time = datetime.datetime.fromtimestamp(unix_time).strftime('%Y-%m-%d_%H-%M-%S')
        # rename the xml file name using the human-readable time
        os.rename(xml_folder + xml_name, human_readable_time + '_' + xml_name)

        column_names = df.columns.values
        uniques, counts = np.unique(column_names, return_counts=True)

        if len(column_names) == len(uniques):
            ## save csv file
            xml_name_new = human_readable_time + '_' + xml_name
            df.to_csv(xml_folder + xml_name_new[:-4] + '.csv')
        else:
            for unique, count in zip(uniques, counts):
                if count > 1:
                    print(f"There is a duplicating column\n column name: {unique}, count: {count}")

            print("There is at least one duplicating column." \
                  f"len of columns_names: {len(column_names)}. " \
                  f"len of unique_columns: {len(uniques)}.")
            sys.exit(0)

    # # plot the spectra
    for column in df.columns:
        # print(df[column])
        plt.plot(df.index[30:500], df[column][30:500], label=column)
    plt.legend(loc="upper right")
    plt.xlabel('wavelength (nm)')
    plt.ylabel('abs')
    # labelLines(plt.gca().get_lines(), zorder=2.5)
    plt.show()

    # save plot
    # plt.savefig(xml_folder + xml_name[:-4] + '.png')

    return df


if __name__ == "__main__":

    data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

    folders = [
        #            '2023-08-21-run01',
        #            '2023-08-22-run01',
        #            '2023-08-22-run02',
        #            '2023-08-28-run01',
        #            '2023-08-29-run01',
        #            '2023-08-29-run02',
        #            'reproducibility_test',
        #            'reference_for_simple_SN1',
        #            'background_test',
        #            'reference_for_E1'
        #            '2023-09-06-run01',
        #              '2023-09-07-run01',
        #               'Yanqiu'
        #               'versatility_test/Suzuki_JC'
        #               'versatility_test/Claisen_WaiShing'
        #               'versatility_test/Ullman_BP'
        #               'versatility_test/Glaser_WaiShing'
        #             '2023-11-28-run01',
        #             '2023-11-29-run01',
        #             '2023-11-29-run02',
        #             '2023-12-02-run01'
                    '2023-12-04-run02'

    ]

    for folder in folders:

        if 'run' in folder:
            xml_folder = data_folder + '/simple-reactions/' + folder + '/nanodrop_spectra/'
        else:
            xml_folder = data_folder + '/nanodrop-spectrometer-measurements/' + folder + '/'

        os.chdir(xml_folder)
        xml_files = glob.glob('*.xml')

        for xml_name in xml_files:
            print("Working on the following file: ")
            print(f'xml_folder: {xml_folder}')
            print(f'xml_name: {xml_name}')


            df = treat_one_file(xml_folder=xml_folder,
                            xml_name=xml_name)

            # if "2023" not in xml_name:
            #     df = treat_one_file(xml_folder=xml_folder,
            #                         xml_name=xml_name)
            # else:
            #     print('This file has been processed already!')
            #
            # if 'EtOH' in xml_name:
            #     df = treat_one_file(xml_folder = xml_folder,
            #                    xml_name = xml_name)
