'''Takes in the run_info.csv file, adds certain time delay to the finishing time and makes a table.
This comvenient for planning the next operation.'''

import pandas as pd
import datetime
import os

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'
experiment_name = 'multicomp-reactions/2023-07-04-run01/'

# df = pd.read_csv(data_folder + experiment_name + 'pipetter_io/run_info_0322.csv', delimiter=', ',
#                  usecols=['plate_code', 'experiment_name', 'finish_time_unix', 'note'])

# df = pd.read_csv(data_folder + experiment_name + 'pipetter_io/run_info.csv', delimiter=', ',
#                  usecols=['plate_code', 'experiment_name', 'finish_time_unix', 'note'])

df = pd.read_csv(data_folder + experiment_name + 'pipetter_io/run_info.csv', delimiter=',', header=0,
                 names=['plate_code', 'experiment_name', 'start_time_unix',
         'start_time_string', 'finish_time_unix', 'finish_time_string', 'note'])


delay_in_hours = 16
df['new_start_unixtime'] = pd.to_datetime(df['finish_time_unix'] + delay_in_hours * 60 * 60, unit='s')
df.new_start_unixtime = df.new_start_unixtime.dt.tz_localize('UTC').dt.tz_convert('Asia/Seoul')
# df = df.loc[df.new_start_unixtime.dt.date == datetime.datetime.now().date()]
df_out = df[['plate_code', 'new_start_unixtime']]

# df.new_start_unixtime = df.new_start_unixtime + pd.Timedelta('06:00:00')
print(df_out.to_string(index=False))

run_id = experiment_name.split('/')[1]
for row in df_out.itertuples():
    # format datetime to show only hour-minute-second
    print(run_id, f'**{row.plate_code}**', row.new_start_unixtime.strftime('%H:%M'))

print('5 minutes before:')
for row in df_out.itertuples():
    # format datetime to show only hour-minute-second
    print(run_id, f'**{row.plate_code}**', (row.new_start_unixtime - datetime.timedelta(hours=0, minutes=5)).strftime('%H:%M'))
