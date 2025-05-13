from robowski.settings import *
import os
import numpy as np
import pandas as pd

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

minc = 0.12
maxc = 0.3

unique_ptsas = np.loadtxt(data_folder + 'multicomp-reactions/2023-06-19-run01/outVandC/unique_cats.txt')
column_names = ['ald001', 'ptsa', 'am001', 'ic001']

ics = []
alds = []
amines = []
ptsas = []
for ic in np.linspace(minc, maxc, 4):
    for ald in np.linspace(minc, maxc, 4):
        for amine in np.linspace(minc, maxc, 4):
            for ptsa in unique_ptsas:
                ics.append(ic)
                alds.append(ald)
                amines.append(amine)
                ptsas.append(ptsa)

df = pd.DataFrame({'ic001': ics, 'ald001': alds, 'am001': amines, 'ptsa': ptsas}, dtype=object)

# remove the cube face where ic001 is 0.30 (maximum):
df.drop(df[df.ic001 == 0.3].index, inplace=True)

# add the plane where ic001 is 0.3 and there are 7 points across the ald and amine axes:
ics = []
alds = []
amines = []
ptsas = []
ic = 0.3
for ald in np.linspace(minc, maxc, 7):
    for amine in np.linspace(minc, maxc, 7):
        for ptsa in unique_ptsas:
            ics.append(ic)
            alds.append(ald)
            amines.append(amine)
            ptsas.append(ptsa)
# append this face to the dataframe:
df = df.append(pd.DataFrame({'ic001': ics, 'ald001': alds, 'am001': amines, 'ptsa': ptsas}), ignore_index=True, dtype=object)

# remove the cube face where amine is 0.12 (minimum):
df.drop(df[df.am001 == 0.12].index, inplace=True)

# add the plane where ic001 is 0.3 and there are 7 points across the ald and amine axes:
ics = []
alds = []
amines = []
ptsas = []
amine = 0.12
for ald in np.linspace(minc, maxc, 7):
    for ic in np.linspace(minc, maxc, 7):
        for ptsa in unique_ptsas:
            ics.append(ic)
            alds.append(ald)
            amines.append(amine)
            ptsas.append(ptsa)
# append this face to the dataframe:
df = df.append(pd.DataFrame({'ic001': ics, 'ald001': alds, 'am001': amines, 'ptsa': ptsas}), ignore_index=True, dtype=object)

# remove the cube face where ald001 is 0.12 (minimum):
df.drop(df[df.ald001 == 0.12].index, inplace=True)

# add the plane where ic001 is 0.3 and there are 7 points across the ald and amine axes:
ics = []
alds = []
amines = []
ptsas = []
ald = 0.12
for amine in np.linspace(minc, maxc, 7):
    for ic in np.linspace(minc, maxc, 7):
        for ptsa in unique_ptsas:
            ics.append(ic)
            alds.append(ald)
            amines.append(amine)
            ptsas.append(ptsa)
# append this face to the dataframe:
df = df.append(pd.DataFrame({'ic001': ics, 'ald001': alds, 'am001': amines, 'ptsa': ptsas}), ignore_index=True, dtype=object)

# randomize the order of rows
df = df.sample(frac=1).reset_index(drop=True)

df.to_csv(data_folder + 'multicomp-reactions/2023-06-19-run01/' + 'outVandC/outC.csv')
