import os
import numpy as np
import pandas as pd

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'


def load_outVandC(run_folder, suffix=''):
    true_c_df = pd.read_csv(run_folder + f'outVandC/outC{suffix}.csv')
    true_c_df.drop('Unnamed: 0', inplace=True, axis=1)
    true_v_df = pd.read_csv(run_folder + f'outVandC/outV{suffix}.csv')
    true_v_df.drop('Unnamed: 0', inplace=True, axis=1)
    return true_c_df, true_v_df


# These files are made by Rafal using his own code. So they are the "true" concentrations and volumes.
true_c_df, true_v_df = load_outVandC(run_folder=data_folder + 'multicomp-reactions/2023-03-20-run01/',
                                     suffix='RF038202303201421')

ptsa_c_thresh = 0.08

xs = true_c_df['ptsa'].to_numpy()
ys = true_v_df['ptsa_dil_x_5'].to_numpy()
mask = np.logical_and(xs < ptsa_c_thresh, xs > 0.016)
slope1_ptsa, _ = np.polyfit(xs[mask], ys[mask], 1)


def ptsa_concentration_to_ptsa_dil_x_5_volume(ptsa):
    if ptsa >= ptsa_c_thresh:
        return 0
    else:
        return slope1_ptsa * ptsa


assert np.isclose(np.array([ptsa_concentration_to_ptsa_dil_x_5_volume(c) for c in xs[mask]]),
                  ys[mask]).all()

xs = true_c_df['ptsa'].to_numpy()
ys = true_v_df['ptsa'].to_numpy()
slope2_ptsa, _ = np.polyfit(xs[xs > ptsa_c_thresh], ys[xs > ptsa_c_thresh], 1)


def ptsa_concentration_to_ptsa_volume(ptsa):
    if ptsa < ptsa_c_thresh:
        return 0
    else:
        return slope2_ptsa * ptsa


assert np.isclose(np.array([ptsa_concentration_to_ptsa_volume(c) for c in true_c_df['ptsa'].to_numpy()]),
                  true_v_df['ptsa'].to_numpy()).all()

slopes_by_substrate = {substrate: np.polyfit(true_c_df[substrate].to_numpy(), true_v_df[substrate].to_numpy(), 1)[0]
                       for substrate in ['ald001', 'ic001', 'am001']}


def substrate_concentration_to_volume(substrate, substrate_concentration):
    return slopes_by_substrate[substrate] * substrate_concentration


for substrate in ['ald001', 'ic001', 'am001']:
    assert np.isclose(
        np.array([substrate_concentration_to_volume(substrate, c) for c in true_c_df[substrate].to_numpy()]),
        true_v_df[substrate].to_numpy()).all()

if __name__ == '__main__':
    experiment_name = 'multicomp-reactions/2023-06-19-run01/'

    concentrations_df = pd.read_csv(data_folder + experiment_name + 'outVandC/' + 'outC.csv')
    concentrations_df.drop('Unnamed: 0', inplace=True, axis=1)
    volumes_df = concentrations_df.copy()

    # convert concentrations to volumes
    for substrate in ['ald001', 'ic001', 'am001']:
        volumes_df[substrate] = concentrations_df[substrate].apply(
            lambda x: substrate_concentration_to_volume(substrate, x))
    volumes_df['ptsa'] = concentrations_df['ptsa'].apply(lambda x: ptsa_concentration_to_ptsa_volume(x))
    volumes_df['ptsa_dil_x_5'] = concentrations_df['ptsa'].apply(lambda x: ptsa_concentration_to_ptsa_dil_x_5_volume(x))
    volumes_df['DMF'] = 200 - volumes_df['ald001'] - volumes_df['ic001'] - volumes_df['am001'] - volumes_df['ptsa'] - \
                        volumes_df['ptsa_dil_x_5']

    # rearrange columns to reflect the addition order
    volumes_df = volumes_df[['DMF', 'ald001', 'ptsa', 'ptsa_dil_x_5', 'am001', 'ic001']]

    volumes_df.to_csv(data_folder + experiment_name + 'outVandC/outV.csv')
