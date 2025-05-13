from robowski.settings import *
import numpy as np
import pandas as pd
import os

from matplotlib import pyplot as plt

data_folder = os.environ['ROBOCHEM_DATA_PATH'].replace('\\', '/') + '/'

def exrobotocrudes_old():
    run_name = 'BPRF/2024-03-06-run02/'
    df_results_2 = pd.read_csv(data_folder + run_name + f'results/product_concentration.csv')
    df_results = pd.read_excel(data_folder + run_name + f'NMR/hantzschexnmroyes.xlsx', sheet_name=0)


    # xs = df_results.index.to_numpy()
    # ys = df_results['pc#HRP01'].to_numpy()
    # xs = df_results.index.to_numpy()
    # colstoplot = ['NMR HE C', 'OYES HE C']

    colstoplot = ['NMR HA C', 'OYES HA C']
    xs = df_results[colstoplot[0]].to_numpy()
    # ys = df_results[colstoplot[1]].to_numpy()
    ys = df_results_2['pc#bb017']
    # ys = df_results_2['pc#HRP01']
    plt.xlabel(colstoplot[0] + 'oncentration, mol/L')
    plt.ylabel(colstoplot[1] + 'oncentration, mol/L')

    maxval = max([np.max(xs), np.max(ys)])
    plt.plot([0, maxval], [0, maxval], 'k--')

    # set same color foe each three points
    labels = ['Conditions for HE max', 'Conditions for HA max', 'Conditions for HE max, repeat', 'Conditions for HA max, repeat']
    for i in range(4):
        plt.scatter(xs[i*3:i*3+3], ys[i*3:i*3+3], c=f'C{i}', marker='o', label=labels[i])

    # annotate points by index
    for i in range(len(xs)):
        plt.annotate(i, (xs[i], ys[i]))

    # plt.xlim(0, np.max(xs))
    # plt.ylim(0, np.max(ys))
    plt.legend()
    plt.show()

def inrobotocrudes_old(name):
    run_name = 'BPRF/2024-03-06-run01/'
    df_results_2 = pd.read_csv(data_folder + run_name + f'results/product_concentration.csv')
    df_results = pd.read_excel(data_folder + run_name + f'NMR/hantzschrobotnmroyes.xlsx', sheet_name=0)

    f1 = plt.figure(figsize=(7,7))

    if name == 'Hantzsch ester':
        ys = df_results['NMR HE C'].to_numpy()
        xs = df_results_2['pc#HRP01'].to_numpy()
        xs_err = df_results_2['pcerr#HRP01'].to_numpy()

    if name == 'Piperidone':
        ys = df_results['NMR HA C'].to_numpy()
        xs = df_results_2['pc#bb017'].to_numpy()
        xs_err = df_results_2['pcerr#bb017'].to_numpy()

    # plot with error bars in x
    plt.title('In roboto crudes, 2024-03-06-run01, $\lambda_{min}$=235 nm')
    dataset_dividing_indes = 9
    plt.errorbar(xs[:dataset_dividing_indes], ys[:dataset_dividing_indes], xerr=xs_err[:dataset_dividing_indes], fmt='o',
                 capsize=5, capthick=2, alpha=0.5, label='Repeated condition A\n(for Hantzsch ester maximum yield)', color='C0')
    plt.errorbar(xs[dataset_dividing_indes:], ys[dataset_dividing_indes:], xerr=xs_err[dataset_dividing_indes:], fmt='o',
                 capsize=5, capthick=2, alpha=0.5, label='Repeated condition B\n(for Piperidone maximum yield)', color='C2')
    maxval = max([np.max(xs), np.max(ys)])
    plt.plot([0, maxval], [0, maxval], 'k--', label='x=y line')
    for i in range(len(xs)):
        plt.annotate(i, (xs[i], ys[i]))
    plt.ylabel(f'{name} concentration by NMR, mol/L')
    plt.xlabel(f'{name} concentration by UV-VIS, mol/L')
    # plt.xlim(-0.05 * maxval, 1.1 * maxval)
    # plt.ylim(-0.05 * maxval, 1.1 * maxval)

    plt.xlim(-0.001, 0.012)
    plt.ylim(-0.001, 0.012)

    plt.legend()
    plt.show()
    # xs = df_results.index.to_numpy()
    # ys = df_results['pc#HRP01'].to_numpy()
    # xs = df_results.index.to_numpy()
    # colstoplot = ['NMR HE C', 'OYES HE C']


def inrobotocrudes(name, ax, runshortname):
    run_name = f'BPRF/{runshortname}/'
    df_results_2 = pd.read_csv(data_folder + run_name + f'results/product_concentration.csv')
    df_results = pd.read_excel(data_folder + run_name + f'NMR/hantzschrobotnmroyes.xlsx', sheet_name=0)
    # sort df_results by index column
    df_results = df_results.sort_values(by='index')

    # f1 = plt.figure(figsize=(7,7))

    if name == 'Hantzsch ester':
        ys = df_results['NMR HE C'].to_numpy()
        xs = df_results_2['pc#HRP01'].to_numpy()
        xs_err = df_results_2['pcerr#HRP01'].to_numpy()

    if name == 'Piperidone':
        ys = df_results['NMR HA C'].to_numpy()
        xs = df_results_2['pc#bb017'].to_numpy()
        xs_err = df_results_2['pcerr#bb017'].to_numpy()

    if name == 'Vinylamine':
        ys = df_results['dm88_4 C'].to_numpy()
        xs = df_results_2['pc#dm088_4'].to_numpy()
        xs_err = df_results_2['pcerr#dm088_4'].to_numpy()

    if name == 'dm070':
        ys = df_results['dm070 C'].to_numpy()
        # xs = df_results_2['pc#dm70'].to_numpy() + df_results_2['pc#dm37'].to_numpy()
        # xs_err = df_results_2['pcerr#dm70'].to_numpy() + df_results_2['pcerr#dm37'].to_numpy()

        xs = df_results_2['pc#dm70'].to_numpy()
        xs_err = df_results_2['pcerr#dm70'].to_numpy()
        # xs = ys
        # xs_err = np.ones_like(xs)*0.0001

    if name == 'dm053':
        ys = df_results['dm053 C'].to_numpy()
        # xs = ys
        # xs_err = np.ones_like(xs)*0.0001
        xs = df_results_2['pc#dm053'].to_numpy()
        xs_err = df_results_2['pcerr#dm053'].to_numpy()


    xs  *= 1000
    xs_err *= 1000
    ys *= 1000

    # plot with error bars in x
    # plt.title('In roboto crudes, 2024-03-12-run01, $\lambda_{min}$=235 nm')
    odd_indices = np.where(df_results['NMR HE C'].to_numpy() > 0.0005)
    even_indices = np.where(df_results['NMR HE C'].to_numpy() <= 0.0005)

    markers, caps, bars = ax.errorbar(xs[even_indices], ys[even_indices], xerr=xs_err[even_indices], fmt='o',
                 alpha=0.5, label='Repeated condition A\n(for Hantzsch ester maximum yield)', color='C0')
    [bar.set_alpha(0.3) for bar in bars]
    [cap.set_alpha(0.3) for cap in caps]
    markers, caps, bars = ax.errorbar(xs[odd_indices], ys[odd_indices], xerr=xs_err[odd_indices], fmt='o',
                 alpha=0.5, label='Repeated condition B\n(for Piperidone maximum yield)', color='C2')
    [bar.set_alpha(0.3) for bar in bars]
    [cap.set_alpha(0.3) for cap in caps]
    max_xs = max(xs)
    max_ys = max(ys)
    maxval = max([max_xs, max_ys])
    print(f'maxval: {maxval}')
    ax.plot([0, maxval], [0, maxval], 'k--', label='x=y line')
    # for i in range(len(xs)):
    #     plt.annotate(i, (xs[i], ys[i]))
    ax.set_ylabel(f'{name} concentration by NMR, mM')
    ax.set_xlabel(f'{name} concentration by UV-VIS, mM')
    ax.set_xlim(-0.05 * maxval, 1.1 * maxval)
    ax.set_ylim(-0.05 * maxval, 1.1 * maxval)

    # plt.xlim(-0.001, 0.012)
    # plt.ylim(-0.001, 0.012)

    # plt.legend()
    # plt.show()
    # xs = df_results.index.to_numpy()
    # ys = df_results['pc#HRP01'].to_numpy()
    # xs = df_results.index.to_numpy()
    # colstoplot = ['NMR HE C', 'OYES HE C']

if __name__ == '__main__':
    # inrobotocrudes(name='Hantzsch ester')
    # inrobotocrudes(name='Piperidone')
    # inrobotocrudes(name='dm88_4')
    # inrobotocrudes(name='dm070')
    # Vinylamine
    # f1 = plt.figure(10, 8)
    # make a gris of four subplots
    runshortname = '2024-03-12-run01'
    # runshortname = '2024-03-20-run01'
    # runshortname = '2024-03-06-run01'
    fig, axs = plt.subplots(2, 3, figsize=(12, 8))
    fig.suptitle(f'In roboto crudes, {runshortname}, ' + '$\lambda_{min}$=221 nm')
    for i, name in enumerate(['Hantzsch ester', 'Piperidone', 'Vinylamine', 'dm070', 'dm053']):
        inrobotocrudes(name, axs[i//3, i%3], runshortname=runshortname)

    plt.show()