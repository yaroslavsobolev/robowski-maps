from robowski.settings import *
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.ticker as mtick


df_results = pd.read_csv('interpolated_product_concentration.csv')
# df_results['yield'] = df_results['yield'].apply(lambda x: 0 if x < 0 else x)
substances = ['ic001','am001','ald001','ptsa']
substance_titles = ['Isocyanide', 'Amine', 'Aldehyde', 'p-TSA']

# df_results = df_results[(df_results['ic001'] == 0.3) & (df_results['ald001'] == 0.15)]
# # unique amine

unique_amine = list(df_results['am001'].unique())
sorted_amine = np.sort(unique_amine)
print(f'sorted_amine: {sorted_amine}')

# unique aldehyde
unique_aldehyde = df_results['ald001'].unique()
print(f'unique_aldehyde: {unique_aldehyde}')

def division_figure():
    yields = []
    prods = []
    lims = []
    indices_to_plot = []
    for i, amine in enumerate(sorted_amine):
        # if amine = 0.24
        df_local = df_results[np.isclose(df_results['ic001'], 0.3) & np.isclose(df_results['am001'], amine) & (np.isclose(df_results['ald001'], 0.3 + 0.12 - amine))]

        # find the value of ptsa closest to 'amine" value
        ptsa_all = df_local['ptsa'].unique()
        ptsa_closest_to_am = ptsa_all[np.argmin(np.abs(ptsa_all - amine))]
        df_local = df_local[df_local['ptsa'] == ptsa_closest_to_am]

        aldehyde_here = 0.3 + 0.12 - amine
        limiting_reagent = min(amine, aldehyde_here)
        max_yield = df_local['yield'].max()
        product_concentration = max_yield * limiting_reagent
        lims.append(limiting_reagent)
        prods.append(product_concentration)

        yields.append(max_yield)
        # find the ptsa value that corresponds to the highest yield
        ptsa_value = df_local[df_local['yield'] == max_yield]['ptsa'].values[0]
        print(f'amine: {amine}, ptsa: {ptsa_closest_to_am}, yield: {max_yield}')
        # find ptsa at max yield
        print(f'BY MAX: amine: {amine}, ptsa: {ptsa_value}, yield: {max_yield}')

    # model results

    df_model = pd.read_csv('ugi_v3_REV_outputs/ugi_model_yields_on_diagonal_85c_2025-02-15.csv')
    # df_model = df_model[df_model['yield']>-1]
    xs_model = df_model['am001']
    # aldehyde_model = 0.3 + 0.12 - xs_model
    ys_model = df_model['yield']
    limiting_reagent_2 = np.minimum(xs_model, 0.3 + 0.12 - xs_model)
    ys_model_conc = ys_model * limiting_reagent_2

    factor = 0.85
    fig, axarr = plt.subplots(3, figsize=(5*factor, 10*factor), sharex=True)
    axarr[2].plot(sorted_amine, 100*np.array(yields), 'o', color='k')
    axarr[2].plot(xs_model, 100*ys_model, '-', color='k')
    axarr[2].set_ylabel('Yield')
    # fmt = '%.0f%%' # Format you want the ticks, e.g. '40%'
    # xticks = mtick.FormatStrFormatter(fmt)
    # axarr[2].yaxis.set_major_formatter(xticks)
    axarr[2].yaxis.set_major_formatter(mtick.PercentFormatter())

    axarr[0].plot(sorted_amine, 1000*np.array(prods), 'o', color='m')
    axarr[0].plot(xs_model, 1000*ys_model_conc, '-', color='m')
    axarr[0].set_ylim(0, 1000*0.019)
    axarr[0].set_ylabel('Product concentration (mM)\nfor optimal p-TSA concentration\n([p-TSA]$\\approx$[amine])')

    axarr[1].plot(sorted_amine, lims, 'o-', color='b')
    axarr[1].set_ylim(0, 0.221)
    axarr[1].set_ylabel('Limiting reagent concentration, M')
    # set ylim of ax2 to (0, 0.19)

    axarr[2].set_xlabel('Amine concentration, M')
    axarr[2].set_xticks(sorted_amine)
    axarr[2].set_ylim(0, 20)

    plt.tight_layout()

    ax1 = axarr[2]
    ax2 = ax1.twiny()

    # Add some extra space for the second axis at the bottom
    fig.subplots_adjust(bottom=0.15)

    new_tick_locations = np.array(sorted_amine)

    def tick_function(X):
        V = 0.3 + 0.12 - X
        return ["%.2f" % z for z in V]

    # Move twinned axis ticks and label from top to bottom
    ax2.xaxis.set_ticks_position("bottom")
    ax2.xaxis.set_label_position("bottom")

    # Offset the twin axis below the host
    ax2.spines["bottom"].set_position(("axes", -0.31))

    # Turn on the frame for the twin axis, but then hide all
    # but the bottom spine
    ax2.set_frame_on(True)
    ax2.patch.set_visible(False)

    # as @ali14 pointed out, for python3, use this
    # for sp in ax2.spines.values():
    # and for python2, use this
    # for sp in ax2.spines.itervalues():
    #     sp.set_visible(False)
    ax2.spines["bottom"].set_visible(True)
    ax2.set_xlim(ax1.get_xlim())

    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(tick_function(new_tick_locations))
    ax2.set_xlabel(r"Aldehyde concentration, M")
    fig.savefig('ugi_diagonal_max_v2.png', dpi=300)
    plt.show()


def p3p2_figure():
    df_model = pd.read_csv('ugi_v3_REV_outputs/ugi_model_yields_on_diagonal_p3p2_2025-02-16.csv')
    # df_model = df_model[df_model['yield']>-1]
    xs_model = df_model['am001']
    # aldehyde_model = 0.3 + 0.12 - xs_model
    ys_model = df_model['p3/p2']
    limiting_reagent_2 = np.minimum(xs_model, 0.3 + 0.12 - xs_model)
    ys_model_conc = ys_model * limiting_reagent_2

    factor = 0.85
    fig, ax0 = plt.subplots(1, figsize=(5 * factor, 4 * factor), sharex=True)
    ax0.plot(xs_model, ys_model, '-', color='k')

    ax0.set_xlabel('Amine concentration, M')
    ax0.set_xticks(sorted_amine)
    # ax0.set_ylim(0, 20)
    ax0.set_ylabel('[P3] / [P2]')

    plt.tight_layout()

    ax1 = ax0
    ax2 = ax1.twiny()

    # Add some extra space for the second axis at the bottom
    fig.subplots_adjust(bottom=0.35)

    new_tick_locations = np.array(sorted_amine)

    def tick_function(X):
        V = 0.3 + 0.12 - X
        return ["%.2f" % z for z in V]

    # Move twinned axis ticks and label from top to bottom
    ax2.xaxis.set_ticks_position("bottom")
    ax2.xaxis.set_label_position("bottom")

    # Offset the twin axis below the host
    ax2.spines["bottom"].set_position(("axes", -0.31))

    # Turn on the frame for the twin axis, but then hide all
    # but the bottom spine
    ax2.set_frame_on(True)
    ax2.patch.set_visible(False)

    # as @ali14 pointed out, for python3, use this
    # for sp in ax2.spines.values():
    # and for python2, use this
    # for sp in ax2.spines.itervalues():
    #     sp.set_visible(False)
    ax2.spines["bottom"].set_visible(True)
    ax2.set_xlim(ax1.get_xlim())

    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(tick_function(new_tick_locations))
    ax2.set_xlabel(r"Aldehyde concentration, M")
    fig.savefig('ugi_diagonal_p3p2.png', dpi=300)
    plt.show()

if __name__ == '__main__':
    # division_figure()
    p3p2_figure()
    pass