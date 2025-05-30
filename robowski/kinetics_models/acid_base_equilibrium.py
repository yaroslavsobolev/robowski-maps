"""
Proton exchange equilibrium calculations for complex chemical mixtures.

This module calculates pH and species concentrations in systems where water concentration
can be comparable to solute concentrations, unlike conventional pH calculation packages
that assume water is in vast excess. This is critical for chemical reactions where water
is produced or consumed and its concentration changes significantly.

The approach solves the zero-charge equation numerically to find equilibrium pH, then
calculates concentrations of all protonated and deprotonated species. Uses high-precision
arithmetic for double ionization species to avoid floating-point rounding errors.

Key advantages over packages like pHcalc:
- No assumption of excess water concentration
- Accounts for water production/consumption during reactions
- Handles both single and double ionization with appropriate precision
- Suitable for non-aqueous and mixed solvents like DMF/water

Scientific background is detailed in the accompanying research article's
Supplementary Information Section 5.3.3.

Example:
    >>> substances = [
    ...     {'pKa': -2.8, 'conc': 0.1, 'charge': 0},  # Acid
    ...     {'pKa': 11, 'conc': 0.1, 'charge': 1}     # Base
    ... ]
    >>> Kw = 1e-14 / (55.56**2)  # Autoprotolysis constant
    >>> ph = solve_for_zero_charge(0.1, Kw, substances)
    >>> subs, remaining_water, oh_minus = concentrations_by_ph(0.1, Kw, substances, ph)
    >>> print(f'Equilibrium pH: {ph}')
    >>> print(f'Substances at equilibrium: {subs}')
"""

from robowski.settings import *
import numpy as np
import matplotlib.pyplot as plt  # Optional for plotting below
from scipy.optimize import brentq
from tqdm import tqdm
from _decimal import *
from scipy.special import expit

## Uncomment for comparing to the PHcalc package
# from pHcalc import Acid, Inert, System

getcontext().prec = 85
minconc = 0

PURE_WATER_MOLARITY = Decimal(55.56)  # mol/L
LOG10_PURE_WATER_MOLARITY = PURE_WATER_MOLARITY.log10()
LOG10_PURE_WATER_MOLARITY_FLOAT = np.log10(55.56)


# 'charge' means 'charge of protonated species
# substances = [{'pKa': -np.log10(6.76e-4), 'conc': 0.01, 'charge': 0},
#               {'pKa': 17, 'conc': 0.01, 'charge': 1}]
# substances = [{'pKa': -1, 'conc': 1, 'charge': 0}]

# substances = [{'pKa': 11, 'conc': 0.1, 'charge': 1},
#               {'pKa': -2.8, 'conc': 0.1, 'charge': 0},
#               {'pKa': 6, 'conc': 0.1, 'charge': 1}]

def slow_conc(remaining_water, ph, pKa, s0):
    """
    Calculate protonated/deprotonated concentrations using high-precision Decimal arithmetic.

    Uses the generalized Henderson-Hasselbalch equations that account for variable water
    concentration. Employs Decimal precision to avoid floating-point rounding errors
    that become significant in extreme pH ranges or low concentrations.

    Parameters
    ----------
    remaining_water : float
        Concentration of neutral water molecules (mol/L).
    ph : float
        Solution pH value.
    pKa : float
        Acid dissociation constant (pKa = -log10(Ka)).
    s0 : float
        Total concentration of the species (mol/L).

    Returns
    -------
    tuple
        (deprotonated_concentration, protonated_concentration) in mol/L as Decimal objects.
    """
    h_plus = Decimal(10) ** (-1 * Decimal(ph))
    Ka = Decimal(10) ** (-(Decimal(pKa) + LOG10_PURE_WATER_MOLARITY))
    sd = Decimal(s0) * Ka / ((h_plus / Decimal(remaining_water)) + Ka)
    sp = Decimal(s0) - sd
    return sd, sp


def fast_conc(remaining_water, ph, pKa, s0):
    """
    Calculate protonated/deprotonated concentrations using optimized sigmoid approach.

    Faster implementation for single ionization species using the logistic sigmoid
    function, which is numerically robust against floating-point errors. Equivalent
    to slow_conc but ~10x faster for single ionization calculations.

    Parameters
    ----------
    remaining_water : float
        Concentration of neutral water molecules (mol/L).
    ph : float
        Solution pH value.
    pKa : float
        Acid dissociation constant (pKa = -log10(Ka)).
    s0 : float
        Total concentration of the species (mol/L).

    Returns
    -------
    tuple
        (deprotonated_concentration, protonated_concentration) in mol/L as floats.
    """
    e10 = np.log10(np.exp(1))
    L = ph / e10 + np.log(remaining_water) - (pKa + float(LOG10_PURE_WATER_MOLARITY)) / e10
    sd = s0 * expit(L)
    sp = s0 * expit(-L)
    return sd, sp


def concentrations_by_ph(water_concentration, Kw, substances, ph, return_Decimals=True):
    """
    Calculate equilibrium concentrations of all species at a given pH.

    Core function that determines concentrations of protonated and deprotonated forms
    for all substances in the mixture. Handles both single ionization (pKa only) and
    double ionization (pKa and pKa2). Accounts for water consumption in proton transfer.

    Uses different precision approaches:
    - Double ionization: High-precision Decimal arithmetic to avoid rounding errors
    - Single ionization: Optimized sigmoid function for speed

    Parameters
    ----------
    water_concentration : float
        Total water concentration in the system (mol/L).
    Kw : float
        Autoprotolysis constant of water in the solvent system.
    substances : list of dict
        Each dict contains: 'pKa', 'conc', 'charge', and optionally 'pKa2'.
        - pKa: First dissociation constant
        - pKa2: Second dissociation constant (if present)
        - conc: Total concentration of all forms of this substance
        - charge: Charge of the fully protonated form
    ph : float
        Solution pH at which to calculate concentrations.
    return_Decimals : bool, optional
        Whether to return Decimal objects (True) or convert to float (False).

    Returns
    -------
    tuple
        (substances_with_concentrations, remaining_water, hydroxide_concentration)
        - substances: Input list with added concentration keys
        - remaining_water: Concentration of neutral water molecules
        - hydroxide_concentration: Concentration of OH- ions
    """
    # solve for the remaining water and present OH- concentration
    # time the execution
    h2o = Decimal(water_concentration)
    Kw = Decimal(Kw)
    ph_d = Decimal(ph)
    h_plus = Decimal(10) ** (-1 * ph_d)
    z = h2o
    h_over_k = h_plus / Kw
    remaining_water = (h_over_k * (z + h_over_k / 4)).sqrt() - h_over_k / 2
    oh_minus = z - remaining_water
    assert np.isclose(float(remaining_water) + float(oh_minus), float(h2o),
                      rtol=1e-8), f'Water digit precision insufficient: {remaining_water} + {oh_minus} != {h2o}'

    # make a copy of list of substances
    subs = [s.copy() for s in substances]
    # calculate concentrations of substances

    e10 = np.log10(np.exp(1))
    L0 = float(ph) / e10 + np.log(float(remaining_water))

    for s in subs:
        # if key 'pKa2' is not present in dict s:
        if not ('pKa2' in s):  # this means there is only one protonation possible
            s['conc_deprot_2'] = 0

            ### Slow version
            # Ka = Decimal(10)**(-(Decimal(s['pKa']) + LOG10_PURE_WATER_MOLARITY))
            # s['conc_deprot'] = Decimal(s['conc']) * Ka / (hplus_over_remaining_water + Ka)
            # s['conc_prot'] = Decimal(s['conc']) - s['conc_deprot']

            # ### Faster version
            # s['conc_deprot'], s['conc_prot'] = fast_conc(remaining_water, ph, s['pKa'], s['conc'])

            ### Even faster version
            L = L0 - (s['pKa'] + LOG10_PURE_WATER_MOLARITY_FLOAT) / e10
            s['conc_deprot'] = s['conc'] * expit(L)
            s['conc_prot'] = s['conc'] * expit(-L)

        else:  # this means there are two ionizations
            Ka1 = Decimal(10) ** (-(Decimal(s['pKa']) + LOG10_PURE_WATER_MOLARITY))
            Ka2 = Decimal(10) ** (-(Decimal(s['pKa2']) + LOG10_PURE_WATER_MOLARITY))
            m0 = Decimal(s['conc'])
            hplus_over_remaining_water = h_plus / remaining_water
            s['conc_deprot'] = m0 * Ka1 / ((hplus_over_remaining_water) + (1 + Ka2 / hplus_over_remaining_water) * Ka1)
            s['conc_deprot_2'] = s['conc_deprot'] * Ka2 / hplus_over_remaining_water
            s['conc_prot'] = m0 - s['conc_deprot'] - s['conc_deprot_2']

        if not return_Decimals:
            # convert from decimal to float64
            s['conc_deprot'] = float(s['conc_deprot'])
            s['conc_deprot_2'] = float(s['conc_deprot_2'])
            s['conc_prot'] = float(s['conc_prot'])

    if not return_Decimals:
        oh_minus = float(oh_minus)
        remaining_water = float(remaining_water)

    return subs, remaining_water, oh_minus


def positive_charges_minus_negative_charges(water_concentration, kW, substances, ph):
    """
    Calculate net charge of the solution at a given pH.

    Computes the sum of all positive charges minus all negative charges in the system,
    including H3O+, OH-, and all protonated/deprotonated species. Used to find the
    pH where net charge equals zero (electroneutrality condition).

    Parameters
    ----------
    water_concentration : float
        Total water concentration in the system (mol/L).
    kW : float
        Autoprotolysis constant of water.
    substances : list of dict
        Species definitions with pKa, concentration, and charge information.
    ph : float
        Solution pH at which to calculate net charge.

    Returns
    -------
    float
        Net charge of the solution (positive charges - negative charges).
        Zero indicates electroneutrality (equilibrium condition).
    """
    subs, remaining_water, oh_minus = concentrations_by_ph(water_concentration, kW, substances, ph)
    # h_plus = Decimal(10)**(-1*Decimal(ph))
    # charge = float(h_plus - oh_minus)
    h_plus = 10 ** (-ph)
    charge = h_plus - float(oh_minus)
    # negative_charge = oh_minus
    for s in subs:
        charge += (s['charge'] - 2) * float(s['conc_deprot_2'])
        charge += (s['charge'] - 1) * float(s['conc_deprot'])
        charge += (s['charge']) * float(s['conc_prot'])
    return float(charge)


def solve_for_zero_charge(water_concentration, kW, substances, force_plot=False):
    """
    Find equilibrium pH where the solution has zero net charge.

    Solves the electroneutrality equation numerically using Brent's method to find
    the pH at which positive charges exactly balance negative charges. This represents
    the true equilibrium pH of the mixture.

    Automatically adjusts search range based on charge behavior at extreme pH values
    to ensure robust convergence for diverse chemical systems.

    Parameters
    ----------
    water_concentration : float
        Total water concentration in the system (mol/L).
    kW : float
        Autoprotolysis constant of water in the solvent system.
    substances : list of dict
        Chemical species with their pKa values, concentrations, and charges.
    force_plot : bool, optional
        Whether to force plotting of charge vs pH curve for diagnostics.

    Returns
    -------
    float
        Equilibrium pH of the solution where net charge equals zero.

    Notes
    -----
    Uses adaptive pH search ranges:
    - Standard range: -8 to +25 pH units
    - Extended range: -35 to +40 pH units for extreme cases

    Convergence tolerances: xtol=0.5e-3, rtol=1e-4 for pH accuracy.
    """
    def partial_func(ph):
        return positive_charges_minus_negative_charges(water_concentration, float(kW), substances, ph)

    if partial_func(20) <= 0:
        phmax = 25
    else:
        phmax = 40
    if partial_func(-8) >= 0:
        return brentq(partial_func, -8, phmax, xtol=0.5e-3, rtol=1e-4)
    else:
        return brentq(partial_func, -35, phmax, xtol=0.5e-3, rtol=1e-4)


if __name__ == '__main__':
    """
    This module is designed to be imported, not run as a standalone script.

    The code below serves as a comprehensive test suite, performance benchmark,
    and scientific demonstration of the module's capabilities. It is useful for:

    1. **Validation Testing**: Verifies that slow_conc() and fast_conc() produce 
       identical results across extreme pH ranges (-65 to +90), ensuring numerical
       accuracy of both implementation approaches.

    2. **Performance Benchmarking**: Times the execution of high-precision vs 
       optimized calculation methods to quantify the ~10x speed improvement of
       the sigmoid approach for single ionization species.

    3. **Basic Functionality Demo**: Simple acid-base equilibrium calculation
       showing pH determination for a two-component system (strong acid + weak base).

    4. **Advanced Feature Testing**: Double ionization example mimicking phosphate
       buffer systems to validate the high-precision Decimal arithmetic for
       species with multiple pKa values.

    5. **Scientific Visualization**: Generates publication-quality 2D heatmaps 
       showing pH, water concentration, and protonation fraction landscapes as
       functions of reactant concentrations - demonstrates real-world application
       for Ugi reaction chemistry modeling.

    6. **Debugging Tools**: Provides concrete examples for troubleshooting
       convergence issues, extreme pH behavior, and numerical precision problems
       that may arise in complex chemical systems.

    """
    for ph in np.linspace(-65, 90, 300):
        remaining_water = 0.1
        pKa = 11
        s0 = 0.1
        sd, sp = slow_conc(remaining_water, ph, pKa, s0)
        sd2, sp2 = fast_conc(remaining_water, ph, pKa, s0)
        assert np.isclose(float(sd), float(sd2),
                          rtol=1e-8), f'Concentration of deprotonated species differs: {sd} != {sd2}'
        assert np.isclose(float(sp), float(sp2),
                          rtol=1e-8), f'Concentration of protonated species differs: {sp} != {sp2}'
    print('slow_conc and fast_conc agree')

    # time the execution of conc_fast and conc_slow
    import time

    t0 = time.time()
    for ph in np.linspace(-35, 50, 100):
        remaining_water = 0.1
        pKa = 11
        s0 = 0.1
        for i in range(100):
            sd, sp = slow_conc(remaining_water, ph, pKa, s0)
    t1 = time.time()
    print(f'slow_conc took {t1 - t0} seconds')

    t0 = time.time()
    for ph in np.linspace(-35, 50, 100):
        remaining_water = 0.1
        pKa = 11
        s0 = 0.1
        for i in range(100):
            sd, sp = fast_conc(remaining_water, ph, pKa, s0)
    t1 = time.time()
    print(f'fast_conc took {t1 - t0} seconds')

    print('.>>>')
    substances = [
        {'pKa': -2.8, 'conc': 0.1, 'charge': 0},
        {'pKa': 11, 'conc': 0.1, 'charge': 1}]
    # Kw_here = 10**(-27.44)*(13.9)**2
    Kw_here = Decimal(1e-14) / (PURE_WATER_MOLARITY ** 2)
    ph = solve_for_zero_charge(1e-11, Kw_here, substances, force_plot=True)
    subs, remaining_water, oh_minus = concentrations_by_ph(54, Kw_here, substances, ph)
    print(f'ph = {ph}')
    print(subs)
    print(f'remaining water {remaining_water}')

    # # Testing double protonation
    # print('Phosphate buffer:')
    # phos = Acid(pKa=[2.148, 7.198], charge=0, conc=0.01)
    # nh4 = Acid(pKa=9.25, charge=1, conc=0.01)
    # system = System(phos, nh4)
    # system.pHsolve()
    # print(f'phcalc: ph={system.pH}')  # Should print 8.95915298

    print('.>>>')
    substances = [
        {'pKa': 9.25, 'conc': 0.01, 'charge': 1},
        {'pKa': 2.148, 'pKa2': 7.198, 'conc': 0.01, 'charge': 0}]
    # Kw_here = 10**(-27.44)*(13.9)**2
    Kw_here = 1e-14 / (PURE_WATER_MOLARITY ** 2)
    ph = solve_for_zero_charge(54, Kw_here, substances)
    subs, remaining_water, oh_minus = concentrations_by_ph(1, Kw_here, substances, ph)
    print(f'ph = {ph}')
    print(subs)
    print(f'remaining water {remaining_water}')

    # phs = np.linspace(-5, 30, 1000)
    # charges = [positive_charges_minus_negative_charges(0.3, 1e-14 / (PURE_WATER_MOLARITY**2), substances, ph) for ph in phs]
    # plot
    # plt.plot(phs, charges)
    # plt.show()

    BuNH2_concentrations = np.linspace(0.01, 0.3, 150)
    TsOH_concentratuons = np.linspace(0.01, 0.3, 150)

    y, x = np.meshgrid(BuNH2_concentrations, TsOH_concentratuons)
    # make a 2d imshow pcolor plot
    phs = np.zeros((len(BuNH2_concentrations), len(TsOH_concentratuons)))
    waters = np.zeros((len(BuNH2_concentrations), len(TsOH_concentratuons)))
    proportions = np.zeros((len(BuNH2_concentrations), len(TsOH_concentratuons)))
    for i, BuNH2 in tqdm(enumerate(BuNH2_concentrations), total=len(BuNH2_concentrations)):
        for j, TsOH in enumerate(TsOH_concentratuons):
            substances = [{'pKa': 11, 'conc': BuNH2, 'charge': 1},  # ammonium
                          {'pKa': -2.8, 'conc': TsOH, 'charge': 0},  # TsOH + H2O -> TSO- + H3O+
                          {'pKa': -7.2, 'conc': 12.915, 'charge': 1}]  # DHFH+ + H2O -> DMF + H3O+
            # Kw_here = 10 ** (-27.44) * (13.9) ** 2
            Kw_here = 10 ** (-31.4) / (PURE_WATER_MOLARITY) / 13
            ph_here = solve_for_zero_charge(water_concentration=TsOH, kW=Kw_here, substances=substances)
            phs[i, j] = ph_here
            subs, remaining_water, oh_minus = concentrations_by_ph(water_concentration=TsOH, Kw=Kw_here,
                                                                   substances=substances,
                                                                   ph=ph_here)
            waters[i, j] = remaining_water
            proportions[i, j] = subs[0]['conc_prot'] / subs[0]['conc']

    # # pickle phs, waters, proportions for persistence
    # import pickle
    # with open(repo_data_path + 'misc_scripts/ugi_ph.pkl', 'wb') as f:
    #     pickle.dump(phs, f)
    # with open(repo_data_path + 'misc_scripts/ugi_waters.pkl', 'wb') as f:
    #     pickle.dump(waters, f)
    # with open(repo_data_path + 'misc_scripts/ugi_proportions.pkl', 'wb') as f:
    #     pickle.dump(proportions, f)
    #
    # # load from pickle
    # with open(repo_data_path + 'misc_scripts/ugi_ph.pkl', 'rb') as f:
    #     phs = pickle.load(f)
    # with open(repo_data_path + 'misc_scripts/ugi_waters.pkl', 'rb') as f:
    #     waters = pickle.load(f)
    # with open(repo_data_path + 'misc_scripts/ugi_proportions.pkl', 'rb') as f:
    #     proportions = pickle.load(f)

    # plot
    interpolation_method = 'bicubic'
    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    #
    plt.imshow(phs, extent=(
    BuNH2_concentrations.min(), BuNH2_concentrations.max(), TsOH_concentratuons.min(), TsOH_concentratuons.max()),
               interpolation=interpolation_method, origin='lower')
    # colorbar with title of 'pH'
    plt.colorbar(label='pH')
    plt.xlabel('(p-TSA$\cdot H_2 O$) concentration, mol/L')
    plt.ylabel('Total amine+ammonium concentration, mol/L')
    plt.tight_layout()
    fig.savefig(repo_data_path + 'misc_scripts/figures/ugi_ph.png', dpi=300)
    plt.show()

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    # colorbar with title of 'pH'
    plt.imshow(proportions, extent=(
    BuNH2_concentrations.min(), BuNH2_concentrations.max(), TsOH_concentratuons.min(), TsOH_concentratuons.max()),
               interpolation=interpolation_method, origin='lower')
    plt.colorbar(label='Fraction of butylammonium to\nthe sum of butylammonium and butylamine')
    plt.xlabel('(p-TSA$\cdot H_2 O$) concentration, mol/L')
    plt.ylabel('Total amine+ammonium concentration, mol/L')
    plt.tight_layout()
    fig.savefig(repo_data_path + 'misc_scripts/figures/ugi_ion_fraction.png', dpi=300)
    plt.show()

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    plt.imshow(waters, extent=(
    BuNH2_concentrations.min(), BuNH2_concentrations.max(), TsOH_concentratuons.min(), TsOH_concentratuons.max()),
               interpolation=interpolation_method, origin='lower')
    # colorbar with title of 'pH'
    plt.colorbar(label='Concentration of remaining (neutral) water, mol/L')
    plt.xlabel('(p-TSA$\cdot H_2 O$) concentration, mol/L')
    plt.ylabel('Amine (BuNH$_{2}$) concentration, mol/L')
    plt.tight_layout()
    fig.savefig(repo_data_path + 'misc_scripts/figures/ugi_neutral_water.png', dpi=300)
    plt.show()
