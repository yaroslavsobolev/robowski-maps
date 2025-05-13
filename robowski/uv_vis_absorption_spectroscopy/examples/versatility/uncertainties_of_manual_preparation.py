from robowski.settings import *
import numpy as np
import uncertainties as unc

# Friedel-Crafts. Author: Galimzhan
powder_weight = unc.ufloat(186, 0.6)
solvent_volume = unc.ufloat(20, 0.02)

aliquot = unc.ufloat(1, 0.01)
net_volume = unc.ufloat(100, 0.1)

result = (powder_weight / solvent_volume) * (aliquot / net_volume)

relative_std_of_result = result.std_dev / result.nominal_value
print(f'Relative standard deviation of result: {relative_std_of_result:.2%}, and -0.5% must be added due to neglect of volume')


# Imine. Author: Rafal
print('Imine_RF')
# calibrations
powder_weight = unc.ufloat(200, 0.6)
solvent_volume = unc.ufloat(20, 0.02)
aliquit_by_hamilton_syringe = unc.ufloat(1, 0.01)
net_volume = unc.ufloat(100, 0.1)
result = (powder_weight / solvent_volume) * (aliquit_by_hamilton_syringe / net_volume)
relative_std_of_result = result.std_dev / result.nominal_value
print(f'Relative standard deviation of result: {relative_std_of_result:.2%}')

# manual preparation of reaction mixture
powder_weight = unc.ufloat(400, 0.6)
solvent_volume = unc.ufloat(25, 0.025) + unc.ufloat(250, 2.5) + unc.ufloat(11, 25/100)
result = powder_weight / solvent_volume
relative_std_of_result = result.std_dev / result.nominal_value
print(f'Relative standard deviation of result: {relative_std_of_result:.2%}')

# Suzuki. Author: Juan Carlos
print('Suzuki_JC')
# stock
powder_weight = unc.ufloat(90, 0.6)
solvent_volume = unc.ufloat(5, 5 * np.sqrt(0.003**2 + 0.008**2))
stock_concentration = powder_weight / solvent_volume

# calibrations
aliquot = unc.ufloat(50, np.sqrt(0.8**2 + 0.3**2))
net_volume = unc.ufloat(3, 3 * np.sqrt(0.006**2 + 0.016**2))
result = (stock_concentration * aliquot) / net_volume
relative_std_of_result = result.std_dev / result.nominal_value
print(f'Relative standard deviation of calibration: {relative_std_of_result:.2%}')

# manual preparation of reaction mixture
aliquot2 = unc.ufloat(250, 250*np.sqrt(0.013**2 + 0.005**2))
net_volume2 = aliquot2 + unc.ufloat(250, 250*np.sqrt(0.013**2 + 0.005**2)) + unc.ufloat(250, 250*np.sqrt(0.013**2 + 0.005**2))
reaction_result = (stock_concentration * aliquot2) / net_volume2
relative_std_of_result = reaction_result.std_dev / reaction_result.nominal_value
print(f'Relative standard deviation of reaction mixture: {relative_std_of_result:.2%}')
# dilution
aliquot3 = unc.ufloat(50, np.sqrt(0.8**2 + 0.3**2))
net_volume3 = unc.ufloat(3, 3 * np.sqrt(0.006**2 + 0.016**2))
diluted_result = (reaction_result * aliquot3) / net_volume3
relative_std_of_result = diluted_result.std_dev / diluted_result.nominal_value
print(f'Relative standard deviation of diluted reaction mixture: {relative_std_of_result:.2%}')

# BP_Ullman
print('BP_Ullman')

powder_weight = unc.ufloat(500, 0.6)
solvent_volume = unc.ufloat(3, 1/20)

aliquot = unc.ufloat(102, np.sqrt(1**2 + 0.3**2))
net_volume = unc.ufloat(20, 0.2)

result = (powder_weight / solvent_volume) * (aliquot / net_volume)

aliquot2 = unc.ufloat(1000, np.sqrt(6**2 + 2**2))
net_volume2 = unc.ufloat(10, 0.1)

result2 = result * (aliquot2 / net_volume2)
relative_std_of_result2 = result2.std_dev / result2.nominal_value
print(f'Relative standard deviation of initial concentrations after dilution: {relative_std_of_result2:.2%}')


