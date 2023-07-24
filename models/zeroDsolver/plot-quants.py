# Plot each column in quants.dat.0
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

quants_file = 'quants.dat.0'
# columns in quants
col_names = [
	'time','Ar','CH4','C2H6','C2H4','C2H2','C3H8',
	'C3H6','O2','H2O','H2O2','CH2O','CO','CO2',
	'H2','HO2','O','OH','CH3','C2H5','C2H3','C3H7',
	'CH3O','CHO','H','S1'
	]

quants_df = pd.read_csv(quants_file, sep="\t", index_col=False, header=None, names=col_names)
Ar_data = np.array(quants_df['Ar'])
CH4_data = np.array(quants_df['CH4'])
C2H6_data = np.array(quants_df['C2H6'])
C2H4_data = np.array(quants_df['C2H4'])
O2_data = np.array(quants_df['O2'])
S1_data = np.array(quants_df['S1'])
time_data = np.array(quants_df['time'])

def plot_species(time, data, pltname):
	end = -1
	plt.plot(time[:end], data[:end])
	plt.xlim(0,max(time[:end]))
	plt.xlabel('Time (s)', fontsize=20)
	plt.ylabel('Concentration $mol/m^3$', fontsize=20)
	plt.savefig(f'{pltname}.png',bbox_inches='tight')
	plt.close()

plot_species(time_data, CH4_data, f'CH4concentration')
plot_species(time_data, O2_data, f'O2concentration')
plot_species(time_data, C2H6_data, f'C2H6concentration')
plot_species(time_data, C2H4_data, f'C2H4concentration')
