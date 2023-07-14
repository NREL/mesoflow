# Plot each column in quants.dat.0
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

quants_file = 'quants.dat.0'
# columns in quants
col_names = ['time','H2','H2ADS1','S1']

quants_df = pd.read_csv(quants_file, sep="\t", index_col=False, header=None, names=col_names)
H2_data = np.array(quants_df['H2'])
H2ADS1_data = np.array(quants_df['H2ADS1'])
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

plot_species(time_data, H2_data, f'H2concentration')
plot_species(time_data, H2ADS1_data, f'H2ADS1concentration')
plot_species(time_data, S1_data, f'S1concentration')
