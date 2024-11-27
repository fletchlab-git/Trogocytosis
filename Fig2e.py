""" This script imports processed data from micropipette measurements to plot the tension
of cells."""

## Computing and data processing libraries
import os 
import pandas as pd
import numpy as np
# import scipy.stats as stats

## Plotting libraries
import seaborn as sns
# from statannotations.Annotator import Annotator
import matplotlib.pyplot as plt


## Set global matplotlib parameters
plt.rcParams["font.family"] = "Calibri"
plt.rcParams["xtick.labelsize"] = 10 
plt.rcParams["ytick.labelsize"] = 10
plt.rcParams["axes.linewidth"] = 1.0
plt.rcParams["xtick.major.width"] = 1.0
plt.rcParams["ytick.major.width"] = 1.0
plt.rcParams["errorbar.capsize"] = 1.0

## Set some colors as well, used in the plots
## throughout.
sage = "#96c8a2"
orange = "#d2682b"
blue = "#43599c"
green = "#41757F"


if __name__ == '__main__':

	## Root directory to search for files
	root_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\"
	save_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\Figure 2\\"

	## Load in the data
	df_tension = pd.read_csv(root_dir + "Figure 5\\CellTension_Susp.csv")
	print(df_tension[df_tension["Cell Type"]=="Jurkat"]["Tension"].mean())

	##### Plot it! ######

	## Set up  plotting and statistical annotation

	## Figure size needs to be in inches, conversion factor
	mm = 1/25.4

	plotting_parameters = {'data': df_tension[df_tension['Treatment'] == 'WT'],
						   'x' : 'Cell Type',
						   'y' : 'Tension',
						   'order' : ['HL60', 'Raji B', 'Jurkat'],
						   'orient': 'v',
						   'palette': [sage, orange, blue]}

	figs, ax = plt.subplots(1,1, figsize=(55*mm,55*mm))

	ax = sns.boxplot(**plotting_parameters, boxprops=dict(alpha=.8))
	ax = sns.swarmplot(**plotting_parameters)

	# pairs = [('Jurkat', 'HL60'), ('Jurkat', 'Raji B'), ('HL60', 'Raji B')]
	# annotator = Annotator(ax, pairs, **plotting_parameters)
	# annotator.configure(test='t-test_ind', text_format='star', loc='inside')
	# annotator.apply_and_annotate()

	ax.set_ylabel('tension (mN/m)', fontsize=10)
	ax.set_xlabel('', fontsize=10)
	ax.set_ylim(-0.1,2.0)
	# plt.tight_layout()

	# Save a .pdf of the figure
	plt.savefig(save_dir + "fig2a_susp.pdf")
	plt.show()