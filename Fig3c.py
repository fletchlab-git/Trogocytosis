""" This script imports GUV tension data and plots it for phagocytosis and trogocytosis."""

## Computing and data processing libraries
import os 
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import linregress
# from statsmodels.formula.api import ols

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
blue = "#43599c"
red = "#D24C2B"

if __name__ == '__main__':
	
## Root directory to search for files
	root_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\"
	save_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\Figure 3\\"
	
	## Load in a dataframe and define the fractions of cells that have
	## phagocytosed or trogocytosed
	df_sample = pd.read_csv(root_dir + "Figure 3\\GUVtensions.csv", header=0)
	t_df = df_sample[df_sample['Behavior'] == 'Trogo']['Tension']
	p_df = df_sample[df_sample['Behavior'] == 'Phago']['Tension']

	## Plot it!

	## Set up  plotting and statistical annotation

	## Figure size needs to be in inches, conversion factor
	mm = 1/25.4

	plotting_parameters_box = {'data': df_sample,
						   'x' : 'Behavior',
						   'y' : 'Tension',
						   'order' : ['Phago', 'Trogo'],
						   'orient':'v',
						   'width': 0.5,
						   'linewidth': 1.5,
						   'saturation': 0.3,
						   'palette' : [red, blue]}
	plotting_parameters_swarm = {'data': df_sample,
						   'x' : 'Behavior',
						   'y' : 'Tension',
						   'order' : ['Phago', 'Trogo'],
						   'palette' : [red, blue]}



	figs, ax = plt.subplots(1,1, figsize=(55*mm,45*mm))
	
	ax = sns.swarmplot(**plotting_parameters_swarm)
	ax = sns.boxplot(**plotting_parameters_box)
	ax.set_ylim(0.001, 200)
	ax.set_yscale('log')
	ax.set_xlabel("")
	ax.set_ylabel("tension (mN/m)", fontsize=10)

	# pairs = [('Phago', 'Trogo')]
	# annotator = Annotator(ax, pairs, **plotting_parameters_swarm)
	# annotator.configure(test='t-test_ind', text_format='star', loc='inside')
	# annotator.apply_and_annotate()

	plt.tight_layout()

	# plt.savefig(save_dir + "fig3c v2.pdf")
	plt.show()