""" This script imports processed data from micropipette measurements to plot the tension
of cells with gluteraldehyde treatment for figure 4."""

## Computing and data processing libraries
import os 
import pandas as pd
import numpy as np
import scipy.stats as stats
# from statsmodels.stats.multicomp import pairwise_tukeyhsd


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
orange = "#d2682b"
light_orange = "#fc8947"
dark_orange = "#c45210"
blue = "#43599c"
light_blue = "#6e91fa"
dark_blue = "#2e4078"


if __name__ == '__main__':

	## Root directory to search for files
	root_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\"
	save_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\Figure 4\\"

	## Load in the data
	df_tension = pd.read_csv(root_dir + "Figure 5\\CellTension_Susp.csv")

	##### Plot it! ######

	## Set up  plotting and statistical annotation
	## Perform a one-way ANOVA
	data = df_tension[df_tension['Cell Type']=='HL60']
	# a = data[data['Treatment']=='WT']['Tension'].values
	# b = data[data['Treatment']=='0.002']['Tension'].values
	# c = data[data['Treatment']=='0.0025']['Tension'].values
	# anova_res = stats.f_oneway(a, b, c)
	# print(anova_res)

	# ## Perform Tukey's pairwise test
	# tukey = pairwise_tukeyhsd(endog=data['Tension'], groups=data['Treatment'], alpha=0.05)

	# ## Create a dataframe with the corresponding results from the Tukey HSD
	# tukey_df = pd.DataFrame(data=tukey._results_table.data[1:], columns=tukey._results_table.data[0])
	# print(tukey_df)


	## Figure size needs to be in inches, conversion factor
	mm = 1/25.4

	plotting_parameters = {'data': data,
						   'x' : 'Treatment',
						   'y' : 'Tension',
						   'orient': 'v',
						   'palette': [light_orange, orange, dark_orange]}

	figs, ax = plt.subplots(1,1, figsize=(50*mm,50*mm))

	ax = sns.boxplot(**plotting_parameters, boxprops=dict(alpha=.8))
	ax = sns.swarmplot(**plotting_parameters)

	# pairs = [(tukey_df['group1'][0], tukey_df['group2'][0]),
	# 		 (tukey_df['group1'][1], tukey_df['group2'][1]),
	# 		 (tukey_df['group1'][2], tukey_df['group2'][2])]

	# annotator = Annotator(ax, pairs, **plotting_parameters)
	# annotator.set_pvalues(tukey_df['p-adj'].values)
	# annotator.annotate()

	ax.set_ylabel('tension (mN/m)', fontsize=10)
	ax.set_xlabel('', fontsize=10)
	ax.set_ylim(-0.2,4.3)
	ax.tick_params(axis='x', labelrotation=45)
	plt.tight_layout()

	# Save a .pdf of the figure
	plt.savefig(save_dir + "fig4a_HL60s.pdf")
	plt.show()