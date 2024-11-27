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
green = "#41757F"
blue = "#43599c"

if __name__ == '__main__':

	## Root directory to search for files
	root_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\"
	save_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\Figure 4\\"

	## Load in a dataframe and define the fractions of cells that have
	## phagocytosed or trogocytosed
	df = pd.read_csv(root_dir + "Figure 4\\FcBlock_titration.csv", header=0)

	## Compute the average values for each cell type at each concentration of opsonin
	df['Phago_frac'] = df['Phago'].values/df['Count'].values.astype(int)*100
	df['Trogo_frac'] = df['Trogo'].values/df['Count'].values.astype(int)*100

	# Split dataframe by phago vs trogo, then concatenate (to assist with seaborn plotting)
	p_df = df[['Cell Type', 'Antibody Conc', 'Phago_frac']]
	p_df = p_df.rename(columns = {'Phago_frac': 'Frac'})
	p_df['Behavior'] = 'Phago'

	t_df = df[['Cell Type', 'Antibody Conc', 'Trogo_frac']]
	t_df = t_df.rename(columns = {'Trogo_frac': 'Frac'})
	t_df['Behavior'] = 'Trogo'

	plot_df = p_df

	##### Plot it! ######

	## Figure size needs to be in inches, conversion factor
	mm = 1/25.4

	plotting_parameters = {'data': plot_df,
						   'x' : 'Antibody Conc',
						   'y' : 'Frac',
						   'hue': 'Cell Type',
						   'palette': [blue, orange, green]}

	figs, ax = plt.subplots(1,1, figsize=(65*mm,50*mm))

	ax = sns.barplot(**plotting_parameters, alpha=0.7)
	ax.set_xlabel("antibody solution conc (uM)", fontsize=10)
	ax.set_ylabel("% macrophages", fontsize=10)
	# ax.set_ylim(0,75)
	plt.legend()
	plt.legend([],[], frameon=False)
	plt.tight_layout()

	# Save a .png of the figure
	plt.savefig(save_dir + "figS4_Phagocytosis.pdf")
	plt.show()