""" This script imports processed data from FloJo to plot the fraction of phagocytosis versus trogocytosis
of GUVs opsonized with anti-biotin and osmotically manipulated. This script is specifically for plotting Figure 3b."""
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
	df = pd.read_csv(root_dir + "Figure 3\\GUV_flowResults.csv", header=0)

	## Compute the average values for each cell type at each concentration of opsonin
	df['Phago_frac'] = df['Phago'].values/df['Count'].values.astype(int)*100
	df['Trogo_frac'] = df['Trogo'].values/df['Count'].values.astype(int)*100

	# Split dataframe by phago vs trogo, then concatenate (to assist with seaborn plotting)
	p_df = df[['Osmolarity', 'Phago_frac']]
	p_df = p_df.rename(columns = {'Phago_frac': 'Frac'})
	p_df['Behavior'] = 'Phago'

	t_df = df[['Osmolarity', 'Trogo_frac']]
	t_df = t_df.rename(columns = {'Trogo_frac': 'Frac'})
	t_df['Behavior'] = 'Trogo'

	plot_df = pd.concat([p_df, t_df], axis=0)
	plot_df = plot_df[plot_df['Osmolarity']== 'isotonic']
	

	##### Plot it! ######

	## Set up  plotting and statistical annotation

	## Figure size needs to be in inches, conversion factor
	mm = 1/25.4

	plotting_parameters = {'data': plot_df,
						   'x' : 'Behavior',
						   'y' : 'Frac',
						   'order' : ['Phago', 'Trogo'],
						   # 'hue' : 'Behavior',
						   # 'hue_order': ['Phago', 'Trogo'],
						   'palette': [red, blue]}


	figs, ax = plt.subplots(1,1, figsize=(50*mm,45*mm))

	ax = sns.barplot(**plotting_parameters)
	# pairs = [('Phago', 'Trogo')]
	# annotator = Annotator(ax, pairs, **plotting_parameters)
	# annotator.configure(test='t-test_ind', text_format='star', loc='inside')
	# annotator.apply_and_annotate()

	ax.set_xlabel("")
	ax.set_ylabel("% macrophages", fontsize=10)
	plt.tight_layout()
	plt.legend([],[], frameon=False)

	# Save a .png of the figure
	plt.savefig(save_dir + "fig3b.pdf")
	plt.show()