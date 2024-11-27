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
blue = "#43599c"
red = "#D24C2B"

if __name__ == '__main__':

	## Root directory to search for files
	root_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\"
	save_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\Figure 4\\"

	## Load in a dataframe and define the fractions of cells that have
	## phagocytosed or trogocytosed
	df = pd.read_csv(root_dir + "Figure 4\\Glutaraldehyde_titration.csv", header=0)

	## Drop the 0.001 concentration because we don't have tension measurements
	df = df.drop(df[df['Drug Conc'] == '0.001%'].index)

	## Compute the average values for each cell type at each concentration of opsonin
	df['Phago_frac'] = df['Phago'].values/df['Count'].values.astype(int)*100
	df['Trogo_frac'] = df['Trogo'].values/df['Count'].values.astype(int)*100

	## Create a sub dataframes for each cell type
	sub_df_HL60 = df[df['Cell Type']=='HL60']
	sub_df_jurkat = df[df['Cell Type']=='Jurkat']

	# Split dataframe by phago vs trogo, then concatenate (to assist with seaborn plotting)
	p_df = sub_df_jurkat[['Cell Type', 'Drug Conc', 'Phago_frac']]
	p_df = p_df.rename(columns = {'Phago_frac': 'Frac'})
	p_df['Behavior'] = 'Phago'

	t_df = sub_df_jurkat[['Cell Type', 'Drug Conc', 'Trogo_frac']]
	t_df = t_df.rename(columns = {'Trogo_frac': 'Frac'})
	t_df['Behavior'] = 'Trogo'

	plot_df = t_df
	
	
	##### Plot it! ######

	## Set up  plotting and statistical annotation
	## Perform a one-way ANOVA
	# data = plot_df
	# a = data[data['Drug Conc']=='WT']['Frac'].values
	# b = data[data['Drug Conc']=='0.002']['Frac'].values
	# c = data[data['Drug Conc']=='0.0025']['Frac'].values
	# anova_res = stats.f_oneway(a, b, c)
	# print(anova_res)

	# ## Perform Tukey's pairwise test
	# tukey = pairwise_tukeyhsd(endog=data['Frac'], groups=data['Drug Conc'], alpha=0.1)

	# ## Create a dataframe with the corresponding results from the Tukey HSD
	# tukey_df = pd.DataFrame(data=tukey._results_table.data[1:], columns=tukey._results_table.data[0])
	# print(tukey_df)


	## Figure size needs to be in inches, conversion factor
	mm = 1/25.4

	plotting_parameters = {'data': plot_df,
						   'x' : 'Drug Conc',
						   'y' : 'Frac',
						   'palette': [blue]}

	figs, ax = plt.subplots(1,1, figsize=(50*mm,50*mm))

	ax = sns.barplot(**plotting_parameters)

	# pairs = [(tukey_df['group1'][0], tukey_df['group2'][0]),
	# 		 (tukey_df['group1'][1], tukey_df['group2'][1]),
	# 		 (tukey_df['group1'][2], tukey_df['group2'][2])]

	# annotator = Annotator(ax, pairs, **plotting_parameters)
	# annotator.set_pvalues(tukey_df['p-adj'].values)
	# annotator.annotate()

	ax.set_xlabel("", fontsize=10)
	ax.set_ylabel("% macrophages", fontsize=10)
	ax.set_ylim(0,90)
	ax.tick_params(axis='x', labelrotation=40)
	plt.legend([],[], frameon=False)
	plt.tight_layout()

	# Save a .png of the figure
	plt.savefig(save_dir + "fig4b_jurkat_trogocytosis.pdf")
	plt.show()
