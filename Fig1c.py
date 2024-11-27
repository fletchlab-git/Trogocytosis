""" This script imports processed data from FloJo to plot the fraction of phagocytosis versus trogocytosis
of target cells opsonized with anti-CD47. This script is specifically for plotting Figure 1c."""

## Computing and data processing libraries
import os 
import pandas as pd
import numpy as np


## Plotting libraries
import seaborn as sns
import matplotlib.pyplot as plt


## Set global matplotlib parameters
plt.rcParams["font.family"] = "Calibri"
plt.rcParams["xtick.labelsize"] = 10 
plt.rcParams["ytick.labelsize"] = 10
plt.rcParams["axes.linewidth"] = 0.6
plt.rcParams["xtick.major.width"] = 1.0
plt.rcParams["ytick.major.width"] = 1.0
plt.rcParams["errorbar.capsize"] = 1.0

## Set some colors as well, used in the plots
## throughout.
sage = "#96c8a2"
light_sage = "#d3e1d0"
yellow = "#EBBC4E"
dark_yellow = "#d49c18"
red = "#D24C2B"
orange = "#d2682b"
blue = "#43599c"
dark_blue = "#273253"
green = "#41757F"
dark_green = "#3C6571"


if __name__ == '__main__':
	
## Root directory to search for files
	root_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\"
	save_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\Figure 1\\"
	
	## Load in a dataframe and filter the data by cell type and concentration of antibody
	df_J774 = pd.read_csv(root_dir + "Figure 1\\J774_FlowResults.csv", header=0)
	df_J774 = df_J774[(df_J774['Cell Type'] == 'Jurkat') & (df_J774['Conc']==1.0)].replace('Jurkat', 'J774')

	df_BMDM = pd.read_csv(root_dir + "Figure 1\\BMDMs_FlowResults.csv", header=0)
	df_BMDM = df_BMDM[(df_BMDM['Cell Type'] == 'Jurkat') & (df_BMDM['Conc']==1.0)].replace('Jurkat', 'BMDM')

	df_RAW = pd.read_csv(root_dir + "Figure 1\\CD47titration_FlowResults.csv", header=0)
	df_RAW = df_RAW[(df_RAW['Cell Type'] == 'Jurkat') & (df_RAW['Conc']==1.0)].replace('Jurkat', 'RAW')

	df_LPS = pd.read_csv(root_dir + "Figure 1\\LPS_FlowResults.csv", header=0)

	## Combine all the cell types into one dataframe
	df = pd.concat([df_J774, df_BMDM, df_RAW, df_LPS])

	## Compute the average values for each cell type at each concentration of opsonin
	df['Phago_frac'] = df['Phago'].values/df['Count'].values.astype(int)*100
	df['Trogo_frac'] = df['Trogo'].values/df['Count'].values.astype(int)*100

	# Split dataframe by phago vs trogo, then concatenate (to assist with seaborn plotting)
	p_df = df[['Cell Type', 'Phago_frac']]
	p_df = p_df.rename(columns = {'Phago_frac': 'Frac'})
	p_df['Behavior'] = 'Phago'

	t_df = df[['Cell Type', 'Trogo_frac']]
	t_df = t_df.rename(columns = {'Trogo_frac': 'Frac'})
	t_df['Behavior'] = 'Trogo'

	plot_df = pd.concat([p_df, t_df])
	
	##### Plot it! ######

	## Figure size needs to be in inches, conversion factor
	mm = 1/25.4

	## Set up  plotting and statistical annotation

	plotting_parameters = {'data': plot_df,
						   'x' : 'Cell Type',
						   'y' : 'Frac',
						   'order' : ['RAW', 'RAW + LPS', 'J774', 'BMDM'],
						   'hue' : 'Behavior',
						   'hue_order': ['Phago', 'Trogo'],
						   'errwidth':1.5,
						   'palette': [red, blue]}

	figs, ax = plt.subplots(1,1, figsize=(75*mm, 55*mm))

	ax = sns.barplot(**plotting_parameters)
	ax.set_xlabel("")
	ax.set_ylabel("% of macrophages", fontsize=12)
	plt.tight_layout()
	plt.legend([],[], frameon=False)


	# Save a .png of the figure
	plt.savefig(save_dir + 'fig1c_v2.pdf')
	plt.show()
