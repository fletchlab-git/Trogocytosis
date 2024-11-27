""" This script imports processed data from FloJo to plot the fraction of phagocytosis versus trogocytosis
of target cells opsonized with anti-CD47 or anti-biotin. This script is specifically for plotting Figures 1b-c."""

## Computing and data processing libraries
import os 
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import linregress

## Plotting libraries
import seaborn as sns
from statannotations.Annotator import Annotator
import matplotlib.pyplot as plt


## Set global matplotlib parameters
plt.rcParams["font.family"] = "Calibri"
plt.rcParams["xtick.labelsize"] = 10 
plt.rcParams["ytick.labelsize"] = 10
plt.rcParams["axes.linewidth"] = 0.6
plt.rcParams["xtick.major.width"] = 1.0
plt.rcParams["ytick.major.width"] = 1.0
plt.rcParams["errorbar.capsize"] = 0.5

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

def MESFcalc(df_cal, df_sample):
	''' This function fits AF647 intensity values of samples to the AF647 calibration
	curve and returns an array of MESF values for the sample. The MESF value is scaled by
	3-7x as the AF647 anti-CD47 antibody is labeled with 3-7 individual fluorophores.'''

	MESF_cal = df_cal['MESF'].values
	AF647_med = df_cal['AF647_med'].values

	# Fit the data to a linear curve_fit
	fit = linregress(MESF_cal, AF647_med)

	# Create a sub dataframe of the sample data with just concentration and intensity values
	sub_df = df_sample[['Conc', 'Intensity']]
	print(sub_df)

	## Find the MESF values for sample intensities
	for (conc, sample) in sub_df.items():

		## Compute MESF from the fit. The values span between 3-7x less than the recorded value
		upper = ((sample - fit[1]) / fit[0])/3
		lower = ((sample - fit[1]) / fit[0])/7
		df_sample[conc + '_' + 'MESFavg'] = (upper + lower)/2


	return df_sample

if __name__ == '__main__':
	
## Root directory to search for files
	root_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\"
	save_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\Figure 1\\"
	
	## Load in a dataframe and define the fractions of cells that have
	## phagocytosed or trogocytosed
	df_cd47 = pd.read_csv(root_dir + "Figure 1\\CD47titration_FlowResults_Susp.csv", header=0)
	df_biotin = pd.read_csv(root_dir + "Figure 1\\BiotinTitration_FlowResults_Susp.csv", header=0)
	df_cal = pd.read_csv(root_dir + "Calibration_values.csv", header=0)
	df_size = pd.read_csv(root_dir + "Figure 5\\CellTension_susp.csv", header=0)

	## Using the calibration data from AF647 beads, compute the MESF for each AF647 intensity
	df_c = MESFcalc(df_cal, df_cd47).drop(columns=['Conc_MESFavg'])
	df_b = MESFcalc(df_cal, df_biotin).drop(columns=['Conc_MESFavg'])

	## Compute the number of antibodies per square micron using the average size of the cell
	## I already divided the MESF by the number of fluorophores per antibody, I just need to 
	## divide by the average surface area
	avg_SA = []
	cell_types = ['HL60', 'Raji B', 'Jurkat']
	for i in cell_types:
		avg_SA = 4*np.pi*(df_size[df_size['Cell Type']==i]['Radius'].mean())
		df_c['Surface Density'] = df_c['Intensity_MESFavg'] / avg_SA
		df_b['Surface Density'] = df_b['Intensity_MESFavg'] / avg_SA

	## Compute the average values for each cell type at each concentration of opsonin
	df_c['Phago_frac'] = df_c['Phago'].values/df_c['Count'].values.astype(int)*100
	df_c['Trogo_frac'] = df_c['Trogo'].values/df_c['Count'].values.astype(int)*100

	df_b['Phago_frac'] = df_b['Phago'].values/df_b['Count'].values.astype(int)*100
	df_b['Trogo_frac'] = df_b['Trogo'].values/df_b['Count'].values.astype(int)*100

	## Create a sub dataframe containing only values for concentration 1.0
	df_c = df_c[df_c['Conc']==1.0]
	df_b = df_b[df_b['Conc']==1.0]

	# Split dataframe by phago vs trogo, then concatenate (to assist with seaborn plotting)
	
	# CD47

	t_df_c = df_c[['Cell Type', 'Surface Density', 'Trogo_frac']]
	t_df_c = t_df_c.rename(columns = {'Trogo_frac': 'Frac'})
	t_df_c['Behavior'] = 'Trogo'
	t_df_c['Opsonin'] = 'CD47'

	p_df_c = df_c[['Cell Type', 'Surface Density', 'Phago_frac']]
	p_df_c = p_df_c.rename(columns = {'Phago_frac': 'Frac'})
	p_df_c['Behavior'] = 'Phago'
	p_df_c['Opsonin'] = 'CD47'

	# Biotin

	t_df_b = df_b[['Cell Type', 'Surface Density', 'Trogo_frac']]
	t_df_b = t_df_b.rename(columns = {'Trogo_frac': 'Frac'})
	t_df_b['Behavior'] = 'Trogo'
	t_df_b['Opsonin'] = 'Biotin'

	p_df_b = df_b[['Cell Type', 'Surface Density', 'Phago_frac']]
	p_df_b = p_df_b.rename(columns = {'Phago_frac': 'Frac'})
	p_df_b['Behavior'] = 'Phago'
	p_df_b['Opsonin'] = 'Biotin'

	plot_df_t = pd.concat([t_df_c, t_df_b], axis=0).sort_values('Surface Density')
	plot_df_p = pd.concat([p_df_c, p_df_b], axis=0).sort_values('Surface Density')

	plot_df = pd.concat([plot_df_t, plot_df_p], axis=0)

	## For this figure, we're just plotting values for HL60s
	plot_df = plot_df[plot_df['Cell Type']=='Jurkat']


	##### Plot it! ######

	## Set up  plotting and statistical annotation

	## Figure size needs to be in inches, conversion factor
	mm = 1/25.4

	plotting_parameters = {'data': plot_df,
						   'x' : 'Opsonin',
						   'y' : 'Frac',
						   'order' : ['CD47', 'Biotin'],
						   'hue' : 'Behavior',
						   'hue_order': ['Phago', 'Trogo'],
						   'ci': 'sd',
						   'errwidth': 1.5, 
						   'palette': [red, blue]}


	figs, ax = plt.subplots(1,1, figsize=(75*mm,55*mm))

	# plt.style.use('dark_background')

	ax = sns.barplot(**plotting_parameters)
	# pairs = [(('CD47', 'Phago'), ('CD47', 'Trogo')), (('Biotin', 'Phago'), ('Biotin', 'Trogo'))]
	# annotator = Annotator(ax, pairs, **plotting_parameters)
	# annotator.configure(test='t-test_ind', text_format='star', loc='inside')
	# annotator.apply_and_annotate()

	ax.set_xlabel("")
	ax.set_ylabel("% macrophages", fontsize=12)
	plt.tight_layout()
	plt.legend([],[], frameon=False)

	# Save a .png of the figure
	plt.savefig(save_dir + "fig1b_v2_JurkatOnly.pdf")
	plt.show()
