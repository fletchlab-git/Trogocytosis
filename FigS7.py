""" This script plots antibody titration data obtained from flow."""

### Computing and data processing libraries
import os 
import pandas as pd
import numpy as np
from scipy.interpolate import CubicSpline

## For regression
from scipy.stats import linregress
from scipy.optimize import minimize

## Plotting libraries
import seaborn as sns
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
light_blue = "#6e91fa"
dark_blue = "#2e4078"

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
	save_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\Figure 2\\"

	## Load in a dataframe and define the fractions of cells that have
	## phagocytosed or trogocytosed
	df_cd47 = pd.read_csv(root_dir + "Figure 1\\CD47titration_FlowResults_Susp.csv", header=0)
	df_biotin = pd.read_csv(root_dir + "Figure 1\\BiotinTitration_FlowResults_Susp.csv", header=0)
	df_glut = pd.read_csv(root_dir + "Figure 5\\AntibodyTitration_treatments_FlowResults_JURKATS.csv", header=0)
	df_cal = pd.read_csv(root_dir + "Calibration_values.csv", header=0)
	df_size = pd.read_csv(root_dir + "Figure 5\\CellTension_susp.csv", header=0)

	## Using the calibration data from AF647 beads, compute the MESF for each AF647 intensity
	df_c = MESFcalc(df_cal, df_cd47).drop(columns=['Conc_MESFavg'])
	df_b = MESFcalc(df_cal, df_biotin).drop(columns=['Conc_MESFavg'])
	df_g = MESFcalc(df_cal, df_glut).drop(columns=['Conc_MESFavg'])

	## Compute the number of antibodies per square micron using the average size of the cell
	## I already divided the MESF by the number of fluorophores per antibody, I just need to 
	## divide by the average surface area
	avg_SA = []
	cell_types = ['HL60', 'Raji B', 'Jurkat']
	for i in cell_types:
		avg_SA = 4*np.pi*(df_size[df_size['Cell Type']==i]['Radius'].mean())
		df_c['Surface Density'] = df_c['Intensity_MESFavg'] / avg_SA
		df_b['Surface Density'] = df_b['Intensity_MESFavg'] / avg_SA
		df_g['Surface Density'] = df_g['Intensity_MESFavg'] / avg_SA

	## Compute the average values for each cell type at each concentration of opsonin
	df_c['Phago_frac'] = df_c['Phago'].values/df_c['Count'].values.astype(int)*100
	df_c['Trogo_frac'] = df_c['Trogo'].values/df_c['Count'].values.astype(int)*100

	df_b['Phago_frac'] = df_b['Phago'].values/df_b['Count'].values.astype(int)*100
	df_b['Trogo_frac'] = df_b['Trogo'].values/df_b['Count'].values.astype(int)*100

	df_g['Phago_frac'] = df_g['Phago'].values/df_g['Count'].values.astype(int)*100
	df_g['Trogo_frac'] = df_g['Trogo'].values/df_g['Count'].values.astype(int)*100

	# Split dataframe by phago vs trogo, then concatenate (to assist with seaborn plotting)
	
	# CD47

	t_df_c = df_c[['Cell Type', 'Surface Density', 'Trogo_frac']]
	t_df_c = t_df_c.rename(columns = {'Trogo_frac': 'Trogo_Frac'})

	p_df_c = df_c[['Cell Type', 'Surface Density', 'Phago_frac']]
	p_df_c = p_df_c.rename(columns = {'Phago_frac': 'Phago_Frac'})

	# Biotin

	t_df_b = df_b[['Cell Type', 'Surface Density', 'Trogo_frac']]
	t_df_b = t_df_b.rename(columns = {'Trogo_frac': 'Trogo_Frac'})

	p_df_b = df_b[['Cell Type', 'Surface Density', 'Phago_frac']]
	p_df_b = p_df_b.rename(columns = {'Phago_frac': 'Phago_Frac'})

	# Glutaraldehyde

	t_df_g = df_g[['Treatment', 'Surface Density', 'Trogo_frac']]
	t_df_g = t_df_g.rename(columns = {'Trogo_frac': 'Trogo_Frac'}).sort_values('Surface Density')

	p_df_g = df_g[['Treatment', 'Surface Density', 'Phago_frac']]
	p_df_g = p_df_g.rename(columns = {'Phago_frac': 'Phago_Frac'}).sort_values('Surface Density')

	plot_df_t = pd.concat([t_df_c, t_df_b], axis=0).sort_values('Surface Density')
	plot_df_p = pd.concat([p_df_c, p_df_b], axis=0).sort_values('Surface Density')

	##### Plot it! ######

	## Set up  plotting and statistical annotation

	## Figure size needs to be in inches, conversion factor
	mm = 1/25.4

	plotting_parameters_t1 = {'data': plot_df_t,
						   'x' : 'Surface Density',
						   'y' : 'Trogo_Frac',
						   'hue' : 'Cell Type',
						   'hue_order': ['HL60', 'Raji B', 'Jurkat'],
						   'palette': [orange, sage, blue]}

	plotting_parameters_t2 = {'data': plot_df_t[plot_df_t['Cell Type']=='Jurkat'],
						   'x' : 'Surface Density',
						   'y' : 'Trogo_Frac',
						   'palette': ['k']}

	plotting_parameters_g = {'data': t_df_g,
						   'x' : 'Surface Density',
						   'y' : 'Trogo_Frac',
						   'hue' : 'Treatment',
						   'hue_order': ['WT', '0.002% glut', '0.0025% glut'],
						   'palette': ['grey', light_blue, blue]}


	figs, axes = plt.subplots(1,2, figsize=(190*mm,60*mm))

	# ax = sns.lineplot(**plotting_parameters_t, alpha=0.5,  lw=4)
	sns.scatterplot(**plotting_parameters_t1, s=60, ax=axes[0])
	sns.scatterplot(**plotting_parameters_t2, s=60, ax=axes[1])
	sns.scatterplot(**plotting_parameters_g, s=60, ax=axes[1])
	axes[0].set_xscale('log')
	axes[1].set_xscale('log')
	axes[0].set_xlabel('# antibodies / '+ r'$\mu m^2$', fontsize=10)
	axes[0].set_ylabel("% trogocytosis", fontsize=10)
	axes[1].set_xlabel('# antibodies / '+ r'$\mu m^2$', fontsize=10)
	axes[1].set_ylabel("% trogocytosis", fontsize=10)
	plt.tight_layout()
	# plt.legend([],[], frameon=False)

	# Save a .png of the figure
	plt.savefig(save_dir + "FigS8.jpg")
	plt.show()
