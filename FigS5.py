""" This script plots AF647 antiCD47 titration data for three different cell types"""

## Computing and data processing libraries
import os 
import pandas as pd
import numpy as np
from scipy.stats import linregress
from scipy.optimize import minimize
# from statsmodels.formula.api import ols

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
orange = "#d2682b"
green = "#41757F"
blue = "#43599c"

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
	save_dir = "C:\\Users\\caitl\\Documents\\Fletcher lab\\Data\\APMAP collaboration\\Flow data\\06Aug2024\\"

	## Load in a dataframe and define the fractions of cells that have
	## phagocytosed or trogocytosed
	df = pd.read_csv(save_dir + "CD47_titration.csv", header=0)
	df_cal = pd.read_csv(root_dir + "Calibration_values.csv", header=0)
	# df_tension = pd.read_csv(root_dir + "Figure 5\\CellTension_Susp.csv")

	## Using the calibration data from AF647 beads, compute the MESF for each AF647 intensity
	# df_titration = MESFcalc(df_cal, df).drop(columns=['Conc_MESFavg'])

	## Compute the number of antibodies per square micron using the average size of the cell
	## I already divided the MESF by the number of fluorophores per antibody, I just need to 
	## divide by the average surface area
	# avg_SA = []
	# cell_types = ['HL60', 'Raji B', 'Jurkat']
	# for i in cell_types:
	# 	avg_SA = 4*np.pi*(df_tension[df_tension['Cell Type']==i]['Radius'].mean())
	# 	df_titration['Surface Density'] = df_titration['Intensity_MESFavg'] / avg_SA

	##### Plot it! ######

	## Figure size needs to be in inches, conversion factor
	mm = 1/25.4

	plotting_parameters = {'data': df,
						   'x' : 'Conc',
						   'y' : 'Intensity',
						   'hue': 'Cell Type',
						   'palette': [blue, orange, green]}

	figs, ax = plt.subplots(1,1, figsize=(75*mm,60*mm))

	ax = sns.barplot(**plotting_parameters)
	ax.set_ylabel("AF647 intensity (a.u.)", fontsize=12)
	ax.set_xlabel("antibody solution conc (uM)", fontsize=12)
	# ax.set_ylim(0,75)
	plt.legend()
	# plt.legend([],[], frameon=False)
	plt.tight_layout()

	# Save a .png of the figure
	plt.savefig(save_dir + "CD47_staining.png")
	plt.show()