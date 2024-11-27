## Computing and data processing libraries
import os 
import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import linregress
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
	save_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\Figure 2\\"

	## Load in a dataframe and define the fractions of cells that have
	## phagocytosed or trogocytosed
	## Load in a dataframe and define the fractions of cells that have
	## phagocytosed or trogocytosed
	df_biotin = pd.read_csv(root_dir + "Figure 1\\BiotinTitration_FlowResults_Susp.csv", header=0)
	df_cd47 = pd.read_csv(root_dir + "Figure 1\\CD47Titration_FlowResults_Susp.csv", header=0)
	df_cal = pd.read_csv(root_dir + "Calibration_values.csv", header=0)
	df_size = pd.read_csv(root_dir + "Figure 5\\CellTension_susp.csv", header=0)
	df_tension = pd.read_csv(root_dir + "Figure 5\\CellTension_Susp.csv")
	df_tension = df_tension[df_tension["Treatment"]=="WT"]

	## Using the calibration data from AF647 beads, compute the MESF for each AF647 intensity
	df_combine = pd.concat([df_biotin, df_cd47])
	df_flow = MESFcalc(df_cal, df_combine).drop(columns=['Conc_MESFavg'])

	## Compute the number of antibodies per square micron using the average size of the cell
	## I already divided the MESF by the number of fluorophores per antibody, I just need to 
	## divide by the average surface area
	avg_SA = []
	cell_types = ['HL60', 'Raji B', 'Jurkat']
	for i in cell_types:
		avg_SA = 4*np.pi*(df_size[df_size['Cell Type']==i]['Radius'].mean())
		df_flow['Surface Density'] = df_flow['Intensity_MESFavg'] / avg_SA

	## Compute the average values for each cell type at each concentration of opsonin
	df_flow['Phago_frac'] = df_flow['Phago'].values/df_flow['Count'].values.astype(int)*100
	df_flow['Trogo_frac'] = df_flow['Trogo'].values/df_flow['Count'].values.astype(int)*100


	## Select values where the 'intensity' falls between 500-1000
	df_slice = df_flow[(df_flow['Surface Density'] > 300) & (df_flow['Surface Density'] < 500)]
	df_tension = df_tension[df_tension["Treatment"]=="WT"]

	## Using the calibration data from AF647 beads, compute the MESF for each AF647 intensity
	df_combine = pd.concat([df_biotin, df_cd47])
	df_flow = MESFcalc(df_cal, df_combine).drop(columns=['Conc_MESFavg'])

	## Compute the number of antibodies per square micron using the average size of the cell
	## I already divided the MESF by the number of fluorophores per antibody, I just need to 
	## divide by the average surface area
	avg_SA = []
	cell_types = ['HL60', 'Raji B', 'Jurkat']
	for i in cell_types:
		avg_SA = 4*np.pi*(df_size[df_size['Cell Type']==i]['Radius'].mean())
		df_flow['Surface Density'] = df_flow['Intensity_MESFavg'] / avg_SA

	## Compute the average values for each cell type at each concentration of opsonin
	df_flow['Phago_frac'] = df_flow['Phago'].values/df_flow['Count'].values.astype(int)*100
	df_flow['Trogo_frac'] = df_flow['Trogo'].values/df_flow['Count'].values.astype(int)*100


	## Compute the averages of tension and trogo fraction for each cell type,
	## then concatenate the dataframes
	trogo_avg = df_slice.filter(items=['Trogo_frac', 'Cell Type']).groupby(['Cell Type']).mean()
	trogo_std = df_slice.filter(items=['Trogo_frac', 'Cell Type']).groupby(['Cell Type']).std()
	trogo_avg = trogo_avg[['Trogo_frac']]

	tension_avg = df_tension.filter(items=['Tension', 'Cell Type']).groupby(['Cell Type']).mean()
	tension_std = df_tension.filter(items=['Tension', 'Cell Type']).groupby(['Cell Type']).std()

	plot_df = trogo_avg.join(tension_avg).reset_index()
	error_df = trogo_std.join(tension_std).reset_index()

	## Compute a Pearson's correlation coefficient
	res = stats.pearsonr(plot_df['Trogo_frac'].values, plot_df['Tension'].values)
	p_cor = round(res[0], 2)
	print(p_cor)

	##### Plot it! #####

	## Figure size needs to be in inches, conversion factor
	mm = 1/25.4

	figs, ax = plt.subplots(1,1, figsize=(55*mm,55*mm))

	ax.plot(plot_df['Tension'],plot_df['Trogo_frac'], marker='o', markersize=5, color='grey', linestyle="None")
	ax.errorbar(plot_df['Tension'], plot_df['Trogo_frac'], xerr = tension_std['Tension'] , yerr=trogo_std['Trogo_frac'], linestyle="None", color='grey', elinewidth=1, capsize=0.0)
	ax.set_ylabel('trogocytic efficiency (%)', fontsize=10)
	ax.set_xlabel('tension (mN/m)', fontsize=10)

	# plt.tight_layout()
	# plt.savefig(save_dir + "Fig2f.pdf")
	plt.show()