""" This script fits a sigmoid to antibody titration data for cells treated with different
concentrations of gluteraldehyde. Data for Figure 4c."""

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
light_orange = "#fc8947"
dark_orange = "#c45210"
blue = "#43599c"
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

def Sigmoid(x, theta):

	"""This function computes a sigmoid"""

	## a = asymptote 
	## c = inflection point of the curve
	## d = rate of the sigmoid decay

	a,c,d = theta

	return a * (1 - 1 / (1 + np.exp(-(x - c) / d)))
	
def L2_cost_function(dataset, theta):
	'''This function computes a weighted least squares loss (or cost) function
	for a fit to data. The loss function has the form L2 = sum((y_true - y_pred)**2)'''

	y_pred = Sigmoid(dataset[:,0], theta) 
	return (((dataset[:,1] - y_pred)**2).sum())

def Fit(dataset, theta0):
	'''This function fits the data to a model and minimizes the cost function to 
	give the best fit parameters for the model to the data.'''

	f = lambda x: L2_cost_function(dataset, x)
	result = minimize(f, x0=theta0)
	return result

def BootStrap(ModelResults, dataset):
	''' Compute the confidence intervals for a particular model of the dataset
	using bootstrap resampling. Here, random samples are drawn from the results
	of the "Fit" function. Each of these samples is taken from the parameters and 
	covariance matrix of the fit. Each of these samples is then fit to a line. 
	The percentiles for the sample lines are computed to give confidence intervals.'''

	## Draw random samples of parameters from the fit. 
	#  ModelResults.x  = parameters
	#  ModelResults.hess_inv = covariance matrix
	p_samples = np.random.multivariate_normal(mean=ModelResults.x, cov=ModelResults.hess_inv, size=(5000,))

	## Make a series of x values that cover 0.8 of the minimum value of X in the dataset
	## and 1.15 of the maximum value of x in the dataset
	x = np.linspace(0.8*dataset[:,0].min(), 1.15*dataset[:,0].max(), 100)

	## Fit lines to the random samples
	samples = np.zeros((5000,100))
	for i, p_i in enumerate(p_samples):
		samples[i] = Sigmoid(x, p_i)


	## Compute the percentiles
	low = np.percentile(samples, 10, axis=0)
	high = np.percentile(samples, 90, axis=0)
	mid = Sigmoid(x, np.array([ModelResults.x[0], ModelResults.x[1], ModelResults.x[2]]))


	return low, high, mid, x

if __name__ == '__main__':

	## Root directory to search for files
	root_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\"
	save_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\Figure 5\\"

	## Load in a dataframe and define the fractions of cells that have
	## phagocytosed or trogocytosed
	df_antibody = pd.read_csv(root_dir + "Figure 5\\AntibodyTitration_treatments_FlowResults_Jurkats.csv", header=0)
	df_tension = pd.read_csv(root_dir + "Figure 5\\CellTension_Susp.csv")
	df_cal = pd.read_csv(root_dir + "Calibration_values.csv", header=0)

	## Using the calibration data from AF647 beads, compute the MESF for each AF647 intensity
	df_titration = MESFcalc(df_cal, df_antibody).drop(columns=['Conc_MESFavg'])

	## Compute the number of antibodies per square micron using the average size of the cell
	## I already divided the MESF by the number of fluorophores per antibody, I just need to 
	## divide by the average surface area
	avg_SA = []
	cell_types = ['HL60', 'Raji B', 'Jurkat']
	for i in cell_types:
		avg_SA = 4*np.pi*(df_tension[df_tension['Cell Type']==i]['Radius'].mean())
		df_titration['Surface Density'] = df_titration['Intensity_MESFavg'] / avg_SA


	## Compute the average values for each cell type at each concentration of opsonin
	df_titration['Phago_frac'] = df_titration['Phago'].values/df_titration['Count'].values.astype(int)*100
	df_titration['Trogo_frac'] = df_titration['Trogo'].values/df_titration['Count'].values.astype(int)*100

	# Split dataframe by phago vs trogo, then concatenate (to assist with seaborn plotting)

	plot_df = df_titration[['Cell Type', 'Surface Density', 'Trogo_frac','Treatment']]
	plot_df = plot_df.rename(columns = {'Trogo_frac': 'Trogo_Frac'})

	# ## Fit a sigmoid to each of the cell type trends
	cell_types = plot_df['Cell Type'].values
	treatments = plot_df['Treatment'].values

	#### Comment out this section if you've already saved a .csv of the fits (it speeds up plotting) #####
	# dfs = []

	# for i, j in zip(cell_types, treatments):

	# 	df = plot_df[(plot_df['Cell Type'] == i) & (plot_df['Treatment']== j)]
	# 	dataset = df[['Surface Density', 'Trogo_Frac']].values
	# 	Result = Fit(dataset, np.array([20, 80, -100]))
	# 	low, high, mid, x_array = BootStrap(Result, dataset)

	# 	## Grab the inflection point value for antibody density
	# 	inflec = Result.x[1]

	# 	## Grab the uncertainty in the fit as well
	# 	fit_param_uncertainty = np.diag(Result.hess_inv)[1] ## The second element corresponds to the antibody density param

	# 	## Create a dataframe for the model and bootstrap results
	# 	dfs.append(pd.DataFrame({'Cell Type': i, 'Mid': mid, 'Intensity': x_array, 'Low': low, 'High' : high,
	# 							 'Inflection': inflec, 'fit_param_uncertainty': fit_param_uncertainty, 'Treatment': j}))

	# model_df = pd.concat(dfs)
	# model_df.to_csv(save_dir + "sigmoid_fit_results_Drug_Jurkats.csv")

	##### End comment-able section #####
	
	model_df = pd.read_csv(root_dir + "\\Figure 4\\sigmoid_fit_results_Drug_Jurkats.csv")
	model_avg = model_df.groupby(['Cell Type','Treatment']).mean()

	##### Plot it! ######

	# Set up  plotting and statistical annotation

	plotting_parameters_sample = {'data': plot_df,
						   'x' : 'Surface Density',
						   'y' : 'Trogo_Frac',
						   'hue' : 'Treatment',
						   'hue_order': ['WT', '0.002% glut', '0.0025% glut'],
						   'palette': [light_blue, blue, dark_blue]}


	plotting_parameters_fit = {'data': model_df,
						   'x' : 'Intensity',
						   'y' : 'Mid',
						   'hue' : 'Treatment',
						   'hue_order': ['WT', '0.002% glut', '0.0025% glut'],
						   'palette': [light_blue, blue, dark_blue]}

	## Figure size needs to be in inches, conversion factor
	mm = 1/25.4

	figs, ax = plt.subplots(1,1, figsize=(55*mm,50*mm))

	ax = sns.scatterplot(**plotting_parameters_sample, linestyle="None", alpha=0.8)
	ax = sns.lineplot(**plotting_parameters_fit, linewidth=2)

	cell_types = ['WT', '0.002% glut', '0.0025% glut']
	for i in cell_types:
		inflection = model_df[model_df['Treatment']==i]['Inflection'].unique()
		ax.axvline(x = inflection, color='black', linestyle='--', linewidth=1)

	ax.set_xscale('log')
	ax.set_xlabel('# antibodies / '+ r'$\mu m^2$', fontsize=10)
	ax.set_ylabel("% trogocytosis", fontsize=10)
	plt.tight_layout()
	plt.legend([],[], frameon=False)

	# Save a .pdf of the figure
	# plt.savefig(save_dir + "fig5b_Jurkats.pdf")
	plt.show()