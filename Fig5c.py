""" This script imports sigmoid fit data, along with the inflection point of those fits (defined as
the 'critical antibody concentration'), and cell tension. The plot shows a linear correlation between
the critical antibody concentration and cell tension. This is for figure 5."""

## Computing and data processing libraries
import os 
import pandas as pd
import numpy as np

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

# Set some colors as well, used in the plots
## throughout.
blue = "#43599c"

def Linear(x, theta):

	"""This function computes a linear regression"""

	m, b = theta

	return (m*x) + b

def L2_cost_function(dataset, theta):
	'''This function computes a weighted least squares loss (or cost) function
	for a fit to data. The loss function has the form L2 = sum((y_true - y_pred)**2)'''

	y_pred = Linear(dataset[:,0], theta) 
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
		samples[i] = Linear(x, p_i)


	## Compute the percentiles
	low = np.percentile(samples, 10, axis=0)
	high = np.percentile(samples, 90, axis=0)
	mid = Linear(x, np.array([ModelResults.x[0], ModelResults.x[1]]))

	return low, high, mid, x

if __name__ == '__main__':

	## Root directory to search for files
	root_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\"
	save_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\Figure 5\\"

	df_tension = pd.read_csv(root_dir + "Figure 5\\CellTension_Susp.csv")
	model_df = pd.read_csv(root_dir + "\\Figure 5\\sigmoid_fit_results_Drug_HL60s.csv")
	model_avg = model_df.groupby(['Cell Type','Treatment']).mean()

	## Group the tension dataframe by cell type and condition; compute the mean and standard dev
	df_avg = df_tension.groupby(['Cell Type','Treatment']).mean()
	df_std = df_tension.groupby(['Cell Type','Treatment']).std()
	df_std = df_std.rename(columns={'Tension':"tension_std", 'Radius':"radius_std"})

	# Make a dataframe of the average values (and standard deviations) of tension and the
	# inflection point from the sigmoid curves. This is to plot a correlation between these
	# two parameters
	corr_df = pd.concat([model_avg, df_avg, df_std], axis=1).reset_index().sort_values('Inflection')
	corr_df = corr_df[corr_df['Cell Type'] == 'HL60']

	# ## Linear regression of corr_df curve
	Result = Fit(corr_df[['Inflection', 'Tension']].values, np.array([400, 0]))
	low, high, mid, x_array = BootStrap(Result, corr_df[['Inflection', 'Tension']].values)

	##### Plot it! ######

	# Set up  plotting and statistical annotation

	## Figure size needs to be in inches, conversion factor
	mm = 1/25.4

	figs, ax = plt.subplots(1,1, figsize=(50*mm,40*mm))

	ax.plot(x_array, mid, color=blue, marker='None', linestyle='-', linewidth=2)
	ax.plot(corr_df['Inflection'], corr_df['Tension'], color=blue, marker='o', markersize=5, linestyle='None', alpha=0.7)
	ax.errorbar(corr_df['Inflection'], corr_df['Tension'],yerr=corr_df['tension_std'], xerr=corr_df['fit_param_uncertainty'],linestyle='None', color=blue, elinewidth=1)#, xerr=model_avg['fit_param_uncertainty']
	ax.fill_between(x_array, low, high, alpha=0.2, color=blue, linestyle='None')
	ax.set_xscale('linear')
	ax.set_xlabel(r"$\rho_{crit}$" +" " + "(# antibodies / " + "$\\mu m^{2}$)", fontsize=8)
	ax.set_ylabel("cell tension (mN/m)", fontsize=8)

	plt.tight_layout()
	plt.legend([],[], frameon=False)

	# Save a .pdf of the figure
	plt.savefig(save_dir + "fig5b_HL60s_fits.pdf")
	plt.show()
