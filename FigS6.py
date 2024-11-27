""" This plots a histogram of GUV tensions for Figure S6."""

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

if __name__ == '__main__':

	## Root directory to search for files
	root_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\"
	save_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\Figure 2\\"

	## Load in a dataframe and define the fractions of cells that have
	## phagocytosed or trogocytosed
	df = pd.read_csv(root_dir + "Figure 2\\GUVtensions_only.csv", header=0)
	print(df)

	##### Plot it! ######

	## Figure size needs to be in inches, conversion factor
	mm = 1/25.4

	plotting_parameters = {'data': df,
						   'x' : 'Tension',
						   'bins':27}

	figs, ax = plt.subplots(1,1, figsize=(75*mm,60*mm))

	ax = sns.histplot(**plotting_parameters, color='k', alpha=0.7)
	ax.set_xlabel("tension (mN/m)", fontsize=12)
	ax.set_ylabel("count", fontsize=12)
	plt.legend([],[], frameon=False)
	plt.tight_layout()

	# Save a .png of the figure
	plt.savefig(save_dir + "figS6.png")
	plt.show()