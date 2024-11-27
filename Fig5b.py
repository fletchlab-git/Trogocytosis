""" This script takes .csv files output from CellProfiler with information on the size of trogocytic bites
in cells. We want to see the distribution of sizes of bites under different cell tensions"""

## Computing and data processing libraries
import os 
import pandas as pd
import numpy as np
from scipy.stats import linregress
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
	save_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\Figure 5\\"
	
	## Load in a dataframe and add a column with the cell type
	df_j = pd.read_csv(root_dir + "Figure 5\\Jurkats_19Sept2023Bites.csv", header=0)
	df_j.insert(1, "Cell Type", "Jurkat")

	## Extract the drug condition from the filename and create a new dataframe with the columns of interest
	condition = []
	for filename in df_j['FileName_ColorImage']:

		## Extract the condition
		condition.append(filename[filename.rfind("MAX_")+4:filename.rfind("jurkats")])

	df_j['condition'] = condition

	## Create a new dataframe with columns of interest
	df_j = df_j[['Cell Type', 'condition', 'AreaShape_Area', 'Parent_Cells']]


	df_h = pd.read_csv(root_dir + "Figure 5\\HL60s_26Sept2023Bites.csv", header=0)
	df_h.insert(1, "Cell Type", "HL60")

	condition = []
	for filename in df_h['FileName_ColorImage']:

		## Extract the condition
		condition.append(filename[filename.rfind("MAX_")+4:filename.rfind("_100x")])

	df_h['condition'] = condition

	## Create a new dataframe with columns of interest
	df_h = df_h[['Cell Type', 'condition', 'AreaShape_Area', 'Parent_Cells']]

	## Combine the dataframes
	df = pd.concat([df_j, df_h])

	## Create a new dataframe containing only the columns of interest. We only want data where the column
	## 'Parent_Cells' is greater than 0. We can also extract the drug condition from the file name.
	df_cells = df[df['Parent_Cells'].values > 0]

	## Convert the area values to microns from pixels
	df_cells['AreaShape_Area'] = np.sqrt((df_cells['AreaShape_Area'].values / 15.3846)/np.pi)

	##### Plot it! ######

	## Set up  plotting and statistical annotation

	## Figure size needs to be in inches, conversion factor
	mm = 1/25.4

	plotting_parameters = {'data': df_cells[df_cells['condition']=='WT'],
						   'x' : 'AreaShape_Area',
						   'hue' : 'Cell Type',
						   'palette': ['k', 'k']}

	figs, ax = plt.subplots(1,1, figsize=(95*mm,73*mm))

	ax = sns.histplot(**plotting_parameters)
	ax.set_xlabel("diameter of bites (um)", fontsize=12)
	ax.set_ylabel(" ", fontsize=12)
	plt.tight_layout()

	# Save a .png of the figure
	plt.savefig(save_dir + "FigS9.pdf")
	plt.show()





	