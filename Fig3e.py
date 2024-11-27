""" This script loads three channel images of macrophages with GUVs as arrays. It finds where
the three channels overlap and outputs an image highlighting colocalized pixels. The script
is subdivided as follows:

	1. 3 channel images are loaded and each channel is separated
	2. Channels are thresholded to obtain masks of 0's and 1's
	3. The masks of lipid/antibody are multiplied by the mask of the macrophages
	4. The product of this multiplication is summed to determine overlap
	5. The original image is plotted multiplied by some alpha value to de-emphasize
	6. The new colocalization mask is plotted over the original image to highlight colocalization
	7. This new image is saved as a .tif """

## Computing and dataprocessing libraries
import os
import pandas as pd
import numpy as np
import seaborn as sns

## Image processing libraries
import skimage.filters as sf
import skimage.io as io
import scipy.ndimage as snd
import skimage.morphology as smo
from skimage import color

## Plotting libraries
import matplotlib.pyplot as plt

## Set global matplotlib parameters
plt.rcParams["font.family"] = "Calibri"
plt.rcParams["xtick.labelsize"] = 10 
plt.rcParams["ytick.labelsize"] = 10
plt.rcParams["axes.linewidth"] = 1.0
plt.rcParams["xtick.major.width"] = 1.0
plt.rcParams["ytick.major.width"] = 1.0
plt.rcParams["errorbar.capsize"] = 1.0

def FileSifter(data_dir):
	"""This function searches the data directory for .tiff files."""

	fnames = []

	for root, dirs, files in os.walk(data_dir):

		## Loop through the files
		## and find .tiff extensions
		## Append these to fnames.

		for file in files:
			if file.endswith(".tif"):
				fnames.append(os.path.join(root, file))

	return fnames

def ImageMask(fnames):
	""" This function loops through the filenames, imports the image, splits each image into
	the three color channels and creates a binary mask for each channel. It then computes the
	area of each mask and stores it in an array."""

	areas = []
	lipid_m = []
	cell_m = []
	antibody_m = []
	combined_m = []
	time = []

	for name in fnames:

		## Load in an image
		im = io.imread(name)

		time.append(int(name[name.rfind("time")+4:name.rfind(".tif")]))

		## Separate the image into three separate channels
		cell = im[:,:,0]
		lipid = im[:,:,1]
		antibody = im[:,:,2]

		## Create masks
		cell_smooth = sf.gaussian(cell, sigma=1)
		cell_thresh = snd.binary_fill_holes(cell_smooth > sf.threshold_triangle(cell_smooth))
		cell_mask = smo.remove_small_objects(cell_thresh, min_size=10)

		lipid_smooth = sf.gaussian(lipid, sigma=1)
		lipid_thresh = lipid_smooth > sf.threshold_triangle(lipid_smooth)
		lipid_mask = smo.remove_small_objects(lipid_thresh, min_size=5)

		antibody_smooth = sf.gaussian(antibody, sigma=2)
		antibody_thresh = antibody_smooth > sf.threshold_triangle(antibody_smooth)
		antibody_mask = smo.remove_small_objects(antibody_thresh, min_size=10)

		## Compute the colocalization masks (include only lipid)
		combined = cell_mask * lipid_mask

		## Store all the masks
		lipid_m.append(lipid_mask)
		antibody_m.append(antibody_mask)
		cell_m.append(cell_mask)
		combined_m.append(combined)

		## compute the area of the co-localization mask
		areas.append(np.sum(combined) / 9.2) # this is the resolution

	return areas, lipid_m, cell_m, antibody_m, combined_m, time


if __name__ == '__main__':

	## Retrieve data
	print("\nGetting the data...")
	data_dir = "C:\\Users\\caitl\\Documents\\Fletcher lab\\Data\\Cell Phagocytosis\\Data for April 19 presentation\\Dynamic Pipette Nov 4\\"
	save_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\Figure 3\\"
	fnames = FileSifter(data_dir)

	mask_areas, lipid_mask, cell_mask, antibody_mask, combined, time = ImageMask(fnames)

	df = pd.DataFrame({'time':time, 'areas':mask_areas}).sort_values(by='time')

	## Convert the time axis to seconds (time loop is 34 s/frame)
	df['time'] = df['time']*34

	# ## Multiply the mask by each channel and make a new image
	# cell_co = cell * combined
	# lipid_co = lipid * combined
	# antibody_co = antibody * combined

	# ## Create a new image from the masks
	# mask_im = np.dstack((cell_co, lipid_co, antibody_co))

	# ## Create a compressed gray-scale version of the original image
	# im_gray = color.rgb2gray(test_im)

	# figs, axes = plt.subplots()
	# axes.imshow(im_gray, cmap='gray')
	# axes.imshow(mask_im, alpha=0.5)

	# plt.show()

	# Plot the co-localization mask area over the time

	## Figure size needs to be in inches, conversion factor
	mm = 1/25.4

	plotting_parameters = {'data': df,
						   'x' : 'time',
						   'y' : 'areas',
						   'linewidth' : 3,
						   'color': "grey"}


	figs, ax = plt.subplots(1,1, figsize=(52*mm,50*mm))

	ax = sns.lineplot(**plotting_parameters)
	ax.set_xscale('linear')
	ax.locator_params(axis='x', nbins=4)
	ax.locator_params(axis='y', nbins=4)
	ax.set_xlabel('time (s)', fontsize=10)
	ax.set_ylabel(r'overlap area $(\mu m^{2}$)', fontsize=10)
	plt.tight_layout()


	figs.savefig(save_dir +  'Fig3e.pdf')
	plt.show()




