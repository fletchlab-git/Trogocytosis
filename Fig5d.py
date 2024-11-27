"""This script plots a phase diagram of trogocytosis and phagocytosis based on
the antibody density of a target and the tension of the target."""

## Computing and data processing libraries
import os 
import pandas as pd
import numpy as np
import scipy


## Plotting libraries
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors


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
	save_dir = "C:\\Users\\caitl\\Documents\\Fletcher Lab\\Data\\Trogocytosis Paper\\Figure 5\\"

	## Define variables measured from the literature
	# kT in Newton-meters
	kT = 4.1*1e-12*1e-9

	# Magnitude of normal stresses at the phagocytic interface in Pa (from Theriot ref)
	sigma_normal = 150

	# Effective target tension in mN/m
	gamma_t = 1e-1

	# Membrane bending stiffness in Nm
	kappa_t = 400*kT

	# Size of phagocytic interface in m
	interface_size = 1e-6

	# Target cell size in m
	cell_size = 10e-6

	# Define a critical antibody density in which we start to see trogocytosis in #/um^2
	rho_AB_crit = 100

	# Define the maximum normal stress as measured in Theriot ref (in Pa)
	sigma_norm_max = 150

	## Define some functions to compute trogocytic bite size from these parameters
	def tension_bite_size(gamma_t, sigma_normal):
		'''This function is the scaling relationship for tension'''
		return gamma_t/sigma_normal

	def bending_bite_size(kappa_t, sigma_normal):
		'''This function is the scaling relationship for bending rigidity (at vanshingly small gamma)'''
		return (kappa_t/(sigma_normal)**(1/3))

	def sigma_normal(rho_AB):
		'''This function includes the effects of antibody density into the normal stress'''

		if rho_AB <= rho_AB_crit:
			return 0
		else:
			return sigma_norm_max*(rho_AB/(rho_AB + rho_AB_crit))

	def phago_trogo_phase_diagram(rho_AB, gamma_t):
		'''This function returns the phase diagram of phagocytic behaviors:
		No interaction: 0
		Trogocytosis: 1
		Phagocytosis: 2 '''

		# Calculate the AB density dependent normal stress magnitude
		sigma_norm_mag = sigma_normal(rho_AB)

		# No interaction
		if sigma_norm_mag==0:
			return 0, np.nan

		# Trogocytosis and Phagocytosis (something happens)
		else:
			# Tension bite size (at high enough gamma to dominate)
			R_tension = tension_bite_size(gamma_t, sigma_norm_mag)

			# Bending bite size (at low gamma, bending dominates)
			R_bend = bending_bite_size(kappa_t, sigma_norm_mag)

			if R_tension>=R_bend:
				R_bite = R_tension
			else:
				R_bite = R_bend

			# Phagocytosis
			if R_bite>=cell_size/10:
				return 2, R_bite
			# Trogocytosis
			else:
				return 1, R_bite

	## Set up parameters for plotting the phase diagram
	gamma_t_array = 1e-3*np.logspace(-2,0.2,100)
	rho_AB_array = np.linspace(0, 4000, 100)
	gamma_t_matrix, rho_AB_matrix = np.meshgrid(gamma_t_array, rho_AB_array)

	# Set up an empty matrix to contain the outcomes and a bite size matrix for R values
	outcome_matrix = np.zeros_like(gamma_t_matrix)
	bite_size_matrix = np.zeros_like(gamma_t_matrix)

	## Compute the phase diagram!
	for ii in range(len(gamma_t_array)):
		for jj in range(len(rho_AB_array)):

			gamma_t, rho_AB = gamma_t_matrix[jj, ii], rho_AB_matrix[jj, ii]

			outcome_matrix[jj, ii], bite_size_matrix[jj,ii] = phago_trogo_phase_diagram(rho_AB, gamma_t)

	#### Plot it! ####

	
	## Figure size needs to be in inches, conversion factor
	mm = 1/25.4

	## Get the colormap
	cmap = plt.get_cmap('bone')

	figs, ax = plt.subplots(1,1, figsize=(95*mm,73*mm))

	CS = ax.contourf(gamma_t_array*1e3, rho_AB_array, outcome_matrix, cmap = cmap, alpha=0.3)

	ax.plot(0.93, 352, marker='o', markersize=5, color='k')
	ax.plot(0.09, 2102, marker='o', markersize=5, color='k')
	ax.plot(0.09, 14, marker='o', markersize=5, color='k')


	ax.set_xlabel('cortical tension (mN/m)', fontsize=12)
	ax.set_ylabel(r"density of antibody" +" " + "(#/ " + "$\\mu m^{2}$)", fontsize=12)
	ax.set_xscale('log')
	ax.set_ylim(0, 2500)

	plt.tight_layout()
	plt.savefig(save_dir + 'phase_diagram.pdf')
	plt.show()
















