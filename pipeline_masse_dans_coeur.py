import numpy as np

import pymses
import sys
import glob as glob
import os 

from pymses.filters import CellsToPoints
from pymses.utils import constants as cst
from pymses.analysis import Camera, raytracing, slicing
from pymses.analysis import ScalarOperator, FractionOperator, MaxLevelOperator
import matplotlib.pyplot as plt
plt.ion()
from matplotlib.colors import Normalize

#from astropy.io import fits
import h5py

import os



if __name__ == '__main__':
	#-------------------------------------------------------
	#Entree le nom de la simulation et le numero de l'output
	#-------------------------------------------------------
	simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire10_hr'
	owner = 'averliat'
	output_min = 1
	output_max = 99

	save = True

	single = True
	num_output = 25

	rho_seuil = 1e11 #particules par cm^3

	name_save = 'Pourcent_mass_in_'+str(rho_seuil)+'_Hcc_output_'+str(output_min)+'_to_'+str(output_max)+'.hdf5'




def masse_dans_coeur_output_unique(simu, owner, num_output, rho_seuil):
	#--------------------------------------------------------------
	#Chemin de la simulation et du dossier de sauvegarde des images
	#--------------------------------------------------------------
	if owner=='averliat':
		path='/gpfs/data1/'+owner+'/'+simu+'/'
	if owner=='phennebe':
		path='/drf/projets/capucine/'+owner+'/'+simu+'/'
	path_save='/gpfs/data1/averliat/analyses/'+simu+'/'



	#-------------------
	#Lecture de l'output
	#-------------------
	ro=pymses.RamsesOutput(path,num_output)
	lbox_pc = ro.info['unit_length'].express(cst.pc)

	amr = ro.amr_source(["rho","vel","P","phi","g"])

	cell_source = CellsToPoints(amr)
	cells = cell_source.flatten(verbose=False)
	dx = cells.get_sizes()

	pos = cells.points
	rho = cells["rho"]



	#------------------------------------------------------------
	#Facteur de conversion des unites de code en unites physiques
	#------------------------------------------------------------
	lbox=ro.info['boxlen']
	lbox_m = ro.info['unit_length'].express(cst.m)
	factor_dist_m = ro.info['unit_length'].express(cst.m)
	factor_dist_cm = ro.info['unit_length'].express(cst.cm)
	lbox_au = ro.info['unit_length'].express(cst.au)
	lbox_cm = ro.info['unit_length'].express(cst.cm)
	factor_time_yr = ro.info['unit_time'].express(cst.year)

	factor_rho_Hcc = ro.info['unit_density'].express(cst.H_cc)
	factor_rho = ro.info['unit_density'].express(cst.kg/(cst.m)**3)



	#------------------------------
	#Conversion en unites physiques
	#------------------------------
	simulation_time = ro.info['time']*factor_time_yr
	dx *= factor_dist_m
	rho_Hcc = rho * factor_rho_Hcc
	rho *= factor_rho



	#---------------------------------------------------------
	#Definition du centre des images et de leur niveau de zoom
	#---------------------------------------------------------
	mask = rho_Hcc > rho_seuil

	rho_masked = rho[mask]
	dx_masked = dx[mask]

	mass_core = np.sum( rho_masked * dx_masked**3 ) / 1.9889e+30  #En Msun
	mass_tot = np.sum( rho * dx**3 ) / 1.9889e+30  #En Msun

	pourcent_mass_core = mass_core / mass_tot *100

	print
	print("==========================================================================")
	print("    M_central_core = "+str(mass_core)+ " Msun")
	print("    M_tot = "+str(mass_tot)+ " Msun")
	print("    Pourcent_M_central_core = "+str(pourcent_mass_core)+ " %")
	print("==========================================================================")
	print
	
	return (mass_core, mass_tot, pourcent_mass_core)



if __name__=='__main__':
	if single==True:
		masse_dans_coeur_output_unique(simu, owner, num_output, rho_seuil)



	else:
		#--------------------------------------------------------------
		#Chemin de la simulation et du dossier de sauvegarde des images
		#--------------------------------------------------------------
		if owner=='averliat':
			path='/gpfs/data1/'+owner+'/'+simu+'/'
		if owner=='phennebe':
			path='/drf/projets/capucine/'+owner+'/'+simu+'/'
		path_save='/gpfs/data1/averliat/analyses/'+simu+'/'


		nbre_output=output_max-output_min+1

		pourcent_mass_core = np.zeros(nbre_output)
		simulation_time = np.zeros(nbre_output)

		output_manquant=0
		for l in range(nbre_output):
			current_output = (l+1)
			try:
				ro=pymses.RamsesOutput(path, current_output)
			except:
				output_manquant += 1
				continue

			#-------------------
			#Lecture de l'output
			#-------------------
			lbox_pc = ro.info['unit_length'].express(cst.pc)

			amr = ro.amr_source(["rho","vel","P","phi","g"])

			cell_source = CellsToPoints(amr)
			cells = cell_source.flatten(verbose=False)
			dx = cells.get_sizes()

			pos = cells.points
			rho = cells["rho"]



			#------------------------------------------------------------
			#Facteur de conversion des unites de code en unites physiques
			#------------------------------------------------------------
			lbox=ro.info['boxlen']
			lbox_m = ro.info['unit_length'].express(cst.m)
			factor_dist_m = ro.info['unit_length'].express(cst.m)
			factor_dist_cm = ro.info['unit_length'].express(cst.cm)
			lbox_au = ro.info['unit_length'].express(cst.au)
			lbox_cm = ro.info['unit_length'].express(cst.cm)
			factor_time_yr = ro.info['unit_time'].express(cst.year)

			factor_rho_Hcc = ro.info['unit_density'].express(cst.H_cc)
			factor_rho = ro.info['unit_density'].express(cst.kg/(cst.m)**3)



			#------------------------------
			#Conversion en unites physiques
			#------------------------------
			simulation_time[l] = ro.info['time']*factor_time_yr
			dx *= factor_dist_m
			rho_Hcc = rho * factor_rho_Hcc
			rho *= factor_rho



			#---------------
			#Mask et calculs
			#---------------
			output_manquant=0

			mask = rho_Hcc > rho_seuil

			rho_masked = rho[mask]
			dx_masked = dx[mask]

			mass_core = np.sum( rho_masked * dx_masked**3 ) / 1.9889e+30  #En Msun
			mass_tot = np.sum( rho * dx**3 ) / 1.9889e+30  #En Msun
		
			pourcent_mass_core[l] = mass_core / mass_tot *100

			print
			print("==========================================================================")
			print("    M_central_core = "+str(mass_core)+ " Msun")
			print("    M_tot = "+str(mass_tot)+ " Msun")
			print("    Pourcent_M_central_core = "+str(pourcent_mass_core[l])+ " %")
			print("==========================================================================")
			print

			if save == True:
				h5f = h5py.File(path_save+name_save, 'w')

				h5f.create_dataset('pourcent_mass_core', data=pourcent_mass_core)
				h5f.create_dataset('simulation_time',data=simulation_time)

				h5f.close()






