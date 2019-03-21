# -*- coding: utf-8 -*-

#----------------------------------------------------------------------------------------
#Script lisant les outputs d'une simulation et sortant le numero de l'output 
#(et le temps associe) ou le seuil de masse acretee est depasse. On peut 
#alors considere ce temps comme le t0 de la simulation, permettant de comparer
#les simulations entre elles
#----------------------------------------------------------------------------------------


import numpy as np

import pymses
import sys
import glob as glob
import os 

from pymses.filters import CellsToPoints
from pymses.utils import constants as cst
import matplotlib.pyplot as plt
plt.ion()

#from astropy.io import fits



def temps_0_simu(path_simu, seuil_rho=1e11, critere_mass_seuil=10, num_min=1, num_max=150, sortie_output=1):
	#sortie_output=1 nouveau parametre pour prendre en compte si les output sont sortie tous les pas ou tous les n pas de temps (dans ce cas mettre sortie_output=n avec n=2, 3, ...)

	#--------------------------------------------
	#Caracteristiques de la simulation consideree
	#--------------------------------------------
	#simu='/gpfs/data1/averliat/B335_noturb_rot3_hydro/'
	#num_min=1
	#num_max=150  #Pas obligatoirement l'output max, peut depasser

	#seuil_rho = 1e-10  #kg.m^-3

	nmbr_output = num_max-num_min+1



	#-------------------------------------------------
	#Lecture du temps et des densites de chaque output
	#-------------------------------------------------
	simulation_time = np.zeros(nmbr_output)
	accreted_mass = np.zeros(nmbr_output)

	output_seuil = 0
	time_seuil = 0

	for l in range(nmbr_output):
		if os.path.isdir(path_simu + 'output_' + str(sortie_output*(l+1)).zfill(5)) == False:
			print("===========================================================")
			print("===========================================================")
			print("      /!\  Output manquante ou seuil trop haut  /!\ ")
			print("===========================================================")
			print("===========================================================")
			break

		ro=pymses.RamsesOutput(path_simu,sortie_output*(l+1))

		factor_time_Myr = ro.info['unit_time'].express(cst.Myr)
		simulation_time[l] = ro.info['time']*factor_time_Myr

		amr = ro.amr_source(["rho"])
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
		simulation_time_yr = ro.info['time']*factor_time_yr
		dx *= factor_dist_m
		rho_Hcc = rho * factor_rho_Hcc
		rho *= factor_rho



		#---------------------------------------------------------
		#Definition du centre des images et de leur niveau de zoom
		#---------------------------------------------------------
		mask = rho_Hcc > seuil_rho

		rho_masked = rho[mask]
		dx_masked = dx[mask]

		mass_core = np.sum( rho_masked * dx_masked**3 ) / 1.9889e+30  #En Msun
		mass_tot = np.sum( rho * dx**3 ) / 1.9889e+30  #En Msun

		pourcent_mass_core = mass_core / mass_tot *100

		print
		print
		print("==========================================================================")
		print("    M_central_core = "+str(mass_core)+ " Msun")
		print("    M_tot = "+str(mass_tot)+ " Msun")
		print("    Pourcent_M_central_core = "+str(pourcent_mass_core)+ " %")
		print("==========================================================================")
		print
		print
	


		accreted_mass[l] = pourcent_mass_core



		if accreted_mass[l] >= critere_mass_seuil:
			output_seuil = sortie_output*(l+1)
			time_seuil = simulation_time[l]
			print("===========================================================")
			print("===========================================================")
			print("         output_seuil = " +str(output_seuil))
			print("         time_seuil = " +str(time_seuil) +" Myr")
			print("===========================================================")
			print("===========================================================")
			break


	if output_seuil == 0:
		print("===========================================================")
		print("===========================================================")
		print("          /!\  Seuil de masse accretee non atteint  /!\ ")
		print("===========================================================")
		print("===========================================================")


	'''
	plt.semilogy(simulation_time, max_rho, linestyle=None, marker='.')
	plt.show()
	'''

	return (output_seuil, time_seuil)

