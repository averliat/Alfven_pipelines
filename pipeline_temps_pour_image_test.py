# -*- coding: utf-8 -*-

#----------------------------------------------------------------------------------------
#Script lisant les outputs d'une simulation et sortant le numero de l'output 
#(et le temps associe) ou le seuil de densite maximal est depasse. On peut 
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

from astropy.io import fits



#--------------------------------------------
#Caracteristiques de la simulation consideree
#--------------------------------------------
path_simu='/gpfs/data1/averliat/B335_noturb_rot3_hydro/'
num_min=1
num_max=100  #Pas obligatoirement l'output max, peut depasser

seuil_rho = 1e-13  #kg.m^-3

nmbr_output = num_max-num_min+1



#-------------------------------------------------
#Lecture du temps et des densites de chaque output
#-------------------------------------------------
simulation_time = np.zeros(nmbr_output)
max_rho = np.zeros(nmbr_output)

output_seuil = 0
time_seuil = 0

for l in range(nmbr_output):
	if os.path.isdir(path_simu + 'output_' + str(l+1).zfill(5)) == False:
		print("===========================================================")
		print("===========================================================")
		print("      /!\  Output manquante ou seuil trop haut  /!\ ")
		print("===========================================================")
		print("===========================================================")
		break

	ro=pymses.RamsesOutput(path_simu,l+1)

	factor_time_Myr = ro.info['unit_time'].express(cst.Myr)
	simulation_time[l] = ro.info['time']*factor_time_Myr

	amr = ro.amr_source(["rho"])
	cell_source = CellsToPoints(amr)
	cells = cell_source.flatten()

	rho = cells["rho"]
	factor_rho = ro.info['unit_density'].express(cst.kg/(cst.m)**3)
	rho *= factor_rho

	max_rho[l] = np.max(rho)
	
	print
	print
	print("===========================================================")
	print('max(rho) =  ' +str(np.max(rho)) + '  kg.m^-3')
	print("===========================================================")
	print
	print
	print

	if max_rho[l] >= seuil_rho:
		output_seuil = l+1
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
	print("          /!\  Seuil de densité non atteint  /!\ ")
	print("===========================================================")
	print("===========================================================")


'''
plt.semilogy(simulation_time, max_rho, linestyle=None, marker='.')
plt.show()
'''



