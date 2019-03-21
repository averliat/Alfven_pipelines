# -*- coding: utf-8 -*-
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



#-------------------------------------------------------------------
#Chemin de la sauvegarde du tableau contenant les moments cinetiques
#-------------------------------------------------------------------
path_1 = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_hr2/Moment_cinetique_tot_output_1_102_par_densite_normalise.fits'
path_2 = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_hr2/Integrale_valeurs_absolues_moments_cinetiques_output_1_102_par_densite_normalise.fits'

path_simu='/gpfs/data1/phennebe/B335_noturb_norot_hydro_hr2/'
num_min=1
num_max=102

class_rho=np.array([1e-21, 1e-18, 1e-15, 1e-12, 1e-9, 1e-6])



#---------------------------
#Calcul de l'abscisse en Myr
#---------------------------
nmbr_output = num_max - num_min + 1
simulation_time = np.zeros(nmbr_output)

for l in range(nmbr_output):
	ro=pymses.RamsesOutput(path_simu,l+1)
	factor_time_Myr = ro.info['unit_time'].express(cst.Myr)
	simulation_time[l] = ro.info['time']*factor_time_Myr



#--------------------------
#Ouverture de la sauvegarde
#--------------------------
J_tot_tab = fits.open(path_1)[0].data
J_abs_tot_tab = fits.open(path_2)[0].data



#------------------
#Debut de la figure
#------------------
plt.clf()
#output_tab = np.linspace(num_min,num_max, num_max-num_min+1

#------------------
#Couleurs des plots
#------------------
#red = 255 /255.
#green = 0 /255.
#blue_scatter = np.linspace(0,255,len(class_rho)) /255.

cmappts = plt.get_cmap('viridis')
colorspts_J_tot = [cmappts(i) for i in np.linspace(0,1,(len(class_rho)))]
cmappts = plt.get_cmap('magma')
colorspts_J_abs_tot = [cmappts(i) for i in np.linspace(0,0.8,(len(class_rho)))]


#-----------------
#Plot de J_tot_tab
#-----------------
#'''
for i in range(len(class_rho)):
	if i==0:
		label=ur'$\rho$ < ' + str(class_rho[i]) + ur' $kg.m^{-3}$'
	elif i != len(class_rho)-1:
		label=str(class_rho[i-1]) + ur'$ < \rho$ < ' + str(class_rho[i]) + ur' $kg.m^{-3}$'
	else:
		label=ur'$\rho$ > ' + str(class_rho[i]) + ur' $kg.m^{-3}$'

	plt.semilogy(simulation_time, J_tot_tab[:,i], marker='.', color=colorspts_J_tot[i], label=label)
#plt.ylabel(ur"Moment cinétique total ($kg.m^2.s^-1$)")
plt.ylabel(ur"Moment cinétique total par unité de masse ($m^2.s^-1$)")


#---------------------
#Plot de J_abs_tot_tab
#---------------------
'''
for i in range(len(class_rho)):
	if i==0:
		label=ur'$\rho$ < ' + str(class_rho[i]) + ur'$kg.m^{-3}$'
	elif i != len(class_rho)-1:
		label=str(class_rho[i-1]) + ur'$ < \rho$ < ' + str(class_rho[i]) + ur'$kg.m^{-3}$'
	else:
		label=ur'$\rho$ > ' + str(class_rho[i]) + ur'$kg.m^{-3}$'

	plt.semilogy(simulation_time, J_abs_tot_tab[:,i], marker = '.', color= colorspts_J_abs_tot[i], label=label)
#plt.ylabel(ur"Moment cinétique absolu ($kg.m^2.s^-1$)")
plt.ylabel(ur"Moment cinétique absolu par unité de masse ($m^2.s^-1$)")
'''

#-----------------------
#Legendes, titre et axes
#-----------------------
plt.xlabel(ur"Temps ($Myr$)")
plt.title(ur'"B335_noturb_norot_hydro_hr2"')

plt.legend(loc='best')

plt.show()

#plt.savefig('../B335_noturb_norot_hydro_hr2/Moment_cinetique_tot_output_1_102_par_densite_normalise_fig.pdf', bbox_inches='tight')
#Ouvrir avec evince
