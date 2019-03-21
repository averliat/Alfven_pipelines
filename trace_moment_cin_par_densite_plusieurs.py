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



class_rho=np.array([1e-21, 1e-18, 1e-15, 1e-12, 1e-9, 1e-6])

#------------------------------------------------------------------
#Chargement moments cinetiques de la simulation perturbee a la main
#------------------------------------------------------------------
path_1_pert_llf = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_pert_llf/Moment_cinetique_tot_output_1_48_par_densite.fits'
path_2_pert_llf = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_pert_llf/Integrale_valeurs_absolues_moments_cinetiques_output_1_48_par_densite.fits'

path_simu_pert_llf='/gpfs/data1/phennebe/B335_noturb_norot_hydro_pert_llf/'
num_min_pert_llf=1
num_max_pert_llf=48

nmbr_output_pert_llf = num_max_pert_llf-num_min_pert_llf+1

simulation_time_pert_llf = np.zeros(nmbr_output_pert_llf)

for l in range(nmbr_output_pert_llf):
	ro=pymses.RamsesOutput(path_simu_pert_llf,l+1)
	factor_time_pert_llf_Myr = ro.info['unit_time'].express(cst.Myr)
	simulation_time_pert_llf[l] = ro.info['time']*factor_time_pert_llf_Myr



#--------------------------------------------------
#Chargement moments cinetiques de la simulation HR2
#--------------------------------------------------
path_1_hr2 = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_hr2/Moment_cinetique_tot_output_1_78_par_densite.fits'
path_2_hr2 = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_hr2/Integrale_valeurs_absolues_moments_cinetiques_output_1_78_par_densite.fits'

path_simu_hr2='/gpfs/data1/phennebe/B335_noturb_norot_hydro_hr2/'
num_min_hr2=1
num_max_hr2=78

nmbr_output_hr2 = num_max_hr2-num_min_hr2+1

simulation_time_hr2 = np.zeros(nmbr_output_hr2)

for l in range(nmbr_output_hr2):
	ro=pymses.RamsesOutput(path_simu_hr2,l+1)
	factor_time_hr2_Myr = ro.info['unit_time'].express(cst.Myr)
	simulation_time_hr2[l] = ro.info['time']*factor_time_hr2_Myr



#--------------------------
#Ouverture de la sauvegarde
#--------------------------
J_tot_tab_pert_llf = fits.open(path_1_pert_llf)[0].data
J_abs_tot_tab_pert_llf = fits.open(path_2_pert_llf)[0].data

J_tot_tab_hr2 = fits.open(path_1_hr2)[0].data
J_abs_tot_tab_hr2 = fits.open(path_2_hr2)[0].data



#------------------
#Debut de la figure
#------------------
plt.clf()


#------------------
#Couleurs des plots
#------------------
#red = 255 /255.
#green = 0 /255.
#blue_scatter = np.linspace(0,255,len(class_rho)) /255.

cmappts = plt.get_cmap('viridis')
colorspts_J_tot_1 = [cmappts(i) for i in np.linspace(0.3,1,(len(class_rho)))]

cmappts = plt.get_cmap('magma')
colorspts_J_tot_2 = [cmappts(i) for i in np.linspace(0,0.8,(len(class_rho)))]


#-----------------
#Plot de J_tot_tab
#-----------------
#'''
for i in range(len(class_rho)):
	if i==0:
		label=ur'$\rho_{hr2}$ < ' + str(class_rho[i]) + ur' $kg.m^{-3}$'
	elif i != len(class_rho)-1:
		label=str(class_rho[i-1]) + ur'$ < \rho_{hr2}$ < ' + str(class_rho[i]) + ur' $kg.m^{-3}$'
	else:
		label=ur'$\rho_{hr2}$ > ' + str(class_rho[i]) + ur' $kg.m^{-3}$'

	plt.semilogy(simulation_time_hr2, J_tot_tab_hr2[:,i], marker='.', color=colorspts_J_tot_1[i], label=label)


for i in range(len(class_rho)):
	if i==0:
		label=ur'$\rho_{pert\_ llf}$ < ' + str(class_rho[i]) + ur' $kg.m^{-3}$'
	elif i != len(class_rho)-1:
		label=str(class_rho[i-1]) + ur'$ < \rho_{pert\_ llf}$ < ' + str(class_rho[i]) + ur' $kg.m^{-3}$'
	else:
		label=ur'$\rho_{pert\_ llf}$ > ' + str(class_rho[i]) + ur' $kg.m^{-3}$'

	plt.semilogy(simulation_time_pert_llf, J_tot_tab_pert_llf[:,i], marker='.', color=colorspts_J_tot_2[i], label=label)

#plt.ylabel(ur"Moment cinétique total ($kg.m^2.s^-1$)")
plt.ylabel(ur"Moment cinétique total ($kg.m^2.s^-1$)")


#---------------------
#Plot de J_abs_tot_tab
#---------------------
'''
for i in range(len(class_rho)):
	if i==0:
		label=ur'$\rho_{hr2}$ < ' + str(class_rho[i]) + ur' $kg.m^{-3}$'
	elif i != len(class_rho)-1:
		label=str(class_rho[i-1]) + ur'$ < \rho_{hr2}$ < ' + str(class_rho[i]) + ur' $kg.m^{-3}$'
	else:
		label=ur'$\rho_{hr2}$ > ' + str(class_rho[i]) + ur' $kg.m^{-3}$'

	plt.semilogy(simulation_time_hr2, J_abs_tot_tab_hr2[:,i], marker='.', color=colorspts_J_tot_1[i], label=label)


for i in range(len(class_rho)):
	if i==0:
		label=ur'$\rho_{pert\_ llf}$ < ' + str(class_rho[i]) + ur' $kg.m^{-3}$'
	elif i != len(class_rho)-1:
		label=str(class_rho[i-1]) + ur'$ < \rho_{pert\_ llf}$ < ' + str(class_rho[i]) + ur' $kg.m^{-3}$'
	else:
		label=ur'$\rho_{pert\_ llf}$ > ' + str(class_rho[i]) + ur' $kg.m^{-3}$'

	plt.semilogy(simulation_time_pert_llf, J_abs_tot_tab_pert_llf[:,i], marker='.', color=colorspts_J_tot_2[i], label=label)

#plt.ylabel(ur"Moment cinétique absolu ($kg.m^2.s^-1$)")
plt.ylabel(ur"Moment cinétique absolu par unité de masse ($m^2.s^-1$)")
'''

#-----------------------
#Legendes, titre et axes
#-----------------------
plt.xlabel(ur"Temps ($Myr$)")
plt.title(ur'"B335_noturb_norot_hydro HR2 et PERT_LLF"')

plt.legend(loc='best')

plt.show()

#plt.savefig('../Comparaison_classe_densité/Moment_cinetique_tot_PERT_LLF_1_48_HR2_1_78_par_densite_fig.pdf', bbox_inches='tight')
#Ouvrir avec evince
