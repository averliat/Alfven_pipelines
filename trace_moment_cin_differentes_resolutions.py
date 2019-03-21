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



#-------------------------------------------------
#Chargement moments cinetiques de la simulation HR
#-------------------------------------------------
path_1_hr = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_hr/Moment_cinetique_tot_output_1_44_barycentre.fits'
path_2_hr = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_hr/Integrale_valeurs_absolues_moments_cinetiques_output_1_44_barycentre.fits'

path_simu_hr='/gpfs/data1/phennebe/B335_noturb_norot_hydro_hr/'
num_min_hr=1
num_max_hr=44

nmbr_output_hr = num_max_hr-num_min_hr+1

simulation_time_hr = np.zeros(nmbr_output_hr)

for l in range(nmbr_output_hr):
	ro=pymses.RamsesOutput(path_simu_hr,l+1)
	factor_time_hr_Myr = ro.info['unit_time'].express(cst.Myr)
	simulation_time_hr[l] = ro.info['time']*factor_time_hr_Myr



#--------------------------------------------------
#Chargement moments cinetiques de la simulation HR2
#--------------------------------------------------
path_1_hr2 = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_hr2/Moment_cinetique_tot_output_1_102_barycentre.fits'
path_2_hr2 = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_hr2/Integrale_valeurs_absolues_moments_cinetiques_output_1_102_barycentre.fits'

path_simu_hr2='/gpfs/data1/phennebe/B335_noturb_norot_hydro_hr2/'
num_min_hr2=1
num_max_hr2=102

nmbr_output_hr2 = num_max_hr2-num_min_hr2+1

simulation_time_hr2 = np.zeros(nmbr_output_hr2)

for l in range(nmbr_output_hr2):
	ro=pymses.RamsesOutput(path_simu_hr2,l+1)
	factor_time_hr2_Myr = ro.info['unit_time'].express(cst.Myr)
	simulation_time_hr2[l] = ro.info['time']*factor_time_hr2_Myr







#--------------------------
#Ouverture de la sauvegarde
#--------------------------
J_tot_tab_hr = fits.open(path_1_hr)[0].data
J_abs_tot_tab_hr = fits.open(path_2_hr)[0].data

J_tot_tab_hr2 = fits.open(path_1_hr2)[0].data
J_abs_tot_tab_hr2 = fits.open(path_2_hr2)[0].data



#------------------
#Debut de la figure
#------------------
plt.clf()
#'''
plt.semilogy(simulation_time_hr, J_tot_tab_hr, color='midnightblue', label=ur'Moment cinétique total HR',marker='.')
plt.semilogy(simulation_time_hr2, J_tot_tab_hr2, color='deepskyblue', label=ur'Moment cinétique total HR2',marker='.') 
'''
plt.semilogy(simulation_time_hr, J_abs_tot_tab_hr, color= 'saddlebrown', label=ur'$\int |$ moments cinétiques $| $ HR',marker='.')
plt.semilogy(simulation_time_hr2, J_abs_tot_tab_hr2, color= 'chocolate', label=ur'$\int |$ moments cinétiques $| $ HR2',marker='.')
'''

plt.xlabel(ur"Temps ($Myr$)")
plt.ylabel(ur"Moment cinétique total ($kg.m^2.s^-1$)")
plt.title(ur'"B335_noturb_norot" a differentes resolutions')

plt.legend(loc='best')
plt.tightlayout()

plt.show()

#plt.savefig('../Comparaison_resolutions_B335_noturb_norot_hydro/Log_Moment_cinetique_tot_HR_1_44_HR2_1_102_barycentre_fig.pdf', bbox_inches='tight')
#Ouvrir avec evince
