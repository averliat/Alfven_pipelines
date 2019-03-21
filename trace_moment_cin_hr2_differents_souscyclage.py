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




#--------------------------------------------------
#Chargement moments cinetiques de la simulation hr2
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



#-------------------------------------------------------------------
#Chargement moments cinetiques de la simulation hr2 moins souscyclee
#-------------------------------------------------------------------
path_1_hr2_m_s = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_hr2_moins_souscyclee/Moment_cinetique_tot_output_1_103_barycentre.fits'
path_2_hr2_m_s = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_hr2_moins_souscyclee/Integrale_valeurs_absolues_moments_cinetiques_output_1_103_barycentre.fits'

path_simu_hr2_m_s='/gpfs/data1/averliat/B335_noturb_norot_hydro_hr2_moins_souscyclee'
num_min_hr2_m_s=1
num_max_hr2_m_s=103

nmbr_output_hr2_m_s = num_max_hr2_m_s-num_min_hr2_m_s+1

simulation_time_hr2_m_s = np.zeros(nmbr_output_hr2_m_s)

for l in range(nmbr_output_hr2_m_s):
	ro=pymses.RamsesOutput(path_simu_hr2_m_s,l+1)
	factor_time_hr2_m_s_Myr = ro.info['unit_time'].express(cst.Myr)
	simulation_time_hr2_m_s[l] = ro.info['time']*factor_time_hr2_m_s_Myr



#--------------------------
#Ouverture de la sauvegarde
#--------------------------
J_tot_tab_hr2 = fits.open(path_1_hr2)[0].data
J_abs_tot_tab_hr2 = fits.open(path_2_hr2)[0].data

J_tot_tab_hr2_m_s = fits.open(path_1_hr2_m_s)[0].data
J_abs_tot_tab_hr2_m_s = fits.open(path_2_hr2_m_s)[0].data



#------------------
#Debut de la figure
#------------------
plt.clf()
'''
plt.semilogy(simulation_time_hr2_m_s, J_tot_tab_hr2_m_s, color='deepskyblue', label=ur'Moment cinétique total HR2 moins sous-cyclée',marker='.')
plt.semilogy(simulation_time_hr2, J_tot_tab_hr2, color='mediumorchid', label=ur'Moment cinétique total HR2',marker='.') 


'''
plt.semilogy(simulation_time_hr2_m_s, J_abs_tot_tab_hr2_m_s, color= 'saddlebrown', label=ur'$\int |$ moments cinétiques $| $ HR2 moins sous-cyclée',marker='.')
plt.semilogy(simulation_time_hr2, J_abs_tot_tab_hr2, color='chocolate', label=ur'$\int |$ moments cinétiques $| $ HR2',marker='.')
#'''

plt.xlabel(ur"Temps ($Myr$)")
plt.ylabel(ur"Moment cinétique total ($kg.m^2.s^-1$)")
plt.title(ur'"B335_noturb_norot_hydro_pert" pour différents sous-cyclage')

plt.legend(loc='best')

plt.show()

#plt.savefig('../Comparaison_sous_cyclage_B335_noturb_norot_hydro_hr2/Log_Moment_cinetique_tot_HR2_1_102_HR2moins_souscyclee_1_103_barycentre_zoom_fig.pdf', bbox_inches='tight')
#Ouvrir avec evince
