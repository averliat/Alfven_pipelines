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



#----------------------------------------------------------------------
#Chargement moments cinetiques de la simulation et calcul de l'abscisse
#----------------------------------------------------------------------
path_1_hr2 = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_hr2/Moment_cinetique_tot_output_1_36.fits'
path_2_hr2 = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_hr2/Integrale_valeurs_absolues_moments_cinetiques_output_1_36.fits'

path_3_hr2 = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_hr2/Moment_cinetique_tot_output_1_55_barycentre.fits'
path_4_hr2 = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_hr2/Integrale_valeurs_absolues_moments_cinetiques_output_1_55_barycentre.fits'

path_simu_hr2='/gpfs/data1/phennebe/B335_noturb_norot_hydro_hr2/'
num_min_hr2=1
num_max_hr2=36


nmbr_output_hr2 = num_max_hr2-num_min_hr2+1

simulation_time_hr2 = np.zeros(nmbr_output_hr2)

for l in range(nmbr_output_hr2):
	ro=pymses.RamsesOutput(path_simu_hr2,l+1)
	factor_time_hr2_Myr = ro.info['unit_time'].express(cst.Myr)
	simulation_time_hr2[l] = ro.info['time']*factor_time_hr2_Myr



#--------------------------
#Ouverture de la sauvegarde
#--------------------------
J_tot_tab_hr2 = fits.open(path_1_hr2)[0].data[0:nmbr_output_hr2]
J_abs_tot_tab_hr2 = fits.open(path_2_hr2)[0].data[0:nmbr_output_hr2]

J_tot_tab_hr2_barycentre = fits.open(path_3_hr2)[0].data[0:nmbr_output_hr2]
J_abs_tot_tab_hr2_barycentre = fits.open(path_4_hr2)[0].data[0:nmbr_output_hr2]



#-----------------------------------------------------------
#Difference relative entre calculs avec centre et barycentre
#-----------------------------------------------------------
diff_abs = np.nan_to_num( (J_abs_tot_tab_hr2 - J_abs_tot_tab_hr2_barycentre)/J_abs_tot_tab_hr2_barycentre )
diff_tot = np.nan_to_num( (J_tot_tab_hr2 - J_tot_tab_hr2_barycentre)/J_tot_tab_hr2_barycentre )


#------------------
#Debut de la figure
#------------------
'''
plt.plot(simulation_time_hr2, J_tot_tab_hr2, color='midnightblue', label=ur'Moment cinétique total HR2 centre',marker='.')
plt.plot(simulation_time_hr2, J_tot_tab_hr2_barycentre, color='deepskyblue', label=ur'Moment cinétique total HR2 barycentre',marker='.')


plt.plot(simulation_time_hr2, J_abs_tot_tab_hr2, color= 'chocolate', label=ur'$\int |$ moments cinétiques $| $ HR2 centre',marker='.')
plt.plot(simulation_time_hr2, J_abs_tot_tab_hr2_barycentre, color= 'darkorange', label=ur'$\int |$ moments cinétiques $| $ HR2 barycentre',marker='.')
'''

#plt.ylim(np.min(diff_abs)-0.00002,np.max(diff_abs)+0.00002) #zoom
plt.ylim(-0.000009,0.000004) #zoom2
plt.xlim(0.0915,0.0945) #zoom2

plt.plot(simulation_time_hr2, diff_abs, color= 'chocolate', label=ur'$ \frac{\int | \sigma_{centre} | - \int | \sigma_{bar} | }{\int | \sigma_{bar} | } $',marker='.')

plt.plot(simulation_time_hr2, diff_tot, color='midnightblue', label=ur'$ \frac{\sigma_{centre} - \sigma_{bar} }{\sigma_{bar}} $',marker='.')


plt.xlabel(ur"Temps ($Myr$)")
plt.ylabel(ur"Moment cinétique total ($kg.m^2.s^-1$)")
plt.title(ur'"B335_noturb_norot_hr2" pour différents points de calcul')

plt.legend(loc='best')

plt.show()

#plt.savefig('../Comparaison_centres_B335_noturb_norot_hydro/Moment_cinetique_et_abs_HR2_1_55_centre_barycentre_fig_zoom2.pdf', bbox_inches='tight')
#Ouvrir avec evince
