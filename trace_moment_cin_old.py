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
path_1 = '/home/users/mnt/averliat/analyses/B335_noturb_rot_hydro/Moment_cinetique_tot_output_1_68_barycentre.fits'
path_2 = '/home/users/mnt/averliat/analyses/B335_noturb_rot_hydro/Integrale_valeurs_absolues_moments_cinetiques_output_1_68_barycentre.fits'

path_simu='/gpfs/data1/averliat/B335_noturb_rot_hydro/'
num_min=1
num_max=68



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
#output_tab = np.linspace(num_min,num_max, num_max-num_min+1)

#plt.plot(simulation_time, J_tot_tab, marker='.', color='midnightblue', label=ur'Moment cinétique total')
#plt.plot(simulation_time, J_abs_tot_tab, marker = '.', color= 'chocolate', label=ur'$\int |$ moments cinétiques $| $')

plt.plot(simulation_time, J_tot_tab, marker='.', color='midnightblue', label=ur'Moment cinétique total')
plt.plot(simulation_time, J_abs_tot_tab, marker = '.', color= 'chocolate', label=ur'$\int |$ moments cinétiques $| $')

plt.xlabel(ur"Temps ($Myr$)")
#plt.xlabel(ur"Numéro de l'output")
plt.ylabel(ur"Moment cinétique total ($kg.m^2.s^-1$)")
plt.title(ur'"B335_noturb_norot_hydro_pert_llf"')

plt.legend(loc='best')

plt.show()

#plt.savefig('../B335_noturb_norot_hydro_hr/Moment_cinetique_et_abs_output_1_44_fig.pdf', bbox_inches='tight')
#Ouvrir avec evince
