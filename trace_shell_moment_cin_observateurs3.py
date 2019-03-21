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

import h5py



#-----------------------------------------------------------
#Noms des simulations et caracteristiques du calcul du rayon
#-----------------------------------------------------------
simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire50_vhr'
tag='hr2'

couleur=6
output_min = 1
output_max = 145


R_min = 0
R_max = 'all'
dr = 10
bar=True
#tag = 'F'



#------------------
#Differents chemins
#------------------
path_analyses='/mnt/magmist/magmist/simu_B335_averliat/analyses/'+simu+'/'
file_save = 'Moment_cin_observateurs_output'+str(output_min)+'_to_'+str(output_max)+'.hdf5'



#------------------------------------
#Chargement des differentes quantites
#------------------------------------
h5f = h5py.File(path_analyses+file_save, 'r')

norm_mom_cin_obs=h5f['norm_mom_cin_obs'][:]
simulation_time=h5f['simulation_time'][:]
norm_GC_AU=h5f['norm_GC_AU'][:]
norm_center_C_AU=h5f['norm_center_C_AU'][:]
norm_center_G_AU=h5f['norm_center_G_AU'][:]
vel_GC_tab=h5f['vel_GC_tab'][:]
norm_mom_cin_obs_par_integ=h5f['norm_mom_cin_obs_par_integ'][:]
erreur_mom=h5f['erreur_mom'][:]

h5f.close()



#------------------
#Debut de la figure
#------------------
cmappts = plt.get_cmap('magma')
colorspts = [cmappts(i) for i in np.linspace(0.1,0.9,7)]

tab_x = np.linspace(2,10,5) #Utile pour hr2
#tab_x = np.linspace(output_min,output_max,output_max-output_min+1)

#Moment cinetique integre, le bon !

#plt.figure()
'''
plt.semilogy(tab_x, norm_mom_cin_obs_par_integ, marker='.',color=colorspts[3], label=ur'Norme du moment cinetique'+tag)
plt.semilogy(tab_x, erreur_mom, marker='.',color=colorspts[5], label=ur'Erreur sur le moment cinetique'+tag)
plt.xlabel(ur'Output')
plt.ylabel(ur"Moment cinetique ($kg.m^2.s^{-1}$)")
plt.legend(loc='best')
'''

#Histogramme moins moche
plt.figure()
mask=(norm_mom_cin_obs_par_integ-norm_mom_cin_obs)>0
mask_inv=np.logical_not(mask)
plt.bar(tab_x[mask], norm_mom_cin_obs_par_integ[mask], width=1., color=colorspts[5], label=ur'Norme du moment cinetique reel '+tag)
plt.bar(tab_x[mask_inv], norm_mom_cin_obs[mask_inv], width=1., color=colorspts[1], label=ur'Norme du moment cinetique analytique '+tag)
plt.bar(tab_x[mask_inv], norm_mom_cin_obs_par_integ[mask_inv], width=1., color='w', edgecolor='black', linewidth=0.4)
plt.bar(tab_x[mask], norm_mom_cin_obs[mask], width=1., color='w', edgecolor='black', linewidth=0.4)
plt.xlabel(ur'Output')
plt.ylabel(ur"Moment cinetique ($kg.m^2.s^{-1}$)")
plt.legend(loc='best')


#plt.yscale('log')
#plt.xlim((0,53))
#plt.ylim((0,1.2e47))


plt.show()

#plt.savefig('../Comparaisons/Verification_moment_physique/Moment_cinetique_analytique_et_reel_hr2.pdf', bbox_inches='tight')

#plt.savefig('../Comparaisons/Verification_moment_physique/Difference_absolue_moment_cinetique_analytique_et_reel_hr2.pdf', bbox_inches='tight')

#plt.savefig('../Comparaisons/Verification_moment_physique/Difference_relative_moment_cinetique_analytique_et_reel_hr2.pdf', bbox_inches='tight')

#plt.savefig('../Comparaisons/Verification_moment_physique/Moment_cinetique_analytique_et_reel_hr2_output'+str(output_min)+'_'+str(output_max)+'.pdf', bbox_inches='tight')
