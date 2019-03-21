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



if 'lecture_time' not in locals():
	lecture_time = True  #Ne pas modifier, relancer le code avec "run -i trace_moment_cin_differentes_rotation" pour pas reouvrir les output


#---------------------------------------------------------------
#Chargement moments cinetiques de la simulation HR sans rotation
#---------------------------------------------------------------
path_1_norot = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_hr/Moment_cinetique_tot_output_1_44_barycentre.fits'
path_2_norot = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_hr/Integrale_valeurs_absolues_moments_cinetiques_output_1_44_barycentre.fits'

path_simu_norot='/gpfs/data1/phennebe/B335_noturb_norot_hydro_hr/'
num_min_norot=1
num_max_norot=44

nmbr_output_norot = num_max_norot-num_min_norot+1


if lecture_time == True:
	simulation_time_norot = np.zeros(nmbr_output_norot)
	for l in range(nmbr_output_norot):
		ro=pymses.RamsesOutput(path_simu_norot,l+1)
		factor_time_norot_Myr = ro.info['unit_time'].express(cst.Myr)
		simulation_time_norot[l] = ro.info['time']*factor_time_norot_Myr



#----------------------------------------------------------------
#Chargement moments cinetiques de la simulation HR2 sans rotation
#----------------------------------------------------------------
path_1_norot2 = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_hr2/Moment_cinetique_tot_output_1_102_barycentre.fits'
path_2_norot2 = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_hr2/Integrale_valeurs_absolues_moments_cinetiques_output_1_102_barycentre.fits'

path_simu_norot2='/gpfs/data1/phennebe/B335_noturb_norot_hydro_hr2/'
num_min_norot2=1
num_max_norot2=102

nmbr_output_norot2 = num_max_norot2-num_min_norot2+1


if lecture_time == True:
	simulation_time_norot2 = np.zeros(nmbr_output_norot2)
	for l in range(nmbr_output_norot2):
		ro=pymses.RamsesOutput(path_simu_norot2,l+1)
		factor_time_norot2_Myr = ro.info['unit_time'].express(cst.Myr)
		simulation_time_norot2[l] = ro.info['time']*factor_time_norot2_Myr



#-----------------------------------------------------------------------
#Chargement moments cinetiques de la simulation perturbee symetriquement
#-----------------------------------------------------------------------
path_1_pert_llf = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_pert_llf/Moment_cinetique_tot_output_1_48_barycentre.fits'
path_2_pert_llf = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_pert_llf/Integrale_valeurs_absolues_moments_cinetiques_output_1_48_barycentre.fits'

path_simu_pert_llf='/gpfs/data1/phennebe/B335_noturb_norot_hydro_pert_llf/'
num_min_pert_llf=1
num_max_pert_llf=48

nmbr_output_pert_llf = num_max_pert_llf-num_min_pert_llf+1


if lecture_time == True:
	simulation_time_pert_llf = np.zeros(nmbr_output_pert_llf)
	for l in range(nmbr_output_pert_llf):
		ro=pymses.RamsesOutput(path_simu_pert_llf,l+1)
		factor_time_pert_llf_Myr = ro.info['unit_time'].express(cst.Myr)
		simulation_time_pert_llf[l] = ro.info['time']*factor_time_pert_llf_Myr



#----------------------------------------------------------------------------------
#Chargement moments cinetiques de la simulation perturbee asymetriquement a la main
#----------------------------------------------------------------------------------
path_1_man = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_pert_asym_manuelle/Moment_cinetique_tot_output_1_60_barycentre.fits'
path_2_man = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_pert_asym_manuelle/Integrale_valeurs_absolues_moments_cinetiques_output_1_60_barycentre.fits'

path_simu_man='/gpfs/data1/averliat/B335_noturb_norot_hydro_pert_asym_manuelle/'
num_min_man=1
num_max_man=60

nmbr_output_man = num_max_man-num_min_man+1

if lecture_time == True:
	simulation_time_man = np.zeros(nmbr_output_man)
	for l in range(nmbr_output_man):
		ro=pymses.RamsesOutput(path_simu_man,l+1)
		factor_time_man_Myr = ro.info['unit_time'].express(cst.Myr)
		simulation_time_man[l] = ro.info['time']*factor_time_man_Myr



#--------------------------------------------------------------------------------------
#Chargement moments cinetiques de la simulation perturbee asymetriquement aleatoirement
#--------------------------------------------------------------------------------------
path_1_aleat = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_pert_asym_aleatoire/Moment_cinetique_tot_output_1_75_barycentre.fits'
path_2_aleat = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_pert_asym_aleatoire/Integrale_valeurs_absolues_moments_cinetiques_output_1_75_barycentre.fits'

path_simu_aleat='/gpfs/data1/averliat/B335_noturb_norot_hydro_pert_asym_aleatoire/'
num_min_aleat=1
num_max_aleat=75

nmbr_output_aleat = num_max_aleat-num_min_aleat+1

if lecture_time == True:
	simulation_time_aleat = np.zeros(nmbr_output_aleat)
	for l in range(nmbr_output_aleat):
		ro=pymses.RamsesOutput(path_simu_aleat,l+1)
		factor_time_aleat_Myr = ro.info['unit_time'].express(cst.Myr)
		simulation_time_aleat[l] = ro.info['time']*factor_time_aleat_Myr



#-------------------------------------------------------------------------------------------
#Chargement moments cinetiques de la simulation perturbee asymetriquement aleatoirement 0.99
#-------------------------------------------------------------------------------------------
path_1_aleat099 = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_pert_asym_aleatoire0.99/Moment_cinetique_tot_output_1_61_barycentre.fits'
path_2_aleat099 = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_pert_asym_aleatoire0.99/Integrale_valeurs_absolues_moments_cinetiques_output_1_61_barycentre.fits'

path_simu_aleat099='/gpfs/data1/averliat/B335_noturb_norot_hydro_pert_asym_aleatoire0.99/'
num_min_aleat099=1
num_max_aleat099=61

nmbr_output_aleat099 = num_max_aleat099-num_min_aleat099+1

if lecture_time == True:
	simulation_time_aleat099 = np.zeros(nmbr_output_aleat099)
	for l in range(nmbr_output_aleat099):
		ro=pymses.RamsesOutput(path_simu_aleat099,l+1)
		factor_time_aleat099_Myr = ro.info['unit_time'].express(cst.Myr)
		simulation_time_aleat099[l] = ro.info['time']*factor_time_aleat099_Myr



#--------------------------
#Ouverture de la sauvegarde
#--------------------------
if lecture_time == True:
	J_tot_tab_norot = fits.open(path_1_norot)[0].data
	J_abs_tot_tab_norot = fits.open(path_2_norot)[0].data

	J_tot_tab_norot2 = fits.open(path_1_norot2)[0].data
	J_abs_tot_tab_norot2 = fits.open(path_2_norot2)[0].data

	J_tot_tab_pert_llf = fits.open(path_1_pert_llf)[0].data
	J_abs_tot_tab_pert_llf = fits.open(path_2_pert_llf)[0].data

	J_tot_tab_man = fits.open(path_1_man)[0].data
	J_abs_tot_tab_man = fits.open(path_2_man)[0].data

	J_tot_tab_aleat = fits.open(path_1_aleat)[0].data
	J_abs_tot_tab_aleat = fits.open(path_2_aleat)[0].data

	J_tot_tab_aleat099 = fits.open(path_1_aleat099)[0].data
	J_abs_tot_tab_aleat099 = fits.open(path_2_aleat099)[0].data



lecture_time = False  #Ne pas modifier

#------------------
#Debut de la figure
#------------------
plt.clf()
#'''
#plt.semilogy(simulation_time_norot, J_tot_tab_norot, color='lightslategrey', label=ur'Moment cinétique total HR',marker='.')
plt.semilogy(simulation_time_norot2, J_tot_tab_norot2, color='lightslategrey', label=ur'Moment cinétique total HR2',marker='.')
#plt.semilogy(simulation_time_aleat, J_tot_tab_aleat, color='black', label=ur'Moment cinétique total ALEAT 0.1',marker='.')
plt.semilogy(simulation_time_aleat099, J_tot_tab_aleat099, color='midnightblue', label=ur'Moment cinétique total ALEAT 0.99',marker='.')
#plt.semilogy(simulation_time_man, J_tot_tab_man, color='deepskyblue', label=ur'Moment cinétique total MAN',marker='.') 
#plt.semilogy(simulation_time_pert_llf, J_tot_tab_pert_llf, color='mediumorchid', label=ur'Moment cinétique total PERT_LLF',marker='.') 


'''
#plt.semilogy(simulation_time_norot, J_abs_tot_tab_norot, color= 'lightslategrey', label=ur'$\int |$ moments cinétiques $| $ HR',marker='.')
plt.semilogy(simulation_time_norot2, J_abs_tot_tab_norot2, color= 'lightslategrey', label=ur'$\int |$ moments cinétiques $| $ HR2',marker='.')
#plt.semilogy(simulation_time_aleat, J_abs_tot_tab_aleat, color= 'saddlebrown', label=ur'$\int |$ moments cinétiques $| $ ALEAT 0.1',marker='.')
plt.semilogy(simulation_time_aleat099, J_abs_tot_tab_aleat099, color= 'chocolate', label=ur'$\int |$ moments cinétiques $| $ ALEAT 0.99',marker='.')
#plt.semilogy(simulation_time_man, J_abs_tot_tab_man, color= 'darkorange', label=ur'$\int |$ moments cinétiques $| $ MAN',marker='.')
#plt.semilogy(simulation_time_pert_llf, J_abs_tot_tab_pert_llf, color='darkgreen', label=ur'$\int |$ moments cinétiques $| $ PERT_LLF',marker='.')
'''

plt.xlabel(ur"Temps ($Myr$)")
plt.ylabel(ur"Moment cinétique total ($kg.m^2.s^-1$)")
plt.title(ur'"B335_noturb_norot_hydro_pert" pour différentes perturbations')

plt.legend(loc='best')

plt.show()

#plt.savefig('../Comparaison_pert_B335_noturb_norot_hydro/Log_Moment_cinetique_tot_PERTLLF_1_48_MAN_1_60_ALEA0.1_1_75_ALEA0.99_1_61_barycentre_fig.pdf', bbox_inches='tight')
#plt.savefig('../Comparaison_pert_B335_noturb_norot_hydro/Log_Moment_cinetique_abs_HR_1_44_ALEA0.1_1_75_ALEA0.99_1_61_barycentre_fig.pdf', bbox_inches='tight')
#plt.savefig('../Comparaison_pert_B335_noturb_norot_hydro/Log_Moment_cinetique_tot_HR2_1_102_ALEA0.99_1_61_barycentre_fig.pdf', bbox_inches='tight')
#Ouvrir avec evince
