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


#------------------------------------------------------------
#Chargement moments cinetiques de la simulation sans rotation
#------------------------------------------------------------
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



#-------------------------------------------------------------------
#Chargement moments cinetiques de la simulation incluant 1% rotation
#-------------------------------------------------------------------
path_1_rot = '/home/users/mnt/averliat/analyses/B335_noturb_rot_hydro/Moment_cinetique_tot_output_1_98_barycentre.fits'
path_2_rot = '/home/users/mnt/averliat/analyses/B335_noturb_rot_hydro/Integrale_valeurs_absolues_moments_cinetiques_output_1_98_barycentre.fits'

path_simu_rot='/gpfs/data1/averliat/B335_noturb_rot_hydro/'
num_min_rot=1
num_max_rot=98

nmbr_output_rot = num_max_rot-num_min_rot+1


if lecture_time == True:
	simulation_time_rot = np.zeros(nmbr_output_rot)
	for l in range(nmbr_output_rot):
		ro=pymses.RamsesOutput(path_simu_rot,l+1)
		factor_time_rot_Myr = ro.info['unit_time'].express(cst.Myr)
		simulation_time_rot[l] = ro.info['time']*factor_time_rot_Myr



#---------------------------------------------------------------------
#Chargement moments cinetiques de la simulation incluant 0.1% rotation
#---------------------------------------------------------------------
path_1_01 = '/home/users/mnt/averliat/analyses/B335_noturb_rot0.1_hydro/Moment_cinetique_tot_output_1_78_barycentre.fits'
path_2_01 = '/home/users/mnt/averliat/analyses/B335_noturb_rot0.1_hydro/Integrale_valeurs_absolues_moments_cinetiques_output_1_78_barycentre.fits'

path_simu_01='/gpfs/data1/averliat/B335_noturb_rot0.1_hydro/'
num_min_01=1
num_max_01=78

nmbr_output_01 = num_max_01-num_min_01+1


if lecture_time == True:
	simulation_time_01 = np.zeros(nmbr_output_01)
	for l in range(nmbr_output_01):
		ro=pymses.RamsesOutput(path_simu_01,l+1)
		factor_time_01_Myr = ro.info['unit_time'].express(cst.Myr)
		simulation_time_01[l] = ro.info['time']*factor_time_01_Myr



#----------------------------------------------------------------------
#Chargement moments cinetiques de la simulation incluant 0.01% rotation
#----------------------------------------------------------------------
path_1_001 = '/home/users/mnt/averliat/analyses/B335_noturb_rot0.01_hydro/Moment_cinetique_tot_output_1_58_barycentre.fits'
path_2_001 = '/home/users/mnt/averliat/analyses/B335_noturb_rot0.01_hydro/Integrale_valeurs_absolues_moments_cinetiques_output_1_58_barycentre.fits'

path_simu_001='/gpfs/data1/averliat/B335_noturb_rot0.01_hydro/'
num_min_001=1
num_max_001=58

nmbr_output_001 = num_max_001-num_min_001+1


if lecture_time == True:
	simulation_time_001 = np.zeros(nmbr_output_001)
	for l in range(nmbr_output_001):
		ro=pymses.RamsesOutput(path_simu_001,l+1)
		factor_time_001_Myr = ro.info['unit_time'].express(cst.Myr)
		simulation_time_001[l] = ro.info['time']*factor_time_001_Myr



#-------------------------------------------------------------------
#Chargement moments cinetiques de la simulation incluant 2% rotation
#-------------------------------------------------------------------
path_1_2 = '/home/users/mnt/averliat/analyses/B335_noturb_rot2_hydro/Moment_cinetique_tot_output_1_20_barycentre.fits'
path_2_2 = '/home/users/mnt/averliat/analyses/B335_noturb_rot2_hydro/Integrale_valeurs_absolues_moments_cinetiques_output_1_20_barycentre.fits'

path_simu_2='/gpfs/data1/averliat/B335_noturb_rot2_hydro/'
num_min_2=1
num_max_2=20

nmbr_output_2 = num_max_2-num_min_2+1


if lecture_time == True:
	simulation_time_2 = np.zeros(nmbr_output_2)
	for l in range(nmbr_output_2):
		ro=pymses.RamsesOutput(path_simu_2,l+1)
		factor_time_2_Myr = ro.info['unit_time'].express(cst.Myr)
		simulation_time_2[l] = ro.info['time']*factor_time_2_Myr



#-------------------------------------------------------------------
#Chargement moments cinetiques de la simulation incluant 3% rotation
#-------------------------------------------------------------------
path_1_3 = '/home/users/mnt/averliat/analyses/B335_noturb_rot3_hydro/Moment_cinetique_tot_output_1_20_barycentre.fits'
path_2_3 = '/home/users/mnt/averliat/analyses/B335_noturb_rot3_hydro/Integrale_valeurs_absolues_moments_cinetiques_output_1_20_barycentre.fits'

path_simu_3='/gpfs/data1/averliat/B335_noturb_rot3_hydro/'
num_min_3=1
num_max_3=20

nmbr_output_3 = num_max_3-num_min_3+1


if lecture_time == True:
	simulation_time_3 = np.zeros(nmbr_output_3)
	for l in range(nmbr_output_3):
		ro=pymses.RamsesOutput(path_simu_3,l+1)
		factor_time_3_Myr = ro.info['unit_time'].express(cst.Myr)
		simulation_time_3[l] = ro.info['time']*factor_time_3_Myr



#--------------------------
#Ouverture de la sauvegarde
#--------------------------
if lecture_time == True:
	J_tot_tab_rot = fits.open(path_1_rot)[0].data
	J_abs_tot_tab_rot = fits.open(path_2_rot)[0].data

	J_tot_tab_01 = fits.open(path_1_01)[0].data
	J_abs_tot_tab_01 = fits.open(path_2_01)[0].data

	J_tot_tab_001 = fits.open(path_1_001)[0].data
	J_abs_tot_tab_001 = fits.open(path_2_001)[0].data

	J_tot_tab_2 = fits.open(path_1_2)[0].data
	J_abs_tot_tab_2 = fits.open(path_2_2)[0].data

	J_tot_tab_3 = fits.open(path_1_3)[0].data
	J_abs_tot_tab_3 = fits.open(path_2_3)[0].data

	J_tot_tab_norot = fits.open(path_1_norot)[0].data
	J_abs_tot_tab_norot = fits.open(path_2_norot)[0].data



lecture_time = False  #Ne pas modifier

#------------------
#Debut de la figure
#------------------
plt.clf()
'''
plt.semilogy(simulation_time_norot, J_tot_tab_norot, color='lightslategrey', label=ur'Moment cinétique total NOROT',marker='.')
plt.semilogy(simulation_time_001, J_tot_tab_001, color='black', label=ur'Moment cinétique total 0.01%',marker='.')
plt.semilogy(simulation_time_01, J_tot_tab_01, color='midnightblue', label=ur'Moment cinétique total 0.1%',marker='.') 
plt.semilogy(simulation_time_rot, J_tot_tab_rot, color='deepskyblue', label=ur'Moment cinétique total 1%',marker='.') 
plt.semilogy(simulation_time_2, J_tot_tab_2, color='mediumorchid', label=ur'Moment cinétique total 2%',marker='.') 
plt.semilogy(simulation_time_3, J_tot_tab_3, color='purple', label=ur'Moment cinétique total 3%',marker='.') 
'''
plt.semilogy(simulation_time_norot, J_abs_tot_tab_norot, color= 'lightslategrey', label=ur'$\int |$ moments cinétiques $| $ NOROT',marker='.')
plt.semilogy(simulation_time_001, J_abs_tot_tab_001, color= 'saddlebrown', label=ur'$\int |$ moments cinétiques $| $ 0.01%',marker='.')
plt.semilogy(simulation_time_01, J_abs_tot_tab_01, color= 'chocolate', label=ur'$\int |$ moments cinétiques $| $ 0.1%',marker='.')
plt.semilogy(simulation_time_rot, J_abs_tot_tab_rot, color= 'darkorange', label=ur'$\int |$ moments cinétiques $| $ 1%',marker='.')
plt.semilogy(simulation_time_2, J_abs_tot_tab_2, color='yellowgreen', label=ur'$\int |$ moments cinétiques $| $ 2%',marker='.') 
plt.semilogy(simulation_time_3, J_abs_tot_tab_3, color='darkgreen', label=ur'$\int |$ moments cinétiques $| $ 3%',marker='.') 
#'''

plt.xlabel(ur"Temps ($Myr$)")
plt.ylabel(ur"Moment cinétique total ($kg.m^2.s^-1$)")
plt.title(ur'"B335_noturb_norot" a differentes resolutions')

plt.legend(loc='best')

plt.show()

#plt.savefig('../Comparaison_rotation_B335_noturb_rot_hydro/Log_Moment_cinetique_total_hr1_44_rot0.01_1_58_rot0.1_1_78_rot_1_98_rot2_1_20_rot3_1_20_fig.pdf', bbox_inches='tight')
#Ouvrir avec evince
