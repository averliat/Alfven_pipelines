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



#-----------------------------------------------------------
#Noms des simulations et caracteristiques du calcul du rayon
#-----------------------------------------------------------
simu_1 = 'B335_noturb_norot_hydro_pert_asym_aleatoire_forte'
output_max_1 = 66
dr_1 = 1

simu_2 = 'B335_noturb_norot_hydro_pert_asym_aleatoire_forte'
output_max_2 = 66
dr_2 = 1


t_cor = True
comparaison_moments_cumules = False
critere_r_disk = 0.2



#------------------------
#Construction des chemins
#------------------------
path_r_1 = '/gpfs/data1/averliat/analyses/'+simu_1+'/rayon_disque_'+str(output_max_1)+'_dr'+str(dr_1)+'_cumul='+str(comparaison_moments_cumules)+'_critere'+str(critere_r_disk)+'.fits'
path_t_1 = '/gpfs/data1/averliat/analyses/'+simu_1+'/simulation_time_'+str(output_max_1)+'_dr'+str(dr_1)+'_cumul='+str(comparaison_moments_cumules)+'_critere'+str(critere_r_disk)+'.fits'
path_tcor_1 = '/gpfs/data1/averliat/analyses/'+simu_1+'/simulation_time_cor_'+str(output_max_1)+'_dr'+str(dr_1)+'_cumul='+str(comparaison_moments_cumules)+'_critere'+str(critere_r_disk)+'.fits'

path_r_2 = '/gpfs/data1/averliat/analyses/'+simu_2+'/rayon_disque_'+str(output_max_2)+'_dr'+str(dr_2)+'_cumul='+str(comparaison_moments_cumules)+'_critere'+str(critere_r_disk)+'.fits'
path_t_2 = '/gpfs/data1/averliat/analyses/'+simu_2+'/simulation_time_'+str(output_max_2)+'_dr'+str(dr_2)+'_cumul='+str(comparaison_moments_cumules)+'_critere'+str(critere_r_disk)+'.fits'
path_tcor_2 = '/gpfs/data1/averliat/analyses/'+simu_2+'/simulation_time_cor_'+str(output_max_2)+'_dr'+str(dr_2)+'_cumul='+str(comparaison_moments_cumules)+'_critere'+str(critere_r_disk)+'.fits'



#nmbr_output_norot = num_max_norot-num_min_norot+1





#-------------------------
#Ouverture des sauvegardes
#-------------------------
r_1 = fits.open(path_r_1)[0].data
t_1 = fits.open(path_t_1)[0].data
tcor_1 = fits.open(path_tcor_1)[0].data

r_2 = fits.open(path_r_2)[0].data
t_2 = fits.open(path_t_2)[0].data
tcor_2 = fits.open(path_tcor_2)[0].data



#------------------
#Debut de la figure
#------------------
plt.figure()
if t_cor == False:
	plt.plot(t_1,r_1, marker='.',color='midnightblue', label=ur'$R_{disk}$')
	plt.plot(t_2,r_2, marker='.',color='chocolate', label=ur'$R_{disk}$')
	plt.xlabel(ur'Temps ($Myr$)')
	plt.ylabel(ur'Rayon du disque ($AU$)')
	plt.legend(loc='best')

if t_cor == True:
	plt.plot(tcor_1,r_1, marker='.',color='midnightblue', label=ur'$R_{disk}$')
	plt.plot(tcor_2,r_2, marker='.',color='chocolate', label=ur'$R_{disk}$')
	plt.xlabel(ur'Temps corrigé ($Myr$)')
	plt.ylabel(ur'Rayon du disque ($AU$)')
	plt.legend(loc='best')

plt.show()

#plt.savefig('../Comparaison_pert_B335_noturb_norot_hydro/Log_Moment_cinetique_tot_PERTLLF_1_48_MAN_1_60_ALEA0.1_1_75_ALEA0.99_1_61_barycentre_fig.pdf', bbox_inches='tight')
#plt.savefig('../Comparaison_pert_B335_noturb_norot_hydro/Log_Moment_cinetique_abs_HR_1_44_ALEA0.1_1_75_ALEA0.99_1_61_barycentre_fig.pdf', bbox_inches='tight')
#plt.savefig('../Comparaison_pert_B335_noturb_norot_hydro/Log_Moment_cinetique_tot_HR2_1_102_ALEA0.99_1_61_barycentre_fig.pdf', bbox_inches='tight')
#Ouvrir avec evince
