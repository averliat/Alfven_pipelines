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
nmbre_simu = 5
simu=[0 for i in range(nmbre_simu)]
output=[0 for i in range(nmbre_simu)]
tag=[0 for i in range(nmbre_simu)]


simu[0] = 'B335_noturb_norot_hydro_pert_asym_aleatoire_forte'
output[0] = 66
tag[0] = 'F'

simu[1] = 'B335_noturb_norot_hydro_pert_asym_aleatoire_forte_lr'
output[1] = 45
tag[1] = 'F_lr'

simu[2] = 'B335_noturb_norot_hydro_pert_asym_aleatoire_tres_forte'
output[2] = 62
tag[2] = 'TF'

simu[3] = 'B335_noturb_norot_hydro_pert_asym_aleatoire30'
output[3] = 60
tag[3] = '30'

simu[4] = 'B335_noturb_norot_hydro_pert_asym_aleatoire70'
output[4] = 60
tag[4] = '70'


trace_t_corrige = True



#------------------------
#Construction des chemins
#------------------------
path_r=[0 for i in range(nmbre_simu)]
path_t=[0 for i in range(nmbre_simu)]
path_tcor=[0 for i in range(nmbre_simu)]


for i in range(nmbre_simu):
	path_r[i] = '/gpfs/data1/averliat/analyses/'+simu[i]+'/rayon_disque_par_pdf'+str(output[i])+'.fits'
	path_t[i] = '/gpfs/data1/averliat/analyses/'+simu[i]+'/simulation_time_par_pdf'+str(output[i])+'.fits'
	path_tcor[i] = '/gpfs/data1/averliat/analyses/'+simu[i]+'/simulation_time_cor_par_pdf'+str(output[i])+'.fits'



#-------------------------
#Ouverture des sauvegardes
#-------------------------
r=[0 for i in range(nmbre_simu)]
t=[0 for i in range(nmbre_simu)]
tcor=[0 for i in range(nmbre_simu)]

for i in range(nmbre_simu):
	r[i]=fits.open(path_r[i])[0].data
	t[i]=fits.open(path_t[i])[0].data
	tcor[i]=fits.open(path_tcor[i])[0].data



#------------------
#Debut de la figure
#------------------
cmappts = plt.get_cmap('magma')
colorspts = [cmappts(i) for i in np.linspace(0,0.9,nmbre_simu)]

plt.figure()
if trace_t_corrige == False:
	for i in range(nmbre_simu):
		plt.plot(t[i],r[i], marker='.',color=colorspts[i], label=tag[i])
	plt.xlabel(ur'Temps ($Myr$)')
	plt.ylabel(ur'Rayon du disque ($AU$)')
	plt.legend(loc='best')

if trace_t_corrige == True:
	for i in range(nmbre_simu):
		plt.plot(tcor[i],r[i], marker='.',color=colorspts[i], label=tag[i])
	plt.xlabel(ur'Temps corrigé ($Myr$)')
	plt.ylabel(ur'Rayon du disque ($AU$)')
	plt.legend(loc='best')

plt.show()

#plt.savefig('../Comparaison_pert_B335_noturb_norot_hydro/Log_Moment_cinetique_tot_PERTLLF_1_48_MAN_1_60_ALEA0.1_1_75_ALEA0.99_1_61_barycentre_fig.pdf', bbox_inches='tight')
#plt.savefig('../Comparaison_pert_B335_noturb_norot_hydro/Log_Moment_cinetique_abs_HR_1_44_ALEA0.1_1_75_ALEA0.99_1_61_barycentre_fig.pdf', bbox_inches='tight')
#plt.savefig('../Comparaison_pert_B335_noturb_norot_hydro/Log_Moment_cinetique_tot_HR2_1_102_ALEA0.99_1_61_barycentre_fig.pdf', bbox_inches='tight')
#Ouvrir avec evince
