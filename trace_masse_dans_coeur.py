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
simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire50_hr'
owner = 'averliat'
output_min = 1
output_max = 99

tag='HR'
couleur=1

rho_seuil = 1e11 #particules par cm^3

name_save = 'Pourcent_mass_in_'+str(rho_seuil)+'_Hcc_output_'+str(output_min)+'_to_'+str(output_max)+'.hdf5'



#------------------
#Differents chemins
#------------------
path_analyses='/gpfs/data1/averliat/analyses/'+simu+'/'



#------------------------------------
#Chargement des differentes quantites
#------------------------------------
h5f = h5py.File(path_analyses+name_save, 'r')

pourcent_mass_core=h5f['pourcent_mass_core'][:]
simulation_time=h5f['simulation_time'][:]

h5f.close()



#------------------
#Debut de la figure
#------------------
cmappts = plt.get_cmap('magma')
colorspts = [cmappts(i) for i in np.linspace(0.1,0.9,7)]



#Moment cinetique integre, le bon !
#plt.figure()
plt.plot(simulation_time, pourcent_mass_core, marker='.',color=colorspts[couleur], label=ur'Masse du coeur '+tag, linestyle='None')
plt.xlabel(ur'Temps ($year$)')
#plt.ylabel(ur"Masse du coeur ($M_\odot$)")
plt.ylabel(ur"Masse du coeur / Masse de l'enveloppe (%)")
plt.legend(loc='best')

plt.show()

#plt.savefig('/gpfs/data1/averliat/analyses/Comparaison_masse_core/Masse_dans_coeur_sur_enveloppe_over_time_50hr_1to99_50vhr_1to99.pdf', bbox_inches='tight')

#Ouvrir avec evince
