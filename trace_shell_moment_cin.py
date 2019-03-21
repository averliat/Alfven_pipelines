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
simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire_forte'
output = 43
R_min = 0
R_max = 'all'
dr = 10
bar=True
tag = 'F'



#------------------
#Differents chemins
#------------------
path_analyses='/gpfs/data1/averliat/analyses/'+simu+'/'
file_save = 'Shell_moment_cin_et_cumul_output'+str(output)+'_rmin='+str(R_min)+'AU_rmax='+str(R_max)+'_dr='+str(dr)+'_bar='+str(bar)+'.hdf5'



#------------------------------------
#Chargement des differentes quantites
#------------------------------------
h5f = h5py.File(path_analyses+file_save, 'r')

radius_shell=h5f['radius_shell'][:]
J_shell_tab=h5f['J_shell_tab'][:]
J_abs_shell_tab=h5f['J_abs_shell_tab'][:]
J_int_tab=h5f['J_int_tab'][:]
J_abs_int_tab=h5f['J_abs_int_tab'][:]
J_x_shell=h5f['J_x_shell'][:]
J_y_shell=h5f['J_y_shell'][:]
J_z_shell=h5f['J_z_shell'][:]
J_x_int_tab=h5f['J_x_int_tab'][:]
J_y_int_tab=h5f['J_y_int_tab'][:]
J_z_int_tab=h5f['J_z_int_tab'][:]
J_x_y_z_int_tab=h5f['J_x_y_z_int_tab'][:]

h5f.close()



#------------------
#Debut de la figure
#------------------
cmappts = plt.get_cmap('magma')
colorspts = [cmappts(i) for i in np.linspace(0.1,0.9,7)]

'''
#Moments cinetique total (norme) et absolu
plt.figure()
plt.plot(radius_shell,J_shell_tab, marker='.',color=colorspts[0], label=ur'Norme du moment cinétique')
plt.plot(radius_shell,J_abs_shell_tab, marker='.',color=colorspts[1], label=ur'Moment cinetique absolu')
plt.xlabel(ur'Rayon ($AU$)')
plt.ylabel(ur'Moment cinetique ($kg.m^2.s^{-1}$)')
plt.legend(loc='best')
'''
'''
#Comme au-dessus mais cumules
plt.figure()
plt.plot(radius_shell,J_int_tab, marker='.',color=colorspts[2], label=ur'Norme du moment cinetique cumulée')
plt.plot(radius_shell,J_abs_int_tab, marker='.',color=colorspts[3], label=ur'Moment cinetique absolu cumulé')
plt.xlabel(ur'Rayon ($AU$)')
plt.ylabel(ur'Norme du moment cinetique cumulé ($kg.m^2.s^{-1}$)')
plt.legend(loc='best')
'''

#Moment cinetique suivant x, y et z
plt.figure()
plt.plot(radius_shell,J_x_shell, marker='.',color=colorspts[4], label=ur'x')
plt.plot(radius_shell,J_y_shell, marker='.',color=colorspts[5], label=ur'y')
plt.plot(radius_shell,J_z_shell, marker='.',color=colorspts[6], label=ur'z')
plt.xlabel(ur'Rayon ($AU$)')
plt.ylabel(ur"Moment cinetique selon l'axe ($kg.m^2.s^{-1}$)")
plt.legend(loc='best')

#Moment cinetique suivant x, y et z - zoom suivant ordonnees
'''
plt.figure()
plt.plot(radius_shell,J_x_shell, marker='.',color=colorspts[4], label=ur'x', linestyle='None')
plt.ylim(-2.3e43,2.3e43)
plt.xlabel(ur'Rayon ($AU$)')
plt.ylabel(ur"Moment cinetique selon l'axe ($kg.m^2.s^{-1}$)")
plt.legend(loc='best')

plt.figure()
plt.plot(radius_shell,J_y_shell, marker='.',color=colorspts[5], label=ur'y', linestyle='None')
plt.ylim(-2.3e43,2.3e43)
plt.xlabel(ur'Rayon ($AU$)')
plt.ylabel(ur"Moment cinetique selon l'axe ($kg.m^2.s^{-1}$)")
plt.legend(loc='best')

plt.figure()
plt.plot(radius_shell,J_z_shell, marker='.',color=colorspts[6], label=ur'z', linestyle='None')
plt.ylim(-2.3e43,2.3e43)
plt.xlabel(ur'Rayon ($AU$)')
plt.ylabel(ur"Moment cinetique selon l'axe ($kg.m^2.s^{-1}$)")
plt.legend(loc='best')
'''

#Moment cinetique cumule suivant x, y et z
plt.figure()
plt.plot(radius_shell,J_x_int_tab, marker='.',color=colorspts[4], label=ur'x')#, linestyle='None')
plt.plot(radius_shell,J_y_int_tab, marker='.',color=colorspts[5], label=ur'y')#, linestyle='None')
plt.plot(radius_shell,J_z_int_tab, marker='.',color=colorspts[6], label=ur'z')#, linestyle='None')
plt.xlabel(ur'Rayon ($AU$)')
plt.ylabel(ur"Moment cinetique cumulé selon l'axe ($kg.m^2.s^{-1}$)")
plt.legend(loc='best')

#Moment cinetique integre, le bon !
plt.figure()
plt.plot(radius_shell,J_x_y_z_int_tab, marker='.',color=colorspts[4], label=ur'Norme du moment cinétique cumulé')
plt.xlabel(ur'Rayon ($AU$)')
plt.ylabel(ur"Moment cinetique ($kg.m^2.s^{-1}$)")
plt.legend(loc='best')

plt.show()

#plt.savefig('../Comparaison_pert_B335_noturb_norot_hydro/Log_Moment_cinetique_tot_PERTLLF_1_48_MAN_1_60_ALEA0.1_1_75_ALEA0.99_1_61_barycentre_fig.pdf', bbox_inches='tight')
#plt.savefig('../Comparaison_pert_B335_noturb_norot_hydro/Log_Moment_cinetique_abs_HR_1_44_ALEA0.1_1_75_ALEA0.99_1_61_barycentre_fig.pdf', bbox_inches='tight')
#plt.savefig('../Comparaison_pert_B335_noturb_norot_hydro/Log_Moment_cinetique_tot_HR2_1_102_ALEA0.99_1_61_barycentre_fig.pdf', bbox_inches='tight')
#Ouvrir avec evince
