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
import scipy.optimize as opt



if __name__=='__main__':
#-----------------------------------------------------------
#Noms des simulations et caracteristiques du calcul du rayon
#-----------------------------------------------------------
    tag='50_shr'
    legend='7,  10*40'  #Inutile pour le moment
    marker='.'

    simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire'+tag

    output = 685

    R_min = 50
    R_max = 17000
    dr = 500


    save_fig = True




def trace_V_gradient_vitesse_kolmogorov(simu,tag,output,R_min,R_max,dr,legend,marker,save_fig):
#------------------
#Differents chemins
#------------------
    path_analyses='/home/averliat/these/analyses/'+simu+'/'
    file_save_vinfall = 'J_turb_log_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.hdf5'



#------------------------------------
#Chargement des differentes quantites
#------------------------------------
    h5f = h5py.File(path_analyses+file_save_vinfall, 'r')

    shells_au=h5f['shells_au'][:]
    j_turb_shell_tab=h5f['j_turb_shell_tab'][:]
    j_turb_shell_tab_2=h5f['j_turb_shell_tab_2'][:]

    h5f.close()


    veltest=6.09e-7*shells_au**1.6



#------------------
#Debut de la figure
#------------------
    cmappts = plt.get_cmap('magma')
    colorspts = [cmappts(i) for i in np.linspace(0.1,0.9,4)]


    plt.figure()
    plt.loglog(shells_au,j_turb_shell_tab, marker='.',color=colorspts[1],label=ur'$<\Delta V_{los}>.R_{shell}$')
    plt.loglog(shells_au,j_turb_shell_tab_2, marker='.',color=colorspts[2],label=ur'$<\Delta V_{los}.R>$')

    plt.loglog(shells_au,veltest,color=colorspts[3],label=ur'$v\propto R^{1.6}$')
    plt.ylim((0.00051009,0.004))

    plt.xlabel(ur'Rayon ($AU$)')
    plt.ylabel(ur"Moment cinétique spécifique  $(km.s^{-1}.pc)$")
    plt.legend(loc='best')

    if save_fig==True:
        plt.savefig('/home/averliat/these/analyses/Gradients_vitesse/Analyses_champ_vitesse/j_turb_kolmogorov/j_turb_kolmogorov_'+simu+'_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.pdf', bbox_inches='tight')




if __name__=='__main__':
    trace_V_gradient_vitesse_kolmogorov(simu,tag,output,R_min,R_max,dr,legend,marker,save_fig)

