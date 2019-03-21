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

    output = 277

    R_min = 50
    R_max = 10000
    dr = 50


    save_fig = True




def trace_Vinfall_gradient_vitesse(simu,tag,output,R_min,R_max,dr,legend,marker,save_fig):
#------------------
#Differents chemins
#------------------
    path_analyses='/home/averliat/these/analyses/'+simu+'/'
    file_save_vinfall = 'Vinfall_gradient_vitesse_log_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.hdf5'
    #file_save_ang = 'Angle_omega_J_gradient_vitesse_log_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.hdf5'
    #file_save_ang = 'Angle_gradient_vitesse_output'+str(output)+'comparaison_mathilde.hdf5'



#------------------------------------
#Chargement des differentes quantites
#------------------------------------
    h5f = h5py.File(path_analyses+file_save_vinfall, 'r')

    shells_au=h5f['shells_au'][:]
    V_r_shell_moy_tab=h5f['V_r_shell_moy_tab'][:]
    V_r_rms_shell_moy_tab=h5f['V_r_rms_shell_moy_tab'][:]
    V_paral_shell_moy_tab=h5f['V_paral_shell_moy_tab'][:]
    sigma_shell_moy_tab=h5f['sigma_shell_moy_tab'][:]

    h5f.close()



#------------------
#Debut de la figure
#------------------
    cmappts = plt.get_cmap('magma')
    colorspts = [cmappts(i) for i in np.linspace(0.1,0.9,4)]


    plt.figure()
    plt.semilogx(shells_au,np.abs(V_r_shell_moy_tab), marker='.',color=colorspts[0],label=ur'$<V_r>$')
    plt.semilogx(shells_au,V_r_rms_shell_moy_tab, marker='.',color=colorspts[1],label=ur'$V_{r, rms}$')
    plt.semilogx(shells_au,V_paral_shell_moy_tab, marker='.',color=colorspts[2],label=ur'$V_{par}$')
    plt.semilogx(shells_au, sigma_shell_moy_tab, marker='.',color=colorspts[3],label=ur'$\sigma$')
    plt.xlabel(ur'Rayon ($AU$)')
    plt.ylabel(ur"V  $(km.s^{-1})$")
    plt.legend(loc='best')

    if save_fig==True:
        plt.savefig('/home/averliat/these/analyses/Gradients_vitesse/Analyses_champ_vitesse/V_infall/V_infall_'+simu+'_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.pdf', bbox_inches='tight')



    plt.figure()
    plt.semilogx(shells_au, np.abs(V_r_shell_moy_tab)/sigma_shell_moy_tab, marker='.',color=colorspts[0],label=ur'$<V_r>$ / $\sigma$')
    plt.semilogx(shells_au, np.abs(V_r_shell_moy_tab)/V_paral_shell_moy_tab, marker='.',color=colorspts[1],label=ur'$<V_r>$ / $V_{par}$')
    plt.semilogx(shells_au, V_r_rms_shell_moy_tab/sigma_shell_moy_tab, marker='.',color=colorspts[2],label=ur'$V_{r, rms}$ / $\sigma$')
    plt.semilogx(shells_au, V_r_rms_shell_moy_tab/V_paral_shell_moy_tab, marker='.',color=colorspts[3],label=ur'$V_{r, rms}$ / $V_{par}$')
    plt.xlabel(ur'Rayon ($AU$)')
    plt.ylabel(ur"Rapport")
    plt.legend(loc='best')

    if save_fig==True:
        plt.savefig('/home/averliat/these/analyses/Gradients_vitesse/Analyses_champ_vitesse/V_infall/Rapport_v_infall_'+simu+'_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.pdf', bbox_inches='tight')



if __name__=='__main__':
    trace_Vinfall_gradient_vitesse(simu,tag,output,R_min,R_max,dr,legend,marker,save_fig)

