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
    tag='lr'
    simu = 'M30_mu7_mach2_'+tag

    output = 60

    R_min = 1
    R_max = 'all'#10000
    dr = 100


    legend='7,  10*40'  #Inutile pour le moment
    marker='.'
    color='black'

    save_fig = False





    tag2='lr_jets_crit0.17'
    simu2 = 'M30_mu7_mach2_'+tag2

    output2 = 60

    R_min2 = 1
    R_max2 = 'all'#10000
    dr2 = 100


    legend2='7,  10*40'  #Inutile pour le moment
    marker2='.'
    color2='grey'





    tag3='lr_jets'
    simu3 = 'M30_mu7_mach2_'+tag3

    output3 = 60

    R_min3 = 1
    R_max3 = 'all'#10000
    dr3 = 100


    legend3='7,  10*40'  #Inutile pour le moment
    marker3='.'
    color3='b'
    




    tag4='lr_jets_co30'
    simu4 = 'M30_mu7_mach2_'+tag4

    output4 = 60

    R_min4 = 1
    R_max4 = 'all'#10000
    dr4 = 100


    legend4='7,  10*40'  #Inutile pour le moment
    marker4='.'
    color4='orange'





    tag5='lr_jets_vjets66'
    simu5 = 'M30_mu7_mach2_'+tag5

    output5 = 60

    R_min5 = 1
    R_max5 = 'all'#10000
    dr5 = 100


    legend5='7,  10*40'  #Inutile pour le moment
    marker5='.'
    color5='y'





    tag6='lr_jets_co30_vjets66'
    simu6 = 'M30_mu7_mach2_'+tag6

    output6 = 60

    R_min6 = 1
    R_max6 = 'all'#10000
    dr6 = 100


    legend6='7,  10*40'  #Inutile pour le moment
    marker6='.'
    color6='purple'




def trace_Vinfall_gradient_vitesse(simu,tag,output,R_min,R_max,dr,legend,marker,save_fig,color='b'):
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


    #plt.figure()
    #plt.semilogx(shells_au,np.abs(V_r_shell_moy_tab), marker='.',color=colorspts[0],label=ur'$<V_r>$')
    #plt.semilogx(shells_au,V_r_rms_shell_moy_tab, marker='.',color=colorspts[1],label=ur'$V_{r, rms}$')
    #plt.semilogx(shells_au,V_paral_shell_moy_tab, marker='.',color=colorspts[2],label=ur'$V_{par}$')
    #plt.semilogx(shells_au, sigma_shell_moy_tab, marker='.',color=colorspts[3],label=ur'$\sigma$')
    #plt.xlabel(ur'Rayon ($AU$)')
    #plt.ylabel(ur"V  $(km.s^{-1})$")
    #plt.legend(loc='best')

    #if save_fig==True:
    #    plt.savefig('/home/averliat/these/analyses/Gradients_vitesse/Analyses_champ_vitesse/V_infall/V_infall_'+simu+'_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.pdf', bbox_inches='tight')


    '''
    plt.figure(10)
    plt.semilogx(shells_au, np.abs(V_r_shell_moy_tab)/sigma_shell_moy_tab, marker='.',color=colorspts[0],label=ur'$<V_r>$ / $\sigma$')
    plt.semilogx(shells_au, np.abs(V_r_shell_moy_tab)/V_paral_shell_moy_tab, marker='.',color=colorspts[1],label=ur'$<V_r>$ / $V_{par}$')
    plt.semilogx(shells_au, V_r_rms_shell_moy_tab/sigma_shell_moy_tab, marker='.',color=colorspts[2],label=ur'$V_{r, rms}$ / $\sigma$')
    plt.semilogx(shells_au, V_r_rms_shell_moy_tab/V_paral_shell_moy_tab, marker='.',color=colorspts[3],label=ur'$V_{r, rms}$ / $V_{par}$')
    plt.xlabel(ur'Rayon ($AU$)')
    plt.ylabel(ur"Rapport")
    plt.legend(loc='best')
    '''

    plt.figure(10)
    plt.semilogx(shells_au, sigma_shell_moy_tab, marker='.',color=color,label=tag)
    plt.xlabel(ur'Rayon ($AU$)')
    plt.ylabel(ur"$\sigma_v$ ($km.s^{-1}$)")
    plt.legend(loc='best')

    if save_fig==True:
        plt.savefig('/home/averliat/these/analyses/Gradients_vitesse/Analyses_champ_vitesse/V_infall/Rapport_v_infall_'+simu+'_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.pdf', bbox_inches='tight')



if __name__=='__main__':
    trace_Vinfall_gradient_vitesse(simu,tag,output,R_min,R_max,dr,legend,marker,save_fig,color)
    trace_Vinfall_gradient_vitesse(simu2,tag2,output2,R_min2,R_max2,dr2,legend2,marker2,save_fig,color2)
    trace_Vinfall_gradient_vitesse(simu3,tag3,output3,R_min3,R_max3,dr3,legend3,marker3,save_fig,color3)
    trace_Vinfall_gradient_vitesse(simu4,tag4,output4,R_min4,R_max4,dr4,legend4,marker4,save_fig,color4)
    trace_Vinfall_gradient_vitesse(simu5,tag5,output5,R_min5,R_max5,dr5,legend5,marker5,save_fig,color5)
    trace_Vinfall_gradient_vitesse(simu6,tag6,output6,R_min6,R_max6,dr6,legend6,marker6,save_fig,color6)

