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



if __name__=='__main__':
#-----------------------------------------------------------
#Noms des simulations et caracteristiques du calcul du rayon
#-----------------------------------------------------------
    tag='_shr_bigbox_0pourc'
    legend='8,  10*40'
    marker='.'

    simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire'+tag

    output_min = 1
    output_max = 41

    width_bar = 1

    hist=True
    diff_abs=True
    diff_relat=True
    moyenne_glissante=True

    savefig=False



def trace_moment_cin_all(simu,tag,output_min,output_max,width_bar,hist,diff_abs,diff_relat,legend,marker,t1=0.,moyenne_glissante=True,savefig=False):
#------------------
#Differents chemins
#------------------
    path_analyses='/home/averliat/these/analyses/'+simu+'/'
    file_save = 'Moment_cin_observateurs_output'+str(output_min)+'_to_'+str(output_max)+'.hdf5'
    path_data = '/drf/projets/alfven-data/averliat/'+simu+'/'



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
    numero_output=h5f['numero_output'][:]

    h5f.close()



#------------------
#Debut de la figure
#------------------
    cmappts = plt.get_cmap('magma')
    colorspts = [cmappts(i) for i in np.linspace(0.1,0.9,6)]


    '''
    plt.figure()
    plt.plot(tab_x,norm_center_C_AU, marker='.',color=colorspts[2], label=tag)
    plt.xlabel(ur'Output')
    plt.ylabel(ur"Distance ($AU$)")
    plt.legend(loc='best')
    '''
    '''
    plt.figure()
    plt.semilogy(tab_x,norm_GC_AU, marker='.',color=colorspts[2])#, label=tag)
    plt.xlabel(ur'Output')
    plt.ylabel(ur"Distance ($AU$)")
    plt.legend(loc='best')

    plt.figure()
    plt.plot(tab_x,norm_center_G_AU, marker='.',color=colorspts[2])#, label=tag)
    plt.xlabel(ur'Output')
    plt.ylabel(ur"Distance ($AU$)")
    plt.legend(loc='best')
    '''


    if hist==True:
#Histogramme
        plt.figure()
        mask=(norm_mom_cin_obs_par_integ-norm_mom_cin_obs)>0
        mask_inv=np.logical_not(mask)
        plt.bar(numero_output[mask], norm_mom_cin_obs_par_integ[mask], width=width_bar, color=colorspts[5], label=ur'Norme du moment cinétique réel '+tag)
        plt.bar(numero_output[mask_inv], norm_mom_cin_obs[mask_inv], width=width_bar, color=colorspts[1], label=ur'Norme du moment cinétique analytique '+tag)
        plt.bar(numero_output[mask_inv], norm_mom_cin_obs_par_integ[mask_inv], width=width_bar, color='w', edgecolor='black', linewidth=0.4)
        plt.bar(numero_output[mask], norm_mom_cin_obs[mask], width=width_bar, color='w', edgecolor='black', linewidth=0.4)
        plt.xlabel(ur'Output')
        plt.ylabel(ur"Moment cinetique ($kg.m^2.s^{-1}$)")
        plt.legend(loc='best')

        #plt.yscale('log')
        #plt.xlim((0,53))
        #plt.ylim((0,1.2e47))

    if diff_abs==True:
        plt.figure(10)
        plt.plot(simulation_time-t1, np.abs(norm_mom_cin_obs_par_integ-norm_mom_cin_obs), marker=marker, label=legend+'  - '+tag, linestyle='None')
        plt.xlabel(ur'Temps (Myr)')
        plt.ylabel(ur"Différence absolue ($kg.m^2.s^{-1}$)")
        plt.legend(loc='best')

    if diff_relat==True:
        plt.figure(11)
        plt.plot(simulation_time-t1, np.abs((norm_mom_cin_obs_par_integ-norm_mom_cin_obs)/norm_mom_cin_obs_par_integ), marker=marker, label=legend+'  - '+tag, linestyle='None')
        plt.xlabel(ur'Temps (Myr)')
        plt.ylabel(ur"Différence relative")
        plt.xlim((0,0.009))
        plt.legend(loc='best')


    erreur=np.abs((norm_mom_cin_obs_par_integ-norm_mom_cin_obs)/norm_mom_cin_obs_par_integ)
    erreur_moy=np.mean(np.nan_to_num(erreur[np.where(simulation_time-t1>0)]))
    #plt.figure(12)
    #plt.plot(tag, erreur_moy,marker=marker)
    #plt.ylabel(ur"Moyenne de l'erreur")




    if moyenne_glissante==True:
        erreur=np.nan_to_num(erreur[np.where(simulation_time-t1>0)])
        time_pour_liss=np.nan_to_num((simulation_time-t1)[np.where(simulation_time-t1>0)])
        nbr_pts=len(erreur)/10
        grp=len(erreur)/nbr_pts
        erreur_lisse=np.zeros(nbr_pts)
        time_lisse=np.zeros(nbr_pts)
        for i in range(nbr_pts):
            erreur_lisse[i]=np.mean(erreur[i*grp:(i+1)*grp])
            time_lisse[i]=np.mean(time_pour_liss[i*grp:(i+1)*grp])
        plt.figure(13)
        plt.plot(time_lisse*1e3,erreur_lisse, marker=marker,color=colorspts[int(tag)/10-1], label=ur'$\varepsilon=$ '+tag+ur'\%')
        plt.xlabel(ur'Time (kyr)')
        plt.ylabel(ur"$\left| \frac{\Delta \sigma}{\sigma} \right|$")
        plt.legend(loc='upper left')
        plt.xlim((0,9))
        plt.ylim((0,0.15))

        if savefig==True:
            plt.tight_layout(pad=0.1) #pad en inch si besoin
            plt.savefig(path_analyses+file_save+'.pdf')#, bbox_inches='tight')



    plt.show()



    '''
    plt.figure(1)
    plt.savefig('../Comparaisons/Verification_moment_physique_propre/Moment_cinetique_analytique_et_reel_50vhr_synchro_output1_123.pdf', bbox_inches='tight')

    plt.figure(2)
    plt.savefig('../Comparaisons/Verification_moment_physique_propre/Moment_cinetique_analytique_et_reel_50hr_synchro_output1_92.pdf', bbox_inches='tight')

    plt.figure(3)
    plt.savefig('../Comparaisons/Verification_moment_physique_propre/Moment_cinetique_analytique_et_reel_10vhr_synchro_output1_166.pdf', bbox_inches='tight')

    plt.figure(4)
    plt.savefig('../Comparaisons/Verification_moment_physique_propre/Moment_cinetique_analytique_et_reel_10hr_synchro_output1_118.pdf', bbox_inches='tight')

    plt.figure(5)
    plt.savefig('../Comparaisons/Verification_moment_physique_propre/Moment_cinetique_analytique_et_reel_hr2_synchro_output1_102.pdf', bbox_inches='tight')
    '''



if __name__=='__main__':
    trace_moment_cin_all(simu,tag,output_min,output_max,width_bar,hist,diff_abs,diff_relat,legend,marker,moyenne_glissante)

