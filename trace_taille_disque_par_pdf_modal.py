# -*- coding: utf-8 -*-
import numpy as np

import pymses
import sys
import glob as glob
import os

from pymses.filters import CellsToPoints
from pymses.utils import constants as cst
import matplotlib.pyplot as plt
plt.ion() #Mis 2 lignes plus bas

import h5py


#Remise a zero des parametres de matplotlib car ils sont modifies par pipeline_taille_disque_par_pdf.py ! Inutile si on quitte ipython entre les 2...
#plt.rcParams.update(plt.rcParamsDefault)
#plt.ion()

if __name__ == '__main__':
#-----------------------------------------------------------
#Noms des simulations et caracteristiques du calcul du rayon
#-----------------------------------------------------------
    simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire_bigbox_50pourc_sink_seuil_haut_MHD_lr'
    tag = '50_shr_bigbox'

    output_min = 'None'
    output_max = 291

    seuil_rho = 1e-10




def trace_taille_disque(simu,tag,tag2,output_min,output_max,color,t1=0.,legend='',marker='.',seuil_rho=1e-10,output_frag='None'):
#------------------
#Differents chemins
#------------------
    path_analyses='/home/averliat/these/analyses/'+simu+'/'

    if output_min=='None':
        ref=np.loadtxt(path_analyses+'t0_seuil_rho_'+str(seuil_rho)+'.txt')
        output_0 = int(ref[0])
        time_0 = ref[1]
        output_min = output_0
        file_save = 'Rayon_disque_par_pdf_modal_bin4_'+str(output_min)+'_'+str(output_max)+'.hdf5'



#------------------------------------
#Chargement des differentes quantites
#------------------------------------
    h5f = h5py.File(path_analyses+file_save, 'r')

    mean_rad_loc_tab=h5f['mean_rad_loc_tab'][:]
    time_tab=h5f['time_tab'][:]
    time_cor_tab=h5f['time_cor_tab'][:]
    min_rad_loc_tab=h5f['min_rad_loc_tab'][:]
    max_rad_loc_tab=h5f['max_rad_loc_tab'][:]
    numero_output=h5f['numero_output'][:]
    mass_disk_tab=h5f['mass_disk_tab'][:]
    log_rho_disk_tab=h5f['log_rho_disk_tab'][:]
    mag_mean_broad_tab=h5f['mag_mean_broad_tab'][:]
    median_rad_loc_tab=h5f['median_rad_loc_tab'][:]
    rad_estim_tab=h5f['rad_estim_tab'][:]

    h5f.close()



#------------------------------------------------------------------
#Redimensionnement tableaux pour la simu 50_mhr (disque trop grand)
#------------------------------------------------------------------
    if tag=='50_mhr':
        time_tab=time_tab[0:len(mean_rad_loc_tab)-7]
        numero_output=numero_output[0:len(mean_rad_loc_tab)-7]
        time_cor_tab=time_cor_tab[0:len(mean_rad_loc_tab)-7]
        mean_rad_loc_tab=mean_rad_loc_tab[0:len(mean_rad_loc_tab)-7]

    if tag=='50_hr':
        time_tab=time_tab[0:84]
        numero_output=numero_output[0:84]
        time_cor_tab=time_cor_tab[0:84]
        mean_rad_loc_tab=mean_rad_loc_tab[0:84]



#--------------------------------------------------
#Filtre mediane glissante pour avoir moins de bruit
#--------------------------------------------------
    half_window=20
    if 'rot' in tag2:#tag2=='20_rot1':
        half_window=2
        if 'MHD' in tag2:#tag2=='20_rot1':
            half_window=5
    def filtre_median(tab):
        nbr_pnts=len(tab)-2*half_window
        tab_filtre=np.zeros(nbr_pnts)
        for i in range(nbr_pnts):
            ind=i+half_window+1
            tab_filtre[i]=np.median(tab[ind-half_window-1:ind+half_window])
        return tab_filtre

    def filtre_mean(tab):
        nbr_pnts=len(tab)-2*half_window
        tab_filtre=np.zeros(nbr_pnts)
        for i in range(nbr_pnts):
            ind=i+half_window+1
            tab_filtre[i]=np.mean(tab[ind-half_window-1:ind+half_window])
        return tab_filtre

    def filtre_time(tab):
        tab_filtre=tab[half_window:-half_window]
        return tab_filtre



#------------------
#Debut de la figure
#------------------

    if output_frag != 'None':
        frag=output_frag-output_min+1
    else:
        frag=output_max-output_min+1

    mean_rad_loc_tab=mean_rad_loc_tab[0:frag]
    time_tab=time_tab[0:frag]
    time_cor_tab=time_cor_tab[0:frag]
    min_rad_loc_tab=min_rad_loc_tab[0:frag]
    max_rad_loc_tab=max_rad_loc_tab[0:frag]
    numero_output=numero_output[0:frag]
    mass_disk_tab=mass_disk_tab[0:frag]
    log_rho_disk_tab=log_rho_disk_tab[0:frag]
    mag_mean_broad_tab=mag_mean_broad_tab[0:frag]
    median_rad_loc_tab=median_rad_loc_tab[0:frag]
    rad_estim_tab=rad_estim_tab[0:frag]
#cmappts = plt.get_cmap('magma')
#colorspts = [cmappts(i) for i in np.linspace(0.1,0.9,7)]

    '''
    plt.figure(1)
    plt.plot(numero_output,mean_rad_loc_tab, marker=marker, label=legend+'  - '+tag)
    plt.xlabel(ur'Numero output')
    plt.ylabel(ur'Rayon moyen du disque ($AU$)')
    plt.legend(loc='best')

    plt.figure(2)
    plt.plot(time_tab-t1,mean_rad_loc_tab, marker=marker, label=legend+'  - '+tag)
    plt.xlabel(ur'Temps ($Myr$)')
    plt.ylabel(ur'Rayon moyen du disque ($AU$)')
    plt.legend(loc='best')

    plt.figure(3)
    plt.plot(time_tab,mean_rad_loc_tab, marker=marker, label=legend+'  - '+tag)
    plt.xlabel(ur'Temps ($Myr$)')
    plt.ylabel(ur'Rayon moyen du disque ($AU$)')
    plt.legend(loc='best')
    '''
    plt.figure(4)
    #plt.plot(time_cor_tab,rad_estim_tab,linestyle='None',marker='.', label=tag2)
    #plt.plot(filtre_time(time_cor_tab),filtre_median(rad_estim_tab),linestyle='None',marker='.', label=tag2)
    #plt.plot(filtre_time(time_cor_tab),filtre_mean(rad_estim_tab), color=color, label=tag2)
    #plt.plot(filtre_time(numero_output),filtre_median(rad_estim_tab), label=tag2)
    plt.plot(numero_output,rad_estim_tab,marker='.', label=tag2)
    #plt.plot(time_cor_tab,median_rad_loc_tab, linestyle='None',color='r',marker='+')
    #plt.plot(time_cor_tab,mean_rad_loc_tab, linestyle='None',color='g',marker='x')
    plt.xlabel(ur'Temps corrigé ($Myr$)')
    plt.ylabel(ur'Rayon moyen du disque ($AU$)')
    plt.legend(loc='best')

    plt.figure(5)
    plt.plot(filtre_time(time_cor_tab),filtre_median(rad_estim_tab), color=color, label=tag2)
    plt.xlabel(ur'Temps corrigé ($Myr$)')
    plt.ylabel(ur'Rayon moyen du disque ($AU$)')
    plt.legend(loc='best')

    '''
    plt.figure()
    plt.plot(time_tab,mass_disk_tab, marker='.',color='green', label=ur'$R_{disk}$')
    plt.xlabel(ur'Temps ($Myr$)')
    plt.ylabel(ur'Masse du disque ($Msun$)')
    plt.legend(loc='best')

    plt.figure()
    plt.plot(time_tab,min_rad_loc_tab, marker='.',color='red', label=ur'$R_{disk}$')
    plt.xlabel(ur'Temps ($Myr$)')
    plt.ylabel(ur'Rayon min du disque ($AU$)')
    plt.legend(loc='best')

    plt.figure()
    plt.plot(time_tab,max_rad_loc_tab, marker='.',color='darkorange', label=ur'$R_{disk}$')
    plt.xlabel(ur'Temps ($Myr$)')
    plt.ylabel(ur'Rayon max du disque ($AU$)')
    plt.legend(loc='best')

    plt.figure()
    plt.plot(time_tab,log_rho_disk_tab, marker='.',color='black', label=ur'$R_{disk}$')
    plt.xlabel(ur'Temps ($Myr$)')
    plt.ylabel(ur'log_rho_tab_disk')
    plt.legend(loc='best')
    '''
    plt.show()



if __name__ == '__main__':
    trace_taille_disque(simu,tag,tag,output_min,output_max,'b')







