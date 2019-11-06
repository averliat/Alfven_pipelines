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
    simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire_bigbox_20pourc_sink_seuil_haut_rot0.25'
    tag = '50_shr_bigbox'

    output_min = 'None'
    output_max = 40

    seuil_rho = 1e-10




def trace_taille_disque(simu,tag,tag2,output_min,output_max,color,indice,tagfull=[1],t1=0.,legend='',marker='.',seuil_rho=1e-10,output_frag='None'):
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



    fig=plt.figure(5)
    ax=plt.gca()
    ax.plot(filtre_time(time_cor_tab)*1e3,filtre_median(rad_estim_tab), color=color, label=legend)
    plt.xlabel(r'Time (kyr)')
    plt.ylabel(r'Disk radius (AU)')
    #plt.legend(loc='best')

    
    if indice==len(tagfull)-1:
        #Get artists and labels for legend and chose which ones to display
        handles, labels = ax.get_legend_handles_labels()
        #display_rot = [5,6,7,8]
        display_rot = [3,4,5,6]
        #display_MHD = [3,4]
        display_norot = [0,1,2]


        #Create custom artists
        #simArtist = plt.Line2D((0,1),(0,0), color='dimgrey', linestyle='solid')
        #anyArtist = plt.Line2D((0,1),(0,0), color='dimgrey',linestyle='dashed')

        #Create custom artists for color
        #handles=[]
        #for j in range(len(tag)):
        #    handles.append(plt.Line2D((0,1),(0,0), color=colors, linestyle='solid'))#, label='The red data')



        #Create legend from custom artist/label lists
        norot_legend=ax.legend([handle for i,handle in enumerate(handles) if i in display_norot],
                [label for i,label in enumerate(labels) if i in display_norot],loc='lower right',ncol=1,fancybox=True)#,framealpha=0.1)
        ax.add_artist(norot_legend)

        rot_legend=ax.legend([handle for i,handle in enumerate(handles) if i in display_rot],
                [label for i,label in enumerate(labels) if i in display_rot],loc='upper left',ncol=1,fancybox=True)#,framealpha=0.1)
        #ax.add_artist(rot_legend)
        #ax.add_artist(rot_legend)
        #rot_legend.get_frame().set_alpha(0.1)
        #rot_legend.set_frame_on(True)
        #rot_legend.get_frame().set_edgecolor('black')

        #MHD_legend=ax.legend([handle for i,handle in enumerate(handles) if i in display_MHD],
                #[label for i,label in enumerate(labels) if i in display_MHD],loc='lower right',ncol=1)#,fancybox=True)#,framealpha=0.1)
    

    #plt.xlim((-2,40))
    #plt.xlim((0,8))
    #plt.ylim((0,100))


    plt.show()



if __name__ == '__main__':
    trace_taille_disque(simu,tag,tag,output_min,output_max,'b')







