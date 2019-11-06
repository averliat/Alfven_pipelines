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
import matplotlib.patches as mpatches

import h5py

plt.style.use("pdf")
plt.style.use("aanda_modif")

#Remise a zero des parametres de matplotlib car ils sont modifies par pipeline_taille_disque_par_pdf.py ! Inutile si on quitte ipython entre les 2...
#plt.rcParams.update(plt.rcParamsDefault)
#plt.ion()

if __name__ == '__main__':
#-----------------------------------------------------------
#Noms des simulations et caracteristiques du calcul du rayon
#-----------------------------------------------------------
    simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire_bigbox_10pourc_sink_seuil_haut'
    tag = '50_shr_bigbox'

    output_min = 'None'
    output_max = 440

    seuil_rho = 1e-10




def trace_taille_disque(simu,tag,output_min,output_max,t1=0.,legend='',marker='.',seuil_rho=1e-10,ind=0,tagtot=1):
#------------------
#Differents chemins
#------------------
    path_analyses='/home/averliat/these/analyses/'+simu+'/'

    if output_min=='None':
        ref=np.loadtxt(path_analyses+'t0_seuil_rho_'+str(seuil_rho)+'.txt')
        output_0 = int(ref[0])
        time_0 = ref[1]
        output_min = output_0
        file_save = 'Rayon_disque_par_pdf'+str(output_min)+'_'+str(output_max)+'.hdf5'



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



#------------------
#Debut de la figure
#------------------
    cmappts = plt.get_cmap('plasma')
    #colorspts = [cmappts(i) for i in np.linspace(0.1,0.9,7)]
    #colorspts = [cmappts(i) for i in np.linspace(0.1,0.82,tagtot)]#0.1-0.82 pour 4 simus
    colorspts = [cmappts(i) for i in np.linspace(0,1,tagtot)]#0.1-0.82 pour 4 simus

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
    

    plt.figure(4)
    plt.plot(time_cor_tab*1e3,mean_rad_loc_tab, marker='.', label=tag,color=colorspts[ind])
    plt.plot(time_cor_tab*1e3,max_rad_loc_tab, marker='.', linestyle='dashed',label=tag,color=colorspts[ind])
    plt.plot(time_cor_tab*1e3,min_rad_loc_tab, marker='.', linestyle='dotted',label=tag,color=colorspts[ind])
    plt.xlabel(ur'Time (kyr)')
    plt.ylabel(ur'Disk radius (AU)')
    plt.legend(loc='best')
    '''




    plt.figure(5)
    ax=plt.gca()

    nbr_pts=len(mean_rad_loc_tab)/20
    grp=len(mean_rad_loc_tab)/nbr_pts

    mean_rad_loc_tab_lisse=np.zeros(nbr_pts)
    max_rad_loc_tab_lisse=np.zeros(nbr_pts)
    min_rad_loc_tab_lisse=np.zeros(nbr_pts)
    time_cor_tab_lisse=np.zeros(nbr_pts)
    for i in range(nbr_pts):
        mean_rad_loc_tab_lisse[i]=np.mean(mean_rad_loc_tab[i*grp:(i+1)*grp])
        max_rad_loc_tab_lisse[i]=np.mean(max_rad_loc_tab[i*grp:(i+1)*grp])
        min_rad_loc_tab_lisse[i]=np.mean(min_rad_loc_tab[i*grp:(i+1)*grp])
        time_cor_tab_lisse[i]=np.mean(time_cor_tab[i*grp:(i+1)*grp])



    ax.plot(time_cor_tab_lisse*1e3,mean_rad_loc_tab_lisse, label=tag,color=colorspts[ind])
    #ax.plot(time_cor_tab_lisse*1e3,max_rad_loc_tab_lisse, linestyle='dashed',color=colorspts[ind])
    #ax.plot(time_cor_tab_lisse*1e3,min_rad_loc_tab_lisse, linestyle='dotted',color=colorspts[ind])

    if ind==tagtot-1:
        #Get artists and labels for legend and chose which ones to display
        handles, labels = ax.get_legend_handles_labels()
        display = range(0,tagtot)#(0,1,2,3)

        #Create custom artists
        meanArtist = plt.Line2D((0,1),(0,0), color='dimgrey', linestyle='solid')
        maxArtist = plt.Line2D((0,1),(0,0), color='dimgrey',linestyle='dashed')
        minArtist = plt.Line2D((0,1),(0,0), color='dimgrey',linestyle='dotted')

        #Create custom artists for colors
        handles=[]
        for j in range(tagtot):
            handles.append(mpatches.Patch(color=colorspts[j]))#, label='The red data')


        plt.xlabel(ur'Time (kyr)')
        plt.ylabel(ur'Disk radius (AU)')
        plt.legend(loc='best')


        #Create legend from custom artist/label lists
        color_legend=ax.legend([handle for i,handle in enumerate(handles) if i in display],
                [label for i,label in enumerate(labels) if i in display],loc='upper left')#,ncol=2)

        ax.add_artist(color_legend)

        ax.legend([meanArtist,maxArtist,minArtist],
                [ur'$R_{\text{mean}}$',ur'$R_{\text{max}}$',ur'$R_{\text{min}}$'],loc='upper right')




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
    trace_taille_disque(simu,tag,output_min,output_max)







