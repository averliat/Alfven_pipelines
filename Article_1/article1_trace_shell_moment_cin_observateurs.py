# -*- coding: utf-8 -*-
import numpy as np

import pymses
import sys
import glob as glob
import os

from pymses.filters import CellsToPoints
from pymses.utils import constants as cst
import matplotlib.pyplot as plt
from matplotlib import gridspec

plt.ion()

import h5py
#import matplotlib_pdf

plt.style.use("pdf")
plt.style.use("aanda_modif")




#params = {
#        'axes.titlesize': 10, 
#        'axes.labelsize': 9, 
#        'font.size': 10, 
#        'figure.titlesize': 10, 
#        'legend.fontsize': 8, 
#        'xtick.labelsize': 8, 
#        'ytick.labelsize': 8}
#plt.rcParams.update(params)
#plt.rc('text',usetex=True)
#font = {'family':'serif'}
#plt.rc('font',**font)
# #plt.rc('legend',**{'fontsize':10})
#plt.rcParams['text.latex.preamble']=[r'\usepackage{amsmath}']



if __name__=='__main__':
#-----------------------------------------------------------
#Noms des simulations et caracteristiques du calcul du rayon
#-----------------------------------------------------------
    tag=''#_shr_bigbox_50pourc'
    legend='8,  10*40'
    marker='.'

    simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire_shr_bigbox_50pourc'#+tag

    output_min = 1
    output_max = 110

    width_bar = 1

    hist=False
    diff_abs=False
    diff_relat=False
    moyenne_glissante=False
    article1=True

    save=True
    save_plusieurs=False



def trace_moment_cin_all(simu,tag,output_min,output_max,width_bar,hist,diff_abs,diff_relat,legend,marker,t1=0.,moyenne_glissante=True,save=False,article1=True,save_plusieurs=False):
#------------------
#Differents chemins
#------------------
    path_analyses='/home/averliat/these/analyses/'+simu+'/'
    file_save = 'Moment_cin_observateurs_output'+str(output_min)+'_to_'+str(output_max)+'.hdf5'
    path_data = '/drf/projets/alfven-data/averliat/'+simu+'/'

    name_save = simu+'_moment_cinetique_analytique_et_reel_output_'+str(output_min)+'_'+str(output_max)+'.pdf'
    path_save = path_analyses+name_save
    path_save_plusieurs = '/home/averliat/these/analyses/article1_figures/'



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
    colorspts = [cmappts(i) for i in np.linspace(0.1,0.9,7)]


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
        plt.bar(numero_output[mask], norm_mom_cin_obs_par_integ[mask], width=width_bar, color=colorspts[5], label=ur'Norm of real momentum')
        plt.bar(numero_output[mask_inv], norm_mom_cin_obs[mask_inv], width=width_bar, color=colorspts[1], label=ur'Norm of analytical momentum')
        plt.bar(numero_output[mask_inv], norm_mom_cin_obs_par_integ[mask_inv], width=width_bar, color='w', edgecolor='black', linewidth=0.4)
        plt.bar(numero_output[mask], norm_mom_cin_obs[mask], width=width_bar, color='w', edgecolor='black', linewidth=0.4)
        plt.xlabel(ur'Output')
        plt.ylabel(ur"Angular momentum ($kg.m^2.s^{-1}$)")
        plt.legend(loc='best')
        plt.tight_layout(pad=0.1) #pad en inch si besoin
        if save==True:
            plt.savefig(path_save, bbox_inches='tight')


        #plt.yscale('log')
        #plt.xlim((0,53))
        #plt.ylim((0,1.2e47))

    if diff_abs==True:
        plt.figure(10)
        plt.plot(simulation_time-t1, np.abs(norm_mom_cin_obs_par_integ-norm_mom_cin_obs), marker=marker, label=legend+'  - '+tag, linestyle='None')
        plt.xlabel(ur'Temps (Myr)')
        plt.ylabel(ur"Différence absolue ($kg.m^2.s^{-1}$)")
        plt.legend(loc='best')
        plt.tight_layout(pad=0.1) #pad en inch si besoin

    if diff_relat==True:
        plt.figure(11)
        plt.plot(numero_output, np.abs((norm_mom_cin_obs_par_integ-norm_mom_cin_obs)/norm_mom_cin_obs_par_integ), marker=marker, label=legend, linestyle='None')
        plt.xlabel(ur'Temps (Myr)')
        plt.ylabel(ur"Différence relative")
        plt.legend(loc='best')
        plt.tight_layout(pad=0.1) #pad en inch si besoin
        if save==True and save_plusieurs==True:
            plt.tight_layout(pad=0.1) #pad en inch si besoin
            plt.savefig(path_save_plusieurs+'test_diff_relat.pdf')


    erreur=np.abs((norm_mom_cin_obs_par_integ-norm_mom_cin_obs)/norm_mom_cin_obs_par_integ)
    erreur_moy=np.mean(np.nan_to_num(erreur[np.where(simulation_time-t1>0)]))
    #plt.figure(12)
    #plt.plot(tag, erreur_moy,marker=marker)
    #plt.ylabel(ur"Moyenne de l'erreur")




    if moyenne_glissante==True:
        erreur=np.nan_to_num(erreur[np.where(simulation_time-t1>0)])
        time_pour_liss=np.nan_to_num((simulation_time-t1)[np.where(simulation_time-t1>0)])
        #grp=10
        #nbr_pts=len(erreur)/10
        #erreur_lisse=np.zeros(nbr_pts)
        #time_lisse=np.zeros(nbr_pts)
        #for i in range(nbr_pts):
        #    erreur_lisse[i]=np.mean(erreur[i*grp:(i+1)*grp])
        #    time_lisse[i]=np.mean(time_pour_liss[i*grp:(i+1)*grp])
        nbr_liss=21 #doit etre impair
        nbr_pts_effectif=len(erreur)-2*nbr_liss/2 #les premiers et derniers points
        erreur_lisse=np.zeros(nbr_pts_effectif)
        time_lisse=np.zeros(nbr_pts_effectif)
        for i in range(nbr_pts_effectif):
            pt=i+2
            erreur_lisse[i] = np.mean(erreur[pt-nbr_liss/2:pt+nbr_liss/2+1])
            time_lisse[i] = time_pour_liss[pt]

        plt.figure(13)
        plt.semilogy(time_lisse,erreur_lisse, marker=marker, label=legend)
        plt.xlabel(ur'Temps corrigé lissé (Myr)')
        plt.ylabel(ur"Différence relative lissée")
        plt.legend(loc='upper left')
        plt.tight_layout(pad=0.1) #pad en inch si besoin
        #plt.ylim((0,1.3))
        if save==True and save_plusieurs==True:
            #plt.tight_layout(pad=0.1) #pad en inch si besoin
            plt.savefig(path_save_plusieurs+'test_moy_gliss2.pdf')






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





    if article1==True:

        f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=False)#, figsize=(3.5433,2.6575))

        #plt.rc('text', usetex=True)
        gs  = gridspec.GridSpec(2, 1, height_ratios=[7, 2])
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])

        #mask=(norm_mom_cin_obs_par_integ-norm_mom_cin_obs)>0
        #mask_inv=np.logical_not(mask)
        #ax1.bar(numero_output[mask], norm_mom_cin_obs_par_integ[mask], width=width_bar, color=colorspts[5], label=ur'Norm of real momentum')
        #ax1.bar(numero_output[mask_inv], norm_mom_cin_obs[mask_inv], width=width_bar, color=colorspts[1], label=ur'Norm of analytical momentum')
        #ax1.bar(numero_output[mask_inv], norm_mom_cin_obs_par_integ[mask_inv], width=width_bar, color='w', edgecolor='black', linewidth=0.4)
        #ax1.bar(numero_output[mask], norm_mom_cin_obs[mask], width=width_bar, color='w', edgecolor='black', linewidth=0.4)
        #ax1.yaxis.set_ticks(np.array([0.2,0.4,0.6,0.8,1.0,1.2,1.4])*1e48)
        #ax1.set_xlabel(ur'Output')
        #ax1.set_ylabel(ur"Angular momentum ($kg.m^2.s^{-1}$)",labelpad = 28)
        #ax1.legend(loc='best')
        

        ax1.plot(numero_output, norm_mom_cin_obs_par_integ, marker=marker, color='Midnightblue',linestyle='None')
        ax1.set_xlabel(ur'Output')
        ax1.set_ylabel(ur"\begin{center} Angular momentum \\ \indent ($kg.m^2.s^{-1}$) \end{center}",labelpad = 6.8)
        #ax1.legend(loc='best')


        error_relat=(norm_mom_cin_obs_par_integ-norm_mom_cin_obs)/norm_mom_cin_obs_par_integ

        ax2.plot(numero_output, error_relat, marker=marker, linestyle='None',color='Darkmagenta')
        ax2.set_xlabel(ur'Output')
        ax2.set_ylabel(ur'$\frac{\Delta \sigma}{\sigma}$',rotation='horizontal',labelpad=10)
        ax2.set_ylim((1.3*np.sort(np.nan_to_num(error_relat))[0],-1.3*np.sort(np.nan_to_num(error_relat))[0]))

        plt.tight_layout(pad=0.1) #pad en inch si besoin
        f.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        ax2.yaxis.set_label_coords(-0.2, 0.25)#, 0.28)



        #if save==True:
            #plt.savefig(path_save, bbox_inches='tight',dpi=300)
            #plt.tight_layout() #pad en inch si besoin
            #plt.savefig(path_save,dpi=300)
        #    plt.savefig(path_save)









    plt.show()




if __name__=='__main__':
    trace_moment_cin_all(simu,tag,output_min,output_max,width_bar,hist,diff_abs,diff_relat,legend,marker,0.,moyenne_glissante,save)

