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

plt.style.use("pdf")
plt.style.use("aanda_modif")



if __name__=='__main__':
#-----------------------------------------------------------
#Noms des simulations et caracteristiques du calcul du rayon
#-----------------------------------------------------------
    #simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire_shr_bigbox_50pourc'
    tag=''
    legend='8,  10*40'
    marker='.'

    simu='B335_noturb_norot_hydro_pert_asym_aleatoire_shr_bigbox_50pourc'#+tag

    output_min_disque = 'None'
    output_min_global = 1
    output_max = 110
    output_frag = 80  #output au-dela duquel le disque fragmente : algorithme faux, 'None' si fragmente pas

    seuil_rho = 1e-10
    t1 = 0.

    width_bar = 1

    hist=True
    diff_abs=False
    diff_relat=False
    moyenne_glissante=False

    save=True



def trace_moment_cin_disque(simu,tag,output_min_disque,output_min_global,output_max,width_bar,hist,diff_abs,diff_relat,legend,marker,t1=0.,seuil_rho=1e-10,moyenne_glissante=True,save=False,output_frag='None'):
#------------------
#Differents chemins
#------------------
    path_analyses='/home/averliat/these/analyses/'+simu+'/'

    if output_min_disque=='None':
            ref=np.loadtxt(path_analyses+'t0_seuil_rho_'+str(seuil_rho)+'.txt')
            output_0 = int(ref[0])
            time_0 = ref[1]
            output_min_disque = output_0

    file_save = 'Moment_cin_observateurs_output'+str(output_min_global)+'_to_'+str(output_max)+'.hdf5'
    file_save2 = 'Moment_cin_observateurs_dans_disque_par_pdf_output'+str(output_min_disque)+'_to_'+str(output_max)+'.hdf5'
#file_save2 = 'Moment_cin_observateurs_dans_coquille_de_densite1e+11_output'+str(output_min)+'_to_'+str(output_max)+'.hdf5'

    if tag=='hr2':
        output_min_disque = output_min_disque -4  #Hack pour hr2


    name_save = simu+'_moment_dans_disque_vs_erreur_output_'+str(output_min_global)+'_'+str(output_max)+'.pdf'
    path_save = path_analyses+'Article_1/'+name_save


#------------------------------------
#Chargement des differentes quantites
#------------------------------------
    h5f = h5py.File(path_analyses+file_save, 'r')
    norm_mom_cin_obs=h5f['norm_mom_cin_obs'][output_min_disque-1:]
    simulation_time=h5f['simulation_time'][output_min_disque-1:]
    norm_GC_AU=h5f['norm_GC_AU'][output_min_disque-1:]
    norm_center_C_AU=h5f['norm_center_C_AU'][output_min_disque-1:]
    norm_center_G_AU=h5f['norm_center_G_AU'][output_min_disque-1:]
    vel_GC_tab=h5f['vel_GC_tab'][output_min_disque-1:]
    norm_mom_cin_obs_par_integ=h5f['norm_mom_cin_obs_par_integ'][output_min_disque-1:]
    numero_output=h5f['numero_output'][output_min_disque-1:]
    h5f.close()

    h5f = h5py.File(path_analyses+file_save2, 'r')
#norm_mom_cin_obs2=h5f['norm_mom_cin_obs'][:]
    simulation_time2=h5f['simulation_time'][:]
    norm_GC_AU2=h5f['norm_GC_AU'][:]
    norm_center_C_AU2=h5f['norm_center_C_AU'][:]
    norm_center_G_AU2=h5f['norm_center_G_AU'][:]
    vel_GC_tab2=h5f['vel_GC_tab'][:]
    norm_mom_cin_obs_par_integ2=h5f['norm_mom_cin_obs_par_integ'][:]
    numero_output2=h5f['numero_output'][:]
    h5f.close()



#------------------
#Debut de la figure
#------------------
    cmappts = plt.get_cmap('magma')
    colorspts = [cmappts(i) for i in np.linspace(0.1,0.9,7)]

#tab_x = np.linspace(2,102,51) #Utile pour hr2
#tab_x = np.linspace(output_min,output_max,output_max-output_min+1)
#tab_x = np.linspace(1,output_max,len(norm_mom_cin_obs_par_integ))

    if output_frag != 'None':
        frag=output_frag-output_min_disque+1
    else:
        frag=output_max-output_min_disque+1

    if hist==True:
        plt.figure()
        #mask=np.abs(norm_mom_cin_obs_par_integ-norm_mom_cin_obs)-norm_mom_cin_obs_par_integ2>0
        #mask_inv=np.logical_not(mask)

        #plt.bar(numero_output[mask], np.abs(norm_mom_cin_obs_par_integ[mask] - norm_mom_cin_obs[mask]), width=width_bar, color=colorspts[5], label=ur'| Moment cinétique réel - analytique | '+tag)

        #plt.bar(numero_output2[mask_inv], norm_mom_cin_obs_par_integ2[mask_inv], width=width_bar, color=colorspts[1], label=ur'Moment dans le disque '+tag)

        #plt.bar(numero_output[mask_inv], np.abs(norm_mom_cin_obs_par_integ[mask_inv] - norm_mom_cin_obs[mask_inv]), width=width_bar, color='w', edgecolor='black', linewidth=0.4)

        #plt.bar(numero_output2[mask], norm_mom_cin_obs_par_integ2[mask], width=width_bar, color='w', edgecolor='black', linewidth=0.4)



        plt.bar(numero_output2[0:frag], norm_mom_cin_obs_par_integ2[0:frag], width=width_bar, color=colorspts[1], label=ur'Momentum in the disk', edgecolor='black', linewidth=0.05)
        plt.bar(numero_output[0:frag], np.abs(norm_mom_cin_obs_par_integ - norm_mom_cin_obs)[0:frag], width=width_bar, color=colorspts[6], label=ur'$\left| \Delta \sigma  \right|$', edgecolor='black', linewidth=0.05)
        plt.xlabel(ur'Output')
        plt.ylabel(ur"Angular momentum ($kg.m^2.s^{-1}$)")
        plt.legend(loc='best')

        plt.tight_layout(pad=0.1) #pad en inch si besoin
        if save==True:
            plt.savefig(path_save)

        #plt.legend(loc='best')

        #plt.yscale('log')
        #plt.xlim((0,53))
        #plt.ylim((0,1.2e47))


    if diff_abs==True:
        plt.figure(20)
        plt.plot(simulation_time-t1, (norm_mom_cin_obs_par_integ2-np.abs((norm_mom_cin_obs_par_integ-norm_mom_cin_obs))), marker=marker, label=legend+'  - '+tag, linestyle='None')
        plt.xlabel(ur'Temps (Myr)')
        plt.ylabel(ur"Différence absolue entre erreur et moment disque($kg.m^2.s^{-1}$)")
        plt.legend(loc='best')

    if diff_relat==True:
        plt.figure(21)
        plt.plot(simulation_time-t1, (norm_mom_cin_obs_par_integ2-np.abs((norm_mom_cin_obs_par_integ-norm_mom_cin_obs)))/norm_mom_cin_obs_par_integ2, marker=marker, label=legend+'  - '+tag, linestyle='None')
        plt.xlabel(ur'Temps (Myr)')
        plt.ylabel(ur"Différence relative entre erreur et moment disque")
        plt.legend(loc='upper left')
        plt.ylim((0,1.3))







    if moyenne_glissante==True:
        erreur=(norm_mom_cin_obs_par_integ2-np.abs((norm_mom_cin_obs_par_integ-norm_mom_cin_obs)))/norm_mom_cin_obs_par_integ2
        erreur_moy=np.mean(np.nan_to_num(erreur[np.where(simulation_time-t1>0)]))
        plt.figure(22)
        plt.plot(tag, erreur_moy,marker=marker)
        plt.ylabel(ur"Moyenne de l'ecart entre l'erreur et le moment dans le disque")

        erreur=np.nan_to_num(erreur[np.where(simulation_time-t1>0)])
        time_pour_liss=np.nan_to_num((simulation_time-t1)[np.where(simulation_time-t1>0)])
        grp=10
        nbr_pts=len(erreur)/10
        erreur_lisse=np.zeros(nbr_pts)
        time_lisse=np.zeros(nbr_pts)
        for i in range(nbr_pts):
            erreur_lisse[i]=np.mean(erreur[i*grp:(i+1)*grp])
            time_lisse[i]=np.mean(time_pour_liss[i*grp:(i+1)*grp])
        plt.figure(23)
        plt.plot(time_lisse,erreur_lisse, marker=marker, label=legend+'  - '+tag)
        plt.xlabel(ur'Temps corrigé lissé (Myr)')
        plt.ylabel(ur"Différence relative entre erreur et moment disque, lissée")
        plt.legend(loc='upper left')
        plt.ylim((0,1.3))







    plt.show()


#plt.savefig('../Comparaisons/Verification_moment_physique/Moment_cinetique_analytique_et_reel_hr2.pdf', bbox_inches='tight')

#plt.savefig('../Comparaisons/Verification_moment_physique/Difference_absolue_moment_cinetique_analytique_et_reel_hr2.pdf', bbox_inches='tight')

#plt.savefig('../Comparaisons/Verification_moment_physique/Difference_relative_moment_cinetique_analytique_et_reel_hr2.pdf', bbox_inches='tight')

#plt.savefig('../Comparaisons/Verification_moment_physique/Difference_moment_cinetique_analytique_et_reel_vs_moment_dans_disque_50vhr_output'+str(output_min)+'_'+str(output_max)+'.pdf', bbox_inches='tight')
if __name__=='__main__':
        trace_moment_cin_disque(simu,tag,output_min_disque,output_min_global,output_max,width_bar,hist,diff_abs,diff_relat,legend,marker,t1,seuil_rho,moyenne_glissante,save,output_frag)
