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
    tag='_shr_bigbox_20pourc'
    legend='8,  10*40'  #Inutile pour le moment
    marker='.'

    simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire'+tag

    output = 40

    R_min = 50
    R_max = 10000
    dr = 50

    ang_absolu = True
    ang_relat = True
    omega = True
    moment_spec = True

    save_fig = False




def trace_angle_gradient_vitesse(simu,tag,output,R_min,R_max,dr,legend,marker,ang_absolu,ang_relat,omega,moment_spec,save_fig):
#------------------
#Differents chemins
#------------------
    path_analyses='/home/averliat/these/analyses/'+simu+'/'
    file_save_ang = 'Angle_omega_J_gradient_vitesse_log_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.hdf5'
    #file_save_ang = 'Angle_gradient_vitesse_output'+str(output)+'comparaison_mathilde.hdf5'



#------------------------------------
#Chargement des differentes quantites
#------------------------------------
    h5f = h5py.File(path_analyses+file_save_ang, 'r')

    shells_au=h5f['shells_au'][:]
    ang_x_tab=h5f['ang_x_tab'][:]
    ang_y_tab=h5f['ang_y_tab'][:]
    ang_z_tab=h5f['ang_z_tab'][:]
    if omega==True:
        omega_x_tab=h5f['omega_x_tab'][:]
        omega_y_tab=h5f['omega_y_tab'][:]
        omega_z_tab=h5f['omega_z_tab'][:]
    if moment_spec==True:
        moment_spec_x_tab=h5f['moment_spec_x_tab'][:]
        moment_spec_y_tab=h5f['moment_spec_y_tab'][:]
        moment_spec_z_tab=h5f['moment_spec_z_tab'][:]

    h5f.close()



#---------------------------------------------------
#Recalage des angles par rapport à l'angle initial
#---------------------------------------------------
    ang_x_tab[np.where(ang_x_tab<0)] = ang_x_tab[np.where(ang_x_tab<0)]+360
    ang_y_tab[np.where(ang_y_tab<0)] = ang_y_tab[np.where(ang_y_tab<0)]+360
    ang_z_tab[np.where(ang_z_tab<0)] = ang_z_tab[np.where(ang_z_tab<0)]+360

    
    seuil_ang=200

    ang_0=ang_x_tab[0]
    for i in range(len(ang_x_tab)):
        if np.abs(ang_0-ang_x_tab[i])>seuil_ang:
            ang_x_tab[i]=ang_x_tab[i]+360
            if np.abs(ang_0-ang_x_tab[i])>seuil_ang:
                ang_x_tab[i]=ang_x_tab[i]-720
        ang_0=ang_x_tab[i]

    ang_0=ang_y_tab[0]
    for i in range(len(ang_y_tab)):
        if np.abs(ang_0-ang_y_tab[i])>seuil_ang:
            ang_y_tab[i]=ang_y_tab[i]+360
            if np.abs(ang_0-ang_y_tab[i])>seuil_ang:
                ang_y_tab[i]=ang_y_tab[i]-720
        ang_0=ang_y_tab[i]

    ang_0=ang_z_tab[0]
    for i in range(len(ang_z_tab)):
        if np.abs(ang_0-ang_z_tab[i])>seuil_ang:
            ang_z_tab[i]=ang_z_tab[i]+360
            if np.abs(ang_0-ang_z_tab[i])>seuil_ang:
                ang_z_tab[i]=ang_z_tab[i]-720
        ang_0=ang_z_tab[i]
    

    ang_x_tab_recal = ang_x_tab - ang_x_tab[0]
    ang_y_tab_recal = ang_y_tab - ang_y_tab[0]
    ang_z_tab_recal = ang_z_tab - ang_z_tab[0]
    
    

#------------------
#Debut de la figure
#------------------
    cmappts = plt.get_cmap('magma')
    colorspts = [cmappts(i) for i in np.linspace(0.1,0.9,7)]

    if ang_absolu==True:
        plt.figure()
        plt.semilogx(shells_au, ang_x_tab, marker='.',color=colorspts[2], label='x,  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"Angle (°)")
        plt.legend(loc='best')

        #plt.figure()
        plt.semilogx(shells_au, ang_y_tab, marker='.',color=colorspts[4], label='y,  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"Angle (°)")
        plt.legend(loc='best')

        #plt.figure()
        plt.semilogx(shells_au, ang_z_tab, marker='.',color=colorspts[6], label='z,  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"Angle (°)")
        plt.legend(loc='best')

        if save_fig==True:
            plt.savefig('/home/averliat/these/analyses/Gradients_vitesse/Angle_gradient_vitesse_'+simu+'_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.pdf', bbox_inches='tight')


    if ang_relat==True:
        plt.figure()
        plt.semilogx(shells_au, ang_x_tab_recal, marker='.',color=colorspts[2], label='x,  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"Angle (°)")
        plt.legend(loc='best')

        #plt.figure()
        plt.semilogx(shells_au, ang_y_tab_recal, marker='.',color=colorspts[4], label='y,  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"Angle (°)")
        plt.legend(loc='best')

        #plt.figure()
        plt.semilogx(shells_au, ang_z_tab_recal, marker='.',color=colorspts[6], label='z,  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"Angle (°)")
        plt.legend(loc='best')

        if save_fig==True:
            plt.savefig('/home/averliat/these/analyses/Gradients_vitesse/Angle_relatif_gradient_vitesse_'+simu+'_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.pdf', bbox_inches='tight')



    if omega==True:
        plt.figure()
        plt.loglog(shells_au, omega_x_tab, marker='.',color=colorspts[2], label='x,  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"$\Omega$  $(km.s^{-1}.pc^{-1})")
        plt.legend(loc='best')

        #plt.figure()
        plt.loglog(shells_au, omega_y_tab, marker='.',color=colorspts[4], label='y,  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"$\Omega$  $(km.s^{-1}.pc^{-1})$")
        plt.legend(loc='best')

        #plt.figure()
        plt.loglog(shells_au, omega_z_tab, marker='.',color=colorspts[6], label='z,  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"$\Omega$  $(km.s^{-1}.pc^{-1})$")
        plt.legend(loc='best')

        if save_fig==True:
            plt.savefig('/home/averliat/these/analyses/Gradients_vitesse/Amplitude_gradient_vitesse_'+simu+'_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.pdf', bbox_inches='tight')



    if moment_spec==True:
        plt.figure()
        plt.loglog(shells_au, moment_spec_x_tab, marker='.',color=colorspts[2], label='x,  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"$j$  $(km.s^{-1}.pc)$")
        plt.legend(loc='best')

        #plt.figure()
        plt.loglog(shells_au, moment_spec_y_tab, marker='.',color=colorspts[4], label='y,  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"$j$  $(km.s^{-1}.pc)$")
        plt.legend(loc='best')

        #plt.figure()
        plt.loglog(shells_au, moment_spec_z_tab, marker='.',color=colorspts[6], label='z,  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"$j$  $(km.s^{-1}.pc)$")
        plt.legend(loc='best')

        if save_fig==True:
            plt.savefig('/home/averliat/these/analyses/Gradients_vitesse/Moment_specifique_gradient_vitesse_'+simu+'_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.pdf', bbox_inches='tight')









    #------------------------------------------------------------------------
    #Pour ajouter la droite de pente -1.5 sur le graph log(omega) = f(log(r))
    #------------------------------------------------------------------------
    y=-1.8*np.log10(shells_au)+6.8

    plt.figure(3)
    plt.loglog(shells_au,10**y,color='g')
    #------------------------------------------------------------------------







if __name__=='__main__':
    trace_angle_gradient_vitesse(simu,tag,output,R_min,R_max,dr,legend,marker,ang_absolu,ang_relat,omega,moment_spec,save_fig)

