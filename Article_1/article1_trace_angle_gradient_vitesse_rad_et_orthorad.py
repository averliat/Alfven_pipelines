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
    tag='_bigbox_50pourc_sink_seuil_haut_rot1'#'_shr_bigbox_50pourc'
    legend='8,  10*40'  #Inutile pour le moment
    marker='.'

    simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire'+tag

    output = 31


    view='x'

    R_min = 50
    R_max = 10000
    dr = 50

    ang_absolu = True
    ang_relat = True
    omega = True
    moment_spec = True

    save_fig = True

    reposition_fig=True



def trace_angle_gradient_vitesse(simu,tag,output,R_min,R_max,dr,legend,marker,ang_absolu,ang_relat,omega,moment_spec,save_fig):
#------------------
#Differents chemins
#------------------
    path='/drf/projets/alfven-data/averliat/'+simu+'/'
    path_analyses='/home/averliat/these/analyses/'+simu+'/'
    file_save_ang_rad_et_orthorad = 'Angle_omega_J_gradient_vitesse_rad_et_orthorad_log_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.hdf5'
    file_save_ang = 'Angle_omega_J_gradient_vitesse_log_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.hdf5'
    #file_save_ang = 'Angle_gradient_vitesse_output'+str(output)+'comparaison_mathilde.hdf5'




#------------------------------------
#Chargement des differentes quantites
#------------------------------------
    h5f = h5py.File(path_analyses+file_save_ang, 'r')

    shells_au=h5f['shells_au'][:]
    ang_tab=h5f['ang_'+view+'_tab'][:]
    if omega==True:
        omega_tab=h5f['omega_'+view+'_tab'][:]
    if moment_spec==True:
        moment_spec_tab=h5f['moment_spec_'+view+'_tab'][:]

    h5f.close()

    h5f = h5py.File(path_analyses+file_save_ang_rad_et_orthorad, 'r')

    shells_au2=h5f['shells_au'][:]
    if not (shells_au == shells_au2).all():
        print('==================================================')
        print('Attention, gradients pas calcules aux memes rayons')
        print('==================================================')
        return
    ang_rad_tab=h5f['ang_'+view+'_rad_tab'][:]
    ang_orthorad_tab=h5f['ang_'+view+'_orthorad_tab'][:]
    if omega==True:
        omega_rad_tab=h5f['omega_'+view+'_rad_tab'][:]
        omega_orthorad_tab=h5f['omega_'+view+'_orthorad_tab'][:]
    if moment_spec==True:
        moment_spec_rad_tab=h5f['moment_spec_'+view+'_rad_tab'][:]
        moment_spec_orthorad_tab=h5f['moment_spec_'+view+'_orthorad_tab'][:]

    h5f.close()



#---------------------------------------------------
#Recalage des angles par rapport à l'angle initial
#---------------------------------------------------
    ang_tab[np.where(ang_tab<0)] = ang_tab[np.where(ang_tab<0)]+360
    ang_rad_tab[np.where(ang_rad_tab<0)] = ang_rad_tab[np.where(ang_rad_tab<0)]+360
    ang_orthorad_tab[np.where(ang_orthorad_tab<0)] = ang_orthorad_tab[np.where(ang_orthorad_tab<0)]+360

    
    seuil_ang=190

    ang_0=ang_tab[0]
    for i in range(len(ang_tab)):
        if np.abs(ang_0-ang_tab[i])>seuil_ang:
            ang_tab[i]=ang_tab[i]+360
            if np.abs(ang_0-ang_tab[i])>seuil_ang:
                ang_tab[i]=ang_tab[i]-720
        ang_0=ang_tab[i]

    ang_0=ang_rad_tab[0]
    for i in range(len(ang_rad_tab)):
        if np.abs(ang_0-ang_rad_tab[i])>seuil_ang:
            ang_rad_tab[i]=ang_rad_tab[i]+360
            if np.abs(ang_0-ang_rad_tab[i])>seuil_ang:
                ang_rad_tab[i]=ang_rad_tab[i]-720
        ang_0=ang_rad_tab[i]

    ang_0=ang_orthorad_tab[0]
    for i in range(len(ang_orthorad_tab)):
        if np.abs(ang_0-ang_orthorad_tab[i])>seuil_ang:
            ang_orthorad_tab[i]=ang_orthorad_tab[i]+360
            if np.abs(ang_0-ang_orthorad_tab[i])>seuil_ang:
                ang_orthorad_tab[i]=ang_orthorad_tab[i]-720
        ang_0=ang_orthorad_tab[i]
    

    ang_tab_recal = ang_tab - ang_tab[0]
    ang_rad_tab_recal = ang_rad_tab - ang_rad_tab[0]
    ang_orthorad_tab_recal = ang_orthorad_tab - ang_orthorad_tab[0]
    

    seuil_ang_2 = 250

    if ang_tab[0] > seuil_ang_2:
        ang_tab -= 360
    elif ang_tab[0] < -seuil_ang_2:
        ang_tab += 360

    if ang_rad_tab[0] > seuil_ang_2:
        ang_rad_tab -= 360
    elif ang_rad_tab[0] < -seuil_ang_2:
        ang_rad_tab += 360

    if ang_orthorad_tab[0] > seuil_ang_2:
        ang_orthorad_tab -= 360
    elif ang_orthorad_tab[0] < -seuil_ang_2:
        ang_orthorad_tab += 360



#---------------------------------------------------------------------------------
#Calcul des t_0 de chaque simulation ou lecture des fichiers si calculs deja faits
#---------------------------------------------------------------------------------
    ro=pymses.RamsesOutput(path,output)
    factor_time_yr = ro.info['unit_time'].express(cst.year)
    simulation_time = ro.info['time']*factor_time_yr
    seuil_rho = 1e-10
    if os.path.isfile(path_analyses+'t0_seuil_rho_'+str(seuil_rho)+'.txt') == False:
        ref = np.array(t_0.temps_0_simu(path_t0, seuil_rho, sortie_output=1))
        np.savetxt(path_analyses+'t0_seuil_rho_'+str(seuil_rho)+'.txt', ref)
    else:
        ref=np.loadtxt(path_analyses+'t0_seuil_rho_'+str(seuil_rho)+'.txt')



#------------------
#Debut de la figure
#------------------
    cmappts = plt.get_cmap('cool')
    colorspts = [cmappts(i) for i in np.linspace(0.0,1,3)]
    colorspts = ['darkviolet','deepskyblue','magenta']

    if ang_absolu==True:
        plt.figure()
        plt.semilogx(shells_au, ang_tab, linestyle='dashdot',marker='None',color=colorspts[0], label='Full velocity')
        plt.xlabel(ur'Radius (AU)')
        plt.ylabel(ur"Angle (°)")
        plt.legend(loc='best')
        plt.title('Time = '+str(int(simulation_time - ref[1]*1e6))+' years  ~~~~~ $ \\varepsilon= 50 \%$ ~~~~~ $\\beta = 1 \%$')

        #plt.figure()
        plt.semilogx(shells_au, ang_rad_tab, marker='None',color=colorspts[1], linestyle='dashed',label='Radial component')
        plt.xlabel(ur'Radius (AU)')
        plt.ylabel(ur"Angle (°)")
        plt.legend(loc='best')

        #plt.figure()
        plt.semilogx(shells_au, ang_orthorad_tab, marker='None',color=colorspts[2], linestyle='dotted',label='Orthoradial component')
        plt.xlabel(ur'Radius (AU)')
        plt.ylabel(ur"Angle (°)")
        plt.legend(loc='best')


        if save_fig==True:
            plt.tight_layout(pad=0.1) #pad en inch si besoin
            plt.savefig(path_analyses+'Article_1/grad_rad_ortho/Angle_gradient_vitesse_'+simu+'_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.pdf')#, bbox_inches='tight')


    if ang_relat==True:
        plt.figure()
        plt.semilogx(shells_au, ang_tab_recal, linestyle='dashdot',marker='None',color=colorspts[0], label='Full velocity')
        plt.xlabel(ur'Radius (AU)')
        plt.ylabel(ur"Angle (°)")
        plt.legend(loc='best')
        plt.title('Time = '+str(int(simulation_time - ref[1]*1e6))+' years  ~~~~~ $ \\varepsilon= 50 \%$ ~~~~~ $\\beta = 1 \%$')

        #plt.figure()
        plt.semilogx(shells_au, ang_rad_tab_recal, marker='None',color=colorspts[1], linestyle='dashed',label='Radial component')
        plt.xlabel(ur'Radius (AU)')
        plt.ylabel(ur"Angle (°)")
        plt.legend(loc='best')

        #plt.figure()
        plt.semilogx(shells_au, ang_orthorad_tab_recal, marker='None',color=colorspts[2], linestyle='dotted',label='Orthoradial component')
        plt.xlabel(ur'Radius (AU)')
        plt.ylabel(ur"Angle (°)")
        plt.legend(loc='best')


        if save_fig==True:
            plt.tight_layout(pad=0.1) #pad en inch si besoin
            plt.savefig(path_analyses+'Article_1/grad_rad_ortho/Angle_relatif_gradient_vitesse_'+simu+'_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.pdf')#, bbox_inches='tight')



    if omega==True:
        plt.figure()
        plt.loglog(shells_au, omega_tab, linestyle='dashdot',marker='None',color=colorspts[0], label='Full velocity')
        plt.xlabel(ur'Radius (AU)')
        plt.ylabel(ur"$\Omega$  $(km.s^{-1}.pc^{-1})")
        plt.legend(loc='best')
        plt.title('Time = '+str(int(simulation_time - ref[1]*1e6))+' years  ~~~~~ $ \\varepsilon= 50 \%$ ~~~~~ $\\beta = 1 \%$')

        #plt.figure()
        plt.loglog(shells_au, omega_rad_tab, marker='None',color=colorspts[1], linestyle='dashed',label='Radial component')
        plt.xlabel(ur'Radius (AU)')
        plt.ylabel(ur"$\Omega$  $(km.s^{-1}.pc^{-1})$")
        plt.legend(loc='best')

        #plt.figure()
        plt.loglog(shells_au, omega_orthorad_tab, marker='None',color=colorspts[2], linestyle='dotted',label='Orthoradial component')
        plt.xlabel(ur'Radius (AU)')
        plt.ylabel(ur"$\Omega$  $(km.s^{-1}.pc^{-1})$")
        plt.legend(loc='best')


        if save_fig==True:
            plt.tight_layout(pad=0.1) #pad en inch si besoin
            plt.savefig(path_analyses+'Article_1/grad_rad_ortho/Amplitude_gradient_vitesse_'+simu+'_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.pdf')#, bbox_inches='tight')



    if moment_spec==True:
        plt.figure()
        plt.loglog(shells_au, moment_spec_tab, linestyle='dashdot',marker='None',color=colorspts[0], label='Full velocity')
        plt.xlabel(ur'Radius (AU)')
        plt.ylabel(ur"$j$  $(km.s^{-1}.pc)$")
        plt.legend(loc='best')
        plt.title('Time = '+str(int(simulation_time - ref[1]*1e6))+' years  ~~~~~ $ \\varepsilon= 50 \%$ ~~~~~ $\\beta = 1 \%$')

        #plt.figure()
        plt.loglog(shells_au, moment_spec_rad_tab, marker='None',color=colorspts[1], linestyle='dashed',label='Radial component')
        plt.xlabel(ur'Radius (AU)')
        plt.ylabel(ur"$j$  $(km.s^{-1}.pc)$")
        plt.legend(loc='best')

        #plt.figure()
        plt.loglog(shells_au, moment_spec_orthorad_tab, marker='None',color=colorspts[2], linestyle='dotted',label='Orthoradial component')
        plt.xlabel(ur'Radius (AU)')
        plt.ylabel(ur"$j$  $(km.s^{-1}.pc)$")
        plt.legend(loc='best')


        if save_fig==True:
            plt.tight_layout(pad=0.1) #pad en inch si besoin
            plt.savefig(path_analyses+'Article_1/grad_rad_ortho/Moment_specifique_gradient_vitesse_'+simu+'_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.pdf')#, bbox_inches='tight')









    #------------------------------------------------------------------------
    #Pour ajouter la droite de pente -1.5 sur le graph log(omega) = f(log(r))
    #------------------------------------------------------------------------
    #y=-1.8*np.log10(shells_au)+6.8

    #plt.figure(3)
    #plt.loglog(shells_au,10**y,color='g')
    #------------------------------------------------------------------------







if __name__=='__main__':
    trace_angle_gradient_vitesse(simu,tag,output,R_min,R_max,dr,legend,marker,ang_absolu,ang_relat,omega,moment_spec,save_fig)

    if reposition_fig==True:
        for i in range(2):
            i *= 2
            for j in range(2):
                plt.figure(i+j+1)
                mngr = plt.get_current_fig_manager()
                geom = mngr.window.geometry()
                x,y,dx,dy = geom.getRect()
                mngr.window.setGeometry(1920/3*(i/2),1080/2*j,dx,dy)

        plt.figure(1)  #Oblige sinon figure 1 mal placee...
        mngr = plt.get_current_fig_manager()
        geom = mngr.window.geometry()
        x,y,dx,dy = geom.getRect()
        #mngr.window.setGeometry(65,52,dx,dy)
        mngr.window.setGeometry(1,1,dx,dy)
