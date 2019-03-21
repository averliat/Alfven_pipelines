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

    output = 685

    R_min = 50
    R_max = 10000
    dr = 50


    V_rms = True
    superpose_omega = True

    V_abs_moy = False
    V_abs_moy_fit = False

    save_fig = True
    save_with_fit = False




def trace_Vrms_gradient_vitesse(simu,tag,output,R_min,R_max,dr,legend,marker,V_rms,V_abs_moy,V_abs_moy_fit,superpose_omega,save_fig,save_with_fit):
#------------------
#Differents chemins
#------------------
    path_analyses='/home/averliat/these/analyses/'+simu+'/'
    file_save_vrms = 'Vrms_gradient_vitesse_log_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.hdf5'
    file_save_ang = 'Angle_omega_J_gradient_vitesse_log_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.hdf5'
    #file_save_ang = 'Angle_gradient_vitesse_output'+str(output)+'comparaison_mathilde.hdf5'



#------------------------------------
#Chargement des differentes quantites
#------------------------------------
    h5f = h5py.File(path_analyses+file_save_vrms, 'r')

    shells_au=h5f['shells_au'][:]
    grad_v=h5f['grad_v'][:]
    Vrms=h5f['Vrms'][:]
    v_abs_moy=h5f['v_abs_moy'][:]

    h5f.close()


    if superpose_omega==True:
        h5f = h5py.File(path_analyses+file_save_ang, 'r')

        omega_x_tab=h5f['omega_x_tab'][:]
        omega_y_tab=h5f['omega_y_tab'][:]
        omega_z_tab=h5f['omega_z_tab'][:]

        h5f.close()



#------------------
#Debut de la figure
#------------------
    cmappts = plt.get_cmap('magma')
    colorspts = [cmappts(i) for i in np.linspace(0.1,0.9,7)]


    if V_rms==True:
        #plt.figure()
        #plt.plot(shells_au,Vrms, marker='.',color=colorspts[2])
        #plt.xlabel(ur'Rayon ($AU$)')
        #plt.ylabel(ur"$V_{RMS}$  $(km.s^{-1})$")
        #plt.legend(loc='best')

        plt.figure()
        plt.loglog(shells_au,Vrms, marker='.',color=colorspts[2])
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"$V_{RMS}$  $(km.s^{-1})$")
        #plt.legend(loc='best')

        #plt.figure()
        #plt.loglog(shells_au, grad_v, marker='.',color=colorspts[4])
        #plt.xlabel(ur'Rayon ($AU$)')
        #plt.ylabel(ur"$\Omega_{estimation}=\frac{V_{RMS}}{R}$  $(km.s^{-1}.pc^{-1})$")
        #plt.legend(loc='best')

        if save_fig==True:
            plt.savefig('/home/averliat/these/analyses/Gradients_vitesse/Analyses_champ_vitesse/V_RMS/V_RMS_'+simu+'_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.pdf', bbox_inches='tight')



    if superpose_omega==True:
        plt.figure()
        plt.loglog(shells_au, omega_x_tab, marker='.',color=colorspts[2], label='x,  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"$\Omega$  $(km.s^{-1}.pc^{-1})")
        plt.legend(loc='best')

        plt.loglog(shells_au, omega_y_tab, marker='.',color=colorspts[4], label='y,  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"$\Omega$  $(km.s^{-1}.pc^{-1})$")
        plt.legend(loc='best')

        plt.loglog(shells_au, omega_z_tab, marker='.',color=colorspts[6], label='z,  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"$\Omega$  $(km.s^{-1}.pc^{-1})$")
        plt.legend(loc='best')

        plt.loglog(shells_au, grad_v, marker='.',color='g', label=ur'Estimation avec $\frac{V_{RMS}}{R}$,  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"$\Omega$  $(km.s^{-1}.pc^{-1})$")
        plt.legend(loc='best')

        if save_fig==True:
            plt.savefig('/home/averliat/these/analyses/Gradients_vitesse/Analyses_champ_vitesse/V_RMS/Estimation_amplitude_gradient_vitesse_'+simu+'_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.pdf', bbox_inches='tight')



    if V_abs_moy==True:
        plt.figure(20)
        plt.plot(shells_au,v_abs_moy, marker='.',color=colorspts[2])
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"$<|V|>$  $(km.s^{-1})$")
        #plt.legend(loc='best')

        plt.figure(21)
        plt.loglog(shells_au,v_abs_moy, marker='.',color=colorspts[2])
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"$<|V|>$  $(km.s^{-1})$")
        #plt.legend(loc='best')

        if V_abs_moy_fit==True:
            #-----------------------
            #1er fit du diagramme PV
            #-----------------------
            def func_2(x,b):
                y=-0.5*x+b
                return y

            abs_min=100  #np.min(shells_au) #AU
            abs_max=1000  #np.max(shells_au) #AU

            abs_fit=np.log10(shells_au[np.min(np.where(shells_au>abs_min)):np.max(np.where(shells_au<abs_max))])
            ord_fit=np.log10(v_abs_moy[np.min(np.where(shells_au>abs_min)):np.max(np.where(shells_au<abs_max))])

            #Ajustement
            p02= 0.
            popt,pcov=opt.curve_fit(func_2,abs_fit,ord_fit,p0=p02) #On fixe r**0.5

            G_grav = 6.67408e-11  #m3.kg-1.s-2
            Msun = 1.9884e30  #kg
            AU = 149597870700  #m

            M_kepler = 1e6*AU/(G_grav*Msun) *10**(2.*popt[0])
            print('Masse keplerienne = '+str(M_kepler)+' Msun')

            y1=-0.5*np.log10(shells_au) + popt[0]
            y1=10**y1
            x1=shells_au

            plt.figure(20)
            plt.plot(x1,y1,color='g',label=ur'Fit ['+str(abs_min)+','+str(abs_max)+ur']AU, $v=\left( \frac{G.M}{R} \right) ^{-0.5}$ : $M='+str(int(100*M_kepler)/100.)+'M_\odot$')

            plt.figure(21)
            plt.loglog(x1,y1,color='g',label=ur'Fit ['+str(abs_min)+','+str(abs_max)+ur']AU, $v=\left( \frac{G.M}{R} \right) ^{-0.5}$ : $M='+str(int(100*M_kepler)/100.)+'M_\odot$')



            #------------------------
            #2eme fit du diagramme PV
            #------------------------
            def func_fit(x,a,b):
                y=a*x+b
                return y

            abs_min=100  #np.min(shells_au) #AU
            abs_max=1000  #np.max(shells_au) #AU

            abs_fit=np.log10(shells_au[np.min(np.where(shells_au>abs_min)):np.max(np.where(shells_au<abs_max))])
            ord_fit=np.log10(v_abs_moy[np.min(np.where(shells_au>abs_min)):np.max(np.where(shells_au<abs_max))])

            #Ajustement
            p0= 0., 1.
            popt,pcov=opt.curve_fit(func_fit,abs_fit,ord_fit,p0=p0)
            print(popt, pcov)

            G_grav = 6.67408e-11  #m3.kg-1.s-2
            Msun = 1.9884e30  #kg
            AU = 149597870700  #m

            M_kepler = 1e6*AU/(G_grav*Msun) *10**(2.*popt[1])
            print('Masse keplerienne = '+str(M_kepler)+' Msun')


            #y1=popt[0]*abs_fit + popt[1]
            y1=popt[0]*np.log10(shells_au) + popt[1]
            y1=10**y1
            #x1=10**abs_fit
            x1=shells_au

            plt.figure(20)
            plt.plot(x1,y1,color='b',label=ur'Fit ['+str(abs_min)+','+str(abs_max)+ur']AU, $v=\left( \frac{G.M}{R} \right) ^{\alpha}$ : $\alpha ='+str(int(1000*popt[0])/1000.)+'$ ,   $M='+str(int(100*M_kepler)/100.)+'M_\odot$')
            plt.legend(loc='best')

            plt.figure(21)
            plt.loglog(x1,y1,color='b',label=ur'Fit ['+str(abs_min)+','+str(abs_max)+ur']AU, $v=\left( \frac{G.M}{R} \right) ^{\alpha}$ : $\alpha ='+str(int(1000*popt[0])/1000.)+'$,  $M='+str(int(100*M_kepler)/100.)+'M_\odot$')
            plt.legend(loc='best')

            #----------
            #Sauvegarde
            #----------
            if save_with_fit==True:
                plt.figure(20)
                plt.savefig('/home/averliat/these/analyses/Gradients_vitesse/Analyses_champ_vitesse/Moy_V_abs_et_fit_'+simu+'_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.pdf', bbox_inches='tight')
                plt.figure(21)
                plt.savefig('/home/averliat/these/analyses/Gradients_vitesse/Analyses_champ_vitesse/Moy_V_abs_et_fit_log_'+simu+'_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.pdf', bbox_inches='tight')



if __name__=='__main__':
    trace_Vrms_gradient_vitesse(simu,tag,output,R_min,R_max,dr,legend,marker,V_rms,V_abs_moy,V_abs_moy_fit,superpose_omega,save_fig,save_with_fit)

