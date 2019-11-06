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
import matplotlib.patches as mpatches

import h5py


plt.style.use("pdf")
plt.style.use("aanda_modif")




if __name__=='__main__':
#-----------------------------------------------------------
#Noms des simulations et caracteristiques du calcul du rayon
#-----------------------------------------------------------
    tag=[10,20,30,40,50,60]
    outputs=[23,40,50,70,59,60]
    #tag=[10,20,30,50]
    #outputs=[23,40,50,59]


    legend='8,  10*40'  #Inutile pour le moment
    marker='.'


    view='z'

    R_min = 50
    R_max = 10000
    dr = 50

    ang_absolu = False
    ang_relat = False
    omega = False
    moment_spec = False

    save_fig = False

    reposition_fig=False



def trace_angle_gradient_vitesse(simu,tag,output,R_min,R_max,dr,legend,marker,ang_absolu,ang_relat,omega,moment_spec,save_fig,ind):
#------------------
#Differents chemins
#------------------
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


    if tag[ind]==20:
        ang_orthorad_tab+=360


#------------------
#Debut de la figure
#------------------
    cmappts = plt.get_cmap('plasma')
    #colorspts = [cmappts(i) for i in np.linspace(0.1,0.9,7)]
    colorspts = [cmappts(i) for i in np.linspace(0.,0.91,len(tag))]#0.1-0.82 pour 4 simus

    if ang_absolu==True:
        plt.figure(10)
        plt.semilogx(shells_au, ang_tab,linestyle='solid',color=colorspts[ind], label=view+',  '+str(tag[ind])+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"Angle (°)")
        plt.legend(loc='best')

        #plt.figure()
        plt.semilogx(shells_au, ang_rad_tab,linestyle='dashed',color=colorspts[ind], label='rad,  '+str(tag[ind])+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"Angle (°)")
        plt.legend(loc='best')

        #plt.figure()
        plt.semilogx(shells_au, ang_orthorad_tab,linestyle='dotted',color=colorspts[ind], label='orthorad,  '+str(tag[ind])+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"Angle (°)")
        plt.legend(loc='best')


        if save_fig==True:
            plt.savefig('/home/averliat/these/analyses/Gradients_vitesse_rad_et_orthorad/Angle_gradient_vitesse_'+simu+'_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.pdf', bbox_inches='tight')


    if ang_relat==True:
        plt.figure(11)
        plt.semilogx(shells_au, ang_tab_recal, marker='.',color=colorspts[ind], label=view+',  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"Angle (°)")
        plt.legend(loc='best')

        #plt.figure()
        plt.semilogx(shells_au, ang_rad_tab_recal, marker='.',color=colorspts[ind], label='rad,  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"Angle (°)")
        plt.legend(loc='best')

        #plt.figure()
        plt.semilogx(shells_au, ang_orthorad_tab_recal, marker='.',color=colorspts[ind], label='orthorad,  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"Angle (°)")
        plt.legend(loc='best')


        if save_fig==True:
            plt.savefig('/home/averliat/these/analyses/Gradients_vitesse_rad_et_orthorad/Angle_relatif_gradient_vitesse_'+simu+'_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.pdf', bbox_inches='tight')



    if omega==True:
        plt.figure(12)
        plt.loglog(shells_au, omega_tab, marker='.',color=colorspts[ind], label=view+',  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"$\Omega$  $(km.s^{-1}.pc^{-1})")
        plt.legend(loc='best')

        #plt.figure()
        plt.loglog(shells_au, omega_rad_tab, marker='.',color=colorspts[ind], label='rad,  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"$\Omega$  $(km.s^{-1}.pc^{-1})$")
        plt.legend(loc='best')

        #plt.figure()
        plt.loglog(shells_au, omega_orthorad_tab, marker='.',color=colorspts[ind], label='orthorad,  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"$\Omega$  $(km.s^{-1}.pc^{-1})$")
        plt.legend(loc='best')


        if save_fig==True:
            plt.savefig('/home/averliat/these/analyses/Gradients_vitesse_rad_et_orthorad/Amplitude_gradient_vitesse_'+simu+'_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.pdf', bbox_inches='tight')



    if moment_spec==True:
        plt.figure(13)
        plt.loglog(shells_au, moment_spec_tab, marker='.',color=colorspts[ind], label=view+',  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"$j$  $(km.s^{-1}.pc)$")
        plt.legend(loc='best')

        #plt.figure()
        plt.loglog(shells_au, moment_spec_rad_tab, marker='.',color=colorspts[ind], label='rad,  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"$j$  $(km.s^{-1}.pc)$")
        plt.legend(loc='best')

        #plt.figure()
        plt.loglog(shells_au, moment_spec_orthorad_tab, marker='.',color=colorspts[ind], label='orthorad,  '+tag+',  output '+str(output))
        plt.xlabel(ur'Rayon ($AU$)')
        plt.ylabel(ur"$j$  $(km.s^{-1}.pc)$")
        plt.legend(loc='best')


        if save_fig==True:
            plt.savefig('/home/averliat/these/analyses/Gradients_vitesse_rad_et_orthorad/Moment_specifique_gradient_vitesse_'+simu+'_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.pdf', bbox_inches='tight')




    fig=plt.figure(17)
    ax=plt.gca()
    ax.semilogx(shells_au, ang_rad_tab-ang_tab,linestyle='solid',color=colorspts[ind], label=ur'$\varepsilon = $'+str(tag[ind])+ur'$\%$')#, label=ur'$\theta_{\text{rad}} - \theta~~~~$ $\epsilon = $'+str(tag[ind])+ur'$\%$')
    plt.xlabel(ur'Radius (AU)')
    plt.ylabel(ur"Angle ($^\circ$)")
    #plt.legend(loc='best')

    #plt.figure()
    ax.semilogx(shells_au, ang_orthorad_tab-ang_tab,linestyle='dashed',color=colorspts[ind])#, label=ur'$\theta_{\text{orthorad}} - \theta~~~~$ $\epsilon = $'+str(tag[ind])+ur'$\%$')
    plt.xlabel(ur'Radius (AU)')
    plt.ylabel(ur"Angle ($^\circ$)")
    #plt.legend(loc='best')


    if ind==len(tag)-1:
        #Get artists and labels for legend and chose which ones to display
        handles, labels = ax.get_legend_handles_labels()
        display = range(0,len(tag))#(0,1,2,3)

        #Create custom artists
        simArtist = plt.Line2D((0,1),(0,0), color='dimgrey', linestyle='solid')
        anyArtist = plt.Line2D((0,1),(0,0), color='dimgrey',linestyle='dashed')

        #Create custom artists for color
        handles=[]
        for j in range(len(tag)):
            handles.append(mpatches.Patch(color=colorspts[j]))#, label='The red data')




        #Create legend from custom artist/label lists
        #ax.legend([simArtist,anyArtist]+[handle for i,handle in enumerate(handles) if i in display],
        #        [ur'$\theta_{\text{rad}} - \theta$', ur'$\theta_{\text{orthorad}} - \theta$']+[label for i,label in enumerate(labels) if i in display])
        color_legend=ax.legend([handle for i,handle in enumerate(handles) if i in display],
                [label for i,label in enumerate(labels) if i in display],loc='best',ncol=2)

        ax.add_artist(color_legend)

        ax.legend([simArtist,anyArtist],
                #[ur'$\theta_{\text{rad}} - \theta$', ur'$\theta_{\text{orthorad}} - \theta$'],loc='upper right')
                [ur'$\Delta \theta_{\text{rad}}$', ur'$\Delta \theta_{\text{orthorad}}$'],loc='upper right')


    #------------------------------------------------------------------------
    #Pour ajouter la droite de pente -1.5 sur le graph log(omega) = f(log(r))
    #------------------------------------------------------------------------
    #y=-1.8*np.log10(shells_au)+6.8

    #plt.figure(3)
    #plt.loglog(shells_au,10**y,color='g')
    #------------------------------------------------------------------------







if __name__=='__main__':

    for ind in range(len(tag)):
        simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire_shr_bigbox_'+str(tag[ind])+'pourc'
        output=outputs[ind]


        trace_angle_gradient_vitesse(simu,tag,output,R_min,R_max,dr,legend,marker,ang_absolu,ang_relat,omega,moment_spec,save_fig,ind)










    #plt.tight_layout(pad=0.1) #pad en inch si besoin
    #plt.savefig('/home/averliat/these/analyses/article1_figures/Angles_gradient_rad_orthorad_plusieurs_allsimus.pdf')


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
