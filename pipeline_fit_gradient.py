import numpy as np

import pymses
import sys
import glob as glob
import os

from pymses.filters import CellsToPoints
from pymses.utils import constants as cst
from pymses.analysis import Camera, raytracing, slicing
from pymses.analysis import ScalarOperator, FractionOperator, MaxLevelOperator
import matplotlib.pyplot as plt
plt.ion()
from matplotlib.colors import Normalize

#from astropy.io import fits

import os
import pipeline_temps_0_simulation as t_0

import scipy.optimize as opt




#-------------------------------------------------------
#Entree le nom de la simulation et le numero de l'output
#-------------------------------------------------------
simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire50_vhr'

owner = 'sapanais'#'averliat_alfven'
num_output = 254

save = False
dir_save = 'Coupe_vitesse_integree'

radius_zoom = 3

v_proj = True

title_time=True
title_time_cor=False
seuil_rho = 1e-10

fleche_vel = False
nbre_fleche = 20  #En fait plus c'est faible plus y'a de fleches...


selon_x = True
selon_y = True
selon_z = True

vmin_vel = None#-0.25
vmax_vel = None#-vmin_vel

vmin_dens = None
vmax_dens = None


reposition_fig = False #Pour repositionner les figures ouvertes par matplotlib


verif_fit = True  #Pour verifier si le fit du gradient de vitesse est juste

#Pour fermer toutes les figures avec matplotlib : plt.close('all')



#--------------------------------------------------------------
#Chemin de la simulation et du dossier de sauvegarde des images
#--------------------------------------------------------------
if owner=='averliat':
    path='/mnt/magmist/magmist/simu_B335_averliat/'+simu+'/'
if owner=='phennebe':
    path='/drf/projets/capucine/'+owner+'/'+simu+'/'
if owner=='sapanais':
    path='/dsm/anais/storageA/magmist/'+simu+'/'
if owner=='averliat_alfven':
    path='/drf/projets/alfven-data/averliat/'+simu+'/'
    
path_save='/home/averliat/these/analyses/'+simu+'/'+dir_save+'/'
path_analyse='/home/averliat/these/analyses/'+simu+'/'

#if simu == 'B335_noturb_norot_hydro_pert_asym_aleatoire50_vhr':
#   path_t0='/mnt/magmist/magmist/simu_B335_averliat/'+simu+'/'
#else:
#   path_t0=path
path_t0=path

if save==True:
    if os.path.isdir(path_save) == False:
        os.mkdir(path_save)



#-------------------
#Lecture de l'output
#-------------------
ro=pymses.RamsesOutput(path,num_output)
lbox_pc = ro.info['unit_length'].express(cst.pc)

amr = ro.amr_source(["rho","vel","P","phi","g"])

cell_source = CellsToPoints(amr)
cells = cell_source.flatten()

pos = cells.points
rho = cells["rho"]



#------------------------------------------------------------
#Facteur de conversion des unites de code en unites physiques
#------------------------------------------------------------
lbox=ro.info['boxlen']
lbox_m = ro.info['unit_length'].express(cst.m)
lbox_au = ro.info['unit_length'].express(cst.au)
lbox_cm = ro.info['unit_length'].express(cst.cm)
factor_time_yr = ro.info['unit_time'].express(cst.year)
factor_vel_km_s = ro.info['unit_velocity'].express(cst.km/cst.s)

simulation_time = ro.info['time']*factor_time_yr



#---------------------------------------------------------------------------------
#Calcul des t_0 de chaque simulation ou lecture des fichiers si calculs deja faits
#---------------------------------------------------------------------------------
if title_time_cor == True:
    if os.path.isfile(path_analyse+'t0_seuil_rho_'+str(seuil_rho)+'.txt') == False:
        ref = np.array(t_0.temps_0_simu(path_t0, seuil_rho, sortie_output=1))
        np.savetxt(path_analyse+'t0_seuil_rho_'+str(seuil_rho)+'.txt', ref)
    else:
        ref=np.loadtxt(path_analyse+'t0_seuil_rho_'+str(seuil_rho)+'.txt')


    #Erreur si on cherche a etudier un fichier inferieur au t0 de la simulation de reference
    if num_output < ref[0]:
        print
        print
        print("=================================================")
        print("=================================================")
        print("/!\   output_ref  <  output_ref_t0   /!\ ")
        print("=================================================")
        print("=================================================")
        title_time_cor=0



#---------------------------------------------------------
#Definition du centre des images et de leur niveau de zoom
#---------------------------------------------------------
#Position du "centre" de la simulation = endroit le plus dense
if radius_zoom==5:
    center = [0.5,0.5,0.5]
else:
    arg_centre = np.argmax(rho)
    center = [pos[:,0][arg_centre],pos[:,1][arg_centre],pos[:,2][arg_centre]]
    

zoom_v=[0.045, 0.015, 0.005, 0.005/3., lbox_pc*4]
radius=zoom_v[radius_zoom-1]
#radius=float(zoom_v[np.where(zoom_v==radius_zoom)])     #0.015#0.005 #Niveau de zoom correspondant au niveau '3' des images de "pipeline_image_unique.py"


#---------------------------------------------------------------------------------------------------------------
#Definition prise sur https://stackoverflow.com/questions/20144529/shifted-colorbar-matplotlib/20146989#20146989 
#pour avoir le centre de la colorbar a 0
#---------------------------------------------------------------------------------------------------------------
class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))




#=============================================================================
#----------------------------------------
#Fonction de fit des gradients de vitesse
#----------------------------------------
def func_vel((x,y),v0,a,b):
    v=v0 + a*x + b*y
    return v.ravel()


def fit_gradient_vel((map_V,extent),verif_fit=verif_fit):
    #Initialisation des differents tableaux
    tab_abs=np.linspace(extent[0],extent[1],len(map_V)) - center[0]*lbox_au
    tab_ord=np.linspace(extent[2],extent[3],len(map_V)) - center[1]*lbox_au
    tab_abs,tab_ord=np.meshgrid(tab_abs,tab_ord)

    #Valeurs approchees des parametres a ajuster
    p0= 0., 1., 1.

    #Ajustement
    popt,pcov=opt.curve_fit(func_vel,(tab_abs,tab_ord),map_V.ravel(),p0=p0)
    v0=popt[0]
    a=popt[1]
    b=popt[2]

    theta=np.arctan(b/a)  #orientation de la ligne de vitesse nulle
    vel=v0+a*tab_abs+b*tab_ord



    #Ajout de la fleche representant le gradient de vitesse
    x_fleche=np.linspace((center[0]-radius)*lbox_au, (center[0]+radius)*lbox_au,100)
    C=(center[1]-np.tan(theta)*center[0])*lbox_au
    y_fleche=np.tan(theta)*x_fleche+C

    taille_fleche_x=np.sqrt( 16*(radius*lbox_au)**2/25. /(1+np.tan(theta)**2) )  #Fleche comprise entre 1/10eme et 9/10eme de la boite
    mask_fleche=(x_fleche>center[0]*lbox_au-taille_fleche_x) & (x_fleche<center[0]*lbox_au+taille_fleche_x)



    #Necessaire pour determiner si la fleche va de droite a gauche ou de gauche a droite
    tab_abs=np.linspace(extent[0],extent[1],len(map_V))
    tab_ord=np.linspace(extent[2],extent[3],len(map_V))

    ind_abs_max=np.where(tab_abs==tab_abs.flat[np.abs(tab_abs - x_fleche[mask_fleche][-1]).argmin()])
    ind_ord_max=np.where(tab_ord==tab_ord.flat[np.abs(tab_ord - y_fleche[mask_fleche][-1]).argmin()])

    ind_abs_min=np.where(tab_abs==tab_abs.flat[np.abs(tab_abs - x_fleche[mask_fleche][0]).argmin()])
    ind_ord_min=np.where(tab_ord==tab_ord.flat[np.abs(tab_ord - y_fleche[mask_fleche][0]).argmin()])

    vel_droite = vel[ind_ord_max,ind_abs_max]
    vel_gauche = vel[ind_ord_min,ind_abs_min]


    if vel_droite > vel_gauche:
        plt.arrow(x_fleche[mask_fleche][0],y_fleche[mask_fleche][0],x_fleche[mask_fleche][-1]-x_fleche[mask_fleche][0],y_fleche[mask_fleche][-1]-y_fleche[mask_fleche][0],head_width=10, head_length=10,linestyle="--",linewidth=1.1,color='g',length_includes_head=True,joinstyle='round',overhang=0.7)
    elif vel_droite < vel_gauche:
        plt.arrow(x_fleche[mask_fleche][-1],y_fleche[mask_fleche][-1],x_fleche[mask_fleche][0]-x_fleche[mask_fleche][-1],y_fleche[mask_fleche][0]-y_fleche[mask_fleche][-1],head_width=10, head_length=10,linestyle="--",linewidth=1.1,color='g',length_includes_head=True,joinstyle='round',overhang=0.7)



    #Trace de la carte du gradient de vitesse ajuste pour verification
    if verif_fit==True:
        plt.figure()
        plt.imshow(vel,extent=extent,origin='lower',cmap='RdBu_r')
        if vel_droite > vel_gauche:
            plt.arrow(x_fleche[mask_fleche][0],y_fleche[mask_fleche][0],x_fleche[mask_fleche][-1]-x_fleche[mask_fleche][0],y_fleche[mask_fleche][-1]-y_fleche[mask_fleche][0],head_width=10, head_length=10,linestyle="--",linewidth=1.1,color='g',length_includes_head=True,joinstyle='round',overhang=0.7)
        elif vel_droite < vel_gauche:
            plt.arrow(x_fleche[mask_fleche][-1],y_fleche[mask_fleche][-1],x_fleche[mask_fleche][0]-x_fleche[mask_fleche][-1],y_fleche[mask_fleche][0]-y_fleche[mask_fleche][-1],head_width=10, head_length=10,linestyle="--",linewidth=1.1,color='g',length_includes_head=True,joinstyle='round',overhang=0.7)

        plt.xlabel('(AU)')     
        plt.ylabel('(AU)')
        cbar = plt.colorbar()         
        cbar.set_label(r'$v \, (km.s^{-1})$')  
        if title_time==True:
            plt.title('Time = '+str(int(simulation_time))+' years')
        if title_time_cor==True:
            plt.title('Time = '+str(int(simulation_time))+' years \n Corrected time = '+str(int(simulation_time - ref[1]*1e6))+' years')

    return
#=============================================================================














#-----------------
#-----------------
#Calcul des cartes
#-----------------
#-----------------

if selon_x==True:
    #--------------------------------------------
    #Calcul de la carte ou l'on regarde suivant x
    #--------------------------------------------
    cam_x = Camera(center=center,line_of_sight_axis='x',region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector='z',map_max_size=512)

    rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_density"])
    rt = raytracing.RayTracer(amr,ro.info,rho_op)
    datamap = rt.process(cam_x, surf_qty=True)
    map_col = np.log10(datamap.map.T*lbox_cm)

    if v_proj == True:
        Vx_op = ScalarOperator(lambda dset: dset["vel"][...,0]*dset["rho"] ,  ro.info["unit_velocity"])
        rt = raytracing.RayTracer(amr,ro.info,Vx_op)
        datamap_vx = rt.process(cam_x, surf_qty=True)
        map_Vx = datamap_vx.map.T / datamap.map.T * factor_vel_km_s



        plt.figure()
        plt.imshow(map_col,extent=[(-radius+center[0])*lbox_au,(radius+center[0])*lbox_au,(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au],origin='lower', vmin=vmin_dens, vmax=vmax_dens)   
        plt.xlabel('$y$ (AU)')     
        plt.ylabel('$z$ (AU)')
        cbar = plt.colorbar()
        if title_time==True:
            plt.title('Time = '+str(int(simulation_time))+' years')
        if title_time_cor==True:
            plt.title('Time = '+str(int(simulation_time))+' years \n Corrected time = '+str(int(simulation_time - ref[1]*1e6))+' years')
        cbar.set_label(r'$log(N) \, \, (cm^{-2})$')
        if radius_zoom==5:
            plt.xlim([0,lbox_au])
            plt.ylim([0,lbox_au])
        if save==True:
            plt.savefig(path_save+'dens_x_'+str(radius_zoom)+'_'+str(num_output)+'.pdf', bbox_inches='tight')


    if v_proj == True:
        plt.figure()
        norm = MidpointNormalize(midpoint=0)  #Pour avoir le centre de la colormap a 0
        extent=[(-radius+center[0])*lbox_au,(radius+center[0])*lbox_au,(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au]
        plt.imshow(map_Vx,extent=[(-radius+center[0])*lbox_au,(radius+center[0])*lbox_au,(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au],origin='lower',cmap='RdBu_r',norm=norm, vmin=vmin_vel, vmax=vmax_vel)
        plt.xlabel('$y$ (AU)')     
        plt.ylabel('$z$ (AU)')
        cbar = plt.colorbar()         
        cbar.set_label(r'$v_x \, (km.s^{-1})$')  
        if title_time==True:
            plt.title('Time = '+str(int(simulation_time))+' years')
        if title_time_cor==True:
            plt.title('Time = '+str(int(simulation_time))+' years \n Corrected time = '+str(int(simulation_time - ref[1]*1e6))+' years')
        if radius_zoom==5:
            plt.xlim([0,lbox_au])
            plt.ylim([0,lbox_au])
        if save==True:
            plt.savefig(path_save+'vel_x_'+str(radius_zoom)+'_'+str(num_output)+'.pdf', bbox_inches='tight')


        #Fit du gradient de vitesse
        fit_gradient_vel((map_Vx,extent))





if selon_y==True:
    #--------------------------------------------
    #Calcul de la carte ou l'on regarde suivant y
    #--------------------------------------------
    cam_y = Camera(center=center,line_of_sight_axis='y',region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector='x',map_max_size=512)

    rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_density"])
    rt = raytracing.RayTracer(amr,ro.info,rho_op)
    datamap = rt.process(cam_y, surf_qty=True)
    map_col = np.log10(datamap.map.T*lbox_cm)

    if v_proj == True:
        Vy_op = ScalarOperator(lambda dset: dset["vel"][...,1]*dset["rho"] ,  ro.info["unit_velocity"])
        rt = raytracing.RayTracer(amr,ro.info,Vy_op)
        datamap_vy = rt.process(cam_y, surf_qty=True)
        map_Vy = datamap_vy.map.T / datamap.map.T * factor_vel_km_s



        plt.figure()
        im = plt.imshow(map_col,extent=[(-radius+center[0])*lbox_au,(radius+center[0])*lbox_au,(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au],origin='lower', vmin=vmin_dens, vmax=vmax_dens)   
        plt.xlabel('$z$ (AU)')     
        plt.ylabel('$x$ (AU)')
        cbar=plt.colorbar()                                                                                   
        cbar.set_label(r'$log(N) \, \, (cm^{-2})$')
        if title_time==True:
            plt.title('Time = '+str(int(simulation_time))+' years')
        if title_time_cor==True:
            plt.title('Time = '+str(int(simulation_time))+' years \n Corrected time = '+str(int(simulation_time - ref[1]*1e6))+' years')
        if radius_zoom==5:
            plt.xlim([0,lbox_au])
            plt.ylim([0,lbox_au])
        if save==True:
            plt.savefig(path_save+'dens_y_'+str(radius_zoom)+'_'+str(num_output)+'.pdf', bbox_inches='tight') 


    if v_proj == True:
        plt.figure()
        norm = MidpointNormalize(midpoint=0)  #Pour avoir le centre de la colormap a 0
        extent=[(-radius+center[0])*lbox_au,(radius+center[0])*lbox_au,(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au]
        plt.imshow(map_Vy,extent=[(-radius+center[0])*lbox_au,(radius+center[0])*lbox_au,(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au],origin='lower',cmap='RdBu_r',norm=norm, vmin=vmin_vel, vmax=vmax_vel)
        plt.xlabel('$z$ (AU)')     
        plt.ylabel('$x$ (AU)')
        cbar = plt.colorbar()          
        cbar.set_label(r'$v_y \, (km.s^{-1})$')
        if title_time==True:
            plt.title('Time = '+str(int(simulation_time))+' years')
        if title_time_cor==True:
            plt.title('Time = '+str(int(simulation_time))+' years \n Corrected time = '+str(int(simulation_time - ref[1]*1e6))+' years')  
        if radius_zoom==5:
            plt.xlim([0,lbox_au])
            plt.ylim([0,lbox_au])
        if save==True:
            plt.savefig(path_save+'vel_y_'+str(radius_zoom)+'_'+str(num_output)+'.pdf', bbox_inches='tight')


        #Fit du gradient de vitesse
        fit_gradient_vel((map_Vy,extent))





if selon_z==True:
    #--------------------------------------------
    #Calcul de la carte ou l'on regarde suivant z
    #--------------------------------------------
    cam_z = Camera(center=center,line_of_sight_axis='z',region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector='y',map_max_size=512)

    rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_density"])
    rt = raytracing.RayTracer(amr,ro.info,rho_op)
    datamap = rt.process(cam_z, surf_qty=True)
    map_col = np.log10(datamap.map.T*lbox_cm)

    if v_proj == True:
        Vz_op = ScalarOperator(lambda dset: dset["vel"][...,2]*dset["rho"] ,  ro.info["unit_velocity"])
        rt = raytracing.RayTracer(amr,ro.info,Vz_op)
        datamap_vz = rt.process(cam_z, surf_qty=True)
        map_Vz = datamap_vz.map.T / datamap.map.T * factor_vel_km_s

        plt.figure()
        im = plt.imshow(map_col,extent=[(-radius+center[0])*lbox_au,(radius+center[0])*lbox_au,(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au],origin='lower', vmin=vmin_dens, vmax=vmax_dens)   
        plt.xlabel('$x$ (AU)')     
        plt.ylabel('$y$ (AU)')
        cbar=plt.colorbar()                                                                  
        cbar.set_label(r'$log(N) \, \, (cm^{-2})$')  
        if title_time==True:
            plt.title('Time = '+str(int(simulation_time))+' years')
        if title_time_cor==True:
            plt.title('Time = '+str(int(simulation_time))+' years \n Corrected time = '+str(int(simulation_time - ref[1]*1e6))+' years')
        if radius_zoom==5:
            plt.xlim([0,lbox_au])
            plt.ylim([0,lbox_au])
        if save==True:
            plt.savefig(path_save+'dens_z_'+str(radius_zoom)+'_'+str(num_output)+'.pdf', bbox_inches='tight')


    if v_proj == True:
        plt.figure()
        norm = MidpointNormalize(midpoint=0)  #Pour avoir le centre de la colormap a 0
        extent=[(-radius+center[0])*lbox_au,(radius+center[0])*lbox_au,(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au]
        plt.imshow(map_Vz,extent=extent,origin='lower',cmap='RdBu_r',norm=norm, vmin=vmin_vel, vmax=vmax_vel)
        plt.xlabel('$x$ (AU)')     
        plt.ylabel('$y$ (AU)')
        cbar = plt.colorbar()          
        cbar.set_label(r'$v_z \, (km.s^{-1})$')  
        if title_time==True:
            plt.title('Time = '+str(int(simulation_time))+' years')
        if title_time_cor==True:
            plt.title('Time = '+str(int(simulation_time))+' years \n Corrected time = '+str(int(simulation_time - ref[1]*1e6))+' years')
        if radius_zoom==5:
            plt.xlim([0,lbox_au])
            plt.ylim([0,lbox_au])
        if save==True:
            plt.savefig(path_save+'vel_z_'+str(radius_zoom)+'_'+str(num_output)+'.pdf', bbox_inches='tight')


        #Fit du gradient de vitesse
        fit_gradient_vel((map_Vz,extent))





if reposition_fig==True:
    for i in range(3):
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



