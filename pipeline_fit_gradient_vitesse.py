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

import h5py

import os
import pipeline_temps_0_simulation as t_0

import scipy.optimize as opt




#-------------------------------------------------------
#Entree le nom de la simulation et le numero de l'output
#-------------------------------------------------------
simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire_shr_bigbox_20pourc'

owner = 'averliat_alfven'
num_output = 40


dir_save = 'Gradient_vitesse'  #Repertoire de sauvegarde des figures

save_ang = True  #Pour sauvegarder les tableaux contenant les angles des gradients


#Caracteristiques des coquilles en AU:
R_min = 50
R_max = 10000  #='all' pour prendre toute la boite
dr = 50


file_save_ang = 'Angle_omega_J_gradient_vitesse_log_output'+str(num_output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.hdf5'
#file_save_ang = 'Angle_gradient_vitesse_output'+str(num_output)+'comparaison_mathilde.hdf5'


title_time=True
title_time_cor=False
seuil_rho = 1e-10


selon_x = True
selon_y = True
selon_z = True


fig_coldens = False   #Pour afficher et sauvegarder les images de densite colonne
save_coldens = True

fig_vel = False   #Pour afficher et sauvegarder les images de vitesses projetees
save_vel = True

verif_fit = False   #Pour verifier si le fit du gradient de vitesse est juste
save_fit = True



save = np.bool(save_coldens*fig_coldens+save_vel*fig_vel+save_fit*verif_fit)


coef_vel=5
vmin_vel = None#-0.25
vmax_vel = None#-vmin_vel
vmin_dens = None
vmax_dens = None




reposition_fig = False #Pour repositionner les figures ouvertes par matplotlib

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
lbox_pc = ro.info['unit_length'].express(cst.pc)
factor_time_yr = ro.info['unit_time'].express(cst.year)
factor_vel_km_s = ro.info['unit_velocity'].express(cst.km/cst.s)

simulation_time = ro.info['time']*factor_time_yr



#---------------------------------------------------------------------------------
#Calcul des t_0 de chaque simulation ou lecture des fichiers si calculs deja faits
#---------------------------------------------------------------------------------
if title_time_cor == True:
    if os.path.isfile(path_analyse+'t0_seuil_rho_'+str(seuil_rho)+'.txt') == False:
        ref = np.array(t_0.temps_0_simu(path, seuil_rho, sortie_output=1))
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
#if radius_zoom==5:
#    center = [0.5,0.5,0.5]
#else:
arg_centre = np.argmax(rho)
center = [pos[:,0][arg_centre],pos[:,1][arg_centre],pos[:,2][arg_centre]]

nbr_shell = int((R_max-R_min)/dr)  #nombre de coquilles
#shells_au=np.linspace(R_min,R_max,nbr_shell)  #rayon des coquilles en AU
shells_au=np.logspace(np.log10(R_min),np.log10(R_max),nbr_shell)  #rayon des coquilles en AU
shells=shells_au/lbox_au  #rayon des coquilles en unite de code (correspondant au zoom_v ci-dessous)
#nbr_shell = 5
#shells_au=np.array([50.,300.,800.,2000.,6000.])
#shells=shells_au/lbox_au


#zoom_v=[0.045, 0.015, 0.005, 0.005/3., lbox_pc*4]
#radius=zoom_v[radius_zoom-1]
#radius=float(zoom_v[np.where(zoom_v==radius_zoom)])     #0.015#0.005 #Niveau de zoom correspondant au niveau '3' des images de "pipeline_image_unique.py"


#==========================================================================================
#Definition prise sur 
#https://stackoverflow.com/questions/20144529/shifted-colorbar-matplotlib/20146989#20146989 
#pour avoir le centre de la colorbar a 0
#==========================================================================================
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


def fit_gradient_vel((map_V,extent), tag, radius=0, verif_fit=verif_fit, save_fit=save_fit, fig_vel=fig_vel, save_vel=save_vel):
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

    omega=np.sqrt(a**2+b**2)  #en km.s-1.AU-1
    omega=omega*lbox_au/lbox_pc  #en km.s-1.pc-1
    moment_spec = omega*(radius*lbox_pc)**2
    theta=np.arctan(b/a)  #orientation de la ligne de vitesse nulle
    vel=v0+a*tab_abs+b*tab_ord



    #Ajout de la fleche representant le gradient de vitesse
    pts_fleche=512
    x_fleche=np.linspace((center[0]-radius)*lbox_au, (center[0]+radius)*lbox_au,pts_fleche)
    C=(center[1]-np.tan(theta)*center[0])*lbox_au
    y_fleche=np.tan(theta)*x_fleche+C

    taille_fleche_x=np.sqrt( 16*(radius*lbox_au)**2/25. /(1+np.tan(theta)**2) )  #Fleche comprise entre 1/10eme et 9/10eme de la boite, ici demi taille
    mask_fleche=(x_fleche>center[0]*lbox_au-taille_fleche_x) & (x_fleche<center[0]*lbox_au+taille_fleche_x)

    #Necessaire pour determiner si la fleche va de droite a gauche ou de gauche a droite
    tab_abs=np.linspace(extent[0],extent[1],len(map_V))
    tab_ord=np.linspace(extent[2],extent[3],len(map_V))



    if np.sum(mask_fleche)<=1:  #Cas ou la fleche est trop verticale et donc le mask_fleche est vide...
        x_fleche_c=x_fleche[np.where(x_fleche==x_fleche.flat[np.abs(x_fleche-center[0]*lbox_au).argmin()])]
        x_fleche_c_ind=np.where(x_fleche==x_fleche_c)

        y_fleche=np.linspace((center[1]-radius)*lbox_au, (center[1]+radius)*lbox_au,pts_fleche)
        taille_fleche_y=np.sqrt( 16*(radius*lbox_au)**2/25. - taille_fleche_x**2 )
        mask_fleche=(y_fleche>center[1]*lbox_au-taille_fleche_y) & (y_fleche<center[1]*lbox_au+taille_fleche_y)

        ind_ord_max=np.where(tab_ord==tab_ord.flat[np.abs(tab_ord - y_fleche[mask_fleche][-1]).argmin()])
        ind_ord_min=np.where(tab_ord==tab_ord.flat[np.abs(tab_ord - y_fleche[mask_fleche][0]).argmin()])


        vel_haut=vel[x_fleche_c_ind,ind_ord_max]
        vel_bas=vel[x_fleche_c_ind,ind_ord_min]
        if vel_haut > vel_bas:
            if fig_vel==True:
                plt.arrow(x_fleche_c,y_fleche[mask_fleche][0],0.,y_fleche[mask_fleche][-1]-y_fleche[mask_fleche][0],head_width=10*radius*lbox_au/172., head_length=10*radius*lbox_au/172.,linestyle="-",linewidth=1.1,color='g',length_includes_head=True,joinstyle='round',overhang=0.7)
                if save_vel==True:
                    plt.savefig(path_save+'vel_'+tag+'_'+str(int(shells_au[ind]))+'_'+str(num_output)+'.pdf', bbox_inches='tight')
            ang=theta*180/np.pi

        elif vel_haut < vel_bas:
            if fig_vel==True:
                plt.arrow(x_fleche_c,y_fleche[mask_fleche][-1],0.,y_fleche[mask_fleche][0]-y_fleche[mask_fleche][-1],head_width=10*radius*lbox_au/172., head_length=10*radius*lbox_au/172.,linestyle="-",linewidth=1.1,color='g',length_includes_head=True,joinstyle='round',overhang=0.7)
                if save_vel==True:
                    plt.savefig(path_save+'vel_'+tag+'_'+str(int(shells_au[ind]))+'_'+str(num_output)+'.pdf', bbox_inches='tight')
            if theta>0:
                ang=theta*180/np.pi - 180
            elif theta<0:
                ang=theta*180/np.pi + 180           

        if verif_fit==True:
            plt.figure()
            plt.imshow(vel,extent=extent,origin='lower',cmap='RdBu_r')
            if vel_haut > vel_bas:
                plt.arrow(x_fleche_c,y_fleche[mask_fleche][0],0.,y_fleche[mask_fleche][-1]-y_fleche[mask_fleche][0],head_width=10*radius*lbox_au/172., head_length=10*radius*lbox_au/172.,linestyle="-",linewidth=1.1,color='g',length_includes_head=True,joinstyle='round',overhang=0.7)
            elif vel_haut < vel_bas:
                plt.arrow(x_fleche_c,y_fleche[mask_fleche][-1],0.,y_fleche[mask_fleche][0]-y_fleche[mask_fleche][-1],head_width=10*radius*lbox_au/172., head_length=10*radius*lbox_au/172.,linestyle="-",linewidth=1.1,color='g',length_includes_head=True,joinstyle='round',overhang=0.7)
            plt.xlabel('(AU)')
            plt.ylabel('(AU)')
            cbar = plt.colorbar()
            cbar.set_label(r'$v \, (km.s^{-1})$')
            if title_time==True:
                plt.title('Time = '+str(int(simulation_time))+' years')
            if title_time_cor==True:
                plt.title('Time = '+str(int(simulation_time))+' years \n Corrected time = '+str(int(simulation_time - ref[1]*1e6))+' years')
            if save_fit==True:
                plt.savefig(path_save+'fit_'+tag+'_'+str(int(shells_au[ind]))+'_'+str(num_output)+'.pdf', bbox_inches='tight')
    


    else:            
        ind_abs_max=np.where(tab_abs==tab_abs.flat[np.abs(tab_abs - x_fleche[mask_fleche][-1]).argmin()])
        ind_ord_max=np.where(tab_ord==tab_ord.flat[np.abs(tab_ord - y_fleche[mask_fleche][-1]).argmin()])

        ind_abs_min=np.where(tab_abs==tab_abs.flat[np.abs(tab_abs - x_fleche[mask_fleche][0]).argmin()])
        ind_ord_min=np.where(tab_ord==tab_ord.flat[np.abs(tab_ord - y_fleche[mask_fleche][0]).argmin()])

        vel_droite = vel[ind_ord_max,ind_abs_max]
        vel_gauche = vel[ind_ord_min,ind_abs_min]


        if vel_droite > vel_gauche:
            if fig_vel==True:
                plt.arrow(x_fleche[mask_fleche][0],y_fleche[mask_fleche][0],x_fleche[mask_fleche][-1]-x_fleche[mask_fleche][0],y_fleche[mask_fleche][-1]-y_fleche[mask_fleche][0],head_width=10*radius*lbox_au/172., head_length=10*radius*lbox_au/172.,linestyle="-",linewidth=1.1,color='g',length_includes_head=True,joinstyle='round',overhang=0.7)
                if save_vel==True:
                    plt.savefig(path_save+'vel_'+tag+'_'+str(int(shells_au[ind]))+'_'+str(num_output)+'.pdf', bbox_inches='tight')
            ang=theta*180/np.pi

        elif vel_droite < vel_gauche:
            if fig_vel==True:
                plt.arrow(x_fleche[mask_fleche][-1],y_fleche[mask_fleche][-1],x_fleche[mask_fleche][0]-x_fleche[mask_fleche][-1],y_fleche[mask_fleche][0]-y_fleche[mask_fleche][-1],head_width=10*radius*lbox_au/172., head_length=10*radius*lbox_au/172.,linestyle="-",linewidth=1.1,color='g',length_includes_head=True,joinstyle='round',overhang=0.7)
                if save_vel==True:
                    plt.savefig(path_save+'vel_'+tag+'_'+str(int(shells_au[ind]))+'_'+str(num_output)+'.pdf', bbox_inches='tight')
            if theta>0:
                ang=theta*180/np.pi - 180
            elif theta<0:
                ang=theta*180/np.pi + 180



        #Trace de la carte du gradient de vitesse ajuste pour verification
        if verif_fit==True:
            plt.figure()
            plt.imshow(vel,extent=extent,origin='lower',cmap='RdBu_r')
            if vel_droite > vel_gauche:
                plt.arrow(x_fleche[mask_fleche][0],y_fleche[mask_fleche][0],x_fleche[mask_fleche][-1]-x_fleche[mask_fleche][0],y_fleche[mask_fleche][-1]-y_fleche[mask_fleche][0],head_width=10*radius*lbox_au/172., head_length=10*radius*lbox_au/172.,linestyle="-",linewidth=1.1,color='g',length_includes_head=True,joinstyle='round',overhang=0.7)
            elif vel_droite < vel_gauche:
                plt.arrow(x_fleche[mask_fleche][-1],y_fleche[mask_fleche][-1],x_fleche[mask_fleche][0]-x_fleche[mask_fleche][-1],y_fleche[mask_fleche][0]-y_fleche[mask_fleche][-1],head_width=10*radius*lbox_au/172., head_length=10*radius*lbox_au/172.,linestyle="-",linewidth=1.1,color='g',length_includes_head=True,joinstyle='round',overhang=0.7)

            plt.xlabel('(AU)')
            plt.ylabel('(AU)')
            cbar = plt.colorbar()
            cbar.set_label(r'$v \, (km.s^{-1})$')
            if title_time==True:
                plt.title('Time = '+str(int(simulation_time))+' years')
            if title_time_cor==True:
                plt.title('Time = '+str(int(simulation_time))+' years \n Corrected time = '+str(int(simulation_time - ref[1]*1e6))+' years')

            if save_fit==True:
                plt.savefig(path_save+'fit_'+tag+'_'+str(int(shells_au[ind]))+'_'+str(num_output)+'.pdf', bbox_inches='tight')

    
    return (ang, omega, moment_spec)  #en degre
#=============================================================================











#=============================================================================
#=============================================================================
#Boucle sur les coquilles de rayons croissants
#=============================================================================
#=============================================================================

#Initialisation des tableaux de stockage
ang_x_tab=np.zeros(nbr_shell)
ang_y_tab=np.zeros(nbr_shell)
ang_z_tab=np.zeros(nbr_shell)
omega_x_tab=np.zeros(nbr_shell)
omega_y_tab=np.zeros(nbr_shell)
omega_z_tab=np.zeros(nbr_shell)
moment_spec_x_tab=np.zeros(nbr_shell)
moment_spec_y_tab=np.zeros(nbr_shell)
moment_spec_z_tab=np.zeros(nbr_shell)
ind=0
verif=0

#Debut de la boucle sur les coquilles
for radius in shells:
    if selon_x==True:
        #--------------------------------------------
        #Calcul de la carte ou l'on regarde suivant x
        #--------------------------------------------
        cam_x = Camera(center=center,line_of_sight_axis='x',region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector='z',map_max_size=512)

        rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_density"])
        rt = raytracing.RayTracer(amr,ro.info,rho_op)
        datamap = rt.process(cam_x, surf_qty=True)

        if fig_coldens==True:
            map_col = np.log10(datamap.map.T*lbox_cm)

            plt.figure()
            plt.imshow(map_col,extent=[(-radius+center[0])*lbox_au,(radius+center[0])*lbox_au,(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au],origin='lower', vmin=vmin_dens, vmax=vmax_dens)
            plt.xlabel('$y$ (AU)')
            plt.ylabel('$z$ (AU)')
            cbar = plt.colorbar()
            cbar.set_label(r'$log(N) \, \, (cm^{-2})$')
            if title_time==True:
                plt.title('Time = '+str(int(simulation_time))+' years')
            if title_time_cor==True:
                plt.title('Time = '+str(int(simulation_time))+' years \n Corrected time = '+str(int(simulation_time - ref[1]*1e6))+' years')
            if save_coldens==True:
                plt.savefig(path_save+'dens_x_'+str(int(shells_au[ind]))+'_'+str(num_output)+'.pdf', bbox_inches='tight')



        Vx_op = ScalarOperator(lambda dset: dset["vel"][...,0]*dset["rho"] ,  ro.info["unit_velocity"])
        rt = raytracing.RayTracer(amr,ro.info,Vx_op)
        datamap_vx = rt.process(cam_x, surf_qty=True)
        map_Vx = datamap_vx.map.T / datamap.map.T * factor_vel_km_s

        extent=[(-radius+center[0])*lbox_au,(radius+center[0])*lbox_au,(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au]


        if fig_vel==True:
            plt.figure()
            norm = MidpointNormalize(midpoint=0)  #Pour avoir le centre de la colormap a 0
            mask=np.zeros((len(map_Vx),)*2)
            for i in range(len(map_Vx)):
                for j in range(len(map_Vx)):
                    mask[i,j]=map_Vx[i,j]>0
            r=map_Vx*mask
            vmax_vel=coef_vel*np.median(r[np.where(r>0)])
            mask=-(mask-1)
            r=map_Vx*mask
            vmin_vel=coef_vel*np.median(r[np.where(r<0)])
            plt.imshow(map_Vx,extent=extent,origin='lower',cmap='RdBu_r',norm=norm, vmin=vmin_vel, vmax=vmax_vel)
            plt.xlabel('$y$ (AU)')
            plt.ylabel('$z$ (AU)')
            cbar = plt.colorbar()
            cbar.set_label(r'$v_x \, (km.s^{-1})$')
            if title_time==True:
                plt.title('Time = '+str(int(simulation_time))+' years')
            if title_time_cor==True:
                plt.title('Time = '+str(int(simulation_time))+' years \n Corrected time = '+str(int(simulation_time - ref[1]*1e6))+' years')


        #Fit du gradient de vitesse
        (ang_x_tab[ind],omega_x_tab[ind],moment_spec_x_tab[ind])=fit_gradient_vel((map_Vx,extent), 'x',radius=radius)





    if selon_y==True:
        #--------------------------------------------
        #Calcul de la carte ou l'on regarde suivant y
        #--------------------------------------------
        cam_y = Camera(center=center,line_of_sight_axis='y',region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector='x',map_max_size=512)

        rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_density"])
        rt = raytracing.RayTracer(amr,ro.info,rho_op)
        datamap = rt.process(cam_y, surf_qty=True)

        if fig_coldens==True:
            map_col = np.log10(datamap.map.T*lbox_cm)

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
            if save_coldens==True:
                plt.savefig(path_save+'dens_y_'+str(int(shells_au[ind]))+'_'+str(num_output)+'.pdf', bbox_inches='tight')



        Vy_op = ScalarOperator(lambda dset: dset["vel"][...,1]*dset["rho"] ,  ro.info["unit_velocity"])
        rt = raytracing.RayTracer(amr,ro.info,Vy_op)
        datamap_vy = rt.process(cam_y, surf_qty=True)
        map_Vy = datamap_vy.map.T / datamap.map.T * factor_vel_km_s

        extent=[(-radius+center[0])*lbox_au,(radius+center[0])*lbox_au,(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au]


        if fig_vel==True:
            plt.figure()
            norm = MidpointNormalize(midpoint=0)  #Pour avoir le centre de la colormap a 0
            mask=np.zeros((len(map_Vy),)*2)
            for i in range(len(map_Vy)):
                for j in range(len(map_Vy)):
                    mask[i,j]=map_Vy[i,j]>0
            r=map_Vy*mask
            vmax_vel=coef_vel*np.median(r[np.where(r>0)])
            mask=-(mask-1)
            r=map_Vy*mask
            vmin_vel=coef_vel*np.median(r[np.where(r<0)])
            plt.imshow(map_Vy,extent=extent,origin='lower',cmap='RdBu_r',norm=norm, vmin=vmin_vel, vmax=vmax_vel)
            plt.xlabel('$z$ (AU)')
            plt.ylabel('$x$ (AU)')
            cbar = plt.colorbar()
            cbar.set_label(r'$v_y \, (km.s^{-1})$')
            if title_time==True:
                plt.title('Time = '+str(int(simulation_time))+' years')
            if title_time_cor==True:
                plt.title('Time = '+str(int(simulation_time))+' years \n Corrected time = '+str(int(simulation_time - ref[1]*1e6))+' years')  


        #Fit du gradient de vitesse
        (ang_y_tab[ind],omega_y_tab[ind],moment_spec_y_tab[ind])=fit_gradient_vel((map_Vy,extent), 'y',radius=radius)





    if selon_z==True:
        #--------------------------------------------
        #Calcul de la carte ou l'on regarde suivant z
        #--------------------------------------------
        cam_z = Camera(center=center,line_of_sight_axis='z',region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector='y',map_max_size=512)

        rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_density"])
        rt = raytracing.RayTracer(amr,ro.info,rho_op)
        datamap = rt.process(cam_z, surf_qty=True)

        if fig_coldens==True:
            map_col = np.log10(datamap.map.T*lbox_cm)

            plt.figure()
            plt.imshow(map_col,extent=[(-radius+center[0])*lbox_au,(radius+center[0])*lbox_au,(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au],origin='lower', vmin=vmin_dens, vmax=vmax_dens)
            plt.xlabel('$x$ (AU)')
            plt.ylabel('$y$ (AU)')
            cbar=plt.colorbar()
            cbar.set_label(r'$log(N) \, \, (cm^{-2})$')
            if title_time==True:
                plt.title('Time = '+str(int(simulation_time))+' years')
            if title_time_cor==True:
                plt.title('Time = '+str(int(simulation_time))+' years \n Corrected time = '+str(int(simulation_time - ref[1]*1e6))+' years')
            if save_coldens==True:
                plt.savefig(path_save+'dens_z_'+str(int(shells_au[ind]))+'_'+str(num_output)+'.pdf', bbox_inches='tight')



        Vz_op = ScalarOperator(lambda dset: dset["vel"][...,2]*dset["rho"] ,  ro.info["unit_velocity"])
        rt = raytracing.RayTracer(amr,ro.info,Vz_op)
        datamap_vz = rt.process(cam_z, surf_qty=True)
        map_Vz = datamap_vz.map.T / datamap.map.T * factor_vel_km_s

        extent=[(-radius+center[0])*lbox_au,(radius+center[0])*lbox_au,(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au]


        if fig_vel==True:
            plt.figure()
            norm = MidpointNormalize(midpoint=0)  #Pour avoir le centre de la colormap a 0
            mask=np.zeros((len(map_Vz),)*2)
            for i in range(len(map_Vz)):
                for j in range(len(map_Vz)):
                    mask[i,j]=map_Vz[i,j]>0
            r=map_Vz*mask
            vmax_vel=coef_vel*np.median(r[np.where(r>0)])
            mask=-(mask-1)
            r=map_Vz*mask
            vmin_vel=coef_vel*np.median(r[np.where(r<0)])
            plt.imshow(map_Vz,extent=extent,origin='lower',cmap='RdBu_r',norm=norm, vmin=vmin_vel, vmax=vmax_vel)
            plt.xlabel('$x$ (AU)')
            plt.ylabel('$y$ (AU)')
            cbar = plt.colorbar()
            cbar.set_label(r'$v_z \, (km.s^{-1})$')
            if title_time==True:
                plt.title('Time = '+str(int(simulation_time))+' years')
            if title_time_cor==True:
                plt.title('Time = '+str(int(simulation_time))+' years \n Corrected time = '+str(int(simulation_time - ref[1]*1e6))+' years')


        #Fit du gradient de vitesse
        (ang_z_tab[ind],omega_z_tab[ind],moment_spec_z_tab[ind])=fit_gradient_vel((map_Vz,extent), 'z',radius=radius)


    ind += 1
    verif=1





if (save_ang == True and verif==1):
    h5f = h5py.File(path_analyse+file_save_ang, 'w')

    h5f.create_dataset('shells_au', data=shells_au)
    h5f.create_dataset('ang_x_tab',data=ang_x_tab)
    h5f.create_dataset('ang_y_tab',data=ang_y_tab)
    h5f.create_dataset('ang_z_tab',data=ang_z_tab)
    h5f.create_dataset('omega_x_tab',data=omega_x_tab)
    h5f.create_dataset('omega_y_tab',data=omega_y_tab)
    h5f.create_dataset('omega_z_tab',data=omega_z_tab)
    h5f.create_dataset('moment_spec_x_tab',data=moment_spec_x_tab)
    h5f.create_dataset('moment_spec_y_tab',data=moment_spec_y_tab)
    h5f.create_dataset('moment_spec_z_tab',data=moment_spec_z_tab)

    h5f.close()





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









