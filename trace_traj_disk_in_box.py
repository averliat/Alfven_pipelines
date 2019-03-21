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
import h5py




#-------------------------------------------------------
#Entree le nom de la simulation et le numero de l'output
#-------------------------------------------------------
simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire_shr_bigbox_50pourc'
owner = 'averliat_alfven'

output_min = 'None'
output_max = 110

num_output = 1  #Output par dessus lequel sera overplot la trajectoire du disque

file_save = 'Traj_disk_in_box_output_'+str(output_min)+'_'+str(output_max)+'.hdf5'

selon_x = True
selon_y = True
selon_z = True



#-------------------
#Ne pas toucher a ca
#-------------------
save = False
dir_save = 'Trajectoire_disque_dans_boite'

radius_zoom = 5

v_proj = False

title_time=False
title_time_cor=False
seuil_rho = 1e-10

fleche_vel = False
nbre_fleche = 20  #En fait plus c'est faible plus y'a de fleches...

vmin_vel = None#-0.25
vmax_vel = None#-vmin_vel

vmin_dens = None
vmax_dens = None

color_sink_colmap = 'firebrick'
transparence_sink_colmap = 0.4
color_sink_velmap = 'limegreen'
transparence_sink_velmap = 0.7

reposition_fig = False  #Pour repositionner les figures ouvertes par matplotlib

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
path_analyses='/home/averliat/these/analyses/'+simu+'/'

#if simu == 'B335_noturb_norot_hydro_pert_asym_aleatoire50_vhr':
#   path_t0='/mnt/magmist/magmist/simu_B335_averliat/'+simu+'/'
#else:
#   path_t0=path
path_t0=path

if save==True:
    if os.path.isdir(path_save) == False:
        os.mkdir(path_save)

if output_min=='None':
    ref=np.loadtxt(path_analyses+'t0_seuil_rho_'+str(seuil_rho)+'.txt')
    output_0 = int(ref[0])
    time_0 = ref[1]
    output_min = output_0
    file_save = 'Traj_disk_in_box_output_'+str(output_min)+'_'+str(output_max)+'.hdf5'



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
    

zoom_v=[0.045, 0.015, 0.005, 0.005/3., 0.5]


if 'bigbox' or 'jet' in simu:
    zoom_v=[0.045/2, 0.015/2, 0.005/2, 0.005/3./2, 0.5, 0.5/2, 0.005/20.]
    if radius_zoom==6:
        center = [0.5,0.5,0.5]

if 'hugebox' in simu:
    zoom_v=[0.045/4, 0.015/4, 0.005/4, 0.005/3./4, 0.5, 0.5/4]
    if radius_zoom==6:
        center = [0.5,0.5,0.5]

radius=zoom_v[radius_zoom-1]
#radius=float(zoom_v[np.where(zoom_v==radius_zoom)])     #0.015#0.005 #Niveau de zoom correspondant au niveau '3' des images de "pipeline_image_unique.py"



#------------------------------------------------------------------------
#Get the properties of the particules with just the particules' positions
#Copie depuis module_extract.py
#------------------------------------------------------------------------
def read_sink_csv(num,directory,no_cvs=None):

    name = directory+'output_'+str(num).zfill(5)+'/sink_'+str(num).zfill(5)+'.csv'
    print 'read sinks from ', name

    if(no_cvs is None):
        sinks = np.loadtxt(name,delimiter=',',ndmin=2,usecols=(0,1,2,3,4,5,6,7,8)) #,9,10,11,12))
    else:
        sinks = np.loadtxt(name,ndmin=2,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12))

    return sinks



#------------------------------
#Pour visualiser les sinks
#Adapted from module_extract.py
#------------------------------
def visualise_sink_AU(path, num_output):
    sinks=read_sink_csv(num_output,path)
    if len(sinks) != 0:
        mass_sinks=sinks[:,1]  #masse des sinks en masses solaires
        x_sinks=sinks[:,3]*lbox_au/lbox  #coordonnees x des sinks en AU
        y_sinks=sinks[:,4]*lbox_au/lbox  #coordonnees y des sinks en AU
        z_sinks=sinks[:,5]*lbox_au/lbox  #coordonnees z des sinks en AU

        return mass_sinks, x_sinks, y_sinks, z_sinks



#---------------------------------------
#Lecture du fichier sink.csv s'il existe
#---------------------------------------
sink_plot=False
if os.path.isfile(path+'output_'+str(num_output).zfill(5)+'/sink_'+str(num_output).zfill(5)+'.csv') == True:
    #sink_plot=True 
    sinks=read_sink_csv(num_output,path)
    if len(sinks) != 0:
        sink_plot=True
    if sink_plot==True:
        m_sinks,x_sinks,y_sinks,z_sinks=visualise_sink_AU(path,num_output)
        #size_sinks=20*np.tanh(25*m_sinks)
        size_sinks=20*np.sqrt(m_sinks)

    if radius_zoom != 5:  #centrage sur la plus grosse sink au lieu de rho max
        arg_biggest_sink=np.argmax(m_sinks)
        x_biggest_sink=x_sinks[arg_biggest_sink]/lbox_au
        y_biggest_sink=y_sinks[arg_biggest_sink]/lbox_au
        z_biggest_sink=z_sinks[arg_biggest_sink]/lbox_au
        center=[x_biggest_sink,y_biggest_sink,z_biggest_sink]




#------------------------------------
#Chargement des differentes quantites
#------------------------------------
h5f = h5py.File(path_analyses+file_save, 'r')

time_tab=h5f['time_tab'][:]
time_cor_tab=h5f['time_cor_tab'][:]
numero_output=h5f['numero_output'][:]
max_rho_tab=h5f['max_rho_tab'][:]*lbox_au

h5f.close()



#--------------------------------------------
#Determination de la fenetre de visualisation
#--------------------------------------------
extent_max_x = np.max(max_rho_tab[:,0])
extent_min_x = np.min(max_rho_tab[:,0])
extent_x = extent_max_x - extent_min_x 

extent_max_y = np.max(max_rho_tab[:,1])
extent_min_y = np.min(max_rho_tab[:,1])
extent_y = extent_max_y - extent_min_y 

extent_max_z = np.max(max_rho_tab[:,2])
extent_min_z = np.min(max_rho_tab[:,2])
extent_z = extent_max_z - extent_min_z 


extent_max_x += 0.05*extent_x 
extent_min_x -= 0.05*extent_x

extent_max_y += 0.05*extent_y 
extent_min_y -= 0.05*extent_y

extent_max_z += 0.05*extent_z 
extent_min_z -= 0.05*extent_z

cmappts = plt.get_cmap('magma')
colorspts = [cmappts(i) for i in np.linspace(0,1,len(time_tab))]



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


    if fleche_vel == True:
        Vy_depuis_x_op = ScalarOperator(lambda dset: dset["vel"][...,1]*dset["rho"] ,  ro.info["unit_velocity"])
        rt = raytracing.RayTracer(amr,ro.info,Vy_depuis_x_op)
        datamap_vy_depuis_x = rt.process(cam_x, surf_qty=True)
        map_vy_depuis_x = datamap_vy_depuis_x.map.T / datamap.map.T * factor_vel_km_s
        map_vy_depuis_x_red = map_vy_depuis_x[::nbre_fleche, ::nbre_fleche]

        Vz_depuis_x_op = ScalarOperator(lambda dset: dset["vel"][...,2]*dset["rho"] ,  ro.info["unit_velocity"])
        rt = raytracing.RayTracer(amr,ro.info,Vz_depuis_x_op)
        datamap_vz_depuis_x = rt.process(cam_x, surf_qty=True)
        map_vz_depuis_x = datamap_vz_depuis_x.map.T / datamap.map.T * factor_vel_km_s
        map_vz_depuis_x_red = map_vz_depuis_x[::nbre_fleche, ::nbre_fleche]


        plt.figure()
        im = plt.imshow(map_col,extent=[(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au,(-radius+center[2])*lbox_au,(radius+center[2])*lbox_au],origin='lower')
        plt.xlabel('$y$ (AU)')
        plt.ylabel('$z$ (AU)')
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
            plt.savefig(path_save+'dens_x_'+str(radius_zoom)+'_'+str(num_output)+'.pdf', bbox_inches='tight')  

        nx = map_vy_depuis_x_red.shape[0]
        ny = map_vy_depuis_x_red.shape[1]
        vec_x = (np.arange(nx)*2./nx*radius - radius + center[0] + radius/nx)*lbox_au
        vec_y = (np.arange(ny)*2./ny*radius - radius + center[1] + radius/nx)*lbox_au
        xx,yy = np.meshgrid(vec_x,vec_y)

        plt.quiver(xx,yy,map_vy_depuis_x_red,map_vz_depuis_x_red)


    else:
        plt.figure()
        plt.imshow(map_col,extent=[(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au,(-radius+center[2])*lbox_au,(radius+center[2])*lbox_au],origin='lower', vmin=vmin_dens, vmax=vmax_dens)   
        if sink_plot==True:
            for i in range(len(m_sinks)):
                if (y_sinks[i]>(-radius+center[1])*lbox_au)&(y_sinks[i]<(radius+center[1])*lbox_au)&(z_sinks[i]>(-radius+center[2])*lbox_au)&(z_sinks[i]<(radius+center[2])*lbox_au):
                    plt.plot(y_sinks[i],z_sinks[i],'.',color=color_sink_colmap,markersize=size_sinks[i],alpha=transparence_sink_colmap)

        for i in range(len(time_tab)):
            plt.plot(max_rho_tab[i][1],max_rho_tab[i][2],'.',color=colorspts[i])


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
        plt.imshow(map_Vx,extent=[(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au,(-radius+center[2])*lbox_au,(radius+center[2])*lbox_au],origin='lower',cmap='RdBu_r',norm=norm, vmin=vmin_vel, vmax=vmax_vel)
        if sink_plot==True:
            for i in range(len(m_sinks)):
                if (y_sinks[i]>(-radius+center[1])*lbox_au)&(y_sinks[i]<(radius+center[1])*lbox_au)&(z_sinks[i]>(-radius+center[2])*lbox_au)&(z_sinks[i]<(radius+center[2])*lbox_au):
                    plt.plot(y_sinks[i],z_sinks[i],'.',color=color_sink_velmap,markersize=size_sinks[i],alpha=transparence_sink_velmap)
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


    if fleche_vel == True:
        Vx_depuis_y_op = ScalarOperator(lambda dset: dset["vel"][...,0]*dset["rho"] ,  ro.info["unit_velocity"])
        rt = raytracing.RayTracer(amr,ro.info,Vx_depuis_y_op)
        datamap_vx_depuis_y = rt.process(cam_y, surf_qty=True)
        map_vx_depuis_y = datamap_vx_depuis_y.map.T / datamap.map.T * factor_vel_km_s
        map_vx_depuis_y_red = map_vx_depuis_y[::nbre_fleche, ::nbre_fleche]

        Vz_depuis_y_op = ScalarOperator(lambda dset: dset["vel"][...,2]*dset["rho"] ,  ro.info["unit_velocity"])
        rt = raytracing.RayTracer(amr,ro.info,Vz_depuis_y_op)
        datamap_vz_depuis_y = rt.process(cam_y, surf_qty=True)
        map_vz_depuis_y = datamap_vz_depuis_y.map.T / datamap.map.T * factor_vel_km_s
        map_vz_depuis_y_red = map_vz_depuis_y[::nbre_fleche, ::nbre_fleche]


        plt.figure()
        im = plt.imshow(map_col,extent=[(-radius+center[2])*lbox_au,(radius+center[2])*lbox_au,(-radius+center[0])*lbox_au,(radius+center[0])*lbox_au],origin='lower')
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

        nx = map_vx_depuis_y_red.shape[0]
        ny = map_vx_depuis_y_red.shape[1]
        vec_x = (np.arange(nx)*2./nx*radius - radius + center[0] + radius/nx)*lbox_au
        vec_y = (np.arange(ny)*2./ny*radius - radius + center[1] + radius/nx)*lbox_au
        xx,yy = np.meshgrid(vec_x,vec_y)

        plt.quiver(xx,yy,map_vz_depuis_y_red,map_vx_depuis_y_red)


    else:
        plt.figure()
        im = plt.imshow(map_col,extent=[(-radius+center[2])*lbox_au,(radius+center[2])*lbox_au,(-radius+center[0])*lbox_au,(radius+center[0])*lbox_au],origin='lower', vmin=vmin_dens, vmax=vmax_dens)   
        if sink_plot==True:
            for i in range(len(m_sinks)):
                if (z_sinks[i]>(-radius+center[2])*lbox_au)&(z_sinks[i]<(radius+center[2])*lbox_au)&(x_sinks[i]>(-radius+center[0])*lbox_au)&(x_sinks[i]<(radius+center[0])*lbox_au):
                    plt.plot(z_sinks[i],x_sinks[i],'.',color=color_sink_colmap,markersize=size_sinks[i],alpha=transparence_sink_colmap)

        for i in range(len(time_tab)):
            plt.plot(max_rho_tab[i][2],max_rho_tab[i][0],'.',color=colorspts[i])


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
        plt.imshow(map_Vy,extent=[(-radius+center[2])*lbox_au,(radius+center[2])*lbox_au,(-radius+center[0])*lbox_au,(radius+center[0])*lbox_au],origin='lower',cmap='RdBu_r',norm=norm, vmin=vmin_vel, vmax=vmax_vel)
        if sink_plot==True:
            for i in range(len(m_sinks)):
                if (z_sinks[i]>(-radius+center[2])*lbox_au)&(z_sinks[i]<(radius+center[2])*lbox_au)&(x_sinks[i]>(-radius+center[0])*lbox_au)&(x_sinks[i]<(radius+center[0])*lbox_au):
                    plt.plot(z_sinks[i],x_sinks[i],'.',color=color_sink_velmap,markersize=size_sinks[i],alpha=transparence_sink_velmap)
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


    if fleche_vel == True:
        Vx_depuis_z_op = ScalarOperator(lambda dset: dset["vel"][...,0]*dset["rho"] ,  ro.info["unit_velocity"])
        rt = raytracing.RayTracer(amr,ro.info,Vx_depuis_z_op)
        datamap_vx_depuis_z = rt.process(cam_z, surf_qty=True)
        map_vx_depuis_z = datamap_vx_depuis_z.map.T / datamap.map.T * factor_vel_km_s
        map_vx_depuis_z_red = map_vx_depuis_z[::nbre_fleche, ::nbre_fleche]

        Vy_depuis_z_op = ScalarOperator(lambda dset: dset["vel"][...,1]*dset["rho"] ,  ro.info["unit_velocity"])
        rt = raytracing.RayTracer(amr,ro.info,Vy_depuis_z_op)
        datamap_vy_depuis_z = rt.process(cam_z, surf_qty=True)
        map_vy_depuis_z = datamap_vy_depuis_z.map.T / datamap.map.T * factor_vel_km_s
        map_vy_depuis_z_red = map_vy_depuis_z[::nbre_fleche, ::nbre_fleche]


        plt.figure()
        im = plt.imshow(map_col,extent=[(-radius+center[0])*lbox_au,(radius+center[0])*lbox_au,(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au],origin='lower')
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

        nx = map_vx_depuis_z_red.shape[0]
        ny = map_vx_depuis_z_red.shape[1]
        vec_x = (np.arange(nx)*2./nx*radius - radius + center[0] + radius/nx)*lbox_au
        vec_y = (np.arange(ny)*2./ny*radius - radius + center[1] + radius/nx)*lbox_au
        xx,yy = np.meshgrid(vec_x,vec_y)

        plt.quiver(xx,yy,map_vx_depuis_z_red,map_vy_depuis_z_red)


    else:
        plt.figure()
        im = plt.imshow(map_col,extent=[(-radius+center[0])*lbox_au,(radius+center[0])*lbox_au,(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au],origin='lower', vmin=vmin_dens, vmax=vmax_dens)   
        if sink_plot==True:
            for i in range(len(m_sinks)):
                if (x_sinks[i]>(-radius+center[0])*lbox_au)&(x_sinks[i]<(radius+center[0])*lbox_au)&(y_sinks[i]>(-radius+center[1])*lbox_au)&(y_sinks[i]<(radius+center[1])*lbox_au):
                    plt.plot(x_sinks[i],y_sinks[i],'.',color=color_sink_colmap,markersize=size_sinks[i],alpha=transparence_sink_colmap)

        for i in range(len(time_tab)):
            plt.plot(max_rho_tab[i][0],max_rho_tab[i][1],'.',color=colorspts[i])


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
        plt.imshow(map_Vz,extent=[(-radius+center[0])*lbox_au,(radius+center[0])*lbox_au,(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au],origin='lower',cmap='RdBu_r',norm=norm, vmin=vmin_vel, vmax=vmax_vel)
        if sink_plot==True:
            for i in range(len(m_sinks)):
                if (x_sinks[i]>(-radius+center[0])*lbox_au)&(x_sinks[i]<(radius+center[0])*lbox_au)&(y_sinks[i]>(-radius+center[1])*lbox_au)&(y_sinks[i]<(radius+center[1])*lbox_au):
                    plt.plot(x_sinks[i],y_sinks[i],'.',color=color_sink_velmap,markersize=size_sinks[i],alpha=transparence_sink_velmap)
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






if v_proj==True and reposition_fig==True:
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



plt.figure(1)
plt.xlim((extent_min_y,extent_max_y))
plt.ylim((extent_min_z,extent_max_z))

plt.figure(2)
plt.xlim((extent_min_z,extent_max_z))
plt.ylim((extent_min_x,extent_max_x))

plt.figure(3)
plt.xlim((extent_min_x,extent_max_x))
plt.ylim((extent_min_y,extent_max_y))


