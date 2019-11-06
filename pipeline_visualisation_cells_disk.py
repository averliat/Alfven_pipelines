import numpy as np
import extract_disk as ed

import pymses
import sys
import glob as glob
import os

from pymses.filters import CellsToPoints
from pymses.utils import constants as cst
from pymses.analysis import Camera, raytracing, slicing
from pymses.analysis import ScalarOperator, FractionOperator, MaxLevelOperator
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rc
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.ion()
from matplotlib.colors import Normalize

#from astropy.io import fits

import os
import pipeline_temps_0_simulation as t_0

plt.style.use("pdf")
plt.style.use("aanda_modif")


#-------------------------------------------------------
#Entree le nom de la simulation et le numero de l'output
#-------------------------------------------------------
simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire_bigbox_50pourc_sink_seuil_haut_MHD_lr'

owner = 'averliat_alfven'
num_output = 114

#-------------------------------------------------------
#-------------------------------------------------------
save = True
dir_save = 'Masque_identification_disque'

radius_zoom = 3

v_proj = False

title_time=False
title_time_cor=True
seuil_rho = 1e-10

fleche_vel = False
nbre_fleche = 20  #En fait plus c'est faible plus y'a de fleches...


selon_x = False
selon_y = False
selon_z = False
selon_coord = True
face_on=True
edge_on=True

#=========================================
#Calcul des directions face-on et edge-on:
#=========================================
#Pour recherche du edge-on : mettre view='edge-on' et changer les coord ci-dessous
#Ensuite faire la figure face-on
view_diagram_edgeon=np.array([0,1,0]) #edge-on
up_vector_edgeon=np.array([1,0,0]) #edge-on

norm_view_diagram_edgeon_2=np.sum(view_diagram_edgeon**2)
view_diagram_faceon=up_vector_edgeon-np.sum(up_vector_edgeon*view_diagram_edgeon)*view_diagram_edgeon/norm_view_diagram_edgeon_2
up_vector_faceon=view_diagram_edgeon

#vd(output59)=[1,0,2]
#uv(output59)=[0,-0.15,1]

#vd(output22)=[1.5,0,2]
#uv(output22)=[0,-0.04,1]
#=========================================


vmin_vel = None#-0.25
vmax_vel = None#-vmin_vel

vmin_dens = None
vmax_dens = None

color_sink_colmap = 'firebrick'
transparence_sink_colmap = 0.4
color_sink_velmap = 'limegreen'
transparence_sink_velmap = 0.7

reposition_fig = False#Pour repositionner les figures ouvertes par matplotlib

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



#--------------------
#Definition du masque
#--------------------
#mask_rho_disk, mask = ed.pdf_to_singledisc_cells(path=path, num=num_output)
#mask_tot=np.copy(mask_rho_disk)
#mask_tot[mask_tot]=mask
rho_min_disk, pos_max_x, pos_max_y, pos_max_z, vel_max_x, vel_max_y, vel_max_z=ed.pdf_to_singledisc_for_plot(path=path, num=num_output)



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

rho_max_codeunits=np.max(rho)

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
    zoom_v=[0.045/2, 0.015/2, 0.005/2, 0.005/3./2, 0.5, 0.5/2, 0.005/2/1.5]
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



#-------------------------------------------
#-------------------------------------------
#Fonction qui calcul le radio entre Vrad et 
#Vpar comme dans "pdf_to_singledisc_cells" 
#de extract_disk.py
#-------------------------------------------
#-------------------------------------------
def compute_vel_ratio(pos_max_x,pos_max_y,pos_max_z,vel_max_x,vel_max_y,vel_max_z,pos,vel):
    #position of disk cell
    pos_disk_x = pos[:,0] - pos_max_x
    pos_disk_y = pos[:,1] - pos_max_y
    pos_disk_z = pos[:,2] - pos_max_z

    #position of disk cell
    vel_disk_x = vel[:,0] - vel_max_x
    vel_disk_y = vel[:,1] - vel_max_y
    vel_disk_z = vel[:,2] - vel_max_z

    #radial component of V
    norm_pos = np.sqrt(pos_disk_x**2+pos_disk_y**2+pos_disk_z**2)
    mask = norm_pos == 0.
    norm_pos[mask] = 1.e-10
    Vrad = (vel_disk_x*pos_disk_x + vel_disk_y*pos_disk_y + vel_disk_z*pos_disk_z) / norm_pos

    #non radial component of V
    Vpar_x = (vel_disk_z*pos_disk_y - vel_disk_y*pos_disk_z) / norm_pos
    Vpar_y = (vel_disk_x*pos_disk_z - vel_disk_z*pos_disk_x) / norm_pos
    Vpar_z = (vel_disk_y*pos_disk_x - vel_disk_x*pos_disk_y) / norm_pos
    Vpar = np.sqrt(Vpar_x**2+Vpar_y**2+Vpar_z**2)
    Vpar[mask] = 1.e-10

    #Ratio
    vrad_orth_ratio = -Vrad / Vpar

    return vrad_orth_ratio



#-------------------------------------------
#Definition de la fonction (densite masquee)
#-------------------------------------------
def disk_op(dset):
    rho = dset["rho"]
    vel = dset["vel"]
    pos = dset.get_cell_centers()
    flag = np.zeros_like(rho)
    mask = rho > rho_min_disk
    vrad_orth_ratio = compute_vel_ratio(pos_max_x,pos_max_y,pos_max_z,vel_max_x,vel_max_y,vel_max_z,pos[mask],vel[mask])
    #mask *= vrad_orth_ratio < 0.5
    mask2 = vrad_orth_ratio < 0.5
    mask[mask]=mask2
    flag[mask] = 1.0
    flag[~mask] = 0.0 
    return flag
#-------------------------------------------
#-------------------------------------------



#-----------------
#-----------------
#Calcul des cartes
#-----------------
#-----------------
def carte_dens_vel(view_diagram,up_vector,pos_colorbar_dens,pos_colorbar_vel,tag_compl):
    #-----------------------------------------------------------
    #Calcul de la carte ou l'on regarde suivant un vecteur donne
    #-----------------------------------------------------------
    cam_x = Camera(center=center,line_of_sight_axis=view_diagram,region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector=up_vector,map_max_size=512)

    #rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_density"])
    rho_op=ScalarOperator(disk_op, ro.info["unit_density"])
    rt = raytracing.RayTracer(amr,ro.info,rho_op)
    datamap = rt.process(cam_x)#, surf_qty=True)
    #map_col = np.log10(datamap.map.T*lbox_cm)
    map_col = datamap.map.T

    if v_proj == True:
        Vx_op = ScalarOperator(lambda dset:  (  (dset["vel"][...,0]*np.array(view_diagram)[0]) + (dset["vel"][...,1]*np.array(view_diagram)[1]) + (dset["vel"][...,2]*np.array(view_diagram)[2]) )   /np.sqrt(np.sum(np.array(view_diagram)**2))  *dset["rho"] ,  ro.info["unit_velocity"])
        rt = raytracing.RayTracer(amr,ro.info,Vx_op)
        datamap_vx = rt.process(cam_x, surf_qty=True)
        map_Vx = datamap_vx.map.T / datamap.map.T * factor_vel_km_s


    #--------------------
    #Debut figure coldens
    #--------------------
    plt.figure()
    ax=plt.gca()
    im1=plt.imshow(map_col,extent=[(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au,(-radius+center[2])*lbox_au,(radius+center[2])*lbox_au],origin='lower', vmin=vmin_dens, vmax=vmax_dens)   
    if pos_colorbar_dens=='left':
        ax.yaxis.set_ticks_position('right')
        ax.yaxis.set_label_position("right")

    if sink_plot==True:
        for i in range(len(m_sinks)):
            if (y_sinks[i]>(-radius+center[1])*lbox_au)&(y_sinks[i]<(radius+center[1])*lbox_au)&(z_sinks[i]>(-radius+center[2])*lbox_au)&(z_sinks[i]<(radius+center[2])*lbox_au):
                plt.plot(y_sinks[i],z_sinks[i],'.',color=color_sink_colmap,markersize=size_sinks[i],alpha=transparence_sink_colmap)
    ax.set_xlabel(r'$y$ (AU)')     
    ax.set_ylabel(r'$z$ (AU)')
    plt.xticks(rotation=90)
    plt.locator_params(axis='x', nbins=7)

    #if title_time==True:
        #plt.title('Time = '+str(int(simulation_time))+' years')
    if title_time_cor==True:
        plt.tight_layout(pad=0.1) #pad en inch si besoin
        plt.title('Time = '+str(int(simulation_time - ref[1]*1e6))+' years')

    # Adding the colorbar
    divider1 = make_axes_locatable(ax)
    if pos_colorbar_dens=='left':
        cax = divider1.append_axes("left", size="6%", pad=0.15)
        cbar1 = plt.colorbar(im1, cax=cax)#, **kw)
        cbar1.ax.yaxis.set_ticks_position('left')
        cbar1.ax.yaxis.set_label_position('left')
    else:
        cax = divider1.append_axes("right", size="6%", pad=0.15)
        cbar1 = plt.colorbar(im1, cax=cax)#, **kw)

    cbar1.set_label(r'$\text{log} \left( N \right) \, \, \left( cm^{-2} \right)$')
    if radius_zoom==5:
        plt.xlim([0,lbox_au])
        plt.ylim([0,lbox_au])
    if save==True:
        plt.tight_layout(pad=0.1) #pad en inch si besoin
        plt.savefig(path_save+simu+'_maskdisk_'+tag_compl+str(radius_zoom)+'_'+str(num_output)+'.pdf')#, bbox_inches='tight')







    #--------------------
    #Debut figure velproj
    #--------------------
    if v_proj == True:
        plt.figure()
        plt.xticks(rotation=90)
        plt.locator_params(axis='x', nbins=7)
        ax=plt.gca()
        norm = MidpointNormalize(midpoint=0)  #Pour avoir le centre de la colormap a 0
        im2=plt.imshow(map_Vx,extent=[(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au,(-radius+center[2])*lbox_au,(radius+center[2])*lbox_au],origin='lower',cmap='RdBu_r',norm=norm, vmin=vmin_vel, vmax=vmax_vel)
        if pos_colorbar_vel=='left':
            ax.yaxis.set_ticks_position('right')
            ax.yaxis.set_label_position("right")

        if sink_plot==True:
            for i in range(len(m_sinks)):
                if (y_sinks[i]>(-radius+center[1])*lbox_au)&(y_sinks[i]<(radius+center[1])*lbox_au)&(z_sinks[i]>(-radius+center[2])*lbox_au)&(z_sinks[i]<(radius+center[2])*lbox_au):
                    plt.plot(y_sinks[i],z_sinks[i],'.',color=color_sink_velmap,markersize=size_sinks[i],alpha=transparence_sink_velmap)
        plt.xlabel(r'$y$ (AU)')     
        plt.ylabel(r'$z$ (AU)')
        
        #if title_time==True:
            #plt.title('Time = '+str(int(simulation_time))+' years')
        if title_time_cor==True:
            plt.tight_layout(pad=0.1) #pad en inch si besoin
            plt.title('Time = '+str(int(simulation_time - ref[1]*1e6))+' years')

        # Adding the colorbar
        divider2 = make_axes_locatable(ax)
        if pos_colorbar_vel=='left':
            cax = divider2.append_axes("left", size="6%", pad=0.15)
            cbar2 = plt.colorbar(im2, cax=cax)#, **kw)
            cbar2.ax.yaxis.set_ticks_position('left')
            cbar2.ax.yaxis.set_label_position('left')
        else:
            cax = divider2.append_axes("right", size="6%", pad=0.15)
            cbar2 = plt.colorbar(im2, cax=cax)#, **kw)

        cbar2.set_label(r'$v_x \, \left( km.s^{-1} \right) $')  
        if radius_zoom==5:
            plt.xlim([0,lbox_au])
            plt.ylim([0,lbox_au])
        if save==True: 
            plt.tight_layout(pad=0.1) #pad en inch si besoin
            plt.savefig(path_save+simu+'_vel_'+tag_compl+str(radius_zoom)+'_'+str(num_output)+'.pdf')#, bbox_inches='tight')

    return




#=================================================
#Trace des cartes pour les vues edge_on et face-on 
#=================================================
#Edge-on:
if edge_on==True:
    pos_colorbar_dens='right'
    pos_colorbar_vel='left'
    tag_compl='edgeon_'
    carte_dens_vel(view_diagram_edgeon,up_vector_edgeon,pos_colorbar_dens,pos_colorbar_vel,tag_compl)
#Face-on:
if face_on==True:
    pos_colorbar_dens='left'
    pos_colorbar_vel='right'
    tag_compl='faceon_'
    carte_dens_vel(view_diagram_faceon,up_vector_faceon,pos_colorbar_dens,pos_colorbar_vel,tag_compl)
#=================================================
#=================================================




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



