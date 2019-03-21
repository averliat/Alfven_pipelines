import numpy as np

import pymses
import sys
import glob as glob
import os

from pymses.filters import CellsToPoints
from pymses.utils import constants as cst
from pymses.analysis import Camera, raytracing, slicing
from pymses.analysis import ScalarOperator, FractionOperator, MaxLevelOperator
from pymses.analysis.cube3d import CubeExtractor

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
simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire50_shr'

owner = 'averliat_alfven'
num_output = 685

level_max = 8+10



dir_save = 'Kolmogorov_gradient_vitesse'  #Repertoire de sauvegarde des figures

save_j_turb = True  #Pour sauvegarder les tableaux contenant les angles des gradients


#Caracteristiques des coquilles en AU:
R_min = 50
R_max = 17000  #='all' pour prendre toute la boite
dr = 500


file_save_j_turb = 'J_turb_log_output'+str(num_output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.hdf5'


title_time=True
title_time_cor=False
seuil_rho = 1e-10




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


if save_j_turb==True:
    if os.path.isdir(path_save) == False:
        os.mkdir(path_save)



#-------------------
#Lecture de l'output
#-------------------
ro=pymses.RamsesOutput(path,num_output)
lbox=ro.info['boxlen']
lbox_pc = ro.info['unit_length'].express(cst.pc)
lbox_au = ro.info['unit_length'].express(cst.au)
factor_vel_km_s = ro.info['unit_velocity'].express(cst.km/cst.s)



amr = ro.amr_source(["rho","vel","P","phi","g"])
cell_source = CellsToPoints(amr)
cells = cell_source.flatten()
pos = cells.points
rho = cells["rho"]


#Centre du disque
arg_centre = np.argmax(rho)
center = [pos[:,0][arg_centre],pos[:,1][arg_centre],pos[:,2][arg_centre]]

radius = lbox * R_max / lbox_au
#center = [0.5,0.5,0.5]
#radius = lbox_pc*4
los='x'

if los=='x':
    upvect='z'
    v_ind=0
elif los=='y':
    upvect='x'
    v_ind=1
elif los=='z':
    upvect='y'
    v_ind=2


#level_max=8
map_size=2048
map_size=map_size/2*2+1



cam = Camera(center=center,line_of_sight_axis=los,region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector=upvect,map_max_size=map_size)


rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_density"])
rt = raytracing.RayTracer(amr,ro.info,rho_op)
datamap = rt.process(cam, surf_qty=True)
#map_col = np.log10(datamap.map.T*lbox_cm)


V_op = ScalarOperator(lambda dset: dset["vel"][...,v_ind]*dset["rho"] ,  ro.info["unit_velocity"])
rt = raytracing.RayTracer(amr,ro.info,V_op)
datamap_v = rt.process(cam, surf_qty=True)
map_V = datamap_v.map.T / datamap.map.T * factor_vel_km_s


dist=np.zeros(shape=(map_size,map_size))
for i in range(map_size):
    d_abs=(i-(map_size-1)/2)* R_max/((map_size-1)/2)
    for j in range(map_size):
        d_ord=(j-(map_size-1)/2)* R_max/((map_size-1)/2)
        dist[i,j]=np.sqrt(d_abs**2+d_ord**2)







nbr_shell = int((R_max-R_min)/dr)  #nombre de coquilles
#shells_au=np.linspace(R_min,R_max,nbr_shell)  #rayon des coquilles en AU
shells_au=np.logspace(np.log10(R_min),np.log10(R_max),nbr_shell)  #rayon des coquilles en AU
shells_pc=shells_au/lbox_au*lbox_pc
#shells=shells_au/lbox_au  #rayon des coquilles en unite de code (correspondant au zoom_v ci-dessous)
#nbr_shell = 5
#shells_au=np.array([50.,300.,800.,2000.,6000.])
#shells=shells_au/lbox_au





#=============================================================================
#=============================================================================
#Boucle sur les coquilles de rayons croissants
#=============================================================================
#=============================================================================

#Initialisation des tableaux de stockage
j_turb_shell_tab = np.zeros(nbr_shell)
j_turb_shell_tab_2 = np.zeros(nbr_shell)
radius_int = -0.1

count=0
ind=0
verif=0


#Debut de la boucle sur les coquilles
for radius in shells_au:
    if count==10:
        print radius
        count=0
    count += 1

    mask = (dist > radius_int) & (dist <= radius)

    #Quantites dans la coquille
    shell_V = map_V[mask]
    shell_dist = dist[mask]
    

    #Calcul de j_turb
    moy_delta_V_los = np.mean(np.abs(shell_V))
    j_turb_shell = moy_delta_V_los*radius/lbox_au*lbox_pc #km.s-1.pc

    j_turb_shell_2 = np.mean(np.abs(shell_V*shell_dist/lbox_au*lbox_pc))


    #Sauvegarde dans les tableaux
    j_turb_shell_tab[ind] = j_turb_shell
    j_turb_shell_tab_2[ind] = j_turb_shell_2


    radius_int = radius
    ind += 1
    verif=1


if (save_j_turb == True and verif==1):
    h5f = h5py.File(path_analyse+file_save_j_turb, 'w')

    h5f.create_dataset('shells_au', data=shells_au)
    h5f.create_dataset('j_turb_shell_tab', data=j_turb_shell_tab)
    h5f.create_dataset('j_turb_shell_tab_2', data=j_turb_shell_tab_2)

    h5f.close()



#V_op = ScalarOperator(lambda dset: dset["vel"] ,  ro.info["unit_velocity"])
#dcube = CubeExtractor(amr,V_op)
#vel = dcube.process(cam,cube_size=2.*radius,resolution=2**level_max)
#vel = vel.data

#rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_density"])
#dcube = CubeExtractor(amr,rho_op)
#rho = dcube.process(cam,cube_size=2.*radius,resolution=2**level_max)
#rho = rho.data

#radius_cell_op = ScalarOperator(lambda dset: np.sqrt((dset["pos"][...,0]-center[0])**2+(dset["pos"][...,1]-center[1])**2+(dset["pos"][...,2]-center[2])**2) ,  ro.info["unit_length"])
#dcube = CubeExtractor(amr,radius_cell_op)
#radius_cell = dcube.process(cam,cube_size=2.*radius,resolution=2**level_max)
#radius_cell = radius_cell.data

