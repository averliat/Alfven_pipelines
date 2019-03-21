# -*- coding: utf-8 -*-
import numpy as np
import extract_disk as ed
import module_extract as me

import pymses
import sys
import glob as glob
import os 

from pymses.filters import CellsToPoints
from pymses.utils import constants as cst
import matplotlib.pyplot as plt
plt.ion()

import h5py

import os
import pipeline_temps_0_simulation as t_0




#-------------------------------------------------------
#Entrer le nom de la simulation et le numero de l'output    
#-------------------------------------------------------
simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire_shr_bigbox_50pourc'
owner = 'averliat_alfven'

output = 59

seuil_rho = 1e-10


#Caracteristiques des coquilles en AU:
R_min = 1  #Devrait etre egale a dr en dessous
R_max = 'all_disk'  #='all_disk' pour prendre tout le disque
dr = 1


#Sauvegarde des quantites finales:
save = True
file_save = 'Structure_dominee_par_rotation_output'+str(output)+'_dr'+str(dr)+'AU.hdf5'




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

path_analyses='/home/averliat/these/analyses/'+simu+'/'



#-------------------
#Lecture de l'output
#-------------------
ro=pymses.RamsesOutput(path,output)
lbox_pc = ro.info['unit_length'].express(cst.pc)

amr = ro.amr_source(["rho","vel","P","phi","g"])

cell_source = CellsToPoints(amr)
cells = cell_source.flatten(verbose=False)
dx = cells.get_sizes()

pos = cells.points
vel = cells["vel"]
rho = cells["rho"]



#-------------------
#Lecture de l'output
#-------------------
if 'sink' in simu:
    sink = me.read_sink_cvs(output,path)[0]
    Msink = sink[1]  #en Masse solaire


#------------------------------
#Conversion en unites physiques
#------------------------------
lbox_au = ro.info['unit_length'].express(cst.au)
factor_dist= ro.info['unit_length'].express(cst.m)
factor_vel = ro.info['unit_velocity'].express(cst.m/cst.s)
factor_rho = ro.info['unit_density'].express(cst.kg/(cst.m)**3)
factor_time_Myr = ro.info['unit_time'].express(cst.Myr)
factor_vel_km_s = ro.info['unit_velocity'].express(cst.km/cst.s)

simulation_time = ro.info['time']*factor_time_Myr


dx *= factor_dist
pos *= factor_dist
vel *= factor_vel
rho *= factor_rho



#-------------------------
#Definition de la coquille
#-------------------------
lbox_m = factor_dist
#center = np.array(([lbox_m/2., lbox_m/2., lbox_m/2.]))

arg_centre = np.argmax(rho)
pos_center = [pos[:,0][arg_centre],pos[:,1][arg_centre],pos[:,2][arg_centre]]
vel_center = [vel[:,0][arg_centre],vel[:,1][arg_centre],vel[:,2][arg_centre]]



#--------------------
#Definition du masque
#--------------------
mask_rho_disk, mask = ed.pdf_to_singledisc_cells(path=path, num=output)

if (np.sum(mask)<=1):
    sys.exit("Erreur masque de selection du disque vide")
    


rho = rho[mask_rho_disk][mask]
dx = dx[mask_rho_disk][mask]
vel = vel[mask_rho_disk][mask]
pos = pos[mask_rho_disk][mask]



#------------------
#Rayon des cellules
#------------------
radius_cell = np.sqrt( np.sum( (pos[:,:]-pos_center)**2 , axis=1 ) )  #Tableau contenant le rayon de chaque cellule
pos_rel = pos[:,:] - pos_center #vecteur centre --> point
vel_norm = np.sqrt( np.sum( vel[:,:]**2 , axis=1 ) )  #Tableau contenant la norme de la vitesse de chaque cellule



#------------------------
#Definition des coquilles
#------------------------
if R_max=='all_disk':
    R_max=np.max(radius_cell)/factor_dist*lbox_au  #rayon max des cellules du disque en AU
nbr_shell = int((R_max-R_min)/dr)  #nombre de coquilles
shells_au=np.linspace(R_min,R_max,nbr_shell)  #rayon des coquilles en AU
#shells_au=np.logspace(np.log10(R_min),np.log10(R_max),nbr_shell)  #rayon des coquilles en AU
shells=shells_au/lbox_au  #rayon des coquilles en unite de code (correspondant au zoom_v ci-dessous)
shells_m=shells*factor_dist  #rayon des coquilles en m

#nbr_shell = 5
#shells_au=np.array([50.,300.,800.,2000.,6000.])
#shells=shells_au/lbox_au



#------------------------------------
#Debut de la boucle sur les coquilles
#------------------------------------
sum_m_int = 0.
radius_int = 0.

V_r_shell_moy_tab = np.zeros(nbr_shell)
V_r_rms_shell_moy_tab = np.zeros(nbr_shell)
V_paral_shell_moy_tab = np.zeros(nbr_shell)
sigma_shell_moy_tab = np.zeros(nbr_shell)
M_cumul_tab = np.zeros(nbr_shell)
M_int_tab = np.zeros(nbr_shell)

count=0
ind=0
verif=0
M_cumul=0

for radius in shells_m:
    if count==10:
        print radius/factor_dist*lbox_au
        count=0
    count += 1

    mask_shell = (radius_cell > radius_int) & (radius_cell <= radius)

    #Quantites dans la coquille
    shell_pos = pos[mask_shell]
    shell_vel = vel[mask_shell]
    shell_rho = rho[mask_shell]
    shell_dx = dx[mask_shell]
    shell_radius_cell = radius_cell[mask_shell]
    shell_pos_rel = pos_rel[mask_shell]
    shell_vel_norm = vel_norm[mask_shell]

    #Masse contenue dans la coquille
    sum_m_int = np.sum( shell_rho* shell_dx**3 ) #mettre += pour prendre en compte l'interieur
    M_cumul += sum_m_int
    
    #V_r moyen par shell
    V_r_shell = np.sum(shell_vel*shell_pos_rel,axis=1) / shell_radius_cell
    V_r_shell_moy = np.sum(V_r_shell*shell_rho*shell_dx**3)/sum_m_int
    V_r_rms_shell_moy = np.sqrt( np.sum(V_r_shell**2*shell_rho*shell_dx**3)/sum_m_int )

    #V_parallele moyen par shell
    V_paral_shell = np.sqrt( shell_vel_norm**2 - V_r_shell**2 )
    V_paral_shell_moy = np.sum(V_paral_shell*shell_rho*shell_dx**3)/sum_m_int 

    #sigma = sum_sur_shell( (v_vecteur - V_r_shell_moy.er)**2 ) pondere par masse
    sigma_shell = np.sqrt( (V_r_shell-V_r_shell_moy)**2 + V_paral_shell**2 )
    sigma_shell_moy = np.sum(sigma_shell*shell_rho*shell_dx**3)/sum_m_int

    #Sauvegarde dans les tableaux
    V_r_shell_moy_tab[ind] = V_r_shell_moy
    V_r_rms_shell_moy_tab[ind] = V_r_rms_shell_moy
    V_paral_shell_moy_tab[ind] = V_paral_shell_moy
    sigma_shell_moy_tab[ind] = sigma_shell_moy
    M_cumul_tab[ind] = M_cumul
    M_int_tab[ind] = sum_m_int


    radius_int = radius
    ind += 1
    verif=1




#-----------------------------------------------------------------------
#Sauvegarde du tableau contenant les moments cinetiques de chaque output
#-----------------------------------------------------------------------
if (save == True and verif==1):
    h5f = h5py.File(path_analyses+file_save, 'w')
    h5f.create_dataset('shells_au', data=shells_au)
    h5f.create_dataset('V_r_shell_moy_tab', data=V_r_shell_moy_tab)
    h5f.create_dataset('V_r_rms_shell_moy_tab', data=V_r_rms_shell_moy_tab)
    h5f.create_dataset('V_paral_shell_moy_tab', data=V_paral_shell_moy_tab)
    h5f.create_dataset('sigma_shell_moy_tab', data=sigma_shell_moy_tab)
    h5f.create_dataset('M_cumul_tab',data=M_cumul_tab)
    h5f.create_dataset('M_int_tab',data=M_int_tab)
    h5f.create_dataset('simulation_time', data=simulation_time)
    #h5f.create_dataset('Msink',data=Msink)


    h5f.close()

plt.close('all')


