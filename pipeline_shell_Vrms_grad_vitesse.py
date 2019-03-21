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

import os




#-------------------------------------------------------
#Entrer le nom de la simulation et le numero de l'output    
#-------------------------------------------------------
simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire50_shr'

owner = 'averliat_alfven'
output = 685


dir_save = 'Gradient_vitesse'  #Repertoire de sauvegarde des figures

save = True


#Caracteristiques des coquilles en AU:
R_min = 50
R_max = 10000  #='all' pour prendre toute la boite
dr = 50


#Sauvegarde des quantites finales:
file_save = 'Vrms_gradient_vitesse_log_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.hdf5'
#file_save = 'Vrms_gradient_vitesse_log_boxref_output'+str(output)+'_Rmin'+str(R_min)+'_Rmax'+str(R_max)+'_dr'+str(dr)+'.hdf5'  #Si vitesse par rapport a la boite, et pas par raport au disque



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

#path_save='/home/averliat/these/analyses/'+simu+'/'+dir_save+'/'
path_analyse='/home/averliat/these/analyses/'+simu+'/'



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



#------------------------------
#Conversion en unites physiques
#------------------------------
lbox_au = ro.info['unit_length'].express(cst.au)
lbox_pc = ro.info['unit_length'].express(cst.pc)
factor_dist= ro.info['unit_length'].express(cst.m)
factor_vel = ro.info['unit_velocity'].express(cst.m/cst.s)
factor_rho = ro.info['unit_density'].express(cst.kg/(cst.m)**3)
factor_time_Myr = ro.info['unit_time'].express(cst.Myr)
factor_vel_km_s = ro.info['unit_velocity'].express(cst.km/cst.s)

simulation_time = ro.info['time']*factor_time_Myr


dx *= factor_dist
#pos *= factor_dist
vel *= factor_vel_km_s
rho *= factor_rho



#---------------------------------------------------------
#Definition du centre des images et de leur niveau de zoom
#---------------------------------------------------------
#Position du "centre" de la simulation = endroit le plus dense
arg_centre = np.argmax(rho)
center = [pos[:,0][arg_centre],pos[:,1][arg_centre],pos[:,2][arg_centre]]
vel_center = [vel[:,0][arg_centre],vel[:,1][arg_centre],vel[:,2][arg_centre]**2]

nbr_shell = int((R_max-R_min)/dr)  #nombre de coquilles
#shells_au=np.linspace(R_min,R_max,nbr_shell)  #rayon des coquilles en AU
shells_au=np.logspace(np.log10(R_min),np.log10(R_max),nbr_shell)  #rayon des coquilles en AU
shells=shells_au/lbox_au  #rayon des coquilles en unite de code (correspondant au zoom_v ci-dessous)
#nbr_shell = 5
#shells_au=np.array([50.,300.,800.,2000.,6000.])
#shells=shells_au/lbox_au



#------------------
#Rayon des cellules
#------------------
radius_cell = np.sqrt( np.sum( (pos[:,:]-center)**2 , axis=1 ) )  #Tableau contenant le rayon de chaque cellule



#-----------------------------------
#Initialisation des moments integres
#-----------------------------------
sum_v2_int = 0.
sum_m_int = 0.
sum_v_abs_int = 0.
radius_int = 0.

sum_v2_int_tab = np.zeros(nbr_shell)
sum_m_int_tab = np.zeros(nbr_shell)
sum_v_abs_int_tab = np.zeros(nbr_shell)
shells_pc_tab = np.array(shells*lbox_pc)



#------------------------------------
#Debut de la boucle sur les coquilles
#------------------------------------
count=0
ind=0
verif=0

for radius in shells:
    if count==10:
        print radius*lbox_au
        count=0
    count += 1

    mask = (radius_cell > radius_int) & (radius_cell <= radius)

    #Quantites dans la coquille
    shell_pos = pos[mask]
    shell_vel = vel[mask]
    shell_rho = rho[mask]
    shell_dx = dx[mask]
    shell_radius_cell = radius_cell[mask]

    #Masse contenue dans la coquille
    sum_m_int = np.sum( shell_rho* shell_dx**3 ) #mettre += pour prendre en compte l'interieur
    sum_v2_int = np.sum( (  np.sqrt( np.sum( (shell_vel[:,:]-vel_center)**2 , axis=1 ) )    )**2 *shell_rho* shell_dx**3 ) #pareil
    #sum_v2_int = np.sum( (  np.sqrt( np.sum( (shell_vel[:,:])**2 , axis=1 ) )    )**2 *shell_rho* shell_dx**3 ) #Vitesse par rapport a la boite, et pas par rapport au disque
    sum_v_abs_int = np.sum(  np.sqrt( np.sum( shell_vel[:,:]**2 , axis=1 ) )     *shell_rho* shell_dx**3 ) #pareil

    sum_v2_int_tab[ind]=sum_v2_int
    sum_m_int_tab[ind]=sum_m_int
    sum_v_abs_int_tab[ind]=sum_v_abs_int

    radius_int = radius
    ind += 1
    verif=1



#-----------------------------------
#Calcul du gradient a chaque echelle
#-----------------------------------
Vrms = np.sqrt( sum_v2_int_tab / sum_m_int_tab )
grad_v = Vrms / shells_pc_tab
v_abs_moy = sum_v_abs_int_tab / sum_m_int_tab



#-----------------------------------------------------------------------
#Sauvegarde du tableau contenant les moments cinetiques de chaque output
#-----------------------------------------------------------------------
if (save == True and verif==1):
    h5f = h5py.File(path_analyse+file_save, 'w')
    h5f.create_dataset('shells_au', data=shells_au)
    h5f.create_dataset('grad_v',data=grad_v)
    h5f.create_dataset('Vrms',data=Vrms)
    h5f.create_dataset('v_abs_moy',data=v_abs_moy)

    h5f.close()



