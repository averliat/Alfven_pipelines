import numpy as np

import pymses
import sys
import glob as glob
import os 

from pymses.filters import CellsToPoints
from pymses.utils import constants as cst
import matplotlib.pyplot as plt
plt.ion()

from astropy.io import fits
import h5py
import os



#------------------------------------------
#Entree pour la simulation et la sauvegarde
#------------------------------------------
#Chemin de la simulation et numeros des outputs
simu='B335_noturb_norot_hydro_pert_asym_aleatoire50_sshr_hugebox'
owner = 'averliat_alfven'

num_min = 1
num_max = 102


#Pour la sauvegarde
save_final = True



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

dir_save = '/home/averliat/these/analyses/'+simu+'/'
file_save = 'Moment_cinetique_tot_et_abs_output_'+str(num_min)+'_'+str(num_max)+'.hdf5'




#-----------------------------------------------------
#Initialisation des tableaux de stockage des resultats
#-----------------------------------------------------
num_tot = num_max-num_min+1

J_tot_tab = np.zeros(num_tot)
J_abs_tot_tab = np.zeros(num_tot)
numero_output = np.zeros(num_tot)
time_tab = np.zeros(num_tot)



#---------------------------------------------
#Debut de la boucle sur les differents outputs
#---------------------------------------------
output_manquant=0
nbre_output_effectif=0
verif=0
for l in range(num_tot):
    num = num_min + l
    try:
        ro=pymses.RamsesOutput(path,num)
    except:
        output_manquant += 1
        if output_manquant==100:
            print
            print
            print("=========================================================================================")
            print("=========================================================================================")
            print("Pas assez d'output pour la simulation : "+path )
            print("=========================================================================================")
            print("=========================================================================================")
            break
        continue

    verif=1
    output_manquant=0

    #-------------------
    #Lecture des donnees
    #-------------------
    ro=pymses.RamsesOutput(path,num)
    #lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)
    lbox_pc = ro.info['unit_length'].express(cst.pc)

    amr = ro.amr_source(["rho","vel","P","phi","g"])

    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    dx = cells.get_sizes()

    pos = cells.points
    vel = cells["vel"]
    rho = cells["rho"]

    centre = np.array([0.5,0.5,0.5])



    #------------------------------
    #Conversion en unites physiques
    #------------------------------
    #lbox_m = ro.info['unit_length'].express(cst.m)
    factor_dist= ro.info['unit_length'].express(cst.m)
    factor_vel = ro.info['unit_velocity'].express(cst.m/cst.s)
    factor_rho = ro.info['unit_density'].express(cst.kg/(cst.m)**3)
    factor_time_Myr = ro.info['unit_time'].express(cst.Myr)

    simulation_time = ro.info['time']*factor_time_Myr

    dx *= factor_dist
    pos *= factor_dist
    centre *= factor_dist
    vel *= factor_vel
    rho *= factor_rho

    pos_moins_centre = pos - centre


    #------------------------------------------------------------------
    #Calcul du moment cinetique de chaque cellule par rapport au centre
    #------------------------------------------------------------------
    #radius_cell_tab = np.sqrt( np.sum( (pos[:,:]-centre)**2 , axis=1 ) )  #Tableau contenant le rayon de chaque cellule

    J_tot_abs = 0
    moment_tot_x = 0
    moment_tot_y = 0
    moment_tot_z = 0

    for i in range(len(dx)):  #Boucle sur les differentes cellules
        #Termes du produit vectoriel du rayon par la vitesse
        premier_terme = pos_moins_centre[i,1]*vel[i,2] - pos_moins_centre[i,2]*vel[i,1]
        deuxieme_terme = pos_moins_centre[i,2]*vel[i,0] - pos_moins_centre[i,0]*vel[i,2]
        troisieme_terme = pos_moins_centre[i,0]*vel[i,1] - pos_moins_centre[i,1]*vel[i,0]

        #Composantes du moment cinetique total
        moment_tot_x += rho[i] * (dx[i]**3) * premier_terme
        moment_tot_y += rho[i] * (dx[i]**3) * deuxieme_terme
        moment_tot_z += rho[i] * (dx[i]**3) * troisieme_terme
        
        #Norme du produit vectoriel
        norm_rad_vect_vel = np.sqrt( premier_terme**2 + deuxieme_terme**2 + troisieme_terme**2 )

        #Moment cinetique absolue de la cellule consideree
        J_tot_abs += rho[i] * (dx[i]**3) * norm_rad_vect_vel



    #--------------------------------------------------
    #Moment cinetique par rapport au centre et stockage
    #--------------------------------------------------
    J_tot = np.sqrt( moment_tot_x**2 + moment_tot_y**2 + moment_tot_z**2 )
    J_tot_tab[nbre_output_effectif] = J_tot

    #-----------------------------------------------------------------------
    #Valeur absolue des moments cinetiques par rapport au centre et stockage
    #-----------------------------------------------------------------------
    J_abs_tot_tab[nbre_output_effectif] = J_tot_abs


    numero_output[nbre_output_effectif] = num
    time_tab[nbre_output_effectif] = simulation_time

    nbre_output_effectif += 1



#-----------------------------------------------------------------------
#Sauvegarde du tableau contenant les moments cinetiques de chaque output
#-----------------------------------------------------------------------
J_tot_tab = J_tot_tab[0:nbre_output_effectif]
J_abs_tot_tab = J_abs_tot_tab[0:nbre_output_effectif]
numero_output = numero_output[0:nbre_output_effectif]
time_tab = time_tab[0:nbre_output_effectif]



#-----------------------------------------------------------------------
#Sauvegarde du tableau contenant les moments cinetiques de chaque output
#-----------------------------------------------------------------------
if (save_final == True and verif==1):
    h5f = h5py.File(dir_save+file_save,'w') 

    h5f.create_dataset('J_tot_tab', data=J_tot_tab)
    h5f.create_dataset('J_abs_tot_tab', data=J_abs_tot_tab)
    h5f.create_dataset('numero_output', data=numero_output)
    h5f.create_dataset('time_tab',data=time_tab)

    h5f.close()




#print("Avancement : ",round((n+1)/N_cluster*100.,2),"%", end="\r")
