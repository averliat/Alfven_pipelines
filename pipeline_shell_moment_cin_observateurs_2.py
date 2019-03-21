# -*- coding: utf-8 -*-
import numpy as np
import extract_disk as ed

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
simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire_shr_bigbox_60pourc'
owner = 'averliat_alfven'

output_min = 'None'
output_max = 97

seuil_rho = 1e-10


#rho_seuil = 1e11 #Hcc, pour la selection des cellules

#Sauvegarde des quantites finales:
save_tab = True
#file_save = 'Moment_cin_observateurs_dans_coquille_de_densite'+str(rho_seuil)+'_output'+str(output_min)+'_to_'+str(output_max)+'.hdf5'




#------------------------
#Construction des chemins   
#------------------------
#Chemin de la simulation et numeros des outputs
if owner=='averliat':
    path='/mnt/magmist/magmist/simu_B335_averliat/'+simu+'/'
    path_analyses='/home/averliat/these/analyses/'+simu+'/'
if owner=='phennebe':
    path='/drf/projets/capucine/phennebe/'+simu+'/'
    path_analyses='/home/averliat/these/analyses/'+simu+'/'
if owner=='sapanais':
    path='/dsm/anais/storageA/magmist/'+simu+'/'
    path_analyses='/home/averliat/these/analyses/'+simu+'/'
if owner=='averliat_alfven':
    path='/drf/projets/alfven-data/averliat/'+simu+'/'
    path_analyses='/home/averliat/these/analyses/'+simu+'/'



#---------------------------------------------------------------------
#Lecture du temps et de l'output initial recale si output_min = 'None'
#---------------------------------------------------------------------
if os.path.isfile(path_analyses+'t0_seuil_rho_'+str(seuil_rho)+'.txt') == False:
    ref = np.array(t_0.temps_0_simu(path, seuil_rho))
    np.savetxt(path_analyses+'t0_seuil_rho_'+str(seuil_rho)+'.txt', ref)
else:
    ref=np.loadtxt(path_analyses+'t0_seuil_rho_'+str(seuil_rho)+'.txt')

output_0 = int(ref[0])
time_0 = ref[1]


if output_min == 'None':
    output_min = output_0


nbre_output = output_max - output_min + 1



#-----------------------------------
#Creation des tableaux de sauvegarde
#-----------------------------------
if save_tab == True:
    file_save = 'Moment_cin_observateurs_dans_disque_par_pdf_output'+str(output_min)+'_to_'+str(output_max)+'.hdf5'



#----------------------------------------------------
#Initialisation des differents tableaux de sauvegarde	
#----------------------------------------------------
#norm_mom_cin_obs = np.zeros(nbre_output)
simulation_time = np.zeros(nbre_output)
norm_GC_AU = np.zeros(nbre_output)
norm_center_C_AU = np.zeros(nbre_output)
norm_center_G_AU = np.zeros(nbre_output)
vel_GC_tab = np.zeros(nbre_output)
norm_mom_cin_obs_par_integ = np.zeros(nbre_output)
numero_output = np.zeros(nbre_output)



#---------------------------------------------------------------------------------------
#Boucle sur les differents output pour extraire les moments des disques de chaque output   
#---------------------------------------------------------------------------------------
output_manquant=0
nbre_output_effectif=0
verif=0
for l in range(nbre_output):
    output=output_min + l
    try:
        ro=pymses.RamsesOutput(path,output)
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
    numero_output[nbre_output_effectif] = output
    #-------------------
    #Lecture de l'output
    #-------------------
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
    factor_dist= ro.info['unit_length'].express(cst.m)
    factor_vel = ro.info['unit_velocity'].express(cst.m/cst.s)
    factor_rho = ro.info['unit_density'].express(cst.kg/(cst.m)**3)
    factor_rho_Hcc = ro.info['unit_density'].express(cst.H_cc)
    factor_time_Myr = ro.info['unit_time'].express(cst.Myr)

    simulation_time[nbre_output_effectif] = ro.info['time']*factor_time_Myr


    dx *= factor_dist
    pos *= factor_dist
    vel *= factor_vel
    rho_Hcc = rho * factor_rho_Hcc
    rho *= factor_rho



    #-------------------------
    #Definition de la coquille
    #-------------------------
    AU = 149597870700 #m

    #Position des differents points
    lbox_m = factor_dist
    center = np.array(([lbox_m/2., lbox_m/2., lbox_m/2.]))

    arg_centre = np.argmax(rho)
    pos_C = [pos[:,0][arg_centre],pos[:,1][arg_centre],pos[:,2][arg_centre]]
    vel_C = [vel[:,0][arg_centre],vel[:,1][arg_centre],vel[:,2][arg_centre]]


    mass_tot = np.sum( rho * dx**3 )
    coord_G_x = 1./mass_tot * np.sum( rho * dx**3 * pos[:,0] )
    coord_G_y = 1./mass_tot * np.sum( rho * dx**3 * pos[:,1] )
    coord_G_z = 1./mass_tot * np.sum( rho * dx**3 * pos[:,2] )
    pos_G = np.array([coord_G_x,coord_G_y,coord_G_z])


    norm_GC_AU[nbre_output_effectif] = np.sqrt( np.sum( ((pos_G - pos_C)/AU)**2 ) )
    norm_center_C_AU[nbre_output_effectif] = np.sqrt( np.sum( ((center - pos_C)/AU)**2 ) )
    norm_center_G_AU[nbre_output_effectif] = np.sqrt( np.sum( ((center - pos_G)/AU)**2 ) )


    print
    print
    print('==========================================================================')
    print("Difference barycentre - point le plus dense = "+str(norm_GC_AU[nbre_output_effectif])+' AU')
    print('==========================================================================')
    print

    print
    print('==========================================================================')
    print("Difference barycentre - centre de la boite = "+str(norm_center_G_AU[nbre_output_effectif])+' AU')
    print('==========================================================================')
    print

    print
    print('==========================================================================')
    print("Difference centre de la boite - point le plus dense = "+str(norm_center_C_AU[nbre_output_effectif])+' AU')
    print('==========================================================================')
    print
    print   



    #Recuperation de l'indice du barycentre pour pouvoir recuperer sa vitesse
    #radius_cell_from_G = np.sqrt( np.sum( (pos[:,:]-pos_G)**2 , axis=1 ) )  #Tableau contenant la distance de chaque cellule par rapport au barycentre
    #arg_G = np.argmin(radius_cell_from_G)
    #pos_G_cor = [pos[:,0][arg_G],pos[:,1][arg_G],pos[:,2][arg_G]]
    #vel_G_cor = [vel[:,0][arg_G],vel[:,1][arg_G],vel[:,2][arg_G]]
    #norm_GGcor_AU = np.sqrt( np.sum( ((pos_G - pos_G_cor)/AU)**2 ) )

    #print
    #print('==========================================================================')
    #print("Difference barycentre calcule - barycentre cellule = "+str(norm_GGcor_AU)+' AU')
    #print('==========================================================================')
    #print
    #print



    #--------------------
    #Definition du masque
    #--------------------
    #mask0 = rho_Hcc >= rho_seuil
    mask_rho_disk, mask = ed.pdf_to_singledisc_cells(path=path, num=output)
    
    if (np.sum(mask)<=1):
        norm_mom_cin_obs_par_integ[nbre_output_effectif]
        nbre_output_effectif += 1
        continue





    rho = rho[mask_rho_disk][mask]
    dx = dx[mask_rho_disk][mask]
    vel = vel[mask_rho_disk][mask]
    pos = pos[mask_rho_disk][mask]




    #Calcul du moment cinetique observateurs
    radius_cell_C = np.sqrt( np.sum( (pos[:,:]-pos_C)**2 , axis=1 ) )  #Tableau contenant le rayon de chaque cellule par rapport a C
    pos_cell_from_C = pos[:,:]-pos_C
    vel_cell_from_C = vel[:,:]-vel_C

    '''
    GC = pos_C - pos_G
    vel_G_C = vel_C

    vel_GC_tab[nbre_output_effectif] = np.sqrt( np.sum( np.array(vel_G_C)**2 ))

    #Termes du produit vectoriel du rayon par la vitesse
    premier_terme = GC[1]*vel_G_C[2] - GC[2]*vel_G_C[1]
    deuxieme_terme = GC[2]*vel_G_C[0] - GC[0]*vel_G_C[2]
    troisieme_terme = GC[0]*vel_G_C[1] - GC[1]*vel_G_C[0]

    #Norme du moment cinetique observateur
    norm_mom_cin_obs[nbre_output_effectif] = mass_tot * np.sqrt( premier_terme**2 + deuxieme_terme**2 + troisieme_terme**2 )
    '''







    #------------------------------------------------------------------
    #Calcul du moment cinetique de chaque cellule par rapport au centre
    #------------------------------------------------------------------ 
    moment_x = 0
    moment_y = 0
    moment_z = 0

    for i in range(len(dx)):  #Boucle sur les differentes cellules dans la coquille
        #Termes du produit vectoriel du rayon par la vitesse
        premier_terme = pos_cell_from_C[i,1]*vel_cell_from_C[i,2] - pos_cell_from_C[i,2]*vel_cell_from_C[i,1]
        deuxieme_terme = pos_cell_from_C[i,2]*vel_cell_from_C[i,0] - pos_cell_from_C[i,0]*vel_cell_from_C[i,2]
        troisieme_terme = pos_cell_from_C[i,0]*vel_cell_from_C[i,1] - pos_cell_from_C[i,1]*vel_cell_from_C[i,0]

        #Composantes du moment cinetique total
        moment_x += rho[i] * (dx[i]**3) * premier_terme
        moment_y += rho[i] * (dx[i]**3) * deuxieme_terme
        moment_z += rho[i] * (dx[i]**3) * troisieme_terme


    #--------------------------------------------------
    #Moment cinetique par rapport au centre et stockage
    #--------------------------------------------------
    norm_mom_cin_obs_par_integ[nbre_output_effectif] = np.sqrt( moment_x**2 + moment_y**2 + moment_z**2 )

    nbre_output_effectif += 1



#norm_mom_cin_obs = norm_mom_cin_obs[0:nbre_output_effectif]
simulation_time = simulation_time[0:nbre_output_effectif]
norm_GC_AU = norm_GC_AU[0:nbre_output_effectif]
norm_center_C_AU = norm_center_C_AU[0:nbre_output_effectif]
norm_center_G_AU = norm_center_G_AU[0:nbre_output_effectif]
vel_GC_tab = vel_GC_tab[0:nbre_output_effectif]
norm_mom_cin_obs_par_integ = norm_mom_cin_obs_par_integ[0:nbre_output_effectif]
numero_output = numero_output[0:nbre_output_effectif]




#-----------------------------------------------------------------------
#Sauvegarde du tableau contenant les moments cinetiques de chaque output
#-----------------------------------------------------------------------
if (save_tab == True and verif==1):
    h5f = h5py.File(path_analyses+file_save, 'w')

    #h5f.create_dataset('norm_mom_cin_obs', data=norm_mom_cin_obs)
    h5f.create_dataset('simulation_time',data=simulation_time)
    h5f.create_dataset('norm_GC_AU',data=norm_GC_AU)
    h5f.create_dataset('norm_center_C_AU',data=norm_center_C_AU)
    h5f.create_dataset('norm_center_G_AU',data=norm_center_G_AU)
    h5f.create_dataset('vel_GC_tab',data=vel_GC_tab)
    h5f.create_dataset('norm_mom_cin_obs_par_integ',data=norm_mom_cin_obs_par_integ)
    h5f.create_dataset('numero_output',data=numero_output)

    h5f.close()




