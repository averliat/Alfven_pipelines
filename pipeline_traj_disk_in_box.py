# -*- coding: utf-8 -*-
import os 
import numpy as np
import matplotlib.pyplot as plt
import h5py
import pymses
from pymses.filters import CellsToPoints
from pymses.utils import constants as cst
plt.ion()
import pipeline_temps_0_simulation as t_0



#-------------------------------------------------------
#Entrer le nom de la simulation et le numero de l'output	
#-------------------------------------------------------
simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire_shr_bigbox_50pourc'
owner = 'averliat_alfven'

output_min = 'None'
output_max = 110

seuil_rho = 1e-10

save_tab = True



#------------------------
#Construction des chemins 	
#------------------------
#Chemin de la simulation et numeros des outputs
if owner == 'averliat_alfven':
    path='/drf/projets/alfven-data/averliat/'+simu+'/'
if owner=='sapanais':
	path='/dsm/anais/storageA/magmist/'+simu+'/'
if owner=='phennebe':
    path='/drf/projets/capucine/phennebe/'+simu+'/'
    
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


nmbr_output = output_max - output_min + 1



#-----------------------------------
#Creation des tableaux de sauvegarde
#-----------------------------------
if save_tab == True:
	file_save = 'Traj_disk_in_box_output_'+str(output_min)+'_'+str(output_max)+'.hdf5'



#----------------------------------------------------
#Initialisation des differents tableaux de sauvegarde	
#----------------------------------------------------
numero_output = np.zeros(nmbr_output)
time_tab = np.zeros(nmbr_output)
time_cor_tab = np.zeros(nmbr_output)
max_rho_tab = np.zeros((nmbr_output,3))


#------------------------------------------------------------------------------------
#Boucle sur les differents output pour extraire le rayon des disques de chaque output	
#------------------------------------------------------------------------------------
output_manquant=0
nbre_output_effectif=0
verif=0
for i in range(nmbr_output):
    num = i+output_min
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
    #Lecture de l'output
    #-------------------
    ro=pymses.RamsesOutput(path,num)
    #lbox_pc = ro.info['unit_length'].express(cst.pc)
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
    factor_time_Myr = ro.info['unit_time'].express(cst.Myr)
    factor_vel_km_s = ro.info['unit_velocity'].express(cst.km/cst.s)

    simulation_time = ro.info['time']*factor_time_Myr

    arg_centre = np.argmax(rho)
    center = [pos[:,0][arg_centre],pos[:,1][arg_centre],pos[:,2][arg_centre]]

    numero_output[nbre_output_effectif] = num
    time_tab[nbre_output_effectif] = simulation_time
    time_cor_tab[nbre_output_effectif] = time_tab[nbre_output_effectif]-time_0
    max_rho_tab[nbre_output_effectif] = center

    nbre_output_effectif += 1



#-----------------------------------------------------------------------
#Sauvegarde du tableau contenant les moments cinetiques de chaque output
#-----------------------------------------------------------------------
time_tab = time_tab[0:nbre_output_effectif]
time_cor_tab = time_cor_tab[0:nbre_output_effectif]
numero_output = numero_output[0:nbre_output_effectif]
max_rho_tab = max_rho_tab[0:nbre_output_effectif]



if (save_tab == True and verif==1):
    h5f = h5py.File(path_analyses+file_save,'w') 

    h5f.create_dataset('time_tab', data=time_tab)
    h5f.create_dataset('time_cor_tab', data=time_cor_tab)
    h5f.create_dataset('numero_output', data=numero_output)
    h5f.create_dataset('max_rho_tab', data=max_rho_tab)


    h5f.close()




'''
#-------------------------------
#Trace des differentes quantites	
#-------------------------------

plt.plot(time_tab,mean_rad_loc_tab, marker='.',color='midnightblue', label=ur'$R_{disk}$')
plt.xlabel(ur'Temps ($Myr$)')
plt.ylabel(ur'Rayon moyen du disque ($AU$)')
plt.legend(loc='best')

plt.figure()
plt.plot(time_tab-time_0,mean_rad_loc_tab, marker='.',color='chocolate', label=ur'$R_{disk}$')
plt.xlabel(ur'Temps corrigé ($Myr$)')
plt.ylabel(ur'Rayon moyen du disque ($AU$)')
plt.legend(loc='best')

plt.figure()
plt.plot(time_tab,mass_disk_tab, marker='.',color='green', label=ur'$R_{disk}$')
plt.xlabel(ur'Temps ($Myr$)')
plt.ylabel(ur'Masse du disque ($Msun$)')
plt.legend(loc='best')

plt.figure()
plt.plot(time_tab,min_rad_loc_tab, marker='.',color='red', label=ur'$R_{disk}$')
plt.xlabel(ur'Temps ($Myr$)')
plt.ylabel(ur'Rayon min du disque ($AU$)')
plt.legend(loc='best')

plt.figure()
plt.plot(time_tab,max_rad_loc_tab, marker='.',color='darkorange', label=ur'$R_{disk}$')
plt.xlabel(ur'Temps ($Myr$)')
plt.ylabel(ur'Rayon max du disque ($AU$)')
plt.legend(loc='best')

plt.figure()
plt.plot(time_tab,log_rho_disk_tab, marker='.',color='black', label=ur'$R_{disk}$')
plt.xlabel(ur'Temps ($Myr$)')
plt.ylabel(ur'log_rho_tab_disk')
plt.legend(loc='best')

plt.show()
'''



#plt.savefig('/gpfs/data1/averliat/analyses/B335_noturb_norot_hydro_pert_asym_aleatoire_tres_forte/Rayon_disque_vs_time_cor_PDF_vs_J.pdf', bbox_inches='tight')

