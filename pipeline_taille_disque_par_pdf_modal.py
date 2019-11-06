# -*- coding: utf-8 -*-
import extract_disk as ed
import os 
import numpy as np
import matplotlib.pyplot as plt
import h5py
import pymses
plt.ion()
import pipeline_temps_0_simulation as t_0



#-------------------------------------------------------
#Entrer le nom de la simulation et le numero de l'output	
#-------------------------------------------------------
simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire_bigbox_50pourc_sink_seuil_haut_MHD_lr'
owner = 'averliat_alfven'

output_min = 292#'None'
output_max = 367

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
	file_save = 'Rayon_disque_par_pdf_modal_bin4_'+str(output_min)+'_'+str(output_max)+'.hdf5'



#----------------------------------------------------
#Initialisation des differents tableaux de sauvegarde	
#----------------------------------------------------
numero_output = np.zeros(nmbr_output)
time_tab = np.zeros(nmbr_output)
time_cor_tab = np.zeros(nmbr_output)
mass_disk_tab = np.zeros(nmbr_output)
max_rad_loc_tab = np.zeros(nmbr_output)
min_rad_loc_tab = np.zeros(nmbr_output)
mean_rad_loc_tab = np.zeros(nmbr_output)
log_rho_disk_tab = np.zeros(nmbr_output)
mag_mean_broad_tab = np.zeros(nmbr_output)
median_rad_loc_tab = np.zeros(nmbr_output)
rad_estim_tab = np.zeros(nmbr_output)



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

    time_tab[nbre_output_effectif], mass_disk_tab[nbre_output_effectif], max_rad_loc_tab[nbre_output_effectif], min_rad_loc_tab[nbre_output_effectif], mean_rad_loc_tab[nbre_output_effectif], log_rho_disk_tab[nbre_output_effectif], mag_mean_broad_tab[nbre_output_effectif], median_rad_loc_tab[nbre_output_effectif], rad_estim_tab[nbre_output_effectif] = ed.pdf_to_singledisc_rad(path=path, num=num)
    numero_output[nbre_output_effectif] = num
    time_cor_tab[nbre_output_effectif] = time_tab[nbre_output_effectif]-time_0

    nbre_output_effectif += 1



#-----------------------------------------------------------------------
#Sauvegarde du tableau contenant les moments cinetiques de chaque output
#-----------------------------------------------------------------------
mean_rad_loc_tab = mean_rad_loc_tab[0:nbre_output_effectif]
time_tab = time_tab[0:nbre_output_effectif]
time_cor_tab = time_cor_tab[0:nbre_output_effectif]
min_rad_loc_tab = min_rad_loc_tab[0:nbre_output_effectif]
max_rad_loc_tab = max_rad_loc_tab[0:nbre_output_effectif]
numero_output = numero_output[0:nbre_output_effectif]
mass_disk_tab = mass_disk_tab[0:nbre_output_effectif]
log_rho_disk_tab = log_rho_disk_tab[0:nbre_output_effectif]
mag_mean_broad_tab = mag_mean_broad_tab[0:nbre_output_effectif]
median_rad_loc_tab = median_rad_loc_tab[0:nbre_output_effectif]
rad_estim_tab = rad_estim_tab[0:nbre_output_effectif]




if (save_tab == True and verif==1):
    h5f = h5py.File(path_analyses+file_save,'w') 

    h5f.create_dataset('mean_rad_loc_tab', data=mean_rad_loc_tab)
    h5f.create_dataset('time_tab', data=time_tab)
    h5f.create_dataset('time_cor_tab', data=time_cor_tab)
    h5f.create_dataset('min_rad_loc_tab', data=min_rad_loc_tab)
    h5f.create_dataset('max_rad_loc_tab', data=max_rad_loc_tab)
    h5f.create_dataset('numero_output', data=numero_output)
    h5f.create_dataset('mass_disk_tab', data=mass_disk_tab)
    h5f.create_dataset('log_rho_disk_tab', data=log_rho_disk_tab)
    h5f.create_dataset('mag_mean_broad_tab', data=mag_mean_broad_tab)
    h5f.create_dataset('median_rad_loc_tab', data=median_rad_loc_tab)
    h5f.create_dataset('rad_estim_tab', data=rad_estim_tab)


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

