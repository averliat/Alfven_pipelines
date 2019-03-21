# -*- coding: utf-8 -*-
import numpy as np
import h5py
import sys
import os



#---------------------
#Parametres a modifier
#---------------------
simu='B335_noturb_norot_hydro_pert_asym_aleatoire50_shr'

num_min1 = 10
num_max1 = 685

output_depart=13
modulo=2


dir_name = '/home/averliat/these/analyses/'+simu+'/'
tab_name_pref = 'Moment_cin_observateurs_dans_disque_par_pdf_output'

path1=dir_name+tab_name_pref+str(num_min1)+'_to_'+str(num_max1)+'_complet.hdf5'

path_fin=dir_name+tab_name_pref+str(num_min1)+'_to_'+str(num_max1)+'.hdf5'



#------------------------------------------------------------------
#Verification que le fichier original a bien ete copie au prealable
#------------------------------------------------------------------
if os.path.isfile(path1)==False:
    print("=================================================================================")
    print("Attention, le fichier n'a pas été sauvegardé au nom '*_complet.hdf5' au préalable")
    print("=================================================================================")
    sys.exit()



#----------------------------
#Pour moments cinetique total
#----------------------------
h5f = h5py.File(path1, 'r')
#norm_mom_cin_obs1=h5f['norm_mom_cin_obs'][:]
simulation_time1=h5f['simulation_time'][:]
norm_GC_AU1=h5f['norm_GC_AU'][:]
norm_center_C_AU1=h5f['norm_center_C_AU'][:]
norm_center_G_AU1=h5f['norm_center_G_AU'][:]
vel_GC_tab1=h5f['vel_GC_tab'][:]
norm_mom_cin_obs_par_integ1=h5f['norm_mom_cin_obs_par_integ'][:]
numero_output1=h5f['numero_output'][:]
h5f.close()


#norm_mom_cin_obs=np.zeros(num_max1)
simulation_time=np.zeros(num_max1)
norm_GC_AU=np.zeros(num_max1)
norm_center_C_AU=np.zeros(num_max1)
norm_center_G_AU=np.zeros(num_max1)
vel_GC_tab=np.zeros(num_max1)
norm_mom_cin_obs_par_integ=np.zeros(num_max1)
numero_output=np.zeros(num_max1)




#---------------
#Tri des outputs
#---------------
num_tot=len(numero_output1)

output_effectif=0
count=1
for i in range(num_tot):
    i=i+num_min1
    if (i<=output_depart):
        #norm_mom_cin_obs[i-1]=norm_mom_cin_obs1[i-1]
        simulation_time[output_effectif]=simulation_time1[i-num_min1]
        norm_GC_AU[output_effectif]=norm_GC_AU1[i-num_min1]
        norm_center_C_AU[output_effectif]=norm_center_C_AU1[i-num_min1]
        norm_center_G_AU[output_effectif]=norm_center_G_AU1[i-num_min1]
        vel_GC_tab[output_effectif]=vel_GC_tab1[i-num_min1]
        norm_mom_cin_obs_par_integ[output_effectif]=norm_mom_cin_obs_par_integ1[i-num_min1]
        numero_output[output_effectif]=numero_output1[i-num_min1]
        output_effectif+=1

    else:
        if count==modulo:
            #norm_mom_cin_obs[output_effectif]=norm_mom_cin_obs1[i-1]
            simulation_time[output_effectif]=simulation_time1[i-num_min1]
            norm_GC_AU[output_effectif]=norm_GC_AU1[i-num_min1]
            norm_center_C_AU[output_effectif]=norm_center_C_AU1[i-num_min1]
            norm_center_G_AU[output_effectif]=norm_center_G_AU1[i-num_min1]
            vel_GC_tab[output_effectif]=vel_GC_tab1[i-num_min1]
            norm_mom_cin_obs_par_integ[output_effectif]=norm_mom_cin_obs_par_integ1[i-num_min1]
            numero_output[output_effectif]=numero_output1[i-num_min1]
            output_effectif+=1
            count=1
        else:
            count+=1
        

if len(np.where(numero_output==num_max1)[0])==0:
            #norm_mom_cin_obs[output_effectif]=norm_mom_cin_obs1[num_max1-1]
            simulation_time[output_effectif]=simulation_time1[num_tot-1]
            norm_GC_AU[output_effectif]=norm_GC_AU1[num_tot-1]
            norm_center_C_AU[output_effectif]=norm_center_C_AU1[num_tot-1]
            norm_center_G_AU[output_effectif]=norm_center_G_AU1[num_tot-1]
            vel_GC_tab[output_effectif]=vel_GC_tab1[num_tot-1]
            norm_mom_cin_obs_par_integ[output_effectif]=norm_mom_cin_obs_par_integ1[num_tot-1]
            numero_output[output_effectif]=numero_output1[num_tot-1]
            output_effectif+=1


#norm_mom_cin_obs = norm_mom_cin_obs[0:output_effectif]
simulation_time = simulation_time[0:output_effectif]
norm_GC_AU = norm_GC_AU[0:output_effectif]
norm_center_C_AU = norm_center_C_AU[0:output_effectif]
norm_center_G_AU = norm_center_G_AU[0:output_effectif]
vel_GC_tab = vel_GC_tab[0:output_effectif]
norm_mom_cin_obs_par_integ = norm_mom_cin_obs_par_integ[0:output_effectif]
numero_output = numero_output[0:output_effectif]










h5f = h5py.File(path_fin, 'w')
#h5f.create_dataset('norm_mom_cin_obs', data=norm_mom_cin_obs)
h5f.create_dataset('simulation_time',data=simulation_time)
h5f.create_dataset('norm_GC_AU',data=norm_GC_AU)
h5f.create_dataset('norm_center_C_AU',data=norm_center_C_AU)
h5f.create_dataset('norm_center_G_AU',data=norm_center_G_AU)
h5f.create_dataset('vel_GC_tab',data=vel_GC_tab)
h5f.create_dataset('norm_mom_cin_obs_par_integ',data=norm_mom_cin_obs_par_integ)
h5f.create_dataset('numero_output',data=numero_output)
h5f.close()







