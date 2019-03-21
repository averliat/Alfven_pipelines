import numpy as np
import h5py



#---------------------
#Parametres a modifier
#---------------------
simu='B335_noturb_norot_hydro_pert_asym_aleatoire_shr_bigbox_50pourc'

num_min1 = 1
num_max1 = 59

num_min2 = num_max1 + 1
num_max2 = 110


dir_name = '/home/averliat/these/analyses/'+simu+'/'
tab_name_pref = 'Moment_cin_observateurs_output'

path1=dir_name+tab_name_pref+str(num_min1)+'_to_'+str(num_max1)+'.hdf5'
path2=dir_name+tab_name_pref+str(num_min2)+'_to_'+str(num_max2)+'.hdf5'

path_fin=dir_name+tab_name_pref+str(num_min1)+'_to_'+str(num_max2)+'.hdf5'




#----------------------------
#Pour moments cinetique total
#----------------------------
h5f = h5py.File(path1, 'r')
norm_mom_cin_obs1=h5f['norm_mom_cin_obs'][:]
simulation_time1=h5f['simulation_time'][:]
norm_GC_AU1=h5f['norm_GC_AU'][:]
norm_center_C_AU1=h5f['norm_center_C_AU'][:]
norm_center_G_AU1=h5f['norm_center_G_AU'][:]
vel_GC_tab1=h5f['vel_GC_tab'][:]
norm_mom_cin_obs_par_integ1=h5f['norm_mom_cin_obs_par_integ'][:]
numero_output1=h5f['numero_output'][:]
h5f.close()

h5f = h5py.File(path2, 'r')
norm_mom_cin_obs2=h5f['norm_mom_cin_obs'][:]
simulation_time2=h5f['simulation_time'][:]
norm_GC_AU2=h5f['norm_GC_AU'][:]
norm_center_C_AU2=h5f['norm_center_C_AU'][:]
norm_center_G_AU2=h5f['norm_center_G_AU'][:]
vel_GC_tab2=h5f['vel_GC_tab'][:]
norm_mom_cin_obs_par_integ2=h5f['norm_mom_cin_obs_par_integ'][:]
numero_output2=h5f['numero_output'][:]
h5f.close()


norm_mom_cin_obs = np.concatenate((norm_mom_cin_obs1,norm_mom_cin_obs2))
simulation_time = np.concatenate((simulation_time1,simulation_time2))
norm_GC_AU = np.concatenate((norm_GC_AU1,norm_GC_AU2))
norm_center_C_AU = np.concatenate((norm_center_C_AU1,norm_center_C_AU2))
norm_center_G_AU = np.concatenate((norm_center_G_AU1,norm_center_G_AU2))
vel_GC_tab = np.concatenate((vel_GC_tab1,vel_GC_tab2))
norm_mom_cin_obs_par_integ = np.concatenate((norm_mom_cin_obs_par_integ1,norm_mom_cin_obs_par_integ2))
numero_output = np.concatenate((numero_output1,numero_output2))



h5f = h5py.File(path_fin, 'w')
h5f.create_dataset('norm_mom_cin_obs', data=norm_mom_cin_obs)
h5f.create_dataset('simulation_time',data=simulation_time)
h5f.create_dataset('norm_GC_AU',data=norm_GC_AU)
h5f.create_dataset('norm_center_C_AU',data=norm_center_C_AU)
h5f.create_dataset('norm_center_G_AU',data=norm_center_G_AU)
h5f.create_dataset('vel_GC_tab',data=vel_GC_tab)
h5f.create_dataset('norm_mom_cin_obs_par_integ',data=norm_mom_cin_obs_par_integ)
h5f.create_dataset('numero_output',data=numero_output)
h5f.close()







