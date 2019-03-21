import numpy as np
import h5py



#---------------------
#Parametres a modifier
#---------------------
simu='B335_noturb_norot_hydro_pert_asym_aleatoire50_emhr'

num_min1 = 6
num_max1 = 225

num_min2 = num_max1 + 1
num_max2 = 690


dir_name = '/home/averliat/these/analyses/'+simu+'/'
tab_name_pref = 'Moment_cinetique_tot_et_abs_output_'

path1=dir_name+tab_name_pref+str(num_min1)+'_'+str(num_max1)+'.hdf5'
path2=dir_name+tab_name_pref+str(num_min2)+'_'+str(num_max2)+'.hdf5'

path_fin=dir_name+tab_name_pref+str(num_min1)+'_'+str(num_max2)+'.hdf5'




#----------------------------
#Pour moments cinetique total
#----------------------------
h5f = h5py.File(path1, 'r')
J_tot_tab1=h5f['J_tot_tab'][:]
J_abs_tot_tab1=h5f['J_abs_tot_tab'][:]
numero_output1=h5f['numero_output'][:]
time_tab1=h5f['time_tab'][:]
h5f.close()

h5f = h5py.File(path2, 'r')
J_tot_tab2=h5f['J_tot_tab'][:]
J_abs_tot_tab2=h5f['J_abs_tot_tab'][:]
numero_output2=h5f['numero_output'][:]
time_tab2=h5f['time_tab'][:]
h5f.close()


J_tot_tab=np.concatenate((J_tot_tab1,J_tot_tab2))
J_abs_tot_tab=np.concatenate((J_abs_tot_tab1,J_abs_tot_tab2))
numero_output=np.concatenate((numero_output1,numero_output2))
time_tab=np.concatenate((time_tab1,time_tab2))


h5f = h5py.File(path_fin, 'w')
h5f.create_dataset('J_tot_tab', data=J_tot_tab)
h5f.create_dataset('J_abs_tot_tab', data=J_abs_tot_tab)
h5f.create_dataset('numero_output', data=numero_output)
h5f.create_dataset('time_tab', data=time_tab)
h5f.close()







