import numpy as np
import h5py



#---------------------
#Parametres a modifier
#---------------------
simu='B335_noturb_norot_hydro_pert_asym_aleatoire_bigbox_50pourc_sink_seuil_haut_MHD_lr'

num_min1 = 32
num_max1 = 291

num_min2 = num_max1 + 1
num_max2 = 367


dir_name = '/home/averliat/these/analyses/'+simu+'/'
tab_name_pref = 'Rayon_disque_par_pdf_modal_bin4_'

path1=dir_name+tab_name_pref+str(num_min1)+'_'+str(num_max1)+'.hdf5'
path2=dir_name+tab_name_pref+str(num_min2)+'_'+str(num_max2)+'.hdf5'

path_fin=dir_name+tab_name_pref+str(num_min1)+'_'+str(num_max2)+'.hdf5'




#----------------------------
#Pour moments cinetique total
#----------------------------
h5f = h5py.File(path1, 'r')
mean_rad_loc_tab1=h5f['mean_rad_loc_tab'][:]
time_tab1=h5f['time_tab'][:]
time_cor_tab1=h5f['time_cor_tab'][:]
min_rad_loc_tab1=h5f['min_rad_loc_tab'][:]
max_rad_loc_tab1=h5f['max_rad_loc_tab'][:]
numero_output1=h5f['numero_output'][:]
mass_disk_tab1=h5f['mass_disk_tab'][:]
log_rho_disk_tab1=h5f['log_rho_disk_tab'][:]
mag_mean_broad_tab1=h5f['mag_mean_broad_tab'][:]
median_rad_loc_tab1=h5f['median_rad_loc_tab'][:]
rad_estim_tab1=h5f['rad_estim_tab'][:]
h5f.close()

h5f = h5py.File(path2, 'r')
mean_rad_loc_tab2=h5f['mean_rad_loc_tab'][:]
time_tab2=h5f['time_tab'][:]
time_cor_tab2=h5f['time_cor_tab'][:]
min_rad_loc_tab2=h5f['min_rad_loc_tab'][:]
max_rad_loc_tab2=h5f['max_rad_loc_tab'][:]
numero_output2=h5f['numero_output'][:]
mass_disk_tab2=h5f['mass_disk_tab'][:]
log_rho_disk_tab2=h5f['log_rho_disk_tab'][:]
mag_mean_broad_tab2=h5f['mag_mean_broad_tab'][:]
median_rad_loc_tab2=h5f['median_rad_loc_tab'][:]
rad_estim_tab2=h5f['rad_estim_tab'][:]
h5f.close()


mean_rad_loc_tab = np.concatenate((mean_rad_loc_tab1,mean_rad_loc_tab2))
time_tab = np.concatenate((time_tab1,time_tab2))
time_cor_tab = np.concatenate((time_cor_tab1,time_cor_tab2))
min_rad_loc_tab = np.concatenate((min_rad_loc_tab1,min_rad_loc_tab2))
max_rad_loc_tab = np.concatenate((max_rad_loc_tab1,max_rad_loc_tab2))
numero_output = np.concatenate((numero_output1,numero_output2))
mass_disk_tab = np.concatenate((mass_disk_tab1,mass_disk_tab2))
log_rho_disk_tab = np.concatenate((log_rho_disk_tab1,log_rho_disk_tab2))
mag_mean_broad_tab = np.concatenate((mag_mean_broad_tab1,mag_mean_broad_tab2))
median_rad_loc_tab = np.concatenate((median_rad_loc_tab1,median_rad_loc_tab2))
rad_estim_tab = np.concatenate((rad_estim_tab1,rad_estim_tab2))



h5f = h5py.File(path_fin, 'w')
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







