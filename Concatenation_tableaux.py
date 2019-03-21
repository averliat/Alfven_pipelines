import numpy as np
from astropy.io import fits



#---------------------
#Parametres a modifier
#---------------------
dir_name = "/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_pert_asym_manuelle"
tab_name_pref = "/Moment_cinetique_tot_output_"
tab_name_suff = "_barycentre.fits"

num_min11 = 1
num_min12 = 36

num_min21 = num_min12 + 1
num_min22 = 60



#----------------------------
#Pour moments cinetique total
#----------------------------
tab1 = fits.open(dir_name+tab_name_pref+str(num_min11)+'_'+str(num_min12)+tab_name_suff)[0].data
tab2 = fits.open(dir_name+tab_name_pref+str(num_min21)+'_'+str(num_min22)+tab_name_suff)[0].data

tab3= np.concatenate((tab1,tab2))

hdu = fits.PrimaryHDU(tab3)
hdu.writeto(dir_name+tab_name_pref+str(num_min11)+'_'+str(num_min22)+tab_name_suff)



#------------------------------------------
#Pour valeurs absolues du moments cinetique
#------------------------------------------
tab4 = fits.open(dir_name+'/Integrale_valeurs_absolues_moments_cinetiques_output_'+str(num_min11)+'_'+str(num_min12)+'_barycentre.fits')[0].data
tab5 = fits.open(dir_name+'/Integrale_valeurs_absolues_moments_cinetiques_output_'+str(num_min21)+'_'+str(num_min22)+'_barycentre.fits')[0].data

tab6= np.concatenate((tab4,tab5))

hdu = fits.PrimaryHDU(tab6)
hdu.writeto(dir_name+'/Integrale_valeurs_absolues_moments_cinetiques_output_'+str(num_min11)+'_'+str(num_min22)+'_barycentre.fits')


