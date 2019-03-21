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
import astropy.coordinates as coord

import os




#-------------------------------------------------------
#Entrer le nom de la simulation et le numero de l'output	
#-------------------------------------------------------
simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire_forte'
owner = 'averliat'
output = 43


#Caracteristiques des coquilles en AU:
R_min = 0
R_max = 'all'  #='all' pour prendre toute la boite
dr = 10

#Calcul par rapport ou barycentre (True) ou au point le plus dense (False)
bar=False


#Sauvegarde des quantites finales:
save_tab = True
file_save = 'Shell_moment_cin_et_cumul_output'+str(output)+'_rmin='+str(R_min)+'AU_rmax='+str(R_max)+'_dr='+str(dr)+'_bar='+str(bar)+'.hdf5'


#Chemin de la simulation et numeros des outputs
path='/gpfs/data1/'+owner+'/'+simu+'/'
path_analyses='/gpfs/data1/averliat/analyses/'+simu+'/'



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
factor_dist= ro.info['unit_length'].express(cst.m)
factor_vel = ro.info['unit_velocity'].express(cst.m/cst.s)
factor_rho = ro.info['unit_density'].express(cst.kg/(cst.m)**3)
factor_time_Myr = ro.info['unit_time'].express(cst.Myr)

simulation_time = ro.info['time']*factor_time_Myr


dx *= factor_dist
pos *= factor_dist
vel *= factor_vel
rho *= factor_rho



#-------------------------
#Definition de la coquille
#-------------------------
AU = 149597870700 #m

#Position du point de calcul des moments
if bar==False:
	arg_centre = np.argmax(rho)
	center = [pos[:,0][arg_centre],pos[:,1][arg_centre],pos[:,2][arg_centre]]
if bar==True:
	mass_tot = np.sum( rho * dx**3 )
	coord_G_x = 1./mass_tot * np.sum( rho * dx**3 * pos[:,0] )
	coord_G_y = 1./mass_tot * np.sum( rho * dx**3 * pos[:,1] )
	coord_G_z = 1./mass_tot * np.sum( rho * dx**3 * pos[:,2] )
	center = np.array([coord_G_x,coord_G_y,coord_G_z])

	arg_centre = np.argmax(rho)
	point_plus_dense = [pos[:,0][arg_centre],pos[:,1][arg_centre],pos[:,2][arg_centre]]

	diff_bar_ptdense_AU = np.sqrt( np.sum( ((center - point_plus_dense)/AU)**2 ) )

	print
	print('==========================================================================')
	print("Difference barycentre - point le plus dense = "+str(diff_bar_ptdense_AU)+' AU')
	print('==========================================================================')
	print
	


radius_cell = np.sqrt( np.sum( (pos[:,:]-center)**2 , axis=1 ) )  #Tableau contenant le rayon de chaque cellule


radius_min = R_min*AU

if R_max == 'all':
	radius_max = int(lbox_au/2.)*AU
else:
	radius_max = R_max*AU

dim_shell = dr*AU
nbre_shell = int((radius_max-radius_min)/dim_shell)

radius_shell = np.zeros(nbre_shell)



#-----------------------------------
#Initialisation des moments integres
#-----------------------------------
J_int = 0.
J_abs_int = 0.
J_x_int = 0.
J_y_int = 0.
J_z_int = 0.
mass_shell_int = 0.


J_shell_tab = np.zeros(nbre_shell)
J_abs_shell_tab = np.zeros(nbre_shell)
J_int_tab = np.zeros(nbre_shell)
J_abs_int_tab = np.zeros(nbre_shell)
J_x_shell = np.zeros(nbre_shell)
J_y_shell = np.zeros(nbre_shell)
J_z_shell = np.zeros(nbre_shell)
J_x_int_tab = np.zeros(nbre_shell)
J_y_int_tab = np.zeros(nbre_shell)
J_z_int_tab = np.zeros(nbre_shell)
J_x_y_z_int_tab = np.zeros(nbre_shell)



#------------------------------------
#Debut de la boucle sur les coquilles
#------------------------------------
count=0

for s in range(nbre_shell):
	if count==10:
		print s
		count=0
	count += 1

	radius_in = radius_min +s*dim_shell
	radius_out = radius_in + dim_shell

	mask = (radius_cell > radius_in) & (radius_cell <= radius_out)

	#Quantites dans la coquille
	shell_pos = pos[mask]
	shell_vel = vel[mask]
	shell_rho = rho[mask]
	shell_dx = dx[mask]
	shell_radius_cell = radius_cell[mask]

	#Masse contenue dans la coquille
	shell_mass = np.sum( shell_rho* shell_dx**3 )

	#Rayon moyen de la coquille
	R_shell = (radius_in + radius_out)/2.
	radius_shell[s] = (radius_in + radius_out)/2. /AU  #Sauvegarde de R_shell en AU pour l'abscisse du trace

	#Vecteurs cellule dans la coquille <-- centre (tableau)
	pos_moins_centre = (shell_pos - center)


	#------------------------------------------------------------------
	#Calcul du moment cinetique de chaque cellule par rapport au centre
	#------------------------------------------------------------------	
	J_shell_abs = 0
	J_shell = 0
	moment_shell_x = 0
	moment_shell_y = 0
	moment_shell_z = 0

	for i in range(len(shell_dx)):  #Boucle sur les differentes cellules dans la coquille
		#Termes du produit vectoriel du rayon par la vitesse
		premier_terme =	pos_moins_centre[i,1]*shell_vel[i,2] - pos_moins_centre[i,2]*shell_vel[i,1]
		deuxieme_terme = pos_moins_centre[i,2]*shell_vel[i,0] - pos_moins_centre[i,0]*shell_vel[i,2]
		troisieme_terme = pos_moins_centre[i,0]*shell_vel[i,1] - pos_moins_centre[i,1]*shell_vel[i,0]

		#Composantes du moment cinetique total
		moment_shell_x += shell_rho[i] * (shell_dx[i]**3) * premier_terme
		moment_shell_y += shell_rho[i] * (shell_dx[i]**3) * deuxieme_terme
		moment_shell_z += shell_rho[i] * (shell_dx[i]**3) * troisieme_terme
	
		#Norme du produit vectoriel
		norm_rad_vect_vel = np.sqrt( premier_terme**2 + deuxieme_terme**2 + troisieme_terme**2 )

		#Moment cinetique absolue de la cellule consideree
		J_shell_abs += shell_rho[i] * (shell_dx[i]**3) * norm_rad_vect_vel

		mass_shell_int += shell_rho[i] * (shell_dx[i]**3)


	#--------------------------------------------------
	#Moment cinetique par rapport au centre et stockage
	#--------------------------------------------------
	J_shell = np.sqrt( moment_shell_x**2 + moment_shell_y**2 + moment_shell_z**2 )


	#------------------------------------------------------------
	#Moments cinetiques cumules par rapport au centre et stockage
	#------------------------------------------------------------
	J_int += J_shell
	J_abs_int += J_shell_abs
	
	J_x_int += moment_shell_x
	J_y_int += moment_shell_y
	J_z_int += moment_shell_z

	J_x_y_z_int = np.sqrt( J_x_int**2 + J_y_int**2 + J_z_int**2 )


	#---------------------------------
	#Sauvegarde des differents moments
	#---------------------------------
	J_shell_tab[s] = J_shell
	J_abs_shell_tab[s] = J_shell_abs
	J_int_tab[s] = J_int
	J_abs_int_tab[s] = J_abs_int
	J_x_shell[s] = moment_shell_x
	J_y_shell[s] = moment_shell_y
	J_z_shell[s] = moment_shell_z
	J_x_int_tab[s] = J_x_int
	J_y_int_tab[s] = J_y_int
	J_z_int_tab[s] = J_z_int
	J_x_y_z_int_tab[s] = J_x_y_z_int




#-----------------------------------------------------------------------
#Sauvegarde du tableau contenant les moments cinetiques de chaque output
#-----------------------------------------------------------------------
if save_tab == True:
	h5f = h5py.File(path_analyses+file_save, 'w')

	h5f.create_dataset('radius_shell', data=radius_shell)
	h5f.create_dataset('J_shell_tab',data=J_shell_tab)
	h5f.create_dataset('J_abs_shell_tab',data=J_abs_shell_tab)
	h5f.create_dataset('J_int_tab',data=J_int_tab)
	h5f.create_dataset('J_abs_int_tab',data=J_abs_int_tab)
	h5f.create_dataset('J_x_shell',data=J_x_shell)
	h5f.create_dataset('J_y_shell',data=J_y_shell)
	h5f.create_dataset('J_z_shell',data=J_z_shell)
	h5f.create_dataset('J_x_int_tab',data=J_x_int_tab)
	h5f.create_dataset('J_y_int_tab',data=J_y_int_tab)
	h5f.create_dataset('J_z_int_tab',data=J_z_int_tab)
	h5f.create_dataset('J_x_y_z_int_tab',data=J_x_y_z_int_tab)

	h5f.close()



#plt.savefig('/gpfs/data1/averliat/analyses/B335_noturb_norot_hydro_pert_asym_aleatoire_tres_forte/Rayon_disque_vs_time_cor.pdf', bbox_inches='tight')

#plt.savefig('/gpfs/data1/averliat/analyses/B335_noturb_norot_hydro_pert_asym_aleatoire_forte/Vitesse_rotation_output12.pdf', bbox_inches='tight')


#print("Avancement : ",round((n+1)/N_cluster*100.,2),"%", end="\r")
