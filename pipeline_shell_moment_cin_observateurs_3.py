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
simu = 'B335_noturb_norot_hydro_hr2_copie_phennebe'#pert_asym_aleatoire50_vhr'
owner = 'averliat'
output_min = 1
output_max = 10




#Sauvegarde des quantites finales:
save_tab = True
#file_save = 'Shell_moment_cin_et_cumul_observateurs_output'+str(output)+'_rmin='+str(R_min)+'AU_rmax='+str(R_max)+'_dr='+str(dr)+'.hdf5'
file_save = 'Moment_cin_fixe_et_err_output'+str(output_min)+'_to_'+str(output_max)+'.hdf5'


#Chemin de la simulation et numeros des outputs
if owner=='averliat':
	path='/mnt/magmist/magmist/simu_B335_averliat/'+simu+'/'
	path_analyses='/mnt/magmist/magmist/simu_B335_averliat/analyses/'+simu+'/'
if owner=='phennebe':
	path='/drf/projets/capucine/phennebe/'+simu+'/'
	path_analyses='/mnt/magmist/magmist/simu_B335_averliat/analyses/'+simu+'/'
if owner=='sapanais':
	path='/dsm/anais/storageA/magmist/'+simu+'/'
	path_analyses='/mnt/magmist/magmist/simu_B335_averliat/analyses/'+simu+'/'



nbre_output = output_max-output_min+1

norm_mom_cin_obs = np.zeros(nbre_output)
simulation_time = np.zeros(nbre_output)
norm_GC_AU = np.zeros(nbre_output)
norm_center_C_AU = np.zeros(nbre_output)
norm_center_G_AU = np.zeros(nbre_output)
vel_GC_tab = np.zeros(nbre_output)
norm_mom_cin_obs_par_integ = np.zeros(nbre_output)
erreur_mom = np.zeros(nbre_output)


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
			print("Pas assez d'output pour la simulation : "+path_simu_a_comparer )
			print("=========================================================================================")
			print("=========================================================================================")
			break
		continue

	verif=1
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
	factor_time_Myr = ro.info['unit_time'].express(cst.Myr)

	simulation_time[nbre_output_effectif] = ro.info['time']*factor_time_Myr


	dx *= factor_dist
	pos *= factor_dist
	vel *= factor_vel
	rho *= factor_rho



	#-------------------------
	#Definition de la coquille
	#-------------------------
	AU = 149597870700 #m

	#Position des differents points
	lbox_m = factor_dist
	center = np.array(([lbox_m/2., lbox_m/2., lbox_m/2.]))


	#arg_centre = np.argmax(rho)
	#pos_C = [pos[:,0][arg_centre],pos[:,1][arg_centre],pos[:,2][arg_centre]]
	#vel_C = [vel[:,0][arg_centre],vel[:,1][arg_centre],vel[:,2][arg_centre]]


	mass_tot = np.sum( rho * dx**3 )
	#coord_G_x = 1./mass_tot * np.sum( rho * dx**3 * pos[:,0] )
	#coord_G_y = 1./mass_tot * np.sum( rho * dx**3 * pos[:,1] )
	#coord_G_z = 1./mass_tot * np.sum( rho * dx**3 * pos[:,2] )
	#pos_G = np.array([coord_G_x,coord_G_y,coord_G_z])


	#norm_GC_AU[nbre_output_effectif] = np.sqrt( np.sum( ((pos_G - pos_C)/AU)**2 ) )
	#norm_center_C_AU[nbre_output_effectif] = np.sqrt( np.sum( ((center - pos_C)/AU)**2 ) )
	#norm_center_G_AU[nbre_output_effectif] = np.sqrt( np.sum( ((center - pos_G)/AU)**2 ) )

	'''
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
	'''


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




	#Calcul du moment cinetique observateurs
	radius_cell_C = np.sqrt( np.sum( (pos[:,:]-center)**2 , axis=1 ) )  #Tableau contenant le rayon de chaque cellule par rapport a C
	pos_cell_from_C = pos[:,:]-center
	vel_cell_from_C = vel[:,:]

	#GC = pos_C - pos_G
	#vel_G_C = vel_C

	#vel_GC_tab[nbre_output_effectif] = np.sqrt( np.sum( np.array(vel_G_C)**2 ))

	#Termes du produit vectoriel du rayon par la vitesse
	#premier_terme =	GC[1]*vel_G_C[2] - GC[2]*vel_G_C[1]
	#deuxieme_terme = GC[2]*vel_G_C[0] - GC[0]*vel_G_C[2]
	#troisieme_terme = GC[0]*vel_G_C[1] - GC[1]*vel_G_C[0]

	#Norme du moment cinetique observateur
	#norm_mom_cin_obs[nbre_output_effectif] = mass_tot * np.sqrt( premier_terme**2 + deuxieme_terme**2 + troisieme_terme**2 )








	#------------------------------------------------------------------
	#Calcul du moment cinetique de chaque cellule par rapport au centre
	#------------------------------------------------------------------	
	moment_x = 0
	moment_y = 0
	moment_z = 0

	erreur_moment = 0

	for i in range(len(dx)):  #Boucle sur les differentes cellules dans la coquille
		#Termes du produit vectoriel du rayon par la vitesse
		premier_terme =	pos_cell_from_C[i,1]*vel_cell_from_C[i,2] - pos_cell_from_C[i,2]*vel_cell_from_C[i,1]
		deuxieme_terme = pos_cell_from_C[i,2]*vel_cell_from_C[i,0] - pos_cell_from_C[i,0]*vel_cell_from_C[i,2]
		troisieme_terme = pos_cell_from_C[i,0]*vel_cell_from_C[i,1] - pos_cell_from_C[i,1]*vel_cell_from_C[i,0]

		#Composantes du moment cinetique total
		moment_x += rho[i] * (dx[i]**3) * premier_terme
		moment_y += rho[i] * (dx[i]**3) * deuxieme_terme
		moment_z += rho[i] * (dx[i]**3) * troisieme_terme

		#Calcul des erreurs
		exp_rho = np.float([c for c in str("%e"%rho[i])][-3]+[c for c in str("%e"%rho[i])][-2]+[c for c in str("%e"%rho[i])][-1])
		delta_rho = 5*10**(-16+exp_rho)

		#exp_pos1 = np.float([c for c in str("%e"%pos_cell_from_C[i,1])][-3]+[c for c in str("%e"%pos_cell_from_C[i,1])][-2]+[c for c in str("%e"%pos_cell_from_C[i,1])][-1])
		#delta_pos = 5*10**(-16+exp_pos1)

		#exp_vel1 = np.float([c for c in str("%e"%vel_cell_from_C[i,1])][-3]+[c for c in str("%e"%vel_cell_from_C[i,1])][-2]+[c for c in str("%e"%vel_cell_from_C[i,1])][-1])
		#delta_vel = 5*10**(-16+exp_vel1

		exp_1 = np.float([c for c in str("%e"%premier_terme)][-3]+[c for c in str("%e"%premier_terme)][-2]+[c for c in str("%e"%premier_terme)][-1])
		delta_1 = 5*10**(-16+exp_1)	

		exp_2 = np.float([c for c in str("%e"%deuxieme_terme)][-3]+[c for c in str("%e"%deuxieme_terme)][-2]+[c for c in str("%e"%deuxieme_terme)][-1])
		delta_2 = 5*10**(-16+exp_2)	

		exp_3 = np.float([c for c in str("%e"%troisieme_terme)][-3]+[c for c in str("%e"%troisieme_terme)][-2]+[c for c in str("%e"%troisieme_terme)][-1])
		delta_3 = 5*10**(-16+exp_3)	
		


		df_drho = dx[i]**3*(premier_terme+deuxieme_terme+troisieme_terme)
		df_d1 = rho[i]*dx[i]**3
		df_d2 = rho[i]*dx[i]**3
		df_d3 = rho[i]*dx[i]**3

		erreur_moment += (df_drho*delta_rho)**2 + (df_d1*delta_1)**2 + (df_d2*delta_2)**2 + (df_d3*delta_3)**2


	erreur_moment = np.sqrt(erreur_moment)

	#--------------------------------------------------
	#Moment cinetique par rapport au centre et stockage
	#--------------------------------------------------
	norm_mom_cin_obs_par_integ[nbre_output_effectif] = np.sqrt( moment_x**2 + moment_y**2 + moment_z**2 )
	erreur_mom[nbre_output_effectif] = erreur_moment

	nbre_output_effectif += 1



#norm_mom_cin_obs = norm_mom_cin_obs[0:nbre_output_effectif]
simulation_time = simulation_time[0:nbre_output_effectif]
#norm_GC_AU = norm_GC_AU[0:nbre_output_effectif]
#norm_center_C_AU = norm_center_C_AU[0:nbre_output_effectif]
#norm_center_G_AU = norm_center_G_AU[0:nbre_output_effectif]
#vel_GC_tab = vel_GC_tab[0:nbre_output_effectif]
norm_mom_cin_obs_par_integ = norm_mom_cin_obs_par_integ[0:nbre_output_effectif]
erreur_mom = erreur_mom[0:nbre_output_effectif]



#-----------------------------------------------------------------------
#Sauvegarde du tableau contenant les moments cinetiques de chaque output
#-----------------------------------------------------------------------
if (save_tab == True and verif==1):
	h5f = h5py.File(path_analyses+file_save, 'w')

	h5f.create_dataset('norm_mom_cin_obs', data=norm_mom_cin_obs)
	h5f.create_dataset('simulation_time',data=simulation_time)
	h5f.create_dataset('norm_GC_AU',data=norm_GC_AU)
	h5f.create_dataset('norm_center_C_AU',data=norm_center_C_AU)
	h5f.create_dataset('norm_center_G_AU',data=norm_center_G_AU)
	h5f.create_dataset('vel_GC_tab',data=vel_GC_tab)
	h5f.create_dataset('norm_mom_cin_obs_par_integ',data=norm_mom_cin_obs_par_integ)
	h5f.create_dataset('erreur_mom',data=erreur_mom)

	h5f.close()



'''
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
'''
