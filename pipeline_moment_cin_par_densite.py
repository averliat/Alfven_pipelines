import numpy as np

import pymses
import sys
import glob as glob
import os 

from pymses.filters import CellsToPoints
from pymses.utils import constants as cst
import matplotlib.pyplot as plt
plt.ion()

from astropy.io import fits

import os



#------------------------------------------
#Entree pour la simulation et la sauvegarde
#------------------------------------------
#Chemin de la simulation et numeros des outputs
path='/gpfs/data1/phennebe/B335_noturb_norot_hydro_pert_llf/'

num_min = 1
num_max=48

num_tot = num_max-num_min+1


#Pour la sauvegarde
save_output = False  #Mettre False pour ne pas sauvegarder les outputs sous forme de tableaux directement lisibles
save_final = True

dir_save = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_pert_llf/'

dir_save_output = 'save_output_barycentre'
file_save_output = '/save_output_barycentre'

file_save = dir_save + 'Moment_cinetique_tot_output_'+str(num_min)+'_'+str(num_max)+'_par_densite.fits'
file_save_2 = dir_save + 'Integrale_valeurs_absolues_moments_cinetiques_output_'+str(num_min)+'_'+str(num_max)+'_par_densite.fits'


class_rho=np.array([1e-21, 1e-18, 1e-15, 1e-12, 1e-9, 1e-6])  #A classer dans l'odre croissant



#---------------------------------------------
#Creation du dossier de sauvegarde des outputs
#---------------------------------------------
if save_output == True:
	if os.path.isdir(dir_save + dir_save_output) == False:
		os.mkdir(dir_save + dir_save_output)



#-----------------------------------------------------
#Initialisation des tableaux de stockage des resultats
#-----------------------------------------------------
J_tot_tab = np.zeros((num_tot, len(class_rho)+1))
J_abs_tot_tab = np.zeros((num_tot, len(class_rho)+1))



#---------------------------------------------
#Debut de la boucle sur les differents outputs
#---------------------------------------------
for l in range(num_tot):
	num = num_min + l


	#-------------------
	#Lecture des donnees
	#-------------------
	#Verification que l'output a deja ete lue
	if os.path.isfile(dir_save + dir_save_output + file_save_output + str(num) + '.npz') == False:
		#Lecture puis sauvegarde des tableaux si l'output n'a jamais ete lue
		ro=pymses.RamsesOutput(path,num)
		#lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)
		lbox_pc = ro.info['unit_length'].express(cst.pc)

		amr = ro.amr_source(["rho","vel","P","phi","g"])

		cell_source = CellsToPoints(amr)
		cells = cell_source.flatten()
		dx = cells.get_sizes()

		pos = cells.points
		vel = cells["vel"]
		rho = cells["rho"]


		#------------------------------
		#Conversion en unites physiques
		#------------------------------
		#lbox_m = ro.info['unit_length'].express(cst.m)
		factor_dist= ro.info['unit_length'].express(cst.m)
		factor_vel = ro.info['unit_velocity'].express(cst.m/cst.s)
		factor_rho = ro.info['unit_density'].express(cst.kg/(cst.m)**3)

		dx *= factor_dist
		pos *= factor_dist
		vel *= factor_vel
		rho *= factor_rho


		#--------------------------------------------
		#Sauvegarde des tableaux dans un fichier .npz
		#--------------------------------------------
		if save_output == True:
			np.savez(dir_save + dir_save_output + file_save_output + str(num), dx=dx, pos=pos, vel=vel, rho=rho)


	#Lecture des tableaux sauvegardes si l'output a deja ete lue
	else:
		output_tab = np.load(dir_save + dir_save_output + file_save_output + str(num) + '.npz')
		print("Lecture output_"+str(num))

		dx = output_tab['dx']
		pos = output_tab['pos']
		vel = output_tab['vel']
		rho = output_tab['rho']

	#-------------------------------------------------------------------------------------
	#A partir d'ici les tableaux sont charges et sauvegardes si ce n'etait pas deja le cas
	#-------------------------------------------------------------------------------------
	#Calcul du barycentre pour chaque output
	mass_tot = np.sum( rho * dx**3 )
	coord_G_x = 1./mass_tot * np.sum( rho * dx**3 * pos[:,0] )
	coord_G_y = 1./mass_tot * np.sum( rho * dx**3 * pos[:,1] )
	coord_G_z = 1./mass_tot * np.sum( rho * dx**3 * pos[:,2] )
	coord_G = np.array([coord_G_x,coord_G_y,coord_G_z])
	
	pos_moins_centre = pos - coord_G



	#------------------------------------------------------------------
	#Calcul du moment cinetique de chaque cellule par rapport au centre
	#------------------------------------------------------------------
	#radius_cell_tab = np.sqrt( np.sum( (pos[:,:]-centre)**2 , axis=1 ) )  #Tableau contenant le rayon de chaque cellule

	for u in range(len(class_rho)+1):
		if u < len(class_rho):
			print("Calcul pour rho = "+str(class_rho[u]))
		else:
			print("Calcul pour rho > "+str(class_rho[u-1]))
		#current_rho = np.ma.masked_where(rho > class_rho[u], rho)  #On masque toutes les valeurs de rho superieures a class_rho[u] (valeurs masquees quand condition vraie)

		#-------------------------------------------------------
		#Creation du masque pour selectionner les valeurs de rho
		#-------------------------------------------------------
		if u==0:
			mask = rho < class_rho[u]
		elif u != len(class_rho):
			mask = (rho >= class_rho[u-1]) & (rho < class_rho[u])
		else: #u==len(class_rho)
			mask = rho >= class_rho[u-1]


		#--------------------------------
		#Masquage des tableaux et calculs
		#--------------------------------
		rho_masked = rho[mask]
		pos_moins_centre_masked = pos_moins_centre[mask]
		vel_masked = vel[mask]
		dx_masked = dx[mask]

		if len(rho_masked)==0:
			J_tot_tab[l,u] = 0
			J_abs_tot_tab[l,u] = 0
		else:
			J_tot_abs = 0
			moment_tot_x = 0
			moment_tot_y = 0
			moment_tot_z = 0

			debug=0
			for i in range(len(rho_masked)):  #Boucle sur les differentes cellules
				#debug+=1
				if debug==100000:
					print(i)
					debug=0
				#ind = indice_cell[mask][i]
				#Termes du produit vectoriel du rayon par la vitesse
				premier_terme =	pos_moins_centre_masked[i,1]*vel_masked[i,2] - pos_moins_centre_masked[i,2]*vel_masked[i,1]
				deuxieme_terme = pos_moins_centre_masked[i,2]*vel_masked[i,0] - pos_moins_centre_masked[i,0]*vel_masked[i,2]
				troisieme_terme = pos_moins_centre_masked[i,0]*vel_masked[i,1] - pos_moins_centre_masked[i,1]*vel_masked[i,0]

				#Composantes du moment cinetique total
				moment_tot_x += rho_masked[i] * (dx_masked[i]**3) * premier_terme
				moment_tot_y += rho_masked[i] * (dx_masked[i]**3) * deuxieme_terme
				moment_tot_z += rho_masked[i] * (dx_masked[i]**3) * troisieme_terme
		
				#Norme du produit vectoriel
				norm_rad_vect_vel = np.sqrt( premier_terme**2 + deuxieme_terme**2 + troisieme_terme**2 )

				#Moment cinetique absolue de la cellule consideree
				J_tot_abs += rho_masked[i] * (dx_masked[i]**3) * norm_rad_vect_vel


			#-------------------------------------------------
			#Masse total dans la tranche de densite consideree
			#-------------------------------------------------
			#M_tot_rho = np.sum(rho_masked * dx_masked**3)


			#--------------------------------------------------
			#Moment cinetique par rapport au centre et stockage
			#--------------------------------------------------
			J_tot = np.sqrt( moment_tot_x**2 + moment_tot_y**2 + moment_tot_z**2 )
			J_tot_tab[l,u] = J_tot #/ M_tot_rho

			#-----------------------------------------------------------------------
			#Valeur absolue des moments cinetiques par rapport au centre et stockage
			#-----------------------------------------------------------------------
			J_abs_tot_tab[l,u] = J_tot_abs #/ M_tot_rho



#-----------------------------------------------------------------------
#Sauvegarde du tableau contenant les moments cinetiques de chaque output
#-----------------------------------------------------------------------
if save_final == True:
	hdu = fits.PrimaryHDU(J_tot_tab)
	hdu.writeto(file_save)

	hdu = fits.PrimaryHDU(J_abs_tot_tab)
	hdu.writeto(file_save_2)




#print("Avancement : ",round((n+1)/N_cluster*100.,2),"%", end="\r")
