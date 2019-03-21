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
path='/gpfs/data1/averliat/B335_noturb_norot_hydro_pert_asym_manuelle/'

num_min = 37
num_max = 60

num_tot = num_max-num_min+1


#Pour la sauvegarde
save_output = False  #Mettre False pour ne pas sauvegarder les outputs sous forme de tableaux directement lisibles
save_final = True

dir_save = '/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_pert_asym_manuelle/'

dir_save_output = 'save_output_barycentre'
file_save_output = '/save_output_barycentre'

file_save = dir_save + 'Moment_cinetique_tot_output_'+str(num_min)+'_'+str(num_max)+'_barycentre.fits'
file_save_2 = dir_save + 'Integrale_valeurs_absolues_moments_cinetiques_output_'+str(num_min)+'_'+str(num_max)+'_barycentre.fits'



#---------------------------------------------
#Creation du dossier de sauvegarde des outputs
#---------------------------------------------
if save_output == True:
	if os.path.isdir(dir_save + dir_save_output) == False:
		os.mkdir(dir_save + dir_save_output)



#-----------------------------------------------------
#Initialisation des tableaux de stockage des resultats
#-----------------------------------------------------
J_tot_tab = np.zeros(num_tot)
J_abs_tot_tab = np.zeros(num_tot)



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


	#A partir d'ici les tableaux sont charges et sauvegardes si ce n'etait pas deja le cas

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

	J_tot_abs = 0
	moment_tot_x = 0
	moment_tot_y = 0
	moment_tot_z = 0

	for i in range(len(dx)):  #Boucle sur les differentes cellules
		#Termes du produit vectoriel du rayon par la vitesse
		premier_terme =	pos_moins_centre[i,1]*vel[i,2] - pos_moins_centre[i,2]*vel[i,1]
		deuxieme_terme = pos_moins_centre[i,2]*vel[i,0] - pos_moins_centre[i,0]*vel[i,2]
		troisieme_terme = pos_moins_centre[i,0]*vel[i,1] - pos_moins_centre[i,1]*vel[i,0]

		#Composantes du moment cinetique total
		moment_tot_x += rho[i] * (dx[i]**3) * premier_terme
		moment_tot_y += rho[i] * (dx[i]**3) * deuxieme_terme
		moment_tot_z += rho[i] * (dx[i]**3) * troisieme_terme
		
		#Norme du produit vectoriel
		norm_rad_vect_vel = np.sqrt( premier_terme**2 + deuxieme_terme**2 + troisieme_terme**2 )

		#Moment cinetique absolue de la cellule consideree
		J_tot_abs += rho[i] * (dx[i]**3) * norm_rad_vect_vel



	#--------------------------------------------------
	#Moment cinetique par rapport au centre et stockage
	#--------------------------------------------------
	J_tot = np.sqrt( moment_tot_x**2 + moment_tot_y**2 + moment_tot_z**2 )
	J_tot_tab[l] = J_tot

	#-----------------------------------------------------------------------
	#Valeur absolue des moments cinetiques par rapport au centre et stockage
	#-----------------------------------------------------------------------
	J_abs_tot_tab[l] = J_tot_abs



#-----------------------------------------------------------------------
#Sauvegarde du tableau contenant les moments cinetiques de chaque output
#-----------------------------------------------------------------------
if save_final == True:
	hdu = fits.PrimaryHDU(J_tot_tab)
	hdu.writeto(file_save)

	hdu = fits.PrimaryHDU(J_abs_tot_tab)
	hdu.writeto(file_save_2)




#print("Avancement : ",round((n+1)/N_cluster*100.,2),"%", end="\r")
