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

from astropy.io import fits
import astropy.coordinates as coord

import os




#-------------------------------------------------------
#Entrer le nom de la simulation et le numero de l'output	
#-------------------------------------------------------
simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire_forte'
owner = 'averliat'
output_max = 66

comparaison_moments_cumules = False

#recalage_temporel = False

#Caracteristiques des coquilles en AU:
R_min = 0
R_max = 20000
dr = 1

#Critere pour delimiter le rayon du disque = difference relative entre le moment total et absolu > critere_r_disk
critere_r_disk = 0.2


#Sauvegarde des quantites finales:
save_tab = True
file_save = 'rayon_disque_'+str(output_max)+'_dr'+str(dr)+'_cumul='+str(comparaison_moments_cumules)+'_critere'+str(critere_r_disk)+'.fits'
file_save2 = 'simulation_time_'+str(output_max)+'_dr'+str(dr)+'_cumul='+str(comparaison_moments_cumules)+'_critere'+str(critere_r_disk)+'.fits'
file_save3 = 'simulation_time_cor_'+str(output_max)+'_dr'+str(dr)+'_cumul='+str(comparaison_moments_cumules)+'_critere'+str(critere_r_disk)+'.fits'


#Chemin de la simulation et numeros des outputs
path='/gpfs/data1/'+owner+'/'+simu+'/'
path_analyses='/gpfs/data1/averliat/analyses/'+simu+'/'


#Output de depart
seuil_rho = 1e-10
if os.path.isfile(path_analyses+'t0_seuil_rho_'+str(seuil_rho)+'.txt') == False:
	ref = np.array(t_0.temps_0_simu(path_simu_ref, seuil_rho))
	np.savetxt(path_simu_ref_analyses+'t0_seuil_rho_'+str(seuil_rho)+'.txt', ref)
else:
	ref=np.loadtxt(path_analyses+'t0_seuil_rho_'+str(seuil_rho)+'.txt')

#if recalage_temporel == True:
output_0 = int(ref[0])
time_0 = ref[1]



#---------------------------------------
#Creation des differents tableaux utiles	
#---------------------------------------
nbre_output=output_max-output_0+1
R_disk = np.zeros(nbre_output)
simulation_time = np.zeros(nbre_output)
simulation_time_cor = np.zeros(nbre_output)



#---------------------------------------------
#Debut de la boucle sur les differents outputs
#---------------------------------------------
for l in range(nbre_output):
	num_output = output_0 + l
	#-------------------
	#Lecture de l'output
	#-------------------
	ro=pymses.RamsesOutput(path,num_output)
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
	factor_time_Myr = ro.info['unit_time'].express(cst.Myr)
	
	simulation_time[l] = ro.info['time']*factor_time_Myr
	simulation_time_cor[l] = ro.info['time']*factor_time_Myr - time_0

	dx *= factor_dist
	pos *= factor_dist
	vel *= factor_vel
	rho *= factor_rho



	#-------------------------
	#Definition de la coquille
	#-------------------------
	#Position du "centre" de la simulation = endroit le plus dense
	arg_centre = np.argmax(rho)
	center = [pos[:,0][arg_centre],pos[:,1][arg_centre],pos[:,2][arg_centre]]

	radius_cell = np.sqrt( np.sum( (pos[:,:]-center)**2 , axis=1 ) )  #Tableau contenant le rayon de chaque cellule
	

	AU = 149597870700 #m

	radius_min = R_min*AU
	radius_max = R_max*AU
	dim_shell = dr*AU
	nbre_shell = int((radius_max-radius_min)/dim_shell)

	radius_shell = np.zeros(nbre_shell)



	#-----------------------------------
	#Initialisation des moments integres
	#-----------------------------------
	J_int = 0.
	J_abs_int = 0.
	


	count=0

	for s in range(nbre_shell):
		if count==100:
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



		#------------------
		#------------------
		#Calcul des moments
		#------------------
		#------------------
		#Vecteurs cellule dans la coquille <-- centre (tableau)
		pos_moins_centre = (shell_pos - center)  #Mettre un moins pour vecteurs cellule --> centre (change les signes de vect_unit, vel_rad, vel_rad_mean, vel_rad_mean_corrigee mais ne change pas V_RMS)


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


		#--------------------------------------------------
		#Moment cinetique par rapport au centre et stockage
		#--------------------------------------------------
		J_shell = np.sqrt( moment_shell_x**2 + moment_shell_y**2 + moment_shell_z**2 )



		#--------------------------------------------------
		#Moment cinetique par rapport au centre et stockage
		#--------------------------------------------------
		J_int += J_shell
		J_abs_int += J_shell_abs



		#--------------------------------------------------
		#Comparaison des moments cinétiques total et absolu
		#--------------------------------------------------
		if comparaison_moments_cumules == False:
			if s>=1:
				diff_prec = diff
			diff = abs((J_shell-J_shell_abs)/J_shell_abs)
		
		
			if s>=1:
				if (diff > diff_prec) and (diff > critere_r_disk):
					R_disk[l] = radius_shell[s]
					break  #Sortie de la boucle sur les coquilles --> passage a l'output suivant

		if comparaison_moments_cumules == True:
			if s>=1:
				diff_prec = diff
			diff = abs((J_int-J_abs_int)/J_abs_int)
		
		
			if s>=1:
				if (diff > diff_prec) and (diff > critere_r_disk):
					R_disk[l] = radius_shell[s]
					break  #Sortie de la boucle sur les coquilles --> passage a l'output suivant



#--------------------------------
#Trace des differentes quantites
#--------------------------------
plt.plot(simulation_time,R_disk, marker='.',color='midnightblue', label=ur'$R_{disk}$')
plt.xlabel(ur'Temps ($Myr$)')
plt.ylabel(ur'Rayon du disque ($AU$)')
plt.legend(loc='best')

plt.figure()
plt.plot(simulation_time_cor,R_disk, marker='.',color='midnightblue', label=ur'$R_{disk}$')
plt.xlabel(ur'Temps corrigé ($Myr$)')
plt.ylabel(ur'Rayon du disque ($AU$)')
plt.legend(loc='best')



#-----------------------------------------------------------------------
#Sauvegarde du tableau contenant les moments cinetiques de chaque output
#-----------------------------------------------------------------------
if save_tab == True:
	hdu = fits.PrimaryHDU(R_disk)
	hdu.writeto(path_analyses+file_save)

	hdu = fits.PrimaryHDU(simulation_time)
	hdu.writeto(path_analyses+file_save2)

	hdu = fits.PrimaryHDU(simulation_time_cor)
	hdu.writeto(path_analyses+file_save3)



#plt.savefig('/gpfs/data1/averliat/analyses/B335_noturb_norot_hydro_pert_asym_aleatoire_tres_forte/Rayon_disque_vs_time_cor.pdf', bbox_inches='tight')

#plt.savefig('/gpfs/data1/averliat/analyses/B335_noturb_norot_hydro_pert_asym_aleatoire_forte/Vitesse_rotation_output12.pdf', bbox_inches='tight')


#print("Avancement : ",round((n+1)/N_cluster*100.,2),"%", end="\r")
