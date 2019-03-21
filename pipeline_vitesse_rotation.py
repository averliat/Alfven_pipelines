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

#from astropy.io import fits
import astropy.coordinates as coord

import os




#-------------------------------------------------------
#Entrer le nom de la simulation et le numero de l'output	
#-------------------------------------------------------
simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire_forte'
owner = 'averliat'
num_output = 20

direction_J = False

#Caracteristiques des coquilles en AU:
R_min = 0
R_max = 800
dr = 1


#Chemin de la simulation et numeros des outputs
path='/gpfs/data1/'+owner+'/'+simu+'/'



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

V_rot_stock = np.zeros(nbre_shell)
J_abs_shell_tab = np.zeros(nbre_shell)
J_shell_tab = np.zeros(nbre_shell)
J_abs_int_tab = np.zeros(nbre_shell)
J_int_tab = np.zeros(nbre_shell)
radius_shell = np.zeros(nbre_shell)
if direction_J == True:
	dir_theta = np.zeros(nbre_shell)
	dir_phi = np.zeros(nbre_shell)

count=0
mass_int_shell=0
J_abs_int=0
J_int=0

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

	#Elargissement du tableau shell_radius_cell pour la division dans la partie suivante
	shell_radius_cell_redim = np.zeros((len(shell_radius_cell),3))
	for i in range(3):
		shell_radius_cell_redim[:,i]=shell_radius_cell



	#-------------------------------------
	#-------------------------------------
	#Evalutation de V_rot dans la coquille
	#-------------------------------------
	#-------------------------------------
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
	J_shell_tab[s] = J_shell

	#-----------------------------------------------------------------------
	#Valeur absolue des moments cinetiques par rapport au centre et stockage
	#-----------------------------------------------------------------------
	J_abs_shell_tab[s] = J_shell_abs


	#-------------------------------
	#Evaluation de V_rot et stockage
	#-------------------------------
	V_rot = J_shell / (R_shell*shell_mass)
	V_rot_stock[s] = V_rot


	#---------------------------------------------------------------
	#Moments cinetiques integre et moment absolu integre et stockage
	#---------------------------------------------------------------
	J_int += J_shell
	J_abs_int += J_shell_abs

	J_int_tab[s] = J_int
	J_abs_int_tab[s] = J_abs_int


	#---------------------------------
	#Directions de moments et stockage
	#---------------------------------
	if direction_J == True:
		(rayon, theta, phi) = coord.cartesian_to_spherical(moment_shell_x,moment_shell_y,moment_shell_z)
		phi = phi.value
		if theta.value >= 0.:
			theta = np.pi/2. - theta.value
		else:
			theta = np.pi/2. + theta.value
		dir_theta[s] = theta *180/np.pi  #theta en degre
		dir_phi[s] = phi *180/np.pi  #phi en degre



#--------------------------------
#Trace des differentes quantites
#--------------------------------
plt.semilogy(radius_shell,V_rot_stock, marker='.',color='midnightblue', label=ur'$V_{rot}$')
plt.xlabel(ur'Position de la coquille ($AU$)')
plt.ylabel(ur'Vitesses ($m.s^{-1}$)')
plt.legend(loc='best')

plt.figure()
plt.semilogy(radius_shell,J_shell_tab, marker='.',color='midnightblue', label=ur'$J$')
plt.semilogy(radius_shell,J_abs_shell_tab, marker='.',color='chocolate', label=ur'$J_{abs}$')
plt.xlabel(ur'Position de la coquille ($AU$)')
plt.ylabel(ur'Moment ($kg.m^2.s^{-1}$)')
plt.legend(loc='best')

plt.figure()
plt.plot(radius_shell,J_int_tab, marker='.',color='midnightblue', label=ur'$J_{tot}$')
plt.plot(radius_shell,J_abs_int_tab, marker='.',color='chocolate', label=ur'$J_{abs, tot}$')
plt.xlabel(ur'Position de la coquille ($AU$)')
plt.ylabel(ur'Moment ($kg.m^2.s^{-1}$)')
plt.legend(loc='best')

if direction_J == True:
	plt.figure()
	plt.plot(radius_shell,dir_theta, marker='.',color='midnightblue', label=ur'$\theta$')
	plt.plot(radius_shell,dir_phi, marker='.',color='chocolate', label=ur'$\phi$')
	plt.xlabel(ur'Position de la coquille ($AU$)')
	plt.ylabel(ur'Angle (degrés)')
	plt.legend(loc='best')


'''
#-----------------------------------------------------------------------
#Sauvegarde du tableau contenant les moments cinetiques de chaque output
#-----------------------------------------------------------------------
if save_final == True:
	hdu = fits.PrimaryHDU(J_tot_tab)
	hdu.writeto(file_save)

	hdu = fits.PrimaryHDU(J_abs_tot_tab)
	hdu.writeto(file_save_2)
'''


#plt.savefig('/gpfs/data1/averliat/analyses/B335_noturb_norot_hydro_pert_asym_aleatoire_forte/Shell_moment_cinetique_tot_et_abs_output12.pdf', bbox_inches='tight')

#plt.savefig('/gpfs/data1/averliat/analyses/B335_noturb_norot_hydro_pert_asym_aleatoire_forte/Vitesse_rotation_output12.pdf', bbox_inches='tight')


#print("Avancement : ",round((n+1)/N_cluster*100.,2),"%", end="\r")
