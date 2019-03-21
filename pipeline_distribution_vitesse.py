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



#-------------------------------------------------------
#Entrer le nom de la simulation et le numero de l'output
#-------------------------------------------------------
simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire_forte'
owner = 'averliat'
num_output = 12


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
pres = cells["P"]



#------------------------------
#Conversion en unites physiques
#------------------------------
#lbox_m = ro.info['unit_length'].express(cst.m)
factor_dist= ro.info['unit_length'].express(cst.m)
factor_vel = ro.info['unit_velocity'].express(cst.m/cst.s)
factor_rho = ro.info['unit_density'].express(cst.kg/(cst.m)**3)
factor_pres = ro.info['unit_pressure'].express(cst.kg/(cst.m)**1/(cst.s)**2)

dx *= factor_dist
pos *= factor_dist
vel *= factor_vel
rho *= factor_rho
pres *= factor_pres



#-------------------------
#Definition de la coquille
#-------------------------
#Position du "centre" de la simulation = endroit le plus dense
arg_centre = np.argmax(rho)
center = [pos[:,0][arg_centre],pos[:,1][arg_centre],pos[:,2][arg_centre]]


radius_cell = np.sqrt( np.sum( (pos[:,:]-center)**2 , axis=1 ) )  #Tableau contenant le rayon de chaque cellule




AU = 149597870700 #m

radius_min = 0*AU
radius_max = 2000*AU
dim_shell = 1*AU
nbre_shell = int((radius_max-radius_min)/dim_shell)

cs_mean_stock = np.zeros(nbre_shell)
V_RMS_stock = np.zeros(nbre_shell)
v_infall_stock = np.zeros(nbre_shell)
radius_shell = np.zeros(nbre_shell)
vel_rad_mean_stock = np.zeros(nbre_shell)

count=0
mass_int_shell=0

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
	shell_pres = pres[mask]
	shell_dx = dx[mask]
	shell_radius_cell = radius_cell[mask]

	#Masse contenue dans la coquille
	shell_mass = np.sum( shell_rho* shell_dx**3 )

	#Elargissement du tableau shell_radius_cell pour la division dans la partie suivante
	shell_radius_cell_redim = np.zeros((len(shell_radius_cell),3))
	for i in range(3):
		shell_radius_cell_redim[:,i]=shell_radius_cell



	#-------------------------------------
	#Evalutation de V_rms dans la coquille
	#-------------------------------------
	#Vecteurs cellule dans la coquille <-- centre (tableau)
	pos_moins_centre = (shell_pos - center)  #Mettre un moins pour vecteurs cellule --> centre (change les signes de vect_unit, vel_rad, vel_rad_mean, vel_rad_mean_corrigee mais ne change pas V_RMS)

	#Vecteurs unitaires cellule <-- centre (tableau)
	vect_unit = pos_moins_centre / shell_radius_cell_redim


	#Vitesses radiales (tableau)
	vel_rad = np.sum( shell_vel*vect_unit, axis = 1 )

	#Vitesse radiale moyenne ponderee par la masse <--> vitesse d'infall (scalaire)
	vel_rad_mean = np.sum(vel_rad *shell_rho*shell_dx**3) /shell_mass 

	#Vitesses radiales corrigees de la vitesse d'infall (tableau)
	vel_rad_corrigee = vel_rad - vel_rad_mean


	#Norme des vitesses au carre (tableau)
	norm_vel_carre = np.sum(shell_vel**2, axis=1)

	#Norme des vitesses orthoradiales au carre (tableau)
	norm_vel_orthorad_carre = norm_vel_carre - vel_rad**2 

	#Norme des vitesses corrigee au carre (tableau)
	v_corrigee_norm_carre = vel_rad_corrigee**2 + norm_vel_orthorad_carre


	#Calcul de V_RMS
	V_RMS = np.sqrt(  np.sum( v_corrigee_norm_carre* shell_rho*shell_dx**3 ) /shell_mass )


	#Masse de tout l'interieur de la coquille
	mass_int_shell += np.sum(shell_rho*shell_dx**3)
	#mask_mass_int_shell = (radius_cell <= radius_out)
	#mass_int_shell = np.sum( rho[mask_mass_int_shell]*dx[mask_mass_int_shell]**3 )

	#Vitesse de chute libre
	v_infall = np.sqrt( 2*6.67e-11*mass_int_shell / ((radius_in+radius_out)/2.) )



	#-------------------------------------------------
	#Evalutation de la vitesse du son dans la coquille
	#-------------------------------------------------
	cs_tab = np.sqrt( 1.*shell_pres/shell_rho )  #gamma=1
	cs_mean = np.sum(cs_tab*shell_rho*shell_dx**3) /shell_mass


	
	#------------------------------------
	#Sauvegarde des differentes quantites
	#------------------------------------
	cs_mean_stock[s] = cs_mean
	V_RMS_stock[s] = V_RMS
	v_infall_stock[s] = v_infall
	radius_shell[s] = (radius_in + radius_out)/2. /AU #En AU pour le trace (abscisse)
	vel_rad_mean_stock[s] = vel_rad_mean



#--------------------------------
#Trace des differentes quantites
#--------------------------------
plt.semilogy(radius_shell,V_RMS_stock, marker='.',color='midnightblue', label=ur'$V_{RMS}$')
plt.semilogy(radius_shell,cs_mean_stock, marker='.',color='chocolate', label=ur'$c_s$')
plt.semilogy(radius_shell,v_infall_stock, marker='.', color= 'mediumorchid', label=ur'$v_{infall}$')
plt.semilogy(radius_shell,np.abs(vel_rad_mean_stock), marker='.', color= 'darkmagenta', label=ur'$v_{rad}$')
plt.xlabel(ur'Position de la coquille ($AU$)')
plt.ylabel(ur'Vitesses ($m.s^{-1}$)')
plt.legend(loc='best')



















