import numpy as np

import pymses
import sys
import glob as glob
import os 

from pymses.filters import CellsToPoints
from pymses.utils import constants as cst
from pymses.analysis import Camera, raytracing, slicing
from pymses.analysis import ScalarOperator, FractionOperator, MaxLevelOperator
import matplotlib.pyplot as plt
plt.ion()

#from astropy.io import fits

import os




#-------------------------------------------------------
#Entree le nom de la simulation et le numero de l'output
#-------------------------------------------------------
simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire_tres_forte'
owner = 'averliat'
num_output = 39

save = True
dir_save = 'Coupe_vitesse_integree'
radius_zoom = 2


selon_x = True
selon_y = True
selon_z = True

#Pour fermer toutes les figures avec matplotlib : plt.close('all')



#--------------------------------------------------------------
#Chemin de la simulation et du dossier de sauvegarde des images
#--------------------------------------------------------------
path='/gpfs/data1/'+owner+'/'+simu+'/'
path_save='/gpfs/data1/averliat/analyses/'+simu+'/'+dir_save+'/'

if save==True:
	if os.path.isdir(path_save) == False:
		os.mkdir(path_save)



#-------------------
#Lecture de l'output
#-------------------
ro=pymses.RamsesOutput(path,num_output)
lbox_pc = ro.info['unit_length'].express(cst.pc)

amr = ro.amr_source(["rho","vel","P","phi","g"])

cell_source = CellsToPoints(amr)
cells = cell_source.flatten()

pos = cells.points
rho = cells["rho"]



#------------------------------------------------------------
#Facteur de conversion des unites de code en unites physiques
#------------------------------------------------------------
lbox=ro.info['boxlen']
lbox_m = ro.info['unit_length'].express(cst.m)
lbox_cm = ro.info['unit_length'].express(cst.cm)

factor_vel_km_s = ro.info['unit_velocity'].express(cst.km/cst.s)



#---------------------------------------------------------
#Definition du centre des images et de leur niveau de zoom
#---------------------------------------------------------
#Position du "centre" de la simulation = endroit le plus dense
arg_centre = np.argmax(rho)
center = [pos[:,0][arg_centre],pos[:,1][arg_centre],pos[:,2][arg_centre]]

zoom_v=[0.045, 0.015, 0.005, 0.005/3.]
radius=zoom_v[radius_zoom-1]
#radius=float(zoom_v[np.where(zoom_v==radius_zoom)])     #0.015#0.005 #Niveau de zoom correspondant au niveau '3' des images de "pipeline_image_unique.py"



#-----------------
#-----------------
#Calcul des cartes
#-----------------
#-----------------

if selon_x==True:
	#--------------------------------------------
	#Calcul de la carte ou l'on regarde suivant x
	#--------------------------------------------
	cam_x = Camera(center=center,line_of_sight_axis='x',region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector='z',map_max_size=512)

	rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_density"])
	rt = raytracing.RayTracer(amr,ro.info,rho_op)
	datamap = rt.process(cam_x, surf_qty=True)
	map_col = np.log10(datamap.map.T*lbox_cm)

	Vx_op = ScalarOperator(lambda dset: dset["vel"][...,0]*dset["rho"] ,  ro.info["unit_velocity"])
	rt = raytracing.RayTracer(amr,ro.info,Vx_op)
	datamap_vx = rt.process(cam_x, surf_qty=True)
	map_Vx = datamap_vx.map.T / datamap.map.T * factor_vel_km_s


	plt.figure()
	plt.imshow(map_col,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')   
	plt.xlabel('$y$ (pc)')     
	plt.ylabel('$z$ (pc)')
	cbar = plt.colorbar()                                                                           
	cbar.set_label(r'$log(N) \, cm^{-2}$')  
	if save==True:
		plt.savefig(path_save+'dens_x_'+str(radius_zoom)+'_'+str(num_output)+'.pdf', bbox_inches='tight')


	plt.figure()
	plt.imshow(map_Vx,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')
	plt.xlabel('$y$ (pc)')     
	plt.ylabel('$z$ (pc)')
	cbar = plt.colorbar()         
	cbar.set_label(r'$v_x \, (km.s^{-1})$')  

	if save==True:
		plt.savefig(path_save+'vel_x_'+str(radius_zoom)+'_'+str(num_output)+'.pdf', bbox_inches='tight')


if selon_y==True:
	#--------------------------------------------
	#Calcul de la carte ou l'on regarde suivant y
	#--------------------------------------------
	cam_y = Camera(center=center,line_of_sight_axis='y',region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector='x',map_max_size=512)

	rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_density"])
	rt = raytracing.RayTracer(amr,ro.info,rho_op)
	datamap = rt.process(cam_y, surf_qty=True)
	map_col = np.log10(datamap.map.T*lbox_cm)

	Vy_op = ScalarOperator(lambda dset: dset["vel"][...,0]*dset["rho"] ,  ro.info["unit_velocity"])
	rt = raytracing.RayTracer(amr,ro.info,Vy_op)
	datamap_vy = rt.process(cam_y, surf_qty=True)
	map_Vy = datamap_vy.map.T / datamap.map.T * factor_vel_km_s


	plt.figure()
	im = plt.imshow(map_col,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')   
	plt.xlabel('$z$ (pc)')     
	plt.ylabel('$x$ (pc)')
	cbar=plt.colorbar()                                                                                   
	cbar.set_label(r'$log(N) \, cm^{-2}$')
	if save==True:
		plt.savefig(path_save+'dens_y_'+str(radius_zoom)+'_'+str(num_output)+'.pdf', bbox_inches='tight')  


	plt.figure()
	plt.imshow(map_Vy,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')
	plt.xlabel('$z$ (pc)')     
	plt.ylabel('$x$ (pc)')
	cbar = plt.colorbar()          
	cbar.set_label(r'$v_y \, (km.s^{-1})$')  
	if save==True:
		plt.savefig(path_save+'vel_y_'+str(radius_zoom)+'_'+str(num_output)+'.pdf', bbox_inches='tight')


if selon_z==True:
	#--------------------------------------------
	#Calcul de la carte ou l'on regarde suivant z
	#--------------------------------------------
	cam_z = Camera(center=center,line_of_sight_axis='z',region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector='y',map_max_size=512)

	rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_density"])
	rt = raytracing.RayTracer(amr,ro.info,rho_op)
	datamap = rt.process(cam_z, surf_qty=True)
	map_col = np.log10(datamap.map.T*lbox_cm)

	Vz_op = ScalarOperator(lambda dset: dset["vel"][...,0]*dset["rho"] ,  ro.info["unit_velocity"])
	rt = raytracing.RayTracer(amr,ro.info,Vz_op)
	datamap_vz = rt.process(cam_z, surf_qty=True)
	map_Vz = datamap_vz.map.T / datamap.map.T * factor_vel_km_s


	plt.figure()
	im = plt.imshow(map_col,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')   
	plt.xlabel('$x$ (pc)')     
	plt.ylabel('$y$ (pc)')
	cbar=plt.colorbar()                                                                                     
	cbar.set_label(r'$log(N) \, cm^{-2}$')  
	if save==True:
		plt.savefig(path_save+'dens_z_'+str(radius_zoom)+'_'+str(num_output)+'.pdf', bbox_inches='tight')


	plt.figure()
	plt.imshow(map_Vz,extent=[(-radius+center[0])*lbox,(radius+center[0])*lbox,(-radius+center[1])*lbox,(radius+center[1])*lbox],origin='lower')
	plt.xlabel('$x$ (pc)')     
	plt.ylabel('$y$ (pc)')
	cbar = plt.colorbar()          
	cbar.set_label(r'$v_z \, (km.s^{-1})$')  
	if save==True:
		plt.savefig(path_save+'vel_z_'+str(radius_zoom)+'_'+str(num_output)+'.pdf', bbox_inches='tight')




















