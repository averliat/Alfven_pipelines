

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
from matplotlib.colors import Normalize

#from astropy.io import fits

import os
import pipeline_temps_0_simulation as t_0

import scipy.optimize as opt




#-------------------------------------------------------
#Entree le nom de la simulation et le numero de l'output
#-------------------------------------------------------
simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire50_hr_sink'
owner = 'averliat_alfven'
num_output = 825


save = False
dir_save = 'Diagramme_PV'

radius_zoom = 4  #6 pour avoir ~1500AU et voir "l'enveloppe" 
view_diagram = 'y'
direction_diagram = 'hori'  #'vert'=verticale, 'hori'=horizontale

taille_pts=1 #AU, pour lissage


v_proj = True

title_time=True
title_time_cor=False
seuil_rho = 1e-10


vmin_vel = None#-0.5
vmax_vel = None#-vmin_vel 

vmin_dens = None
vmax_dens = None

#Pour fermer toutes les figures avec matplotlib : plt.close('all')



#--------------------------------------------------------------
#Chemin de la simulation et du dossier de sauvegarde des images
#--------------------------------------------------------------
if owner=='averliat':
	path='/mnt/magmist/magmist/simu_B335_averliat/'+simu+'/'
if owner=='phennebe':
	path='/drf/projets/capucine/'+owner+'/'+simu+'/'
if owner=='sapanais':
	path='/dsm/anais/storageA/magmist/'+simu+'/'
if owner=='averliat_alfven':
	path='/drf/projets/alfven-data/averliat/'+simu+'/'
	
path_save='/home/averliat/these/analyses/'+simu+'/'+dir_save+'/'
path_analyse='/home/averliat/these/analyses/'+simu+'/'

#if simu == 'B335_noturb_norot_hydro_pert_asym_aleatoire50_vhr':
#	path_t0='/mnt/magmist/magmist/simu_B335_averliat/'+simu+'/'
#else:
#	path_t0=path
path_t0=path

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
dx = cells.get_sizes()

pos = cells.points
rho = cells["rho"]



#------------------------------------------------------------
#Facteur de conversion des unites de code en unites physiques
#------------------------------------------------------------
lbox=ro.info['boxlen']
lbox_m = ro.info['unit_length'].express(cst.m)
lbox_au = ro.info['unit_length'].express(cst.au)
lbox_cm = ro.info['unit_length'].express(cst.cm)
factor_time_yr = ro.info['unit_time'].express(cst.year)
factor_vel_km_s = ro.info['unit_velocity'].express(cst.km/cst.s)
factor_rho = ro.info['unit_density'].express(cst.kg/(cst.m)**3)
factor_dist= ro.info['unit_length'].express(cst.m)

simulation_time = ro.info['time']*factor_time_yr


'''
#---------------------------------------------------------------------------------
#Calcul des t_0 de chaque simulation ou lecture des fichiers si calculs deja faits
#---------------------------------------------------------------------------------
if title_time_cor == True:
	if os.path.isfile(path_analyse+'t0_seuil_rho_'+str(seuil_rho)+'.txt') == False:
		ref = np.array(t_0.temps_0_simu(path_t0, seuil_rho, sortie_output=1))
		np.savetxt(path_analyse+'t0_seuil_rho_'+str(seuil_rho)+'.txt', ref)
	else:
		ref=np.loadtxt(path_analyse+'t0_seuil_rho_'+str(seuil_rho)+'.txt')


	#Erreur si on cherche a etudier un fichier inferieur au t0 de la simulation de reference
	if num_output < ref[0]:
		print
		print
		print("=================================================")
		print("=================================================")
		print("/!\   output_ref  <  output_ref_t0   /!\ ")
		print("=================================================")
		print("=================================================")
		title_time_cor=0
'''


#---------------------------------------------------------
#Definition du centre des images et de leur niveau de zoom
#---------------------------------------------------------
#Position du "centre" de la simulation = endroit le plus dense
if radius_zoom==5:
	center = [0.5,0.5,0.5]
else:
	arg_centre = np.argmax(rho)
	center = [pos[:,0][arg_centre],pos[:,1][arg_centre],pos[:,2][arg_centre]]
	

zoom_v=[0.045, 0.015, 0.005, 0.005/3., lbox_pc*4, 0.025]
radius=zoom_v[radius_zoom-1]
#radius=float(zoom_v[np.where(zoom_v==radius_zoom)])     #0.015#0.005 #Niveau de zoom correspondant au niveau '3' des images de "pipeline_image_unique.py"


#---------------------------------------------------------------------------------------------------------------
#Definition prise sur https://stackoverflow.com/questions/20144529/shifted-colorbar-matplotlib/20146989#20146989 
#pour avoir le centre de la colorbar a 0
#---------------------------------------------------------------------------------------------------------------
class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))



#-----------------
#-----------------
#Calcul des cartes
#-----------------
#-----------------


#------------------
#Calcul de la carte
#------------------
if view_diagram=='x':
	up_vector='z'
	V_op_number=0
elif view_diagram=='y':
	up_vector='x'
	V_op_number=1
elif view_diagram=='z':
	up_vector='y'
	V_op_number=2

#view_diagram=[1,-0.5,0]
#up_vector=[0.8,0.8,1]


cam = Camera(center=center,line_of_sight_axis=view_diagram,region_size=[2.*radius,2.*radius],distance=radius,far_cut_depth=radius,up_vector=up_vector,map_max_size=512)

rho_op = ScalarOperator(lambda dset: dset["rho"] ,  ro.info["unit_density"])
rt = raytracing.RayTracer(amr,ro.info,rho_op)
datamap = rt.process(cam, surf_qty=True)
map_col = np.log10(datamap.map.T*lbox_cm)

if v_proj == True:
	V_op = ScalarOperator(lambda dset: dset["vel"][...,V_op_number]*dset["rho"] ,  ro.info["unit_velocity"])
	rt = raytracing.RayTracer(amr,ro.info,V_op)
	datamap_v = rt.process(cam, surf_qty=True)
	map_V = datamap_v.map.T / datamap.map.T * factor_vel_km_s



plt.figure()
im = plt.imshow(map_col,extent=[(-radius+center[0])*lbox_au,(radius+center[0])*lbox_au,(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au],origin='lower', vmin=vmin_dens, vmax=vmax_dens)   
plt.xlabel('$x$ (AU)')     
plt.ylabel('$y$ (AU)')
cbar=plt.colorbar()                                                                  
cbar.set_label(r'$log(N) \, \, (cm^{-2})$')  
if title_time==True:
	plt.title('Time = '+str(int(simulation_time))+' years')
if title_time_cor==True:
	plt.title('Time = '+str(int(simulation_time))+' years \n Corrected time = '+str(int(simulation_time - ref[1]*1e6))+' years')
if radius_zoom==5:
	plt.xlim([0,lbox_au])
	plt.ylim([0,lbox_au])

if direction_diagram=='hori':  #Pour tracer les lignes horizontales ou verticales sur les cartes
	plt.plot([(-radius+center[0])*lbox_au,(radius+center[0])*lbox_au],[center[1]*lbox_au,center[1]*lbox_au], color ='m', linewidth=1.1, linestyle="--")
elif direction_diagram=='vert':
	plt.plot([center[0]*lbox_au,center[0]*lbox_au],[(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au], color ='m', linewidth=1.1, linestyle="--")

if save==True:
	plt.savefig(path_save+'dens_'+view_diagram+'_'+str(radius_zoom)+'_'+str(num_output)+'_pour_diagPV.pdf', bbox_inches='tight')



if v_proj == True:
	plt.figure()
	norm = MidpointNormalize(midpoint=0)  #Pour avoir le centre de la colormap a 0
	plt.imshow(map_V,extent=[(-radius+center[0])*lbox_au,(radius+center[0])*lbox_au,(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au],origin='lower',cmap='RdBu_r',norm=norm, vmin=vmin_vel, vmax=vmax_vel)
	plt.xlabel('$x$ (AU)')     
	plt.ylabel('$y$ (AU)')
	cbar = plt.colorbar()          
	cbar.set_label(r'$v_{los} \, (km.s^{-1})$')  
	if title_time==True:
		plt.title('Time = '+str(int(simulation_time))+' years')
	if title_time_cor==True:
		plt.title('Time = '+str(int(simulation_time))+' years \n Corrected time = '+str(int(simulation_time - ref[1]*1e6))+' years')
	if radius_zoom==5:
		plt.xlim([0,lbox_au])
		plt.ylim([0,lbox_au])

	if direction_diagram=='hori':
		plt.plot([(-radius+center[0])*lbox_au,(radius+center[0])*lbox_au],[center[1]*lbox_au,center[1]*lbox_au], color ='limegreen', linewidth=1.1, linestyle="--")
	elif direction_diagram=='vert':
		plt.plot([center[0]*lbox_au,center[0]*lbox_au],[(-radius+center[1])*lbox_au,(radius+center[1])*lbox_au], color ='limegreen', linewidth=1.1, linestyle="--")

	if save==True:
		plt.savefig(path_save+'vel_'+view_diagram+'_'+str(radius_zoom)+'_'+str(num_output)+'_pour_diagPV.pdf', bbox_inches='tight')



#Diagramme initial
if direction_diagram=='hori':
	diagpv_ord = map_V[256,:] #ligne horizontale, map_V[:,256] pour ligne verticale
	diagpv_abs = np.linspace((-radius+center[0])*lbox_au,(radius+center[0])*lbox_au,512) #mettre [1] si verticale

if direction_diagram=='vert':
	diagpv_ord = map_V[:,256]
	diagpv_abs = np.linspace((-radius+center[1])*lbox_au,(radius+center[1])*lbox_au,512)

plt.figure()
plt.plot(diagpv_abs,diagpv_ord)
plt.xlabel('Distance (AU)')     
plt.ylabel(r'$v_{los} \, (km.s^{-1})$')
#plt.loglog(diagpv_abs[0:220],diagpv_ord[0:220])



#Diagramme avec branches retournees
plt.figure()
if direction_diagram=='hori':
	test_sens_rot_droite=np.sum(map_V[256,256:256+10])
	test_sens_rot_gauche=np.sum(map_V[256,256-10:256])
	if test_sens_rot_gauche>test_sens_rot_droite:
		diagpv_ord_1 = map_V[256,0:256][::-1]  #Part du centre, vers la gauche 
		diagpv_abs_1 = np.linspace(radius*lbox_au,0,256)[::-1] #mettre [1] si verticale
		diagpv_ord_2 = -map_V[256,257:512] #Part du centre, vers la droite
		diagpv_abs_2 = np.linspace(0,radius*lbox_au,255) #mettre [1] si verticale
	if test_sens_rot_gauche<test_sens_rot_droite:
		diagpv_ord_2 = -map_V[256,0:256][::-1]  #Part du centre, vers la gauche
		diagpv_abs_2 = np.linspace(radius*lbox_au,0,256)[::-1] #mettre [1] si verticale
		diagpv_ord_1 = map_V[256,257:512] #ligne horizontale, map_V[:,256] pour ligne verticale
		diagpv_abs_1 = np.linspace(0,radius*lbox_au,255) #mettre [1] si verticale

if direction_diagram=='vert':
	test_sens_rot_haut=np.sum(map_V[256:256+10,256])
	test_sens_rot_bas=np.sum(map_V[256-10:256,256])
	if test_sens_rot_bas>test_sens_rot_haut:
		diagpv_ord_1 = map_V[0:256,256][::-1]  #Part du centre, vers le bas
		diagpv_abs_1 = np.linspace(radius*lbox_au,0,256)[::-1] #mettre [1] si verticale
		diagpv_ord_2 = -map_V[257:512,256] #ligne horizontale, map_V[:,256] pour ligne verticale
		diagpv_abs_2 = np.linspace(0,radius*lbox_au,255) #mettre [1] si verticale
	if test_sens_rot_bas<test_sens_rot_haut:
		diagpv_ord_2 = -map_V[0:256,256][::-1]  #Part du centre, vers le bas
		diagpv_abs_2 = np.linspace(center[1]*lbox_au,0,256)[::-1] #mettre [1] si verticale
		diagpv_ord_1 = map_V[257:512,256] #ligne horizontale, map_V[:,256] pour ligne verticale
		diagpv_abs_1 = np.linspace(0,center[1]*lbox_au,255) #mettre [1] si verticale

plt.plot(diagpv_abs_1,diagpv_ord_1,color='r')
plt.plot(diagpv_abs_2,diagpv_ord_2,color='b')
plt.xlabel('Distance au centre (AU)')     
plt.ylabel(r'$v_{los} \, (km.s^{-1})$')
		
plt.figure()
plt.loglog(diagpv_abs_1,diagpv_ord_1,color='r')
plt.loglog(diagpv_abs_2,diagpv_ord_2,color='b')
plt.xlabel('Distance au centre (AU)')
plt.ylabel(r'$v_{los} \, (km.s^{-1})$')

plt.figure()
plt.loglog(diagpv_abs_1[np.argmax(diagpv_ord_1):],diagpv_ord_1[np.argmax(diagpv_ord_1):],color='r')
plt.loglog(diagpv_abs_2[np.argmax(diagpv_ord_2):],diagpv_ord_2[np.argmax(diagpv_ord_2):],color='b')
plt.xlabel('Distance au centre (AU)')
plt.ylabel(r'$v_{los} \, (km.s^{-1})$')



#Lissage:
#On veut des pixels plus gros que 10AU
#On a 256 points pour decrire une taille radius*lbox_au
nmbre_points_max = np.int(radius*lbox_au/taille_pts)

for i in range(0,256):
	reste=256%(256-i)
	if (reste==0 and nmbre_points_max>=(256-i)):
		n_final = 256-i
		n_par_group=256/(256-i)
		#print n_final
		#print n_par_group
		break

diagpv_abs_1_moy = np.zeros(n_final)
diagpv_abs_2_moy = np.zeros(n_final)
diagpv_ord_1_moy = np.zeros(n_final)
diagpv_ord_2_moy = np.zeros(n_final)


#REVOIR CA EN HAUT QUAND 1 et 2 SONT INVERSES
test=np.copy(diagpv_abs_2)
for i in range(0,n_final):
	diagpv_abs_1_moy[i] = np.sum(diagpv_abs_1[i*n_par_group:i*n_par_group+n_par_group])/n_par_group
	diagpv_ord_1_moy[i] = np.sum(diagpv_ord_1[i*n_par_group:i*n_par_group+n_par_group])/n_par_group

for i in range(0,n_final-1):
	diagpv_abs_2_moy[i] = np.sum(diagpv_abs_2[i*n_par_group:i*n_par_group+n_par_group])/n_par_group
	diagpv_ord_2_moy[i] = np.sum(diagpv_ord_2[i*n_par_group:i*n_par_group+n_par_group])/n_par_group
	test[i*n_par_group:i*n_par_group+n_par_group]=0.

diagpv_abs_2_moy[n_final-1] = np.sum(diagpv_abs_2[(n_final-1)*n_par_group:(n_final-1)*n_par_group+n_par_group-1])/(n_par_group-1)
diagpv_ord_2_moy[n_final-1] = np.sum(diagpv_ord_2[(n_final-1)*n_par_group:(n_final-1)*n_par_group+n_par_group-1])/(n_par_group-1)
test[(n_final-1)*n_par_group:(n_final-1)*n_par_group+n_par_group]=0.


plt.figure()
plt.plot(diagpv_abs_1_moy,diagpv_ord_1_moy,color='r',marker='.')
plt.plot(diagpv_abs_2_moy,diagpv_ord_2_moy,color='b',marker='.')
plt.xlabel('Distance au centre (AU)')     
plt.ylabel(r'$v_{los} \, (km.s^{-1})$')
plt.title('Diagramme PV lisse')

plt.figure()
plt.loglog(diagpv_abs_1_moy[np.argmax(diagpv_ord_1_moy):],diagpv_ord_1_moy[np.argmax(diagpv_ord_1_moy):],color='r',marker='.')
plt.loglog(diagpv_abs_2_moy[np.argmax(diagpv_ord_2_moy):],diagpv_ord_2_moy[np.argmax(diagpv_ord_2_moy):],color='b',marker='.')
plt.xlabel('Distance au centre (AU)')     
plt.ylabel(r'$v_{los} \, (km.s^{-1})$')
plt.title('Diagramme PV lisse')


'''
if direction_diagram=='vert':
	diagpv_ord = map_V[:,256]
	diagpv_abs = np.linspace((-radius+center[1])*lbox_au,(radius+center[1])*lbox_au,512)

plt.figure()
plt.plot(diagpv_abs,diagpv_ord)
#plt.loglog(diagpv_abs[0:220],diagpv_ord[0:220])
'''



#-------------------
#Fit du diagramme PV
#-------------------
def func_fit(x,a,b):
    y=a*x+b
    return y


abs_max=45 #AU
abs_min=15 #AU
#abs_fit_1=np.log10(diagpv_abs_1[np.argmax(diagpv_ord_1):np.max(np.where(diagpv_abs_1<abs_max))])
#ord_fit_1=np.log10(diagpv_ord_1[np.argmax(diagpv_ord_1):np.max(np.where(diagpv_abs_1<abs_max))])
#abs_fit_2=np.log10(diagpv_abs_2[np.argmax(diagpv_ord_2):np.max(np.where(diagpv_abs_2<abs_max))])
#ord_fit_2=np.log10(diagpv_ord_2[np.argmax(diagpv_ord_2):np.max(np.where(diagpv_abs_2<abs_max))])

abs_fit_1=np.log10(diagpv_abs_1[np.min(np.where(diagpv_abs_1>abs_min)):np.max(np.where(diagpv_abs_1<abs_max))])
ord_fit_1=np.log10(np.abs((diagpv_ord_1[np.min(np.where(diagpv_abs_1>abs_min)):np.max(np.where(diagpv_abs_1<abs_max))])))
abs_fit_2=np.log10(diagpv_abs_2[np.min(np.where(diagpv_abs_2>abs_min)):np.max(np.where(diagpv_abs_2<abs_max))])
ord_fit_2=np.log10(diagpv_ord_2[np.min(np.where(diagpv_abs_2>abs_min)):np.max(np.where(diagpv_abs_2<abs_max))])

p0= 0., 1.

#Ajustement
popt1,pcov1=opt.curve_fit(func_fit,abs_fit_1,ord_fit_1,p0=p0)
popt2,pcov2=opt.curve_fit(func_fit,abs_fit_2,ord_fit_2,p0=p0)

y1=popt1[0]*abs_fit_1 + popt1[1]
y1=10**y1
x1=10**abs_fit_1

y2=popt2[0]*abs_fit_2 + popt2[1]
y2=10**y2
x2=10**abs_fit_2

plt.figure(6)
plt.loglog(x1,y1,color='orangered')
plt.loglog(x2,y2,color='dodgerblue')




#-------------------------
#overplot de np.sqrt(GM/r)
#-------------------------
abs_max_int=diagpv_abs_1[np.argmax(diagpv_ord_1)] #AU
pos *= lbox_au  #AU

radius_cells_au = np.sqrt( np.sum( (pos[:,:]-np.array(center)*lbox_au)**2 , axis=1 ) )  #Tableau contenant le rayon de chaque cellule

mask = (radius_cells_au < abs_max_int/2.)

rho_int = rho[mask]
dx_int = dx[mask]

dx_int *= factor_dist
rho_int *= factor_rho

#Masse contenue dans la partie non keplerienne
M_int = np.sum( rho_int* dx_int**3 ) #mettre = pour prendre en compte l'interieur

G_grav = 6.67408e-11  #m3.kg-1.s-2

diag_pv_theor_1=np.sqrt(G_grav*M_int/(diagpv_abs_1[np.argmax(diagpv_ord_1):]/lbox_au*lbox_m))/1000.  #km/s

plt.figure(12)
plt.plot(diagpv_abs_1[np.argmax(diagpv_ord_1):], diag_pv_theor_1, marker='.', color='g')

plt.figure(14)
plt.loglog(diagpv_abs_1[np.argmax(diagpv_ord_1):], diag_pv_theor_1, marker='.', color='g')

