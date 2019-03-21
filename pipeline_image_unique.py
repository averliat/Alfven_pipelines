##extraction pipeline
import numpy as np
import module_extract as me
import pymses
import sys
import glob as glob
import os 



#-----------------------------
#Simulation et output a sortir
#-----------------------------
path='/gpfs/data1/averliat/B335_noturb_norot_hydro_pert_asym_aleatoire_tres_forte/'
path_out='/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_pert_asym_aleatoire_tres_forte/figures/'

num_output = 14

coldens_only = True  #True pour avoir que l'image de densite colonne



#-----------------------------
#Creation du dossier si besoin
#-----------------------------
if os.path.isdir(path_out) == False:
	os.mkdir(path_out)



#-----------------------------------------------------
#Definition des niveaux de zoom et creation des images
#-----------------------------------------------------
print(num_output)
#zoom_v=[0.5/16.,0.5/64.,0.5/256.,0.5/1024.]
zoom_v=[0.045, 0.015, 0.005, 0.005/3.]
#zoom_v=[0.5] #Pour avoir toute la boite de la simu

print path

me.make_image_zoom(path,num_output,zoom_v,path_out,sinks=False,force=False,center_dmax=True,ps=True,col_dens_only=coldens_only)





