##extraction pipeline
import numpy as np
import module_extract as me
import pymses
import sys
import glob as glob
import os 



path='/gpfs/data1/averliat/B335_noturb_norot_hydro_pert_asym_aleatoire0.99/'
path_out='/home/users/mnt/averliat/analyses/B335_noturb_norot_hydro_pert_asym_aleatoire0.99/figures'

#'''
num_min = 25
num_max = 37

interv = 3

num_tot = num_max-num_min+1

if os.path.isdir(path_out) == False:
	os.mkdir(path_out)


for i in range(int(num_tot/interv)):
	num = num_min + i*interv
	print(num)
	#zoom_v=[0.5/16.,0.5/64.,0.5/256.,0.5/1024.]
	zoom_v=[0.045, 0.015, 0.005, 0.005/3.]  #[0.005/3., 0.005, 0.015, 0.045]
	#zoom_v=[0.005]

	print path

	me.make_image_zoom(path,num,zoom_v,path_out,sinks=False,force=False,center_dmax=True,ps=True)





'''
num = 40

#zoom_v=[0.5/16.,0.5/64.,0.5/256.,0.5/1024.]
zoom_v=[0.005/3., 0.005, 0.015, 0.045] # essayer 0.0005



print path

me.make_image_zoom(path,num,zoom_v,path_out,sinks=False,force=False,center_dmax=True,ps=True)
'''
