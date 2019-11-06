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




#-------------------------------------------------------
#Entree le nom de la simulation et le numero de l'output
#-------------------------------------------------------
simu = 'jet_fullM_rotNS_dirJ'#'B335_noturb_norot_hydro_pert_asym_aleatoire_lllr_bigbox_10pourc_sink'

owner = 'averliat_alfven'
num_output = 403



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

path_analyse='/home/averliat/these/analyses/'+simu+'/'



#-------------------
#Lecture de l'output
#-------------------
ro=pymses.RamsesOutput(path,num_output)
lbox_pc = ro.info['unit_length'].express(cst.pc)

amr = ro.amr_source(["rho","vel","P","phi","g"])

cell_source = CellsToPoints(amr)
cells = cell_source.flatten(verbose=False)
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
factor_rho_Hcc = ro.info['unit_density'].express(cst.H_cc)
factor_dist_m = ro.info['unit_length'].express(cst.m)

simulation_time = ro.info['time']*factor_time_yr



#------------------------------
#Conversion en unites physiques
#------------------------------
dx *= factor_dist_m
rho_CU=np.copy(rho)
rho_hcc = rho*factor_rho_Hcc
rho *= factor_rho



mass_tot = np.sum( rho * dx**3 ) / 1.9889e+30  #En Msun



#------------------------------------------------------------------------
#Get the properties of the particules with just the particules' positions
#Copie depuis module_extract.py
#------------------------------------------------------------------------
def read_sink_csv(num,directory,no_cvs=None):

    name = directory+'output_'+str(num).zfill(5)+'/sink_'+str(num).zfill(5)+'.csv'
    print 'read sinks from ', name

    if(no_cvs is None):
        sinks = np.loadtxt(name,delimiter=',',ndmin=2,usecols=(0,1,2,3,4,5,6,7,8)) #,9,10,11,12))
    else:
        sinks = np.loadtxt(name,ndmin=2,usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12))

    return sinks



#---------------------------------------
#Lecture du fichier sink.csv s'il existe
#---------------------------------------
#mass_sinks_tab=np.zeros(output_max+1-output_min)
#numero_output=np.zeros(output_max+1-output_min)
#ind_output=0
#for num_output in range(output_min,output_max+1):
#    if os.path.isfile(path+'output_'+str(num_output).zfill(5)+'/sink_'+str(num_output).zfill(5)+'.csv') == True:
#        sinks=read_sink_csv(num_output,path)
#        if len(sinks) != 0:
#            m_sinks=sinks[:,1]  #masse des sinks en masses solaires
#        else:
#            m_sinks=0.



#    mass_sinks_tab[ind_output]=m_sinks
#    numero_output[ind_output]=num_output
#    ind_output+=1


if os.path.isfile(path+'output_'+str(num_output).zfill(5)+'/sink_'+str(num_output).zfill(5)+'.csv') == True:
    sinks=read_sink_csv(num_output,path)
    if len(sinks) != 0:
        m_sinks=sinks[:,1]  #masse des sinks en masses solaires
    else:
        m_sinks=0.



print
print("==========================================================================")
print("    M_gaz = "+str(mass_tot)+ " Msun")
print("    M_sink = "+str(m_sinks)+ " Msun")
print("    M_tot = "+str(mass_tot+m_sinks)+ " Msun")
print("==========================================================================")
print
