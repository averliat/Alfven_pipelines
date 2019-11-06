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
simu = 'jet_test_18'#test_18'#'B335_noturb_norot_hydro_pert_asym_aleatoire_lllr_bigbox_10pourc_sink'

owner = 'averliat_alfven'
output_min=1
output_max = 463

color='r'

#save = True
#dir_save = 'mass_sinks'



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
    
#path_save='/home/averliat/these/analyses/'+simu+'/'+dir_save+'/'
path_analyse='/home/averliat/these/analyses/'+simu+'/'

path_t0=path

#if save==True:
#    if os.path.isdir(path_save) == False:
#        os.mkdir(path_save)




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
mass_sinks_tab=np.zeros(output_max+1-output_min)
numero_output=np.zeros(output_max+1-output_min)
ind_output=0
for num_output in range(output_min,output_max+1):
    m_sinks=None
    if os.path.isfile(path+'output_'+str(num_output).zfill(5)+'/sink_'+str(num_output).zfill(5)+'.csv') == True:
        sinks=read_sink_csv(num_output,path)
        if len(sinks) != 0:
            m_sinks=sinks[:,1]  #masse des sinks en masses solaires
        else:
            m_sinks=0.



    mass_sinks_tab[ind_output]=m_sinks
    numero_output[ind_output]=num_output
    ind_output+=1



#---------------------------------------------------
#Trace de la masse des sinks en fonction de l'output
#---------------------------------------------------
plt.figure(20)
plt.plot(numero_output,mass_sinks_tab,linestyle='None',marker='.',color=color)
plt.xlabel(r"Numero de l'output")
plt.ylabel(r"Masse de la sink ($M_\odot$)")







