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




if __name__ == '__main__':
#-------------------------------------------------------
#Entree le nom de la simulation et le numero de l'output
#-------------------------------------------------------
    simu1 = 'cluster_test_ir4_nsink8'

    owner = 'averliat_alfven'
    output_min1=12
    output_max1=51

    colorsink1='red'
    colorgas1='green'
    colortot1='blue'

    legend1='  avec jets'
    


    simu2 = 'cluster_test_ir4_nsink8_ssjet'

    output_min2=12
    output_max2=21

    colorsink2='darkred'
    colorgas2='darkgreen'
    colortot2='darkblue'

    legend2='  sans jets'





def SFR_calculation(simu,owner,output_min,output_max,colorsink,colorgas,colortot,legend):
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

    path_t0=path





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
    mass_sinks_tab=np.zeros(output_max+1-output_min,dtype=object)
    mass_gas_tab=np.zeros(output_max+1-output_min)
    numero_output=np.zeros(output_max+1-output_min)
    simulation_time=np.zeros(output_max+1-output_min)
    ind_output=0

    for num_output in range(output_min,output_max+1):
        m_sinks=None
        if os.path.isfile(path+'output_'+str(num_output).zfill(5)+'/sink_'+str(num_output).zfill(5)+'.csv') == True:
            sinks=read_sink_csv(num_output,path)
            if len(sinks) != 0:
                m_sinks=sinks[:,1]  #masse des sinks en masses solaires
            else:
                m_sinks=0.



        mass_sinks_tab[ind_output]=np.array(m_sinks)
        numero_output[ind_output]=num_output
        

        
        #-------------------
        #lecture de l'output
        #-------------------
        ro=pymses.RamsesOutput(path,num_output)
        
        amr = ro.amr_source(["rho"])

        cell_source = CellsToPoints(amr)
        cells = cell_source.flatten(verbose=False)
        dx = cells.get_sizes()

        rho = cells["rho"]



        #------------------------------
        #conversion en unites physiques
        #------------------------------
        lbox_au = ro.info['unit_length'].express(cst.au)
        factor_dist= ro.info['unit_length'].express(cst.m)
        factor_vel = ro.info['unit_velocity'].express(cst.m/cst.s)
        factor_rho = ro.info['unit_density'].express(cst.kg/(cst.m)**3)
        
        factor_time_Myr = ro.info['unit_time'].express(cst.Myr)

        simulation_time[ind_output] = ro.info['time']*factor_time_Myr
        
        dx *= factor_dist
        rho *= factor_rho


        mass_tot = np.sum( rho * dx**3 ) / 1.9889e+30  #En Msun
        mass_gas_tab[ind_output]=mass_tot
        

        ind_output+=1


    mass_sinks_sum_tab=np.zeros(output_max+1-output_min)
    for i in range(len(mass_sinks_tab)):
        mass_sinks_sum_tab[i]=np.sum(mass_sinks_tab[i])



    #---------------------------------------------------
    #Trace de la masse des sinks en fonction de l'output
    #---------------------------------------------------
    plt.figure(20)
    plt.plot(numero_output,mass_sinks_sum_tab,linestyle='solid',marker='.',color=colorsink)
    plt.xlabel(r"Numero de l'output")
    plt.ylabel(r"Masse des sinks ($M_\odot$)")


    plt.figure(21)
    plt.plot(simulation_time,mass_sinks_sum_tab,linestyle='solid',marker='.',color=colorsink)
    plt.xlabel(r"Temps (Myr)")
    plt.ylabel(r"Masse des sinks ($M_\odot$)")


    plt.figure(22)
    plt.plot(numero_output,mass_sinks_sum_tab+mass_gas_tab,linestyle='solid',marker='.',color=colortot)
    plt.xlabel(r"Numero de l'output")
    plt.ylabel(r"Masse totale ($M_\odot$)")


    plt.figure(23)
    plt.semilogy(simulation_time,mass_gas_tab,linestyle='solid',marker='.',color=colorgas,label='Gas'+legend)
    plt.semilogy(simulation_time,mass_sinks_sum_tab,linestyle='solid',marker='.',color=colorsink,label='Sinks'+legend)
    plt.semilogy(simulation_time,mass_sinks_sum_tab+mass_gas_tab,linestyle='solid',marker='.',color=colortot,label='Gas + sinks'+legend)
    plt.xlabel(r"Temps (Myr)")
    plt.ylabel(r"Masse ($M_\odot$)")
    plt.legend(loc='best')





if __name__ == '__main__':
    SFR_calculation(simu1,owner,output_min1,output_max1,colorsink1,colorgas1,colortot1,legend1)
    SFR_calculation(simu2,owner,output_min2,output_max2,colorsink2,colorgas2,colortot2,legend2)
