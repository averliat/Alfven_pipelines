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
import h5py



if __name__ == '__main__':
#-------------------------------------------------------------------
#Chemin de la sauvegarde du tableau contenant les moments cinetiques
#-------------------------------------------------------------------
    solo=False

    simu='B335_noturb_norot_hydro_pert_asym_aleatoire_shr_bigbox_50pourc'
    tag='50_shr'

    num_min = 1
    num_max = 59



def trace_moment_cin(simu,tag,num_min,num_max,t1=0.,legend='',marker='.'):
#------------------
#Differents chemins
#------------------
    path_analyses='/home/averliat/these/analyses/'+simu+'/'
    file_save = 'Moment_cinetique_tot_et_abs_output_'+str(num_min)+'_'+str(num_max)+'.hdf5'



#--------------------------
#Ouverture de la sauvegarde
#--------------------------
    h5f = h5py.File(path_analyses+file_save, 'r')

    J_tot_tab=h5f['J_tot_tab'][:]
    J_abs_tot_tab=h5f['J_abs_tot_tab'][:]
    numero_output=h5f['numero_output'][:]
    time_tab=h5f['time_tab'][:]

    h5f.close()



#------------------
#Debut de la figure
#------------------
    '''
    plt.figure(1)
    plt.semilogy(numero_output, J_tot_tab, marker=marker, label=ur'Moment cinétique total')
    plt.semilogy(numero_output, J_abs_tot_tab, marker = '.', color= 'chocolate', label=ur'$\int |$ moments cinétiques $| $')

    plt.xlabel(ur"Output")
    plt.ylabel(ur"Moment cinétique ($kg.m^2.s^-1$)")
    plt.legend(loc='best')
    '''

    plt.figure(1)
    plt.semilogy(numero_output, J_tot_tab, marker=marker, label=legend)
    plt.xlabel(ur"Output")
    plt.ylabel(ur"Moment cinétique total ($kg.m^2.s^-1$)")
    plt.legend(loc='best')

    plt.figure(2)
    plt.semilogy(numero_output, J_abs_tot_tab, marker=marker, label=legend)
    plt.xlabel(ur"Output")
    plt.ylabel(ur"$\int |$ moments cinétiques $| $ ($kg.m^2.s^-1$)")
    plt.legend(loc='best')

    plt.figure(3)
    plt.semilogy(time_tab, J_tot_tab, marker=marker, label=legend)
    plt.xlabel(ur"Time ($Myr$)")
    plt.ylabel(ur"Moment cinétique total ($kg.m^2.s^-1$)")
    plt.legend(loc='best')

    plt.figure(4)
    plt.semilogy(time_tab, J_abs_tot_tab, marker=marker, label=legend)
    plt.xlabel(ur"Time ($Myr$)")
    plt.ylabel(ur"$\int |$ moments cinétiques $| $ ($kg.m^2.s^-1$)")
    plt.legend(loc='best')

    plt.figure(5)
    plt.semilogy(time_tab-t1, J_tot_tab, marker=marker, label=legend)
    plt.xlabel(ur"Corrected time ($Myr$)")
    plt.ylabel(ur"Moment cinétique total ($kg.m^2.s^-1$)")
    plt.legend(loc='best')

    plt.figure(6)
    plt.semilogy(time_tab-t1, J_abs_tot_tab, marker=marker, label=legend)
    plt.xlabel(ur"Corrected time ($Myr$)")
    plt.ylabel(ur"$\int |$ moments cinétiques $| $ ($kg.m^2.s^-1$)")
    plt.legend(loc='best')
    
    plt.show()



if __name__ == '__main__':
    if solo==True:
        trace_moment_cin(simu,tag,num_min,num_max)

    elif solo==False:
        #tag=['50_shr','50_hhr','50_vhr','50_emhr','50_mhr','50_hr']
        #legend=['8, 10*40', '7, 10*80', '7, 10*40', '7, 10*20', '7,   8*40', '7,   8*20']
        #marker=['.','v','p','s','P','*']

        #output_max=[721, 285, 254, 420, 250, 144]

        #t1=[0.09912+6.95e-6, 0.09708-5.41e-4, 0.09653+5.72e-4-6.75e-5, 0.09243+2.72e-5, 0.09068-5.16e-4, 0.09068+1.84e-6]

        #base_simu='B335_noturb_norot_hydro_pert_asym_aleatoire'


        #for i in range(len(tag)):
        #    i=len(tag)-i-1
        #    trace_moment_cin(base_simu+tag[i],tag[i],1,output_max[i],t1[i],tag[i]+' ,  '+legend[i],marker[i])

        #trace_moment_cin('B335_noturb_norot_hydro_hr2','hr2',1,346,0.09725-1.16e-4,'hr2 ,  7,  10*40','X')

        #trace_moment_cin('B335_noturb_rot1_hydro_pert_asym_aleatoire50_vhr','rot1',1,58,0,'rot1 ,  7,  10*40','x')

        #trace_moment_cin('B335_noturb_norot_hydro_pert_asym_aleatoire50_shr_bigbox','50_shr_bigbox',1,108,0,'50_shr_bigbox ,  8,  10*40','p')

        #trace_moment_cin('B335_noturb_norot_hydro_pert_asym_aleatoire50_hr_niMHD','50_hr_niMHD',1,2720,0,'50_hr_niMHD ,  7,  8*20','p')




        tag=['0pourc','10pourc','20pourc','50pourc']
        legend=['8, 10*40','8, 10*40','8, 10*40','8, 10*40']
        marker=['.','v','p','s']

        output_max=[41,57,60,59]

        #t1=[0.09912+6.95e-6, 0.09708-5.41e-4, 0.09653+5.72e-4-6.75e-5, 0.09243+2.72e-5, 0.09068-5.16e-4, 0.09068+1.84e-6]
        t1=[0,0,0,0]

        base_simu='B335_noturb_norot_hydro_pert_asym_aleatoire_shr_bigbox_'


        for i in range(len(tag)):
            i=len(tag)-i-1
            trace_moment_cin(base_simu+tag[i],tag[i],1,output_max[i],t1[i],tag[i]+' ,  '+legend[i],marker[i])
