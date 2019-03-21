
# -*- coding: utf-8 -*-

#----------------------------------------------------------------------------------------
#Script lisant un output d'une simulation et sortant le temps associe.
#En faisant cela pour des outputs de simulations differentes, on peut les comparer
#au meme temps
#----------------------------------------------------------------------------------------


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


#import pipeline_temps_0_simulation as t_0
import pipeline_temps_0_simulation_par_masse_coeur as pmdc



#--------------------------------------------
#Caracteristiques de la simulation consideree
#--------------------------------------------
simu_ref = 'B335_noturb_norot_hydro_pert_asym_aleatoire50_vhr'
num_output_ref = 99

simu_a_comparer = 'B335_noturb_norot_hydro_pert_asym_aleatoire50_hr'


recalage_temporel = True
seuil_rho = 1e11   #Utile que si recalage_temporel==True
critere_mass_seuil = 10


#Auteur des simulations:
owner_ref = 'averliat'
owner_a_comparer = 'averliat'

#Construction des chemins:
if owner_ref == 'averliat':
	path_simu_ref = '/gpfs/data1/'+owner_ref+'/'+simu_ref+'/'
if owner_ref =='phennebe':
	path_simu_ref = '/drf/projets/capucine/'+owner_ref+'/'+simu_ref+'/'

if owner_a_comparer == 'averliat':
	path_simu_a_comparer = '/gpfs/data1/'+owner_a_comparer+'/'+simu_a_comparer+'/'
if owner_a_comparer == 'phennebe':
	path_simu_a_comparer = '/drf/projets/capucine/'+owner_a_comparer+'/'+simu_a_comparer+'/'

path_simu_ref_analyses = '/gpfs/data1/averliat/analyses/'+simu_ref+'/'
path_simu_a_comparer_analyses = '/gpfs/data1/averliat/analyses/'+simu_a_comparer+'/'

num_output_max = 1000  #Ne pas modifier a priori, peut depassser le nombre max d'output pour tout lire
sortie_output_ref=1  #Voir dans "pipeline_temps_0_simulation" a quoi ca sert  #ATTENTION, PLUS UTILISER, LAISSER A 1
sortie_output_a_comparer=1  #Voir dans "pipeline_temps_0_simulation" a quoi ca sert  #IDEM



#---------------------------------------------------------------------------------
#Calcul des t_0 de chaque simulation ou lecture des fichiers si calculs deja faits
#---------------------------------------------------------------------------------
if recalage_temporel == True:
	if os.path.isfile(path_simu_ref_analyses+'t0_seuil_mass_core_'+str(seuil_rho)+'_'+str(critere_mass_seuil)+'.txt') == False:
		ref = np.array(pmdc.temps_0_simu(path_simu_ref, seuil_rho, critere_mass_seuil, sortie_output=sortie_output_ref))
		np.savetxt(path_simu_ref_analyses+'t0_seuil_mass_core_'+str(seuil_rho)+'_'+str(critere_mass_seuil)+'.txt', ref)
	else:
		ref=np.loadtxt(path_simu_ref_analyses+'t0_seuil_mass_core_'+str(seuil_rho)+'_'+str(critere_mass_seuil)+'.txt')


	if os.path.isfile(path_simu_a_comparer_analyses+'t0_seuil_mass_core_'+str(seuil_rho)+'_'+str(critere_mass_seuil)+'.txt') == False:
		a_comparer = np.array(pmdc.temps_0_simu(path_simu_a_comparer, seuil_rho, critere_mass_seuil, sortie_output=sortie_output_a_comparer))
		np.savetxt(path_simu_a_comparer_analyses+'t0_seuil_mass_core_'+str(seuil_rho)+'_'+str(critere_mass_seuil)+'.txt', a_comparer)
	else:
		a_comparer=np.loadtxt(path_simu_a_comparer_analyses+'t0_seuil_mass_core_'+str(seuil_rho)+'_'+str(critere_mass_seuil)+'.txt')


	#Erreur si on cherche a etudier un fichier inferieur au t0 de la simulation de reference
	if num_output_ref < ref[0]:
		print
		print
		print("=================================================")
		print("=================================================")
		print("/!\   output_ref  <  output_ref_t0   /!\ ")
		print("=================================================")
		print("=================================================")
		System.exit (0)

else:
	ref=(1,0)
	a_comparer=(1,0)



#-----------------------------------------
#Lecture du temps de l'output de reference
#-----------------------------------------
ro=pymses.RamsesOutput(path_simu_ref, num_output_ref)

factor_time_Myr_ref = ro.info['unit_time'].express(cst.Myr)
simulation_time_ref = ro.info['time']*factor_time_Myr_ref - ref[1]



#----------------------------------------------------------------
#Recherche des output correspondant dans la simulation a comparer
#----------------------------------------------------------------
output_manquant=0
for l in range(num_output_max):
	current_output = (l+1) + int(a_comparer[0])
	try:
		ro=pymses.RamsesOutput(path_simu_a_comparer, current_output)
	except:
		output_manquant += 1
		if output_manquant==100:
			print
			print
			print("=========================================================================================")
			print("=========================================================================================")
			print("Pas assez d'output pour la simulation : "+path_simu_a_comparer )
			print("=========================================================================================")
			print("=========================================================================================")
			break
		continue

	output_prec = current_output - (output_manquant+1)
	output_manquant=0

	factor_time_Myr = ro.info['unit_time'].express(cst.Myr)
	simulation_time = ro.info['time']*factor_time_Myr - a_comparer[1]

	if simulation_time >= simulation_time_ref:
		output_sup = current_output
		simulation_time_sup = simulation_time
		
		output_inf = output_prec #current_output-1
		ro=pymses.RamsesOutput(path_simu_a_comparer, output_inf)
		factor_time_Myr = ro.info['unit_time'].express(cst.Myr)
		simulation_time_inf = ro.info['time']*factor_time_Myr - a_comparer[1]
		#output_inf = output_prec #current_output-1
		
		diff_simulation_time_min_sup = abs(simulation_time_sup - simulation_time_ref)
		diff_simulation_time_min_inf = abs(simulation_time_ref - simulation_time_inf)

		print
		print
		print
		print("==========================================================================")
		print("==========================================================================")
		print("Pour la simulation : "+path_simu_ref )
		print
		print("       output_ref = " +str(num_output_ref))
		print("       time_ref = " +str(simulation_time_ref) +" Myr")
		print("==========================================================================")
		print("==========================================================================")
		print
		
		if diff_simulation_time_min_sup < diff_simulation_time_min_inf:
			print
			print("==========================================================================")
			print("==========================================================================")
			print("Pour la simulation : "+path_simu_a_comparer )
			print
			print("    output_inf = "+str(output_inf) )
			print("    time_inf = "+str(simulation_time_inf)+ " Myr")
			print
			print("    output_sup = "+str(output_sup)+"    <---" )
			print("    time_sup = "+str(simulation_time_sup)+ " Myr")
			print("==========================================================================")
			print("==========================================================================")
			break

		if diff_simulation_time_min_sup >= diff_simulation_time_min_inf:
			print
			print("==========================================================================")
			print("==========================================================================")
			print("Pour la simulation : "+path_simu_a_comparer )
			print
			print("    output_inf = "+str(output_inf)+"    <---" )
			print("    time_inf = "+str(simulation_time_inf)+ " Myr")
			print
			print("    output_sup = "+str(output_sup) )
			print("    time_sup = "+str(simulation_time_sup)+ " Myr")
			print("==========================================================================")
			print("==========================================================================")
			break






