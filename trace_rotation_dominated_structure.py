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

import h5py
import scipy.optimize as opt



simu = 'B335_noturb_norot_hydro_pert_asym_aleatoire_shr_bigbox_50pourc'
output = 59
dr = 1

tag='50pourc'
legend='8,  10*40'
marker='.'


#Pour les fits:
#abs_min=7 #AU
#abs_max=40 #AU


seuil_rho = 1e-10



#------------------
#Differents chemins
#------------------
path_analyses='/home/averliat/these/analyses/'+simu+'/'
file_save = 'Structure_dominee_par_rotation_output'+str(output)+'_dr'+str(dr)+'AU.hdf5'



#------------------------------------
#Chargement des differentes quantites
#------------------------------------
h5f = h5py.File(path_analyses+file_save, 'r')
shells_au=h5f['shells_au'][:]
V_r_shell_moy_tab=h5f['V_r_shell_moy_tab'][:] /1e3  #en km.s^-1
V_r_rms_shell_moy_tab=h5f['V_r_rms_shell_moy_tab'][:] /1e3  #en km.s^-1
V_paral_shell_moy_tab=h5f['V_paral_shell_moy_tab'][:] /1e3  #en km.s^-1
sigma_shell_moy_tab=h5f['sigma_shell_moy_tab'][:] /1e3  #en km.s^-1
M_cumul_tab=h5f['M_cumul_tab'][:]  #en kg
M_int_tab=h5f['M_int_tab'][:]  #en kg
simulation_time=h5f['simulation_time'][()]  #en Myr
#Msink=h5f['Msink'][()]  #en Msun
h5f.close()



'''
#-----------------------------------
#Fit du diagramme PV sans parametres
#-----------------------------------
G_grav = 6.67408e-11  #m3.kg-1.s-2
M_sun = 1.9884e30  #kg
AU = 149597870700  #m


abs_fit=np.log10(shells_au[np.min(np.where(shells_au>abs_min)):np.max(np.where(shells_au<abs_max))])
ord_fit=np.log10(np.abs((V_paral_shell_moy_tab[np.min(np.where(shells_au>abs_min)):np.max(np.where(shells_au<abs_max))])))


def func_fit1(x,a,b):
    y=a*x+b
    return y

p0= 0., 1.

#Ajustement
popt1,pcov1=opt.curve_fit(func_fit1,abs_fit,ord_fit,p0=p0)

y1=popt1[0]*abs_fit + popt1[1]
y1=10**y1
x1=10**abs_fit


M_kep1 = 1e6*AU/(G_grav*M_sun)  * 10**(2*popt1[1])


print("=======================================================")
print("Fit avec tous paramètres libres  -- courbe rouge :")
print("Exposant = "+str(popt1[0]))
print("Masse keplerienne = "+str(M_kep1)+" Msun")
print("=======================================================")




#------------------------------------------
#Fit du diagramme PV avec exposant -0.5 fixe
#------------------------------------------
def func_fit2(x,b):
    y=-0.5*x+b
    return y

p0= 1.

#Ajustement
popt2,pcov2=opt.curve_fit(func_fit2,abs_fit,ord_fit,p0=p0)

y2=-0.5*abs_fit + popt2[0]
y2=10**y2
x2=10**abs_fit


M_kep2 = 1e6*AU/(G_grav*M_sun)  * 10**(2*popt2[0])


print("Fit avec exposant -0.5 fixe  -- courbe verte :")
print("Masse keplerienne = "+str(M_kep2)+" Msun")
print("=======================================================")




#---------------
#Fonction exacte
#---------------
vel_th=np.sqrt(G_grav*Msink/shells_au) *np.sqrt(M_sun/(1e6*AU))

print("Masse de la sink = "+str(Msink)+" Msun  -- courbe mauve")
print("=======================================================")
'''




#------------------
#Debut de la figure
#------------------
cmappts = plt.get_cmap('magma')
colorspts = [cmappts(i) for i in np.linspace(0.2,0.8,4)]



'''
plt.figure()
plt.plot(shells_au, V_paral_shell_moy_tab, marker=marker, label=legend+'  - '+tag, linestyle='None')

plt.plot(x1,y1, color='orangered')
plt.plot(x2,y2, color='forestgreen')
plt.plot(shells_au, vel_th, color='darkmagenta')

plt.ylim((0,np.max(np.nan_to_num(V_paral_shell_moy_tab))*1.2))

plt.xlabel(ur'Rayon (AU)')
plt.ylabel(ur"Vitesse orthoradiale moyenne ($m.s^{-1}$)")
#plt.legend(loc='best')



plt.figure()
plt.loglog(shells_au, V_paral_shell_moy_tab, marker=marker, label=legend+'  - '+tag, linestyle='None')

plt.loglog(x1,y1,color='orangered')
plt.loglog(x2,y2, color='forestgreen')
plt.loglog(shells_au, vel_th, color='darkmagenta')

plt.xlabel(ur'Rayon (AU)')
plt.ylabel(ur"Vitesse orthoradiale moyenne ($m.s^{-1}$)")
#plt.legend(loc='best')
'''




plt.figure()
plt.semilogx(shells_au,np.abs(V_r_shell_moy_tab), marker='.',color=colorspts[0],label=ur'$<V_r>$')
plt.semilogx(shells_au,V_r_rms_shell_moy_tab, marker='.',color=colorspts[1],label=ur'$V_{r, rms}$')
plt.semilogx(shells_au,V_paral_shell_moy_tab, marker='.',color=colorspts[2],label=ur'$V_{par}$')
plt.semilogx(shells_au, sigma_shell_moy_tab, marker='.',color=colorspts[3],label=ur'$\sigma$')
plt.xlabel(ur'Rayon ($AU$)')
plt.ylabel(ur"V  $(km.s^{-1})$")
plt.legend(loc='best')




plt.figure()
plt.semilogx(shells_au, sigma_shell_moy_tab/np.abs(V_r_shell_moy_tab), marker='.',color=colorspts[0],label=ur'$\sigma$ / $<V_r>$')
plt.semilogx(shells_au, V_paral_shell_moy_tab/np.abs(V_r_shell_moy_tab), marker='.',color=colorspts[1],label=ur'$V_{par}$/$<V_r>$')
plt.semilogx(shells_au, sigma_shell_moy_tab/V_r_rms_shell_moy_tab, marker='.',color=colorspts[2],label=ur'$\sigma$/$V_{r, rms}$')
plt.semilogx(shells_au,V_paral_shell_moy_tab/V_r_rms_shell_moy_tab, marker='.',color=colorspts[3],label=ur'$V_{par}$/$V_{r, rms}$')
plt.xlabel(ur'Rayon ($AU$)')
plt.ylabel(ur"Rapport")
plt.legend(loc='best')




plt.show()



#plt.savefig('../Comparaisons/Verification_moment_physique/Difference_relative_moment_cinetique_analytique_et_reel_hr2.pdf', bbox_inches='tight')
