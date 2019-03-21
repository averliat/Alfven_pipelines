# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import article1_trace_shell_moment_cin_observateurs as tsmco
import pipeline_temps_0_simulation as t_0
import os

plt.style.use("pdf")
plt.style.use("aanda_modif")




'''
tag=['50_shr','50_hhr','50_vhr','50_emhr','50_mhr','50_hr','50_shr_bigbox']
legend=['8, 10*40', '7, 10*80', '7, 10*40', '7, 10*20', '7,   8*40', '7,   8*20', '8, 10*40']
marker=['.','v','p','s','P','*','p']

output_max=[685, 275, 254, 225, 130, 144, 108]

t1=[0.09912+6.95e-6, 0.09708-5.41e-4, 0.09653+5.72e-4-6.75e-5, 0.09243+2.72e-5, 0.09068-5.16e-4, 0.09068+1.84e-6, 0.0]

base_simu='B335_noturb_norot_hydro_pert_asym_aleatoire'

width_bar=[1,1,1,1,1,1,1]

hist=False
diff_abs=True
diff_relat=True



for i in range(len(tag)):
    i=len(tag)-i-1
    tsmco.trace_moment_cin_all(base_simu+tag[i],tag[i],1,output_max[i],width_bar[i],hist,diff_abs,diff_relat,legend[i],marker[i],t1[i])
'''

#tsmco.trace_moment_cin_all('B335_noturb_norot_hydro_hr2','hr2',1,210,2,hist,diff_abs,diff_relat,'7,  10*40','X',0.09725-1.16e-4)


#plt.figure(10)
#plt.xlim((-0.02,0.015))
#plt.ylim((-0.1e46,2.6e46))
##plt.title("Différence entre la valeur absolue de l'erreur sur les moments et le moment dans le disque")

#plt.figure(11)
#plt.xlim((-0.02,0.015))
#plt.ylim((-0.05,0.85))
##plt.title("Différence relative entre la valeur absolue de l'erreur sur les moments et le moment dans le disque")


#tsmco.trace_moment_cin_all('B335_noturb_norot_hydro_pert_asym_aleatoire50_shr_bigbox','50_shr_bigbox',1,108,1,False,True,True,'8,  10*40','X',0.10093)

#tsmco.trace_moment_cin_all('B335_noturb_norot_hydro_pert_asym_aleatoire50_vhr','50_vhr',1,254,1,False,True,True,'7,  10*40','p',0.09653+5.72e-4-6.75e-5)

#tsmco.trace_moment_cin_all('B335_noturb_norot_hydro_pert_asym_aleatoire50_hr_niMHD','50_hr_niMHD',1,2720,20,True,True,True,'7,  8*20','.',0.10919)






tag=['0pourc','10pourc','20pourc','30pourc','40pourc','50pourc','60pourc']
legend=[r'0 \%', r'10 \%', r'20 \%', r'30 \%', r'40 \%', r'50 \%', r'60 \%']
marker=['.','v','p','s','*','<','>']

output_max=[41, 99, 93, 100, 97, 110, 97]

#t1=[0.09912+6.95e-6, 0.09708-5.41e-4, 0.09653+5.72e-4-6.75e-5, 0.09243+2.72e-5, 0.09068-5.16e-4, 0.09068+1.84e-6, 0.0]
#t1=[0,0,0,0,0,0,0]
base_simu='B335_noturb_norot_hydro_pert_asym_aleatoire_shr_bigbox_'

t1=[]
seuil_rho = 1e-10
for i in range(len(tag)):
    path_analyses='/home/averliat/these/analyses/'+base_simu+tag[i]+'/'
    if os.path.isfile(path_analyses+'t0_seuil_rho_'+str(seuil_rho)+'.txt') == False:
        ref = np.array(t_0.temps_0_simu(path_t0, seuil_rho, sortie_output=1))
        np.savetxt(path_analyses+'t0_seuil_rho_'+str(seuil_rho)+'.txt', ref)
    else:
        ref=np.loadtxt(path_analyses+'t0_seuil_rho_'+str(seuil_rho)+'.txt')
    t1.append(ref[1])



width_bar=[1,1,1,1,1,1,1]

hist=False
diff_abs=False
diff_relat=True
moyenne_glissante=True
article1=False

save=False #Laisser a False, le save est modifie dans la boucle
save_plusieurs=True


for i in range(len(tag)):
    if (i==len(tag)-1 and save_plusieurs==True):
        save=True
    i=len(tag)-i-1
    tsmco.trace_moment_cin_all(base_simu+tag[i],tag[i],1,output_max[i],width_bar[i],hist,diff_abs,diff_relat,legend[i],marker[i],t1[i],moyenne_glissante=moyenne_glissante,save=save,article1=article1,save_plusieurs=save_plusieurs)
