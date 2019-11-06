# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import article1_trace_shell_moment_cin_observateurs2 as tsmco2
import pipeline_temps_0_simulation as t_0
import os



'''
tag=['50_shr','50_hhr','50_vhr','50_emhr','50_mhr','50_hr']
legend=['8, 10*40', '7, 10*80', '7, 10*40', '7, 10*20', '7,   8*40', '7,   8*20']
marker=['.','v','p','s','P','*']

output_max=[685, 275, 254, 225, 130, 144]

t1=[0.09912+6.95e-6, 0.09708-5.41e-4, 0.09653+5.72e-4-6.75e-5, 0.09243+2.72e-5, 0.09068-5.16e-4, 0.09068+1.84e-6]


base_simu='B335_noturb_norot_hydro_pert_asym_aleatoire'

width_bar=[6,10,1,5,5,1]
output_min_global=1


hist=False
diff_abs=True
diff_relat=True



for i in range(len(tag)):
    i=len(tag)-i-1
    tsmco2.trace_moment_cin_disque(base_simu+tag[i],tag[i],'None',1,output_max[i],width_bar[i],hist,diff_abs,diff_relat,legend[i],marker[i],t1[i])
'''

#tsmco2.trace_moment_cin_disque('B335_noturb_norot_hydro_hr2','hr2','None',1,210,2,hist,diff_abs,diff_relat,'7,  10*40','X',0.09725-1.16e-4)

'''
plt.figure(10)
plt.xlim((-0.02,0.015))
plt.ylim((-0.1e46,2.6e46))

plt.figure(11)
plt.xlim((-0.02,0.015))
plt.ylim((-0.05,0.85))
'''




#tsmco2.trace_moment_cin_disque('B335_noturb_norot_hydro_pert_asym_aleatoire50_shr_bigbox','50_shr_bigbox','None',1,108,1,True,True,True,'8,  10*40','.',0.10093)

#tsmco2.trace_moment_cin_disque('B335_noturb_norot_hydro_pert_asym_aleatoire50_vhr','50_vhr','None',1,254,4,True,True,True,'7,  10*40','.',0.09653+5.72e-4-6.75e-5)

#tsmco2.trace_moment_cin_disque('B335_noturb_norot_hydro_pert_asym_aleatoire50_hr_niMHD','50_hr_niMHD','None',1,2720,20,True,True,True,'7,  8*20','.',0.10919)



tag=['10','20','30','40','50','60']
legend=['8, 10*40','8, 10*40','8, 10*40','8, 10*40','8, 10*40','8, 10*40','8, 10*40']
marker=['v','p','s','*','<','>']

output_max=[99, 93, 100, 97, 110, 97]

output_frag=[45,65,68,76,80,69]

#t1=[0.09912+6.95e-6, 0.09708-5.41e-4, 0.09653+5.72e-4-6.75e-5, 0.09243+2.72e-5, 0.09068-5.16e-4, 0.09068+1.84e-6]
#t1=[0,0,0,0]
base_simu='B335_noturb_norot_hydro_pert_asym_aleatoire_shr_bigbox_'

t1=[]
seuil_rho = 1e-10
for i in range(len(tag)):
    path_analyses='/home/averliat/these/analyses/'+base_simu+tag[i]+'pourc/'
    if os.path.isfile(path_analyses+'t0_seuil_rho_'+str(seuil_rho)+'.txt') == False:
        ref = np.array(t_0.temps_0_simu(path_t0, seuil_rho, sortie_output=1))
        np.savetxt(path_analyses+'t0_seuil_rho_'+str(seuil_rho)+'.txt', ref)
    else:
        ref=np.loadtxt(path_analyses+'t0_seuil_rho_'+str(seuil_rho)+'.txt')
    t1.append(ref[1])



width_bar=[1,1,1,1,1,1]
output_min_global=1


hist=False
diff_abs=False
diff_relat=False
moyenne_glissante=True



for i in range(len(tag)):
    i=len(tag)-i-1
    tsmco2.trace_moment_cin_disque(base_simu+tag[i]+'pourc',tag[i],'None',1,output_max[i],width_bar[i],hist,diff_abs,diff_relat,legend[i],marker[i],t1[i],moyenne_glissante=moyenne_glissante,save=False,output_frag=output_frag[i])



plt.tight_layout(pad=0.25) #pad en inch si besoin
plt.savefig('/home/averliat/these/analyses/article1_figures/Comparaison_difference_moment_cin_reel_analytique_et_moment_dans_disque_legend_a_droite_sans_legend_jusqua_frag.pdf')
