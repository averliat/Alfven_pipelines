import numpy as np
import matplotlib.pyplot as plt
import article1_trace_taille_disque_par_pdf as ttdpp



'''
tag=['50_shr','50_hhr','50_vhr','50_emhr','50_mhr','50_hr']
legend=['8, 10*40', '7, 10*80', '7, 10*40', '7, 10*20', '7,   8*40', '7,   8*20']
marker=['.','v','p','s','P','*']

output_max=[685, 275, 254, 690, 130, 144]

t1=[0.09912+6.95e-6, 0.09708-5.41e-4, 0.09653+5.72e-4-6.75e-5, 0.09243+2.72e-5, 0.09068-5.16e-4, 0.09068+1.84e-6]  #Myr, pour recalage temporel et avoir toutes les courbes a peu pres superposees


base_simu='B335_noturb_norot_hydro_pert_asym_aleatoire'


#simu=[]
#for i in range(len(tag)):
#    simu.append(base_simu+tag[i])


for i in range(len(tag)):
    i=len(tag)-i-1
    ttdpp.trace_taille_disque(base_simu+tag[i],tag[i],'None',output_max[i],t1[i],legend[i],'.')
'''


#ttdpp.trace_taille_disque('B335_noturb_norot_hydro_hr2','hr2','None',346,0.09725-1.16e-4,'7,  10*40','.')
'''


tag=['rot1','rot0.5','rot0.1','rot0.01','rot0.001']
legend=['','','','','']
marker=['.','v','p','s','P']

output_max=[26,88,270,245,266]

t1=[0.0878226,0.0871526,0.0922875,0.0957122,0.0967315]

base_simu1='B335_noturb_'
base_simu2='_hydro_pert_asym_aleatoire50_vhr'


for i in range(len(tag)):
    i=len(tag)-i-1
    ttdpp.trace_taille_disque(base_simu1+tag[i]+base_simu2,tag[i],'None',output_max[i],t1[i],legend[i],'.')
'''



#ttdpp.trace_taille_disque('B335_noturb_norot_hydro_pert_asym_aleatoire50_shr_bigbox','50_shr_bigbox','None',108,0.10093,'8,  10*40','.')

#ttdpp.trace_taille_disque('B335_noturb_norot_hydro_pert_asym_aleatoire50_vhr','50_vhr','None',254,0.09653+5.72e-4-6.75e-5,'7,  10*40','.')

#ttdpp.trace_taille_disque('B335_noturb_norot_hydro_pert_asym_aleatoire50_hr_niMHD','50_hr_niMHD','None',2720,0.10919,'7,  8*20','.')




'''
tag=['10','20','30','40','50','60']
legend=['','','','','','']
marker=['.','v','p','s','P','.']

output_max=[99,93,100,97,110,97]

t1=[0,0,0,0,0,0]

base_simu='B335_noturb_norot_hydro_pert_asym_aleatoire_shr_bigbox_'


for i in range(len(tag)):
    i=len(tag)-i-1
    ttdpp.trace_taille_disque(base_simu+tag[i]+'pourc',tag[i],'None',output_max[i],t1[i],legend[i],'.')
'''




'''
tag=['0','10','20','30','40','50','60']
legend=['','','','','','','']
marker=['+','.','v','p','s','P','.']

output_max=[264,78,202,194,181,343,187]

t1=[0,0,0,0,0,0,0]

base_simu='B335_noturb_norot_hydro_pert_asym_aleatoire_lllr_bigbox_'


for i in range(len(tag)):
    i=len(tag)-i-1
    ttdpp.trace_taille_disque(base_simu+tag[i]+'pourc_sink',tag[i],'None',output_max[i],t1[i],legend[i],'.')
'''

tag=['10','20','50','50','50','50','50']
tag2=['','','','_lr','_MHD_lr','_niMHD_lr','_rot1']
legend=['','','','','','','']
marker=['+','x','.','.','.','.','.']

output_max=[440,400,480,68,114,85,102]

t1=[0,0,0,0,0,0,0]

base_simu='B335_noturb_norot_hydro_pert_asym_aleatoire_bigbox_'


for i in range(len(tag)):
    ind=i
    i=len(tag)-i-1
    ttdpp.trace_taille_disque(base_simu+tag[i]+'pourc_sink_seuil_haut'+tag2[i],tag[i],'None',output_max[i],t1[i],legend[i],'.',ind=ind,tagtot=len(tag))


#plt.tight_layout(pad=0.1) #pad en inch si besoin
#plt.savefig('/home/averliat/these/analyses/article1_figures/Taille_disk_10_20_50.pdf')
