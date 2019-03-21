import numpy as np
import matplotlib.pyplot as plt
import trace_masse_disque_par_pdf as tmdpp



'''
tag=['50_shr','50_hhr','50_vhr','50_emhr','50_mhr','50_hr']
legend=['8, 10*40', '7, 10*80', '7, 10*40', '7, 10*20', '7,   8*40', '7,   8*20']
marker=['.','v','p','s','P','*']

output_max=[685, 275, 254, 225, 130, 144]

t1=[0.09912+6.95e-6, 0.09708-5.41e-4, 0.09653+5.72e-4-6.75e-5, 0.09243+2.72e-5, 0.09068-5.16e-4, 0.09068+1.84e-6]  #Myr, pour recalage temporel et avoir toutes les courbes a peu pres superposees


base_simu='B335_noturb_norot_hydro_pert_asym_aleatoire'



for i in range(len(tag)):
    i=len(tag)-i-1
    tmdpp.trace_masse_disque(base_simu+tag[i],tag[i],'None',output_max[i],t1[i],legend[i],'.')



tmdpp.trace_masse_disque('B335_noturb_norot_hydro_hr2','hr2','None',346,0.09725-1.16e-4,'7,  10*40','.')


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
    tmdpp.trace_masse_disque(base_simu1+tag[i]+base_simu2,tag[i],'None',output_max[i],t1[i],legend[i],'.')



