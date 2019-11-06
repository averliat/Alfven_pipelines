import numpy as np
import matplotlib.pyplot as plt
import trace_taille_disque_par_pdf_modal_sans_MHD as ttdppm



plt.style.use("pdf")
plt.style.use("aanda_modif")



cmappts = plt.get_cmap('hot')#autumn')
colorsrot = [cmappts(i) for i in np.linspace(0.3,0.7,4)]

cmappts = plt.get_cmap('summer')
colorsnorot = [cmappts(i) for i in np.linspace(0.2,0.8,3)]

cmappts = plt.get_cmap(r'cool')
colorsMHD = [cmappts(i) for i in np.linspace(0.3,0.7,2)]

#tag=['50','50','20','20','50','50','50','20','10']
#tag2=['_rot1','_rot0.25','_rot1','_rot0.25','_MHD_lr','_MHD_lr_rot1','','','']
tag=['50','50','20','20','50','20','10']
tag2=['_rot1','_rot0.25','_rot1','_rot0.25','','','']
legend=[r'$\varepsilon=50\%$~,~~$\beta=1\%$',
        r'$\varepsilon=50\%$~,~~$\beta=0.25\%$',
        r'$\varepsilon=20\%$~,~~$\beta=1\%$',
        r'$\varepsilon=20\%$~,~~$\beta=0.25\%$',
        #r'MHD,~$\varepsilon=50\%$~',
        #r'MHD,~$\varepsilon=50\%$~,~~$\beta=1\%$',
        r'$\varepsilon=50\%$',
        r'$\varepsilon=20\%$',
        r'$\varepsilon=10\%$']




marker=['+','.','x','p','P','>','<','.','+','x']

#output_max=[30,70,20,40,367,138,480,400,440]
output_max=[30,70,20,40,480,400,440]

#colors=[colorsrot[0],colorsrot[1],colorsrot[2],colorsrot[3],colorsMHD[0],colorsMHD[1],colorsnorot[0],colorsnorot[1],colorsnorot[2]]
colors=[colorsrot[0],colorsrot[1],colorsrot[2],colorsrot[3],colorsnorot[0],colorsnorot[1],colorsnorot[2]]

#output_frag=[24,63,14,30,'None','None',480,400,440]
output_frag=[24,63,14,30,480,400,440]

t1=[0,0,0,0,0,0,0,0,0]

base_simu='B335_noturb_norot_hydro_pert_asym_aleatoire_bigbox_'


for j in range(len(tag)):
    i=len(tag)-j-1
    ttdppm.trace_taille_disque(base_simu+tag[i]+'pourc_sink_seuil_haut'+tag2[i],tag[i],tag[i]+tag2[i],'None',output_max[i],colors[i],j,tag,t1[i],legend[i],'.',output_frag=output_frag[i])



save=True
path_save='/home/averliat/these/analyses/article1_figures/taille_disk/'
if save==True:
    plt.tight_layout(pad=0.1) #pad en inch si besoin
    plt.savefig(path_save+'Taille_disque_all_sans_MHD.pdf')#, bbox_inches='tight')
