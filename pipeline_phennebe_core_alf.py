##extraction pipeline
import numpy as np
import module_extract as me
import pymses
import sys
import glob as glob
import os 
import module_plot_clump as mpc

import extract_disk as ed
import extract_core_profil as ecp

#path='/gpfs/data1/phennebe/RB_DENSE_CORE2/'
#num_v = [400,350,300,250,200,150,100,50]


path='/drf/projets/alfven-data/phennebe/DC_2_lowres/'
#ed.analyse_disk_all(path,path_out=None,force=True,tag='',order='<',ps=True,pdf=True,do_plot=True)

path='/drf/projets/alfven-data/phennebe/DC_2_lowres_earlysink/'
#ed.analyse_disk_all(path,path_out=None,force=True,tag='',order='<',ps=True,pdf=True,do_plot=True)

path='/drf/projets/alfven-data/phennebe/DC_2_res_vlowacc/'
#ed.analyse_disk_all(path,path_out=None,force=True,tag='',order='<',ps=True,pdf=True,do_plot=True)

path='/drf/projets/alfven-data/phennebe/DC_2_rest2_test/'
#ed.analyse_disk_all(path,path_out=None,force=True,tag='',order='<',ps=True,pdf=True,do_plot=True)

path_in_v=['/drf/projets/alfven-data/phennebe/DC_2/','/drf/projets/alfven-data/phennebe/DC_2_res_vlowacc/','/drf/projets/alfven-data/phennebe/DC_2_lowres/','/drf/projets/alfven-data/phennebe/DC_2_lowres_earlysink/','/drf/projets/alfven-data/phennebe/DC_2_rest2_test/']

#ed.prop_disk_all(path,path_out=None,force=True,tag='',gcomp=True,conv_sink=2,ps=True,pdf=True,do_plot=True,tag_selec='0')

ed.plot_disk_serie_over(path_in_v,path_out=None,force=False,tag='',ps=True,pdf=True)
#STOP



path='/drf/projets/alfven-data/phennebe/DC_2_fld_hres/'
zoom_v=[0.5/256.,0.5/512.]

#ed.prop_disk_all(path,path_out=None,force=True,tag='',gcomp=True,conv_sink=2,ps=True,pdf=True,do_plot=True,tag_selec='741')



path='/drf/projets/alfven-data/phennebe/DC_2_rest4/'
zoom_v=[0.5/256.,0.5/512.]

#num=1546
#ed.single_z_disc_prop(path,num)

#me.make_image_zoom_last(path,zoom_v,sinks=False,force=True,center_dmax=True,ps=True,prob_df=False,cpuamr=True,tag='zoom')

#ed.prop_disk_all(path,path_out=None,force=True,tag='',gcomp=True,conv_sink=2,ps=True,pdf=True,do_plot=True,tag_selec='0')



path='/drf/projets/alfven-data/phennebe/DC_2_fld_hres/'
zoom_v=[0.5/256.,0.5/512.]

#ed.prop_disk_all(path,path_out=None,force=True,tag='',gcomp=True,conv_sink=2,ps=True,pdf=True,do_plot=True,tag_selec='0')



path='/drf/projets/alfven-data/phennebe/DC_2_rest3/'
zoom_v=[0.5/256.,0.5/512.]

#ed.prop_disk_all(path,path_out=None,force=True,tag='',gcomp=True,conv_sink=2,ps=True,pdf=True,do_plot=True,tag_selec='0')



#num=1546
#ed.single_z_disc_prop(path,num)

#me.make_image_zoom_last(path,zoom_v,sinks=False,force=True,center_dmax=True,ps=True,prob_df=False,cpuamr=True,tag='zoom')



#STOP





force_disk= True
#path='/drf/projets/alfven-data/phennebe/DC_7/'
#ed.analyse_disk_all(path,path_out=None,force=force_disk,tag='',order='<',ps=True,pdf=True,do_plot=True)



#num_v=[100]
#ed.analyse_disk_serie(path,num_v,path_out=None,force=True,ps=False,tag='',order='<')

#zoom_v=[0.5/128.]

#for num in num_v:

#    me.make_image_zoom(path,num,zoom_v,sinks=False,center_dmax=True,ps=True,force=True,tag='zoom')




#STOP
path='/drf/projets/alfven-data/phennebe/DC_1/'
#ed.analyse_disk_all(path,path_out=None,force=force_disk,tag='',order='<',ps=True,pdf=True,do_plot=True)

#num=1000

#ed.single_z_disc_prop(path,num,force=True,tag='',order='<',do_plot=True,path_out=path)

#ed.plot_disc_prop(path,num,ps=True,pdf=False)


#STOP


force = True
#force=False

skip_core=True#False
#skip_core=False

movie= False #True

force_mov=False


list_path=[]

path='/drf/projets/capucine/phennebe/B335_weakB_FLD_sink_feed/'
#list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_6/'
#list_path.append((path))

path='/drf/projets/capucine/phennebe/B335_simple_lowrot/'
#list_path.append((path))

path='/drf/projets/capucine/phennebe/B335_simple/'
#list_path.append((path))

path='/drf/projets/capucine/phennebe/B335_simple_turb0.5/'
#list_path.append((path))

path='/drf/projets/capucine/phennebe/B335_simple_turb0.2/'
#list_path.append((path))


path='/drf/projets/capucine/phennebe/B335_FLD_sink/'
#list_path.append((path))

path='/drf/projets/alfven-data/phennebe/B335_fld/'
#list_path.append((path))

path='/drf/projets/alfven-data/phennebe/B335_turb_fld/'
#list_path.append((path))

path='/drf/projets/alfven-data/phennebe/B335_turb_fld2/'
#list_path.append((path))


if( not skip_core):
    for path in list_path:

        print path

#        num_v=[10,14,20,200]

        num_v=[5,9,15,30]

#        num_v=[3,4,5]

        time_l, r_l, rho_l, Vrad_l,Vpar_l,sig_l,mach_l,temp_l = ecp.analyse_core_serie(path,num_v,force=force,ps=True,pdf=True)

#STOP

#path='../Project_CORE/DC_1_occigen_res_early_restart/'
#zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.,0.5/128.]
#me.make_image_zoom_all(path,zoom_v,force=True,center_dmax=True,ps=True,prob_df=True)

force = False

disk_extract = False #True

disk_prop = False #True

force_disk= True

skip=False


image_new = True #False

list_path=[]


path='/drf/projets/alfven-data/phennebe/M30_mu7_mach2/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_2_fld_res/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_2_fld_hres/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/B335_weakB_FLD_sink_feed/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_2_res_vlowacc/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_2_lowres_earlysink/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_2_lowres/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_2_test_coefAD_updt/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_2_test_coefAD/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_2_rest2_test/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_2_rest2_test2/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_2_rest2_test3/'
list_path.append((path))



path='/drf/projets/alfven-data/phennebe/DC_2/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_2_fld/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_3/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_4/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_5/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_6/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_7/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_8/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_9/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_1/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_1_try/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_1_restart2/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_1_res/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_1_res_bsink/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_1_nosink/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_1_res_restart/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_1_res_lowacc/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_1_res_vlowacc/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_1_res_vlowacc_rest/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_1_res_lowacc_hlmax/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_1_nosink_res_sink/'
list_path.append((path))




path='/drf/projets/alfven-data/phennebe/DC_2_rest/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_2_rest2/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_2_rest3/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_2_rest4/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_2_nosink/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_2_res/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_2_res_lowacc/'
list_path.append((path))
path='/drf/projets/alfven-data/phennebe/DC_2_res_vlowacc/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_2_res_vlowacc_rest/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_2_res_vlowacc_rest2/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_2_res_vlowacc_rest2_hlmax/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_2_res_restart/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_2_res_rest_hlmax/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_3/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_4/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_5/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_6/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_7/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/DC_8/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_9/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_6_res_rest_hlmax/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_6_res_rest_vhlmax/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_6_res_restart/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_6_restart2/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_6_restart/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_6_res/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_6_res_lowacc/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_6_res_vlowacc_rest/'
list_path.append((path))

path='/drf/projets/alfven-data/phennebe/DC_10_res/'
list_path.append((path))


path='/drf/projets/alfven-data/phennebe/B335_fld/'
#list_path.append((path))


path='/drf/projets/alfven-data/phennebe/B335_fld_res/'
#list_path.append((path))
 
path='/drf/projets/alfven-data/phennebe/B335_turb_fld/'
#list_path.append((path))

path='/drf/projets/alfven-data/phennebe/B335_turb_fld2/'
#list_path.append((path))




zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.,0.5/128.]



if( not skip):
    for path in list_path:

        print path


        if(image_new):

            me.make_image_zoom_last(path,zoom_v,sinks=False,force=False,center_dmax=True,ps=True,prob_df=True)


        if(disk_prop):
            ed.prop_disk_all(path,path_out=None,force=True,tag='',gcomp=True,conv_sink=2,ps=True,pdf=True,do_plot=True)



        if(disk_extract):

            ed.analyse_disk_all(path,path_out=None,force=force_disk,tag='',order='<',ps=True,pdf=True,do_plot=True)


        if(movie):

            zoomov_v=[0.5/128.]

            me.make_image_zoom_all(path,zoomov_v,force=force_mov,col_dens_only=True,center_dmax=True)



STOP


path='../Project_CORE/DC_1_occigen/'
zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.,0.5/128.]
me.make_image_zoom_all(path,zoom_v,force=force,center_dmax=True,ps=True,prob_df=True)



path='../Project_CORE/DC_1_occigen_res/'
zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.,0.5/128.]
me.make_image_zoom_all(path,zoom_v,force=force,center_dmax=True,ps=True,prob_df=True)


path='../Project_CORE/DC_1_occigen_res2/'
zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.,0.5/128.]
me.make_image_zoom_all(path,zoom_v,force=force,center_dmax=True,ps=True,prob_df=True)






skip=True
skip=False

import extract_disk as ed

Force=False

disk_extract=False #True

if(disk_extract):

    force_disk=False


    path='/gpfs/data1/phennebe/DC_9_ps_baro2/'
    ed.analyse_disk_all(path,path_out=None,force=force_disk,tag='',order='<',ps=True,pdf=True,do_plot=True)


    path='/gpfs/data1/phennebe/DC_0b/'
    ed.analyse_disk_all(path,path_out=None,force=force_disk,tag='',order='<',ps=True,pdf=True,do_plot=True)

    path='/gpfs/data1/phennebe/DC_1_consmom/'
    ed.analyse_disk_all(path,path_out=None,force=force_disk,tag='',order='<',ps=True,pdf=True,do_plot=True)

    path='/gpfs/data1/phennebe/DC_3_ps_baro2/'
    ed.analyse_disk_all(path,path_out=None,force=force_disk,tag='',order='<',ps=True,pdf=True,do_plot=True)

    path='/gpfs/data1/phennebe/DC_4_ps_baro2/'
    ed.analyse_disk_all(path,path_out=None,force=force_disk,tag='',order='<',ps=True,pdf=True,do_plot=True)

    path='/gpfs/data1/phennebe/DC_5_ps_baro2/'
    ed.analyse_disk_all(path,path_out=None,force=force_disk,tag='',order='<',ps=True,pdf=True,do_plot=True)

    path='/gpfs/data1/phennebe/DC_6_ps_baro2/'
    ed.analyse_disk_all(path,path_out=None,force=force_disk,tag='',order='<',ps=True,pdf=True,do_plot=True)

    path='/gpfs/data1/phennebe/DC_7_ps_baro2/'
    ed.analyse_disk_all(path,path_out=None,force=force_disk,tag='',order='<',ps=True,pdf=True,do_plot=True)

    path='/gpfs/data1/phennebe/DC_8_ps_baro2/'
    ed.analyse_disk_all(path,path_out=None,force=force_disk,tag='',order='<',ps=True,pdf=True,do_plot=True)


    path='/gpfs/data1/phennebe/DC_1_consmom_hr/'
    ed.analyse_disk_all(path,path_out=None,force=force_disk,tag='',order='<',ps=True,pdf=True,do_plot=True)

    path='/gpfs/data1/phennebe/DC_1_consmom_hr_restart/'
#    ed.analyse_disk_all(path,path_out=None,force=force_disk,tag='',order='<',ps=True,pdf=True,do_plot=True)
    ed.analyse_disk_all(path,path_out=None,force=force_disk,tag='',order='<',ps=True,pdf=True,do_plot=True)


    path='/drf/projets/capucine/phennebe/B335_vweakB_FLD_sink/'
    ed.analyse_disk_all(path,path_out=None,force=force_disk,tag='',order='<',ps=True,pdf=True,do_plot=True)


    path='/drf/projets/capucine/phennebe/B335_weakB_FLD_sink/'
    ed.analyse_disk_all(path,path_out=None,force=force_disk,tag='',order='<',ps=True,pdf=True,do_plot=True)



    path='/drf/projets/capucine/phennebe/B335_FLD_sink/'
    ed.analyse_disk_all(path,path_out=None,force=force_disk,tag='',order='<',ps=True,pdf=True,do_plot=True)

    path='/drf/projets/capucine/phennebe/B335_highrot_FLD_sink/'
    ed.analyse_disk_all(path,path_out=None,force=force_disk,tag='',order='<',ps=True,pdf=True,do_plot=True)




#STOP

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.,0.5/128.]
path='/gpfs/data1/phennebe/DC_7_ps_baro2/'
num_v=[10,20,30,40,50,60,70,80]
#ed.analyse_disk_serie(path,num_v,path_out=None,force=False,ps=False,tag='',order='<')


#for num in num_v:

#    me.make_image_zoom(path,num,zoom_v,sinks=False,center_dmax=True,ps=True,force=True,tag='zoom')


path='/gpfs/data1/phennebe/DC_BD/'
zoom_v=[0.5/32.]
#me.make_image_zoom_all(path,zoom_v,force=True,center_dmax=True,col_dens_only=True,tag='movie',mag_im=False)


path='/gpfs/data1/phennebe/DC_1_consmom_hr/'
#zoom_v=[0.5/128.]
zoom_v=[0.5/64.]
#me.make_image_zoom_all(path,zoom_v,force=True,center_dmax=True,col_dens_only=True,tag='movie',mag_im=False)
#me.make_image_zoom_all(path,zoom_v,force=True,center_zoom=True,col_dens_only=True,tag='movie')

path='/gpfs/data1/phennebe/DC_1_consmom_hr/'
zoom_v=[0.5/1024.]
#me.make_image_zoom_all(path,zoom_v,force=True,center_zoom=True,col_dens_only=True,tag='movie')




zoom_v=[0.5/256.,0.5/1024.]

num_v = [8,10,12,14,16,18,20]
num_v = [25,26,27,28,29,30,34]

#num_v = [57]





list_path=[]

path='/gpfs/data1/phennebe/DC_6_ps_baro2/'
list_path.append((path))

path='/gpfs/data1/phennebe/DC_7_ps_baro2/'
list_path.append((path))


path='/gpfs/data1/phennebe/DC_8_ps_baro2/'
list_path.append((path))

path='/gpfs/data1/phennebe/DC_9_ps_baro2/'
list_path.append((path))

path='/gpfs/data1/phennebe/DC_3_ps_baro2/'
list_path.append((path))
#mass_v, time_v = me.draw_mass_sink_cvs(path)


path='/gpfs/data1/phennebe/DC_4_ps_baro2/'
list_path.append((path))
#mass_v, time_v = me.draw_mass_sink_cvs(path)


path='/gpfs/data1/phennebe/DC_5_ps_baro2/'
list_path.append((path))
#mass_v, time_v = me.draw_mass_sink_cvs(path)


path='/gpfs/data1/phennebe/DC_1_consmom/'
list_path.append((path))
mass_v, time_v = me.draw_mass_sink_cvs(path)


path='/gpfs/data1/phennebe/DC_1_consmom_rest/'
list_path.append((path))
mass_v, time_v = me.draw_mass_sink_cvs(path)


path='/gpfs/data1/phennebe/DC_1_consmom_hr/'
list_path.append((path))
#mass_v, time_v = me.draw_mass_sink_cvs(path)


path='/gpfs/data1/phennebe/DC_1_consmom_hr_restart/'
list_path.append((path))


path='/gpfs/data1/phennebe/DC_1_FLD/'
list_path.append((path))


path='/gpfs/data1/phennebe/DC_3_ps/'
list_path.append((path))
#mass_v, time_v = me.draw_mass_sink_cvs(path)


path='/gpfs/data1/phennebe/DC_0b/'
list_path.append((path))
#mass_v, time_v = me.draw_mass_sink_cvs(path)


path='/gpfs/data1/phennebe/DC_1_FLD_hr/'
list_path.append((path))


path='/gpfs/data1/phennebe/B335_FLD_sink/'
list_path.append((path))
#mass_v, time_v = me.draw_mass_sink_cvs(path)

path='/gpfs/data1/phennebe/B335_highrot_FLD_sink/'
list_path.append((path))
#mass_v, time_v = me.draw_mass_sink_cvs(path)


path='/gpfs/data1/phennebe/B335_weakB_FLD_sink/'
list_path.append((path))
#mass_v, time_v = me.draw_mass_sink_cvs(path)


path='/gpfs/data1/phennebe/B335_vweakB_FLD_sink/'
list_path.append((path))
#mass_v, time_v = me.draw_mass_sink_cvs(path)






zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.,0.5/128.]



if( not skip):
    for path in list_path:

        print path

        me.make_image_zoom_last(path,zoom_v,sinks=False,force=False,center_dmax=True,ps=True,prob_df=True)



STOP





num_v = [54]
path='/gpfs/data1/phennebe/B335_noturb_norot_hydro_hr4/'

#for num in num_v:

#    me.make_image_zoom(path,num,zoom_v,sinks=False,center_dmax=True,ps=True,force=False,tag='zoom')

#num_v = [6]
num_v = [48]

path='/gpfs/data1/phennebe/B335_noturb_norot_hydro_pert_llf/'

#for num in num_v:

#    me.make_image_zoom(path,num,zoom_v,sinks=False,center_dmax=True,ps=True,force=False,tag='zoom')




#path='/gpfs/data1/phennebe/DC_1_consmom_hr/'
#zoom_v=[0.5/16.,0.5/64.]

zoom_v=[0.5/16.,0.5/64.,0.5/256.,0.5/1024.]

num_v = [20,18,16,14,13,12,11,10]

path='/gpfs/data1/phennebe/DC_hydro_larson_norot/'

#for num in num_v:

#    me.make_image_zoom(path,num,zoom_v,sinks=False,center_dmax=True,ps=True,force=False)




list_path=[]


path='/gpfs/data1/phennebe/DC_hydro_larson2/'
list_path.append((path))

path='/gpfs/data1/phennebe/DC_hydro_larson/'
list_path.append((path))


path='/gpfs/data1/phennebe/DC_hydro_larson_norot/'
list_path.append((path))


path='/gpfs/data1/phennebe/DC_hydro_larson_noeos/'
list_path.append((path))



zoom_v=[0.5/16.,0.5/64.,0.5/256.,0.5/1024.]




if( not skip):
    for path in list_path:

        print path

        me.make_image_zoom_last(path,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)




list_path=[]


path='/gpfs/data1/phennebe/DC_2/'
list_path.append((path))

path='/gpfs/data1/phennebe/DC_4/'
list_path.append((path))


path='/gpfs/data1/phennebe/B335_FLD/'
list_path.append((path))




path='/gpfs/data1/phennebe/B335_noturb_norot_hydro_isosink/'
list_path.append((path))



path='/gpfs/data1/phennebe/B335_noturb_norot_hydro_hr2/'
list_path.append((path))


path='/gpfs/data1/phennebe/B335_noturb_norot_hydro_pert_llf/'
list_path.append((path))


path='/gpfs/data1/phennebe/B335_weakB_FLD_sink_feed/'
list_path.append((path))
#mass_v, time_v = me.draw_mass_sink_cvs(path)




path='/gpfs/data1/phennebe/B335_noturb_norot_hydro_hr/'
list_path.append((path))

path='/gpfs/data1/phennebe/B335_noturb_norot_hydro/'
list_path.append((path))


path='/gpfs/data1/phennebe/B335_turb_norot_hydro/'
list_path.append((path))


path='/gpfs/data1/phennebe/B335_turb_rot_hydro/'
list_path.append((path))


path='/gpfs/data1/phennebe/B335_simple_hydro/'
list_path.append((path))





#path='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint_frag2/'
#list_path.append((path))

path='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint_frag/'
list_path.append((path))

#path='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint_frag2/'
#list_path.append((path))

path='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint_sink/'
list_path.append((path))

path='/gpfs/data1/phennebe/RB_DENSE_CORE4_hydro_sink/'
list_path.append((path))

path='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint_sink_bis/'
list_path.append((path))


path='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint_sink_FLD/'
list_path.append((path))


path='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint_sink_FLD_bis/'
list_path.append((path))


path='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint_sink_FLD_bis_hr/'
list_path.append((path))


path='/gpfs/data1/phennebe/B335_simple_lowrot/'
list_path.append((path))


path='/gpfs/data1/phennebe/B335_simple_mhd_ieal/'
list_path.append((path))


path='/gpfs/data1/phennebe/B335_simple_mhd_ieal_Bweak/'
list_path.append((path))


path='/gpfs/data1/phennebe/B335_simple/'
list_path.append((path))


path='/gpfs/data1/phennebe/B335_simple_mu4/'
list_path.append((path))


path='/gpfs/data1/phennebe/B335_simple_turb0.2/'
list_path.append((path))

path='/gpfs/data1/phennebe/B335_simple_turb0.5/'
list_path.append((path))


path='/gpfs/data1/phennebe/B335_simple_Bweak/'
list_path.append((path))


path='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint/'
list_path.append((path))


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]



if( not skip):
    for path in list_path:

        me.make_image_zoom_last(path,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)




Force=True

#Force=True




path='/gpfs/data1/phennebe/B335_noturb_norot_hydro/'
zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

#zoom_v=[0.5/4.,0.5/16.,0.5/64.]

#num_v = [10,16]

num_v = [20]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=Force,center_dmax=True,ps=True,savetxt=True,mag_abs=False,mag_im=False)

    i_im=2
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)
    i_im=3
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)
    i_im=4
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)



path='/gpfs/data1/phennebe/B335_turb_norot_hydro/'
zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

#zoom_v=[0.5/4.,0.5/16.,0.5/64.]

#num_v = [8]
num_v = [15]

#num_v = [64]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=Force,center_dmax=True,ps=True,savetxt=True,mag_abs=False,mag_im=False)

    i_im=2
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)
    i_im=3
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)
    i_im=4
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)



path='/gpfs/data1/phennebe/B335_turb_rot_hydro/'
zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

#zoom_v=[0.5/4.,0.5/16.,0.5/64.]

#num_v = [6]
num_v = [14]

#num_v = [64]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=Force,center_dmax=True,ps=True,savetxt=True,mag_abs=False,mag_im=False)

    i_im=2
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)
    i_im=3
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)
    i_im=4
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)



path='/gpfs/data1/phennebe/B335_simple_lowrot/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

#zoom_v=[0.5/4.,0.5/16.,0.5/64.]

num_v = [64,120,300]

#num_v = [64]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=Force,center_dmax=True,ps=True,savetxt=True,mag_abs=True)

    i_im=2
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)
    i_im=3
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)
    i_im=4
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)




path='/gpfs/data1/phennebe/B335_simple/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [64,150,300]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=Force,center_dmax=True,ps=True,savetxt=True,mag_abs=True)

    i_im=2
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)
    i_im=3
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)
    i_im=4
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)



path='/gpfs/data1/phennebe/B335_simple_turb0.5/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [20,80,150]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=Force,center_dmax=True,ps=True,savetxt=True,mag_abs=True)

    i_im=2
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)
    i_im=3
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)
    i_im=4
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)



path='/gpfs/data1/phennebe/B335_simple_turb0.2/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [20,80,160]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=Force,center_dmax=True,ps=True,savetxt=True,mag_abs=True)

    i_im=2
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)
    i_im=3
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)
    i_im=4
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)



path='/gpfs/data1/phennebe/B335/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [50,100,200,800]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=Force,center_dmax=True,ps=True,savetxt=True,mag_abs=True)

    i_im=2
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)
    i_im=3
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)
    i_im=4
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)



path='/gpfs/data1/phennebe/B335_simple_Bweak/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [64,128,248]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=Force,center_dmax=True,ps=True,savetxt=True,mag_abs=True)

    i_im=2
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)
    i_im=3
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)
    i_im=4
    mpc.do_cut(num,path,i_im=i_im,ps=True,pdf=False)





STOP


path='/gpfs/data1/phennebe/B335_simple_mhd_ieal_Bweak/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [14]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=True,center_dmax=True,ps=True,savetxt=True)






path='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint_frag/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [20]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)





path='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint_part2/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [10,20,30,40,50,60,70,80,100,117,128,140,150,160,166,176,188,210,248]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)





path='/gpfs/data1/phennebe/RB_DENSE_CORE4_hydro_lowturb_part2/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [8,10,12,14,20,24,30,34,40,42,46,50,54,56]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)




path='/gpfs/data1/phennebe/RB_DENSE_CORE4_hydro_lowrot_turb_part2/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [8,10,12,14,16,18,20,24,28,32,36,40,42,50,60,70,90,100,110]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)




path='/gpfs/data1/phennebe/RB_DENSE_CORE4_hydro_lowrot_part2/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [8,10,12,14,16,18,20,24,28,30,34,40,50,60,70,78,82,86,90,92]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)




path='/gpfs/data1/phennebe/RB_DENSE_CORE4_hydro_part2/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [8,10,12,14,16,18,20,24,28,30,34,40,42,50,60,70,80,90,110,110]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)




path='/gpfs/data1/phennebe/RB_DENSE_CORE4_imhd_part2/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [8,10,16,20,30,40,50,60,70,80,90,100,110,118,140,160,180,200]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)









path='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [1600,1500,1400,1300,1200,1120,1060,1000,960,900,860,808]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)





path='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint_sink/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [294,295,296,297,298,299,300,301,310,320,330,340,286,287,288,289,290,291,292,293]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=True,force=False,center_dmax=True,ps=True)






path='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint_part/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [60,80,100,120,140,150,180,210,230,250,270,290,400,500,600,700]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)





path='/gpfs/data1/phennebe/RB_DENSE_CORE4_pint2_part/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [60,80,100,150,170,190,210,230,250,300,400,500,600,700]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)





path='/gpfs/data1/phennebe/RB_DENSE_CORE4_hydro_part/'
#path_out='/gpfs/data1/phennebe/RB_DENSE_CORE4_imhd_sink/'
#if ( len(glob.glob(path_out)) == 0 ):
#    str_mkdir = 'mkdir ' + path_out
#    os.system(str_mkdir)
#path_out=path_out+'/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]


num_v = [10,20,30,40,60,80,100,120,140,160]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True,conv_sink=2)





path='/gpfs/data1/phennebe/RB_DENSE_CORE4_imhd_part/'
#path_out='/gpfs/data1/phennebe/RB_DENSE_CORE4_imhd_sink/'
#if ( len(glob.glob(path_out)) == 0 ):
#    str_mkdir = 'mkdir ' + path_out
#    os.system(str_mkdir)
#path_out=path_out+'/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]


num_v = [10,20,30,40,50,100,150,200,250,300]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True,conv_sink=2)





path='/gpfs/data1/phennebe/RB_DENSE_CORE4_hydro_sink/'
#path_out='/gpfs/data1/phennebe/RB_DENSE_CORE4_imhd_sink/'
#if ( len(glob.glob(path_out)) == 0 ):
#    str_mkdir = 'mkdir ' + path_out
#    os.system(str_mkdir)
#path_out=path_out+'/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]


num_v = [40,50,60,70,80,115,150,200,250,300,350,400]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=True,force=False,center_dmax=True,ps=True,conv_sink=2)





path='/gpfs/data1/phennebe/RB_DENSE_CORE4_imhd_sink/'
#path_out='/gpfs/data1/phennebe/RB_DENSE_CORE4_imhd_sink/'
#if ( len(glob.glob(path_out)) == 0 ):
#    str_mkdir = 'mkdir ' + path_out
#    os.system(str_mkdir)
#path_out=path_out+'/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]


num_v = [40,50,60,70,80,115,150,200,250,300,350,400]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=True,force=False,center_dmax=True,ps=True,conv_sink=2)







#me.make_image_zoom_all(path_in,zoom_v,path_out=path_out)




path='/gpfs/data1/phennebe/B335/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [820,800,750,700,650,620,600,580,530,500,450,400,350,300]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)




path='/gpfs/data1/phennebe/RB_DENSE_CORE4_hex/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [386]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)









path='/gpfs/data1/phennebe/RB_DENSE_CORE4_qua/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [984]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=True,center_dmax=True,ps=True)





path='/gpfs/data1/phennebe/RB_DENSE_CORE4_imhd_part/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [300,400,500,600]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)





path='/gpfs/data1/phennebe/RB_DENSE_CORE4_hydro_part/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

num_v = [2,4,8,10,20,30]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,center_dmax=True,ps=True)





path='/gpfs/data1/phennebe/RB_DENSE_CORE4_ter/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

#zoom_v=[0.5/16.,0.5/64.]

num_v = [1000]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=True,center_dmax=True,ps=True)






path='/gpfs/data1/bcommerc/patrick/'

path_out='/gpfs/data1/phennebe/COEUR_BENOIT/'

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

#num_v = [50,100,150,198]

#num_v = [17,19,24,38,44,55,61,70,98]

#num_v = [180,150,100,50]
num_v = [4001,3001,2001]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=True,center_dmax=True,path_out=path_out,ps=True)





path='/gpfs/data1/phennebe/RB_DENSE_CORE4_bis/'


zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

#num_v = [50,100,150,198]

#num_v = [17,19,24,38,44,55,61,70,98]

#num_v = [180,150,100,50]
num_v = [50,40,30]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=True,center_dmax=True)



path='/gpfs/data1/phennebe/RB_DENSE_CORE4/'
#path='/gpfs/data1/phennebe/RB_DENSE_CORE3/'
#path_out='/gpfs/data1/phennebe/RB_DENSE_CORE/'


#mass_v, time_v = me.draw_mass_sink_cvs(path)

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

#num_v = [50,100,150,198]

#num_v = [17,19,24,38,44,55,61,70,98]

#num_v = [180,150,100,50]
num_v = [70,60,50,40,30]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=True,center_dmax=True)


STOP



path='/gpfs/data1/phennebe/RB_DENSE_CORE_hydro/'
#path_out='/gpfs/data1/phennebe/RB_DENSE_CORE/'


#mass_v, time_v = me.draw_mass_sink_cvs(path)

zoom_v=[0.5,0.5/4.,0.5/16.,0.5/64.]

#num_v = [50,100,150,198]

#num_v = [17,19,24,38,44,55,61,70,98]

num_v = [20,40,60,100,180]

for num in num_v:

    me.make_image_zoom(path,num,zoom_v,sinks=False,force=False,mag_im=False)


