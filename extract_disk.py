#module to extract structure and visualise simulations
import sys
import numpy as np
import os 

import pymses
from pymses.filters import CellsToPoints
from pymses.sources.ramses.tree_utils import octree_compute_neighbors

from pymses.sources.ramses import output

from pymses.analysis import Camera, raytracing, slicing
from pymses.analysis import ScalarOperator, FractionOperator, MaxLevelOperator
from pymses.analysis.point_sampling import PointSamplingProcessor
from pymses.utils.regions import Sphere
from pymses.utils.regions import Box
from pymses.filters import RegionFilter
from pymses.filters import PointFunctionFilter


#from matplotlib import pyplot as P
import pylab as P
import glob as glob
import pickle as pickle
from struct import pack


import module_extract as me



#tick_size = 20
#fontlabel_size = 15
#axislabel_size = 20


#params = {'backend': 'wxAgg', 'lines.markersize' : 3, 'axes.labelsize': axislabel_size, 'font.size': fontlabel_size, 'legend.fontsize': fontlabel_size, 'xtick.labelsize': tick_size, 'ytick.labelsize': tick_size, 'axes.linewidth' : 3}
#P.rcParams.update(params)



#fig_width=8                                                                                                \
#fig_height=4                                                                                               \
#fig_size = [fig_width,fig_height]
P.figure(figsize=(8,6))


P.gcf().subplots_adjust(bottom=0.2)
P.gcf().subplots_adjust(left=0.2)


##########################################################
##do the density PDF
##########################################################
def get_rhopdf(path_in,num,path_out=None,force=False,ps=False,xmin=0,xmax=1.,ymin=0,ymax=1.,zmin=0.,zmax=1.,tag='',order='<',do_plot=False):

    if( path_out is not None):
        directory_out = path_out
    else:
        directory_out = path_in

    if(tag != ''):
        tag = tag+'_'

    name_save_pdf=directory_out+'pdf_' + tag + str(num).zfill(5) + '.save'


    #a list of dictionnary to describe the plots
    list_plot = []

    ro=pymses.RamsesOutput(path_in,num,order=order)
        
    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)  

    amr = ro.amr_source(["rho"])

    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    dx = cells.get_sizes()


    mask1 = cells.points[:,0] > xmin
    mask2 = cells.points[:,0] < xmax
    mask3 = cells.points[:,1] > ymin
    mask4 = cells.points[:,1] < ymax
    mask5 = cells.points[:,2] > zmin
    mask6 = cells.points[:,2] < zmax
    mask = mask1*mask2*mask3*mask4*mask5*mask6


    AU,pc,Ms,Myr,scale_n,scale_d,scale_t,scale_l,scale_v,scale_T2,scale_ener,scale_mag,microG,km_s,Cwnm,scale_mass,unit_col,lbox_pc = me.normalisation_averliat(ro)

    vol  = dx**3
    mass = cells["rho"]*vol * scale_mass / Ms

    log_rho = np.log10(cells['rho'])
    log_min=min(log_rho)
    log_max=max(log_rho)
    nbin=50
    width=(log_max-log_min)/nbin

    hist_logrho , hist_rho_edges = P.histogram(log_rho[mask],bins=nbin,range=(log_min,log_max),weights=mass[mask])
#    pdf_logrho , histrhomv_edges = P.histogram(log_rho[mask],bins=nbin,range=(log_min,log_max),weights=vol[mask])



    if(do_plot):
        #mass weighted histogram of density 
        P.clf()
        P.bar(hist_rho_edges[:-1],np.log10(hist_logrho)+3.,width,-3.,color='none',edgecolor='k')

        xlabel='$log(n)$ (cm$^{-3}$)'
        ylabel='$log(M)$ (M$_\odot$)'

        P.xlabel(r'$log(n)$ (cm$^{-3}$)')
        P.ylabel(r'$log(M)$ (M$_\odot$)')

        name_im = directory_out+'pdfmw_logrho_bd_'+tag+str(num).zfill(5)

        P.savefig(name_im+'.ps')
        P.savefig(name_im+'.png')




    return hist_logrho, hist_rho_edges

##########################################################
#calculate disk radius by looking at the miminim of the mass weighted PDF 
#then ask that Vrad / Vpar < 0.5 (i.e. rotation dominates infall)
#here assume a single disk - results are not correct if more than one disk
#it also assumes that the density maximum belongs to the disk and is close to its center
##########################################################
def pdf_to_singledisc_rad(path,num,path_out=None,force=False,ps=False,xmin=0,xmax=1.,ymin=0,ymax=1.,zmin=0.,zmax=1.,tag='',order='<',do_plot=False):

    print(num)


    #get the mass weighted density prob density function
    hist_logrho,hist_rho_edges = get_rhopdf(path,num,path_out=path_out,force=force,ps=ps,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,zmin=zmin,zmax=zmax,tag=tag,order=order,do_plot=do_plot)


    nbin = len(hist_logrho)

    mask_min = hist_rho_edges[0:nbin] > 8. 
    mask_max = hist_rho_edges[0:nbin] < 11.

    mask = mask_min*mask_max


    amin = np.argmin(hist_logrho[mask])

    #density at the disk edge
    log_rho_disk = hist_rho_edges[0:nbin][mask][amin]

    print 'log_rho_disk', log_rho_disk





#    pdb.set_trace() # debug mode     

    ro=pymses.RamsesOutput(path,num)
    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)

    AU,pc,Ms,Myr,scale_n,scale_d,scale_t,scale_l,scale_v,scale_T2,scale_ener,scale_mag,microG,km_s,Cwnm,scale_mass,unit_col,lbox_pc = me.normalisation_averliat(ro)

    time = ro.info['time'] * scale_t / Myr

    amr = ro.amr_source(["rho","vel","Br"])

    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    dx = cells.get_sizes()

    pos = cells.points


    #cells belonging to the disk
    mask_rho_disk = np.log10(cells['rho']) > log_rho_disk
    
    print 'disk cell number',np.sum(mask_rho_disk)

    #mean magnetic field in the whole selected region
    mrd=mask_rho_disk

    #volume of the selected region 
    vol = np.sum(dx[mrd]**3)

    mag_mean_broad=np.sum(np.sqrt((cells['Br'][:,0][mrd]**2+cells['Br'][:,1][mrd]**2+cells['Br'][:,2][mrd]**2))*dx[mrd]**3)/vol

    #normalise it in cgs
    mag_mean_broad = mag_mean_broad*scale_mag


    #maximum density cell
    mask_rho_max = cells['rho'] == np.max(cells['rho'])

    #mass disk
    mass_disk = np.sum(dx[mask_rho_disk]**3 * cells['rho'][mask_rho_disk])
    mass_disk = mass_disk * scale_mass / Ms

    #position of maximum density cell
    pos_max_x = pos[:,0][mask_rho_max][0]
    pos_max_y = pos[:,1][mask_rho_max][0]
    pos_max_z = pos[:,2][mask_rho_max][0]

    #position of disk cell
    pos_disk_x = pos[:,0][mask_rho_disk] - pos_max_x
    pos_disk_y = pos[:,1][mask_rho_disk] - pos_max_y
    pos_disk_z = pos[:,2][mask_rho_disk] - pos_max_z

    #position of maximum density cell
    vel_max_x = cells['vel'][:,0][mask_rho_max][0]
    vel_max_y = cells['vel'][:,1][mask_rho_max][0]
    vel_max_z = cells['vel'][:,2][mask_rho_max][0]

    #position of disk cell
    vel_disk_x = cells['vel'][:,0][mask_rho_disk] - vel_max_x
    vel_disk_y = cells['vel'][:,1][mask_rho_disk] - vel_max_y
    vel_disk_z = cells['vel'][:,2][mask_rho_disk] - vel_max_z

    #radial component of V
    norm_pos = np.sqrt(pos_disk_x**2+pos_disk_y**2+pos_disk_z**2)
    mask = norm_pos == 0.
    norm_pos[mask] = 1.e-10
    Vrad = (vel_disk_x*pos_disk_x + vel_disk_y*pos_disk_y + vel_disk_z*pos_disk_z) / norm_pos

    #non radial component of V
    Vpar_x = (vel_disk_z*pos_disk_y - vel_disk_y*pos_disk_z) / norm_pos
    Vpar_y = (vel_disk_x*pos_disk_z - vel_disk_z*pos_disk_x) / norm_pos
    Vpar_z = (vel_disk_y*pos_disk_x - vel_disk_x*pos_disk_y) / norm_pos
    Vpar = np.sqrt(Vpar_x**2+Vpar_y**2+Vpar_z**2)
    Vpar[mask] = 1.e-10


    #now select the point that have a parallel velocity larger than the radial one 
    mask = (-Vrad) / Vpar < 0.5

    

    if(np.sum(mask) <= 1):
        return time, 0., 0., 0., 0. , log_rho_disk, mag_mean_broad, 0., 0. #, max_rad_disk, mean_rad_disk, mean_out_rad_disk


    pos_disk_x = pos_disk_x[mask]
    pos_disk_y = pos_disk_y[mask]
    pos_disk_z = pos_disk_z[mask]

    vel_disk_x = vel_disk_x[mask]
    vel_disk_y = vel_disk_y[mask]
    vel_disk_z = vel_disk_z[mask]

#    print 'sum mask',np.sum(mask),len(pos_disk_x)

    #determine the direction of the mean angular momentum
    Jx = np.mean(pos_disk_y*vel_disk_z-pos_disk_z*vel_disk_y)
    Jy = np.mean(pos_disk_z*vel_disk_x-pos_disk_x*vel_disk_z)
    Jz = np.mean(pos_disk_x*vel_disk_y-pos_disk_y*vel_disk_x)

    Jnorm = np.sqrt(Jx**2+Jy**2+Jz**2)
    Jx=Jx/Jnorm
    Jy=Jy/Jnorm
    Jz=Jz/Jnorm

    #project the radial vector in the disk plane
    Jpos = Jx*pos_disk_x + Jy*pos_disk_y + Jz*pos_disk_z
    pos_dpl_x = pos_disk_x - Jpos*Jx
    pos_dpl_y = pos_disk_y - Jpos*Jy
    pos_dpl_z = pos_disk_z - Jpos*Jz

    if( pos_dpl_x[0]**2+pos_dpl_x[1]**2+pos_dpl_x[2]**2 != 0.):
        pos_x = pos_dpl_x[0]
        pos_y = pos_dpl_y[0]
        pos_z = pos_dpl_z[0]
    else:
        pos_x = pos_dpl_x[1]
        pos_y = pos_dpl_y[1]
        pos_z = pos_dpl_z[1]

    #calculate the cosinus 
    cos_pos = (pos_dpl_x*pos_x+pos_dpl_y*pos_y+pos_dpl_z*pos_z)

    mask = pos_dpl_x*pos_dpl_x+pos_dpl_y*pos_dpl_y+pos_dpl_z*pos_dpl_z > 0.

    cos_pos = cos_pos  / np.sqrt(pos_x*pos_x+pos_y*pos_y+pos_z*pos_z)
    cos_pos[mask] = cos_pos[mask]  / np.sqrt(pos_dpl_x*pos_dpl_x+pos_dpl_y*pos_dpl_y+pos_dpl_z*pos_dpl_z)[mask]


    #radius of disk cells
#    rad_disk = np.sqrt( (pos_disk[:,0]-pos_max_x)**2 + (pos_disk[:,1]-pos_max_y)**2 + (pos_disk[:,2]-pos_max_z)**2 )
    rad_disk = np.sqrt( (pos_disk_x)**2 + (pos_disk_y)**2 + (pos_disk_z)**2 )



    #mean disk radius
    mean_rad_disk = np.mean(rad_disk)
    median_rad_disk = np.median(rad_disk)

    
    arg_d = np.argsort(rad_disk)


    ndisk = len(arg_d)

    print 'ndisk',ndisk

    mask_out_disk = arg_d[np.int(ndisk*0.9):ndisk-1]

    #mean radius of the 90% point having higher radius
    mean_out_rad_disk = np.mean(rad_disk[mask_out_disk])


    #now look at the maximum in a series of direction
    ndir = 500
    rad_loc_v = np.zeros(ndir)
    cos_min=-1.
    cos_max = cos_min + 2./(ndir+1)
    for ld in range(ndir):
        
        mask1 = cos_pos >= cos_min
        mask2 = cos_pos <  cos_max
        mask = mask1*mask2


        print 'len mask, sum mask',len(mask),np.sum(mask)

        if(np.sum(mask) !=0):
#            pdb.set_trace()
            rad_loc_v[ld] = np.max(rad_disk[mask])
        else:
            rad_loc_v[ld] = 0.

        cos_min = cos_min + 2./(ndir+1)
        cos_max = cos_max + 2./(ndir+1)


    max_rad_loc = np.max(rad_loc_v[rad_loc_v<500/lbox/pc*AU]) * lbox * pc / AU

    min_rad_loc = np.min(rad_loc_v[rad_loc_v<500/lbox/pc*AU]) * lbox * pc / AU

    mean_rad_loc = np.mean(rad_loc_v[rad_loc_v<500/lbox/pc*AU]) * lbox * pc / AU

    median_rad_loc = np.median(rad_loc_v[rad_loc_v<500/lbox/pc*AU]) * lbox * pc / AU


    #--------------------------
    #Ajout par AV pour essayer de mieux estimer la taille des disques
    rad_loc_v = rad_loc_v*lbox*pc/AU  #listes des rayons en AU
    
    taille_bin = 4 #AU
    largeur_modale = 2

    min_serie = np.floor(np.min(rad_loc_v))
    min_serie -= min_serie%taille_bin
    max_serie = np.floor(np.max(rad_loc_v))+1
    max_serie += (taille_bin-max_serie%taille_bin)

    nbr_bin=(max_serie-min_serie)/taille_bin
    bins=np.arange(nbr_bin)*taille_bin+min_serie+taille_bin/2.

    count=np.zeros(len(bins))
    
    for i in range(len(rad_loc_v)):
        ind=np.int(np.floor((rad_loc_v[i]-min_serie)/taille_bin))
        count[ind] += 1

    bin_central = bins[np.argmax(count)]
    lim_inf = bin_central -taille_bin/2. -largeur_modale*taille_bin
    lim_sup = bin_central +taille_bin/2. +largeur_modale*taille_bin

    rad_new=rad_loc_v[rad_loc_v>lim_inf]
    rad_new=rad_new[rad_new<lim_sup]

    rad_estim = np.median(rad_new)
    #--------------------------



    mean_rad_disk =  mean_rad_disk * lbox * pc / AU

    mean_out_rad_disk = mean_out_rad_disk * lbox * pc / AU

#    pdb.set_trace() # debug mode     

    return time, mass_disk, max_rad_loc, min_rad_loc, mean_rad_loc, log_rho_disk, mag_mean_broad, median_rad_loc, rad_estim #, max_rad_disk, mean_rad_disk, mean_out_rad_disk
##########################################################
##########################################################
##do series of analysis
##########################################################
##########################################################
def analyse_disk_serie(path_in,num_v,path_out=None,force=False,ps=False,tag='',order='<',do_plot=False):


    if( path_out is not None):
        directory_out = path_out
    else:
        directory_out = path_in





    #store the results
    name_save=directory_out+'/res_disk.save'
    if ( len(glob.glob(name_save)) != 0 and not force ):

        f = open(name_save,'r')
        time_v = pickle.load(f)
        mdisk_v = pickle.load(f)
        rmax_v = pickle.load(f)
        rmin_v = pickle.load(f)
        rmean_v = pickle.load(f)
        num_v = pickle.load(f)
        lrhod_v = pickle.load(f)
        Bmean_v = pickle.load(f)
        f.close()

        return time_v, mdisk_v,rmax_v,rmin_v,rmean_v,lrhod_v


    time_v = np.zeros(len(num_v))
    mdisk_v = np.zeros(len(num_v))
    rmax_v = np.zeros(len(num_v))
    rmin_v = np.zeros(len(num_v))
    rmean_v = np.zeros(len(num_v))
    lrhod_v = np.zeros(len(num_v))
    Bmean_v = np.zeros(len(num_v))


    ii = 0
    for num in num_v:


        if(num >= 1900):
            continue

        time,mass_disk, max_rad_loc, min_rad_loc, mean_rad_loc , log_rho_disk , mag_mean_broad = pdf_to_singledisc_rad(path_in,num,xmin=0,xmax=1.,ymin=0,ymax=1.,zmin=0.,zmax=1.,tag='',order='<',do_plot=do_plot)

        time_v[ii] = time
        mdisk_v[ii] = mass_disk
        rmax_v[ii]  = max_rad_loc
        rmin_v[ii]  = min_rad_loc
        rmean_v[ii] = mean_rad_loc
        lrhod_v[ii] = log_rho_disk
        Bmean_v[ii] = mag_mean_broad

        ii = ii + 1

        print 'mass disk', mass_disk

        print 'max_rad_loc, min_rad_loc, mean_rad_loc',max_rad_loc, min_rad_loc, mean_rad_loc

#        print 'max_rad, mean_rad, mean_out_rad', max_rad, mean_rad, mean_out_rad



    #store the results
    name_save=directory_out+'/res_disk.save'
    f = open(name_save,'w')
    pickle.dump(time_v,f)
    pickle.dump(mdisk_v,f)
    pickle.dump(rmax_v,f)
    pickle.dump(rmin_v,f)
    pickle.dump(rmean_v,f)
    pickle.dump(num_v,f)
    pickle.dump(lrhod_v,f)
    pickle.dump(Bmean_v,f)
    f.close()

#    pdb.set_trace() # debug mode     

    return time_v, mdisk_v,rmax_v,rmin_v,rmean_v,lrhod_v,Bmean_v
##############################################################################
##############################################################################
### do images for all output ending with a "0"in path
### path_in : the path
### num the output number
##############################################################################
##############################################################################
def analyse_disk_all(path_in,path_out=None,sinks=False,force=True,no_cvs=None,tag='',gcomp=True,conv_sink=1,center_dmax=False,order='<',ps=False,pdf=False,descrip=None,do_plot=False):

    num_v = list()

    def getkey2(item):
        siz=len(item)
        key=int(item[siz-5:siz])
        return key


    name_search= path_in+'output*0'
    for name  in sorted(glob.glob(name_search),key=getkey2):

        print 'find file',name

        siz=len(name)
        num=int(name[siz-5:siz])

        num_v.append(num)

#    time_v = np.zeros(len(num_v))
#    mdisk_v = np.zeros(len(num_v))
#    rmax_v = np.zeros(len(num_v))
#    rmin_v = np.zeros(len(num_v))
#    rmean_v = np.zeros(len(num_v))


    time_v,mdisk_v,rmax_v,rmin_v,rmean_v,lrhod_v,Bmean_v = analyse_disk_serie(path_in,num_v,path_out=path_out,force=force,ps=ps,tag=tag,order=order,do_plot=do_plot)


#    pdb.set_trace() # debug mode     

    plot_disk_serie(path_in,time_v,mdisk_v,rmax_v,rmin_v,rmean_v,lrhod_v,Bmean_v,path_out=path_out,ps=ps,tag=tag)



    return
##############################################################################
##############################################################################
### do images for all output ending with a "0"in path
### path_in : the path
### num the output number
##############################################################################
##############################################################################
def plot_disk_serie(path_in,time_v,mdisk_v,rmax_v,rmin_v,rmean_v,lrhod_v,Bmean_v,path_out=None,ps=False,pdf=False,tag=''):


#    pdb.set_trace() # debug mode     

    if( path_out is not None):
        directory_out = path_out
    else:
        directory_out = path_in


    P.clf()


#    P.plot([0.,3.],[1.e4,1.e1],color='k')

    P.plot(time_v,np.log10(rmax_v),color='k',linewidth=2)

    P.plot(time_v,np.log10(rmin_v),color='r',linewidth=2)

    P.plot(time_v,np.log10(rmean_v),color='b',linewidth=2)


    P.xlabel(r'$t \, (Myr)$')
    P.ylabel('disk radius (AU)')

#    P.legend( (r'$M^{-1}$','all structures')  ,loc='upper right')
#    P.title=('all structures')

    if(ps):
        P.savefig(directory_out+'rad_disk_t.ps')
    if(pdf):
        P.savefig(directory_out+'rad_disk_t.pdf')
    P.savefig(directory_out+'rad_disk_t.jpeg')



    mass_sink_v, time_sink_v = me.draw_mass_sink_cvs(path_in,file_exclude=None,path_out=directory_out,no_cvs=None)

    P.clf()

    P.plot(time_v,np.log10(mdisk_v),color='k',linewidth=2)

    P.plot(time_sink_v,np.log10(mass_sink_v),color='r',linewidth=2)

    P.xlabel(r'$t \, (Myr)$')
    P.ylabel('mass (Ms)')

#    P.legend( (r'$M^{-1}$','all structures')  ,loc='upper right')
#    P.title=('all structures')

    if(ps):
        P.savefig(directory_out+'mass_disk_t.ps')
    if(pdf):
        P.savefig(directory_out+'mass_disk_t.pdf')
    P.savefig(directory_out+'mass_disk_t.jpeg')




    P.clf()

    P.plot(time_v,lrhod_v,color='k',linewidth=2)

    P.xlabel(r'$t \, (Myr)$')
    P.ylabel(r'$\rho$ (cm$^{-3}$)')

#    P.legend( (r'$M^{-1}$','all structures')  ,loc='upper right')
#    P.title=('all structures')

    if(ps):
        P.savefig(directory_out+'rho_disk_t.ps')
    if(pdf):
        P.savefig(directory_out+'rho_disk_t.pdf')
    P.savefig(directory_out+'rho_disk_t.jpeg')




    P.clf()

    P.plot(time_v,np.log10(Bmean_v),color='k',linewidth=2)

    P.xlabel(r'$t \, (Myr)$')
    P.ylabel(r'$\log(B)$ (G)')

#    P.legend( (r'$M^{-1}$','all structures')  ,loc='upper right')
#    P.title=('all structures')

    if(ps):
        P.savefig(directory_out+'B_disk_t.ps')
    if(pdf):
        P.savefig(directory_out+'B_disk_t.pdf')
    P.savefig(directory_out+'B_disk_t.jpeg')



    return












##########################################################
#Modification of the pdf_to_singledisc_rad function :
#it identifies the disk by looking at the miminim of the mass weighted PDF 
#then ask that Vrad / Vpar < 0.5 (i.e. rotation dominates infall)
#here assume a single disk - results are not correct if more than one disk
#it also assumes that the density maximum belongs to the disk and is close to its center
#it return the cells that belong to the disk
##########################################################
def pdf_to_singledisc_cells(path,num,path_out=None,force=False,ps=False,xmin=0,xmax=1.,ymin=0,ymax=1.,zmin=0.,zmax=1.,tag='',order='<',do_plot=False):

    print(num)


    #get the mass weighted density prob density function
    hist_logrho,hist_rho_edges = get_rhopdf(path,num,path_out=path_out,force=force,ps=ps,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,zmin=zmin,zmax=zmax,tag=tag,order=order,do_plot=do_plot)


    nbin = len(hist_logrho)

    mask_min = hist_rho_edges[0:nbin] > 8. 
    mask_max = hist_rho_edges[0:nbin] < 11.

    mask = mask_min*mask_max


    amin = np.argmin(hist_logrho[mask])

    #density at the disk edge
    log_rho_disk = hist_rho_edges[0:nbin][mask][amin]

    print 'log_rho_disk', log_rho_disk





#    pdb.set_trace() # debug mode     

    ro=pymses.RamsesOutput(path,num)
    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)

    AU,pc,Ms,Myr,scale_n,scale_d,scale_t,scale_l,scale_v,scale_T2,scale_ener,scale_mag,microG,km_s,Cwnm,scale_mass,unit_col,lbox_pc = me.normalisation_averliat(ro)

    #time = ro.info['time'] * scale_t / Myr

    amr = ro.amr_source(["rho","vel"])#,"Br"])

    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    dx = cells.get_sizes()

    pos = cells.points


    #cells belonging to the disk
    mask_rho_disk = np.log10(cells['rho']) > log_rho_disk
    
    print 'disk cell number',np.sum(mask_rho_disk)

    #mean magnetic field in the whole selected region
    #mrd=mask_rho_disk

    #volume of the selected region 
    #vol = np.sum(dx[mrd]**3)

    #mag_mean_broad=np.sum(np.sqrt((cells['Br'][:,0][mrd]**2+cells['Br'][:,1][mrd]**2+cells['Br'][:,2][mrd]**2))*dx[mrd]**3)/vol

    #normalise it in cgs
    #mag_mean_broad = mag_mean_broad*scale_mag


    #maximum density cell
    mask_rho_max = cells['rho'] == np.max(cells['rho'])

    #mass disk
    #mass_disk = np.sum(dx[mask_rho_disk]**3 * cells['rho'][mask_rho_disk])
    #mass_disk = mass_disk * scale_mass / Ms

    #position of maximum density cell
    pos_max_x = pos[:,0][mask_rho_max][0]
    pos_max_y = pos[:,1][mask_rho_max][0]
    pos_max_z = pos[:,2][mask_rho_max][0]

    #position of disk cell
    pos_disk_x = pos[:,0][mask_rho_disk] - pos_max_x
    pos_disk_y = pos[:,1][mask_rho_disk] - pos_max_y
    pos_disk_z = pos[:,2][mask_rho_disk] - pos_max_z

    #velocity of maximum density cell
    vel_max_x = cells['vel'][:,0][mask_rho_max][0]
    vel_max_y = cells['vel'][:,1][mask_rho_max][0]
    vel_max_z = cells['vel'][:,2][mask_rho_max][0]

    #velocity of disk cell
    vel_disk_x = cells['vel'][:,0][mask_rho_disk] - vel_max_x
    vel_disk_y = cells['vel'][:,1][mask_rho_disk] - vel_max_y
    vel_disk_z = cells['vel'][:,2][mask_rho_disk] - vel_max_z

    #radial component of V
    norm_pos = np.sqrt(pos_disk_x**2+pos_disk_y**2+pos_disk_z**2)
    mask = norm_pos == 0.
    norm_pos[mask] = 1.e-10
    Vrad = (vel_disk_x*pos_disk_x + vel_disk_y*pos_disk_y + vel_disk_z*pos_disk_z) / norm_pos

    #non radial component of V
    Vpar_x = (vel_disk_z*pos_disk_y - vel_disk_y*pos_disk_z) / norm_pos
    Vpar_y = (vel_disk_x*pos_disk_z - vel_disk_z*pos_disk_x) / norm_pos
    Vpar_z = (vel_disk_y*pos_disk_x - vel_disk_x*pos_disk_y) / norm_pos
    Vpar = np.sqrt(Vpar_x**2+Vpar_y**2+Vpar_z**2)
    Vpar[mask] = 1.e-10


    #now select the point that have a parallel velocity larger than the radial one 
    mask = (-Vrad) / Vpar < 0.5

    
    
    #if(np.sum(mask) <= 1):
    #    return time, 0., 0., 0., 0. , log_rho_disk, mag_mean_broad #, max_rad_disk, mean_rad_disk, mean_out_rad_disk


    #pos_disk_x = pos_disk_x[mask]
    #pos_disk_y = pos_disk_y[mask]
    #pos_disk_z = pos_disk_z[mask]

    #vel_disk_x = vel_disk_x[mask]
    #vel_disk_y = vel_disk_y[mask]
    #vel_disk_z = vel_disk_z[mask]





    #radius of disk cells
#    rad_disk = np.sqrt( (pos_disk[:,0]-pos_max_x)**2 + (pos_disk[:,1]-pos_max_y)**2 + (pos_disk[:,2]-pos_max_z)**2 )
    #rad_disk = np.sqrt( (pos_disk_x)**2 + (pos_disk_y)**2 + (pos_disk_z)**2 )



    #mean disk radius
    #mean_rad_disk = np.mean(rad_disk)

    return mask_rho_disk, mask





##########################################################
#Modification of the pdf_to_singledisc_cells function :
#it identifies the disk by looking at the miminim of the mass weighted PDF 
#then ask that Vrad / Vpar < 0.5 (i.e. rotation dominates infall)
#here assume a single disk - results are not correct if more than one disk
#it also assumes that the density maximum belongs to the disk and is close to its center
#it return the cells that belong to the disk
##########################################################
def pdf_to_singledisc_for_plot(path,num,path_out=None,force=False,ps=False,xmin=0,xmax=1.,ymin=0,ymax=1.,zmin=0.,zmax=1.,tag='',order='<',do_plot=False):

    print(num)


    #get the mass weighted density prob density function
    hist_logrho,hist_rho_edges = get_rhopdf(path,num,path_out=path_out,force=force,ps=ps,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,zmin=zmin,zmax=zmax,tag=tag,order=order,do_plot=do_plot)


    nbin = len(hist_logrho)

    mask_min = hist_rho_edges[0:nbin] > 8. 
    mask_max = hist_rho_edges[0:nbin] < 11.

    mask = mask_min*mask_max


    amin = np.argmin(hist_logrho[mask])

    #density at the disk edge
    log_rho_disk = hist_rho_edges[0:nbin][mask][amin]

    print 'log_rho_disk', log_rho_disk

    rho_min_disk=10**log_rho_disk





#    pdb.set_trace() # debug mode     

    ro=pymses.RamsesOutput(path,num)
    lbox=ro.info['boxlen'] #boxlen in codeunits (=>pc)

    AU,pc,Ms,Myr,scale_n,scale_d,scale_t,scale_l,scale_v,scale_T2,scale_ener,scale_mag,microG,km_s,Cwnm,scale_mass,unit_col,lbox_pc = me.normalisation_averliat(ro)

    #time = ro.info['time'] * scale_t / Myr

    amr = ro.amr_source(["rho","vel"])#,"Br"])

    cell_source = CellsToPoints(amr)
    cells = cell_source.flatten()
    dx = cells.get_sizes()

    pos = cells.points


    #cells belonging to the disk
    mask_rho_disk = np.log10(cells['rho']) > log_rho_disk
    
    print 'disk cell number',np.sum(mask_rho_disk)

    #mean magnetic field in the whole selected region
    #mrd=mask_rho_disk

    #volume of the selected region 
    #vol = np.sum(dx[mrd]**3)

    #mag_mean_broad=np.sum(np.sqrt((cells['Br'][:,0][mrd]**2+cells['Br'][:,1][mrd]**2+cells['Br'][:,2][mrd]**2))*dx[mrd]**3)/vol

    #normalise it in cgs
    #mag_mean_broad = mag_mean_broad*scale_mag


    #maximum density cell
    mask_rho_max = cells['rho'] == np.max(cells['rho'])

    #mass disk
    #mass_disk = np.sum(dx[mask_rho_disk]**3 * cells['rho'][mask_rho_disk])
    #mass_disk = mass_disk * scale_mass / Ms

    #position of maximum density cell
    pos_max_x = pos[:,0][mask_rho_max][0]
    pos_max_y = pos[:,1][mask_rho_max][0]
    pos_max_z = pos[:,2][mask_rho_max][0]


    #velocity of maximum density cell
    vel_max_x = cells['vel'][:,0][mask_rho_max][0]
    vel_max_y = cells['vel'][:,1][mask_rho_max][0]
    vel_max_z = cells['vel'][:,2][mask_rho_max][0]


    
    
    #if(np.sum(mask) <= 1):
    #    return time, 0., 0., 0., 0. , log_rho_disk, mag_mean_broad #, max_rad_disk, mean_rad_disk, mean_out_rad_disk


    #pos_disk_x = pos_disk_x[mask]
    #pos_disk_y = pos_disk_y[mask]
    #pos_disk_z = pos_disk_z[mask]

    #vel_disk_x = vel_disk_x[mask]
    #vel_disk_y = vel_disk_y[mask]
    #vel_disk_z = vel_disk_z[mask]





    #radius of disk cells
#    rad_disk = np.sqrt( (pos_disk[:,0]-pos_max_x)**2 + (pos_disk[:,1]-pos_max_y)**2 + (pos_disk[:,2]-pos_max_z)**2 )
    #rad_disk = np.sqrt( (pos_disk_x)**2 + (pos_disk_y)**2 + (pos_disk_z)**2 )



    #mean disk radius
    #mean_rad_disk = np.mean(rad_disk)

    return rho_min_disk,pos_max_x,pos_max_y,pos_max_z,vel_max_x,vel_max_y,vel_max_z
