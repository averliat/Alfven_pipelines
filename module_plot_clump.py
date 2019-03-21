import sys
import numpy as np
import pylab as P
import glob as glob
import pdb
import pickle as pickle
from struct import pack

from numpy.linalg import eigh,eig



tick_size = 20
fontlabel_size = 15
axislabel_size = 20


params = {'backend': 'wxAgg', 'lines.markersize' : 3, 'axes.labelsize': axislabel_size, 'font.size': fontlabel_size, 'legend.fontsize': fontlabel_size, 'xtick.labelsize': tick_size, 'ytick.labelsize': tick_size, 'axes.linewidth' : 3}
P.rcParams.update(params)



#fig_width=8                                                                                                \
#fig_height=4                                                                                               \
#fig_size = [fig_width,fig_height]
P.figure(figsize=(8,6))


P.gcf().subplots_adjust(bottom=0.2)
P.gcf().subplots_adjust(left=0.2)

#P.style.use("aanda")


#path_in='/mnt/coast/phennebe/FRIG/FRIG_4/ZOOM4/'
#path_out='/mnt/coast/phennebe/FRIG/FRIG_4/ZOOM4/'

########################################################################################
########################################################################################
## plot the properties of the clumps/cores
########################################################################################
########################################################################################                 
def plot_clump_properties(name_prop_clump,name,path_in,num,path_out=None,ps=False,pdf=False,nbins_glob=25,mark=3,log_min=None,log_max=None,center=None,radius=None,subname=''):

    if( path_out is not None):
        directory_out = path_out
    else:
        directory_out = path_in


    if(subname != ''):
        subname='_'+subname+'_'


    f = open(directory_out+name_prop_clump,'r')
    prop_clump = pickle.load(f)
    thres_v = pickle.load(f)
    f.close()



    #select a subregion for the clumps
    if( radius is not None):
        if(center is None):
            print 'center should be defined as a 3d array'
            return
        if(len(center) != 3 ):
            print 'center should be defined as a 3d array'
            return
        #calculate the position with respect to specified center
        xx = prop_clump['pos_n_max'][:,0]-center[0]
        yy = prop_clump['pos_n_max'][:,1]-center[1]
        zz = prop_clump['pos_n_max'][:,2]-center[2]

        #    rad_v = np.sqrt(xx**2+yy**2+zz**2)
        # use  projection along z axis for now
        rad_v = np.sqrt(xx**2+yy**2)

        mask = rad_v < radius

        #select the clumps of interest
        for key in prop_clump.keys():
            prop_clump[key] = prop_clump[key][mask]




    mask = prop_clump['mass'] > 0.
    mass_vlog = np.log10(prop_clump['mass'][mask])
    if(log_min is None):
        log_min=min(mass_vlog)
    if(log_max is None):
        log_max=max(mass_vlog)



    #histogram of cell number
    ncell_vlog = np.log10(prop_clump['ncell'])
    log_cmin=min(ncell_vlog)
    log_cmax=max(ncell_vlog)
    nbin=50.
    width=(log_cmax-log_cmin)/nbin
    P.clf()
    hist_ncell , hist_edges = P.histogram(ncell_vlog,bins=nbin,range=(log_cmin,log_cmax))
    P.bar(hist_edges[:-1],hist_ncell,width)
    P.yscale('log', nonposy='clip')
    P.xlabel(r'$log(ncell)$')
    P.ylabel(r'$log(N_c)$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist_ncell_'+str(num).zfill(5)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_hist_ncell_'+str(num).zfill(5)+'.pdf')
    P.savefig(directory_out+name+subname+'_hist_ncell_'+str(num).zfill(5)+'.png')



    #series of mass spectra (various thresholds)

    #cut on the number of points
    mask = prop_clump['mass'] > 0.
    mask_n = prop_clump['ncell'] > 100
    mask = mask * mask_n
    mass_vlog = np.log10(prop_clump['mass'][mask])
    nbin=50.
    width=(log_max-log_min)/nbin
    P.clf()
    hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
    P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='k',linewidth=2)


    P.plot([0.,3.],[1.e4,1.e1],color='k')
#    P.plot([0.,3.],[1.e4,1.e0],color='k')

#    P.plot([0.,3.],[1.e4,1.e1],color='k')
    P.yscale('log', nonposy='clip')
    P.xlabel(r'$log(M) \, (M_ \odot)$')
    P.ylabel(r'$dN/dlogM$')

    P.legend( (r'$M^{-1}$','all structures')  ,loc='upper right')
#    P.title=('all structures')

    if(ps):
        P.savefig(directory_out+name+subname+'_hist_mass_'+str(num).zfill(5)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_hist_mass_'+str(num).zfill(5)+'.pdf')
    P.savefig(directory_out+name+subname+'_hist_mass_'+str(num).zfill(5)+'.png')



    #cut on the mass to flux ratio, small mu
    mask = prop_clump['mass'] > 0.
    mask_n = prop_clump['mu'] < 1.
    mask_mag = prop_clump['PB'] != 0.
    mask = mask * mask_n * mask_mag
    mask_n = prop_clump['ncell'] > 100
    mask = mask_n * mask
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass'][mask])
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='k',linewidth=2)

        P.plot([0.,3.],[1.e4,1.e1],color='k')
        P.plot([0.,3.],[1.e4,1.e0],color='k')

        P.yscale('log', nonposy='clip')
        P.xlabel(r'$log(M) \, (M_ \odot)$')
        P.ylabel(r'$dN/dlogM$')

        P.legend( (r'$M^{-1}$',r'$M^{-1.3}$',r'$\mu<1$')  ,loc='upper right')



#        P.title('subcritical cores')

        if(ps):
            P.savefig(directory_out+name+subname+'_hist_mass_mulow_'+str(num).zfill(5)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_hist_mass_mulow_'+str(num).zfill(5)+'.pdf')
        P.savefig(directory_out+name+subname+'_hist_mass_mulow_'+str(num).zfill(5)+'.png')



    #cut on the mass to flux ratio, big mu
    mask = prop_clump['mass'] > 0.
    mask_n = prop_clump['mu'] > 1.
    mask_mag = prop_clump['PB'] != 0.
    mask = mask * mask_n * mask_mag
    mass_vlog = np.log10(prop_clump['mass'][mask])
    mask_n = prop_clump['ncell'] > 100
    mask = mask_n * mask
    nbin=50.
    width=(log_max-log_min)/nbin
    P.clf()
    hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
    P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='k',linewidth=2)
#    P.plot([0.,3.],[1.e4,1.e1],color='k')
    P.yscale('log', nonposy='clip')
    P.xlabel(r'$log(M) \, (M_ \odot)$')
    P.ylabel(r'$dN/dlogM$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist_mass_mucut_'+str(num).zfill(5)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_hist_mass_mucut_'+str(num).zfill(5)+'.pdf')
    P.savefig(directory_out+name+subname+'_hist_mass_mucut_'+str(num).zfill(5)+'.png')


    #add to the previous plot cut on the mass to smallest flux ratio, mu
    mask = prop_clump['mass'] > 0.
    mask_n = prop_clump['mu'] > 1./3.
    mask_mag = prop_clump['PB'] != 0.
    mask = mask * mask_n * mask_mag
    mask_n = prop_clump['ncell'] > 100
    mask = mask_n * mask
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass'][mask])
        nbin=50.
        width=(log_max-log_min)/nbin
#        P.clf()
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width,color='none',linewidth=2,edgecolor='r',linestyle='dashed')

        P.legend((r'$\mu>1$',r'$\mu>0.3$')  ,loc='upper right')

        P.plot([0.,3.],[1.e4,1.e1],color='k')
    
#        P.yscale('log', nonposy='clip')
#        P.xlabel(r'$log(M)$ (M_\odot)')
#        P.ylabel(r'$log(dN/dlogM)$')
        if(ps):
            P.savefig(directory_out+name+subname+'_hist_mass_mucut_0.3_'+str(num).zfill(5)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_hist_mass_mucut_0.3_'+str(num).zfill(5)+'.pdf')
        P.savefig(directory_out+name+subname+'_hist_mass_mucut_0.3_'+str(num).zfill(5)+'.png')


#    pdb.set_trace() # debug mode     




    #cut the high virial parameter
    mask = prop_clump['mass'] > 0.
    mask_vir = np.abs(prop_clump['grav_vir']) > 1.e-10
    mask = mask * mask_vir
    mask_n = np.abs(prop_clump['ener_kin']/prop_clump['grav_vir']) < 1./2.
    mask = mask * mask_n
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass'][mask])
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='k',linewidth=2)
        P.yscale('log', nonposy='clip')
        P.xlabel(r'$log(M) \, (M_ \odot)$')
        P.ylabel(r'$dN/dlogM$')
        if(ps):
            P.savefig(directory_out+name+subname+'_hist_mass_vircut_'+str(num).zfill(5)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_hist_mass_vircut_'+str(num).zfill(5)+'.pdf')
        P.savefig(directory_out+name+subname+'_hist_mass_vircut_'+str(num).zfill(5)+'.png')


    #cut the low virial parameter
    mask = prop_clump['mass'] > 0.
    mask_vir = np.abs(prop_clump['grav_vir']) > 1.e-10
    mask = mask * mask_vir
    mask_n = np.abs(prop_clump['ener_kin']/prop_clump['grav_vir']) > 1./2.
    mask = mask * mask_n
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass'][mask])
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='k',linewidth=2)
        P.yscale('log', nonposy='clip')
        P.xlabel(r'$log(M) \, (M_ \odot)$')
        P.ylabel(r'$dN/dlogM$')
        if(ps):
            P.savefig(directory_out+name+subname+'_hist_mass_virhigh_'+str(num).zfill(5)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_hist_mass_virhigh_'+str(num).zfill(5)+'.pdf')
        P.savefig(directory_out+name+subname+'_hist_mass_virhigh_'+str(num).zfill(5)+'.png')



    #cut the low thermal virial parameter
    mask = prop_clump['mass'] > 0.
    mask_vir = np.abs(prop_clump['grav_vir']) > 1.e-10
    mask = mask * mask_vir
    mask_n = np.abs(prop_clump['therm_vir']/prop_clump['grav_vir']) > 1./2.
    mask = mask * mask_n
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass'][mask])
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='k',linewidth=2)
        P.yscale('log', nonposy='clip')
        P.xlabel(r'$log(M) \, (M_ \odot)$')
        P.ylabel(r'$dN/dlogM$')
        if(ps):
            P.savefig(directory_out+name+subname+'_hist_mass_virthhigh_'+str(num).zfill(5)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_hist_mass_virthhigh_'+str(num).zfill(5)+'.pdf')
        P.savefig(directory_out+name+subname+'_hist_mass_virthhigh_'+str(num).zfill(5)+'.png')



    #cut the high thermal virial parameter
    mask = prop_clump['mass'] > 0.
    mask_vir = np.abs(prop_clump['grav_vir']) > 1.e-10
    mask = mask * mask_vir
    mask_n = np.abs(prop_clump['therm_vir']/prop_clump['grav_vir']) < 1./2.
    mask = mask * mask_n
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass'][mask])
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='k',linewidth=2)
        P.plot([0.,3.],[1.e4,1.e1],color='k')
        P.yscale('log', nonposy='clip')
        P.xlabel(r'$log(M) \, (M_ \odot)$')
        P.ylabel(r'$dN/dlogM$')
        if(ps):
            P.savefig(directory_out+name+subname+'_hist_mass_virthcut_'+str(num).zfill(5)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_hist_mass_virthcut_'+str(num).zfill(5)+'.pdf')
        P.savefig(directory_out+name+subname+'_hist_mass_virthcut_'+str(num).zfill(5)+'.png')


    #add to the previous plot the mass weighted histogram selecting points only over certain density thresholds
    mask = prop_clump['mass_0'] > 0.
    mask_vir = np.abs(prop_clump['grav_vir']) > 1.e-10
    mask = mask * mask_vir
    mask_n = np.abs(prop_clump['therm_vir']/prop_clump['grav_vir']) < 1./2.
    mask_n1 = prop_clump['ncell'] > 100
    mask = mask * mask_n * mask_n1
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass_0'][mask])
        nbin=50.
        width=(log_max-log_min)/nbin
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='k',linestyle='solid',linewidth=2)


    mask = prop_clump['mass_1'] > 0.
    mask_vir = np.abs(prop_clump['grav_vir']) > 1.e-10
    mask = mask * mask_vir
    mask_n = np.abs(prop_clump['therm_vir']/prop_clump['grav_vir']) < 1./2.
    mask_n1 = prop_clump['ncell'] > 100
    mask = mask * mask_n * mask_n1
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass_1'][mask])
        nbin=50.
        width=(log_max-log_min)/nbin
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='r',linestyle='solid',linewidth=2)

    mask = prop_clump['mass_2'] > 0.
    mask_vir = np.abs(prop_clump['grav_vir']) > 1.e-10
    mask = mask * mask_vir
    mask_n = np.abs(prop_clump['therm_vir']/prop_clump['grav_vir']) < 1./2.
    mask_n1 = prop_clump['ncell'] > 100
    mask = mask * mask_n * mask_n1
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass_2'][mask])
        nbin=50.
        width=(log_max-log_min)/nbin
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='m',linestyle='solid',linewidth=2)


    mask = prop_clump['mass_3'] > 0.
    mask_vir = np.abs(prop_clump['grav_vir']) > 1.e-10
    mask = mask * mask_vir
    mask_n = np.abs(prop_clump['therm_vir']/prop_clump['grav_vir']) < 1./2.
    mask_n1 = prop_clump['ncell'] > 100
    mask = mask * mask_n * mask_n1
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass_3'][mask])
        nbin=50.
        width=(log_max-log_min)/nbin
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='b',linestyle='solid',linewidth=2)


        P.yscale('log', nonposy='clip')
        P.xlabel(r'$log(M) \, (M_\odot)$')
        P.ylabel(r'$dN/dlogM$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist_mass_virthcut_mult_'+str(num).zfill(5)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_hist_mass_virthcut_mult_'+str(num).zfill(5)+'.pdf')
    P.savefig(directory_out+name+subname+'_hist_mass_virthcut_mult_'+str(num).zfill(5)+'.png')




    #cut the high thermal virial parameter weighted by freefall time
    #display the mean density vs mass
    mask = prop_clump['mass'] > 0.
    mask_vir = np.abs(prop_clump['grav_vir']) > 1.e-10
    mask = mask * mask_vir
    mask_n = np.abs(prop_clump['therm_vir']/prop_clump['grav_vir']) < 1./2.
    mask = mask * mask_n
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass'][mask])
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()

        rad = (prop_clump['size_iner'][:,0]*prop_clump['size_iner'][:,1]*prop_clump['size_iner'][:,2])**(0.3333) #an estimate of the radius
        rho_mean = prop_clump['mass'][mask] / rad[mask]**3  
        o_ff_time = np.sqrt(rho_mean) #invert of the free fall time

        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max),weights=o_ff_time)
        P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='k',linewidth=2)
        P.plot([0.,3.],[1.e4,1.e1],color='k')
        P.yscale('log', nonposy='clip')
        P.xlabel(r'$log(M) \, (M_ \odot)$')
        P.ylabel(r'$dN/dlogM$')
        if(ps):
            P.savefig(directory_out+name+subname+'_hist_mass_fftw_virthcut_'+str(num).zfill(5)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_hist_mass_fftw_virthcut_'+str(num).zfill(5)+'.pdf')
        P.savefig(directory_out+name+subname+'_hist_mass_fftw_virthcut_'+str(num).zfill(5)+'.png')


#        P.clf()
#        P.plot(np.log10(prop_clump['mass'][mask]),np.log10(rho_mean),'.',markersize=mark)
#        P.xlabel(r'$log(M) \, (M_ \odot)$')
#        P.ylabel(r'$log(n) \, (cm^{-3})$')
#        if(ps):
#            P.savefig(directory_out+name+subname+'_mass_nmean_'+str(num).zfill(5)+'.ps')
#        if(pdf):
#            P.savefig(directory_out+name+subname+'_mass_nmean_'+str(num).zfill(5)+'.pdf')
#        P.savefig(directory_out+name+subname+'_mass_nmean_'+str(num).zfill(5)+'.png')





        #take care with mask; mask_n and rho_mean

    #2D distogram of J.B/|J||B| vs mass, cut the low mu parameter and the large mean densities
    mask = prop_clump['mass'] > 0.
    mask_mu = np.abs(prop_clump['mu']) > 1.
    mask = mask * mask_mu
    mask_n = np.abs(prop_clump['therm_vir']/prop_clump['grav_vir']) < 1./2.
    mask = mask * mask_n
    mask_mag = prop_clump['PB'] > 0.
    mask = mask *  mask_mag
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass'][mask])
        nbin=50.
        width=(log_max-log_min)/nbin

        rad = (prop_clump['size_iner'][:,0]*prop_clump['size_iner'][:,1]*prop_clump['size_iner'][:,2])**(0.3333) #an estimate of the radius
        rho_mean = prop_clump['mass'][mask] / rad[mask]**3  

        alpha = (prop_clump['j_v'][:,0]*prop_clump['B'][:,0]+prop_clump['j_v'][:,1]*prop_clump['B'][:,1]+prop_clump['j_v'][:,2]*prop_clump['B'][:,2])
        alpha = alpha /  np.sqrt(prop_clump['j_v'][:,0]**2+prop_clump['j_v'][:,1]**2+prop_clump['j_v'][:,2]**2)
        alpha = alpha /  np.sqrt(prop_clump['B'][:,0]**2+prop_clump['B'][:,1]**2+prop_clump['B'][:,2]**2)
        alpha = alpha[mask]

        mask_n = rho_mean < 1.e6
        if (np.sum(mask_n) > 0):
            nbins = nbins_glob 
            mass_nmean , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask][mask_n]),(alpha[mask_n]),bins=nbin) #,weights=prop_clump['mass'][mask]) 
            P.clf()
            im = P.imshow(np.log10(mass_nmean.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
            cbar = P.colorbar(im)                                                                                                          
#        cbar.set_label(r'$log(M) \, (M_\odot)$')
            P.xlabel(r'$log(M) \, (M_ \odot)$')
            P.ylabel(r'$\alpha$')                                          
            if(ps):
                P.savefig(directory_out+name+subname+'_hist2D_mass_JB'+'_'+str(num)+'.ps')
            if(pdf):
                P.savefig(directory_out+name+subname+'_hist2D_mass_JB'+'_'+str(num)+'.pdf')
            P.savefig(directory_out+name+subname+'_hist2D_mass_JB'+'_'+str(num)+'.png')



    #2D distogram of |J| vs mass, cut the low mu parameter and the large mean densities
    mask = prop_clump['mass'] > 0.
    mask_mu = np.abs(prop_clump['mu']) > 1.
    mask = mask * mask_mu
    mask_n = np.abs(prop_clump['therm_vir']/prop_clump['grav_vir']) < 1./2.
    mask = mask * mask_n
    mask_mag = prop_clump['PB'] > 0.
    mask = mask *  mask_mag
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass'][mask])
        nbin=50.
        width=(log_max-log_min)/nbin

        rad = (prop_clump['size_iner'][:,0]*prop_clump['size_iner'][:,1]*prop_clump['size_iner'][:,2])**(0.3333) #an estimate of the radius
        rho_mean = prop_clump['mass'][mask] / rad[mask]**3  

        mask_n = rho_mean < 1.e6
        if (np.sum(mask_n) > 0):
            nbins = nbins_glob 
            mass_J , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask][mask_n]),np.log10(prop_clump['j'][mask][mask_n]),bins=nbin) #,weights=prop_clump['mass'][mask]) 
            P.clf()
            im = P.imshow(np.log10(mass_J.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
            cbar = P.colorbar(im)                                                                                                          
#        cbar.set_label(r'$log(M) \, (M_\odot)$')
            cbar.set_label('number of cores')
            P.xlabel(r'$log(M) \, (M_ \odot)$')
            P.ylabel(r'$log(J) \, (pc \, km \, s^{-1})$')                                          
            if(ps):
                P.savefig(directory_out+name+subname+'_hist2D_mass_Jnorm'+'_'+str(num)+'.ps')
            if(pdf):
                P.savefig(directory_out+name+subname+'_hist2D_mass_Jnorm'+'_'+str(num)+'.pdf')
            P.savefig(directory_out+name+subname+'_hist2D_mass_Jnorm'+'_'+str(num)+'.png')




        mask_n = rho_mean < 1.e6
        if (np.sum(mask_n) > 0):
            nbins = nbins_glob 
            sig_J , xedges, yedges = np.histogram2d(np.log10(prop_clump['sig'][mask][mask_n]/np.sqrt(3.)),np.log10(prop_clump['j'][mask][mask_n]),bins=nbin) #,weights=prop_clump['mass'][mask]) 
            P.clf()
            im = P.imshow(np.log10(sig_J.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
            cbar = P.colorbar(im)                                                                                                          
#        cbar.set_label(r'$log(M) \, (M_\odot)$')
            cbar.set_label('number of cores')
            P.xlabel(r'$log(\sigma) \, (km \, s^{-1})$')
            P.ylabel(r'$log(J) \, (pc \, km \, s^{-1})$')                                          
            if(ps):
                P.savefig(directory_out+name+subname+'_hist2D_sig_Jnorm'+'_'+str(num)+'.ps')
            if(pdf):
                P.savefig(directory_out+name+subname+'_hist2D_sig_Jnorm'+'_'+str(num)+'.pdf')
            P.savefig(directory_out+name+subname+'_hist2D_sig_Jnorm'+'_'+str(num)+'.png')





    #cut the high thermal virial parameter and the large mean densities
    mask = prop_clump['mass'] > 0.
    mask_vir = np.abs(prop_clump['grav_vir']) > 1.e-10
    mask = mask * mask_vir
    mask_n = np.abs(prop_clump['therm_vir']/prop_clump['grav_vir']) < 1./2.
    mask = mask * mask_n
    mask_n = prop_clump['ncell'] > 100
    mask = mask_n * mask
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass'][mask])
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()

        rad = (prop_clump['size_iner'][:,0]*prop_clump['size_iner'][:,1]*prop_clump['size_iner'][:,2])**(0.3333) #an estimate of the radius
        rho_mean = prop_clump['mass'][mask] / rad[mask]**3  
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='k',linewidth=2)

        mask_n = rho_mean < 1.e5
        if (np.sum(mask_n) > 0):
            hist_mass , hist_edges = P.histogram(mass_vlog[mask_n],bins=nbin,range=(log_min,log_max))
            P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='m',linewidth=2)

        mask_n = rho_mean < 1.e6
        if (np.sum(mask_n) > 0):
            hist_mass , hist_edges = P.histogram(mass_vlog[mask_n],bins=nbin,range=(log_min,log_max))
            P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='r',linewidth=2)


        mask_n = rho_mean < 1.e7
        if (np.sum(mask_n) > 0):
            hist_mass , hist_edges = P.histogram(mass_vlog[mask_n],bins=nbin,range=(log_min,log_max))
            P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='b',linewidth=2)

        P.plot([0.,3.],[1.e4,1.e1],color='k')
        P.yscale('log', nonposy='clip')
        P.xlabel(r'$log(M) \, (M_ \odot)$')
        P.ylabel(r'$dN/dlogM$')
        if(ps):
            P.savefig(directory_out+name+subname+'_hist_mass_virth_hd_cut_'+str(num).zfill(5)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_hist_mass_virth_hd_cut_'+str(num).zfill(5)+'.pdf')
        P.savefig(directory_out+name+subname+'_hist_mass_virth_hd_cut_'+str(num).zfill(5)+'.png')



        #reduced version of the mass spectrum selected by mean density
        P.clf()
        mask_n = rho_mean < 1.e5
        if (np.sum(mask_n) > 0):
            hist_mass , hist_edges = P.histogram(mass_vlog[mask_n],bins=nbin,range=(log_min,log_max))
            P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='k',linewidth=2)

        mask_n = rho_mean < 1.e6
        if (np.sum(mask_n) > 0):
            hist_mass , hist_edges = P.histogram(mass_vlog[mask_n],bins=nbin,range=(log_min,log_max))
            P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='b',linewidth=2)



        P.yscale('log', nonposy='clip')
        P.xlabel(r'$log(M) \, (M_ \odot)$')
        P.ylabel(r'$dN/dlogM$')

        P.legend((r'$\alpha_{th} < 1 \, n_{mean<10^5}$', r'$\alpha_{th} < 1 \, n_{mean<10^6}$')  ,loc='upper right')

        P.plot([0.,3.],[1.e4,1.e1],color='k')
        P.plot([0.,3.],[1.e4,1.e0],color='k')

        if(ps):
            P.savefig(directory_out+name+subname+'_hist_mass_virth_hd1_cut_'+str(num).zfill(5)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_hist_mass_virth_hd1_cut_'+str(num).zfill(5)+'.pdf')
        P.savefig(directory_out+name+subname+'_hist_mass_virth_hd1_cut_'+str(num).zfill(5)+'.png')




    #cut the low mu parameter and the large mean densities
    mask = prop_clump['mass'] > 0.
    mask_mu = np.abs(prop_clump['mu']) > 1.
    mask = mask * mask_mu
    mask_n = prop_clump['ncell'] > 100
    mask = mask_n * mask
#    mask_n = np.abs(prop_clump['therm_vir']/prop_clump['grav_vir']) < 2.
#    mask = mask * mask_n
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass'][mask])
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()

        rad = (prop_clump['size_iner'][:,0]*prop_clump['size_iner'][:,1]*prop_clump['size_iner'][:,2])**(0.3333) #an estimate of the radius
        rho_mean = prop_clump['mass'][mask] / rad[mask]**3  


        P.clf()
        mask_n = rho_mean < 1.e5
        if (np.sum(mask_n) > 0):
            hist_mass , hist_edges = P.histogram(mass_vlog[mask_n],bins=nbin,range=(log_min,log_max))
            P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='k',linewidth=2)

        mask_n = rho_mean < 1.e6
        if (np.sum(mask_n) > 0):
            hist_mass , hist_edges = P.histogram(mass_vlog[mask_n],bins=nbin,range=(log_min,log_max))
            P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='b',linewidth=2)


        P.yscale('log', nonposy='clip')
        P.xlabel(r'$log(M) \, (M_ \odot)$')
        P.ylabel(r'$dN/dlogM$')


        P.legend((r'$\mu>1 \, n_{mean<10^5}$', r'$\mu>1 \, n_{mean<10^6}$')  ,loc='upper right')

        P.plot([0.,3.],[1.e4,1.e1],color='k')
        P.plot([0.,3.],[1.e4,1.e0],color='k')


        if(ps):
            P.savefig(directory_out+name+subname+'_hist_mass_mu_hd1_cut_'+str(num).zfill(5)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_hist_mass_mu_hd1_cut_'+str(num).zfill(5)+'.pdf')
        P.savefig(directory_out+name+subname+'_hist_mass_mu_hd1_cut_'+str(num).zfill(5)+'.png')




    ### mass histogram of clumps having a maximum density over some threshold
    mask = prop_clump['mass_1'] > 0.
    mask_n = prop_clump['ncell'] > 100
    mask = mask * mask_n
    mask1 = prop_clump['n_max'] > thres_v[2]
    mask = mask*mask1
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass_1'][mask])
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='k',linewidth=2)

    mask = prop_clump['mass_1'] > 0.
    mask_n = prop_clump['ncell'] > 100
    mask = mask * mask_n
    mask1 = prop_clump['n_max'] > thres_v[4]
    mask = mask*mask1
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass_3'][mask])
        nbin=50.
        width=(log_max-log_min)/nbin

        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width,color='none',linewidth=2,edgecolor='b')


    mask = prop_clump['mass_1'] > 0.
    mask_n = prop_clump['ncell'] > 100
    mask = mask * mask_n
    mask1 = prop_clump['n_max'] > thres_v[5]
    mask = mask*mask1
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass_5'][mask])
        nbin=50.
        width=(log_max-log_min)/nbin
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width,color='none',linewidth=2,edgecolor='r')

        P.yscale('log', nonposy='clip')
        P.xlabel(r'$log(M)$ (M_\odot)')
        P.ylabel(r'$dN/dlogM$')

    if(ps):
        P.savefig(directory_out+name+subname+'_hist_mass_'+str(np.log10(thres_v[1]))+'_denscut_'+str(num).zfill(5)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_hist_mass_'+str(np.log10(thres_v[1]))+'_denscut_'+str(num).zfill(5)+'.pdf')
    P.savefig(directory_out+name+subname+'_hist_mass_'+str(np.log10(thres_v[1]))+'_denscut_'+str(num).zfill(5)+'.png')






    #distribution of M as a function of size calculated from inertia matrix
    P.clf()
    rad = (prop_clump['size_iner'][:,0]*prop_clump['size_iner'][:,1]*prop_clump['size_iner'][:,2])**(0.3333) #an estimate of the radius
    mask = rad > 0.
    mask_1 = prop_clump['ncell'] > 10
    mask = mask * mask_1

    P.plot(np.log10(rad[mask]),np.log10(prop_clump['mass'][mask]),'.',markersize=mark)

    P.xlabel(r'$log(R)$ (pc)')
    P.ylabel(r'$log(M) \, (M_ \odot)$')
    if(ps):
        P.savefig(directory_out+name+subname+'_R_mass_'+str(num).zfill(5)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_R_mass_'+str(num).zfill(5)+'.pdf')
    P.savefig(directory_out+name+subname+'_R_mass_'+str(num).zfill(5)+'.png')


    nbins = nbins_glob 
    rad_mass , xedges, yedges = np.histogram2d(np.log10(rad[mask]),np.log10(prop_clump['mass'][mask]),bins=nbin,weights=prop_clump['mass'][mask]) 
    P.clf()
    im = P.imshow(np.log10(rad_mass.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(R)$ (pc)')
    P.ylabel(r'$log(M) \, (M_ \odot)$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_rad_mass'+'_'+str(num)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_hist2D_rad_mass'+'_'+str(num)+'.pdf')
    P.savefig(directory_out+name+subname+'_hist2D_rad_mass'+'_'+str(num)+'.png')




    mask = rad > 0.
    mask_n = prop_clump['rho_mean'] < 1.e5
    mask = mask*mask_n

    nbins = nbins_glob 
    rad_mass , xedges, yedges = np.histogram2d(np.log10(rad[mask]),np.log10(prop_clump['mass'][mask]),bins=nbin,weights=prop_clump['mass'][mask]) 
    P.clf()
    im = P.imshow(np.log10(rad_mass.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(R)$ (pc)')
    P.ylabel(r'$log(M) \, (M_ \odot)$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_rad_mass_nmean'+'_'+str(num)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_hist2D_rad_mass_nmean'+'_'+str(num)+'.pdf')
    P.savefig(directory_out+name+subname+'_hist2D_rad_mass_nmean'+'_'+str(num)+'.png')




    mask = rad > 0.
    mask_n = prop_clump['rho_mean'] < 1.e5
    mask = mask*mask_n
    mask_n = np.abs(prop_clump['therm_vir']/prop_clump['grav_vir']) < 1./2.
    mask = mask*mask_n

    nbins = nbins_glob 
    rad_mass , xedges, yedges = np.histogram2d(np.log10(rad[mask]),np.log10(prop_clump['mass'][mask]),bins=nbin) #,weights=prop_clump['mass'][mask]) 
    P.clf()
    im = P.imshow(np.log10(rad_mass.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
#    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    cbar.set_label('number of cores')
    P.xlabel(r'$log(R)$ (pc)')
    P.ylabel(r'$log(M) \, (M_ \odot)$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_rad_mass_nmean_vircut'+'_'+str(num)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_hist2D_rad_mass_nmean_vircut'+'_'+str(num)+'.pdf')
    P.savefig(directory_out+name+subname+'_hist2D_rad_mass_nmean_vircut'+'_'+str(num)+'.png')



    mask = rad > 0.
    mask_n = np.abs(prop_clump['therm_vir']/prop_clump['grav_vir']) < 1./2.
    mask = mask*mask_n

    nbins = nbins_glob 
    rad_mass , xedges, yedges = np.histogram2d(np.log10(rad[mask]),np.log10(prop_clump['mass'][mask]),bins=nbin) #,weights=prop_clump['mass'][mask]) 
    P.clf()
    im = P.imshow(np.log10(rad_mass.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
#    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    cbar.set_label('number of cores')
    P.xlabel(r'$log(R)$ (pc)')
    P.ylabel(r'$log(M) \, (M_ \odot)$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_rad_mass_vircut'+'_'+str(num)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_hist2D_rad_mass_vircut'+'_'+str(num)+'.pdf')
    P.savefig(directory_out+name+subname+'_hist2D_rad_mass_vircut'+'_'+str(num)+'.png')





    #distribution of sigma as a function of size
    P.clf()
    rad = (3./4./np.pi*prop_clump['vol'])**(0.3333) #an estimate of the radius    

    P.plot(np.log10(rad),np.log10(prop_clump['sig']/np.sqrt(3.)),'.',markersize=mark)

    P.xlabel(r'$log(R)$ (pc)')
    P.ylabel(r'$log(\sigma) \, (km \, s^{-1})$')
    if(ps):
        P.savefig(directory_out+name+subname+'_R_sig_'+str(num).zfill(5)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_R_sig_'+str(num).zfill(5)+'.pdf')
    P.savefig(directory_out+name+subname+'_R_sig_'+str(num).zfill(5)+'.png')


    nbins=nbins_glob
    mask = prop_clump['sig'] > 0.
    rad_sig , xedges, yedges = np.histogram2d(np.log10(rad[mask]),np.log10(prop_clump['sig'][mask]/np.sqrt(3.)),bins=nbin,weights=prop_clump['mass'][mask]) #,range=[[xmin,xmax],[ymin,ymax]])

#    mask = rad_sig == 0.
#    rad_sig[mask] = np.min(rad_sig)
    P.clf()
    im = P.imshow(np.log10(rad_sig.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(R)$ (pc)')
    P.ylabel(r'$log(\sigma) \, (km \, s^{-1})$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_rad_sig'+'_'+str(num)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_hist2D_rad_sig'+'_'+str(num)+'.pdf')
    P.savefig(directory_out+name+subname+'_hist2D_rad_sig'+'_'+str(num)+'.png')





    #distribution of sigma*col_dens as a function of size
    P.clf()
    col_dens = prop_clump['mass'] / np.pi / rad**2  #column density expressed in Ms / pc^2
    P.plot(np.log10(rad*col_dens),np.log10(prop_clump['sig']/np.sqrt(3.)),'.',markersize=mark)
    P.xlabel(r'$log(R \Sigma) (M _\odot pc^{-1})$')
    P.ylabel(r'$log(\sigma) \, (km \, s^{-1})$')
    if(ps):
        P.savefig(directory_out+name+subname+'_RCD_sig_'+str(num).zfill(5)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_RCD_sig_'+str(num).zfill(5)+'.pdf')
    P.savefig(directory_out+name+subname+'_RCD_sig_'+str(num).zfill(5)+'.png')

    nbins=nbins_glob
    mask = prop_clump['sig'] > 0.
    CDrad_sig , xedges, yedges = np.histogram2d(np.log10(rad[mask]*col_dens[mask]),np.log10(prop_clump['sig'][mask]/np.sqrt(3.)),bins=nbin,weights=prop_clump['mass'][mask])
    P.clf()
    im = P.imshow(np.log10(CDrad_sig.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(R \Sigma) (M_\odot pc^{-1})$')
    P.ylabel(r'$log(\sigma) \, (km \, s^{-1})$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_RCD_sig'+'_'+str(num)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_hist2D_RCD_sig'+'_'+str(num)+'.pdf')
    P.savefig(directory_out+name+subname+'_hist2D_RCD_sig'+'_'+str(num)+'.png')




    #distribution of sigma as a function of mass
    P.clf()
    mass_vlog = np.log10(prop_clump['mass'])
    P.plot(mass_vlog,np.log10(prop_clump['sig']/np.sqrt(3.)),'.',markersize=mark)

    P.xlabel(r'$log(M)$ (M_\odot)')
    P.ylabel(r'$log(\sigma) \, (km \, s^{-1})$')
    if(ps):
        P.savefig(directory_out+name+subname+'_mass_sig_'+str(num).zfill(5)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_mass_sig_'+str(num).zfill(5)+'.pdf')
    P.savefig(directory_out+name+subname+'_mass_sig_'+str(num).zfill(5)+'.png')

    nbins=nbins_glob
    mask = prop_clump['sig'] > 0.
    mass_sig , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(prop_clump['sig'][mask]/np.sqrt(3.)),bins=nbin,weights=prop_clump['mass'][mask])
    P.clf()
    im = P.imshow(np.log10(mass_sig.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(M) (M_\odot)$')
    P.ylabel(r'$log(\sigma) \, (km \, s^{-1})$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_mass_sig'+'_'+str(num)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_hist2D_mass_sig'+'_'+str(num)+'.pdf')
    P.savefig(directory_out+name+subname+'_hist2D_mass_sig'+'_'+str(num)+'.png')


    nbins=nbins_glob
    mask = prop_clump['sig'] > 0.

    mask_n = prop_clump['rho_mean'] < 1.e5
    mask = mask*mask_n

    mass_sig , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(prop_clump['sig'][mask]/np.sqrt(3.)),bins=nbin,weights=prop_clump['mass'][mask])
    P.clf()
    im = P.imshow(np.log10(mass_sig.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(M) (M_\odot)$')
    P.ylabel(r'$log(\sigma) \, (km \, s^{-1})$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_mass_sig_nmean_low'+'_'+str(num)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_hist2D_mass_sig_nmean_low'+'_'+str(num)+'.pdf')
    P.savefig(directory_out+name+subname+'_hist2D_mass_sig_nmean_low'+'_'+str(num)+'.png')





    nbins=nbins_glob
    mask = prop_clump['sig'] > 0.

    mask_n = prop_clump['rho_mean'] < 1.e5
    mask = mask*mask_n

    mask_n = np.abs(prop_clump['therm_vir']/prop_clump['grav_vir']) < 1.
    mask = mask*mask_n

    mass_sig , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(prop_clump['sig'][mask]/np.sqrt(3.)),bins=nbin)
    P.clf()
    im = P.imshow(np.log10(mass_sig.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
#    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    cbar.set_label('number of cores')                                          
    P.xlabel(r'$log(M) (M_\odot)$')
    P.ylabel(r'$log(\sigma) \, (km \, s^{-1})$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_mass_sig_nmean_low_virth'+'_'+str(num)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_hist2D_mass_sig_nmean_low_virth'+'_'+str(num)+'.pdf')
    P.savefig(directory_out+name+subname+'_hist2D_mass_sig_nmean_low_virth'+'_'+str(num)+'.png')







    #distribution of volume weighted Mach as a function of mass
    P.clf()
    mass_vlog = np.log10(prop_clump['mass'])
    P.plot(mass_vlog,np.log10(prop_clump['Mach']),'.',markersize=mark)

    P.xlabel(r'$log(M) (M_\odot)$')
    P.ylabel(r'$log({\cal M})$')
    if(ps):
        P.savefig(directory_out+name+subname+'_mass_Mach_'+str(num).zfill(5)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_mass_Mach_'+str(num).zfill(5)+'.pdf')
    P.savefig(directory_out+name+subname+'_mass_Mach_'+str(num).zfill(5)+'.png')



    #distribution of mass weighted Mach as a function of mass
    P.clf()
    mass_vlog = np.log10(prop_clump['mass'])
    P.plot(mass_vlog,np.log10(prop_clump['Mach_mw']),'.',markersize=mark)

    P.xlabel(r'$log(M)$ (M_\odot)')
    P.ylabel(r'$log({\cal M})$')
    if(ps):
        P.savefig(directory_out+name+subname+'_mass_Mach_mw_'+str(num).zfill(5)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_mass_Mach_mw_'+str(num).zfill(5)+'.pdf')
    P.savefig(directory_out+name+subname+'_mass_Mach_mw_'+str(num).zfill(5)+'.png')

    nbins=nbins_glob
    mask = prop_clump['sig'] > 0.
    mass_mach , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(prop_clump['Mach_mw'][mask]),bins=nbin,weights=prop_clump['mass'][mask])
    P.clf()
    im = P.imshow(np.log10(mass_mach.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(M) (M_\odot)$')
    P.ylabel(r'$log({\cal M})$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_mass_Mach_mw'+'_'+str(num)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_hist2D_mass_Mach_mw'+'_'+str(num)+'.pdf')
    P.savefig(directory_out+name+subname+'_hist2D_mass_Mach_mw'+'_'+str(num)+'.png')



    nbins=nbins_glob
    mask = prop_clump['sig'] > 0.

    mask_n = prop_clump['rho_mean'] < 1.e5
    mask = mask*mask_n

    mass_mach , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(prop_clump['Mach_mw'][mask]),bins=nbin,weights=prop_clump['mass'][mask])
    P.clf()
    im = P.imshow(np.log10(mass_mach.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(M) (M_\odot)$')
    P.ylabel(r'$log({\cal M})$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_mass_Mach_mw_nmean_low'+'_'+str(num)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_hist2D_mass_Mach_mw_nmean_low'+'_'+str(num)+'.pdf')
    P.savefig(directory_out+name+subname+'_hist2D_mass_Mach_mw_nmean_low'+'_'+str(num)+'.png')



    #distribution of Alfvenic Mach nmber as a function of mass
    mask = prop_clump['PB'] > 0.
    if(np.sum(mask) > 0):
        P.clf()
        mass_vlog = np.log10(prop_clump['mass'])
        P.plot(mass_vlog,np.log10(prop_clump['Mach_alfv']),'.',markersize=mark)

        P.xlabel(r'$log(M)$ (M_\odot)')
        P.ylabel(r'$log({\cal M}_{alf})$')
        if(ps):
            P.savefig(directory_out+name+subname+'_mass_Mach_alf_'+str(num).zfill(5)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_mass_Mach_alf_'+str(num).zfill(5)+'.pdf')
        P.savefig(directory_out+name+subname+'_mass_Mach_alf_'+str(num).zfill(5)+'.png')


    nbins=nbins_glob
    mask = prop_clump['sig'] > 0.
    mask_mag = prop_clump['PB'] > 0.
    mask = mask *  mask_mag
    if(np.sum(mask) > 0):
        mass_Mach_alf , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(prop_clump['Mach_alfv'][mask]),bins=nbin,weights=prop_clump['mass'][mask])
        P.clf()
        im = P.imshow(np.log10(mass_Mach_alf.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
        cbar = P.colorbar(im)                                                                                                          
        cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
        P.xlabel(r'$log(M) (M_\odot)$')
        P.ylabel(r'$log({\cal M}_{alf})$')
        if(ps):
            P.savefig(directory_out+name+subname+'_hist2D_mass_Mach_alf'+'_'+str(num)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_hist2D_mass_Mach_alf'+'_'+str(num)+'.pdf')
        P.savefig(directory_out+name+subname+'_hist2D_mass_Mach_alf'+'_'+str(num)+'.png')


    nbins=nbins_glob
    mask = prop_clump['sig'] > 0.
    mask_mag = prop_clump['PB'] > 0.
    mask = mask *  mask_mag

    mask_n = prop_clump['rho_mean'] < 1.e5
    mask = mask*mask_n

    if(np.sum(mask) > 0):
        mass_Mach_alf , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(prop_clump['Mach_alfv'][mask]),bins=nbin,weights=prop_clump['mass'][mask])
        P.clf()
        im = P.imshow(np.log10(mass_Mach_alf.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
        cbar = P.colorbar(im)                                                                                                          
        cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
        P.xlabel(r'$log(M) (M_\odot)$')
        P.ylabel(r'$log({\cal M}_{alf})$')
        if(ps):
            P.savefig(directory_out+name+subname+'_hist2D_mass_Mach_alf_nmean_low'+'_'+str(num)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_hist2D_mass_Mach_alf_nmean_low'+'_'+str(num)+'.pdf')
        P.savefig(directory_out+name+subname+'_hist2D_mass_Mach_alf_nmean_low'+'_'+str(num)+'.png')




    #distribution of ncell as a function of mass
    P.clf()
    mass_vlog = np.log10(prop_clump['mass'])
    P.plot(mass_vlog,np.log10(prop_clump['ncell']),'.',markersize=mark)

    P.xlabel(r'$log(M)$ (M_\odot)')
    P.ylabel(r'$log(Ncell)$')
    if(ps):
        P.savefig(directory_out+name+subname+'_mass_Ncell_'+str(num).zfill(5)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_mass_Ncell_'+str(num).zfill(5)+'.pdf')
    P.savefig(directory_out+name+subname+'_mass_Ncell_'+str(num).zfill(5)+'.png')



    #distribution of n_max as a function of mass
    P.clf()
    mass_vlog = np.log10(prop_clump['mass'])
    P.plot(mass_vlog,np.log10(prop_clump['n_max']),'.',markersize=mark)

    P.xlabel(r'$log(M)$ (M_\odot)')
    P.ylabel(r'$log(n_{max})$ (cm$^{-3}$)')
    if(ps):
        P.savefig(directory_out+name+subname+'_mass_n_max_'+str(num).zfill(5)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_mass_n_max_'+str(num).zfill(5)+'.pdf')
    P.savefig(directory_out+name+subname+'_mass_n_max_'+str(num).zfill(5)+'.png')





    #distribution of mu as a function of mass
    mask = prop_clump['PB'] != 0.
    if(np.sum(mask) >0):
        P.clf()
        mass_vlog = np.log10(prop_clump['mass'])
        P.plot(mass_vlog,np.log10(prop_clump['mu']),'.',markersize=mark)

        P.xlabel(r'$log(M)$ (M_\odot)')
        P.ylabel(r'$log( \mu )$')
        if(ps):
            P.savefig(directory_out+name+subname+'_mass_mu_'+str(num).zfill(5)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_mass_mu_'+str(num).zfill(5)+'.pdf')
        P.savefig(directory_out+name+subname+'_mass_mu_'+str(num).zfill(5)+'.png')


        nbins=nbins_glob
        mask = prop_clump['sig'] > 0.
        mass_mu , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(prop_clump['mu'][mask]),bins=nbin,weights=prop_clump['mass'][mask])
        P.clf()
        im = P.imshow(np.log10(mass_mu.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
        cbar = P.colorbar(im)                                                                                                          
        cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
        P.xlabel(r'$log(M) (M_\odot)$')
        P.ylabel(r'$log( \mu)$')
        if(ps):
            P.savefig(directory_out+name+subname+'_hist2D_mass_mu'+'_'+str(num)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_hist2D_mass_mu'+'_'+str(num)+'.pdf')
        P.savefig(directory_out+name+subname+'_hist2D_mass_mu'+'_'+str(num)+'.png')



        nbins=nbins_glob
        mask = prop_clump['sig'] > 0.
        mask_n = prop_clump['rho_mean'] < 1.e5
        mask = mask*mask_n

        mass_mu , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(prop_clump['mu'][mask]),bins=nbin,weights=prop_clump['mass'][mask])
        P.clf()
        im = P.imshow(np.log10(mass_mu.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
        cbar = P.colorbar(im)                                                                                                          
        cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
        P.xlabel(r'$log(M) (M_\odot)$')
        P.ylabel(r'$log( \mu)$')
        if(ps):
            P.savefig(directory_out+name+subname+'_hist2D_mass_mu_nmean_low'+'_'+str(num)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_hist2D_mass_mu_nmean_low'+'_'+str(num)+'.pdf')
        P.savefig(directory_out+name+subname+'_hist2D_mass_mu_nmean_low'+'_'+str(num)+'.png')






    #distribution of mu as a function of mass for the thermally supercritical cores
    mask = prop_clump['PB'] != 0.
    if(np.sum(mask) >0):
        P.clf()
        mass_vlog = np.log10(prop_clump['mass'])

        nbins=nbins_glob
        mask = prop_clump['sig'] > 0.
        mask_n = np.abs(prop_clump['therm_vir']/prop_clump['grav_vir']) < 1.
        mask_n1 = prop_clump['rho_mean'] < 1.e5
        mask = mask*mask_n*mask_n1

        mass_mu , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(prop_clump['mu'][mask]),bins=nbin)
        P.clf()
        im = P.imshow(np.log10(mass_mu.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
        cbar = P.colorbar(im)                                                                                                          
        cbar.set_label('number of cores')                                          
        P.xlabel(r'$log(M) (M_\odot)$')
        P.ylabel(r'$log( \mu)$')
        if(ps):
            P.savefig(directory_out+name+subname+'_hist2D_virth_mass_mu'+'_'+str(num)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_hist2D_virth_mass_mu'+'_'+str(num)+'.pdf')
        P.savefig(directory_out+name+subname+'_hist2D_virth_mass_mu'+'_'+str(num)+'.png')








    #distribution of n_mean as a function of mass
    mask = prop_clump['PB'] != 0.
    if(np.sum(mask) >0):
        P.clf()



        rad = (prop_clump['size_iner'][:,0]*prop_clump['size_iner'][:,1]*prop_clump['size_iner'][:,2])**(0.3333) #an estimate of the radius
        mask = rad > 0.

        rho_mean = prop_clump['mass'][mask] / rad[mask]**3  
        nbins = nbins_glob 
        mass_nmean , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(rho_mean),bins=nbin) #,weights=prop_clump['mass'][mask]) 
        P.clf()
        im = P.imshow(np.log10(mass_nmean.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
        cbar = P.colorbar(im)                                                                                                          
        cbar.set_label('number of objects')
        P.xlabel(r'$log(M) \, (M_ \odot)$')
        P.ylabel(r'$log(n) \, (cm^{-3})$')                                          
        if(ps):
            P.savefig(directory_out+name+subname+'_hist2D_mass_nmean2'+'_'+str(num)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_hist2D_mass_nmean2'+'_'+str(num)+'.pdf')
        P.savefig(directory_out+name+subname+'_hist2D_mass_nmean2'+'_'+str(num)+'.png')



        nbins = nbins_glob 
        mass_nmean , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(prop_clump['rho_mean'][mask]),bins=nbin) #,weights=prop_clump['mass'][mask]) 
        P.clf()
        im = P.imshow(np.log10(mass_nmean.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
        cbar = P.colorbar(im)                                                                                                          
        cbar.set_label('number of objects')
        P.xlabel(r'$log(M) \, (M_ \odot)$')
        P.ylabel(r'$log(n) \, (cm^{-3})$')                                          
        if(ps):
            P.savefig(directory_out+name+subname+'_hist2D_mass_nmean'+'_'+str(num)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_hist2D_mass_nmean'+'_'+str(num)+'.pdf')
        P.savefig(directory_out+name+subname+'_hist2D_mass_nmean'+'_'+str(num)+'.png')



        nbins = nbins_glob 
        nmax_nmean , xedges, yedges = np.histogram2d(np.log10(prop_clump['n_max'][mask]),np.log10(prop_clump['rho_mean'][mask]),bins=nbin) #,weights=prop_clump['mass'][mask]) 
        P.clf()
        im = P.imshow(np.log10(nmax_nmean.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
        cbar = P.colorbar(im)                                                                                                          
        cbar.set_label('number of objects')
        P.xlabel(r'$log(n_{max}) \, (cm^{-3})$') 
        P.ylabel(r'$log(n) \, (cm^{-3})$')                                          
        if(ps):
            P.savefig(directory_out+name+subname+'_hist2D_nmax_nmean'+'_'+str(num)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_hist2D_nmax_nmean'+'_'+str(num)+'.pdf')
        P.savefig(directory_out+name+subname+'_hist2D_nmax_nmean'+'_'+str(num)+'.png')


        nbins = nbins_glob 
        mass_nmax , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(prop_clump['n_max'][mask]),bins=nbin) #,weights=prop_clump['mass'][mask]) 
        P.clf()
        im = P.imshow(np.log10(mass_nmax.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
        cbar = P.colorbar(im)                                                                                                          
        cbar.set_label('number of objects')
        P.xlabel(r'$log(M) \, (M_ \odot)$')
        P.ylabel(r'$log(n_{max}) \, (cm^{-3})$')                                          
        if(ps):
            P.savefig(directory_out+name+subname+'_hist2D_mass_nmax'+'_'+str(num)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_hist2D_mass_nmax'+'_'+str(num)+'.pdf')
        P.savefig(directory_out+name+subname+'_hist2D_mass_nmax'+'_'+str(num)+'.png')



    #distribution of B as a function of mean density
    P.clf()
    rho_mean_vlog = np.log10(prop_clump['rho_mean'])
    B_v = (np.sqrt(prop_clump['B'][:,0]**2+prop_clump['B'][:,1]**2+prop_clump['B'][:,2]**2))
    P.plot(rho_mean_vlog,np.log10(B_v),'.',markersize=mark)

    P.xlabel(r'$log(\bar{n})$ (cm$^{-3}$)')
    P.ylabel(r'$log( B ) ($\mu$G)')
    if(ps):
        P.savefig(directory_out+name+subname+'_n_B_'+str(num).zfill(5)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_n_B_'+str(num).zfill(5)+'.pdf')
    P.savefig(directory_out+name+subname+'_n_B_'+str(num).zfill(5)+'.png')



    #distribution of B as a function of mean density
    P.clf()
    rho_mean_vlog = np.log10(prop_clump['rho_mean'])
    B_vlog = (np.sqrt(prop_clump['PB']*8.*np.pi))*1.e6 #get the field in muG
    mask = prop_clump['PB'] != 0.
    if(np.sum(mask) >0):
        P.plot(rho_mean_vlog,np.log10(B_vlog),'.',markersize=mark)
        P.xlabel(r'$log(\bar{n})$ (cm$^{-3}$)')
        P.ylabel(r'$log( B ) (\mu$G)')
        if(ps):
            P.savefig(directory_out+name+subname+'_n_PB_'+str(num).zfill(5)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_n_PB_'+str(num).zfill(5)+'.pdf')
        P.savefig(directory_out+name+subname+'_n_PB_'+str(num).zfill(5)+'.png')




    #distribution of alpha = 5/3 sig^2 / (GM/R) as a function of mass with radius based on volume
    P.clf()
    P.plot(mass_vlog,np.log10(prop_clump['alpha_vir']),'.',markersize=mark)
    P.xlabel(r'$log(M)$ (M_\odot)')
    P.ylabel(r'$log( \alpha )$')
    P.ylim(ymin=-0.5,ymax=2.)
    if(ps):
        P.savefig(directory_out+name+subname+'_mass_alpha_'+str(num).zfill(5)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_mass_alpha_'+str(num).zfill(5)+'.pdf')
    P.savefig(directory_out+name+subname+'_mass_alpha_'+str(num).zfill(5)+'.png')


    nbins=nbins_glob
    mask = prop_clump['sig'] > 0.
    rad = (prop_clump['size_iner'][:,0]*prop_clump['size_iner'][:,1]*prop_clump['size_iner'][:,2])**(0.3333) #an estimate of the radius
    mask1 = rad > 0.
    mask = mask * mask1
    mass_alf_vir , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(prop_clump['alpha_vir'][mask]),bins=nbin,weights=prop_clump['mass'][mask])
    P.clf()
    im = P.imshow(np.log10(mass_alf_vir.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(M) (M_\odot)$')
    P.ylabel(r'$log(\alpha)$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_mass_alpha'+'_'+str(num)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_hist2D_mass_alpha'+'_'+str(num)+'.pdf')
    P.savefig(directory_out+name+subname+'_hist2D_mass_alpha'+'_'+str(num)+'.png')




    #distribution of alpha = 5/3 sig^2 / (GM/R) as a function of mass with radius based on inertia momentum
    P.clf()
    P.plot(mass_vlog,np.log10(prop_clump['alpha_vir1']),'.',markersize=mark)
    P.xlabel(r'$log(M)$ (M_\odot)')
    P.ylabel(r'$log( \alpha )$')
    P.ylim(ymin=-0.5,ymax=2.)
    if(ps):
        P.savefig(directory_out+name+subname+'_mass_alpha1_'+str(num).zfill(5)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_mass_alpha1_'+str(num).zfill(5)+'.pdf')
    P.savefig(directory_out+name+subname+'_mass_alpha1_'+str(num).zfill(5)+'.png')

    nbins=nbins_glob
    mask = prop_clump['alpha_vir1'] > 0.
    rad = (prop_clump['size_iner'][:,0]*prop_clump['size_iner'][:,1]*prop_clump['size_iner'][:,2])**(0.3333) #an estimate of the radius
    mask1 = rad > 0.
    mask = mask * mask1
    mass_alf_vir1 , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(prop_clump['alpha_vir1'][mask]),bins=nbin,weights=prop_clump['mass'][mask])
    P.clf()
    im = P.imshow(np.log10(mass_alf_vir1.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(M) (M_\odot)$')
    P.ylabel(r'$log(\alpha)$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_mass_alpha1'+'_'+str(num)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_hist2D_mass_alpha1'+'_'+str(num)+'.pdf')
    P.savefig(directory_out+name+subname+'_hist2D_mass_alpha1'+'_'+str(num)+'.png')



    nbins=nbins_glob
    mask = prop_clump['alpha_vir1'] > 0.
    rad = (prop_clump['size_iner'][:,0]*prop_clump['size_iner'][:,1]*prop_clump['size_iner'][:,2])**(0.3333) #an estimate of the radius
    mask1 = rad > 0.
    mask = mask * mask1

    mask_n = prop_clump['rho_mean'] < 1.e5
    mask = mask * mask_n

    mass_alf_vir1 , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(prop_clump['alpha_vir1'][mask]),bins=nbin,weights=prop_clump['mass'][mask])
    P.clf()
    im = P.imshow(np.log10(mass_alf_vir1.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(M) (M_\odot)$')
    P.ylabel(r'$log(\alpha)$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_mass_alpha1_nmean_low'+'_'+str(num)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_hist2D_mass_alpha1_nmean_low'+'_'+str(num)+'.pdf')
    P.savefig(directory_out+name+subname+'_hist2D_mass_alpha1_nmean_low'+'_'+str(num)+'.png')



    #distribution of alpha = Ekin / Egrav
    if( np.sum(prop_clump['grav_vir'] !=0) > 0):
        P.clf()
        P.plot(mass_vlog,np.log10(np.abs(prop_clump['ener_kin']/prop_clump['grav_vir'])),'.',markersize=mark)

        P.xlabel(r'$log(M)$ (M_\odot)')
        P.ylabel(r'$log( \alpha )$')
        if(ps):
            P.savefig(directory_out+name+subname+'_mass_alpha_vir_'+str(num).zfill(5)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_mass_alpha_vir_'+str(num).zfill(5)+'.pdf')
        P.savefig(directory_out+name+subname+'_mass_alpha_vir_'+str(num).zfill(5)+'.png')


        nbins=nbins_glob
        mask = prop_clump['sig'] > 0.
        mask_vir = np.abs(prop_clump['grav_vir']) > 1.e-10
        mask = mask * mask_vir
        mass_alf_vir2 , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(np.abs(prop_clump['ener_kin'][mask]/prop_clump['grav_vir'][mask])),bins=nbin,weights=prop_clump['mass'][mask])
        P.clf()
        im = P.imshow(np.log10(mass_alf_vir2.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
        cbar = P.colorbar(im)                                                                                                          
        cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
        P.xlabel(r'$log(M) (M_\odot)$')
        P.ylabel(r'$log(\alpha_{vir})$')
        if(ps):
            P.savefig(directory_out+name+subname+'_hist2D_mass_alf_vir'+'_'+str(num)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_hist2D_mass_alf_vir'+'_'+str(num)+'.pdf')
        P.savefig(directory_out+name+subname+'_hist2D_mass_alf_vir'+'_'+str(num)+'.png')




        nbins=nbins_glob
        mask = prop_clump['sig'] > 0.
        mask_vir = np.abs(prop_clump['grav_vir']) > 1.e-10
        mask = mask * mask_vir
        mask_n = prop_clump['rho_mean'] > 1.e5
        mask = mask * mask_n

        mass_alf_vir2 , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(np.abs(prop_clump['ener_kin'][mask]/prop_clump['grav_vir'][mask])),bins=nbin,weights=prop_clump['mass'][mask])
        P.clf()
        im = P.imshow(np.log10(mass_alf_vir2.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
        cbar = P.colorbar(im)                                                                                                          
        cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
        P.xlabel(r'$log(M) (M_\odot)$')
        P.ylabel(r'$log(\alpha_{vir})$')
        if(ps):
            P.savefig(directory_out+name+subname+'_hist2D_mass_alf_vir_nmean'+'_'+str(num)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_hist2D_mass_alf_vir_nmean'+'_'+str(num)+'.pdf')
        P.savefig(directory_out+name+subname+'_hist2D_mass_alf_vir_nmean'+'_'+str(num)+'.png')



        nbins=nbins_glob
        mask = prop_clump['sig'] > 0.
        mask_vir = np.abs(prop_clump['grav_vir']) > 1.e-10
        mask = mask * mask_vir
        mask_n = prop_clump['rho_mean'] < 1.e5
        mask = mask * mask_n

        mass_alf_vir2 , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(np.abs(prop_clump['ener_kin'][mask]/prop_clump['grav_vir'][mask])),bins=nbin,weights=prop_clump['mass'][mask])
        P.clf()
        im = P.imshow(np.log10(mass_alf_vir2.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
        cbar = P.colorbar(im)                                                                                                          
        cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
        P.xlabel(r'$log(M) (M_\odot)$')
        P.ylabel(r'$log(\alpha _{vir})$')
        if(ps):
            P.savefig(directory_out+name+subname+'_hist2D_mass_alf_vir_nmean_low'+'_'+str(num)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_hist2D_mass_alf_vir_nmean_low'+'_'+str(num)+'.pdf')
        P.savefig(directory_out+name+subname+'_hist2D_mass_alf_vir_nmean_low'+'_'+str(num)+'.png')




    #distribution of alpha = Etherm / Egrav
    if( np.sum(prop_clump['grav_vir'] !=0) > 0):
        mask_vir = np.abs(prop_clump['grav_vir']) > 1.e-10
        mask =  mask_vir
        P.clf()
        P.plot(mass_vlog[mask],np.log10(np.abs(prop_clump['therm_vir'][mask]/prop_clump['grav_vir'][mask])),'.',markersize=mark)

        P.xlabel(r'$log(M)$ (M_\odot)')
        P.ylabel(r'$log( \alpha_{th} )$')
        if(ps):
            P.savefig(directory_out+name+subname+'_mass_alpha_thvir_'+str(num).zfill(5)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_mass_alpha_thvir_'+str(num).zfill(5)+'.pdf')
        P.savefig(directory_out+name+subname+'_mass_alpha_thvir_'+str(num).zfill(5)+'.png')

        nbins=nbins_glob
        mask = prop_clump['sig'] > 0.
        mask_vir = np.abs(prop_clump['grav_vir']) > 1.e-10
        mask = mask * mask_vir
        mass_alf_thvir , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(np.abs(prop_clump['therm_vir'][mask]/prop_clump['grav_vir'][mask])),bins=nbin,weights=prop_clump['mass'][mask])
        P.clf()
        im = P.imshow(np.log10(mass_alf_thvir.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
        cbar = P.colorbar(im)                                                                                                          
        cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
        P.xlabel(r'$log(M) (M_\odot)$')
        P.ylabel(r'$log(\alpha _{th})$')
        if(ps):
            P.savefig(directory_out+name+subname+'_hist2D_mass_alf_thvir'+'_'+str(num)+'.ps')
        if(pdf):
            P.savefig(directory_out+name+subname+'_hist2D_mass_alf_thvir'+'_'+str(num)+'.pdf')
        P.savefig(directory_out+name+subname+'_hist2D_mass_alf_thvir'+'_'+str(num)+'.png')




    #distribution of aspect ratio 
    mask = prop_clump['ncell'] > 100
    lamb_tot = (prop_clump['size_iner'][:,0]+prop_clump['size_iner'][:,1]+prop_clump['size_iner'][:,2])/2.
    lamb1 = np.sqrt((lamb_tot[mask] - prop_clump['size_iner'][:,2][mask])/prop_clump['mass'][mask])
    lamb2 = np.sqrt((lamb_tot[mask] - prop_clump['size_iner'][:,1][mask])/prop_clump['mass'][mask])
    lamb3 = np.sqrt((lamb_tot[mask] - prop_clump['size_iner'][:,0][mask])/prop_clump['mass'][mask])

#    P.clf()
#    P.plot(lamb2/lamb3,lamb1/lamb2,'.',markersize=mark)

    a_ratio , xedges, yedges = np.histogram2d(lamb2/lamb3,lamb1/lamb2,bins=nbin)
    P.clf()
    im = P.imshow(a_ratio.T,origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$N$')                                          

    P.xlabel(r'$\mu_2 / \mu_3$')
    P.ylabel(r'$\mu_1 / \mu_2$')
    if(ps):
        P.savefig(directory_out+name+subname+'_hist2D_aratio_'+str(num).zfill(5)+'.ps')
    if(pdf):
        P.savefig(directory_out+name+subname+'_hist2D_aratio_'+str(num).zfill(5)+'.pdf')
    P.savefig(directory_out+name+subname+'_hist2D_aratio_'+str(num).zfill(5)+'.png')



    return


########################################################################################
########################################################################################
## plot the properties of the clumps/cores
########################################################################################
########################################################################################                 
def old_plot_clump_properties(name_prop_clump,name,path_in,num,path_out=None):

    if( path_out is not None):
        directory_out = path_out
    else:
        directory_out = path_in



    f = open(directory_out+name_prop_clump,'r')
    prop_clump = pickle.load(f)
    thres_v = pickle.load(f)
    f.close()

    #markersize
    mark = 3



    #histogram of cell number
    ncell_vlog = np.log10(prop_clump['ncell'])
    log_min=min(ncell_vlog)
    log_max=max(ncell_vlog)
    nbin=50.
    width=(log_max-log_min)/nbin
    P.clf()
    hist_ncell , hist_edges = P.histogram(ncell_vlog,bins=nbin,range=(log_min,log_max))
    P.bar(hist_edges[:-1],hist_ncell,width)
    P.yscale('log', nonposy='clip')
    P.xlabel(r'$log(ncell)$')
    P.ylabel(r'$log(N_c)$')
    P.savefig(directory_out+name+'_hist_ncell_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+'_hist_ncell_'+str(num).zfill(5)+'.jpeg')



    #series of mass spectra (various thresholds)

    #cut on the number of points
    mask = prop_clump['mass'] > 0.
    mask_n = prop_clump['ncell'] > 100
    mask = mask * mask_n
    mass_vlog = np.log10(prop_clump['mass'][mask])
    log_min=min(mass_vlog)
    log_max=max(mass_vlog)
    nbin=50.
    width=(log_max-log_min)/nbin
    P.clf()
    hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
    P.bar(hist_edges[:-1],hist_mass,width)
    P.yscale('log', nonposy='clip')
    P.xlabel(r'$log(M)$ (M_\odot)')
    P.ylabel(r'$log(dN/dlogM)$')
    P.savefig(directory_out+name+'_hist_mass_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+'_hist_mass_'+str(num).zfill(5)+'.jpeg')


    #cut on the mass to flux ratio, big mu
    mask = prop_clump['mass'] > 0.
    mask_n = prop_clump['mu'] > 1.
    mask = mask * mask_n
    mass_vlog = np.log10(prop_clump['mass'][mask])
    log_min=min(mass_vlog)
    log_max=max(mass_vlog)
    nbin=50.
    width=(log_max-log_min)/nbin
    P.clf()
    hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
    P.bar(hist_edges[:-1],hist_mass,width)
    P.yscale('log', nonposy='clip')
    P.xlabel(r'$log(M)$ (M_\odot)')
    P.ylabel(r'$log(dN/dlogM)$')
    P.savefig(directory_out+name+'_hist_mass_mucut_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+'_hist_mass_mucut_'+str(num).zfill(5)+'.jpeg')


    #cut on the mass to flux ratio, small mu
    mask = prop_clump['mass'] > 0.
    mask_n = prop_clump['mu'] < 1.
    mask = mask * mask_n
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass'][mask])
        log_min=min(mass_vlog)
        log_max=max(mass_vlog)
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width)
        P.yscale('log', nonposy='clip')
        P.xlabel(r'$log(M)$ (M_\odot)')
        P.ylabel(r'$log(dN/dlogM)$')
        P.savefig(directory_out+name+'_hist_mass_mulow_'+str(num).zfill(5)+'.ps')
        P.savefig(directory_out+name+'_hist_mass_mulow_'+str(num).zfill(5)+'.jpeg')


    #cut on the mass to flux ratio, mu
    mask = prop_clump['mass'] > 0.
    mask_n = prop_clump['mu'] > 1./3.
    mask = mask * mask_n
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass'][mask])
        log_min=min(mass_vlog)
        log_max=max(mass_vlog)
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width)
        P.yscale('log', nonposy='clip')
        P.xlabel(r'$log(M)$ (M_\odot)')
        P.ylabel(r'$log(dN/dlogM)$')
        P.savefig(directory_out+name+'_hist_mass_mucut_0.3_'+str(num).zfill(5)+'.ps')
        P.savefig(directory_out+name+'_hist_mass_mucut_0.3_'+str(num).zfill(5)+'.jpeg')


    #cut the high virial parameter
    mask = prop_clump['mass'] > 0.
    mask_n = np.abs(prop_clump['ener_kin']/prop_clump['grav_vir']) < 1.
    mask = mask * mask_n
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass'][mask])
        log_min=min(mass_vlog)
        log_max=max(mass_vlog)
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width)
        P.yscale('log', nonposy='clip')
        P.xlabel(r'$log(M)$ (M_\odot)')
        P.ylabel(r'$log(dN/dlogM)$')
        P.savefig(directory_out+name+'_hist_mass_vircut30_'+str(num).zfill(5)+'.ps')
        P.savefig(directory_out+name+'_hist_mass_vircut30_'+str(num).zfill(5)+'.jpeg')


    #cut the low virial parameter
    mask = prop_clump['mass'] > 0.
    mask_n = np.abs(prop_clump['ener_kin']/prop_clump['grav_vir']) > 1.
    mask = mask * mask_n
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass'][mask])
        log_min=min(mass_vlog)
        log_max=max(mass_vlog)
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width)
        P.yscale('log', nonposy='clip')
        P.xlabel(r'$log(M)$ (M_\odot)')
        P.ylabel(r'$log(dN/dlogM)$')
        P.savefig(directory_out+name+'_hist_mass_virhigh30_'+str(num).zfill(5)+'.ps')
        P.savefig(directory_out+name+'_hist_mass_virhigh30_'+str(num).zfill(5)+'.jpeg')


    #cut the high thermal virial parameter
    mask = prop_clump['mass'] > 0.
    mask_n = np.abs(prop_clump['therm_vir']/prop_clump['grav_vir']) < 1.
    mask = mask * mask_n
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass'][mask])
        log_min=min(mass_vlog)
        log_max=max(mass_vlog)
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width)
        P.yscale('log', nonposy='clip')
        P.xlabel(r'$log(M)$ (M_\odot)')
        P.ylabel(r'$log(dN/dlogM)$')
        P.savefig(directory_out+name+'_hist_mass_virthcut_'+str(num).zfill(5)+'.ps')
        P.savefig(directory_out+name+'_hist_mass_virthcut_'+str(num).zfill(5)+'.jpeg')


    #cut the low thermal virial parameter
    mask = prop_clump['mass'] > 0.
    mask_n = np.abs(prop_clump['therm_vir']/prop_clump['grav_vir']) > 1.
    mask = mask * mask_n
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass'][mask])
        log_min=min(mass_vlog)
        log_max=max(mass_vlog)
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width)
        P.yscale('log', nonposy='clip')
        P.xlabel(r'$log(M)$ (M_\odot)')
        P.ylabel(r'$log(dN/dlogM)$')
        P.savefig(directory_out+name+'_hist_mass_virthhigh_'+str(num).zfill(5)+'.ps')
        P.savefig(directory_out+name+'_hist_mass_virthhigh_'+str(num).zfill(5)+'.jpeg')


    mask = prop_clump['mass_0'] > 0.
    mask_n = prop_clump['ncell'] > 100
    mask = mask * mask_n
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass_0'][mask])
        log_min=min(mass_vlog)
        log_max=max(mass_vlog)
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width)
        P.yscale('log', nonposy='clip')
        P.xlabel(r'$log(M)$ (M_\odot)')
        P.ylabel(r'$log(dN/dlogM)$')
        P.savefig(directory_out+name+'_hist_mass_'+str(np.log10(thres_v[0]))+'_'+str(num).zfill(5)+'.ps')
        P.savefig(directory_out+name+'_hist_mass_'+str(np.log10(thres_v[0]))+'_'+str(num).zfill(5)+'.jpeg')


    mask = prop_clump['mass_1'] > 0.
    mask_n = prop_clump['ncell'] > 100
    mask = mask * mask_n
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass_1'][mask])
        log_min=min(mass_vlog)
        log_max=max(mass_vlog)
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width)
        P.yscale('log', nonposy='clip')
        P.xlabel(r'$log(M)$ (M_\odot)')
        P.ylabel(r'$log(dN/dlogM)$')
        P.savefig(directory_out+name+'_hist_mass_'+str(np.log10(thres_v[1]))+'_'+str(num).zfill(5)+'.ps')
        P.savefig(directory_out+name+'_hist_mass_'+str(np.log10(thres_v[1]))+'_'+str(num).zfill(5)+'.jpeg')


    mask = prop_clump['mass_1'] > 0.
    mask_n = prop_clump['ncell'] > 100
    mask = mask * mask_n
    mask1 = prop_clump['n_max'] > thres_v[2]
    mask = mask*mask1
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass_1'][mask])
        log_min=min(mass_vlog)
        log_max=max(mass_vlog)
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width)
        P.yscale('log', nonposy='clip')
        P.xlabel(r'$log(M)$ (M_\odot)')
        P.ylabel(r'$log(dN/dlogM)$')
        P.savefig(directory_out+name+'_hist_mass_'+str(np.log10(thres_v[1]))+'_'+str(np.log10(thres_v[2]))+'_'+str(num).zfill(5)+'.ps')
        P.savefig(directory_out+name+'_hist_mass_'+str(np.log10(thres_v[1]))+'_'+str(np.log10(thres_v[2]))+'_'+str(num).zfill(5)+'.jpeg')


    mask = prop_clump['mass_1'] > 0.
    mask_n = prop_clump['ncell'] > 100
    mask = mask * mask_n
    mask1 = prop_clump['n_max'] > thres_v[4]
    mask = mask*mask1
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass_3'][mask])
        log_min=min(mass_vlog)
        log_max=max(mass_vlog)
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width)
        P.yscale('log', nonposy='clip')
        P.xlabel(r'$log(M)$ (M_\odot)')
        P.ylabel(r'$log(dN/dlogM)$')
        P.savefig(directory_out+name+'_hist_mass_'+str(np.log10(thres_v[1]))+'_'+str(np.log10(thres_v[3]))+'_'+str(num).zfill(5)+'.ps')
        P.savefig(directory_out+name+'_hist_mass_'+str(np.log10(thres_v[1]))+'_'+str(np.log10(thres_v[3]))+'_'+str(num).zfill(5)+'.jpeg')


    mask = prop_clump['mass_1'] > 0.
    mask_n = prop_clump['ncell'] > 100
    mask = mask * mask_n
    mask1 = prop_clump['n_max'] > thres_v[5]
    mask = mask*mask1
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass_5'][mask])
        log_min=min(mass_vlog)
        log_max=max(mass_vlog)
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width)
        P.yscale('log', nonposy='clip')
        P.xlabel(r'$log(M)$ (M_\odot)')
        P.ylabel(r'$log(dN/dlogM)$')


    P.savefig(directory_out+name+'_hist_mass_'+str(np.log10(thres_v[1]))+'_'+str(np.log10(thres_v[5]))+'_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+'_hist_mass_'+str(np.log10(thres_v[1]))+'_'+str(np.log10(thres_v[5]))+'_'+str(num).zfill(5)+'.jpeg')


    #distribution of sigma as a function of size
    P.clf()
    rad = (3./4./np.pi*prop_clump['vol'])**(0.3333) #an estimate of the radius    

    P.plot(np.log10(rad),np.log10(prop_clump['sig']/np.sqrt(3.)),'.',markersize=mark)

    P.xlabel(r'$log(R)$ (pc)')
    P.ylabel(r'$log(\sigma) \, (km \, s^{-1})$')
    P.savefig(directory_out+name+'_R_sig_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+'_R_sig_'+str(num).zfill(5)+'.jpeg')


    #distribution of sigma*col_dens as a function of size
    P.clf()
    col_dens = prop_clump['mass'] / np.pi / rad**2  #column density expressed in Ms / pc^2
    P.plot(np.log10(rad*col_dens),np.log10(prop_clump['sig']/np.sqrt(3.)),'.',markersize=mark)
    P.xlabel(r'$log(R \Sigma)$ (M_\odot pc^{-1})')
    P.ylabel(r'$log(\sigma) \, (km \, s^{-1})$')
    P.savefig(directory_out+name+'_RCD_sig_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+'_RCD_sig_'+str(num).zfill(5)+'.jpeg')





    #distribution of sigma as a function of mass
    P.clf()
    mass_vlog = np.log10(prop_clump['mass'])
    P.plot(mass_vlog,np.log10(prop_clump['sig']/np.sqrt(3.)),'.',markersize=mark)

    P.xlabel(r'$log(M)$ (M_\odot)')
    P.ylabel(r'$log(\sigma) \, (km \, s^{-1})$')
    P.savefig(directory_out+name+'_mass_sig_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+'_mass_sig_'+str(num).zfill(5)+'.jpeg')



    #distribution of volume weighted Mach as a function of mass
    P.clf()
    mass_vlog = np.log10(prop_clump['mass'])
    P.plot(mass_vlog,np.log10(prop_clump['Mach']),'.',markersize=mark)

    P.xlabel(r'$log(M)$ (M_\odot)')
    P.ylabel(r'$log({\cal M})$')
    P.savefig(directory_out+name+'_mass_Mach_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+'_mass_Mach_'+str(num).zfill(5)+'.jpeg')



    #distribution of mass weighted Mach as a function of mass
    P.clf()
    mass_vlog = np.log10(prop_clump['mass'])
    P.plot(mass_vlog,np.log10(prop_clump['Mach_mw']),'.',markersize=mark)

    P.xlabel(r'$log(M)$ (M_\odot)')
    P.ylabel(r'$log({\cal M})$')
    P.savefig(directory_out+name+'_mass_Mach_mw_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+'_mass_Mach_mw_'+str(num).zfill(5)+'.jpeg')



    #distribution of Mach as a function of mass
    P.clf()
    mass_vlog = np.log10(prop_clump['mass'])
    P.plot(mass_vlog,np.log10(prop_clump['Mach_alfv']),'.',markersize=mark)

    P.xlabel(r'$log(M)$ (M_\odot)')
    P.ylabel(r'$log({\cal M}_{alf})$')
    P.savefig(directory_out+name+'_mass_Mach_alf_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+'_mass_Mach_alf_'+str(num).zfill(5)+'.jpeg')


    #distribution of ncell as a function of mass
    P.clf()
    mass_vlog = np.log10(prop_clump['mass'])
    P.plot(mass_vlog,np.log10(prop_clump['ncell']),'.',markersize=mark)

    P.xlabel(r'$log(M)$ (M_\odot)')
    P.ylabel(r'$log(Ncell)$')
    P.savefig(directory_out+name+'_mass_Ncell_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+'_mass_Ncell_'+str(num).zfill(5)+'.jpeg')



    #distribution of n_max as a function of mass
    P.clf()
    mass_vlog = np.log10(prop_clump['mass'])
    P.plot(mass_vlog,np.log10(prop_clump['n_max']),'.',markersize=mark)

    P.xlabel(r'$log(M)$ (M_\odot)')
    P.ylabel(r'$log(n_{max})$ (cm$^{-3}$)')
    P.savefig(directory_out+name+'_mass_n_max_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+'_mass_n_max_'+str(num).zfill(5)+'.jpeg')



    #distribution of mu as a function of mass
    P.clf()
    mass_vlog = np.log10(prop_clump['mass'])
    P.plot(mass_vlog,np.log10(prop_clump['mu']),'.',markersize=mark)

    P.xlabel(r'$log(M)$ (M_\odot)')
    P.ylabel(r'$log( \mu )$')
    P.savefig(directory_out+name+'_mass_mu_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+'_mass_mu_'+str(num).zfill(5)+'.jpeg')


    #distribution of B as a function of mean density
    P.clf()
    rho_mean_vlog = np.log10(prop_clump['rho_mean'])
    B_v = (np.sqrt(prop_clump['B'][:,0]**2+prop_clump['B'][:,1]**2+prop_clump['B'][:,2]**2))
    P.plot(rho_mean_vlog,np.log10(B_v),'.',markersize=mark)

    P.xlabel(r'$log(\bar{n})$ (cm$^{-3}$)')
    P.ylabel(r'$log( B ) ($\mu$G)')
    P.savefig(directory_out+name+'_n_B_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+'_n_B_'+str(num).zfill(5)+'.jpeg')



    #distribution of B as a function of mean density
    P.clf()
    rho_mean_vlog = np.log10(prop_clump['rho_mean'])
    B_vlog = (np.sqrt(prop_clump['PB']*8.*np.pi))*1.e6 #get the field in muG
    P.plot(rho_mean_vlog,np.log10(B_vlog),'.',markersize=mark)
    P.xlabel(r'$log(\bar{n})$ (cm$^{-3}$)')
    P.ylabel(r'$log( B ) ($\mu$G)')
    P.savefig(directory_out+name+'_n_PB_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+'_n_PB_'+str(num).zfill(5)+'.jpeg')




    #distribution of alpha = 5/3 sig^2 / (GM/R) as a function of mass with radius based on volume
    P.clf()
    P.plot(mass_vlog,np.log10(prop_clump['alpha_vir']),'.',markersize=mark)
    P.xlabel(r'$log(M)$ (M_\odot)')
    P.ylabel(r'$log( \alpha )$')
    P.ylim(ymin=-0.5,ymax=2.)
    P.savefig(directory_out+name+'_mass_alpha_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+'_mass_alpha_'+str(num).zfill(5)+'.jpeg')



    #distribution of alpha = 5/3 sig^2 / (GM/R) as a function of mass with radius based on inertia momentum
    P.clf()
    P.plot(mass_vlog,np.log10(prop_clump['alpha_vir1']),'.',markersize=mark)
    P.xlabel(r'$log(M)$ (M_\odot)')
    P.ylabel(r'$log( \alpha )$')
    P.ylim(ymin=-0.5,ymax=2.)
    P.savefig(directory_out+name+'_mass_alpha1_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+'_mass_alpha1_'+str(num).zfill(5)+'.jpeg')



    #distribution of alpha = Ekin / Egrav
    P.clf()
    P.plot(mass_vlog,np.log10(np.abs(prop_clump['ener_kin']/prop_clump['grav_vir'])),'.',markersize=mark)

    P.xlabel(r'$log(M)$ (M_\odot)')
    P.ylabel(r'$log( \alpha )$')
    P.savefig(directory_out+name+'_mass_alpha_vir_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+'_mass_alpha_vir_'+str(num).zfill(5)+'.jpeg')



    #distribution of alpha = Etherm / Egrav
    P.clf()
    P.plot(mass_vlog,np.log10(np.abs(prop_clump['therm_vir']/prop_clump['grav_vir'])),'.',markersize=mark)

    P.xlabel(r'$log(M)$ (M_\odot)')
    P.ylabel(r'$log( \alpha_{th} )$')
    P.savefig(directory_out+name+'_mass_alpha_thvir_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+'_mass_alpha_thvir_'+str(num).zfill(5)+'.jpeg')




    #distribution of aspect ratio 
    P.clf()
    P.plot(prop_clump['size_iner'][:,1]/prop_clump['size_iner'][:,2],prop_clump['size_iner'][:,0]/prop_clump['size_iner'][:,1],'.',markersize=mark)

    P.xlabel(r'$\mu_2 / \mu_3$')
    P.ylabel(r'$\mu_1 / \mu_2$')
    P.savefig(directory_out+name+'_aspect_ratio_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+'_aspect_ratio_'+str(num).zfill(5)+'.jpeg')



    #distribution of aspect ratio depending on mu
    P.clf()
    mask_n = prop_clump['mu'] > 1.
    vec1 = prop_clump['size_iner'][:,1][mask_n]/prop_clump['size_iner'][:,2][mask_n]
    vec2 = prop_clump['size_iner'][:,0][mask_n]/prop_clump['size_iner'][:,1][mask_n]
    P.plot(vec1,vec2,'xr',markersize=mark)

    mask_n = prop_clump['mu'] < 1.
    vec1 = prop_clump['size_iner'][:,1][mask_n]/prop_clump['size_iner'][:,2][mask_n]
    vec2 = prop_clump['size_iner'][:,0][mask_n]/prop_clump['size_iner'][:,1][mask_n]
    P.plot(vec1,vec2,'.k',markersize=mark)

    P.xlabel(r'$\mu_2 / \mu_3$')
    P.ylabel(r'$\mu_1 / \mu_2$')
    P.savefig(directory_out+name+'_aspect_ratio1_'+str(num).zfill(5)+'.ps')
    P.savefig(directory_out+name+'_aspect_ratio1_'+str(num).zfill(5)+'.jpeg')




    return

########################################################################################
########################################################################################
## plot the properties of the clumps/cores belonging to a specific region specified by a 
## position and a radius
########################################################################################
########################################################################################                 
def plot_clump_prop_subregion(name_prop_clump,path_in,num,center,radius,name,path_out=None):

    if( path_out is not None):
        directory_out = path_out
    else:
        directory_out = path_in




    f = open(directory_out+name_prop_clump,'r')
    prop_clump = pickle.load(f)
    thres_v = pickle.load(f)
    f.close()


    #calculate the position with respect to specified center
    xx = prop_clump['pos_n_max'][:,0]-center[0]
    yy = prop_clump['pos_n_max'][:,1]-center[1]
    zz = prop_clump['pos_n_max'][:,2]-center[2]

#    rad_v = np.sqrt(xx**2+yy**2+zz**2)
# use  projection along z axis for now
    rad_v = np.sqrt(xx**2+yy**2)


    #cut on the number of points
    mask = prop_clump['mass'] > 0.
    mask_n = prop_clump['ncell'] > 100
    mask = mask * mask_n

    mask_n = rad_v < radius
    mask = mask * mask_n

    mass_vlog = np.log10(prop_clump['mass'][mask])
    log_min=min(mass_vlog)
    log_max=max(mass_vlog)
    nbin=50.
    width=(log_max-log_min)/nbin
    P.clf()
    hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
    P.bar(hist_edges[:-1],hist_mass,width)
    P.yscale('log', nonposy='clip')
    P.xlabel(r'$log(M)$ (M_\odot)')
    P.ylabel(r'$log(dN/dlogM)$')
    P.savefig(directory_out+'hist_mass_'+str(num).zfill(5)+'_'+name+'.ps')
    P.savefig(directory_out+'hist_mass_'+str(num).zfill(5)+'_'+name+'.jpeg')


    #cut on the mass to flux ratio, big mu
    mask = prop_clump['mass'] > 0.
    mask_n = prop_clump['mu'] > 1.
    mask = mask * mask_n

    mask_n = rad_v < radius
    mask = mask * mask_n


    mass_vlog = np.log10(prop_clump['mass'][mask])
    log_min=min(mass_vlog)
    log_max=max(mass_vlog)
    nbin=50.
    width=(log_max-log_min)/nbin
    P.clf()
    hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
    P.bar(hist_edges[:-1],hist_mass,width)
    P.yscale('log', nonposy='clip')
    P.xlabel(r'$log(M)$ (M_\odot)')
    P.ylabel(r'$log(dN/dlogM)$')
    P.savefig(directory_out+'hist_mass_mucut_'+str(num).zfill(5)+'_'+name+'.ps')
    P.savefig(directory_out+'hist_mass_mucut_'+str(num).zfill(5)+'_'+name+'.jpeg')





    #cut on the mass to flux ratio, small mu
    mask = prop_clump['mass'] > 0.
    mask_n = prop_clump['mu'] < 1.
    mask = mask * mask_n

    mask_n = rad_v < radius
    mask = mask * mask_n


    mass_vlog = np.log10(prop_clump['mass'][mask])
    log_min=min(mass_vlog)
    log_max=max(mass_vlog)
    nbin=50.
    width=(log_max-log_min)/nbin
    P.clf()
    hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
    P.bar(hist_edges[:-1],hist_mass,width)
    P.yscale('log', nonposy='clip')
    P.xlabel(r'$log(M)$ (M_\odot)')
    P.ylabel(r'$log(dN/dlogM)$')
    P.savefig(directory_out+'hist_mass_mulow_'+str(num).zfill(5)+'_'+name+'.ps')
    P.savefig(directory_out+'hist_mass_mulow_'+str(num).zfill(5)+'_'+name+'.jpeg')



    return


##############################################################################################
##############################################################################################
##############################################################################################
#do plots for  COASTDB

       
def plot_clump_properties_COASTDB(name_prop_clump,name,path_in,num,path_out=None,ps=False,nbins_glob=25,mark=3,log_min=None,log_max=None,center=None,radius=None,subname=''):

    if( path_out is not None):
        directory_out = path_out
    else:
        directory_out = path_in


    if(subname != ''):
        subname='_'+subname+'_'


    f = open(directory_out+name_prop_clump,'r')
    prop_clump = pickle.load(f)
    thres_v = pickle.load(f)
    f.close()


    #a list of dictionnaries to describe the plots
    list_plot = []


    #select a subregion for the clumps
    if( radius is not None):
        if(center is None):
            print 'center should be defined as a 3d array'
            return
        if(len(center) != 3 ):
            print 'center should be defined as a 3d array'
            return
        #calculate the position with respect to specified center
        xx = prop_clump['pos_n_max'][:,0]-center[0]
        yy = prop_clump['pos_n_max'][:,1]-center[1]
        zz = prop_clump['pos_n_max'][:,2]-center[2]

        #    rad_v = np.sqrt(xx**2+yy**2+zz**2)
        # use  projection along z axis for now
        rad_v = np.sqrt(xx**2+yy**2)

        mask = rad_v < radius

        #select the clumps of interest
        for key in prop_clump.keys():
            prop_clump[key] = prop_clump[key][mask]




    mask = prop_clump['mass'] > 0.
    mass_vlog = np.log10(prop_clump['mass'][mask])
    if(log_min is None):
        log_min=min(mass_vlog)
    if(log_max is None):
        log_max=max(mass_vlog)



    #histogram of cell number
    ncell_vlog = np.log10(prop_clump['ncell'])
    log_cmin=min(ncell_vlog)
    log_cmax=max(ncell_vlog)
    nbin=50.
    width=(log_cmax-log_cmin)/nbin
    P.clf()
    hist_ncell , hist_edges = P.histogram(ncell_vlog,bins=nbin,range=(log_cmin,log_cmax))
    P.bar(hist_edges[:-1],hist_ncell,width,color='none',edgecolor='k')
    P.yscale('log', nonposy='clip')
    xlabel='$log(ncell)$'
    P.xlabel(r'$log(ncell)$')
    ylabel='$log(N_c)$'
    P.ylabel(r'$log(N_c)$')

    name_im = directory_out+name+subname+'_hist_ncell_'+str(num).zfill(5)
    if(ps):
        P.savefig(name_im+'.ps')
    P.savefig(name_im+'.png')


    dd = dict()
    dd.update({'name':'histogram of the cell number'})
    dd.update({'type':'histogram'})
    dd.update({'description':'Histogram of the cell number.','units-x':'none','units-y':'none','xtitle':xlabel,'ytitle':ylabel})
    data={'png':name_im,'x-axis':hist_edges,'y-axis':hist_ncell}
    dd.update({'data':data})
    list_plot.append(dd)






    #series of mass spectra (various thresholds)

    #cut on the number of points
    mask = prop_clump['mass'] > 0.
    mask_n = prop_clump['ncell'] > 100
    mask = mask * mask_n
    mass_vlog = np.log10(prop_clump['mass'][mask])
    nbin=50.
    width=(log_max-log_min)/nbin
    P.clf()
    hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
    P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='k',linewidth=2)

    P.plot([0.,3.],[1.e4,1.e1],color='k')
    P.yscale('log', nonposy='clip')

    xlabel='$log(M) \, (M_ \odot)$'
    ylabel='$dN/dlogM$'
    P.xlabel(r'$log(M) \, (M_ \odot)$')
    P.ylabel(r'$dN/dlogM$')

    name_im = directory_out+name+subname+'_hist_mass_'+str(num).zfill(5)
    if(ps):
        P.savefig(name_im+'.ps')
    P.savefig(name_im+'.png')



    dd = dict()
    dd.update({'name':'mass histogram'})
    dd.update({'type':'histogram'})
    dd.update({'description':'Histogram of the mass. Structures with at least 100 cells are retained.','units-x':'Msun','units-y':'none','xtitle':xlabel,'ytitle':ylabel})
    data={'png':name_im,'x-axis':hist_edges,'y-axis':hist_mass}
    dd.update({'data':data})
    list_plot.append(dd)




    #cut on the mass to flux ratio, small mu
    mask = prop_clump['mass'] > 0.
    mask_n = prop_clump['mu'] < 1.
    mask_mag = prop_clump['PB'] != 0.
    mask = mask * mask_n * mask_mag
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass'][mask])
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='k',linewidth=2)


        P.plot([0.,3.],[1.e4,1.e1],color='k')
        P.yscale('log', nonposy='clip')


        xlabel='$log(M) \, (M_ \odot)$'
        ylabel='$dN/dlogM$'

        P.xlabel(r'$log(M) \, (M_ \odot)$')
        P.ylabel(r'$dN/dlogM$')


        image_im = directory_out+name+subname+'_hist_mass_mulow_'+str(num).zfill(5)
        if(ps):
            P.savefig(image_im+'.ps')
        P.savefig(image_im+'.png')


        dd = dict()
        dd.update({'name':'mass histogram'})
        dd.update({'type':'histogram'})
        dd.update({'description':'Histogram of the mass. Magnetically subcritical structures are selected, i.e. structures having a mass-to-flux ratio smaller than 1. ','units-x':'Msun','units-y':'none','xtitle':xlabel,'ytitle':ylabel})
        data={'png':name_im,'x-axis':hist_edges,'y-axis':hist_mass}
        dd.update({'data':data})
        list_plot.append(dd)



    #cut on the mass to flux ratio, big mu
    mask = prop_clump['mass'] > 0.
    mask_n = prop_clump['mu'] > 1.
    mask_mag = prop_clump['PB'] != 0.
    mask = mask * mask_n * mask_mag
    mass_vlog = np.log10(prop_clump['mass'][mask])
    nbin=50.
    width=(log_max-log_min)/nbin
    P.clf()
    hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
    P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='k',linewidth=2)
    P.plot([0.,3.],[1.e4,1.e1],color='k')
    P.yscale('log', nonposy='clip')

    xlabel='$log(M) \, (M_ \odot)$'
    ylabel='$dN/dlogM$'

    P.xlabel(r'$log(M) \, (M_ \odot)$')
    P.ylabel(r'$dN/dlogM$')


    image_im = directory_out+name+subname+'_hist_mass_mucut_'+str(num).zfill(5)
    if(ps):
        P.savefig(image_im+'.ps')
    P.savefig(image_im+'.png')



    dd = dict()
    dd.update({'name':'mass histogram'})
    dd.update({'type':'histogram'})
    dd.update({'description':'Histogram of the mass. Magnetically supercritical structures are selected, i.e. structures having a mass-to-flux ratio larger than 1. ','units-x':'Msun','units-y':'none','xtitle':xlabel,'ytitle':ylabel})
    data={'png':name_im,'x-axis':hist_edges,'y-axis':hist_mass}
    dd.update({'data':data})
    list_plot.append(dd)





    #cut the high virial parameter
    mask = prop_clump['mass'] > 0.
    mask_vir = np.abs(prop_clump['grav_vir']) > 1.e-10
    mask = mask * mask_vir
    mask_n = np.abs(prop_clump['ener_kin']/prop_clump['grav_vir']) < 1./2.
    mask = mask * mask_n
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass'][mask])
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='k',linewidth=2)
        P.yscale('log', nonposy='clip')

        xlabel='$log(M) \, (M_ \odot)$'
        ylabel='$dN/dlogM$'

        P.xlabel(r'$log(M) \, (M_ \odot)$')
        P.ylabel(r'$dN/dlogM$')


        image_im = directory_out+name+subname+'_hist_mass_vircut_'+str(num).zfill(5)
        if(ps):
            P.savefig(image_im+'.ps')
        P.savefig(image_im+'.png')


        dd = dict()
        dd.update({'name':'mass histogram'})
        dd.update({'type':'histogram'})
        dd.update({'description':'Histogram of the mass.  Structures such that 2 $E_{kin}/E_{grav}$ < 1 are selected. ','units-x':'Msun','units-y':'none','xtitle':xlabel,'ytitle':ylabel})
        data={'png':name_im,'x-axis':hist_edges,'y-axis':hist_mass}
        dd.update({'data':data})
        list_plot.append(dd)




    #cut the low virial parameter
    mask = prop_clump['mass'] > 0.
    mask_vir = np.abs(prop_clump['grav_vir']) > 1.e-10
    mask = mask * mask_vir
    mask_n = np.abs(prop_clump['ener_kin']/prop_clump['grav_vir']) > 1./2.
    mask = mask * mask_n
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass'][mask])
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='k',linewidth=2)
        P.yscale('log', nonposy='clip')


        xlabel='$log(M) \, (M_ \odot)$'
        ylabel='$dN/dlogM$'

        P.xlabel(r'$log(M) \, (M_ \odot)$')
        P.ylabel(r'$dN/dlogM$')

        image_im = directory_out+name+subname+'_hist_mass_virhigh_'+str(num).zfill(5)
        if(ps):
            P.savefig(image_im+'.ps')
        P.savefig(image_im+'.png')


        dd = dict()
        dd.update({'name':'mass histogram'})
        dd.update({'type':'histogram'})
        dd.update({'description':'Histogram of the mass.  Structures such that 2 $E_{kin}/E_{grav}$ > 1 are selected. ','units-x':'Msun','units-y':'none','xtitle':xlabel,'ytitle':ylabel})
        data={'png':name_im,'x-axis':hist_edges,'y-axis':hist_mass}
        dd.update({'data':data})
        list_plot.append(dd)



    #cut the low thermal virial parameter
    mask = prop_clump['mass'] > 0.
    mask_vir = np.abs(prop_clump['grav_vir']) > 1.e-10
    mask = mask * mask_vir
    mask_n = np.abs(prop_clump['therm_vir']/prop_clump['grav_vir']) > 1./2.
    mask = mask * mask_n
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass'][mask])
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='k',linewidth=2)
        P.yscale('log', nonposy='clip')

        xlabel='$log(M) \, (M_ \odot)$'
        ylabel='$dN/dlogM$'

        P.xlabel(r'$log(M) \, (M_ \odot)$')
        P.ylabel(r'$dN/dlogM$')

        image_im = directory_out+name+subname+'_hist_mass_virthhigh_'+str(num).zfill(5)
        if(ps):
            P.savefig(image_im+'.ps')
        P.savefig(image_im+'.png')


        dd = dict()
        dd.update({'name':'mass histogram'})
        dd.update({'type':'histogram'})
        dd.update({'description':'Histogram of the mass.  Structures such that 2 $E_{th}/E_{grav}$ > 1 are selected. ','units-x':'Msun','units-y':'none','xtitle':xlabel,'ytitle':ylabel})
        data={'png':name_im,'x-axis':hist_edges,'y-axis':hist_mass}
        dd.update({'data':data})
        list_plot.append(dd)




    #cut the high thermal virial parameter
    mask = prop_clump['mass'] > 0.
    mask_vir = np.abs(prop_clump['grav_vir']) > 1.e-10
    mask = mask * mask_vir
    mask_n = np.abs(prop_clump['therm_vir']/prop_clump['grav_vir']) < 1./2.
    mask = mask * mask_n
    if(np.sum(mask) >0):
        mass_vlog = np.log10(prop_clump['mass'][mask])
        nbin=50.
        width=(log_max-log_min)/nbin
        P.clf()
        hist_mass , hist_edges = P.histogram(mass_vlog,bins=nbin,range=(log_min,log_max))
        P.bar(hist_edges[:-1],hist_mass,width,color='none',edgecolor='k',linewidth=2)
        P.plot([0.,3.],[1.e4,1.e1],color='k')
        P.yscale('log', nonposy='clip')

        xlabel='$log(M) \, (M_ \odot)$'
        ylabel='$dN/dlogM$'

        P.xlabel(r'$log(M) \, (M_ \odot)$')
        P.ylabel(r'$dN/dlogM$')

        image_im = directory_out+name+subname+'_hist_mass_virthcut_'+str(num).zfill(5)
        if(ps):
            P.savefig(image_im+'.ps')
        P.savefig(image_im+'.png')


        dd = dict()
        dd.update({'name':'mass histogram'})
        dd.update({'type':'histogram'})
        dd.update({'description':'Histogram of the mass.  Structures such that 2 $E_{th}/E_{grav}$ < 1 are selected. ','units-x':'Msun','units-y':'none','xtitle':xlabel,'ytitle':ylabel})
        data={'png':name_im,'x-axis':hist_edges,'y-axis':hist_mass}
        dd.update({'data':data})
        list_plot.append(dd)



        #take care with mask; mask_n and rho_mean


    #distribution of M as a function of size calculated from inertia matrix
    P.clf()
    rad = (prop_clump['size_iner'][:,0]*prop_clump['size_iner'][:,1]*prop_clump['size_iner'][:,2])**(0.3333) #an estimate of the radius
    mask = rad > 0.
    mask_1 = prop_clump['ncell'] > 10
    mask = mask * mask_1

    nbins = nbins_glob 
    rad_mass , xedges, yedges = np.histogram2d(np.log10(rad[mask]),np.log10(prop_clump['mass'][mask]),bins=nbin,weights=prop_clump['mass'][mask]) 
    P.clf()
    im = P.imshow(np.log10(rad_mass.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          



    label='$log(M) \, (M_\odot)$'
    xlabel='$log(R)$ (pc)'
    ylabel='$log(M) \, (M_ \odot)$'


    P.xlabel(r'$log(R)$ (pc)')
    P.ylabel(r'$log(M) \, (M_ \odot)$')

    name_im = directory_out+name+subname+'_hist2D_rad_mass'+'_'+str(num)
    if(ps):
        P.savefig(name_im+'.ps')
    P.savefig(name_im+'.png')


    dd = dict()
    dd.update({'name':'mass-radius 2D-histogram'})
    dd.update({'type':'bidimentional histogram'})
    dd.update({'description':'Mass weighted bidimentional histogram of the mass vs radius for all selected structures.','units-x':'pc','units-y':'Msun','units':'Msun','xtitle':xlabel,'ytitle':ylabel,'title':label})
    data={'png':name_im,'x-axis':xedges,'y-axis':yedges,'values':rad_mass}
    dd.update({'data':data})
    list_plot.append(dd)





    #distribution of M as a function of velocity dispersion
    nbins=nbins_glob
    mask = prop_clump['sig'] > 0.
    rad_sig , xedges, yedges = np.histogram2d(np.log10(rad[mask]),np.log10(prop_clump['sig'][mask]/np.sqrt(3.)),bins=nbin,weights=prop_clump['mass'][mask]) #,range=[[xmin,xmax],[ymin,ymax]])
    P.clf()
    im = P.imshow(np.log10(rad_sig.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                       


    label='$log(M) \, (M_\odot)$'                                          
    xlabel='$log(R)$ (pc)'
    ylabel='$log(\sigma) \, (km \, s^{-1})$'

                                                                                   
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(R)$ (pc)')
    P.ylabel(r'$log(\sigma) \, (km \, s^{-1})$')


    image_im=directory_out+name+subname+'_hist2D_rad_sig'+'_'+str(num)
    if(ps):
        P.savefig(image_im+'.ps')
    P.savefig(image_im+'.png')


    dd = dict()
    dd.update({'name':'velocity dispersion-radius 2D-histogram'})
    dd.update({'type':'bidimentional histogram'})
    dd.update({'description':'Mass weighted bidimentional histogram of velocity dispersion vs radius for all selected structures.','units-x':'pc','units-y':'km/s','units':'Msun','xtitle':xlabel,'ytitle':ylabel,'title':label})
    data={'png':name_im,'x-axis':xedges,'y-axis':yedges,'values':rad_sig}
    dd.update({'data':data})
    list_plot.append(dd)







    #distribution of sigma*col_dens as a function of size
    col_dens = prop_clump['mass'] / np.pi / rad**2  #column density expressed in Ms / pc^2

    nbins=nbins_glob
    mask = prop_clump['sig'] > 0.
    CDrad_sig , xedges, yedges = np.histogram2d(np.log10(rad[mask]*col_dens[mask]),np.log10(prop_clump['sig'][mask]/np.sqrt(3.)),bins=nbin,weights=prop_clump['mass'][mask])
    P.clf()
    im = P.imshow(np.log10(CDrad_sig.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          

    label='$log(M) \, (M_\odot)$'
    xlabel='$log(R \Sigma) (M_\odot pc^{-1})$'
    ylabel='$log(\sigma) \, (km \, s^{-1})$'

    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(R \Sigma) (M_\odot pc^{-1})$')
    P.ylabel(r'$log(\sigma) \, (km \, s^{-1})$')

    image_im = directory_out+name+subname+'_hist2D_RCD_sig'+'_'+str(num)
    if(ps):
        P.savefig(image_im+'.ps')
    P.savefig(image_im+'.png')


    dd = dict()
    dd.update({'name':'velocity dispersion-column density 2D-histogram'})
    dd.update({'type':'bidimentional histogram'})
    dd.update({'description':'Mass weighted bidimentional histogram of velocity dispersion vs mean column density for all selected structures.','units-x':'Msun pc$^{-2}$','units-y':'km/s','units':'Msun','xtitle':xlabel,'ytitle':ylabel,'title':label})
    data={'png':name_im,'x-axis':xedges,'y-axis':yedges,'values':CDrad_sig}
    dd.update({'data':data})
    list_plot.append(dd)







    #distribution of sigma as a function of mass

    nbins=nbins_glob
    mask = prop_clump['sig'] > 0.
    mass_sig , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(prop_clump['sig'][mask]/np.sqrt(3.)),bins=nbin,weights=prop_clump['mass'][mask])
    P.clf()
    im = P.imshow(np.log10(mass_sig.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                             


    label='$log(M) \, (M_\odot)$'
    xlabel='$log(M) \, (M_\odot)$'
    ylabel='$log(\sigma) \, (km \, s^{-1})$'
                                             
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(M) (M_\odot)$')
    P.ylabel(r'$log(\sigma) \, (km \, s^{-1})$')


    image_im = directory_out+name+subname+'_hist2D_mass_sig'+'_'+str(num)
    if(ps):
        P.savefig(image_im+'.ps')
    P.savefig(image_im+'.png')


    dd = dict()
    dd.update({'name':'velocity dispersion-mass 2D-histogram'})
    dd.update({'type':'bidimentional histogram'})
    dd.update({'description':'Mass weighted bidimentional histogram of velocity dispersion vs mass for all selected structures.','units-x':'Msun$','units-y':'km/s','units':'Msun','xtitle':xlabel,'ytitle':ylabel,'title':label})
    data={'png':name_im,'x-axis':xedges,'y-axis':yedges,'values':mass_sig}
    dd.update({'data':data})
    list_plot.append(dd)




    #distribution of mass weighted Mach as a function of mass

    nbins=nbins_glob
    mask = prop_clump['sig'] > 0.
    mass_mach , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(prop_clump['Mach_mw'][mask]),bins=nbin,weights=prop_clump['mass'][mask])
    P.clf()
    im = P.imshow(np.log10(mass_mach.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)


    label='$log(M) \, (M_\odot)$'
    xlabel='$log(M) (M_\odot)$'
    ylabel='$log({\cal M})$'
                                                                                                          
    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(M) (M_\odot)$')
    P.ylabel(r'$log({\cal M})$')

    image_im=directory_out+name+subname+'_hist2D_mass_Mach_mw'+'_'+str(num)
    if(ps):
        P.savefig(image_im+'.ps')
    P.savefig(image_im+'.png')


    dd = dict()
    dd.update({'name':'Mach number-mass 2D-histogram'})
    dd.update({'type':'bidimentional histogram'})
    dd.update({'description':'Mass weighted bidimentional histogram of Mach number vs mass for all selected structures.','units-x':'Msun$','units-y':'none','units':'Msun','xtitle':xlabel,'ytitle':ylabel,'title':label})
    data={'png':name_im,'x-axis':xedges,'y-axis':yedges,'values':mass_mach}
    dd.update({'data':data})
    list_plot.append(dd)




    #distribution of Alfvenic Mach nmber as a function of mass
    nbins=nbins_glob
    mask = prop_clump['sig'] > 0.
    mask_mag = prop_clump['PB'] > 0.
    mask = mask *  mask_mag
    if(np.sum(mask) > 0):
        mass_Mach_alf , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(prop_clump['Mach_alfv'][mask]),bins=nbin,weights=prop_clump['mass'][mask])
        P.clf()
        im = P.imshow(np.log10(mass_Mach_alf.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
        cbar = P.colorbar(im)                                                                                                          



        label='$log(M) \, (M_\odot)$'
        xlabel='$log(M) (M_\odot)$'
        ylabel='$log({\cal M}_{alf})$'

        cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
        P.xlabel(r'$log(M) (M_\odot)$')
        P.ylabel(r'$log({\cal M}_{alf})$')

        image_im=directory_out+name+subname+'_hist2D_mass_Mach_alf'+'_'+str(num)
        if(ps):
            P.savefig(image_im+'.ps')
        P.savefig(image_im+'.png')


        dd = dict()
        dd.update({'name':'Alfvenic Mach number-mass 2D-histogram'})
        dd.update({'type':'bidimentional histogram'})
        dd.update({'description':'Mass weighted bidimentional histogram of Alfvenic Mach number vs mass for all selected structures.','units-x':'Msun$','units-y':'none','units':'Msun','xtitle':xlabel,'ytitle':ylabel,'title':label})
        data={'png':name_im,'x-axis':xedges,'y-axis':yedges,'values':mass_Mach_alf}
        dd.update({'data':data})
        list_plot.append(dd)




    #distribution of mu as a function of mass
    mask = prop_clump['PB'] != 0.
    if(np.sum(mask) >0):
        nbins=nbins_glob
        mask = prop_clump['sig'] > 0.
        mass_mu , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(prop_clump['mu'][mask]),bins=nbin,weights=prop_clump['mass'][mask])
        P.clf()
        im = P.imshow(np.log10(mass_mu.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
        cbar = P.colorbar(im)                                                                                                          



        label='$log(M) \, (M_\odot)$'
        xlabel='$log(M) (M_\odot)$'
        ylabel='$log(\mu)$'

        cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
        P.xlabel(r'$log(M) (M_\odot)$')
        P.ylabel(r'$log(\mu)$')

        image_im=directory_out+name+subname+'_hist2D_mass_mu'+'_'+str(num)
        if(ps):
            P.savefig(image_im+'.ps')
        P.savefig(image_im+'.png')


        dd = dict()
        dd.update({'name':'mass-to-flux ratio - mass 2D-histogram'})
        dd.update({'type':'bidimentional histogram'})
        dd.update({'description':'Mass weighted bidimentional histogram of mass-to-flux over critical mass-to-flux ratio vs mass for all selected structures.','units-x':'Msun$','units-y':'none','units':'Msun','xtitle':xlabel,'ytitle':ylabel,'title':label})
        data={'png':name_im,'x-axis':xedges,'y-axis':yedges,'values':mass_mu}
        dd.update({'data':data})
        list_plot.append(dd)







    #distribution of alpha = simeq 2 Ekin / Egrav
    nbins=nbins_glob
    mask = prop_clump['alpha_vir1'] > 0.
    rad = (prop_clump['size_iner'][:,0]*prop_clump['size_iner'][:,1]*prop_clump['size_iner'][:,2])**(0.3333) #an estimate of the radius
    mask1 = rad > 0.
    mask = mask * mask1
    mass_alf_vir1 , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(prop_clump['alpha_vir1'][mask]),bins=nbin,weights=prop_clump['mass'][mask])
    P.clf()
    im = P.imshow(np.log10(mass_alf_vir1.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          


    label='$log(M) \, (M_\odot)$'
    xlabel='$log(M) (M_\odot)$'
    ylabel='$log(\alpha)$'

    cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
    P.xlabel(r'$log(M) (M_\odot)$')
    P.ylabel(r'$log(\alpha))')

    image_im=directory_out+name+subname+'_hist2D_mass_alpha'+'_'+str(num)
    if(ps):
        P.savefig(image_im+'.ps')
    P.savefig(image_im+'.png')


    dd = dict()
    dd.update({'name':'alpha - mass 2D-histogram'})
    dd.update({'type':'bidimentional histogram'})
    dd.update({'description':'Mass weighted bidimentional histogram of alpha vs mass for all selected structures.','units-x':'Msun$','units-y':'none','units':'Msun','xtitle':xlabel,'ytitle':ylabel,'title':label})
    data={'png':name_im,'x-axis':xedges,'y-axis':yedges,'values':mass_alf_vir1}
    dd.update({'data':data})
    list_plot.append(dd)






    #distribution of alpha = Ekin / Egrav
    if( np.sum(prop_clump['grav_vir'] !=0) > 0):

        nbins=nbins_glob
        mask = prop_clump['sig'] > 0.
        mask_vir = np.abs(prop_clump['grav_vir']) > 1.e-10
        mask = mask * mask_vir
        mass_alf_vir2 , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(np.abs(prop_clump['ener_kin'][mask]/prop_clump['grav_vir'][mask])),bins=nbin,weights=prop_clump['mass'][mask])
        P.clf()
        im = P.imshow(np.log10(mass_alf_vir2.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
        cbar = P.colorbar(im)                                                                                                          


        label='$log(M) \, (M_\odot)$'
        xlabel='$log(M) (M_\odot)$'
        ylabel='$log(E_{kin} / E_{grav})$'

        cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
        P.xlabel(r'$log(M) (M_\odot)$')
        P.ylabel(r'$log(E_{kin} / E_{grav})$')

        image_im=directory_out+name+subname+'_hist2D_mass_Ekin_Egrav'+'_'+str(num)
        if(ps):
            P.savefig(image_im+'.ps')
        P.savefig(image_im+'.png')


        dd = dict()
        dd.update({'name':'kinetic over gravitational energy ratio - mass 2D-histogram'})
        dd.update({'type':'bidimentional histogram'})
        dd.update({'description':'Mass weighted bidimentional histogram of kinetic over gravitational energy ratio vs mass for all selected structures.','units-x':'Msun$','units-y':'none','units':'Msun','xtitle':xlabel,'ytitle':ylabel,'title':label})
        data={'png':name_im,'x-axis':xedges,'y-axis':yedges,'values':mass_alf_vir2}
        dd.update({'data':data})
        list_plot.append(dd)





    #distribution of alpha = Etherm / Egrav
    if( np.sum(prop_clump['grav_vir'] !=0) > 0):

        nbins=nbins_glob
        mask = prop_clump['sig'] > 0.
        mask_vir = np.abs(prop_clump['grav_vir']) > 1.e-10
        mask = mask * mask_vir
        mass_alf_thvir , xedges, yedges = np.histogram2d(np.log10(prop_clump['mass'][mask]),np.log10(np.abs(prop_clump['therm_vir'][mask]/prop_clump['grav_vir'][mask])),bins=nbin,weights=prop_clump['mass'][mask])
        P.clf()
        im = P.imshow(np.log10(mass_alf_thvir.T),origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
        cbar = P.colorbar(im)                                                                                                          

        label='$log(M) \, (M_\odot)$'
        xlabel='$log(M) (M_\odot)$'
        ylabel='$log(E_{th} / E_{grav})$'

        cbar.set_label(r'$log(M) \, (M_\odot)$')                                          
        P.xlabel(r'$log(M) (M_\odot)$')
        P.ylabel(r'$log(E_{th} / E_{grav})')

        image_im=directory_out+name+subname+'_hist2D_mass_Eth_Egrav'+'_'+str(num)
        if(ps):
            P.savefig(image_im+'.ps')
        P.savefig(image_im+'.png')


        dd = dict()
        dd.update({'name':'thermal over gravitational energy ratio - mass 2D-histogram'})
        dd.update({'type':'bidimentional histogram'})
        dd.update({'description':'Mass weighted bidimentional histogram of thermal over gravitational energy ratio vs mass for all selected structures.','units-x':'Msun$','units-y':'none','units':'Msun','xtitle':xlabel,'ytitle':ylabel,'title':label})
        data={'png':name_im,'x-axis':xedges,'y-axis':yedges,'values':mass_alf_thvir}
        dd.update({'data':data})
        list_plot.append(dd)








    #distribution of aspect ratio 
    mask = prop_clump['ncell'] > 100
    lamb_tot = (prop_clump['size_iner'][:,0]+prop_clump['size_iner'][:,1]+prop_clump['size_iner'][:,2])/2.
    lamb1 = np.sqrt((lamb_tot[mask] - prop_clump['size_iner'][:,2][mask])/prop_clump['mass'][mask])
    lamb2 = np.sqrt((lamb_tot[mask] - prop_clump['size_iner'][:,1][mask])/prop_clump['mass'][mask])
    lamb3 = np.sqrt((lamb_tot[mask] - prop_clump['size_iner'][:,0][mask])/prop_clump['mass'][mask])

#    P.clf()
#    P.plot(lamb2/lamb3,lamb1/lamb2,'.',markersize=mark)

    a_ratio , xedges, yedges = np.histogram2d(lamb2/lamb3,lamb1/lamb2,bins=nbin)
    P.clf()
    im = P.imshow(a_ratio.T,origin='lower',extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],aspect='auto')
    cbar = P.colorbar(im)                                                                                                          
    cbar.set_label(r'$N$')                                          

    P.xlabel(r'$\mu_2 / \mu_3$')
    P.ylabel(r'$\mu_1 / \mu_2$')
#    if(ps):
#        P.savefig(directory_out+name+subname+'_hist2D_aratio_'+str(num).zfill(5)+'.ps')
#    P.savefig(directory_out+name+subname+'_hist2D_aratio_'+str(num).zfill(5)+'.png')



    name_im=directory_out+'/list_fig'+'_'+subname+format(num,'05')+'.save'
    f = open(name_im,'w')
    pickle.dump(list_plot,f)
    f.close()


    return




##############################################################################################
##############################################################################################
##############################################################################################
#do cut of radial and azimuthal velocities
def do_cut(num,path,i_im=0,path_out=None,map_size=512,tag='',ps=False,pdf=False):


    directory = path

    if( path_out is not None):
        directory_out = path_out
    else:
        directory_out = path


    name_npz=directory+'/coldens_z'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.npz'

    npzfile = np.load(name_npz)
    map_size = npzfile['map_size']
    length = npzfile['length']
    map_col = npzfile['map_col']


    name_npz=directory+'/V_z'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.npz'
    npzfile = np.load(name_npz)
    map_size = npzfile['map_size']
    length = npzfile['length']
    map_vx_z = npzfile['map_vx']
    map_vy_z = npzfile['map_vy']


    name_npz=directory+'/V_x'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.npz'
    npzfile = np.load(name_npz)
    map_size = npzfile['map_size']
    length = npzfile['length']
    map_vy_x = npzfile['map_vy']
    map_vz_x = npzfile['map_vz']


    name_npz=directory+'/V_y'+'_'+tag+str(i_im)+'_'+format(num,'05')+'.npz'
    npzfile = np.load(name_npz)
    map_size = npzfile['map_size']
    length = npzfile['length']
    map_vx_y = npzfile['map_vx']
    map_vz_y = npzfile['map_vz']


    vect_pos = (np.arange(map_size)-map_size/2.)/map_size * length

    P.clf()

    P.plot(vect_pos,map_vx_z[:,map_size/2],color='k')
    P.plot(vect_pos,map_vy_z[map_size/2,:],color='r')

    P.xlabel(r'$x \, (pc)$')
    P.ylabel(r'$V_r (km s^{-1})$')

#    P.legend( (r'$M^{-1}$','all structures')  ,loc='upper right')
    P.title=('infall')

    if(ps):
        P.savefig(directory_out+'Vr_'+str(i_im)+'_'+str(num).zfill(5)+'.ps')
    if(pdf):
        P.savefig(directory_out+'Vr_'+str(i_im)+'_'+str(num).zfill(5)+'.pdf')
    P.savefig(directory_out+'Vr_'+str(i_im)+'_'+str(num).zfill(5)+'.png')


    P.clf()

    P.plot(vect_pos,-map_vy_z[:,map_size/2],color='k')
    P.plot(vect_pos,map_vx_z[map_size/2,:],color='r')

    P.xlabel(r'$x \, (pc)$')
    P.ylabel(r'$V_ \phi (km s^{-1})$')

#    P.legend( (r'$M^{-1}$','all structures')  ,loc='upper right')
#    P.title=('all structures')

    P.title=('rotation')

    if(ps):
        P.savefig(directory_out+'Vphi_'+str(i_im)+'_'+str(num).zfill(5)+'.ps')
    if(pdf):
        P.savefig(directory_out+'Vphi_'+str(i_im)+'_'+str(num).zfill(5)+'.pdf')
    P.savefig(directory_out+'Vphi_'+str(i_im)+'_'+str(num).zfill(5)+'.png')



    return
