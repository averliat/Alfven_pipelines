import matplotlib.lines as mlines
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg


d=235 #pc


#############################   C18O PdBI   #########################

data=np.loadtxt('C18O/PdBI-selfcal/4A1/data-PVIpix-C18O-PdBI-2pixbeam.txt')
offRA=data[:,0] #arcsec
offDEC=data[:,1]  #arcsec
velo=data[:,2] #km/s
err_velo=data[:,3] #km/s
int_peak=data[:,5]
rms=data[:,6]

offRA_PdBI=[]
offDEC_PdBI=[]
velo_PdBI=[]
err_velo_PdBI=[]
for i in range(0,len(offDEC)):
    if abs(offRA[i])<=5.0 and abs(offDEC[i])<=4.0 and int_peak[i]>=5*rms[i] and abs(err_velo[i])<0.5 and velo[i]>4 and velo[i]<9:
        offRA_PdBI.append(offRA[i]/3600/180*np.pi)
        offDEC_PdBI.append(offDEC[i]/3600/180*np.pi)
        velo_PdBI.append(velo[i])
        err_velo_PdBI.append(err_velo[i])

A = np.column_stack((np.ones(len(offRA_PdBI)), offRA_PdBI, offDEC_PdBI))
coeffs, residuals, rank, sig_val = np.linalg.lstsq(A,velo_PdBI)
coeffs_err, residuals_err, rank_err, sig_val_err = np.linalg.lstsq(A,err_velo_PdBI)
#print coeffs
#print coeffs_err


r2 = 1 - residuals / sum((velo_PdBI-np.mean(velo_PdBI))**2)
#print r2

r2_err = 1 - residuals_err / sum((err_velo_PdBI-np.mean(err_velo_PdBI))**2)
#print r2_err

sys_vel=coeffs[0] #km/s
err_sys_vel=coeffs_err[0] #km/s
coeff_x=coeffs[1] #km/s/rad
err_coeff_x=coeffs_err[1] #km/s/rad
coeff_y=coeffs[2] #km/s/rad
err_coeff_y=coeffs_err[2] #km/s/rad

magn=np.sqrt(coeff_x**2+coeff_y**2)/(180/np.pi*3600*d*4.84814e-6) #km/s/pc
err_magn=((1./np.sqrt(coeff_x**2+coeff_y**2))*np.sqrt(coeff_x**2*err_coeff_x**2+coeff_y**2*err_coeff_y**2))/(180/np.pi*3600*d*4.84814e-6) #km/s/pc
PA=np.arctan(coeff_x/coeff_y)*180/np.pi #deg
err_PA=((1./(1.+(coeff_x/coeff_y)**2))*(1./coeff_y)*np.sqrt(err_coeff_x**2+(coeff_x/coeff_y)**2*err_coeff_y**2))*180/np.pi #deg

print "##   C18O PdBI   ##"
print "vitesse syst = %f +/- %f " %(sys_vel, err_sys_vel)
print "magn = %f +/- %f " %(magn, err_magn)
print "PA = %f +/- %f " %(PA, err_PA)
print "-------"
                                
                                

#############################   C18O PdBI+30m   #########################

data=np.loadtxt('C18O/merging-selfcal/weight10/4A1/data-PVIpix-C18O-merged-2pixbeam.txt')
offRA=data[:,0] #arcsec
offDEC=data[:,1]  #arcsec
velo=data[:,2] #km/s
err_velo=data[:,3] #km/s
int_peak=data[:,5]
rms=data[:,6]

offRA_merged=[]
offDEC_merged=[]
velo_merged=[]
err_velo_merged=[]
for i in range(0,len(offDEC)):
    if abs(offRA[i])<=4.0 and abs(offDEC[i])<=4.0 and int_peak[i]>=5*rms[i] and abs(err_velo[i])<0.5 and velo[i]>4 and velo[i]<9:
        offRA_merged.append(offRA[i]/3600/180*np.pi)
        offDEC_merged.append(offDEC[i]/3600/180*np.pi)
        velo_merged.append(velo[i])
        err_velo_merged.append(err_velo[i])



A = np.column_stack((np.ones(len(offRA_merged)), offRA_merged, offDEC_merged))
coeffs, residuals, rank, sig_val = np.linalg.lstsq(A,velo_merged)
coeffs_err, residuals_err, rank_err, sig_val_err = np.linalg.lstsq(A,err_velo_merged)
#print coeffs
#print coeffs_err


r2 = 1 - residuals / sum((velo_merged-np.mean(velo_merged))**2)
#print r2

r2_err = 1 - residuals_err / sum((err_velo_merged-np.mean(err_velo_merged))**2)
#print r2_err

sys_vel=coeffs[0] #km/s
err_sys_vel=coeffs_err[0] #km/s
coeff_x=coeffs[1] #km/s/rad
err_coeff_x=coeffs_err[1] #km/s/rad
coeff_y=coeffs[2] #km/s/rad
err_coeff_y=coeffs_err[2] #km/s/rad

magn=np.sqrt(coeff_x**2+coeff_y**2)/(180/np.pi*3600*d*4.84814e-6) #km/s/pc
err_magn=((1./np.sqrt(coeff_x**2+coeff_y**2))*np.sqrt(coeff_x**2*err_coeff_x**2+coeff_y**2*err_coeff_y**2))/(180/np.pi*3600*d*4.84814e-6) #km/s/pc
PA=np.arctan(coeff_x/coeff_y)*180/np.pi #deg
err_PA=((1./(1.+(coeff_x/coeff_y)**2))*(1./coeff_y)*np.sqrt(err_coeff_x**2+(coeff_x/coeff_y)**2*err_coeff_y**2))*180/np.pi #deg


print "##   C18O PdBI+30m   ##"
print "vitesse syst = %f +/- %f " %(sys_vel, err_sys_vel)
print "magn = %f +/- %f " %(magn, err_magn)
print "PA = %f +/- %f " %(PA, err_PA)
print "-------"



#############################   C18O 30m   #########################

data=np.loadtxt('C18O/30m/4A1/data-PVIpix-C18O-30m-2pixbeam.txt')
offRA=data[:,0] #arcsec
offDEC=data[:,1]  #arcsec
velo=data[:,2] #km/s
err_velo=data[:,3] #km/s
int_peak=data[:,5]
rms=data[:,6]

offRA_30m=[]
offDEC_30m=[]
velo_30m=[]
err_velo_30m=[]
for i in range(0,len(offDEC)):
    if abs(offRA[i])<=40.0 and abs(offDEC[i])<=40.0 and int_peak[i]>=5*rms[i] and abs(err_velo[i])<0.5 and velo[i]>5 and velo[i]<10:
        offRA_30m.append(offRA[i]/3600/180*np.pi)
        offDEC_30m.append(offDEC[i]/3600/180*np.pi)
        velo_30m.append(velo[i])
        err_velo_30m.append(err_velo[i])



A = np.column_stack((np.ones(len(offRA_30m)), offRA_30m, offDEC_30m))
coeffs, residuals, rank, sig_val = np.linalg.lstsq(A,velo_30m)
coeffs_err, residuals_err, rank_err, sig_val_err = np.linalg.lstsq(A,err_velo_30m)
#print coeffs
#print coeffs_err


r2 = 1 - residuals / sum((velo_30m-np.mean(velo_30m))**2)
#print r2

r2_err = 1 - residuals_err / sum((err_velo_30m-np.mean(err_velo_30m))**2)
#print r2_err

sys_vel=coeffs[0] #km/s
err_sys_vel=coeffs_err[0] #km/s
coeff_x=coeffs[1] #km/s/rad
err_coeff_x=coeffs_err[1] #km/s/rad
coeff_y=coeffs[2] #km/s/rad
err_coeff_y=coeffs_err[2] #km/s/rad

magn=np.sqrt(coeff_x**2+coeff_y**2)/(180/np.pi*3600*d*4.84814e-6) #km/s/pc
err_magn=((1./np.sqrt(coeff_x**2+coeff_y**2))*np.sqrt(coeff_x**2*err_coeff_x**2+coeff_y**2*err_coeff_y**2))/(180/np.pi*3600*d*4.84814e-6) #km/s/pc
PA=np.arctan(coeff_x/coeff_y)*180/np.pi #deg
err_PA=((1./(1.+(coeff_x/coeff_y)**2))*(1./coeff_y)*np.sqrt(err_coeff_x**2+(coeff_x/coeff_y)**2*err_coeff_y**2))*180/np.pi #deg

print "##   C18O 30m   ##"
print "vitesse syst = %f +/- %f " %(sys_vel, err_sys_vel)
print "magn = %f +/- %f " %(magn, err_magn)
print "PA = %f +/- %f " %(PA, err_PA)
print "-------"





#############################   N2H+ PdBI   #########################

data=np.loadtxt('N2H+/PdBI/4A1/data-PVIpix-N2H+-PdBI-2pixbeam.txt')
offRA=data[:,0] #arcsec
offDEC=data[:,1]  #arcsec
velo=data[:,4] #km/s
err_velo=data[:,5] #km/s
int_peak=data[:,11]
rms=data[:,12]

offRA_PdBI=[]
offDEC_PdBI=[]
velo_PdBI=[]
err_velo_PdBI=[]
for i in range(0,len(offDEC)):
    if abs(offRA[i])<=10.0 and abs(offDEC[i])<=20.0 and int_peak[i]>=5*rms[i] and abs(err_velo[i])<0.5 and velo[i]>4 and velo[i]<9:
        offRA_PdBI.append(offRA[i]/3600/180*np.pi)
        offDEC_PdBI.append(offDEC[i]/3600/180*np.pi)
        velo_PdBI.append(velo[i])
        err_velo_PdBI.append(err_velo[i])

A = np.column_stack((np.ones(len(offRA_PdBI)), offRA_PdBI, offDEC_PdBI))
coeffs, residuals, rank, sig_val = np.linalg.lstsq(A,velo_PdBI)
coeffs_err, residuals_err, rank_err, sig_val_err = np.linalg.lstsq(A,err_velo_PdBI)
#print coeffs
#print coeffs_err


r2 = 1 - residuals / sum((velo_PdBI-np.mean(velo_PdBI))**2)
#print r2

r2_err = 1 - residuals_err / sum((err_velo_PdBI-np.mean(err_velo_PdBI))**2)
#print r2_err

sys_vel=coeffs[0] #km/s
err_sys_vel=coeffs_err[0] #km/s
coeff_x=coeffs[1] #km/s/rad
err_coeff_x=coeffs_err[1] #km/s/rad
coeff_y=coeffs[2] #km/s/rad
err_coeff_y=coeffs_err[2] #km/s/rad

magn=np.sqrt(coeff_x**2+coeff_y**2)/(180/np.pi*3600*d*4.84814e-6) #km/s/pc
err_magn=((1./np.sqrt(coeff_x**2+coeff_y**2))*np.sqrt(coeff_x**2*err_coeff_x**2+coeff_y**2*err_coeff_y**2))/(180/np.pi*3600*d*4.84814e-6) #km/s/pc
PA=np.arctan(coeff_x/coeff_y)*180/np.pi #deg
err_PA=((1./(1.+(coeff_x/coeff_y)**2))*(1./coeff_y)*np.sqrt(err_coeff_x**2+(coeff_x/coeff_y)**2*err_coeff_y**2))*180/np.pi #deg


print "##   N2H+ PdBI   ##"
print "vitesse syst = %f +/- %f " %(sys_vel, err_sys_vel)
print "magn = %f +/- %f " %(magn, err_magn)
print "PA = %f +/- %f " %(PA, err_PA)
print "-------"


#############################   N2H+ PdBI+30m   #########################

data=np.loadtxt('N2H+/merging/weight5/4A1/data-PVIpix-N2H+-merging-2pixbeam.txt')
offRA=data[:,0] #arcsec
offDEC=data[:,1]  #arcsec
velo=data[:,4] #km/s
err_velo=data[:,5] #km/s
int_peak=data[:,11]
rms=data[:,12]

offRA_merged=[]
offDEC_merged=[]
velo_merged=[]
err_velo_merged=[]
for i in range(0,len(offDEC)):
    if abs(offRA[i])<=20.0 and abs(offDEC[i])<=20.0 and int_peak[i]>=5*rms[i] and abs(err_velo[i])<0.5 and velo[i]>4 and velo[i]<9:
        offRA_merged.append(offRA[i]/3600/180*np.pi)
        offDEC_merged.append(offDEC[i]/3600/180*np.pi)
        velo_merged.append(velo[i])
        err_velo_merged.append(err_velo[i])



A = np.column_stack((np.ones(len(offRA_merged)), offRA_merged, offDEC_merged))
coeffs, residuals, rank, sig_val = np.linalg.lstsq(A,velo_merged)
coeffs_err, residuals_err, rank_err, sig_val_err = np.linalg.lstsq(A,err_velo_merged)
#print coeffs
#print coeffs_err


r2 = 1 - residuals / sum((velo_merged-np.mean(velo_merged))**2)
#print r2

r2_err = 1 - residuals_err / sum((err_velo_merged-np.mean(err_velo_merged))**2)
#print r2_err

sys_vel=coeffs[0] #km/s
err_sys_vel=coeffs_err[0] #km/s
coeff_x=coeffs[1] #km/s/rad
err_coeff_x=coeffs_err[1] #km/s/rad
coeff_y=coeffs[2] #km/s/rad
err_coeff_y=coeffs_err[2] #km/s/rad


magn=np.sqrt(coeff_x**2+coeff_y**2)/(180/np.pi*3600*d*4.84814e-6) #km/s/pc
err_magn=((1./np.sqrt(coeff_x**2+coeff_y**2))*np.sqrt(coeff_x**2*err_coeff_x**2+coeff_y**2*err_coeff_y**2))/(180/np.pi*3600*d*4.84814e-6) #km/s/pc
PA=np.arctan(coeff_x/coeff_y)*180/np.pi #deg
err_PA=((1./(1.+(coeff_x/coeff_y)**2))*(1./coeff_y)*np.sqrt(err_coeff_x**2+(coeff_x/coeff_y)**2*err_coeff_y**2))*180/np.pi #deg


print "##   N2H+ PdBI+30m   ##"
print "vitesse syst = %f +/- %f " %(sys_vel, err_sys_vel)
print "magn = %f +/- %f " %(magn, err_magn)
print "PA = %f +/- %f " %(PA, err_PA)
print "-------"




#############################   N2H+ 30m   #########################

data=np.loadtxt('N2H+/30m/4A1/data-PVIpix-N2H+-30m-2pixbeam.txt')
offRA=data[:,0] #arcsec
offDEC=data[:,1]  #arcsec
velo=data[:,4] #km/s
err_velo=data[:,5] #km/s
int_peak=data[:,11]
rms=data[:,12]


offRA_30m=[]
offDEC_30m=[]
velo_30m=[]
err_velo_30m=[]
for i in range(0,len(offDEC)):
    if abs(offRA[i])<=80.0 and abs(offDEC[i])<=80.0 and int_peak[i]>=5*rms[i] and abs(err_velo[i])<0.5 and velo[i]>5 and velo[i]<10:
        offRA_30m.append(offRA[i]/3600/180*np.pi)
        offDEC_30m.append(offDEC[i]/3600/180*np.pi)
        velo_30m.append(velo[i])
        err_velo_30m.append(err_velo[i])



A = np.column_stack((np.ones(len(offRA_30m)), offRA_30m, offDEC_30m))
coeffs, residuals, rank, sig_val = np.linalg.lstsq(A,velo_30m)
coeffs_err, residuals_err, rank_err, sig_val_err = np.linalg.lstsq(A,err_velo_30m)
#print coeffs
#print coeffs_err


r2 = 1 - residuals / sum((velo_30m-np.mean(velo_30m))**2)
#print r2

r2_err = 1 - residuals_err / sum((err_velo_30m-np.mean(err_velo_30m))**2)
#print r2_err

sys_vel=coeffs[0] #km/s
err_sys_vel=coeffs_err[0] #km/s
coeff_x=coeffs[1] #km/s/rad
err_coeff_x=coeffs_err[1] #km/s/rad
coeff_y=coeffs[2] #km/s/rad
err_coeff_y=coeffs_err[2] #km/s/rad

magn=np.sqrt(coeff_x**2+coeff_y**2)/(180/np.pi*3600*d*4.84814e-6) #km/s/pc
err_magn=((1./np.sqrt(coeff_x**2+coeff_y**2))*np.sqrt(coeff_x**2*err_coeff_x**2+coeff_y**2*err_coeff_y**2))/(180/np.pi*3600*d*4.84814e-6) #km/s/pc
PA=np.arctan(coeff_x/coeff_y)*180/np.pi #deg
err_PA=((1./(1.+(coeff_x/coeff_y)**2))*(1./coeff_y)*np.sqrt(err_coeff_x**2+(coeff_x/coeff_y)**2*err_coeff_y**2))*180/np.pi #deg


print "##   N2H+ 30m   ##"
print "vitesse syst = %f +/- %f " %(sys_vel, err_sys_vel)
print "magn = %f +/- %f " %(magn, err_magn)
print "PA = %f +/- %f " %(PA, err_PA)
print "-------"


exit()
