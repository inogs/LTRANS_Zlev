import numpy as np
from array import array
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import math
try:
  import netCDF4 as nc4
  netcdf_module='netCDF4'
except:
  import scipy.io.netcdf as NC
  netcdf_module='scipy.io.netcdf'
print('using',netcdf_module)

nt=8
rec_per_netcdf_file=4
nk=10
nij=64
dt=86400
Ext0=0
outdir='fake_fields/'

#############################################################
os.system('mkdir '+outdir+' > /dev/null 2>&1')

fig, ax = plt.subplots()
cax = fig.add_axes([0.78, 0.2, 0.03, 0.6])
U=np.zeros((nt,nk,nij,nij),dtype=float)
for t in range(0,nt):
 for k in range(0,nk):
  for j in range(0,nij):
   for i in range(0,nij):
    U[t,nk-1-k,j,i]=(k+1.0)/(9.0*nk) * math.sin( (i+t*2)/nij *2.0*math.pi)* math.sin( (j+t)/nij *2.0*math.pi)

V=np.zeros((nt,nk,nij,nij),dtype=float)
for t in range(0,nt):
 for k in range(0,nk):
  for j in range(0,nij):
   for i in range(0,nij):
    V[t,nk-1-k,j,i]=(k+1.0)/(9.0*nk) * math.cos( (j)/nij *2.0*math.pi) * math.cos( (i+t*2)/nij *2.0*math.pi)

for t in range(0,nt):
 figs, axs = plt.subplots(5,2)
 figs.set_size_inches(15.0, 22.0)
 for k in range(0,nk):
   (ax_i,ax_j)=(k%5,int((k-k%5)/5))
   axs[ax_i,ax_j].cla()
   axs[ax_i,ax_j].set_title('k='+str(k))
   im=axs[ax_i,ax_j].contourf(np.arange(0,nij),np.arange(0,nij),np.sqrt(U[t,k]**2+V[t,k]**2))
   arrow_scale=np.max((np.abs(U[t,k])+np.abs(V[t,k]))/2)/0.03
   print(k,'->>>',(ax_i,ax_j),arrow_scale*0.03)
   axs[ax_i,ax_j].quiver(np.arange(0,nij,2),np.arange(0,nij,2),U[t,k,::2,::2],V[t,k,::2,::2],pivot='middle',angles='xy',scale=arrow_scale,color='w')
   figs.colorbar(im, ax=axs[ax_i,ax_j], orientation='vertical')
 figs.savefig(outdir+'/'+'UV.'+str(int((Ext0+t*dt)/100)).zfill(10)+'.png')
 print(outdir+'/'+'UV.'+str(int((Ext0+t*dt)/100)).zfill(10)+'.png',arrow_scale)

W=np.zeros((nt,nk,nij,nij),dtype=float)
for t in range(0,nt):
 for k in range(0,nk):
  for j in range(0,nij):
   for i in range(0,nij):
    W[t,nk-1-k,j,i]=(k+10.0)/(3.0*nk)*(1-2*(t%2))*((float(t+nij)/nij)*0.5*math.sin(float(nij/2+j)/float(nij)*2*math.pi)+((t+nij)/nij)*0.5*math.cos(float(nij/2+i)/float(nij)*2*math.pi))/20000.

for t in range(0,nt):
 figs, axs = plt.subplots(5,2)
 figs.set_size_inches(15.0, 22.0)
 for k in range(0,nk):
   (ax_i,ax_j)=(k%5,int((k-k%5)/5))
   axs[ax_i,ax_j].cla()
   axs[ax_i,ax_j].set_title('k='+str(k))
   im=axs[ax_i,ax_j].contourf(np.arange(0,nij),np.arange(0,nij),W[t,k])
   arrow_scale=np.max((np.abs(U[t,k])+np.abs(V[t,k]))/2)/0.03
   axs[ax_i,ax_j].quiver(np.arange(0,nij,2),np.arange(0,nij,2),U[t,k,::2,::2],V[t,k,::2,::2],pivot='middle',angles='xy',scale=arrow_scale,color='w')
   figs.colorbar(im, ax=axs[ax_i,ax_j], orientation='vertical')
 figs.savefig(outdir+'/'+'W.'+str(int((Ext0+t*dt)/100)).zfill(10)+'.png')
 print(outdir+'/'+'W.'+str(int((Ext0+t*dt)/100)).zfill(10)+'.png',arrow_scale)

Eta=np.zeros((nt,nij,nij),dtype=float)
for t in range(0,nt):
 for j in range(0,nij):
  for i in range(0,nij):
   Eta[t,j,i]=(1-2*(t%2))*((float(t+2.)/2.)*0.5*math.sin(float(j)/float(nij)*math.pi)+((t+1.)/2.)*0.5*math.cos(float(i)/float(nij)*math.pi))/50.
 ax.cla()
 cax.cla()
 im=ax.contourf(np.arange(0,nij),np.arange(0,nij),Eta[t])
 fig.colorbar(im, cax=cax, orientation='vertical')
 fig.savefig(outdir+'/'+'Eta.'+str(int((Ext0+t*dt)/100)).zfill(10)+'.png')

Uwind=np.zeros((nt,nij,nij),dtype=float)
for t in range(0,nt):
 for j in range(0,nij):
  for i in range(0,nij):
   Uwind[t,j,i]=(1+t%2)*math.sin(float(i)/float(nij)*math.pi)/5.*3.

Vwind=np.zeros((nt,nij,nij),dtype=float)
for t in range(0,nt):
 for j in range(0,nij):
  for i in range(0,nij):
   Vwind[t,j,i]=math.cos(float(j)/float(nij)*2.0*math.pi)*3.

for t in range(0,nt):
 ax.cla()
 cax.cla()
 im=ax.contourf(np.arange(0,nij),np.arange(0,nij),np.sqrt(Uwind[t]**2+Vwind[t]**2))
 arrow_scale=np.max((np.abs(Uwind[t])+np.abs(Vwind[t]))/2)/0.03
 ax.quiver(np.arange(0,nij,2),np.arange(0,nij,2),Uwind[t,::2,::2],Vwind[t,::2,::2],pivot='middle',angles='xy',scale=arrow_scale,color='w')
 fig.colorbar(im, cax=cax, orientation='vertical')
 fig.savefig(outdir+'/'+'EXFwind.'+str(int((Ext0+t*dt)/100)).zfill(10)+'.png')
 
sustr=np.where(Uwind<0,-(abs(Uwind)/20.659)**(1.0/0.4278),(Uwind/20.659)**(1.0/0.4278))
svstr=np.where(Vwind<0,-(abs(Vwind)/20.659)**(1.0/0.4278),(Vwind/20.659)**(1.0/0.4278))

##################################################################################
for t in range(0,nt):
 output_file = open(outdir+'/'+'U.'+str(int((Ext0+t*dt)/100)).zfill(10)+'.data', 'wb')
 float_array = array('d',U[t].flatten())
 float_array.tofile(output_file)
 output_file.close()
for t in range(0,nt):
 output_file = open(outdir+'/'+'V.'+str(int((Ext0+t*dt)/100)).zfill(10)+'.data', 'wb')
 float_array = array('d',V[t].flatten())
 float_array.tofile(output_file)
 output_file.close()
for t in range(0,nt):
 output_file = open(outdir+'/'+'W.'+str(int((Ext0+t*dt)/100)).zfill(10)+'.data', 'wb')
 float_array = array('d',W[t].flatten())
 float_array.tofile(output_file)
 output_file.close()
for t in range(0,nt):
 output_file = open(outdir+'/'+'Eta.'+str(int((Ext0+t*dt)/100)).zfill(10)+'.data', 'wb')
 float_array = array('d',Eta[t].flatten())
 float_array.tofile(output_file)
 output_file.close()
for t in range(0,nt):
 output_file = open(outdir+'/'+'EXFuwind.'+str(int((Ext0+t*dt)/100)).zfill(10)+'.data', 'wb')
 float_array = array('d',Uwind[t].flatten())
 float_array.tofile(output_file)
 output_file.close()
for t in range(0,nt):
 output_file = open(outdir+'/'+'EXFvwind.'+str(int((Ext0+t*dt)/100)).zfill(10)+'.data', 'wb')
 float_array = array('d',Vwind[t].flatten())
 float_array.tofile(output_file)
 output_file.close()
for t in range(0,nt):
 output_file = open(outdir+'/'+'Ustokesdrift.'+str(int((Ext0+t*dt)/100)).zfill(10)+'.data', 'wb')
 float_array = array('d',Ustokes[t].flatten())
 float_array.tofile(output_file)
 output_file.close()
for t in range(0,nt):
 output_file = open(outdir+'/'+'Vstokesdrift.'+str(int((Ext0+t*dt)/100)).zfill(10)+'.data', 'wb')
 float_array = array('d',Vstokes[t].flatten())
 float_array.tofile(output_file)
 output_file.close()

Stokes_in_U=np.zeros((nt,nk,nij,nij),dtype=float)
for k in range(0,nk):
  Stokes_in_U[:,k,:,:]=Ustokes[:,:,:]
for t in range(0,nt):
 output_file = open(outdir+'/'+'Stokes_in_U.'+str(int((Ext0+t*dt)/100)).zfill(10)+'.data', 'wb')
 float_array = array('d',Stokes_in_U[t].flatten())
 float_array.tofile(output_file)
 output_file.close()

Stokes_in_V=np.zeros((nt,nk,nij,nij),dtype=float)
for k in range(0,nk):
  Stokes_in_V[:,k,:,:]=Vstokes[:,:,:]
for t in range(0,nt):
 output_file = open(outdir+'/'+'Stokes_in_V.'+str(int((Ext0+t*dt)/100)).zfill(10)+'.data', 'wb')
 float_array = array('d',Stokes_in_V[t].flatten())
 float_array.tofile(output_file)
 output_file.close()

###################################################################################
my_s_w=np.linspace(-1.0,0.0,nk+1)
my_s_rho=0.5*(my_s_w[1:]+my_s_w[0:-1])
my_Cs_w =np.sin(0.5*math.pi*np.arange(len(my_s_w))/float(len(my_s_w)-1))**2/1.0 -1.0
my_Cs_r =np.sin(0.5*math.pi*np.arange(len(my_s_rho))/float(len(my_s_rho)-1))**2/1.0 -1.0
my_Cs_r=0.5*(my_Cs_w[1:]+my_Cs_w[0:-1])

Wext=np.zeros((nt,nk+1,nij,nij),dtype=float)
Wext[:,0,:,:]=W[:,0,:,:]
Wext[:,1:,:,:]=W[:,0:,:,:]

for tfile in range(0,nt,rec_per_netcdf_file):
  outfile='fake_hydro_field_'+str(int((Ext0+tfile*dt)/100)).zfill(10)+'.nc'
  t0=tfile
  tF=min(nt,tfile+rec_per_netcdf_file)
  print('-------- WRITING FILE '+str(outfile)+'-----------')
  if(netcdf_module=='netCDF4'):
    ncfile= nc4.Dataset(outdir+'/'+outfile,"w")
  else:
    ncfile=NC.netcdf_file(outdir+'/'+outfile,"w")
  ncfile.createDimension("ocean_time", None)
  ncfile.createDimension("xi_rho",   nij)
  ncfile.createDimension("xi_u",     nij-1)
  ncfile.createDimension("xi_v",     nij)
  ncfile.createDimension("eta_rho",  nij)
  ncfile.createDimension("eta_u",    nij)
  ncfile.createDimension("eta_v",    nij-1)
  ncfile.createDimension("s_rho",    len(my_s_rho))
  ncfile.createDimension("s_w",      len(my_s_w))

  time=ncfile.createVariable("ocean_time", 'i', ('ocean_time',))
  setattr(time,'units'              , 'seconds')
  setattr(time,'long_name'          , 'ocean_time')
  setattr(time,'standard_name'      , 'ocean_time')

  Uvel=ncfile.createVariable("u", 'f', ('ocean_time', 's_rho', 'eta_u', 'xi_u',))
  setattr(Uvel,'units'              , 'm/s')
  setattr(Uvel,'long_name'          , 'Current velocity in X direction')
  setattr(Uvel,'standard_name'      , 'u-vel')
  setattr(Uvel,'axis'               , 'X')
  setattr(Uvel,'valid_min'          , np.min(U))
  setattr(Uvel,'valid_max'          , np.max(U))
  Vvel=ncfile.createVariable("v", 'f', ('ocean_time', 's_rho', 'eta_v', 'xi_v',))
  setattr(Vvel,'units'              , 'm/s')
  setattr(Vvel,'long_name'          , 'Current velocity in Y direction')
  setattr(Vvel,'standard_name'      , 'v-vel')
  setattr(Vvel,'axis'               , 'Y')
  setattr(Vvel,'valid_min'          , np.min(V))
  setattr(Vvel,'valid_max'          , np.max(V))
  Wvel=ncfile.createVariable("w", 'f', ('ocean_time', 's_w', 'eta_rho', 'xi_rho',))
  setattr(Wvel,'units'              , 'm/s')
  setattr(Wvel,'long_name'          , '')
  setattr(Wvel,'standard_name'      , '')
  setattr(Wvel,'axis'               , 'Z')
  setattr(Wvel,'valid_min'          , np.min(W))
  setattr(Wvel,'valid_max'          , np.max(W))
  Zeta=ncfile.createVariable("zeta", 'f', ('ocean_time', 'eta_rho', 'xi_rho',))
  setattr(Zeta,'units'              , 'meters')
  setattr(Zeta,'long_name'          , '')
  setattr(Zeta,'standard_name'      , '')
  setattr(Zeta,'axis'               , 'Z')
  setattr(Zeta,'valid_min'          , np.min(Eta))
  setattr(Zeta,'valid_max'          , np.max(Eta))
  Uwnd=ncfile.createVariable("Uwind", 'f', ('ocean_time', 'eta_u', 'xi_u',))
  setattr(Uwnd,'units'              , 'm/s')
  setattr(Uwnd,'long_name'          , 'surface U-wind component')
  setattr(Uwnd,'standard_name'      , 'Uwind')
  setattr(Uwnd,'axis'               , 'X')
  setattr(Uwnd,'valid_min'          , np.min(Uwind))
  setattr(Uwnd,'valid_max'          , np.max(Uwind))
  Vwnd=ncfile.createVariable("Vwind", 'f', ('ocean_time', 'eta_v', 'xi_v',))
  setattr(Vwnd,'units'              , 'm/s')
  setattr(Vwnd,'long_name'          , 'surface V-wind component')
  setattr(Vwnd,'standard_name'      , 'Vwind')
  setattr(Vwnd,'axis'               , 'Y')
  setattr(Vwnd,'valid_min'          , np.min(Vwind))
  setattr(Vwnd,'valid_max'          , np.max(Vwind))
  SUstr=ncfile.createVariable("sustr", 'f', ('ocean_time', 'eta_u', 'xi_u',))
  setattr(SUstr,'units'              , 'N/m2')
  setattr(SUstr,'long_name'          , 'surface U-momentum wind-stress')
  setattr(SUstr,'standard_name'      , 'sustr')
  setattr(SUstr,'axis'               , 'X')
  setattr(SUstr,'valid_min'          , np.min(sustr))
  setattr(SUstr,'valid_max'          , np.max(sustr))
  SVstr=ncfile.createVariable("svstr", 'f', ('ocean_time', 'eta_v', 'xi_v',))
  setattr(SVstr,'units'              , 'N/m2')
  setattr(SVstr,'long_name'          , 'surface V-momentum wind-stress')
  setattr(SVstr,'standard_name'      , 'svstr')
  setattr(SVstr,'axis'               , 'Y')
  setattr(SVstr,'valid_min'          , np.min(svstr))
  setattr(SVstr,'valid_max'          , np.max(svstr))
  Ustk=ncfile.createVariable("ustokes", 'f', ('ocean_time', 'eta_u', 'xi_u',))
  setattr(Ustk,'long_name'          , 'surface stokes-drift component in x-direction')
  setattr(Ustk,'standard_name'      , 'ustokes')
  setattr(Ustk,'axis'               , 'X')
  setattr(Ustk,'valid_min'          , np.min(Ustokes))
  setattr(Ustk,'valid_max'          , np.max(Ustokes))
  Vstk=ncfile.createVariable("vstokes", 'f', ('ocean_time', 'eta_v', 'xi_v',))
  setattr(Vstk,'long_name'          , 'surface stokes-drift component in y-direction')
  setattr(Vstk,'standard_name'      , 'vstokes')
  setattr(Vstk,'axis'               , 'Y')
  setattr(Vstk,'valid_min'          , np.min(Vstokes))
  setattr(Vstk,'valid_max'          , np.max(Vstokes))
  s_rho=ncfile.createVariable("s_rho", 'f', ('s_rho',))
  setattr(s_rho,'long_name'          , 'S-coordinate at RHO-points')
  setattr(s_rho,'standard_name'      , 's_rho')
  setattr(s_rho,'valid_min'          , np.min(my_s_rho))
  setattr(s_rho,'valid_max'          , np.max(my_s_rho))
  s_rho[:]=my_s_rho[:]
  s_w=ncfile.createVariable("s_w", 'f', ('s_w',))
  setattr(s_w,'long_name'          , 'S-coordinate at W-points')
  setattr(s_w,'standard_name'      , 's_w')
  setattr(s_w,'valid_min'          , np.min(my_s_w))
  setattr(s_w,'valid_max'          , np.max(my_s_w))
  s_w[:]=my_s_w[:]
  Cs_rho=ncfile.createVariable("Cs_r", 'f', ('s_rho',))
  setattr(Cs_rho,'long_name'          , 'S-coordinate stretching curves at RHO-points')
  setattr(Cs_rho,'standard_name'      , 'Cs_r')
  setattr(Cs_rho,'valid_min'          , np.min(my_Cs_r))
  setattr(Cs_rho,'valid_max'          , np.max(my_Cs_r))
  Cs_rho[:]=my_Cs_r[:]
  Cs_w=ncfile.createVariable("Cs_w", 'f', ('s_w',))
  setattr(Cs_w,'long_name'          , 'S-coordinate stretching curves at W-points')
  setattr(Cs_w,'standard_name'      , 'Cs_w')
  setattr(Cs_w,'valid_min'          , np.min(my_Cs_w))
  setattr(Cs_w,'valid_max'          , np.max(my_Cs_w))
  Cs_w[:]=my_Cs_w[:]
  for i,t in enumerate(range(t0,tF)):
    time[i]=Ext0+t*dt
    # add ::-1 to invert k-order of the 3d fields as MITgcm has an inverted k axis respect to ROMs
    Uvel[i,:,:,:]=U[t,::-1,:,1:]
    Vvel[i,:,:,:]=V[t,::-1,1:,:]
    Wvel[i,:,:,:]=Wext[t,::-1,:,:]
    Zeta[i,:,:]=Eta[t,:,:]
    Uwnd[i,:,:]=Uwind[t,:,1:]
    Vwnd[i,:,:]=Vwind[t,1:,:]
    SUstr[i,:,:]=sustr[t,:,1:]
    SVstr[i,:,:]=svstr[t,1:,:]
    Ustk[i,:,:]=Ustokes[t,:,1:]
    Vstk[i,:,:]=Vstokes[t,1:,:]

  ncfile.close()


