import numpy as np
import matplotlib.pyplot as plt
try:
  import netCDF4 as nc4
  netcdf_module='netCDF4'
except:
  import scipy.io as NC
  netcdf_module='scipy.io'

print('using',netcdf_module)
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math
import os
outdir='input/'
os.system('mkdir '+outdir+' > /dev/null 2>&1')
nij=64
nk=10

lon_rho,lat_rho=np.meshgrid( np.linspace(41,41.5,nij) , np.linspace(17,17.5,nij) )
#print('lon_rho:',np.shape(lon_rho),', [',np.min(lon_rho),', ... , ',np.max(lon_rho),']')
lon_u   = 0.5*(lon_rho[:,1:]+lon_rho[:,:-1]) 
#print('lon_u:',np.shape(lon_u),', [',np.min(lon_u),', ... , ',np.max(lon_u),']')
lon_v   =0.5*(lon_rho[1:,:]+lon_rho[:-1,:])
#print('lon_v:',np.shape(lon_v),', [',np.min(lon_v),', ... , ',np.max(lon_v),']')

#print('lat_rho:',np.shape(lat_rho),', [',np.min(lat_rho),', ... , ',np.max(lat_rho),']')
lat_u   = 0.5*(lat_rho[:,1:]+lat_rho[:,:-1])
#print('lat_u:',np.shape(lat_u),', [',np.min(lat_u),', ... , ',np.max(lat_u),']')
lat_v   = 0.5*(lat_rho[1:,:]+lat_rho[:-1,:])
#print('lat_v:',np.shape(lat_v),', [',np.min(lat_v),', ... , ',np.max(lat_v),']')
dlon=lon_rho[0,1]-lon_rho[0,0]
dlat=lat_rho[1,0]-lat_rho[0,0]

Zp1=np.linspace(-4,0.0,nk+1)
Zvec=0.5*(Zp1[1:]+Zp1[:-1])
depth_ROMS=np.zeros((nij,nij),dtype=float)
KBottomRUV=np.zeros((3,nij,nij), dtype=int)
mask_rho=np.zeros((nk,nij,nij),dtype=float)
for j in range(nij):
    for i in range(nij):
        depth_ROMS[j,i]=round(min(0,-3.0*(math.sin(math.pi*j/(nij-1.0))+math.sin(math.pi*i/(nij-1.0)))+2.0), 1)
depth_MITgcm=np.round(depth_ROMS*2.)/2.
for j in range(nij):
    for i in range(nij):
        mask_rho[np.argmin(abs(Zp1-(depth_MITgcm[j,i]))):,j,i]=1
mask_u=np.maximum(np.minimum((mask_rho[:,:,:-1]+mask_rho[:,:,1:])-1,1),0)
mask_v=np.maximum(np.minimum((mask_rho[:,:-1,:]+mask_rho[:,1:,:])-1,1),0)
for j in range(nij):
    for i in range(nij):
        try:KBottomRUV[0,j,i]=np.where(mask_rho[:,j,i]==1)[0][0]+1
        except:KBottomRUV[0,j,i]=nk+1
KBottomRUV[1,:,:-1]=np.maximum(KBottomRUV[0,:,:-1],KBottomRUV[0,:,1:])
KBottomRUV[1,:,nij-1]=KBottomRUV[0,:,nij-1]
KBottomRUV[2,:-1,:]=np.maximum(KBottomRUV[0,:-1,:],KBottomRUV[0,1:,:])
KBottomRUV[2,nij-1,:]=KBottomRUV[0,nij-1,:]
KBottomRUV[:,:,:]=np.maximum(1,KBottomRUV[:,:,:])

#print('depth_ROMS in range',np.min(depth_ROMS),np.max(depth_ROMS),depth_ROMS[np.arange(0,nij),np.arange(0,nij)])
#ax = plt.subplot()
#i_rho,j_rho=np.meshgrid( np.linspace(0,nij,nij) , np.linspace(0,nij,nij) )
#im=plt.pcolor(i_rho,j_rho,depth_MITgcm, edgecolors='k', linewidths=1,shading='nearest',snap=False)
#plt.axis('off')
#dij=0.5
#for ij in range(0,nij):
#    plt.annotate(str(mask_rho[:,ij,ij]),(i_rho[ij,ij],j_rho[ij,ij]))
#for ij in range(0,nij):
#    plt.annotate('R'+str(KBottomRUV[0,nij-1-ij,ij]),(i_rho[nij-1-ij,ij]-0.2*dij,j_rho[nij-1-ij,ij]-0.2*dij))
#for ij in range(0,nij):
#    plt.annotate('U'+str(KBottomRUV[1,int(nij/2),ij]),(i_rho[int(nij/2),ij]+0.4*dij,j_rho[int(nij/2),ij]-0.1*dij))
#for ij in range(0,nij):
#    plt.annotate('V'+str(KBottomRUV[2,ij,int(nij/2)]),(i_rho[ij,int(nij/2)]-0.2*dij,j_rho[ij,int(nij/2)]+0.4*dij))
#divider = make_axes_locatable(ax)
#cax = divider.append_axes("right", size="5%", pad=0.05)
#plt.colorbar(im, cax=cax)
#plt.show()


print(' ')
outfile=outdir+'/'+'GridforLTRANS-Zlevels.nc'
#print('')
print('-------- WRITING GRID FILE '+str(outfile)+'-----------')
if(netcdf_module=='netCDF4'):
  ncgridOUT= nc4.Dataset(outfile,"w")
else:
  ncgridOUT=NC.netcdf_file(outfile,"w")
ncgridOUT.createDimension("xi_rho",   nij)
ncgridOUT.createDimension("xi_u",     nij-1)
ncgridOUT.createDimension("xi_v",     nij)
ncgridOUT.createDimension("eta_rho",  nij)
ncgridOUT.createDimension("eta_u",    nij)
ncgridOUT.createDimension("eta_v",    nij-1)
ncgridOUT.createDimension("Z",        len(Zvec))
ncgridOUT.createDimension("Zi",       len(Zp1))
ncgridOUT.createDimension("RUV",   int(3))


ncvar=ncgridOUT.createVariable("Z"            , 'f', ('Z', ))
setattr(ncvar,'units'              , 'meters depth')
setattr(ncvar,'long_name'          , 'vertical coordinate of cell center'  )
setattr(ncvar,'standard_name'      , 'vertical coordinate of cell center'  )
ncvar[:]=Zvec

ncvar=ncgridOUT.createVariable("Zp1"            , 'f', ('Zi', ))
setattr(ncvar,'units'              , 'meters depth')
setattr(ncvar,'long_name'          , 'vertical coordinate of cell interface'  )
setattr(ncvar,'standard_name'      , 'vertical coordinate of cell interface'  )
ncvar[:]=Zp1

ncvar=ncgridOUT.createVariable("lon_rho", 'f', ('eta_rho', 'xi_rho',))
setattr(ncvar,'units'              , 'degrees_east')
setattr(ncvar,'long_name'          , 'lon_rho')
setattr(ncvar,'standard_name'      , 'lon_rho')
setattr(ncvar,'axis'               , 'X'        )
setattr(ncvar,'valid_min'          , -180.      )
setattr(ncvar,'valid_max'          , +180.      )
setattr(ncvar,'_CoordinateAxisType', 'longitude'      )
ncvar[:]=lon_rho[:,:]

ncvar=ncgridOUT.createVariable("lat_rho", 'f', ('eta_rho', 'xi_rho',))
setattr(ncvar,'units'              , 'degrees_north')
setattr(ncvar,'long_name'          , 'lat_rho')
setattr(ncvar,'standard_name'      , 'lat_rho')
setattr(ncvar,'axis'               , 'Y'        )
setattr(ncvar,'valid_min'          , -90.      )
setattr(ncvar,'valid_max'          , +90.      )
setattr(ncvar,'_CoordinateAxisType', 'latitude'      )
ncvar[:]=lat_rho[:,:]

ncvar=ncgridOUT.createVariable("mask_rho", 'f', ('Z', 'eta_rho', 'xi_rho',))
setattr(ncvar,'units'              , 'binary 0-1')
setattr(ncvar,'long_name'          , 'mask_rho')
setattr(ncvar,'standard_name'      , 'mask_rho')
setattr(ncvar,'axis'               , 'Y'        )
setattr(ncvar,'valid_min'          , 0.      )
setattr(ncvar,'valid_max'          , +1.      )
ncvar[:]=mask_rho

ncvar=ncgridOUT.createVariable("lon_u", 'f', ('eta_u', 'xi_u',))
setattr(ncvar,'units'              , 'degrees_east')
setattr(ncvar,'long_name'          , 'lon_u')
setattr(ncvar,'standard_name'      , 'lon_u')
setattr(ncvar,'axis'               , 'X'        )
setattr(ncvar,'valid_min'          , -180.      )
setattr(ncvar,'valid_max'          , +180.      )
setattr(ncvar,'_CoordinateAxisType', 'longitude'      )
ncvar[:]=lon_u[:,:]

ncvar=ncgridOUT.createVariable("lat_u", 'f', ('eta_u', 'xi_u',))
setattr(ncvar,'units'              , 'degrees_north')
setattr(ncvar,'long_name'          , 'lat_u')
setattr(ncvar,'standard_name'      , 'lat_u')
setattr(ncvar,'axis'               , 'Y'        )
setattr(ncvar,'valid_min'          , -90.      )
setattr(ncvar,'valid_max'          , +90.      )
setattr(ncvar,'_CoordinateAxisType', 'latitude'      )
ncvar[:]=lat_u[:,:]

ncvar=ncgridOUT.createVariable("mask_u", 'f', ('Z', 'eta_u', 'xi_u',))
setattr(ncvar,'units'              , 'binary 0-1')
setattr(ncvar,'long_name'          , 'mask_u')
setattr(ncvar,'standard_name'      , 'mask_u')
setattr(ncvar,'axis'               , 'Y'        )
setattr(ncvar,'valid_min'          , 0.      )
setattr(ncvar,'valid_max'          , +1.      )
ncvar[:]=mask_u

ncvar=ncgridOUT.createVariable("lon_v", 'f', ('eta_v', 'xi_v',))
setattr(ncvar,'units'              , 'degrees_east')
setattr(ncvar,'long_name'          , 'lon_v')
setattr(ncvar,'standard_name'      , 'lon_v')
setattr(ncvar,'axis'               , 'X'        )
setattr(ncvar,'valid_min'          , -180.      )
setattr(ncvar,'valid_max'          , +180.      )
setattr(ncvar,'_CoordinateAxisType', 'longitude'      )
ncvar[:]=lon_v[:,:]

ncvar=ncgridOUT.createVariable("lat_v", 'f', ('eta_v', 'xi_v',))
setattr(ncvar,'units'              , 'degrees_north')
setattr(ncvar,'long_name'          , 'lat_v')
setattr(ncvar,'standard_name'      , 'lat_v')
setattr(ncvar,'axis'               , 'Y'        )
setattr(ncvar,'valid_min'          , -90.      )
setattr(ncvar,'valid_max'          , +90.      )
setattr(ncvar,'_CoordinateAxisType', 'latitude'      )
ncvar[:]=lat_v[:,:]

ncvar=ncgridOUT.createVariable("mask_v", 'f', ('Z', 'eta_v', 'xi_v',))
setattr(ncvar,'units'              , 'binary 0-1')
setattr(ncvar,'long_name'          , 'mask_v')
setattr(ncvar,'standard_name'      , 'mask_v')
setattr(ncvar,'axis'               , 'Y'        )
setattr(ncvar,'valid_min'          , 0.      )
setattr(ncvar,'valid_max'          , +1.      )
ncvar[:]=mask_v

ncvar=ncgridOUT.createVariable("h"      , 'f', ('eta_rho', 'xi_rho',))
setattr(ncvar,'units'              , 'meters' )
setattr(ncvar,'long_name'          , 'depth'  )
setattr(ncvar,'standard_name'      , 'depth'  )
setattr(ncvar,'axis'               , 'Z'      )
setattr(ncvar,'valid_min'          , 0.       )
setattr(ncvar,'positive'           , 'down'   )
setattr(ncvar,'_CoordinateAxisType', 'heigth' )
setattr(ncvar,'_CoordinateAxiszisPositive', 'down')
ncvar[:]=-depth_MITgcm # must be positive

ncvar=ncgridOUT.createVariable("KBottomRUV", 'i', ('RUV', 'eta_rho', 'xi_rho',))
setattr(ncvar,'units'              , 'integer')
setattr(ncvar,'long_name'          , 'Bottom cell index (partial cell) for Rho,U,V nodes')
setattr(ncvar,'standard_name'      , 'Bottom cell index (partial cell)')
setattr(ncvar,'valid_min'          , 0       )
setattr(ncvar,'valid_max'          , nk+1      )
ncvar[:]=KBottomRUV[:,:,:]  # must be negative

ncgridOUT.close()

outfile=outdir+'/'+'GridforLTRANS-Slevels.nc'
#print('')
print('-------- WRITING GRID FILE '+str(outfile)+'-----------')
if(netcdf_module=='netCDF4'):
  ncgridOUT= nc4.Dataset(outfile,"w")
else:
  ncgridOUT=NC.netcdf_file(outfile,"w")
ncgridOUT.createDimension("xi_rho",   nij)
ncgridOUT.createDimension("xi_u",     nij-1)
ncgridOUT.createDimension("xi_v",     nij)
ncgridOUT.createDimension("eta_rho",  nij)
ncgridOUT.createDimension("eta_u",    nij)
ncgridOUT.createDimension("eta_v",    nij-1)

ncvar=ncgridOUT.createVariable("lon_rho", 'f', ('eta_rho', 'xi_rho',))
setattr(ncvar,'units'              , 'degrees_east')
setattr(ncvar,'long_name'          , 'lon_rho')
setattr(ncvar,'standard_name'      , 'lon_rho')
setattr(ncvar,'axis'               , 'X'        )
setattr(ncvar,'valid_min'          , -180.      )
setattr(ncvar,'valid_max'          , +180.      )
setattr(ncvar,'_CoordinateAxisType', 'longitude'      )
ncvar[:]=lon_rho[:,:]

ncvar=ncgridOUT.createVariable("lat_rho", 'f', ('eta_rho', 'xi_rho',))
setattr(ncvar,'units'              , 'degrees_north')
setattr(ncvar,'long_name'          , 'lat_rho')
setattr(ncvar,'standard_name'      , 'lat_rho')
setattr(ncvar,'axis'               , 'Y'        )
setattr(ncvar,'valid_min'          , -90.      )
setattr(ncvar,'valid_max'          , +90.      )
setattr(ncvar,'_CoordinateAxisType', 'latitude'      )
ncvar[:]=lat_rho[:,:]

ncvar=ncgridOUT.createVariable("mask_rho", 'f', ('eta_rho', 'xi_rho',))
setattr(ncvar,'units'              , 'binary 0-1')
setattr(ncvar,'long_name'          , 'mask_rho')
setattr(ncvar,'standard_name'      , 'mask_rho')
setattr(ncvar,'axis'               , 'Y'        )
setattr(ncvar,'valid_min'          , 0.      )
setattr(ncvar,'valid_max'          , +1.      )
ncvar[:]=mask_rho[-1,:,:]

ncvar=ncgridOUT.createVariable("lon_u", 'f', ('eta_u', 'xi_u',))
setattr(ncvar,'units'              , 'degrees_east')
setattr(ncvar,'long_name'          , 'lon_u')
setattr(ncvar,'standard_name'      , 'lon_u')
setattr(ncvar,'axis'               , 'X'        )
setattr(ncvar,'valid_min'          , -180.      )
setattr(ncvar,'valid_max'          , +180.      )
setattr(ncvar,'_CoordinateAxisType', 'longitude'      )
ncvar[:]=lon_u[:,:]

ncvar=ncgridOUT.createVariable("lat_u", 'f', ('eta_u', 'xi_u',))
setattr(ncvar,'units'              , 'degrees_north')
setattr(ncvar,'long_name'          , 'lat_u')
setattr(ncvar,'standard_name'      , 'lat_u')
setattr(ncvar,'axis'               , 'Y'        )
setattr(ncvar,'valid_min'          , -90.      )
setattr(ncvar,'valid_max'          , +90.      )
setattr(ncvar,'_CoordinateAxisType', 'latitude'      )
ncvar[:]=lat_u[:,:]

ncvar=ncgridOUT.createVariable("mask_u", 'f', ('eta_u', 'xi_u',))
setattr(ncvar,'units'              , 'binary 0-1')
setattr(ncvar,'long_name'          , 'mask_u')
setattr(ncvar,'standard_name'      , 'mask_u')
setattr(ncvar,'axis'               , 'Y'        )
setattr(ncvar,'valid_min'          , 0.      )
setattr(ncvar,'valid_max'          , +1.      )
ncvar[:]=mask_u[-1,:,:]

ncvar=ncgridOUT.createVariable("lon_v", 'f', ('eta_v', 'xi_v',))
setattr(ncvar,'units'              , 'degrees_east')
setattr(ncvar,'long_name'          , 'lon_v')
setattr(ncvar,'standard_name'      , 'lon_v')
setattr(ncvar,'axis'               , 'X'        )
setattr(ncvar,'valid_min'          , -180.      )
setattr(ncvar,'valid_max'          , +180.      )
setattr(ncvar,'_CoordinateAxisType', 'longitude'      )
ncvar[:]=lon_v[:,:]

ncvar=ncgridOUT.createVariable("lat_v", 'f', ('eta_v', 'xi_v',))
setattr(ncvar,'units'              , 'degrees_north')
setattr(ncvar,'long_name'          , 'lat_v')
setattr(ncvar,'standard_name'      , 'lat_v')
setattr(ncvar,'axis'               , 'Y'        )
setattr(ncvar,'valid_min'          , -90.      )
setattr(ncvar,'valid_max'          , +90.      )
setattr(ncvar,'_CoordinateAxisType', 'latitude'      )
ncvar[:]=lat_v[:,:]

ncvar=ncgridOUT.createVariable("mask_v", 'f', ('eta_v', 'xi_v',))
setattr(ncvar,'units'              , 'binary 0-1')
setattr(ncvar,'long_name'          , 'mask_v')
setattr(ncvar,'standard_name'      , 'mask_v')
setattr(ncvar,'axis'               , 'Y'        )
setattr(ncvar,'valid_min'          , 0.      )
setattr(ncvar,'valid_max'          , +1.      )
ncvar[:]=mask_v[-1,:,:]

ncvar=ncgridOUT.createVariable("h"      , 'f', ('eta_rho', 'xi_rho',))
setattr(ncvar,'units'              , 'meters' )
setattr(ncvar,'long_name'          , 'depth'  )
setattr(ncvar,'standard_name'      , 'depth'  )
setattr(ncvar,'axis'               , 'Z'      )
setattr(ncvar,'valid_min'          , 0.       )
setattr(ncvar,'positive'           , 'down'   )
setattr(ncvar,'_CoordinateAxisType', 'heigth' )
setattr(ncvar,'_CoordinateAxiszisPositive', 'down')
ncvar[:]=-depth_ROMS # must be positive

ncvar=ncgridOUT.createVariable("angle", 'f', ('eta_rho', 'xi_rho',))
setattr(ncvar,'units'              , 'radians')
setattr(ncvar,'long_name'          , 'angle between xi axis and east')
setattr(ncvar,'standard_name'      , 'angle')
ncvar[:]=np.zeros((nij,nij),dtype=float)

ncgridOUT.close()

def fstr(value):
    return "%.8f" % value
cellskip=6
Parfilename=outdir+'/'+'Iniparloc_every_'+str(cellskip)+'_rhowaternode.csv'
print('writing file '+Parfilename)
fpartloc = open(Parfilename,"w")
count=0
for j in range(1,nij-1,cellskip):
  for i in range(1,nij-1,cellskip):
     if(mask_rho[-1,j,i]==1):
          fpartloc.write(fstr(lon_rho[j,i])+ ', '+fstr(lat_rho[j,i])+', -0.1,  100, 101000 '+'\n')
          count +=1
fpartloc.close()
print(str(count)+' particles written in file ',Parfilename)

