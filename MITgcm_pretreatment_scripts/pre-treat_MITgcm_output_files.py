PACKAGEDIRECTORY="/galileo/home/userexternal/clauren1/LTRANS_Zlev/"
#
writeUniformIniParlocFile=True  
(io,i_f,istep)=(0,0,5)
(jo,j_f,jstep)=(0,0,5)
#
MITgcmdirectory=PACKAGEDIRECTORY+'/MITgcm_outputs/'
f_Eta='Eta.'
f_RHOA='RHOAnoma.'
f_U='U.'
f_V='V.'
f_W='W.'
f_S='S.'
f_T='T.'
f_KPPdiffS='KPPdiffS.'
f_EXFuwind='EXFuwind.'
f_EXFvwind='EXFvwind.'
PRECISION=4
identifier='boxes_NiNj64'
#############################################################
import numpy as np
import os,sys
import glob
import datetime
from numpy import array
import struct
import math
import scipy.io.netcdf as NC
############################
import tools_module as tools
#############################################################
griddirout=PACKAGEDIRECTORY+'/SIM/input/'
LTRANSrundir=PACKAGEDIRECTORY+'/SIM/rundir_'+identifier+'/'
OUTPUTDIR=PACKAGEDIRECTORY+'/SIM/output_'+identifier+'/'
#############################################################
if(True):
        SETUPEXISTS=False
        print 'creating '+PACKAGEDIRECTORY+'plot_scripts/setup_'+identifier+'.py and '+LTRANSrundir+'LTRANS_'+identifier+'.data by calling tools.getgridparamsfromSTDOUT(',MITgcmdirectory,',',identifier,')'
        if(not tools.is_non_zero_file(LTRANSrundir)):os.system('mkdir '+LTRANSrundir)
        #tools.runcommand('cp LTRANS_model_file.batch '+LTRANSrundir+'LTRANS_'+identifier+'.batch')
        #tools.runcommand('sed -i s/PACKAGEDIRECTORY/'+str(PACKAGEDIRECTORY.replace('/','!?'))+'/ '+LTRANSrundir+'LTRANS_'+identifier+'.batch')
        #tools.runcommand('sed -i s/LTRANSRUNDIR/'+str(LTRANSrundir.replace('/','!?'))+'/ '+LTRANSrundir+'LTRANS_'+identifier+'.batch')
        #tools.runcommand('sed -i "s/!?/\//g" '+LTRANSrundir+'LTRANS_'+identifier+'.batch')
        #tools.runcommand('sed -i s/IDENTIFIER/'+identifier+'/ '+LTRANSrundir+'LTRANS_'+identifier+'.batch')
        #tools.runcommand('sed -i s/IDENTIFIER/'+identifier+'/ '+LTRANSrundir+'LTRANS_'+identifier+'.batch')
        #tools.runcommand('sed -i s/SCRATCHDIRECTORY/'+str(OUTPUTDIR.replace('/','!?'))+'/ '+LTRANSrundir+'LTRANS_'+identifier+'.batch')
        #tools.runcommand('sed -i "s/!?/\//g" '+LTRANSrundir+'LTRANS_'+identifier+'.batch')
	tools.getgridparamsfromSTDOUT(MITgcmdirectory,identifier,LTRANSrundir,f_Eta=f_Eta,f_RHOA=f_RHOA,
	f_U=f_U,f_V=f_V,f_W=f_W,f_S=f_S,f_T=f_T,f_KPPdiffS=f_KPPdiffS,f_EXFuwind=f_EXFuwind,f_EXFvwind=f_EXFvwind,plotdir=PACKAGEDIRECTORY+'/plot_scripts/')

	execfile(PACKAGEDIRECTORY+'plot_scripts/setup_'+identifier+'.py', globals())
        setupanddatafilerewritten=True
        setupwriter=open(PACKAGEDIRECTORY+'plot_scripts/setup_'+identifier+'.py','a')
        setupwriter.write('\n')
        setupwriter.write('hydrodir="'+MITgcmdirectory+'" \n')
        setupwriter.write('rundir="'+LTRANSrundir+'" \n')
        setupwriter.write("boundsfile=rundir+'llbounds.bln' \n")
        setupwriter.write('k_rho='+str(nzrho_in)+' \n')


os.system('mkdir '+OUTPUTDIR)
os.system('mkdir '+OUTPUTDIR+'metadata') 
################################################################################
def fstr(value):
    return "%.8f" % value

def eightnumberstring(number):
        if number<10:         completenum='0000000'+str(number)
        elif number<100:      completenum='000000'+str(number)
        elif number<1000:     completenum='00000'+str(number)
        elif number<10000:    completenum='0000'+str(number)
        elif number<100000:   completenum='000'+str(number)
        elif number<1000000:  completenum='00'+str(number)
        elif number<10000000: completenum='0'+str(number)
        elif number<100000000:completenum=str(number)
        else:          print('error')
        #print  completefilename
	return completenum


def reduce_gridvar(gridelement,ny,nx,skip,ARRAY):
  if(len(np.shape(ARRAY))==2):
	nz=0
  else:
	nz=np.shape(ARRAY)[0]
  nx0=0
  ny0=0
  if(gridelement==1):
  	nx0=skip-1
  elif(gridelement==2):
  	ny0=skip-1
  print ny,nx,np.shape(ARRAY)
  if(nz==0):
	REDUCEDARRAY=np.zeros((ny,nx), dtype=ARRAY.dtype)
	for i in range(nx0,skip):
	  for j in range(ny0,skip):
	     REDUCEDARRAY[:ny,:nx] = REDUCEDARRAY[:ny,:nx]+   ARRAY[j::skip,i::skip]
  else:
	REDUCEDARRAY=np.zeros((nz,ny,nx), dtype=ARRAY.dtype)
	for i in range(nx0,skip):
	  for j in range(ny0,skip):
	     REDUCEDARRAY[:,:ny,:nx] = REDUCEDARRAY[:,:ny,:nx]+   ARRAY[:,j::skip,i::skip]
  return REDUCEDARRAY

def emptyfillednumberstring(number,stringsize):  
	outstring=str(number)
	len_num=len(outstring)
	if(len_num<stringsize):
		for count in range(len_num,stringsize): outstring=' '+outstring
	else:
		outstring=outstring[0:stringsize]
	return outstring

def preparemaskforLTRANS(mask):
  (nz,ny,nx)=tools.getnznynx(mask)
  if(np.sum(mask[0])>np.sum(mask[-1])): 
	MODEL='MITgcm'
	korder=-1
	k0=nz-1
	kF=-1
  else: 
	MODEL='LTRANS'  
	korder=1
	k0=0
	kF=nz
  modnodes=np.zeros((ny,nx),dtype=int)
  print "Edit Rho Mask - remove water nodes that don't have at least 2 water node neighbors on a ",MODEL," grid"
  totfounds=0
  founds=0
  for k in range(k0,kF,korder):
    prevfounds=founds
    founds=-1
    while (founds != 0) :
      prevfounds=founds
      founds=0
      for j in range(0,ny):
        for i in range(0,nx):
          #!Edit Rho Mask - remove nodes that don't have at least 2 neighbors
          #!  This removes a situation that may occur where an area that is 
          #!  within the created boundaries, is not covered by either the U 
          #!  or V grid. 
          if(mask[k,j,i]>0):
            count = 0
            numhoriz=0
	    numvert=0
            if(i> 0):
              if(mask[k,j,i-1]>0): 
		count = count + 1
		numhoriz+=1
            if(i<nx-1):
              if(mask[k,j,i+1]>0): 
		count = count + 1
		numhoriz+=1
            if(j> 0):
              if(mask[k,j-1,i]>0): 
		count = count + 1
		numvert+=1
            if(j<ny-1):
              if(mask[k,j+1,i]>0): 
		count = count + 1
		numvert+=1
            if(count < 2 
		or (j==0    and mask[k,j+1,i]==0)
		or (j==ny-1 and mask[k,j-1,i]==0)
		or (i==0    and mask[k,j,i+1]==0)
		or (i==nx-1 and mask[k,j,i-1]==0) 
		or (count==2 and (numvert==0 or numhoriz==0))  ):
              founds=founds+1
              if(founds==1):print 'At level '+str(k+1),': deleting cell Rho(i,j)=',
              if(count>=2) : print ' ('+str(i+1)+','+str(j+1)+')',
              if(MODEL=='MITgcm'):
                 mask[k:,j,i]       = 0                     
                 modnodes[j,i]=1
              else: 
                 mask[:k+1,j,i]     = 0                      
                 modnodes[j,i]=1
      if(founds>0):print ' deleted',founds,'nodes ;'
      totfounds=totfounds+founds
    if(prevfounds>0):print ' '
  print 'in total in 3d '+str(totfounds)+' nodes were deleted'
  print 'in total in 2d '+str(np.sum(modnodes))+' nodes were deleted'
  
  return (mask,modnodes)
######################################
i0=0
j0=0
full_nxrho_in = nxrho_in 
full_nyrho_in = nyrho_in 
full_nzrho_in = nzrho_in 

#################################################################
ZvecMIT=np.zeros((nzrho_in),dtype=float)  # cell-center vertical coordinates       
Zp1MIT=np.zeros((nzrho_in+1),dtype=float)  # cell-interface vertical coordinates
(x,y,rho,u,v)=(0,1,0,1,2)

for k in range(1,nzrho_in+1):
	Zp1MIT[k]=Zp1MIT[k-1]-delZ[k-1]
	ZvecMIT[k-1]=0.5*(Zp1MIT[k-1]+Zp1MIT[k]) 
print 'Zp1MIT=',Zp1MIT
#################################################################

nxu_in=nxrho_in-1
nyu_in=nyrho_in

nxv_in=nxrho_in
nyv_in=nyrho_in-1

nzrho_in_cut=nzrho_in  # can be modified to extract only upper levels of the grid.
nzrho_in_cut=min(nzrho_in_cut,nzrho_in)
nzrho_in_w=nzrho_in+1
nzrho_in_w_cut=nzrho_in_cut+1
print 'nxrho_in= '+str(nxrho_in)
print 'nyrho_in= '+str(nyrho_in)
print 'nzrho_in= '+str(nzrho_in)
tools.runcommand('sed -i "s/HYDROBYTES/'+str(PRECISION)+'/" '+LTRANSrundir+'LTRANS_'+identifier+'.data')
tools.runcommand('sed -i "s/TODAYSDATE/'+datetime.date.today().strftime("%d %B %Y")+'/" '+LTRANSrundir+'LTRANS_'+identifier+'.data')
tools.runcommand('sed -i s/INPUTDIRECTORY/'+str(griddirout.replace('/','!?'))+'/ '+LTRANSrundir+'LTRANS_'+identifier+'.data')
tools.runcommand('sed -i s/SCRATCHDIRECTORY/'+str(OUTPUTDIR.replace('/','!?'))+'/ '+LTRANSrundir+'LTRANS_'+identifier+'.data')
tools.runcommand('sed -i "s/!?/\//g" '+LTRANSrundir+'LTRANS_'+identifier+'.data')
######################################

################### READING COMPUTING  VERTICAL COORDINATES #####################
Zvec=np.zeros((nzrho_in), dtype=np.float32)
for k in range(nzrho_in-1,-1,-1):
	Zvec[k]=ZvecMIT[(nzrho_in-1)-k]
Zp1_mod=np.zeros((nzrho_in+1), dtype=np.float)
Zp1_mod[:]=999999.
Zp1=np.zeros((nzrho_in+1), dtype=np.float)
Zp1[nzrho_in]=0.0
for k in range(nzrho_in,-1,-1):
	Zp1[k] =Zp1MIT[(nzrho_in)-k]
h=np.zeros((nyrho_in,nxrho_in),dtype=float)
KBottomRUV=np.zeros((3,nyrho_in,nxrho_in), dtype=np.int)
#################### READING INVERSING WRITING GRID  ###########################
for nf in range(0,numgridfiles):
        #g_filenames =['XC.data','YC.data','hFacC.data','XG.data','YG.data','Depth.data'] 
	print 'read grid file '+str(nf)+' : '+g_filenames[nf]
	f = open(dirin+g_filenames[nf],"rb")
	#dataread = np.fromfile(f, dtype=np.float32)
	if(nf==nfhFacC):	
          nk=nzrho_in
          full_nk=full_nzrho_in
	else:			
          nk=1
          full_nk=1
        print g_filenames[nf][max(0,len(g_filenames[nf])-5):]
        if(PRECISION==8):  myfmt='d'*full_nk*full_nyrho_in*full_nxrho_in
        elif(PRECISION==4):myfmt='f'*full_nk*full_nyrho_in*full_nxrho_in
        else:
             print 'ERROR PRECISION should be 4 or 8, instead it was set to ',PRECISION
        dataread= struct.unpack(myfmt,f.read(PRECISION*full_nk*full_nyrho_in*full_nxrho_in))
	var=np.zeros((full_nk,full_nyrho_in,full_nxrho_in),dtype=float)
	for k in range(0,full_nk):
		for j in range(0,full_nyrho_in):
			var[k,j,:]=dataread[k*full_nyrho_in*full_nxrho_in+j*full_nxrho_in:k*full_nyrho_in*full_nxrho_in+j*full_nxrho_in+full_nxrho_in]
	f.close()     
	locnxrho_in=int(nxrho_in)
	locnyrho_in=int(nyrho_in)
	locnzrho_in=int(nzrho_in)
	VECz_name='Zmd000027'
	xdimname='xi_rho'
	ydimname='eta_rho'
	zdimname='s_rho'

	if(nf==nfXC) : 
		lon_rho =np.copy(var[0,0:nyrho_in,0:nxrho_in])
                print 'full lon min max=',np.min(var),np.max(var)
                print 'load lon min max=',np.min(lon_rho),np.max(lon_rho)
                LONMIN=np.min(lon_rho)
                tools.runcommand('sed -i s/LONMIN/'+str(LONMIN)+'/ '+LTRANSrundir+'LTRANS_'+identifier+'.data')
		print 'lon_rho:',np.shape(lon_rho),', [',np.min(lon_rho),', ... , ',np.max(lon_rho),']'
		lon_u   = 0.5*(lon_rho[:,1:]+lon_rho[:,:-1]) 
		#lon_u   = lon_rho[0:nyrho_in,0+1:nxrho_in]-1.0/reslon/2.0
		print 'lon_u:',np.shape(lon_u),', [',np.min(lon_u),', ... , ',np.max(lon_u),']'
		lon_v   =np.copy(var[0,0+1:nyrho_in,0:nxrho_in]) 
		print 'lon_v:',np.shape(lon_v),', [',np.min(lon_v),', ... , ',np.max(lon_v),']'
	elif(nf==nfYC) : 
		lat_rho =np.copy(var[0,0:nyrho_in,0:nxrho_in])
                print 'full lat min max=',np.min(var),np.max(var)
                print 'load lat min max=',np.min(lat_rho),np.max(lat_rho)
                LATMIN=np.min(lat_rho)
                tools.runcommand('sed -i s/LATMIN/'+str(LATMIN)+'/ '+LTRANSrundir+'LTRANS_'+identifier+'.data')
		print 'lat_rho:',np.shape(lat_rho),', [',np.min(lat_rho),', ... , ',np.max(lat_rho),']'
		lat_u   =np.copy(var[0,0:nyrho_in,0+1:nxrho_in]) 
		print 'lat_u:',np.shape(lat_u),', [',np.min(lat_u),', ... , ',np.max(lat_u),']'
		lat_v   =0.5*(lat_rho[1:,:]+lat_rho[:-1,:])
		print 'lat_v:',np.shape(lat_v),', [',np.min(lat_v),', ... , ',np.max(lat_v),']'
	elif(nf==nfhFacC):
		hFac   =np.copy(var[0:nzrho_in,0:nyrho_in,0:nxrho_in])
		print 'hfac',np.shape(hFac)
	elif(nf==nfXG) :  
		lon_u=np.copy(var[0,0:nyrho_in,0+1:nxrho_in])
#		lon_u   = 0.5*(lon_rho[:,1:]+lon_rho[:,:-1]) 
		print 'lon_u:',np.shape(lon_u),', [',np.min(lon_u),', ... , ',np.max(lon_u),']'
	elif(nf==nfYG) : 
		lat_v=np.copy(var[0,0+1:nyrho_in,0:nxrho_in])
		print 'lat_v:',np.shape(lat_v),', [',np.min(lat_v),', ... , ',np.max(lat_v),']'
	elif  (nf==nfDepth) : 
		h       =np.copy(var[0,0:nyrho_in,0:nxrho_in])
                print 'full h min max=',np.min(var),np.max(var)
                print 'load h min max=',np.min(h),np.max(h)
		print 'h:',np.shape(h),', min max= [',np.min(h),', ... , ',np.max(h),']'
write_xyz_ASCII=False
if(write_xyz_ASCII):
  fbati = open('Bathymetry_'+identifier+'.xyz',"w")
  fbati.write('   i,   j,         lon,         lat,      depth,\n')
  for j in range(0,len(lat_rho)):
    for i in range(0,len(lat_rho[0])):
      fbati.write(emptyfillednumberstring(i+1,4)+','+emptyfillednumberstring(j+1,4)+', '+fstr(lon_rho[j,i])+', '+fstr(lat_rho[j,i])+', '+fstr(h[j,i])+', \n')
  fbati.close()

if(setupanddatafilerewritten): 
  Polyfilename=griddirout+'Polygon_'+identifier+'_global.csv'
  print 'writing file '+Polyfilename
  tools.runcommand('sed -i s/POLYGONFILE/Polygon_'+identifier+'_global.csv/ '+LTRANSrundir+'LTRANS_'+identifier+'.data')
  fpoly = open(Polyfilename,"w")
  AVLON=(np.min(lon_rho)+np.max(lon_rho))/2.0
  AVLAT=(np.min(lat_rho)+np.max(lat_rho))/2.0
  fpoly.write(' 100000, '+fstr(AVLON)+ ', '+fstr(AVLAT)+ ', '+fstr(np.min(lon_rho))+ ', '+fstr(np.min(lat_rho)) +', \n')
  fpoly.write(' 100000, '+fstr(AVLON)+ ', '+fstr(AVLAT)+ ', '+fstr(np.min(lon_rho))+ ', '+fstr(np.max(lat_rho)) +', \n')
  fpoly.write(' 100000, '+fstr(AVLON)+ ', '+fstr(AVLAT)+ ', '+fstr(np.max(lon_rho))+ ', '+fstr(np.max(lat_rho)) +', \n')
  fpoly.write(' 100000, '+fstr(AVLON)+ ', '+fstr(AVLAT)+ ', '+fstr(np.max(lon_rho))+ ', '+fstr(np.min(lat_rho)) +', \n')
  fpoly.write(' 100000, '+fstr(AVLON)+ ', '+fstr(AVLAT)+ ', '+fstr(np.min(lon_rho))+ ', '+fstr(np.min(lat_rho)) +', \n')
  fpoly.close()

print '        modifying h to fit Zp1 for full cells ...'
for j in range(0,nyrho_in):
	for i in range(0,nxrho_in):
		if(hFac[0,j,i]!=0.0):
			for k in range(1,nzrho_in):
				kinv=(nzrho_in-1)-k
				if((hFac[k,j,i]==0.0 and hFac[k-1,j,i]==1.0) and h[j,i]!=-Zp1MIT[k]):
                                        #if(abs(h[j,i]+Zp1MIT[k])>0.1):  print i,j,' modify h=',h[j,i],' -->> ',-Zp1MIT[k]
                                        #print i,j,' modify h=',h[j,i],' -->> ',-Zp1MIT[k]
                                        h[j,i]=-Zp1MIT[k]
                                        #print h[j,i]
                                        #if(Zp1_mod[kinv+1]==999999.):
                                        #    print 'modify Zp1MIT[',k,']=',Zp1MIT[k],Zp1[kinv+1],' -->> ',-h[j,i]
                                        #    Zp1MIT[k]=-h[j,i]
                                        #    Zp1_mod[kinv+1]=-h[j,i]
                                        #    Zp1[kinv+1]=-h[j,i]
                                        #else:
                                        #    print ' ERROR SECOND MODIF Zp1MIT[',k,']=',Zp1MIT[k],Zp1[kinv+1],' XX>>!! ',h[j,i]
                                            
###############################################################

print '        computing KBottom...'
for j in range(0,nyrho_in):
        print 'At j=',j,', eventual water cells at the bottom of the bathymetry found at i=(', 
	for i in range(0,nxrho_in):
		KBottomRUV[0,j,i]=nzrho_in +(1) # +(1) as LTRANS is fortran, indexes from 1 to nzrho_in
		if(hFac[0,j,i]!=0.0):
			for k in range(1,nzrho_in):
				kinv=(nzrho_in-1)-k
				if(hFac[k,j,i]!=0.0 and hFac[k,j,i]!=1.0):  # PARTIAL BOTTOM CELL
					# MAKE IT FULL CELL
		 			KBottomRUV[0,j,i]=kinv +(1) # +(1) as LTRANS is fortran, indexes from 1 to nzrho_in
                                        if(h[j,i]<-Zp1MIT[k] or h[j,i]>-Zp1MIT[k+1]):
                                         print '--------------------------'
                                         print 'INCONSISTENCY i=',i,' j=',j,' k=',k,' kinv=',kinv +(1),'hfac[k-1,k-k+1]= (',hFac[max(0,k-1),j,i],hFac[k,j,i],hFac[min(nzrho_in-1,k+1),j,i],')'
                                         print 'h=',h[j,i],' Zp1MIT[',max(0,k-1),']=',Zp1MIT[max(0,k-1)],' Zp1MIT[',k,']=',Zp1MIT[k],' Zp1MIT[',k+1,']=',Zp1MIT[k+1]
                                         print '--------------------------'
                                         quit()
                  			break
				elif(hFac[k,j,i]==0.0 and hFac[k-1,j,i]!=0.0):
		 			KBottomRUV[0,j,i]=kinv+1 +(1) # +(1) as LTRANS is fortran, indexes from 1 to nzrho_in
                                        if(h[j,i]<-Zp1MIT[k] or h[j,i]>-Zp1MIT[k+1]):
                                         print '--------------------------'
                                         print 'INCONSISTENCY i=',i,' j=',j,' k=',k,' kinv=',kinv +(1),'hfac[k-1,k-k+1]= (',hFac[max(0,k-1),j,i],hFac[k,j,i],hFac[min(nzrho_in-1,k+1),j,i],')'
                                         print 'h=',h[j,i],' Zp1MIT[',max(0,k-1),']=',Zp1MIT[max(0,k-1)],' Zp1MIT[',k,']=',Zp1MIT[k],' Zp1MIT[',k+1,']=',Zp1MIT[k+1]
                                         print '--------------------------'
                                         quit()
                  			break
				elif(k==nzrho_in-1 and hFac[k,j,i]!=0.0): 
					print i,
					KBottomRUV[0,j,i]=kinv +(1) # +(1) as LTRANS is fortran, indexes from 1 to nzrho_in
                  			break					
				elif(hFac[k,j,i]==1.0):
					continue
				else:
					print ' // '	
					print 'unresolved situation at  i=',i,', j=',j,'  =>> programs STOPS'	
					quit()
        print ')' 
KBottomRUV[1,:,:-1]=np.maximum(KBottomRUV[0,:,:-1],KBottomRUV[0,:,1:])
KBottomRUV[1,:,nxrho_in-1]=KBottomRUV[0,:,nxrho_in-1]
KBottomRUV[2,:-1,:]=np.maximum(KBottomRUV[0,:-1,:],KBottomRUV[0,1:,:])
KBottomRUV[2,nyrho_in-1,:]=KBottomRUV[0,nyrho_in-1,:]
#
#--------------------------------------------------------------------------------------------
#################### COMPUTE  MASKS  (k direction is inversed compared to MITgcm) ######################
print '        computing masks...'

if(nxu_in != nxrho_in-1 or nyv_in !=nyrho_in-1 or nzrho_in_w != nzrho_in+1): stop
LTRANSmask_rho=np.zeros((nzrho_in  ,nyrho_in  ,nxrho_in  ), dtype=np.int)
for j in range(0,nyrho_in):
	for i in range(0,nxrho_in):
		if(KBottomRUV[0,j,i]==0 or KBottomRUV[1,j,i]==0 or KBottomRUV[2,j,i]==0):
			print 'ERROR KBottomRUV[0,j,i]==0 or KBottomRUV[1,j,i]==0 or KBottomRUV[2,j,i]==0'
			quit() 
		ko=KBottomRUV[0,j,i] 
		if( ko<nzrho_in):
			LTRANSmask_rho[ko-1 :nzrho_in   , j , i ] = 1  # WATER

#-------------------------------------------------------------------------------------------------------------
#tools.printmask(LTRANSmask_rho[-1,:70,600:750])#600:750],hFac[0,::12,600:750],hFac[0,::12,600:750])
(LTRANSmask_rho,modnodes1)=preparemaskforLTRANS(LTRANSmask_rho)
maskrho_in=LTRANSmask_rho[::-1]

h_orig=np.copy(h)
#-------------------------------------------------------------------------------------------------------------
#######################################################################################
#######################################################################################
skipindex=1
if(setupanddatafilerewritten):setupwriter.write('DomainMINMAX=['+str(np.min(lon_rho))+','+str(np.max(lon_rho))+','+str(np.min(lat_rho))+','+str(np.max(lat_rho))+'] \n')
sxy=np.zeros((2,3),dtype=int)
nxy=np.zeros((2,3),dtype=int)
(nxy[x,rho],nxy[x,u],nxy[x,v],nxy[y,rho],nxy[y,u],nxy[y,v])=(nxrho_in,nxu_in,nxv_in,nyrho_in,nyu_in,nyv_in)	
if(setupanddatafilerewritten): setupwriter.write('gridfilename="'+griddirout+'GridforLTRANS-'+identifier+'-c1.nc'+'" \n')	

(sxy[x,rho],sxy[x,u],sxy[x,v],sxy[y,rho],sxy[y,u],sxy[y,v])=(nxrho_in,nxu_in,nxv_in,nyrho_in,nyu_in,nyv_in)
maskrho=maskrho_in
masku=np.maximum(np.minimum((maskrho[:,:,:-1]+maskrho[:,:,1:])-1,1),0)
maskv=np.maximum(np.minimum((maskrho[:,:-1,:]+maskrho[:,1:,:])-1,1),0)
#
print '---------- compute sum of merged masks -----------'
print 'rho mask'
summaskrho=reduce_gridvar(rho,sxy[y,rho],sxy[x,rho],skipindex,maskrho)
(summaskrho,modnodes2)=preparemaskforLTRANS(summaskrho)
smaskrho=np.minimum(summaskrho[:,:sxy[y,rho],:sxy[x,rho]],1)
#
print 'u mask'
smasku=np.maximum(np.minimum((smaskrho[:,:,:-1]+smaskrho[:,:,1:])-1,1),0)
summasku=reduce_gridvar(u,sxy[y,u],sxy[x,u],skipindex,masku)
summasku=summasku*smasku
#
print 'v mask',
smaskv=np.maximum(np.minimum((smaskrho[:,:-1,:]+smaskrho[:,1:,:])-1,1),0)
summaskv=reduce_gridvar(v,sxy[y,v],sxy[x,v],skipindex,maskv)
summaskv=summaskv*smaskv
SumMask=[summaskrho,summasku,summaskv]
SMask=[smaskrho,smasku,smaskv]
Mask=[maskrho,masku,maskv]
print 'CHECK THAT PREPARING MASKS IMPACTS CORRECTLY ON SMASKU AND SMASKV'
#
print 'dimensions of new mask are Ny=',len(smaskrho),' and Nx=',len(smaskrho[0])
#
print '-------- merged mask : ---------------------------'
#tools.printmask(smaskrho[-14],smasku[-14],smaskv[-14])
smaskrho_out=smaskrho[::-1]                                # mask inversion 
smasku_out=smasku[::-1]                                    # mask inversion 
smaskv_out=smaskv[::-1]                                    # mask inversion
print '----------------------------------------------------' 
nk=0
ones=np.copy(maskrho[0])
ones[:,:]=1
sumones_rho=reduce_gridvar(rho,sxy[y,rho],sxy[x,rho],skipindex,ones)
ones=np.copy(masku[0])
ones[:,:]=1
sumones_u=reduce_gridvar(u,sxy[y,u],sxy[x,u],skipindex,ones)
ones=np.copy(maskv[0])
ones[:,:]=1
sumones_v=reduce_gridvar(v,sxy[y,v],sxy[x,v],skipindex,ones)
print '--------- DONE WITH MASKS-------------------------'
print ' '
# lon_rh: -----------------------------------------------------------------
var_ext=lon_rho
var_red=reduce_gridvar(rho,sxy[y,rho],sxy[x,rho],skipindex,var_ext)
lon_rho_out=var_red/np.maximum(sumones_rho[0],1)
#tools.printval(lon_rho_out,10,0)
# lat_rho: -----------------------------------------------------------------
var_ext=lat_rho
var_red=reduce_gridvar(rho,sxy[y,rho],sxy[x,rho],skipindex,var_ext)
lat_rho_out=var_red/np.maximum(sumones_rho[0],1)
#tools.printval(lat_rho_out,10,0)
# lon_u: -----------------------------------------------------------------
var_ext=lon_u
var_red=reduce_gridvar(u,sxy[y,u],sxy[x,u],skipindex,var_ext)
lon_u_out=var_red/np.maximum(sumones_u[0],1)
# lat_u: -----------------------------------------------------------------
var_ext=lat_u
var_red=reduce_gridvar(u,sxy[y,u],sxy[x,u],skipindex,var_ext)
lat_u_out=var_red/np.maximum(sumones_u[0],1)
# lon_v: -----------------------------------------------------------------
var_ext=lon_v
var_red=reduce_gridvar(v,sxy[y,v],sxy[x,v],skipindex,var_ext)
lon_v_out=var_red/np.maximum(sumones_v[0],1)
# lat_v: -----------------------------------------------------------------
var_ext=lat_v
var_red=reduce_gridvar(v,sxy[y,v],sxy[x,v],skipindex,var_ext)
lat_v_out=var_red/np.maximum(sumones_v[0],1)
#######################################################################################
print 'correcting h and KBottom to take in account the cells deleted'
count=0
countok=0
countfullcells=0
countpartcells=0
for j in range(0,sxy[y,rho]):
	for i in range(0,sxy[x,rho]):
                #if(modnodes1[j,i]>0):print ' node ij=',i,j,' was modified ',smaskrho_out[:,j,i]
                for k in range(0,nzrho_in):		
			if(smaskrho_out[k,j,i]==1 or k==nzrho_in-1):
                             if((k==nzrho_in-1 and KBottomRUV[0,j,i]!=k+(2) )):
                                  #print 'masks were modified in i,j=',i,j,' KBot=',KBottomRUV[0,j,i],' -->> k+2=',k+(2),'  h=',h[j,i],' -->> Z[k+1]=',-Zp1[k+1],' Zp1[k-1]=',Zp1[k-1],' Zp1[k]=',Zp1[k],' Zp1[k+1]=',Zp1[k+1]
                                  KBottomRUV[0,j,i]=k+(2)
                                  h[j,i]=-Zp1[k+1]
                                  count=count+1
                             elif((KBottomRUV[0,j,i]!=k+(1) and h[j,i]>-Zp1[k])):
                                  #print 'masks were modified in i,j=',i,j,' KBot=',KBottomRUV[0,j,i],' -->> k+1=',k+(1),'  h=',h[j,i],' -->> Z[k  ]=',-Zp1[k],' Zp1[k-1]=',Zp1[k-1],' Zp1[k]=',Zp1[k],' Zp1[k+1]=',Zp1[k+1]
                                  KBottomRUV[0,j,i]=k+(1)
                                  h[j,i]=-Zp1[k]
                                  count=count+1
                             elif((k==nzrho_in-1 and smaskrho_out[k,j,i]==0) and (KBottomRUV[0,j,i]==nzrho_in+1)):
                                  countok=countok+1
                             elif(h[j,i]>-Zp1[k] or h[j,i]<=-Zp1[k+1]):
                                  print 'ERROR    in i,j,k=',i,j,k,' KBot=',KBottomRUV[0,j,i],' != ',k+(1),'  h=',h[j,i],' Zp1[k-1]=',Zp1[k-1],' Zp1[k]=',Zp1[k],' Zp1[k+1]=',Zp1[k+1],' mask=',smaskrho_out[k,j,i],'   ##################'
                                  quit()
                             else:
                                  if(abs(h[j,i]+Zp1[k])<=0.001):
                                    countfullcells+=1
                                    #if(hFac[(nzrho_in-1)-k,j,i]!=1):print ' full    bottom cell found in i,j=',i,j,' KBot=',KBottomRUV[0,j,i],'  h=',h[j,i],' Z[k  ]=',-Zp1[k],' Zp1[k-1]=',Zp1[k-1],' Zp1[k]=',Zp1[k],' Zp1[k+1]=',Zp1[k+1],' diff=',abs(h[j,i]+Zp1[k]),' hFac=',hFac[(nzrho_in-1)-k,j,i]
                                  else:
                                    countpartcells+=1
                                    if(hFac[(nzrho_in-1)-k,j,i]==1):print ' partial bottom cell found in i,j=',i,j,' KBot=',KBottomRUV[0,j,i],'  h=',h[j,i],' Z[k  ]=',-Zp1[k],' Zp1[k-1]=',Zp1[k-1],' Zp1[k]=',Zp1[k],' Zp1[k+1]=',Zp1[k+1],' diff=',abs(h[j,i]+Zp1[k]),' hFac=',hFac[(nzrho_in-1)-k,j,i]
                             break
print count,' nodes have now corrected h and KBottom values, ',countok,' were ok, including ',countfullcells,' full cells and ',countpartcells,' partial cells'
#h=h_orig
KBottomRUV[1,:,:-1]=np.maximum(KBottomRUV[0,:,:-1],KBottomRUV[0,:,1:])
KBottomRUV[1,:,nxrho_in-1]=KBottomRUV[0,:,nxrho_in-1]
KBottomRUV[2,:-1,:]=np.maximum(KBottomRUV[0,:-1,:],KBottomRUV[0,1:,:])
KBottomRUV[2,nyrho_in-1,:]=KBottomRUV[0,nyrho_in-1,:]
for j in range(0,sxy[y,rho]):
	for i in range(0,sxy[x,rho]):
		KBottomRUV[:,j,i]=np.maximum(1,KBottomRUV[:,j,i]-(nzrho_in-nzrho_in_cut))
#######################################################################################
print(' ') 
outfile='GridforLTRANS-'+identifier+'-c'+str(skipindex)+'.nc'
if(setupanddatafilerewritten): setupwriter.write('gridfilename="'+griddirout+outfile+'" \n')
print ''
print '-------- WRITING GRID FILE '+str(outfile)+'-----------'
ncgridOUT= NC.netcdf_file(griddirout+outfile,"w")
ncgridOUT.createDimension("xi_rho",   sxy[x,rho])
ncgridOUT.createDimension("xi_u",     sxy[x,u])
ncgridOUT.createDimension("xi_v",     sxy[x,v])
ncgridOUT.createDimension("eta_rho",  sxy[y,rho])
ncgridOUT.createDimension("eta_u",    sxy[y,u])
ncgridOUT.createDimension("eta_v",    sxy[y,v])
ncgridOUT.createDimension("Z",        nzrho_in_cut) 
ncgridOUT.createDimension("Zi",       nzrho_in_cut+1) 
ncgridOUT.createDimension("RUV",   int(3))


ncvar=ncgridOUT.createVariable("Z"            , 'f', ('Z', ))
setattr(ncvar,'units'              , 'meters depth')
setattr(ncvar,'long_name'          , 'vertical coordinate of cell center'  )
setattr(ncvar,'standard_name'      , 'vertical coordinate of cell center'  )
ncvar[:]=Zvec[nzrho_in-nzrho_in_cut:nzrho_in]

ncvar=ncgridOUT.createVariable("Zp1"            , 'f', ('Zi', ))
setattr(ncvar,'units'              , 'meters depth')
setattr(ncvar,'long_name'          , 'vertical coordinate of cell interface'  )
setattr(ncvar,'standard_name'      , 'vertical coordinate of cell interface'  )
ncvar[:]=Zp1[nzrho_in-nzrho_in_cut:nzrho_in_w]

ncvar=ncgridOUT.createVariable("lon_rho", 'f', ('eta_rho', 'xi_rho',))
setattr(ncvar,'units'              , 'degrees_east')
setattr(ncvar,'long_name'          , 'lon_rho')
setattr(ncvar,'standard_name'      , 'lon_rho')
setattr(ncvar,'axis'               , 'X'        )
setattr(ncvar,'valid_min'          , -180.      )
setattr(ncvar,'valid_max'          , +180.      )
setattr(ncvar,'_CoordinateAxisType', 'longitude'      )
ncvar[:]=lon_rho_out[:,:]

ncvar=ncgridOUT.createVariable("lat_rho", 'f', ('eta_rho', 'xi_rho',))
setattr(ncvar,'units'              , 'degrees_north')
setattr(ncvar,'long_name'          , 'lat_rho')
setattr(ncvar,'standard_name'      , 'lat_rho')
setattr(ncvar,'axis'               , 'Y'        )
setattr(ncvar,'valid_min'          , -90.      )
setattr(ncvar,'valid_max'          , +90.      )
setattr(ncvar,'_CoordinateAxisType', 'latitude'      )
ncvar[:]=lat_rho_out[:,:]

ncvar=ncgridOUT.createVariable("mask_rho", 'f', ('Z', 'eta_rho', 'xi_rho',))
setattr(ncvar,'units'              , 'binary 0-1')
setattr(ncvar,'long_name'          , 'mask_rho')
setattr(ncvar,'standard_name'      , 'mask_rho')
setattr(ncvar,'axis'               , 'Y'        )
setattr(ncvar,'valid_min'          , 0.      )
setattr(ncvar,'valid_max'          , +1.      )
ncvar[:]=smaskrho_out[nzrho_in-nzrho_in_cut:nzrho_in,:,:]

ncvar=ncgridOUT.createVariable("lon_u", 'f', ('eta_u', 'xi_u',))
setattr(ncvar,'units'              , 'degrees_east')
setattr(ncvar,'long_name'          , 'lon_u')
setattr(ncvar,'standard_name'      , 'lon_u')
setattr(ncvar,'axis'               , 'X'        )
setattr(ncvar,'valid_min'          , -180.      )
setattr(ncvar,'valid_max'          , +180.      )
setattr(ncvar,'_CoordinateAxisType', 'longitude'      )
ncvar[:]=lon_u_out[:,:]

ncvar=ncgridOUT.createVariable("lat_u", 'f', ('eta_u', 'xi_u',))
setattr(ncvar,'units'              , 'degrees_north')
setattr(ncvar,'long_name'          , 'lat_u')
setattr(ncvar,'standard_name'      , 'lat_u')
setattr(ncvar,'axis'               , 'Y'        )
setattr(ncvar,'valid_min'          , -90.      )
setattr(ncvar,'valid_max'          , +90.      )
setattr(ncvar,'_CoordinateAxisType', 'latitude'      )
ncvar[:]=lat_u_out[:,:]

ncvar=ncgridOUT.createVariable("mask_u", 'f', ('Z', 'eta_u', 'xi_u',))
setattr(ncvar,'units'              , 'binary 0-1')
setattr(ncvar,'long_name'          , 'mask_u')
setattr(ncvar,'standard_name'      , 'mask_u')
setattr(ncvar,'axis'               , 'Y'        )
setattr(ncvar,'valid_min'          , 0.      )
setattr(ncvar,'valid_max'          , +1.      )
ncvar[:]=smasku_out[nzrho_in-nzrho_in_cut:nzrho_in,:,:]

ncvar=ncgridOUT.createVariable("lon_v", 'f', ('eta_v', 'xi_v',))
setattr(ncvar,'units'              , 'degrees_east')
setattr(ncvar,'long_name'          , 'lon_v')
setattr(ncvar,'standard_name'      , 'lon_v')
setattr(ncvar,'axis'               , 'X'        )
setattr(ncvar,'valid_min'          , -180.      )
setattr(ncvar,'valid_max'          , +180.      )
setattr(ncvar,'_CoordinateAxisType', 'longitude'      )
ncvar[:]=lon_v_out[:,:]

ncvar=ncgridOUT.createVariable("lat_v", 'f', ('eta_v', 'xi_v',))
setattr(ncvar,'units'              , 'degrees_north')
setattr(ncvar,'long_name'          , 'lat_v')
setattr(ncvar,'standard_name'      , 'lat_v')
setattr(ncvar,'axis'               , 'Y'        )
setattr(ncvar,'valid_min'          , -90.      )
setattr(ncvar,'valid_max'          , +90.      )
setattr(ncvar,'_CoordinateAxisType', 'latitude'      )
ncvar[:]=lat_v_out[:,:]

ncvar=ncgridOUT.createVariable("mask_v", 'f', ('Z', 'eta_v', 'xi_v',))
setattr(ncvar,'units'              , 'binary 0-1')
setattr(ncvar,'long_name'          , 'mask_v')
setattr(ncvar,'standard_name'      , 'mask_v')
setattr(ncvar,'axis'               , 'Y'        )
setattr(ncvar,'valid_min'          , 0.      )
setattr(ncvar,'valid_max'          , +1.      )
ncvar[:]=smaskv_out[nzrho_in-nzrho_in_cut:nzrho_in,:,:]

ncvar=ncgridOUT.createVariable("h"      , 'f', ('eta_rho', 'xi_rho',))
setattr(ncvar,'units'              , 'meters' )
setattr(ncvar,'long_name'          , 'depth'  )
setattr(ncvar,'standard_name'      , 'depth'  )
setattr(ncvar,'axis'               , 'Z'      )
setattr(ncvar,'valid_min'          , 0.       )
setattr(ncvar,'positive'           , 'down'   )
setattr(ncvar,'_CoordinateAxisType', 'heigth' )
setattr(ncvar,'_CoordinateAxiszisPositive', 'down')
ncvar[:]=np.minimum(h[:,:],-np.min(Zp1[nzrho_in-nzrho_in_cut:nzrho_in_w]))

ncvar=ncgridOUT.createVariable("KBottomRUV", 'i', ('RUV', 'eta_rho', 'xi_rho',))
setattr(ncvar,'units'              , 'integer')
setattr(ncvar,'long_name'          , 'Bottom cell index (partial cell) for Rho,U,V nodes')
setattr(ncvar,'standard_name'      , 'Bottom cell index (partial cell)')
setattr(ncvar,'valid_min'          , 0       )
setattr(ncvar,'valid_max'          , nzrho_in_cut+1      )
ncvar[:]=KBottomRUV[:,:,:]


ncgridOUT.close()
if(writeUniformIniParlocFile):
   Parfilename=griddirout+'Iniparloc_'+identifier+'_every_'+str(istep)+'-i_'+str(jstep)+'-j_rhowaternode.csv'
   print 'writing file '+Parfilename
   fpartloc = open(Parfilename,"w")
   if(setupanddatafilerewritten): setupwriter.write('driftfilename="'+Parfilename+'" \n')

   count=0
   io=max(1,io)
   if (i_f<io+1):
      i_f=sxy[x,rho]-1
   else:
      i_f=min(i_f,sxy[x,rho]-1)
   jo=max(1,jo)
   if (j_f<jo+1):
      j_f=sxy[y,rho]-1
   else:
      j_f=min(j_f,sxy[y,rho]-1)

   print 'for j in range(',jo,',',j_f,',',jstep,'):'
   print '  for i in range(',io,',',i_f,',',istep,'):'
   print '    smaskrho_out[-1,j,i]=',smaskrho_out[-1,j,i]
   for j in range(jo,j_f,jstep):
     for i in range(io,i_f,istep):
	if(smaskrho_out[-1,j,i]==1):
             fpartloc.write(fstr(lon_rho_out[j,i])+ ', '+fstr(lat_rho_out[j,i])+', -0.1,  100, 101000 '+'\n')
             count +=1
   fpartloc.close()
   print str(count)+' particles written in file ',Parfilename
   if(setupanddatafilerewritten):
    tools.runcommand('sed -i s/NUMPAR/'+str(count)+'/ '+LTRANSrundir+'LTRANS_'+identifier+'.data')
    tools.runcommand('sed -i s/INPUTPARTLOCFILE/Iniparloc_'+identifier+'_every_'+str(istep)+'-i_'+str(jstep)+'-j_rhowaternode/ '+LTRANSrundir+'LTRANS_'+identifier+'.data')
   else:
    print 'WARNING NUMPAR and new PartLoc File to be actualised in '+LTRANSrundir+'LTRANS_'+identifier+'.data'


	
if(setupanddatafilerewritten):setupwriter.close()
print 'DONE'
