MODULE HYDRO_MOD                                                               

#define IMIN 1
#define IMAX 2
#define JMIN 3
#define JMAX 4
#define RNODE 1
#define UNODE 2
#define VNODE 3

!  This module handles all the input from the hydrodynamic NetCDF input files.
!  It is the only module that interacts with NetCDF input files.  It contains
!  all the variables read in from the NetCDF files.  It also contains all the
!  information and variables related to the grid elements.
!
!  Created by:            Zachary Schlag        
!  Created on:            07 Aug 2008
!  Last Modified on:         Feb 2013
  IMPLICIT NONE
  PRIVATE
  SAVE

  INTEGER :: iint, &  !Keeps track of the input file, 0=file 1, 1=file 2, etc.
             stepf, & !Keeps track of the forward time step
             countfilenum
  !Used for reading in NetCDF variables one time step at a time
  INTEGER :: STARTr(4),COUNTr(4),STARTz(3),COUNTz(3)

  !Keeps track of the Rho, U, and V element that each particle is in
  INTEGER, ALLOCATABLE, DIMENSION(:) :: P_r_element,P_u_element,P_v_element
!--- CL-OGS: Keeps track of the k level where are the particles
  INTEGER, ALLOCATABLE, DIMENSION(:) :: P_klev,P_klev_old
!--- CL-OGS: store interpolation coeffficients, nodes lists and number of
!            rho, u and v land nodes per vertical level 
!            where data from the closest water nodes must be copied
  INTEGER, DIMENSION(3) :: NUM_COPNOD   ! Stores the tot number of nodes R,U and V to be updated (ie the number of id-numbers written in Node_COPNOD
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: Node_COPNOD ! Stores the id-number of every updated R,U,V node
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: Nghb_COPNOD ! Stores for every updated R,U,V node its neighbors id-numbers
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: Coef_COPNOD ! Stores for every updated R,U,V node and for every neighbor the averaging coefficient
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: klev_COPNOD ! Stores for every updated R,U,V node the k-level range affected by the averaging proccess

!--- CL-OGS:
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: OMP_ruv

  !These variables keep track of the interpolation method and weights
  INTEGER :: numthreads

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: OMP_Xpar,OMP_Ypar
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: OMP_tOK
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: OMP_t,OMP_u
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: OMP_Wgt
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: OMP_Pint

  !S-Level location variables
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SC,CS,SCW,CSW
!--- CL-OGS: depth in meters of the Z-grid vertical levels at rho and w nodes
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: ZC,ZW 
!--- CL-OGS: vertical level number of the bottom (Z-grid)
  INTEGER         , ALLOCATABLE, DIMENSION(:,:,:) :: BottomK 
!--- CL-OGS: if bottom cell was partial in MITgcm re-interpolate currents at new cell center
!--- CL-OGS: vertical level number of the bottom (Z-grid)
  INTEGER, ALLOCATABLE, DIMENSION(:) :: BottomKatRnode,BottomKatUnode,BottomKatVnode
!--- CL-OGS: form of the 4 nodes elements imported from boundary module  
!--- CL-OGS: to compute Z-grid local depth 
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: eleform 
!--- CL-OGS: save X and Y particles positions when setEle is called to write
!--- CL-OGS: error scripts when particles get out of the domain
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::      &
              Xpar_at_setEle,Ypar_at_setEle,Zpar_at_setEle
!--- CL-OGS: import x and y bounds positions from boundary module to write
!--- CL-OGS: error scripts when particles get out of the domain
   DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: bnd_x,bnd_y
!--- CL-OGS
   INTEGER, ALLOCATABLE, DIMENSION(:):: nbounds

  !Depth at each rho node location
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: depth
  !Depth of the modified Zgrid bathimetry at rho, rho E (element center), U and V nodes
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: depthR,depthE,depthU,depthV
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: GrainSize
  
  !read in zeta,salinity,temperature,vertical diffusivity, and U,V,W velocities 
  !  at hydrodynamic back, center, and forward time
  INTEGER :: t_b,t_c,t_f,t_ijruv(4,3)
  !DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:)  :: t_Swdown
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:)  :: t_zeta
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: t_salt,t_temp,t_KH,       &
                     t_den,t_Uvel,t_Vvel,t_Wvel
!      *****   IMIOM      *****
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: t_hsig                                          ! significant wave height
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: t_tm01                                          ! mean wave period
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: t_uwind, t_vwind,t_iwind                       ! 10m wind components and intensity
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: t_pdir                                          ! principle wave direction
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: t_wlen                                          ! mean wave length
  CHARACTER(len=200) :: swannm
!      ***** END IMIOM *****

!--- CL-OGS: keep track of the nodes where fields were update to make sure that particles do not exit the updated buffer
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: updatenodesbuffer
!--- CL-OGS: vertical level numbers for MITgcm hydro using varying water cell at each level
  INTEGER :: us_tridim, ws_tridim

  !Rho, U, and V grid wet elements(four node numbers that make up the element)
  !  (wet means at least 1 node is masked as water)
!  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: RE,UE,VE
!--- CL-OGS: extension to third dimension of grid wet elements lists
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: RE,UE,VE
!--- CL-OGS: creation of k-level dependant number of wet elements to replace rho_elements,u_elements,v_elements of the original LTRANS v2b version
  INTEGER, ALLOCATABLE, DIMENSION(:) :: rho_kwele ! number of rho elements with at least 1 vertex is water up to plane k-1
  INTEGER, ALLOCATABLE, DIMENSION(:) :: u_kwele   !            u 
  INTEGER, ALLOCATABLE, DIMENSION(:) :: v_kwele   !            v
  !For each element, a list containing itself and all the elements that share a 
  !  node with that element;  used to speed up determining which element the 
  !  particle has moved to, if it has moved at all
!  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: r_Adjacent,u_Adjacent,v_Adjacent
!--- CL-OGS: extension to third dimension of of the Adjacent elements lists
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: r_Adjacent,u_Adjacent,v_Adjacent

  !X/Y location of all the Rho,U,V grid nodes, and the angle between
  !  x-coordinate and true east direction (radian)
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: rho_angle,rx,ry,ux,uy,vx,vy
!  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: r_ele_x,r_ele_y,u_ele_x,   &
!                                                   u_ele_y,v_ele_x,v_ele_y
!--- CL-OGS:  extension to third dimension and modification of the name of the 
!--- CL-OGS:  variables representing positions in x and y of the wet elements
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: r_kwele_x,r_kwele_y,      &
                         u_kwele_x, u_kwele_y,v_kwele_x,v_kwele_y

  !U, and V grid metric node locations, and Rho grid masking
  !  These variables are shared with boundary_module.f90 to create the model 
  !  boundaries
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)  :: x_u,y_u,x_v,y_v
!--- CL-OGS:  extension to third dimension of the masks
!  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:)  :: m_r,m_u,m_v
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:)  :: m_r,m_u,m_v
!--- CL-OGS:  extension to third dimension of the masks
!  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: mask_rho
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: mask_rho

!--- CL-OGS:  adding a dimension to account for the third dimension vertical level for the node-number-dependant-masks
!  INTEGER, ALLOCATABLE, DIMENSION( : ) :: rho_mask,u_mask,v_mask   !ewn.v.2
  INTEGER, ALLOCATABLE, DIMENSION( :,: ) :: rho_mask,u_mask,v_mask   !ewn.v.2
!--- CL-OGS: creating lists storing for each element its distance in Elements from land for stranding processes
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: coast_eleij !0 (land)-1 (border)-2(inner border)-3(touches inner border)...
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: coast_ele !0 (land)-1 (border)-2(inner border)-3(touches inner border)...
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: DeadOrOut


  !Keeps track if the grid has been read in yet or not
  !  If the grid hasn't been read in, the boundaries can't be made
  LOGICAL :: GRD_SET = .FALSE.

  !The concatenated hydrodynamic input file name
  CHARACTER(len=200) :: filenm

  !Counters for NetCDF files
  INTEGER :: NCcount,NCstart,prcount

  INTEGER :: fpy
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: lastplotpos
  !The following procedures have been made public:
  PUBLIC :: initGrid,initHydro,updateHydro,setEle,setEle_all,setInterp,        &
    getInterp,WCTS_ITPI,getSlevel,getWlevel,getMask_Rho,getUVxy,        &
    getR_ele,getP_r_element,finHydro,initNetCDF,createNetCDF,writeNetCDF,      &
    getKRlevel,getDepth,filenm,getP_klev,& !--- CL-OGS
    seteleform,setbounds,setnodesdepth,setDeadOrOut           !--- CL-OGS
     !,getKWlevel,set_closernode,outputdetails_closernode
CONTAINS


  SUBROUTINE initGrid()
    !This subroutine reads in the grid information and with it creates all the 
    !  element variables
    USE PARAM_MOD, ONLY: numpar,ui,vi,uj,vj,us,ws,rho_nodes,u_nodes,v_nodes,   &
        max_rho_elements,max_u_elements,    &
        max_v_elements,NCgridfile,  &
        Zgrid,ADJele_file,ADJele_fname,BoundaryBLNs,                           & !--- CL-OGS
        filestep,Vtransform,Wind,GrainSize_fname,read_GrainSize,         & !--- CL-OGS
        OutDir,NCOutFile,Zgrid_depthinterp,WindIntensity,filenum                         !--- CL-OGS
!    USE CONVERT_MOD, ONLY: lon2x,lat2y                                          !--- CL-OGS
    USE CONVERT_MOD, ONLY: lon2x,lat2y,x2lon,y2lat                               !--- CL-OGS
    USE netcdf
    !$ use OMP_LIB          
#include "VAR_IDs.h"
    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

    INTEGER :: STATUS,NCID,VID

    INTEGER :: i,j,m,count,inele
    INTEGER :: countele                                                          !--- CL-OGS   
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: euleriandepth,  &
                                x_rho,y_rho,angle,GrainSize_tmp
!--- CL-OGS:  extension to third dimension of the masks
!    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: mask_u, mask_v
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: mask_u, mask_v
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: lon_rho,lat_rho,lon_u,    &
                                                     lat_u,lon_v,lat_v
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: r_ele,u_ele,v_ele
    !INTEGER, ALLOCATABLE, DIMENSION( : ) :: rho_mask,u_mask,v_mask  !ewn.v.2
    INTEGER :: k,nfilesin,nf,ios,waiting,nodestocopy,kbot,kmax,maxnodestocopy     !--- CL-OGS 
    INTEGER :: old_i,old_j,old_count                                              !--- CL-OGS 
    DOUBLE PRECISION :: summask 
    DOUBLE PRECISION,DIMENSION(4) :: tmpcoef,oldtmpcoef
    !ALLOCATE MODULE VARIABLES
    !--- CL-OGS :
    ! CL-OGS: Added for OMP
    numthreads=1
    !$OMP PARALLEL
     !$OMP MASTER 
      !$ numthreads=OMP_GET_NUM_THREADS ()
      !$ write(*,*)'in hydrodynamic_module OMP_NUM_THREADS=',numthreads
     !$OMP END MASTER
    !$OMP END PARALLEL
    ALLOCATE(OMP_ruv(4,3,numthreads))
    ALLOCATE(OMP_Xpar(numthreads),OMP_Ypar(numthreads))
    ALLOCATE(OMP_tOK(3,numthreads),OMP_t(3,numthreads),OMP_u(3,numthreads))
    ALLOCATE(OMP_Wgt(4,3,numthreads),OMP_Pint(4,4,3,numthreads))
    ! CL-OGS: Added for OMP - end


    if(Zgrid)then
      us_tridim=us
      ws_tridim=ws
      ALLOCATE(ZC(us))
      ALLOCATE(ZW(ws))
      ALLOCATE(BottomKatRnode(rho_nodes))
      ALLOCATE(BottomKatUnode(  u_nodes))
      ALLOCATE(BottomKatVnode(  v_nodes))
      ALLOCATE(BottomK(vi,uj,3))
    else
      Zgrid_depthinterp=.FALSE.
      us_tridim=1
      ws_tridim=1
      ALLOCATE(SC(us))
      ALLOCATE(CS(us))
      ALLOCATE(SCW(ws))
      ALLOCATE(CSW(ws))
    endif
    ALLOCATE(Xpar_at_setEle(numpar)) !--- CL-OGS
    ALLOCATE(Ypar_at_setEle(numpar)) !--- CL-OGS
    ALLOCATE(Zpar_at_setEle(numpar)) !--- CL-OGS
    ALLOCATE(P_klev(numpar))
    ALLOCATE(P_klev_old(numpar))
    ALLOCATE(DeadOrOut(numpar))
    ALLOCATE(depth(rho_nodes))
    if(read_GrainSize)  ALLOCATE(GrainSize(rho_nodes))  !--- CL-OGS:  for behavior 8
    ALLOCATE(RE(4,max_rho_elements,us_tridim))          !--- CL-OGS:  extension to third dimension 
    ALLOCATE(UE(4,max_u_elements,us_tridim))            !--- CL-OGS:  extension to third dimension  
    ALLOCATE(VE(4,max_v_elements,us_tridim))            !--- CL-OGS:  extension to third dimension 
    ALLOCATE(r_Adjacent(max_rho_elements,10,us_tridim)) !--- CL-OGS:  extension to third dimension 
    ALLOCATE(u_Adjacent(max_u_elements,10,us_tridim))   !--- CL-OGS:  extension to third dimension 
    ALLOCATE(v_Adjacent(max_v_elements,10,us_tridim))   !--- CL-OGS:  extension to third dimension 
    ALLOCATE(rho_angle(rho_nodes))
    ALLOCATE(rx(rho_nodes))
    ALLOCATE(ry(rho_nodes))
    ALLOCATE(ux(u_nodes))
    ALLOCATE(uy(u_nodes))
    ALLOCATE(vx(v_nodes))
    ALLOCATE(vy(v_nodes))
    ALLOCATE(r_kwele_x(4,max_rho_elements,us_tridim))    !--- CL-OGS:  extension to third dimension 
    ALLOCATE(r_kwele_y(4,max_rho_elements,us_tridim))    !--- CL-OGS:  extension to third dimension 
    ALLOCATE(u_kwele_x(4,max_u_elements,us_tridim))      !--- CL-OGS:  extension to third dimension 
    ALLOCATE(u_kwele_y(4,max_u_elements,us_tridim))      !--- CL-OGS:  extension to third dimension 
    ALLOCATE(v_kwele_x(4,max_v_elements,us_tridim))      !--- CL-OGS:  extension to third dimension 
    ALLOCATE(v_kwele_y(4,max_v_elements,us_tridim))      !--- CL-OGS:  extension to third dimension 
    ALLOCATE(P_r_element(numpar))
    ALLOCATE(P_u_element(numpar))
    ALLOCATE(P_v_element(numpar))
    ALLOCATE(mask_rho(vi,uj,us_tridim))                  !--- CL-OGS:  extension to third dimension    
    ALLOCATE(m_r(vi,uj,us_tridim))                       !--- CL-OGS:  extension to third dimension 
    ALLOCATE(m_u(ui,uj,us_tridim))                       !--- CL-OGS:  extension to third dimension 
    ALLOCATE(m_v(vi,vj,us_tridim))                       !--- CL-OGS:  extension to third dimension 
    ALLOCATE(x_u(ui,uj))                                
    ALLOCATE(y_u(ui,uj))                                
    ALLOCATE(x_v(vi,vj))
    ALLOCATE(y_v(vi,vj))
    ALLOCATE(rho_mask(rho_nodes,us_tridim))  !ewn.v.2     !--- CL-OGS:  extension to third dimension   
    ALLOCATE(u_mask(u_nodes,us_tridim))                   !--- CL-OGS:  extension to third dimension 
    ALLOCATE(v_mask(v_nodes,us_tridim))                   !--- CL-OGS:  extension to third dimension 
    ALLOCATE(coast_eleij(vi-1,uj-1))  
    ALLOCATE(coast_ele(max_rho_elements,us_tridim))  

    !ALLOCATE SUBROUTINE VARIABLES
    if(read_GrainSize)  ALLOCATE(GrainSize_tmp(vi,uj))     !--- CL-OGS:  for behavior 8
    ALLOCATE(euleriandepth(vi,uj))
    ALLOCATE(mask_u(ui,uj,us_tridim))                      !--- CL-OGS:  extension to third dimension 
    ALLOCATE(mask_v(vi,vj,us_tridim))                      !--- CL-OGS:  extension to third dimension 
    ALLOCATE(x_rho(vi,uj))
    ALLOCATE(y_rho(vi,uj))
    ALLOCATE(lon_rho(vi,uj))
    ALLOCATE(lat_rho(vi,uj))
    ALLOCATE(lon_u(ui,uj))
    ALLOCATE(lat_u(ui,uj))
    ALLOCATE(lon_v(vi,vj))
    ALLOCATE(lat_v(vi,vj))
    if(.not.Zgrid) ALLOCATE(angle(vi,uj))            !--- CL-OGS: temporary variable not used for MITgcm 
    ALLOCATE(r_ele(4,max_rho_elements))
    ALLOCATE(u_ele(4,max_u_elements))
    ALLOCATE(v_ele(4,max_v_elements))
    !  ALLOCATE(rho_mask(rho_nodes))    !ewn.v.2
    !  ALLOCATE(u_mask(u_nodes))
    !  ALLOCATE(v_mask(v_nodes))
    !--- CL-OGS: k-level dependant number of wet elements replacing rho_elements
    ALLOCATE (rho_kwele(us_tridim),STAT=STATUS)
    if(STATUS /= 0) write(*,*) 'Problem allocating rho_kwele'
    !--- CL-OGS: k-level dependant number of wet elements replacing u_elements
    ALLOCATE (u_kwele(us_tridim),STAT=STATUS)
    if(STATUS /= 0) write(*,*) 'Problem allocating u_kwele'
    !--- CL-OGS: k-level dependant number of wet elements replacing v_elements
    ALLOCATE (v_kwele(us_tridim),STAT=STATUS)
    if(STATUS /= 0) write(*,*) 'Problem allocating v_kwele'
    ALLOCATE(lastplotpos(numpar,5))
    lastplotpos=0 
    rho_angle=0.0                                          !--- CL-OGS: initalized null for MITgcm
    DeadOrOut(:)=.False.

    WRITE(*,*) 'read-in grid information'

    ! *************************** READ IN GRID INFO **************************

    STATUS = NF90_OPEN(TRIM(NCgridfile),NF90_NOWRITE, NCID)
    if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem NF90_OPEN'
    if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! Depth (m)
      STATUS = NF90_INQ_VARID(NCID,'h',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find depth'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,euleriandepth)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read depth'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! longitude at rho (°)
      STATUS = NF90_INQ_VARID(NCID,'lon_rho',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find lon_rho'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,lon_rho)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read lon_rho'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! latitude at rho (°)
      STATUS = NF90_INQ_VARID(NCID,'lat_rho',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find lat_rho'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,lat_rho)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read lat_rho'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! longitude at u (°)
      STATUS = NF90_INQ_VARID(NCID,'lon_u',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find lon_u'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,lon_u)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read lon_u'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! latitude at u (°)
      STATUS = NF90_INQ_VARID(NCID,'lat_u',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find lat_u'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,lat_u)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read lat_u'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! longitude at v (°)
      STATUS = NF90_INQ_VARID(NCID,'lon_v',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find lon_v'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,lon_v)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read lon_v'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! latitude at v (°)
      STATUS = NF90_INQ_VARID(NCID,'lat_v',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find lat_v'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,lat_v)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read lat_v'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! mask on rho grid
      STATUS = NF90_INQ_VARID(NCID,'mask_rho',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find mask_rho'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,mask_rho)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read mask_rho'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! mask on u grid
      STATUS = NF90_INQ_VARID(NCID,'mask_u',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find mask_u'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,mask_u)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read mask_u'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! mask on v grid
      STATUS = NF90_INQ_VARID(NCID,'mask_v',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find mask_v'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,mask_v)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read mask_v'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      if(Zgrid) then !--- CL-OGS: read MITgcm specific grid files
        if(Vtransform.ne.0)then
          write(*,*)'error in input data file: Vtransform set to ',Vtransform
          write(*,*)'while hydroMITgcm uses fixed Z vertical levels'
          stop
        endif
        ! Z-coordinate on rho grid (Z) : cell-centered coordinates
        STATUS = NF90_INQ_VARID(NCID,'Z',VID)
        IF(STATUS /= NF90_NOERR) write(*,*) 'Problem finding Z'
        if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
        STATUS = NF90_GET_VAR(NCID,VID,ZC)
        if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read Z'
        if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

        ! Z-coordinate on w grid (Zp1) : interface-centered coordinates
        STATUS = NF90_INQ_VARID(NCID,'Zp1',VID)
        IF(STATUS /= NF90_NOERR)write(*,*) 'Problem finding Zp1'
        if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
        STATUS = NF90_GET_VAR(NCID,VID,ZW)
        if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read Zp1'
        if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)


        STATUS = NF90_INQ_VARID(NCID,'KBottomRUV',VID)
        IF(STATUS /= NF90_NOERR)write(*,*) 'Problem finding KBottomRUV'
        if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
        STATUS = NF90_GET_VAR(NCID,VID,BottomK)
        if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read KBottomRUV'
        if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      else          !--- CL-OGS: angle read only for ROMS files  
        if(Vtransform.eq.0)then
          write(*,*)'error in input data file: Vtransform set to ',Vtransform
          write(*,*)'while ROMS uses sigma vertical levels'
          stop
        endif
        ! angle between x-coordinate and true east direction (radian)
        STATUS = NF90_INQ_VARID(NCID,'angle',VID)
        if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find angle'
        if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
        STATUS = NF90_GET_VAR(NCID,VID,angle)
        if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read angle'
        if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      endif

    STATUS = NF90_CLOSE(NCID)

    ! *************************** READ IN GRAIN SIZE FILE ********************
    if(read_GrainSize)then
      WRITE(*,*) 'read-in GrainSize file ',GrainSize_fname
      STATUS = NF90_OPEN(TRIM(GrainSize_fname),NF90_NOWRITE, NCID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem NF90_OPEN'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_INQ_VARID(NCID,'GrainSize',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find GrainSize'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,GrainSize_tmp)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read GrainSize'
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_CLOSE(NCID)
    endif



    ! ************************* READ IN S LEVEL INFO *************************
     if(Zgrid)then
       nfilesin=8 ! MITgcm binary outputs without wind files
       if(Wind) nfilesin=nfilesin+2 
       if(WindIntensity) nfilesin=nfilesin+1 
     else 
       nfilesin=1 ! Roms NETcdf outputs
     endif

     call set_filename(VAR_ID_salt,filenum,filenm)

    if(.not.Zgrid) then
      STATUS = NF90_OPEN(TRIM(filenm), NF90_NOWRITE, NCID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem NF90_OPEN:',TRIM(filenm)
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! s-coordinate on rho grid (sc_r)
      STATUS = NF90_INQ_VARID(NCID,'s_rho',VID)
      IF(STATUS /= NF90_NOERR)THEN
        STATUS = NF90_INQ_VARID(NCID,'sc_r',VID)
        if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem finding SC in',TRIM(filenm)
        if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      ENDIF
      STATUS = NF90_GET_VAR(NCID,VID,SC)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read SC in',TRIM(filenm)
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! Cs value on rho grid (Cs_r)
      STATUS = NF90_INQ_VARID(NCID,'Cs_r',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find CS in',TRIM(filenm)
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,CS)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read CS in',TRIM(filenm)
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! s-coordinate on w grid (sc_w)
      STATUS = NF90_INQ_VARID(NCID,'s_w',VID)
      IF(STATUS /= NF90_NOERR)THEN
        STATUS = NF90_INQ_VARID(NCID,'sc_w',VID)
        if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem finding SCW in',TRIM(filenm)
        if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      ENDIF
      STATUS = NF90_GET_VAR(NCID,VID,SCW)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read SCW in',TRIM(filenm)
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      ! Cs value on w grid (Cs_w)
      STATUS = NF90_INQ_VARID(NCID,'Cs_w',VID)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem find CSW in',TRIM(filenm)
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      STATUS = NF90_GET_VAR(NCID,VID,CSW)
      if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem read CSW in',TRIM(filenm)
      if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

      !close the dataset and reassign the NCID
      STATUS = NF90_CLOSE(NCID)
    endif

    ! *************************** CREATE ELEMENTS *****************************

    !Store Mask Values to multiply by...
    m_r = mask_rho
    m_u = mask_u
    m_v = mask_v

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! ~  4B. Prepare Elements (i.e., assign ID numbers to rectangular grids)  ~
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    write(*,*) 'create coast_elements'

    ! Convert rho nodes lon/lat to x/y coordinates
    do j=1,uj
      do i=1,vi
        x_rho(i,j) = lon2x(lon_rho(i,j),lat_rho(i,j))
        y_rho(i,j) = lat2y(lat_rho(i,j))
      enddo
    enddo

    ! Convert u nodes lon/lat to x/y coordinates
    do j=1,uj
      do i=1,ui
        x_u(i,j) = lon2x(lon_u(i,j),lat_u(i,j))
        y_u(i,j) = lat2y(lat_u(i,j))
      enddo
    enddo

    ! Convert v nodes lon/lat to x/y coordinates
    do j=1,vj
      do i=1,vi
        x_v(i,j) = lon2x(lon_v(i,j),lat_v(i,j))
        y_v(i,j) = lat2y(lat_v(i,j))
      enddo
    enddo

   if(Zgrid)then !--- CL-OGS: read number of the bottom - vertical level 
     count=0
     do j=1,uj
       do i=1,vi
         count = count+1
         BottomKatRnode(count)=Bottomk(i,j,1) 
         if(i.eq.1.or.j.eq.1.or.i.eq.vi.or.j.eq.uj)BottomKatRnode(count)=ws
       enddo
     enddo
     count=0
     do j=1,uj
       do i=1,ui
         count = count+1
         BottomKatUnode(count)=Bottomk(i,j,2) 
         if(i.eq.1.or.j.eq.1.or.i.eq.ui.or.j.eq.uj)BottomKatUnode(count)=ws
       enddo
     enddo
     count=0
     do j=1,vj
       do i=1,vi
         count = count+1
         BottomKatVnode(count)=Bottomk(i,j,3) 
         if(i.eq.1.or.j.eq.1.or.i.eq.vi.or.j.eq.vj)BottomKatVnode(count)=ws
       enddo
     enddo
   endif

    !---------------------------
    ! Assign mask values to rho nodes 
    do k=1,us_tridim  !--- CL-OGS: extention to 3d 
    count = 0
    do j=1,uj
      do i=1,vi
        count = count + 1   !move to next node number
            !cycles through each variable replacing the vi,uj part with count
            !  essentially giving it node numbers
        rho_mask(count,k) = mask_rho(i,j,k)
        if(.not.Zgrid) rho_angle(count) = angle(i,j)
      enddo
    enddo
    enddo

    ! Assign mask values to u nodes
    do k=1,us_tridim  !--- CL-OGS: extention to 3d
    count = 0
    do j=1,uj
      do i=1,ui
        count = count + 1   !move to next node number
            !cycles through each variable replacing the vi,uj part with count
            !  essentially giving it node numbers
        u_mask(count,k) = mask_u(i,j,k)
      enddo
    enddo
    enddo

    ! Assign mask values to v nodes
    do k=1,us_tridim  !--- CL-OGS: extention to 3d
    count = 0
    do j=1,vj
      do i=1,vi
        count = count + 1   !move to next node number
            !cycles through each variable replacing the vi,uj part with count
            !  essentially giving it node numbers
        v_mask(count,k) = mask_v(i,j,k)
      enddo
    enddo
    enddo

    IF(Zgrid)THEN   !--- CL-OGS: prepare parameters to copy hydrodynamic fields
                    ! of rho, u and v water nodes in neighbours land nodes
                    ! ----------- rho nodes : --------------------------------
      nodestocopy=0
      do j=1,uj
       do i=1,vi
        kbot=BottomK(i,j,1)
        kmax=max(1,kbot-2)
        if(kbot.eq.ws)kmax=kbot-1
        tmpcoef(:)=-1
        do k=kmax,1,-1
         if(mask_rho(i,j,k).eq.0)then
          summask=float(mask_rho(max( 1,i-1),j,k)+&
                        mask_rho(min(vi,i+1),j,k)+&
                        mask_rho(i,max( 1,j-1),k)+&
                        mask_rho(i,min(uj,j+1),k) )
          if(summask>0)then
           oldtmpcoef(1)=tmpcoef(1)
           oldtmpcoef(2)=tmpcoef(2)
           oldtmpcoef(3)=tmpcoef(3)
           oldtmpcoef(4)=tmpcoef(4)          
           tmpcoef(1)=float(mask_rho(max( 1,i-1),j,k))/summask
           tmpcoef(2)=float(mask_rho(min(vi,i+1),j,k))/summask
           tmpcoef(3)=float(mask_rho(i,max( 1,j-1),k))/summask
           tmpcoef(4)=float(mask_rho(i,min(uj,j+1),k))/summask
           if((k.eq.kmax) .or. &
             (oldtmpcoef(1).ne.tmpcoef(1)).or.&
             (oldtmpcoef(2).ne.tmpcoef(2)).or.&
             (oldtmpcoef(3).ne.tmpcoef(3)).or.&
             (oldtmpcoef(4).ne.tmpcoef(4)))then
               nodestocopy=nodestocopy+1
           elseif(k.ne.kmax)then
               continue
           else
            exit
           endif
          endif
         endif
        enddo
       enddo
      enddo
      write(*,*)'hydro rho node arrays to copy is ',nodestocopy
      maxnodestocopy=nodestocopy
      NUM_COPNOD(RNODE)=nodestocopy
                    ! ----------- u nodes : --------------------------------
      nodestocopy=0
      do j=1,uj
       do i=1,ui
        kbot=BottomK(i,j,2)
        kmax=max(1,kbot-2)
        if(kbot.eq.ws)kmax=kbot-1
        tmpcoef(:)=-1
        do k=kmax,1,-1
         if(mask_u(i,j,k).eq.0)then
          summask=float(mask_u(max( 1,i-1),j,k)+&
                        mask_u(min(ui,i+1),j,k)+&
                        mask_u(i,max( 1,j-1),k)+&
                        mask_u(i,min(uj,j+1),k) )
          if(summask>0)then
           oldtmpcoef(1)=tmpcoef(1)
           oldtmpcoef(2)=tmpcoef(2)
           oldtmpcoef(3)=tmpcoef(3)
           oldtmpcoef(4)=tmpcoef(4)          
           tmpcoef(1)=float(mask_u(max( 1,i-1),j,k))/summask
           tmpcoef(2)=float(mask_u(min(ui,i+1),j,k))/summask
           tmpcoef(3)=float(mask_u(i,max( 1,j-1),k))/summask
           tmpcoef(4)=float(mask_u(i,min(uj,j+1),k))/summask
           if((k.eq.kmax) .or. &
             (oldtmpcoef(1).ne.tmpcoef(1)).or.&
             (oldtmpcoef(2).ne.tmpcoef(2)).or.&
             (oldtmpcoef(3).ne.tmpcoef(3)).or.&
             (oldtmpcoef(4).ne.tmpcoef(4)))then
               nodestocopy=nodestocopy+1
           elseif(k.ne.kmax)then
                continue
           else
            exit
           endif
          endif
         endif
        enddo
       enddo
      enddo
      maxnodestocopy=max(maxnodestocopy,nodestocopy)
      write(*,*)'u nodes to copy is ',nodestocopy
      NUM_COPNOD(UNODE)=nodestocopy
                    ! -----------  v nodes : --------------------------------
      nodestocopy=0
      do j=1,vj
       do i=1,vi
        kbot=BottomK(i,j,3)
        kmax=max(1,kbot-2)
        if(kbot.eq.ws)kmax=kbot-1
        tmpcoef(:)=-1
        do k=kmax,1,-1
         if(mask_v(i,j,k).eq.0)then
          summask=float(mask_v(max( 1,i-1),j,k)+&
                        mask_v(min(vi,i+1),j,k)+&
                        mask_v(i,max( 1,j-1),k)+&
                        mask_v(i,min(vj,j+1),k) )
          if(summask>0)then
           oldtmpcoef(1)=tmpcoef(1)
           oldtmpcoef(2)=tmpcoef(2)
           oldtmpcoef(3)=tmpcoef(3)
           oldtmpcoef(4)=tmpcoef(4)          
           tmpcoef(1)=float(mask_v(max( 1,i-1),j,k))/summask
           tmpcoef(2)=float(mask_v(min(vi,i+1),j,k))/summask
           tmpcoef(3)=float(mask_v(i,max( 1,j-1),k))/summask
           tmpcoef(4)=float(mask_v(i,min(vj,j+1),k))/summask
           if((k.eq.kmax) .or. &
             (oldtmpcoef(1).ne.tmpcoef(1)).or.&
             (oldtmpcoef(2).ne.tmpcoef(2)).or.&
             (oldtmpcoef(3).ne.tmpcoef(3)).or.&
             (oldtmpcoef(4).ne.tmpcoef(4)))then
               nodestocopy=nodestocopy+1
           elseif(k.ne.kmax)then
                continue
           else
            exit
           endif
          endif
         endif
        enddo
       enddo
      enddo
      maxnodestocopy=max(maxnodestocopy,nodestocopy)
      write(*,*)'v nodes to copy is ',nodestocopy
      NUM_COPNOD(VNODE)=nodestocopy

      !-----------------------------------------------------------------------------
      ! ALLOCATE MATRIX STORING COPNOD parameters:  ! ------------------------------
      ALLOCATE(Coef_COPNOD(3,maxnodestocopy,4))
      ALLOCATE(Nghb_COPNOD(3,maxnodestocopy,4))
      ALLOCATE(Node_COPNOD(3,maxnodestocopy))
      ALLOCATE(klev_COPNOD(3,maxnodestocopy,2))

      !---------------------------------------------
      nodestocopy=0
      count=0
      old_i=0
      old_j=0
      old_count=0    
      do j=1,uj
       do i=1,vi
        count = count + 1   !move to next node number
        kbot=BottomK(i,j,1)
        kmax=max(1,kbot-2)
        if(kbot.eq.ws)kmax=kbot-1
        do k=kmax,1,-1
         if(mask_rho(i,j,k).eq.0)then
          summask=float(mask_rho(max( 1,i-1),j,k)+&
                        mask_rho(min(vi,i+1),j,k)+&
                        mask_rho(i,max( 1,j-1),k)+&
                        mask_rho(i,min(uj,j+1),k) )
          if(summask>0)then
           tmpcoef(1)=float(mask_rho(max( 1,i-1),j,k))/summask
           tmpcoef(2)=float(mask_rho(min(vi,i+1),j,k))/summask
           tmpcoef(3)=float(mask_rho(i,max( 1,j-1),k))/summask
           tmpcoef(4)=float(mask_rho(i,min(uj,j+1),k))/summask
           if((k.eq.kmax) .or. &
             (Coef_COPNOD(RNODE,max(1,nodestocopy),1).ne.tmpcoef(1)).or.&
             (Coef_COPNOD(RNODE,max(1,nodestocopy),2).ne.tmpcoef(2)).or.&
             (Coef_COPNOD(RNODE,max(1,nodestocopy),3).ne.tmpcoef(3)).or.&
             (Coef_COPNOD(RNODE,max(1,nodestocopy),4).ne.tmpcoef(4)))then
               nodestocopy=nodestocopy+1
               old_i=i
               old_j=j
               old_count=count
               Coef_COPNOD(RNODE,nodestocopy,1)=tmpcoef(1)
               Coef_COPNOD(RNODE,nodestocopy,2)=tmpcoef(2)
               Coef_COPNOD(RNODE,nodestocopy,3)=tmpcoef(3)
               Coef_COPNOD(RNODE,nodestocopy,4)=tmpcoef(4)
               Nghb_COPNOD(RNODE,nodestocopy,1)=max(1,count-1)
               Nghb_COPNOD(RNODE,nodestocopy,2)=min(uj*vi,count+1)
               Nghb_COPNOD(RNODE,nodestocopy,3)=max(1,count-vi)
               Nghb_COPNOD(RNODE,nodestocopy,4)=min(uj*vi,count+vi)
               Node_COPNOD(RNODE,nodestocopy)=count
               klev_COPNOD(RNODE,nodestocopy,1)=k
               klev_COPNOD(RNODE,nodestocopy,2)=k
           elseif(k.ne.kmax)then
                klev_COPNOD(RNODE,nodestocopy,2)=k
           else
            exit
           endif
          endif
         endif
        enddo
       enddo
      enddo
      !----------------------------------------------------------
      nodestocopy=0
      count=0
      do j=1,uj
       do i=1,ui
        count = count + 1   !move to next node number
        kbot=BottomK(i,j,2)
        kmax=max(1,kbot-2)
        if(kbot.eq.ws)kmax=kbot-1
        do k=kmax,1,-1
         if(mask_u(i,j,k).eq.0)then
          summask=float(mask_u(max( 1,i-1),j,k)+&
                        mask_u(min(ui,i+1),j,k)+&
                        mask_u(i,max( 1,j-1),k)+&
                        mask_u(i,min(uj,j+1),k) )
          if(summask>0)then
           tmpcoef(1)=float(mask_u(max( 1,i-1),j,k))/summask
           tmpcoef(2)=float(mask_u(min(ui,i+1),j,k))/summask
           tmpcoef(3)=float(mask_u(i,max( 1,j-1),k))/summask
           tmpcoef(4)=float(mask_u(i,min(uj,j+1),k))/summask
           if((k.eq.kmax) .or. &
             (Coef_COPNOD(UNODE,max(1,nodestocopy),1).ne.tmpcoef(1)).or.&
             (Coef_COPNOD(UNODE,max(1,nodestocopy),2).ne.tmpcoef(2)).or.&
             (Coef_COPNOD(UNODE,max(1,nodestocopy),3).ne.tmpcoef(3)).or.&
             (Coef_COPNOD(UNODE,max(1,nodestocopy),4).ne.tmpcoef(4)))then
               nodestocopy=nodestocopy+1
               Coef_COPNOD(UNODE,nodestocopy,1)=tmpcoef(1)
               Coef_COPNOD(UNODE,nodestocopy,2)=tmpcoef(2)
               Coef_COPNOD(UNODE,nodestocopy,3)=tmpcoef(3)
               Coef_COPNOD(UNODE,nodestocopy,4)=tmpcoef(4)
               Nghb_COPNOD(UNODE,nodestocopy,1)=max(1,count-1)
               Nghb_COPNOD(UNODE,nodestocopy,2)=min(uj*ui,count+1)
               Nghb_COPNOD(UNODE,nodestocopy,3)=max(1,count-ui)
               Nghb_COPNOD(UNODE,nodestocopy,4)=min(uj*ui,count+ui)
               Node_COPNOD(UNODE,nodestocopy)=count
               klev_COPNOD(UNODE,nodestocopy,1)=k
               klev_COPNOD(UNODE,nodestocopy,2)=k
           elseif(k.ne.kmax)then  ! if(k.ne.kmax ): could be cancelled?
                klev_COPNOD(UNODE,nodestocopy,2)=k
           else   ! could be cancelled?
            exit  ! could be cancelled?
           endif
          endif
         endif
        enddo
       enddo
      enddo
      !-------------------------------------------------------
      nodestocopy=0
      count=0
      do j=1,vj
       do i=1,vi
        count = count + 1   !move to next node number
        kbot=BottomK(i,j,3)
        kmax=max(1,kbot-2)
        if(kbot.eq.ws)kmax=kbot-1
        do k=kmax,1,-1
         if(mask_v(i,j,k).eq.0)then
          summask=float(mask_v(max( 1,i-1),j,k)+&
                        mask_v(min(vi,i+1),j,k)+&
                        mask_v(i,max( 1,j-1),k)+&
                        mask_v(i,min(vj,j+1),k) )
          if(summask>0)then
           tmpcoef(1)=float(mask_v(max( 1,i-1),j,k))/summask
           tmpcoef(2)=float(mask_v(min(vi,i+1),j,k))/summask
           tmpcoef(3)=float(mask_v(i,max( 1,j-1),k))/summask
           tmpcoef(4)=float(mask_v(i,min(vj,j+1),k))/summask
           if((k.eq.kmax) .or. &
             (Coef_COPNOD(VNODE,max(1,nodestocopy),1).ne.tmpcoef(1)).or.&
             (Coef_COPNOD(VNODE,max(1,nodestocopy),2).ne.tmpcoef(2)).or.&
             (Coef_COPNOD(VNODE,max(1,nodestocopy),3).ne.tmpcoef(3)).or.&
             (Coef_COPNOD(VNODE,max(1,nodestocopy),4).ne.tmpcoef(4)))then
               nodestocopy=nodestocopy+1
               Coef_COPNOD(VNODE,nodestocopy,1)=tmpcoef(1)
               Coef_COPNOD(VNODE,nodestocopy,2)=tmpcoef(2)
               Coef_COPNOD(VNODE,nodestocopy,3)=tmpcoef(3)
               Coef_COPNOD(VNODE,nodestocopy,4)=tmpcoef(4)
               Nghb_COPNOD(VNODE,nodestocopy,1)=max(1,count-1)
               Nghb_COPNOD(VNODE,nodestocopy,2)=min(vj*vi,count+1)
               Nghb_COPNOD(VNODE,nodestocopy,3)=max(1,count-vi)
               Nghb_COPNOD(VNODE,nodestocopy,4)=min(vj*vi,count+vi)
               Node_COPNOD(VNODE,nodestocopy)=count
               klev_COPNOD(VNODE,nodestocopy,1)=k
               klev_COPNOD(VNODE,nodestocopy,2)=k
           elseif(k.ne.kmax)then
                klev_COPNOD(VNODE,nodestocopy,2)=k
           else
            exit
           endif
          endif
         endif
        enddo
       enddo
      enddo
    ENDIF  ! (Zgrid) -------- end preparing parameters to copy hydrodynamic fields  -------

    ! Create matrix that contains the node numbers for each rho element
    !  n(4,count)-------n(3,count)       ^ j
    !       |     count     |            |   
    !  n(1,count)-------n(2,count)       |--> i
    ! 
    count = 0
    do j=1,uj-1                         !z2v3.2
      do i=1,vi-1
        count = count + 1
        r_ele(1,count) = i + (j-1)*vi
        r_ele(2,count) = i + 1 + (j-1)*vi
        r_ele(3,count) = i + 1 + j*vi
        r_ele(4,count) = i + j*vi
      enddo
    enddo

    ! Create matrix that contains the node numbers for each u element
    count = 0
    do j=1,uj-1                         !z2v3.2
      do i=1,ui-1
        count = count + 1
        u_ele(1,count) = i + (j-1)*ui
        u_ele(2,count) = i + 1 + (j-1)*ui
        u_ele(3,count) = i + 1 + j*ui
        u_ele(4,count) = i + j*ui
      enddo
    enddo

    ! Create matrix that contains the node numbers for each v element
    count = 0
    do j=1,vj-1                         !z2v3.2
      do i=1,vi-1
        count = count + 1
        v_ele(1,count) = i + (j-1)*vi
        v_ele(2,count) = i + 1 + (j-1)*vi
        v_ele(3,count) = i + 1 + j*vi
        v_ele(4,count) = i + j*vi
      enddo
    enddo

    !--------------------------------------------------------------------------
    !--- CL-OGS: assign number to elements to identify their distance from the coast in elements
    do i=1,vi-1
      do j=1,uj-1
        coast_eleij(i,j)=10 ! water ele
        if(mask_rho(i  ,j+1,us_tridim) == 0) coast_eleij(i,j)=1 ! land/bnd ele
        if(mask_rho(i  ,j  ,us_tridim) == 0) coast_eleij(i,j)=1 ! land/bnd ele
        if(mask_rho(i+1,j+1,us_tridim) == 0) coast_eleij(i,j)=1 ! land/bnd ele
        if(mask_rho(i+1,j  ,us_tridim) == 0) coast_eleij(i,j)=1 ! land/bnd ele                                      
      enddo
    enddo
    ! Identify second water coast_elements which are the ones directly adjacent to boundary coast_elements
    do i=2,vi-2
      do j=2,uj-2
         if(coast_eleij(i,j)==10) then
            if(  ( coast_eleij(i  ,j+1)<=1 ) &
            .or. ( coast_eleij(i  ,j-1)<=1 ) &
            .or. ( coast_eleij(i+1,j  )<=1 ) &
            .or. ( coast_eleij(i-1,j  )<=1 ) ) coast_eleij(i,j)=2
         endif
      enddo
    enddo
    ! Identify third water coast_elements which are the ones directly adjacent to the second water coast_elements
    do i=2,vi-2
      do j=2,uj-2
         if(coast_eleij(i,j)==10) then
            if(  ( coast_eleij(i  ,j+1)<=2 ) &
            .or. ( coast_eleij(i  ,j-1)<=2 ) &
            .or. ( coast_eleij(i+1,j  )<=2 ) &
            .or. ( coast_eleij(i-1,j  )<=2 ) ) coast_eleij(i,j)=3
         endif
      enddo
    enddo
    ! Identify fourth water coast_elements which are the ones directly adjacent to the third water coast_elements
    do i=2,vi-2
      do j=2,uj-2
         if(coast_eleij(i,j)==10) then
            if(  ( coast_eleij(i  ,j+1)<=3 ) &
            .or. ( coast_eleij(i  ,j-1)<=3 ) &
            .or. ( coast_eleij(i+1,j  )<=3 ) &
            .or. ( coast_eleij(i-1,j  )<=3 ) ) coast_eleij(i,j)=4
         endif
      enddo
    enddo
    ! Identify fifth water coast_elements which are the ones directly adjacent to the fourth water coast_elements
    do i=2,vi-2
      do j=2,uj-2
         if(coast_eleij(i,j)==10) then
            if(  ( coast_eleij(i  ,j+1)<=4 ) &
            .or. ( coast_eleij(i  ,j-1)<=4 ) &
            .or. ( coast_eleij(i+1,j  )<=4 ) &
            .or. ( coast_eleij(i-1,j  )<=4 ) ) coast_eleij(i,j)=5
         endif
      enddo
    enddo

    do k=1,us_tridim
     count = 0
     do j=1,uj-1                         !z2v3.2
       do i=1,vi-1
         inele = 0
         !using the mask determine if any of the nodes for the current
         !  element are inbounds, if so set inele to 1
         if(mask_rho(i  ,j+1,k) == 1)  inele=1
         if(mask_rho(i  ,j  ,k) == 1)  inele=1
         if(mask_rho(i+1,j+1,k) == 1)  inele=1
         if(mask_rho(i+1,j  ,k) == 1)  inele=1              !z2v3.2
         !if inele = 1 then at least one of the three nodes for this element
         !  are in bounds.
         if( inele .EQ. 1 ) then
           count = count + 1
           coast_ele(count,k) = coast_eleij(i,j) ! counting only water elements!!!
         endif
       enddo
     enddo
    enddo
   

    ! ------------------------------------------------------------------------

    ! Create matrix that contains only the rho elements that contain a node
    ! whose mask value = 1 (i.e., it has at least one water point).
    RE=0   !--- CL-OGS: initialize because the full RE matrix won't be filled by the following
           ! instructions, as dim(RE(,:,)=max_rho_elements>rho_kwele(k) 
    rho_kwele(:)=0
    do k=1,us_tridim  !--- CL-OGS: extention to 3d
    count = 0
    do i=1,max_rho_elements
      inele = 0
      !using the mask determine if any of the nodes for the current
      !  element are inbounds, if so set inele to 1
      if( rho_mask(r_ele(1,i),k) .EQ. 1) inele=1
      if( rho_mask(r_ele(2,i),k) .EQ. 1) inele=1
      if( rho_mask(r_ele(3,i),k) .EQ. 1) inele=1
      if( rho_mask(r_ele(4,i),k) .EQ. 1) inele=1              !z2v3.2
      !if inele = 1 then at least one of the three nodes for this element
      !  are in bounds.
      if( inele .EQ. 1 ) then
        count = count + 1
        !create array of elements that contain at least one node in bounds
        RE(1,count,k) = r_ele(1,i)
        RE(2,count,k) = r_ele(2,i)
        RE(3,count,k) = r_ele(3,i)
        RE(4,count,k) = r_ele(4,i)                            !z2v3.2
      endif
    enddo
    rho_kwele(k) = count !--- CL-OGS: number of rho elements with at least
                         !--- CL-OGS  one corner with mask_rho=1 at level k
    enddo
    
    ! Create matrix that contains only the u elements that contain a node
    ! whose mask value = 1 (i.e., it has at least one water point).
    UE=0 ! initialize because the full UE matrix won't be filled by the following
         ! instructions, as dim(UE(,:,)=max_u_elements>u_kwele(k)
    u_kwele(:) = 0
    do k=1,us_tridim
    count = 0
    countele=0
    do j=1,uj-1                         !z2v3.2
     do i=1,ui-1
      countele=countele+1
      !do i=1,max_u_elements
      inele = 0
      !using the mask determine if any of the nodes for the current
      !  element are inbounds, if so set inele to 1
      if( u_mask(u_ele(1,countele),k) .EQ. 1) inele=1
      if( u_mask(u_ele(2,countele),k) .EQ. 1) inele=1
      if( u_mask(u_ele(3,countele),k) .EQ. 1) inele=1
      if( u_mask(u_ele(4,countele),k) .EQ. 1) inele=1                    !z2v3.2
      if( mask_rho(i,j,k)==0 .and. mask_rho(i+1,j,k)==1                        &
      .and. mask_rho(i+2,j,k)==0 ) inele=1 ! for 1-cell-width vertical channels
      !if inele = 1 then at least one of the three nodes for this element
      !  are in bounds.
      if( inele .EQ. 1 ) then
        count = count + 1
        !create array of elements that contain at least one node in bounds
        UE(1,count,k) = u_ele(1,countele)
        UE(2,count,k) = u_ele(2,countele)
        UE(3,count,k) = u_ele(3,countele)
        UE(4,count,k) = u_ele(4,countele)                            !z2v3.2
      endif
     enddo 
    enddo
    u_kwele(k) = count !  u  surface elements with at least one corner with mask_u=1
    enddo

    ! Create matrix that contains only the v elements that contain a node
    ! whose mask value = 1 (i.e., it has at least one water point).
    VE=0 ! initialize because the full VE matrix won't be filled by the following
         ! instructions, as dim(VE(,:,)=max_v_elements>v_kwele(k) 
    v_kwele(:) = 0
    do k=1,us_tridim
    count = 0
    countele=0
    do j=1,vj-1                         !z2v3.2
     do i=1,vi-1
      countele=countele+1
      !do i=1,max_v_elements
      inele = 0
      !using the mask determine if any of the nodes for the current
      !  element are inbounds, if so set inele to 1
      if( v_mask(v_ele(1,countele),k) .EQ. 1) inele=1
      if( v_mask(v_ele(2,countele),k) .EQ. 1) inele=1
      if( v_mask(v_ele(3,countele),k) .EQ. 1) inele=1
      if( v_mask(v_ele(4,countele),k) .EQ. 1) inele=1                    !z2v3.2
      if( mask_rho(i,j,k)==0 .and. mask_rho(i,j+1,k)==1                        &
        .and. mask_rho(i,j+2,k)==0) inele=1 ! for 1-cell-heigh horizontal channels
      !if inele = 1 then at least one of the three nodes for this element
      !  are in bounds.
      if( inele .EQ. 1 ) then
        count = count + 1
        !create array of elements that contain at least one node in bounds
        VE(1,count,k) = v_ele(1,countele)
        VE(2,count,k) = v_ele(2,countele)
        VE(3,count,k) = v_ele(3,countele)
        VE(4,count,k) = v_ele(4,countele)                            !z2v3.2
      endif 
     enddo
    enddo
    v_kwele(k) = count !  v  surface elements with at least one corner with mask_v=1
    enddo

    ! Create matrices of  x/y for rho nodes and depth values in rho node number
    !   format 
    count = 0
    do j=1,uj
      do i=1,vi
        count = count + 1   !move to next node number
        !cycles through each variable replacing the vi,uj part with count
        !  essentially giving it node numbers
        rx(count) = x_rho(i,j)
        ry(count) = y_rho(i,j)
        depth(count) = euleriandepth(i,j)
      enddo
    enddo
   
    ! Create matrices of GrainSize in rho node number format
    if(read_GrainSize)then
     count = 0
     do j=1,uj
      do i=1,vi
        count = count + 1   !move to next node number
        GrainSize(count) = GrainSize_tmp(i,j)
      enddo
     enddo
    endif

    ! Create matrices of x/y values for u nodes in u node number format 
    count = 0
    do j=1,uj
      do i=1,ui
        count = count + 1   !move to next node number
        !cycles through each variable replacing the ui,uj part with count
        !  essentially giving it node numbers
        ux(count) = x_u(i,j)
        uy(count) = y_u(i,j)
      enddo
    enddo

    ! Create matrices of x/y values for v nodes in v node number format
    count = 0
    do j=1,vj
      do i=1,vi
        count = count + 1   !move to next node number
        !cycles through each variable replacing the vi,vj part with count
        !  essentially giving it node numbers
        vx(count) = x_v(i,j)
        vy(count) = y_v(i,j)
      enddo
    enddo

    ! Create matrices of x/y node values for each rho, u, and v water element  
    r_kwele_x=0.0
    r_kwele_y=0.0
    do k=1,us_tridim                      !--- CL-OGS: extention to 3d
    do j=1,rho_kwele(k)                   !--- CL-OGS: extention to 3d
      do i=1,4                                          !z2v3.3
        r_kwele_x(i,j,k) = rx(RE(i,j,k))  !--- CL-OGS: extention to 3d
        r_kwele_y(i,j,k) = ry(RE(i,j,k))  !--- CL-OGS: extention to 3d
      enddo
    enddo
    enddo

    if(BoundaryBLNs) then  !--- CL-OGS: write rho_kwele and rho bottom level and mask in csv file
      OPEN(110,FILE='rho_kwele.csv',POSITION='APPEND',status='replace')
      do k=1,us_tridim
          do j=1,rho_kwele(k)
            write(110,"(2(F10.5,','),2(i7,','))") &
                          x2lon(0.25*(r_kwele_x(1,j,k)+r_kwele_x(2,j,k)        &
                                  +r_kwele_x(3,j,k)+r_kwele_x(4,j,k)),         &
                          0.25*(r_kwele_y(1,j,k)+r_kwele_y(2,j,k)+             &
                                  r_kwele_y(3,j,k)+r_kwele_y(4,j,k))),         &
                          y2lat(0.25*(r_kwele_y(1,j,k)+r_kwele_y(2,j,k)+       &
                            r_kwele_y(3,j,k)+r_kwele_y(4,j,k))),k,RE(1,j,k)
          enddo
      enddo
      CLOSE(110)
      OPEN(110,FILE='water_rho_nodes.csv',POSITION='APPEND',status='replace')
      CLOSE(110)
    endif

    u_kwele_x=0.0
    u_kwele_y=0.0
    do k=1,us_tridim                     !--- CL-OGS: extention to 3d
    do j=1,u_kwele(k)                    !--- CL-OGS: extention to 3d
      do i=1,4                                          !z2v3.3
        u_kwele_x(i,j,k) = ux(UE(i,j,k)) !--- CL-OGS: extention to 3d
        u_kwele_y(i,j,k) = uy(UE(i,j,k)) !--- CL-OGS: extention to 3d
      enddo
    enddo
    enddo
 
    if(BoundaryBLNs) then  !--- CL-OGS: write u_kwele and u bottom level and mask in csv file
      k=us_tridim
      OPEN(110,FILE='u_kwele.csv',POSITION='APPEND',status='replace')
          do j=1,u_kwele(k)
            write(110,"(2(F10.5,','))") &
                          x2lon(0.25*(u_kwele_x(1,j,k)+u_kwele_x(2,j,k)+       &
                                   u_kwele_x(3,j,k)+u_kwele_x(4,j,k)),         &
                          0.25*(u_kwele_y(1,j,k)+u_kwele_y(2,j,k)+             &
                                   u_kwele_y(3,j,k)+u_kwele_y(4,j,k))), &
                          y2lat(0.25*(u_kwele_y(1,j,k)+u_kwele_y(2,j,k)+       &
                                   u_kwele_y(3,j,k)+u_kwele_y(4,j,k)))
          enddo
      CLOSE(110)  
    endif 
    v_kwele_x=0.0
    v_kwele_y=0.0
    do k=1,us_tridim                     !--- CL-OGS: extention to 3d
    do j=1,v_kwele(k)                    !--- CL-OGS: extention to 3d
      do i=1,4                                          !z2v3.3
        v_kwele_x(i,j,k) = vx(VE(i,j,k)) !--- CL-OGS: extention to 3d
        v_kwele_y(i,j,k) = vy(VE(i,j,k)) !--- CL-OGS: extention to 3d
      enddo
    enddo
    enddo
 
    if(BoundaryBLNs) then  !--- CL-OGS: write v_kwele and v bottom level and mask in csv file
      k=us_tridim
      OPEN(110,FILE='v_kwele.csv',POSITION='APPEND',status='replace')
          do j=1,v_kwele(k)
            write(110,"(2(F10.5,','))") &
                          x2lon(0.25*(v_kwele_x(1,j,k)+v_kwele_x(2,j,k)+       &
                                   v_kwele_x(3,j,k)+v_kwele_x(4,j,k)),         &
                          0.25*(v_kwele_y(1,j,k)+v_kwele_y(2,j,k)+             &
                                    v_kwele_y(3,j,k)+v_kwele_y(4,j,k))),       &
                          y2lat(0.25*(v_kwele_y(1,j,k)+v_kwele_y(2,j,k)+       &
                                    v_kwele_y(3,j,k)+v_kwele_y(4,j,k)))
          enddo
      CLOSE(110)   
    endif

    ! ************************ FIND ADJACENT ELEMENTS *************************

    ! Create search restriction algorithms
    !if(ADJele_file)then
      write(*,*) 'find adjacent elements, reading or writing in file ',        &
                        TRIM(ADJele_fname)
      open(unit=110,file=TRIM(ADJele_fname),form='unformatted',                &
                 status='unknown', action='readwrite',access='direct',         &
                 recl=4*max_rho_elements*10*us_tridim, iostat=ios)
      if ( ios /= 0 ) then
        do waiting=1,10
         call sleep(30)
         write(*,*)'waiting as opening of file',trim(ADJele_fname),' failed'
         open(unit=110,file=TRIM(ADJele_fname),form='unformatted',             &
               status='unknown',  action='readwrite',access='direct',          &
               recl=4*max_rho_elements*10*us_tridim,iostat=ios)
         if ( ios == 0 ) exit
        enddo
      endif
  
      if ( ios /= 0 ) then
        write(*,*) " ERROR OPENING ",TRIM(ADJele_fname)
        stop
      endif
    !endif
    if(ADJele_file)then
      write(*,*)'read from file ',TRIM(ADJele_fname)
      read(110,rec=1,IOSTAT=ios)r_Adjacent(:,:,:)
      if ( ios /= 0 ) stop " ERROR reading r_Adjacent "
      read(110,rec=2,IOSTAT=ios)u_Adjacent(:,:,:)
      if ( ios /= 0 ) stop " ERROR reading u_Adjacent "
      read(110,rec=3,IOSTAT=ios)v_Adjacent(:,:,:)
      if ( ios /= 0 ) stop " ERROR reading v_Adjacent "
    else
      ! I. For each element, list all elements that are adjacent to it 
      write(*,*) ' - compute rho adjacent elements '
      r_Adjacent=0
      do k=1,us_tridim
        do i=1,rho_kwele(k)
          r_Adjacent(i,1,k) = i
          m=1
          do j=max(i-(vi+2),1),min(i+vi+2,rho_kwele(k))
            if(j.EQ.i) cycle
            if(  (RE(1,i,k).EQ.RE(1,j,k)) .OR. (RE(1,i,k).EQ.RE(2,j,k))        &
            .OR. (RE(1,i,k).EQ.RE(3,j,k)) .OR. (RE(1,i,k).EQ.RE(4,j,k))        &
            .OR. (RE(2,i,k).EQ.RE(1,j,k)) .OR. (RE(2,i,k).EQ.RE(2,j,k))        &
            .OR. (RE(2,i,k).EQ.RE(3,j,k)) .OR. (RE(2,i,k).EQ.RE(4,j,k))        &
            .OR. (RE(3,i,k).EQ.RE(1,j,k)) .OR. (RE(3,i,k).EQ.RE(2,j,k))        &
            .OR. (RE(3,i,k).EQ.RE(3,j,k)) .OR. (RE(3,i,k).EQ.RE(4,j,k))        &
            .OR. (RE(4,i,k).EQ.RE(1,j,k)) .OR. (RE(4,i,k).EQ.RE(2,j,k))        &
            .OR. (RE(4,i,k).EQ.RE(3,j,k)) .OR. (RE(4,i,k).EQ.RE(4,j,k)) )then
              m=m+1
              r_Adjacent(i,m,k) = j
            endif
          enddo 
        enddo 
      enddo 
     
      write(*,*) ' - compute u adjacent elements '
      u_Adjacent=0
      do k=1,us_tridim
        do i=1,u_kwele(k)
          u_Adjacent(i,1,k) = i
          m=1
          do j=max(i-(ui+2),1),min(i+ui+2,u_kwele(k))
            if(j.EQ.i) cycle
            if(  (UE(1,i,k).EQ.UE(1,j,k)) .OR. (UE(1,i,k).EQ.UE(2,j,k))        &
            .OR. (UE(1,i,k).EQ.UE(3,j,k)) .OR. (UE(1,i,k).EQ.UE(4,j,k))        &
            .OR. (UE(2,i,k).EQ.UE(1,j,k)) .OR. (UE(2,i,k).EQ.UE(2,j,k))        &
            .OR. (UE(2,i,k).EQ.UE(3,j,k)) .OR. (UE(2,i,k).EQ.UE(4,j,k))        &
            .OR. (UE(3,i,k).EQ.UE(1,j,k)) .OR. (UE(3,i,k).EQ.UE(2,j,k))        &
            .OR. (UE(3,i,k).EQ.UE(3,j,k)) .OR. (UE(3,i,k).EQ.UE(4,j,k))        &
            .OR. (UE(4,i,k).EQ.UE(1,j,k)) .OR. (UE(4,i,k).EQ.UE(2,j,k))        &
            .OR. (UE(4,i,k).EQ.UE(3,j,k)) .OR. (UE(4,i,k).EQ.UE(4,j,k)) ) then
              m=m+1
              u_Adjacent(i,m,k) = j
            endif
          enddo 
        enddo 
      enddo 
     
      write(*,*) ' - compute v adjacent elements '
      v_Adjacent=0
      do k=1,us_tridim
        do i=1,v_kwele(k)
          v_Adjacent(i,1,k) = i
          m=1
          do j=max(i-(vi+2),1),min(i+vi+2,v_kwele(k))
            if(j.EQ.i) cycle
            if(  (VE(1,i,k).EQ.VE(1,j,k)) .OR. (VE(1,i,k).EQ.VE(2,j,k))        &
            .OR. (VE(1,i,k).EQ.VE(3,j,k)) .OR. (VE(1,i,k).EQ.VE(4,j,k))        &
            .OR. (VE(2,i,k).EQ.VE(1,j,k)) .OR. (VE(2,i,k).EQ.VE(2,j,k))        &
            .OR. (VE(2,i,k).EQ.VE(3,j,k)) .OR. (VE(2,i,k).EQ.VE(4,j,k))        &
            .OR. (VE(3,i,k).EQ.VE(1,j,k)) .OR. (VE(3,i,k).EQ.VE(2,j,k))        &
            .OR. (VE(3,i,k).EQ.VE(3,j,k)) .OR. (VE(3,i,k).EQ.VE(4,j,k))        &
            .OR. (VE(4,i,k).EQ.VE(1,j,k)) .OR. (VE(4,i,k).EQ.VE(2,j,k))        &
            .OR. (VE(4,i,k).EQ.VE(3,j,k)) .OR. (VE(4,i,k).EQ.VE(4,j,k))  )then     
              m=m+1
              v_Adjacent(i,m,k) = j
            endif
          enddo 
        enddo
      enddo 

     write(*,*)'write in rho,u and v adjacent elements in file ',              &
                                                  TRIM(ADJele_fname)
     write(110,rec=1,IOSTAT=ios)r_Adjacent(:,:,:)
     if ( ios /= 0 ) stop " ERROR writing r_Adjacent "
     write(110,rec=2,IOSTAT=ios)u_Adjacent(:,:,:)
     if ( ios /= 0 ) stop " ERROR writing u_Adjacent "
     write(110,rec=3,IOSTAT=ios)v_Adjacent(:,:,:)
     if ( ios /= 0 ) stop " ERROR writing v_Adjacent "
    endif 
    CLOSE(110)
    GRD_SET = .TRUE.

    !DEALLOCATE SUBROUTINE VARIABLES
    DEALLOCATE(euleriandepth,mask_u,mask_v,x_rho,y_rho)
    if(read_GrainSize) DEALLOCATE(GrainSize_tmp) 
    if(.not.Zgrid) DEALLOCATE(angle)
    ! DEALLOCATE(r_ele,u_ele,v_ele,rho_mask,u_mask,v_mask)
    DEALLOCATE(r_ele,u_ele,v_ele)

    fpy=123456
    open (unit=fpy,file=TRIM(OutDir)//'/'//TRIM(NCOutFile)//'PartinEle.py',    &
                 STATUS='REPLACE')
    write(fpy,'(a)')'import matplotlib.pyplot as plt'
    write(fpy,'(a)')'import numpy as np'
    write(fpy,'(a)')'fig,ax=plt.subplots()'
    close(fpy)
    !ff=123456
    !open (unit = ff, file = 'part_not_in_ele.py')
    !write(ff,'(a)')'import matplotlib.pyplot as plt'
  END SUBROUTINE initGrid



  SUBROUTINE initHydro()
    !This Subroutine reads in the hydrodynamic information for the first 
    !  iteration
    USE PARAM_MOD, ONLY: numpar,ui,vi,uj,vj,us,ws,rho_nodes,u_nodes,v_nodes,   &
        filenum,filestep,tdim,numdigits,recordnum,days,dt, &
        readZeta,constZeta,readSalt,constSalt, &
        !readNetcdfSwdown,                &
        readTemp,constTemp,readDens,constDens,readU,constU,readV,constV,readW, &
        constW,readAks,constAks,WindIntensity,readIwind,constIwind,            &
        readUwind,constUwind,readVwind,constVwind,Zgrid,Wind,hydrobytes,       &  !--- CL-OGS:
!      *****   IMIOM      *****
          swan_prefix, swan_suffix,swan_filenum,WindWaveModel,SigWaveHeight,   &
          MeanWavePeriod,UWind_10,VWind_10,PeakDirection,PeakWaveLength,OilOn
!      ***** END IMIOM *****
        
    USE netcdf
    USE RANDOM_MOD, ONLY: genrand_real1 !--- CL-OGS
    USE CONVERT_MOD, ONLY: x2lon,y2lat  !--- CL-OGS
#include "VAR_IDs.h"
    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

    INTEGER :: STATUS,NCID,VID

    INTEGER :: i,j,k,t,count,counter,kmask,kbot
    INTEGER :: nfmax,nfn,nfnn,incrstepf,nf  !--- CL-OGS
    DOUBLE PRECISION :: fac                !--- CL-OGS
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION( :,:,: ) :: romZ !,romSwdown
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: romW,romKH,romS,romT, &
                                                romD,romU, romV
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION( :,:,: ) :: modelUwind,modelVwind,&
                                                         modelIwind  !--- CL-OGS
    !--- CL-OGS: following variables added to handle MITgcm-files
    !--- CL-OGS  (using a different file for every field variable )
    INTEGER :: ios,nvarf,nfilesin,ktlev,waiting,rand15
    INTEGER :: searchnode,nodestocopy,tcopy, k1, k2
    REAL, ALLOCATABLE, DIMENSION(:) :: tmpvec  
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: dbltmpvec  
    DOUBLE PRECISION, DIMENSION(3) :: salt_up,salt_around  
    DOUBLE PRECISION :: temp_up(3),temp_around(3)  
    DOUBLE PRECISION :: den_up(3),den_around(3)  
    DOUBLE PRECISION :: KH_up(3),KH_around(3)  

!      *****   IMIOM      *****
    INTEGER :: scounter                                   !swan file counter
    INTEGER :: startnum                                   !start record to read in file
    INTEGER :: ntloop                                     !numer of time to loop through opening new file each time
    INTEGER :: nloop                                      !loop counter
    INTEGER,PARAMETER :: interpol_from_cell_center_to_CArakawa=1,do_not_interpolate=0
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION( :,:,: ) :: swanHs,swantm01,      &
                                   swanpd,swanwl
!      ***** END IMIOM *****


    !ALLOCATE MODULE VARIABLES
    !ALLOCATE(t_Swdown(3,rho_nodes))
    ALLOCATE(t_zeta(3,rho_nodes))
    ALLOCATE(t_salt(3,rho_nodes,us))
    ALLOCATE(t_temp(3,rho_nodes,us))
    ALLOCATE(t_Wvel(3,rho_nodes,ws))
    ALLOCATE(t_den (3,rho_nodes,us))
    ALLOCATE(t_KH  (3,rho_nodes,ws))
    ALLOCATE(t_Uvel(3,  u_nodes,us))
    ALLOCATE(t_Vvel(3,  v_nodes,us))
    ALLOCATE(updatenodesbuffer(2,uj,3)) !--- CL-OGS

    !t_Swdown = 0.
    t_zeta = 0.
    t_salt = 0.
    t_temp = 0.
    t_den  = 0.
    t_KH   = 0.
    t_Uvel = 0.
    t_Vvel = 0.
    t_Wvel = 0.

    !ALLOCATE SUBROUTINE VARIABLES
    !ALLOCATE(romSwdown(vi,uj,3))
    ALLOCATE(romZ(vi,uj,3))
    ALLOCATE(romW(vi,uj,ws,3))
    ALLOCATE(romS(vi,uj,us,3))
    ALLOCATE(romT(vi,uj,us,3))
    ALLOCATE(romD(vi,uj,us,3))
    ALLOCATE(romU(ui,uj,us,3))
    ALLOCATE(romV(vi,vj,us,3))
    ALLOCATE(romKH(vi,uj,ws,3))
    ALLOCATE(tmpvec(vi)) !--- CL-OGS 
    ALLOCATE(dbltmpvec(vi)) !--- CL-OGS 

    !if(Wind .and. Zgrid)then !--- CL-OGS 
      ALLOCATE(t_uwind(3,  u_nodes))
      ALLOCATE(t_vwind(3,  v_nodes))
      t_uwind = 0
      t_vwind = 0 
      ALLOCATE(modelUwind(ui,uj,3))  
      ALLOCATE(modelVwind(vi,vj,3))  
    !endif
      ALLOCATE(modelIwind(vi,uj,3))
    if(WindIntensity .and. Zgrid)then
      ALLOCATE(t_iwind(3,  rho_nodes))
      t_iwind = 0
    endif
    !-----------------------------------------------------
    IF(OilOn)THEN
        !ALLOCATE IMIOM MODULE AND SUBROUTINE VARIABLES
        ALLOCATE(t_hsig (3,rho_nodes))
        ALLOCATE(t_tm01 (3,rho_nodes))
        ALLOCATE(t_pdir (3,rho_nodes)) 
        ALLOCATE(t_wlen (3,rho_nodes))

        t_hsig  = 0
        t_tm01  = 0
        t_pdir  = 0
        t_wlen  = 0

        ALLOCATE(swanHs(vi,uj,3))  !--- CL-OGS : changed 1 -> 3
        ALLOCATE(swantm01(vi,uj,3))!--- CL-OGS : changed 1 -> 3
        ALLOCATE(swanpd(vi,uj,3))  !--- CL-OGS : changed 1 -> 3
        ALLOCATE(swanwl(vi,uj,3))  !--- CL-OGS : changed 1 -> 3
        !ALLOCATE(modelUwind(ui,uj,1))  !--- CL-OGS : commented out 
        !ALLOCATE(modelVwind(vi,vj,1))  !--- CL-OGS : commented out
    END IF      !OilOn
    iint = 0

    
    !if(not Zgrid)then  !--- CL-OGS : restricting to ROMS hydro files !--- CL-OGS : commented out
    ! if(tdim .lt. 3)then                                                   !--- CL-OGS : commented out
    !      stepf = 1                                                       !--- CL-OGS : commented out
    ! else                                                                  !--- CL-OGS : commented out
    !      stepf = 3                                                       !--- CL-OGS : commented out
    ! end if                                                                !--- CL-OGS : commented out
    !endif                                                                  !--- CL-OGS : commented out
    !--- CL-OGS : commented out the IMIOM following code lines 
    !DO nloop = 1, 3

    ! !Open netCDF file
    ! if(nloop .eq. 1)then
    !       iint = 0
    !       startnum = 1
    ! else
    !       SELECT CASE(tdim)
    !             case(1)
    !                   iint = iint + 1
    !                   startnum = startnum
    !             case(2)
    !                   iint = iint + mod(nloop,tdim)
    !                   if(startnum .eq. 1)then
    !                         startnum = startnum + 1
    !                   else
    !                         startnum = startnum - 1
    !                   end if
    !             case default      !3 or greater
    !                   iint = iint
    !                   startnum = startnum + 1
    !       END SELECT
    ! end if
!    ***** END IMIOM *****

    ! Verfications of consistency of the input parameter for hydro files:
    if  (tdim == 0 .and. filestep.ne.0 ) then
        write(*,*) 'error inconsistency between tdim null=',tdim,              &
                   ' and filestep=',filestep
        stop
    endif

    if  (tdim.ne.0 .and. filestep == 0 )  then
        write(*,*) 'error inconsistency between tdim=',tdim,                   &
                   ' and filestep null=',filestep
        stop
    endif
    if  (tdim == 1 .and. recordnum.ne.1)  then
        write(*,*) 'error inconsistency between tdim unitary=',tdim,          &
                   ' and recordnum=',recordnum
        stop
    endif

    if (tdim==1)then ! each of the three timestep are in different files
      nfmax=3
      nfn=1
      nfnn=1
      incrstepf=1
    elseif(tdim==0 .or. (tdim>=(recordnum+3))) then ! all threetimesteps are in the same file
      nfmax=1
      nfn=1
      nfnn=3
      incrstepf=3
    else
      write(*,*)'case where the first 3 time steps are ',                      &
                 'in 2 different files not yet implemented.'
      write(*,*)'the program will now stop.'
      stop
    endif
    stepf=recordnum-1    !Forward step is (recordnum+3)rd time step of file

    if(Zgrid)then
       nfilesin=8 ! MITgcm binary outputs without wind files
       if(Wind) nfilesin=nfilesin+2 
       if(WindIntensity) nfilesin=nfilesin+1 
       !if(readNetcdfSwdown) nfilesin=nfilesin+1
    else 
       nfilesin=1 ! Roms NETcdf outputs
    endif

    !write(*,*)nfilesin,' variables will be initialised by init_hydro'
      
    !t_ijruv = (/175,195,155,175,175,195,155,175,175,195,155,175/)

    t_b = 1    !Back step is 1st time step in arrays
    t_c = 2    !Center step is 2nd time step in arrays
    t_f = 3    !Forward step is 3rd time step in arrays

    !Get i/j max/min for rho/u/v
    call setijruv()


    DO nf=1,nfmax
      !Open netCDF file
      if (nf>1) then
        iint=iint+filestep
        nfn=nfn+1
        nfnn=nfnn+1
      endif

      counter=iint+filenum  !176 + 1 = 177 --> June 26,1995
      countfilenum=counter
      stepf=stepf+incrstepf
      write(*,*)'reading record ',recordnum,':',recordnum+incrstepf-1

      ! Read in data for first three external time steps
      !------------------------------------
      if(readZeta)then  
        call read_data_from_file(VAR_ID_zeta,vi,uj,1,3,nf,nfn,nfnn,romZ,RNODE,recordnum,incrstepf,do_not_interpolate)
      else
        romZ = constZeta
      endif
      !------------------------------------
      if(readSalt)then
        call read_data_from_file(VAR_ID_salt,vi,uj,us,3,nf,nfn,nfnn,romS,RNODE,recordnum,incrstepf,do_not_interpolate)
      else
        romS = constSalt
      endif
      !------------------------------------
      if(readTemp)then  
        call read_data_from_file(VAR_ID_temp,vi,uj,us,3,nf,nfn,nfnn,romT,RNODE,recordnum,incrstepf,do_not_interpolate)
      else
        romT = constTemp
      endif
      !------------------------------------
      if(readDens)then  
        call read_data_from_file(VAR_ID_den,vi,uj,us,3,nf,nfn,nfnn,romD,RNODE,recordnum,incrstepf,do_not_interpolate)
      else
        romD = constDens
      endif
      !------------------------------------
      if(readU)then  
        call read_data_from_file(VAR_ID_uvel,ui,uj,us,3,nf,nfn,nfnn,romU,UNODE,recordnum,incrstepf,do_not_interpolate)
      else
        romU = constU
      endif
      !------------------------------------
      if(readV)then  
          call read_data_from_file(VAR_ID_vvel,vi,vj,us,3,nf,nfn,nfnn,romV,VNODE,recordnum,incrstepf,do_not_interpolate)
      else
        romV = constV
      endif
      !------------------------------------
      if(readW)then  
        call read_data_from_file(VAR_ID_wvel,vi,uj,ws,3,nf,nfn,nfnn,romW,RNODE,recordnum,incrstepf,do_not_interpolate)
      else
        romW = constW
      endif
      !------------------------------------
      if(readAks)then  
        call read_data_from_file(VAR_ID_kh,vi,uj,us,3,nf,nfn,nfnn,romKH,RNODE,recordnum,incrstepf,do_not_interpolate)
      else
        romKH = constAks
      endif
      !------------------------------------
      if(Wind .and.readUwind)then  
        call read_data_from_file(VAR_ID_uwind,ui,uj,1,3,nf,nfn,nfnn,modelUwind,UNODE,recordnum,incrstepf, &
                    interpol_from_cell_center_to_CArakawa)
      else
        modelUwind = constUwind
      endif       
      !------------------------------------
      if(Wind .and.readVwind)then  
        call read_data_from_file(VAR_ID_vwind,vi,vj,1,3,nf,nfn,nfnn,modelVwind,VNODE,recordnum,incrstepf, &
                    interpol_from_cell_center_to_CArakawa)
      else
        modelVwind = constVwind
      endif
      !------------------------------------
      if(readIwind)then  
        call read_data_from_file(VAR_ID_iwind,vi,uj,1,3,nf,nfn,nfnn,modelIwind,RNODE,recordnum,incrstepf,do_not_interpolate)
      else
        modelIwind = constIwind
      endif
      !------------------------------------


    ENDDO !nf=1,nfmax
      
      ! Store the ranges of nodes that were update
      updatenodesbuffer=0
      ! rho node range
      do j=t_ijruv(JMIN,RNODE),t_ijruv(JMAX,RNODE)
        updatenodesbuffer(1,j,1) = (j-1)*vi + t_ijruv(IMIN,RNODE) ! frst rnode at latitude j
        updatenodesbuffer(2,j,1) = (j-1)*vi + t_ijruv(IMAX,RNODE) ! last rnode at latitude j
      enddo
      write(*,'(2(a,2i5))')'updating rho nodes data in i=',t_ijruv(IMIN,RNODE),t_ijruv(IMAX,RNODE),         &
                ' j=',t_ijruv(JMIN,RNODE),t_ijruv(JMAX,RNODE)

     !write(*,'(a,F10.5,a,F10.5,a)')'lon=[ ',                                  &
     ! x2lon(rx(updatenodesbuffer(1,t_ijruv(JMIN,RNODE),1)),                            &
     !       ry(updatenodesbuffer(1,t_ijruv(JMIN,RNODE),1))),' : ',                     &
     ! x2lon(rx(updatenodesbuffer(2,t_ijruv(JMAX,RNODE),1)),                            &
     !       ry(updatenodesbuffer(2,t_ijruv(JMAX,RNODE),1))),' ]'
     !write(*,'(a,F10.5,a,F10.5,a)')'lat=[ ',                                  &
     ! y2lat(ry(updatenodesbuffer(1,t_ijruv(JMIN,RNODE),1))),' : ' ,                    &
     ! y2lat(ry(updatenodesbuffer(2,t_ijruv(JMAX,RNODE),1))),' ]'

      ! u node range
       
      do j=t_ijruv(JMIN,UNODE),t_ijruv(JMAX,UNODE)
        updatenodesbuffer(1,j,2) = (j-1)*ui + t_ijruv(IMIN,UNODE) ! frst unode at latitude j
        updatenodesbuffer(2,j,2) = (j-1)*ui + t_ijruv(IMAX,UNODE) ! last unode at latitude j
      enddo
      write(*,*)'updating u nodes data in i=',t_ijruv(IMIN,UNODE),t_ijruv(IMAX,UNODE),' j=',     &
                   t_ijruv(JMIN,UNODE),t_ijruv(JMAX,UNODE)
     !write(*,'(a,F10.5,a,F10.5,a)')'lon=[ ',                                  &
     ! x2lon(ux(updatenodesbuffer(1,t_ijruv(JMIN,UNODE),2)),                            &
     !       uy(updatenodesbuffer(1,t_ijruv(JMIN,UNODE),2))),' : ',                     &
     ! x2lon(ux(updatenodesbuffer(2,t_ijruv(JMAX,UNODE),2)),                            &
     !       uy(updatenodesbuffer(2,t_ijruv(JMAX,UNODE),2))),' ]'
     !write(*,'(a,F10.5,a,F10.5,a)')'lat=[ ',                                  &
     ! y2lat(uy(updatenodesbuffer(1,t_ijruv(JMIN,UNODE),2))),' : ',                     &
     ! y2lat(uy(updatenodesbuffer(2,t_ijruv(JMAX,UNODE),2))),' ]'

      ! v node range
      do j=t_ijruv(JMIN,VNODE),t_ijruv(JMAX,VNODE)
        updatenodesbuffer(1,j,3) = (j-1)*vi + t_ijruv(IMIN,VNODE)  ! frst vnode at latitude j
        updatenodesbuffer(2,j,3) = (j-1)*vi + t_ijruv(IMAX,VNODE) ! last vnode at latitude j
      enddo
      write(*,*)'updating v nodes data in i=',t_ijruv(IMIN,VNODE),t_ijruv(IMAX,VNODE),' j=',    &
                t_ijruv(JMIN,VNODE),t_ijruv(JMAX,VNODE)
     !write(*,'(a,F10.5,a,F10.5,a)')'lon=[ ',                                  &
     ! x2lon(vx(updatenodesbuffer(1,t_ijruv(JMIN,VNODE),3)),                           &
     !       vy(updatenodesbuffer(1,t_ijruv(JMIN,VNODE),3))),' : ',                    &
     ! x2lon(vx(updatenodesbuffer(2,t_ijruv(JMAX,VNODE),3)),                           &
     !       vy(updatenodesbuffer(2,t_ijruv(JMAX,VNODE),3))),' ]'                        
     !write(*,'(a,F10.5,a,F10.5,a)')'lat=[ ',                                  &
     ! y2lat(vy(updatenodesbuffer(1,t_ijruv(JMIN,VNODE),3))),' : ',                    &
     ! y2lat(vy(updatenodesbuffer(2,t_ijruv(JMAX,VNODE),3))),' ]'


      !Reshape input to fit node numbers assigned to elements
      do j=t_ijruv(JMIN,RNODE),t_ijruv(JMAX,RNODE)
        !write(*,*)'initHydro nodes',(j-1)*vi + t_ijruv(IMIN,RNODE),':',& 
        !   (j-1)*vi + t_ijruv(IMAX,RNODE)
        do i=t_ijruv(IMIN,RNODE),t_ijruv(IMAX,RNODE)
          count = (j-1)*vi + i
          do k=1,us
            kmask=min(k,us_tridim) 
            t_salt(1,count,k) = romS(i,j,k,1) * m_r(i,j,kmask)
            t_salt(2,count,k) = romS(i,j,k,2) * m_r(i,j,kmask)
            t_salt(3,count,k) = romS(i,j,k,3) * m_r(i,j,kmask)
            t_temp(1,count,k) = romT(i,j,k,1) * m_r(i,j,kmask)
            t_temp(2,count,k) = romT(i,j,k,2) * m_r(i,j,kmask)
            t_temp(3,count,k) = romT(i,j,k,3) * m_r(i,j,kmask)
            t_Wvel(1,count,k+1) = romW(i,j,k+1,1) * m_r(i,j,kmask)
            t_Wvel(2,count,k+1) = romW(i,j,k+1,2) * m_r(i,j,kmask)
            t_Wvel(3,count,k+1) = romW(i,j,k+1,3) * m_r(i,j,kmask)
            t_den(1,count,k) = (romD(i,j,k,1) + DBLE(1000.0)) * m_r(i,j,kmask)
            t_den(2,count,k) = (romD(i,j,k,2) + DBLE(1000.0)) * m_r(i,j,kmask)
            t_den(3,count,k) = (romD(i,j,k,3) + DBLE(1000.0)) * m_r(i,j,kmask)
            t_KH(1,count,k+1) = romKH(i,j,k+1,1) * m_r(i,j,kmask)                
            t_KH(2,count,k+1) = romKH(i,j,k+1,2) * m_r(i,j,kmask)                
            t_KH(3,count,k+1) = romKH(i,j,k+1,3) * m_r(i,j,kmask)                
          enddo                                         
          t_Wvel(1,count,1) = romW(i,j,1,1) * m_r(i,j,1)  ! BEUG? SHOULD BE 0?               
          t_Wvel(2,count,1) = romW(i,j,1,2) * m_r(i,j,1)  ! BEUG? SHOULD BE 0?               
          t_Wvel(3,count,1) = romW(i,j,1,3) * m_r(i,j,1)  ! BEUG? SHOULD BE 0?                
          t_KH(1,count,1) =  romKH(i,j,1,1) * m_r(i,j,1)
          t_KH(2,count,1) =  romKH(i,j,1,2) * m_r(i,j,1)
          t_KH(3,count,1) =  romKH(i,j,1,3) * m_r(i,j,1)
          t_zeta(1,count) =    romZ(i,j,1) *    m_r(i,j,us_tridim)
          t_zeta(2,count) =    romZ(i,j,2) *    m_r(i,j,us_tridim)
          t_zeta(3,count) =    romZ(i,j,3) *    m_r(i,j,us_tridim)
          !t_Swdown(1,count) =    romSwdown(i,j,1) 
          !t_Swdown(2,count) =    romSwdown(i,j,2) 
          !t_Swdown(3,count) =    romSwdown(i,j,3) 
        enddo
      enddo

      do j=t_ijruv(JMIN,UNODE),t_ijruv(JMAX,UNODE)
        do i=t_ijruv(IMIN,UNODE),t_ijruv(IMAX,UNODE)
          count = (j-1)*ui + i
          do k=1,us
            kmask=min(k,us_tridim) 
            t_Uvel(1,count,k) = romU(i,j,k,1) * m_u(i,j,kmask)
            t_Uvel(2,count,k) = romU(i,j,k,2) * m_u(i,j,kmask)
            t_Uvel(3,count,k) = romU(i,j,k,3) * m_u(i,j,kmask)
          enddo
        enddo
      enddo

      do j=t_ijruv(JMIN,VNODE),t_ijruv(JMAX,VNODE)
        do i=t_ijruv(IMIN,VNODE),t_ijruv(IMAX,VNODE)
          count = (j-1)*vi + i
          do k=1,us
            kmask=min(k,us_tridim)
            t_Vvel(1,count,k) = romV(i,j,k,1) * m_v(i,j,kmask)
            t_Vvel(2,count,k) = romV(i,j,k,2) * m_v(i,j,kmask)
            t_Vvel(3,count,k) = romV(i,j,k,3) * m_v(i,j,kmask)
          enddo
        enddo    
      enddo


      if(Zgrid)then
        ! COPY DATA OF FIRST CELL ABOVE THE BOTTOM TO ALL CELLS BELOW 
        nodestocopy=1
        do j=t_ijruv(JMIN,RNODE),t_ijruv(JMAX,RNODE)
         !write(*,'(8(a,i4),a)')&
         !  '---------------------------------------------- j=',j,&
         !' (in [',t_ijruv(JMIN,RNODE),',',t_ijruv(JMAX,RNODE), &
         !']) i=[',t_ijruv(IMIN,RNODE),',',t_ijruv(IMAX,RNODE),'] --- ui=',ui, &
         !' uj=',uj,' us=',us,'----------------------------------------------'
          do i=t_ijruv(IMIN,RNODE),t_ijruv(IMAX,RNODE)
            count = (j-1)*vi + i
            kbot=BottomK(i,j,1)
            if(kbot>1)then
             if(kbot<ws)then
              t_salt(:,count,kbot-1) = t_salt(:,count,kbot)
              t_temp(:,count,kbot-1) = t_temp(:,count,kbot)
              t_den(:,count,kbot-1)  = t_den(:,count,kbot) 
              t_KH(:,count,kbot-1)   = t_KH(:,count,kbot)  
             endif
             do searchnode=nodestocopy,NUM_COPNOD(RNODE)
               if(Node_COPNOD(RNODE,searchnode).le.count)nodestocopy=searchnode
               if(Node_COPNOD(RNODE,searchnode).ge.count)exit
             enddo
             do while(Node_COPNOD(RNODE,nodestocopy).eq.count)
               k1=klev_COPNOD(RNODE,nodestocopy,2)
               k2=klev_COPNOD(RNODE,nodestocopy,1)
               do tcopy=1,3
                 t_salt(tcopy,count,k1:k2)=(                        &
                 t_salt(tcopy,Nghb_COPNOD(RNODE,nodestocopy,1),k1:k2) &
                             *Coef_COPNOD(RNODE,nodestocopy,1) +  &
                 t_salt(tcopy,Nghb_COPNOD(RNODE,nodestocopy,2),k1:k2) &
                             *Coef_COPNOD(RNODE,nodestocopy,2) +  &
                 t_salt(tcopy,Nghb_COPNOD(RNODE,nodestocopy,3),k1:k2) &
                             *Coef_COPNOD(RNODE,nodestocopy,3) +  &
                 t_salt(tcopy,Nghb_COPNOD(RNODE,nodestocopy,4),k1:k2) &
                             *Coef_COPNOD(RNODE,nodestocopy,4) ) 
                 t_temp(tcopy,count,k1:k2)=(                        &
                 t_temp(tcopy,Nghb_COPNOD(RNODE,nodestocopy,1),k1:k2) &
                             *Coef_COPNOD(RNODE,nodestocopy,1) +  &
                 t_temp(tcopy,Nghb_COPNOD(RNODE,nodestocopy,2),k1:k2) &
                             *Coef_COPNOD(RNODE,nodestocopy,2) +  &
                 t_temp(tcopy,Nghb_COPNOD(RNODE,nodestocopy,3),k1:k2) &
                             *Coef_COPNOD(RNODE,nodestocopy,3) +  &
                 t_temp(tcopy,Nghb_COPNOD(RNODE,nodestocopy,4),k1:k2) &
                             *Coef_COPNOD(RNODE,nodestocopy,4) ) 
                 t_den(tcopy,count,k1:k2)=(                        &
                 t_den(tcopy,Nghb_COPNOD(RNODE,nodestocopy,1),k1:k2) &
                            *Coef_COPNOD(RNODE,nodestocopy,1) +  &
                 t_den(tcopy,Nghb_COPNOD(RNODE,nodestocopy,2),k1:k2) &
                            *Coef_COPNOD(RNODE,nodestocopy,2) +  &
                 t_den(tcopy,Nghb_COPNOD(RNODE,nodestocopy,3),k1:k2) &
                            *Coef_COPNOD(RNODE,nodestocopy,3) +  &
                 t_den(tcopy,Nghb_COPNOD(RNODE,nodestocopy,4),k1:k2) &
                            *Coef_COPNOD(RNODE,nodestocopy,4) ) 
                 t_KH(tcopy,count,k1:k2)=(                        &
                 t_KH(tcopy,Nghb_COPNOD(RNODE,nodestocopy,1),k1:k2) &
                           *Coef_COPNOD(RNODE,nodestocopy,1) +  &
                 t_KH(tcopy,Nghb_COPNOD(RNODE,nodestocopy,2),k1:k2) &
                           *Coef_COPNOD(RNODE,nodestocopy,2) +  &
                 t_KH(tcopy,Nghb_COPNOD(RNODE,nodestocopy,3),k1:k2) &
                           *Coef_COPNOD(RNODE,nodestocopy,3) +  &
                 t_KH(tcopy,Nghb_COPNOD(RNODE,nodestocopy,4),k1:k2) &
                           *Coef_COPNOD(RNODE,nodestocopy,4) ) 
               enddo
               nodestocopy = nodestocopy +1
               if(nodestocopy.gt.NUM_COPNOD(RNODE))exit
             enddo !while
            endif
            if(nodestocopy.gt.NUM_COPNOD(RNODE))exit
          enddo    
          if(nodestocopy.gt.NUM_COPNOD(RNODE))exit
        enddo

        nodestocopy=1
        do j=t_ijruv(JMIN,UNODE),t_ijruv(JMAX,UNODE)
          do i=t_ijruv(IMIN,UNODE),t_ijruv(IMAX,UNODE)
           count = (j-1)*ui + i
           kbot=BottomK(i,j,2)
           if( kbot>1)then
             if(kbot<ws)then
               t_Uvel(:,count,kbot-1) =  t_Uvel(:,count,kbot)
             endif
             do searchnode=nodestocopy,NUM_COPNOD(UNODE)
               if(Node_COPNOD(UNODE,searchnode).le.count)nodestocopy=searchnode
               if(Node_COPNOD(UNODE,searchnode).ge.count)exit
             enddo
             do while(Node_COPNOD(UNODE,nodestocopy).eq.count)
               k1=klev_COPNOD(UNODE,nodestocopy,2)
               k2=klev_COPNOD(UNODE,nodestocopy,1)
               do tcopy=1,3
                 t_Uvel(tcopy,count,k1:k2)=(                        &
                 t_Uvel(tcopy,Nghb_COPNOD(UNODE,nodestocopy,1),k1:k2) &
                             *Coef_COPNOD(UNODE,nodestocopy,1) +  &
                 t_Uvel(tcopy,Nghb_COPNOD(UNODE,nodestocopy,2),k1:k2) &
                             *Coef_COPNOD(UNODE,nodestocopy,2) +  &
                 t_Uvel(tcopy,Nghb_COPNOD(UNODE,nodestocopy,3),k1:k2) &
                             *Coef_COPNOD(UNODE,nodestocopy,3) +  &
                 t_Uvel(tcopy,Nghb_COPNOD(UNODE,nodestocopy,4),k1:k2) &
                             *Coef_COPNOD(UNODE,nodestocopy,4) ) 
               enddo  !tcopy
               nodestocopy = nodestocopy +1
               if(nodestocopy.gt.NUM_COPNOD(UNODE))exit
             enddo !While
           endif
           if(nodestocopy.gt.NUM_COPNOD(UNODE))exit
          enddo    
          if(nodestocopy.gt.NUM_COPNOD(UNODE))exit
        enddo


        nodestocopy=1
        do j=t_ijruv(JMIN,VNODE),t_ijruv(JMAX,VNODE)
          do i=t_ijruv(IMIN,VNODE),t_ijruv(IMAX,VNODE)
           count = (j-1)*vi + i
           kbot=BottomK(i,j,3)
           if( kbot>1)then
             if(kbot<ws)then
               t_Vvel(:,count,kbot-1) =  t_Vvel(:,count,kbot)
             endif
             do searchnode=nodestocopy,NUM_COPNOD(VNODE)
               if(Node_COPNOD(VNODE,searchnode).le.count)nodestocopy=searchnode
               if(Node_COPNOD(VNODE,searchnode).ge.count)exit
             enddo
             do while(Node_COPNOD(VNODE,nodestocopy).eq.count)
               k1=klev_COPNOD(VNODE,nodestocopy,2)
               k2=klev_COPNOD(VNODE,nodestocopy,1)
               do tcopy=1,3
                 t_Vvel(tcopy,count,k1:k2)=(                        &
                 t_Vvel(tcopy,Nghb_COPNOD(VNODE,nodestocopy,1),k1:k2) &
                             *Coef_COPNOD(VNODE,nodestocopy,1) +  &
                 t_Vvel(tcopy,Nghb_COPNOD(VNODE,nodestocopy,2),k1:k2) &
                             *Coef_COPNOD(VNODE,nodestocopy,2) +  &
                 t_Vvel(tcopy,Nghb_COPNOD(VNODE,nodestocopy,3),k1:k2) &
                             *Coef_COPNOD(VNODE,nodestocopy,3) +  &
                 t_Vvel(tcopy,Nghb_COPNOD(VNODE,nodestocopy,4),k1:k2) &
                             *Coef_COPNOD(VNODE,nodestocopy,4) ) 
               enddo  !tcopy
               nodestocopy = nodestocopy +1
               if(nodestocopy.gt.NUM_COPNOD(VNODE))exit
             enddo !While
           endif
           if(nodestocopy.gt.NUM_COPNOD(VNODE))exit
          enddo    
          if(nodestocopy.gt.NUM_COPNOD(VNODE))exit
        enddo      
  
        nodestocopy=1
        do j=t_ijruv(JMIN,RNODE),t_ijruv(JMAX,RNODE)
          do i=t_ijruv(IMIN,RNODE),t_ijruv(IMAX,RNODE)
            count = (j-1)*vi + i
            kbot=BottomK(i,j,1)
            if(kbot>1)then
             if(kbot<ws)then
              t_Wvel(:,count,kbot-1) = t_Wvel(:,count,kbot)
             endif
             do searchnode=nodestocopy,NUM_COPNOD(RNODE)
               if(Node_COPNOD(RNODE,searchnode).le.count)nodestocopy=searchnode
               if(Node_COPNOD(RNODE,searchnode).ge.count)exit
             enddo
             do while(Node_COPNOD(RNODE,nodestocopy).eq.count)
               k1=klev_COPNOD(RNODE,nodestocopy,2)
               k2=klev_COPNOD(RNODE,nodestocopy,1)
               do tcopy=1,3
                 t_Wvel(tcopy,count,k1:k2)=(                        &
                 t_Wvel(tcopy,Nghb_COPNOD(RNODE,nodestocopy,1),k1:k2) &
                             *Coef_COPNOD(RNODE,nodestocopy,1) +  &
                 t_Wvel(tcopy,Nghb_COPNOD(RNODE,nodestocopy,2),k1:k2) &
                             *Coef_COPNOD(RNODE,nodestocopy,2) +  &
                 t_Wvel(tcopy,Nghb_COPNOD(RNODE,nodestocopy,3),k1:k2) &
                             *Coef_COPNOD(RNODE,nodestocopy,3) +  &
                 t_Wvel(tcopy,Nghb_COPNOD(RNODE,nodestocopy,4),k1:k2) &
                             *Coef_COPNOD(RNODE,nodestocopy,4) ) 
               enddo
               nodestocopy = nodestocopy +1
               if(nodestocopy.gt.NUM_COPNOD(RNODE))exit
             enddo !while
            endif
            if(nodestocopy.gt.NUM_COPNOD(RNODE))exit
          enddo    
          if(nodestocopy.gt.NUM_COPNOD(RNODE))exit
        enddo

 
     endif ! (Zgrid)

    ! --- CL-OGS : 
    write(*,*)'counter at inithydro=',counter-2*432,counter-432,counter
      ! --- CL-OGS : 
       do j=t_ijruv(JMIN,UNODE),t_ijruv(JMAX,UNODE)
        do i=t_ijruv(IMIN,UNODE),t_ijruv(IMAX,UNODE)
          count = (j-1)*ui + i
          t_uwind(1,count) =    modelUwind(i,j,1) *    m_u(i,j,us_tridim)
          t_uwind(2,count) =    modelUwind(i,j,2) *    m_u(i,j,us_tridim)
          t_uwind(3,count) =    modelUwind(i,j,3) *    m_u(i,j,us_tridim)
        enddo
       enddo
       do j=t_ijruv(JMIN,VNODE),t_ijruv(JMAX,VNODE)
        do i=t_ijruv(IMIN,VNODE),t_ijruv(IMAX,VNODE)
          count = (j-1)*vi + i
          t_vwind(1,count) =    modelVwind(i,j,1) *    m_v(i,j,us_tridim)
          t_vwind(2,count) =    modelVwind(i,j,2) *    m_v(i,j,us_tridim)
          t_vwind(3,count) =    modelVwind(i,j,3) *    m_v(i,j,us_tridim)
        enddo
       enddo
      if(WindIntensity .and. Zgrid)then
       do j=t_ijruv(JMIN,RNODE),t_ijruv(JMAX,RNODE)
        do i=t_ijruv(IMIN,RNODE),t_ijruv(IMAX,RNODE)
          count = (j-1)*vi + i
          t_iwind(1,count) =    modelIwind(i,j,1) *    m_r(i,j,us_tridim)
          t_iwind(2,count) =    modelIwind(i,j,2) *    m_r(i,j,us_tridim)
          t_iwind(3,count) =    modelIwind(i,j,3) *    m_r(i,j,us_tridim)
        enddo
       enddo
      endif
      !  ***    IMIOM *****
      ! WIND WAVE MODEL DATA  ------------------------------------
      IF(OilOn)then
       if(WindWaveModel)THEN
         if (tdim==1)then ! each of the three timestep are in different files
           nfmax=3
           nfn=1
           nfnn=1
           incrstepf=1
         elseif(tdim==0 .or. (tdim>=(recordnum+3))) then ! all threetimesteps are in the same file
           nfmax=1
           nfn=1
           nfnn=3
           incrstepf=3
         else
           write(*,*)'case where the first 3 time steps are in 2 different files not yet implemented.'
           write(*,*)'the program will now stop.'
           stop
         endif
         stepf=recordnum-1    !Forward step is (recordnum+3)rd time step of file
         nfilesin=1 ! Roms NETcdf outputs
         DO nf=1,nfmax
            if (nf>1) then
              iint=iint+filestep
              nfn=nfn+1
              nfnn=nfnn+1
            endif
            counter=iint+filenum  !176 + 1 = 177 --> June 26,1995
            countfilenum=counter

            stepf=stepf+incrstepf
            nvarf=1
            scounter = iint + swan_filenum

            call set_filename(VAR_ID_swan,scounter,swannm)

            ! Read in data for first three external time steps
            STATUS = NF90_OPEN(TRIM(swannm), NF90_NOWRITE, NCID)
            if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem NF90_OPEN'
            if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

                ! **** Hsig ****
                startz(1)=t_ijruv(IMIN,RNODE)
                startz(2)=t_ijruv(JMIN,RNODE)
                startz(3)=recordnum

                countz(1)=t_ijruv(IMAX,RNODE)-t_ijruv(IMIN,RNODE)+1
                countz(2)=t_ijruv(JMAX,RNODE)-t_ijruv(JMIN,RNODE)+1
                countr(3)=incrstepf

                STATUS = NF90_INQ_VARID(NCID,'Hs',VID)
                if (STATUS .NE. NF90_NOERR) then
                  write(*,*) 'Problem find Hs'
                  write(*,*) NF90_STRERROR(STATUS)
                  stop
                endif

                STATUS = NF90_GET_VAR(NCID,VID,swanHs(t_ijruv(IMIN,RNODE):t_ijruv(IMAX,RNODE),   &
                              t_ijruv(JMIN,RNODE):t_ijruv(JMAX,RNODE),nfn:nfnn),STARTz,COUNTz)

                if (STATUS .NE. NF90_NOERR) then
                  write(*,*) 'Problem read SwanHs array'
                  write(*,*) NF90_STRERROR(STATUS)
                  stop
                endif

                ! **** tm01 ****
                startz(1)=t_ijruv(IMIN,RNODE)
                startz(2)=t_ijruv(JMIN,RNODE)
                startz(3)=recordnum

                countz(1)=t_ijruv(IMAX,RNODE)-t_ijruv(IMIN,RNODE)+1
                countz(2)=t_ijruv(JMAX,RNODE)-t_ijruv(JMIN,RNODE)+1
                countr(3)=incrstepf

                STATUS = NF90_INQ_VARID(NCID,'tm01',VID)
                if (STATUS .NE. NF90_NOERR) then
                  write(*,*) 'Problem find tm01'
                  write(*,*) NF90_STRERROR(STATUS)
                  stop
                endif

                STATUS = NF90_GET_VAR(NCID,VID,swantm01(t_ijruv(IMIN,RNODE):t_ijruv(IMAX,RNODE), &
                           t_ijruv(JMIN,RNODE):t_ijruv(JMAX,RNODE),nfn:nfnn),STARTz,COUNTz)

                if (STATUS .NE. NF90_NOERR) then
                  write(*,*) 'Problem read swantm01 array'
                  write(*,*) NF90_STRERROR(STATUS)
                  stop
                endif

                ! **** u10 ****
                startr(1)=t_ijruv(IMIN,UNODE)
                startr(2)=t_ijruv(JMIN,UNODE)
                startz(3)=recordnum

                countr(1)=t_ijruv(IMAX,UNODE)-t_ijruv(IMIN,UNODE)+1
                countr(2)=t_ijruv(JMAX,UNODE)-t_ijruv(JMIN,UNODE)+1
                countr(3)=incrstepf

                STATUS = NF90_INQ_VARID(NCID,'u10',VID)
                if (STATUS .NE. NF90_NOERR) then
                  write(*,*) 'Problem find u10'
                  write(*,*) NF90_STRERROR(STATUS)
                  stop
                endif

                STATUS=NF90_GET_VAR(NCID,VID,modelUwind(t_ijruv(IMIN,UNODE):t_ijruv(IMAX,UNODE), &
                               t_ijruv(JMIN,UNODE):t_ijruv(JMAX,UNODE),nfn:nfnn),STARTz,COUNTz)

                if (STATUS .NE. NF90_NOERR) then
                  write(*,*) 'Problem read swanU array'
                  write(*,*) NF90_STRERROR(STATUS)
                  stop
                endif

                ! **** V10 ****
                startr(1)=t_ijruv(IMIN,VNODE)
                startr(2)=t_ijruv(JMIN,VNODE)
                startz(3)=recordnum

                countr(1)=t_ijruv(IMAX,VNODE)-t_ijruv(IMIN,VNODE)+1
                countr(2)=t_ijruv(JMAX,VNODE)-t_ijruv(JMIN,VNODE)+1
                countr(3)=incrstepf

                STATUS = NF90_INQ_VARID(NCID,'v10',VID)
                if (STATUS .NE. NF90_NOERR) then
                  write(*,*) 'Problem find v10'
                  write(*,*) NF90_STRERROR(STATUS)
                  stop
                endif

                STATUS=NF90_GET_VAR(NCID,VID,modelVwind(t_ijruv(IMIN,VNODE):t_ijruv(IMAX,VNODE),&
                       t_ijruv(JMIN,VNODE):t_ijruv(JMAX,VNODE),nfn:nfnn),STARTz,COUNTz)

                if (STATUS .NE. NF90_NOERR) then
                  write(*,*) 'Problem read swanV array'
                  write(*,*) NF90_STRERROR(STATUS)
                  stop
                endif

                ! **** PeakDir ****
                startz(1)=t_ijruv(IMIN,RNODE)
                startz(2)=t_ijruv(JMIN,RNODE)
                startz(3)=recordnum

                countz(1)=t_ijruv(IMAX,RNODE)-t_ijruv(IMIN,RNODE)+1
                countz(2)=t_ijruv(JMAX,RNODE)-t_ijruv(JMIN,RNODE)+1
                countr(3)=incrstepf

                STATUS = NF90_INQ_VARID(NCID,'Pd',VID)
                if (STATUS .NE. NF90_NOERR) then
                  write(*,*) 'Problem find Pd'
                  write(*,*) NF90_STRERROR(STATUS)
                  stop
                endif

                STATUS=NF90_GET_VAR(NCID,VID,swanpd(t_ijruv(IMIN,RNODE):t_ijruv(IMAX,RNODE),     &
                             t_ijruv(JMIN,RNODE):t_ijruv(JMAX,RNODE),nfn:nfnn),STARTz,COUNTz)

                if (STATUS .NE. NF90_NOERR) then
                  write(*,*) 'Problem read swanpd array'
                  write(*,*) NF90_STRERROR(STATUS)
                  stop
                endif

                ! **** PeakWaveLength ****
                startz(1)=t_ijruv(IMIN,RNODE)
                startz(2)=t_ijruv(JMIN,RNODE)
                startz(3)=recordnum

                countz(1)=t_ijruv(IMAX,RNODE)-t_ijruv(IMIN,RNODE)+1
                countz(2)=t_ijruv(JMAX,RNODE)-t_ijruv(JMIN,RNODE)+1
                countr(3)=incrstepf

                STATUS = NF90_INQ_VARID(NCID,'Pwl',VID)
                if (STATUS .NE. NF90_NOERR) then
                  write(*,*) 'Problem find Pwl'
                  write(*,*) NF90_STRERROR(STATUS)
                  stop
                endif

                STATUS = NF90_GET_VAR(NCID,VID,swanwl(t_ijruv(IMIN,RNODE):t_ijruv(IMAX,RNODE),   &
                                 t_ijruv(JMIN,RNODE):t_ijruv(JMAX,RNODE),nfn:nfnn),STARTz,COUNTz)

                if (STATUS .NE. NF90_NOERR) then
                  write(*,*) 'Problem read swanwl array'
                  write(*,*) NF90_STRERROR(STATUS)
                  stop
                endif

            !close the dataset and reassign the NCID
            STATUS = NF90_CLOSE(NCID)
         ENDDO
        else   ! if WindWaveModel
              write(*,*)'SigWaveHeight=',SigWaveHeight
              swanHs   = SigWaveHeight
              swantm01 = MeanWavePeriod
              if(.not.readUwind)  modelUwind = UWind_10
              if(.not.readVwind)  modelVwind = VWind_10
              swanpd   = PeakDirection
              swanwl   = PeakWaveLength
        endif                ! if WindWaveModel

        !Reshape input to fit node numbers assigned to elements
        DO nloop=1,3
         do j=t_ijruv(JMIN,RNODE),t_ijruv(JMAX,RNODE)
            do i=t_ijruv(IMIN,RNODE),t_ijruv(IMAX,RNODE)
              count = (j-1)*vi + i
              t_hsig(nloop,count) = swanHs(i,j,nloop) * m_r(i,j,us)
              t_tm01(nloop,count) = swantm01(i,j,nloop) * m_r(i,j,us)
              t_pdir(nloop,count) = swanpd(i,j,nloop) * m_r(i,j,us)
              t_wlen(nloop,count) = swanwl(i,j,nloop) * m_r(i,j,us)
            enddo
         enddo
        ENDDO
        if(windwavemodel.or.(UWind_10.ne.0 .and. (.not. Zgrid)))then     !only overwrite ROMS t_uwind (from sustr) if using windwavesmodel OR overwriting by constant U wind
              write(*,*)'overwrite ROMS U Wind'
              DO nloop=1,3
                do j=t_ijruv(JMIN,UNODE),t_ijruv(JMAX,UNODE)
                  do i=t_ijruv(IMIN,UNODE),t_ijruv(IMAX,UNODE)
                    count = (j-1)*ui + i
                    t_uwind(nloop,count) = modelUwind(i,j,nloop) * m_u(i,j,us)
                  enddo
                enddo
              ENDDO
        end if
        if(windwavemodel.or.(VWind_10.ne.0 .and. (.not. Zgrid)))then     !only overwrite ROMS t_uwind (from sustr) if using windwavesmodel OR overwriting by constant U wind
             write(*,*)'overwrite ROMS V Wind'
              DO nloop=1,3
                do j=t_ijruv(JMIN,VNODE),t_ijruv(JMAX,VNODE)
                  do i=t_ijruv(IMIN,VNODE),t_ijruv(IMAX,VNODE)
                    count = (j-1)*vi + i
                t_vwind(nloop,count) = modelVwind(i,j,nloop) * m_v(i,j,us)
                  enddo
                enddo
              ENDDO
        end if
      END IF        ! if OilOn


     !DEALLOCATE SUBROUTINE VARIABLES
     if(OilOn)then
             DEALLOCATE(swanHs,swantm01,swanpd,swanwl)
     END IF        !if OilOn


    !DEALLOCATE SUBROUTINE VARIABLES
    DEALLOCATE(romZ,romW,romD,romKH,romS,romT,romU,romV,&
                 modelUwind,modelVwind,modelIwind)
    !DEALLOCATE(romSwdown)

    write(*,*)'counter at end of inithydro=',counter
  END SUBROUTINE initHydro


  SUBROUTINE updateHydro()
    USE PARAM_MOD, ONLY: ui,vi,uj,vj,us,ws,tdim,rho_nodes,u_nodes,v_nodes,     &
        filenum,numdigits,readZeta,constZeta,readSalt,constSalt, &
        readTemp,constTemp,readDens,constDens,readU,constU,readV,constV,readW, &
        constW,readAks,constAks,&
        !readNetcdfSwdown,                                    &
        startfile,filestep,                                              &
        readUwind,constUwind,readVwind,constVwind,Zgrid,Wind,hydrobytes,       &  !--- CL-OGS:
        WindIntensity,readIwind,constIwind,                                    &
!      *****   IMIOM      *****
          swan_prefix, swan_suffix,swan_filenum,WindWaveModel,SigWaveHeight,   &
          MeanWavePeriod,UWind_10,VWind_10,PeakDirection,PeakWaveLength,OilOn
!      ***** END IMIOM *****
    USE netcdf
    USE RANDOM_MOD, ONLY: genrand_real1
#include "VAR_IDs.h"
    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

    INTEGER :: STATUS,NCID,VID,FID

    INTEGER :: i,j,k,t,count,counter,kmask,kbot

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION( :,:,: ) :: romZf!,romSwdownf
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: romSf,romTf,romDf,    &
                                romUf,romVf,romWf,romKHf
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION( :,:,: ) ::modelUwindf,modelVwindf    !--- CL-OGS
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION( :,:,: ) ::modelIwindf    !--- CL-OGS
    INTEGER :: ios,nvarf,nfilesin,ktlev,waiting,rand15
    INTEGER :: searchnode,nodestocopy, k1, k2
    REAL, ALLOCATABLE, DIMENSION(:) :: tmpvec
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: dbltmpvec  
    INTEGER,PARAMETER :: interpol_from_cell_center_to_CArakawa=1,do_not_interpolate=0
 
    !IMIOM
    INTEGER :: scounter
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: swanHsf,swantm01f,      &
                                        swanpdf,swanwlf

    !ALLOCATE SUBROUTINE VARIABLES
    !ALLOCATE(romSwdownf(vi,uj,1))
    ALLOCATE(romZf(vi,uj,1))
    ALLOCATE(romSf(vi,uj,us,1))
    ALLOCATE(romTf(vi,uj,us,1))
    ALLOCATE(romDf(vi,uj,us,1))
    ALLOCATE(romUf(ui,uj,us,1))
    ALLOCATE(romVf(vi,vj,us,1))
    ALLOCATE(romWf(vi,uj,ws,1))
    ALLOCATE(romKHf(vi,uj,ws,1))
    ALLOCATE(modelUwindf(ui,uj,1))
    ALLOCATE(modelVwindf(vi,vj,1))
    ALLOCATE(modelIwindf(vi,uj,1)) 
    ALLOCATE(tmpvec(vi))
    ALLOCATE(dbltmpvec(vi))
    if(OilOn)then! .and. WindWaveModel)then
        ALLOCATE(swanHsf(vi,uj,1))
        ALLOCATE(swantm01f(vi,uj,1))
        ALLOCATE(swanpdf(vi,uj,1))
        ALLOCATE(swanwlf(vi,uj,1))
    end if
     
    !Rotate Indices
    t_b = mod(t_b,3)+1  ! 1 -> 2 -> 3 -> 1
    t_c = mod(t_c,3)+1  ! 2 -> 3 -> 1 -> 2
    t_f = mod(t_f,3)+1  ! 3 -> 1 -> 2 -> 3


    !if the current input file is not yet finished, just increment stepf to 
    !  the next time step
    IF (((startfile .AND. (iint==0) .AND. (stepf==tdim)) .OR.   &
         (stepf .LT. tdim) )      .OR. filestep==0              ) THEN

      stepf=stepf+1
      counter=countfilenum
      write(*,*)'continue reading in previous hydro file from record ',stepf

    ELSE
    !if the current input file is finished, update filenm to next input file,
    !  and reset stepf to 1

      !Open netCDF file
      iint = iint+filestep
      counter=iint+filenum  !176 + 1 = 177 --> June 26,1995
      countfilenum=counter
      stepf = 1
      write(*,*)'opening new hydro file number',counter
    ENDIF


    !Get i/j max/min for rho/u/v
    call setijruv()
    FID=110 

    if(Zgrid)then
       nfilesin=8 ! MITgcm binary outputs without wind files
       if(Wind) nfilesin=nfilesin+2 
       if(WindIntensity) nfilesin=nfilesin+1
       !if(readNetcdfSwdown) nfilesin=nfilesin+1
    else 
       nfilesin=1 ! Roms NETcdf outputs
    endif
 
      !------------------------------------

      if(readZeta)then  
        call read_data_from_file(VAR_ID_zeta,vi,uj,1,1,1,1,1,romZf,RNODE,stepf,1,do_not_interpolate)
      else
        romZf = constZeta
      endif

       !Reshape input to fit node numbers assigned to elements
       do j=t_ijruv(JMIN,RNODE),t_ijruv(JMAX,RNODE)
         do i=t_ijruv(IMIN,RNODE),t_ijruv(IMAX,RNODE)
           count = (j-1)*vi + i
           t_zeta(t_f,count)     = romZf(i,j,1) * m_r(i,j,us_tridim)
         enddo
       enddo

      !------------------------------------

      if(readSalt)then
        call read_data_from_file(VAR_ID_salt,vi,uj,us,1,1,1,1,romSf,RNODE,stepf,1,do_not_interpolate)
      else
        romSf = constSalt
      endif

       !Reshape input to fit node numbers assigned to elements
       do j=t_ijruv(JMIN,RNODE),t_ijruv(JMAX,RNODE)
         do i=t_ijruv(IMIN,RNODE),t_ijruv(IMAX,RNODE)
           count = (j-1)*vi + i
           do k=1,us    
             kmask=min(k,us_tridim)
             t_salt(t_f,count,k) = romSf(i,j,k,1) * m_r(i,j,kmask)
          enddo
         enddo
       enddo

       ! COPY DATA OF FIRST CELL ABOVE THE BOTTOM TO ALL CELLS BELOW 
       if(Zgrid)then
        nodestocopy=1
        do j=t_ijruv(JMIN,RNODE),t_ijruv(JMAX,RNODE)
          do i=t_ijruv(IMIN,RNODE),t_ijruv(IMAX,RNODE)
            count = (j-1)*vi + i
            kbot=BottomK(i,j,1)
            if(kbot>1)then
              if(kbot<ws )then
               t_salt(t_f,count,kbot-1) =  t_salt(t_f,count,kbot)
              endif
              do searchnode=nodestocopy,NUM_COPNOD(RNODE)
                if(Node_COPNOD(RNODE,searchnode).le.count)nodestocopy=searchnode
                if(Node_COPNOD(RNODE,searchnode).ge.count)exit
              enddo
              do while(Node_COPNOD(RNODE,nodestocopy).eq.count)
                k1=klev_COPNOD(RNODE,nodestocopy,2)
                k2=klev_COPNOD(RNODE,nodestocopy,1)
                t_salt(t_f,count,k1:k2)=(                       &
                  t_salt(t_f,Nghb_COPNOD(RNODE,nodestocopy,1),k1:k2) &
                            *Coef_COPNOD(RNODE,nodestocopy,1) +    &
                  t_salt(t_f,Nghb_COPNOD(RNODE,nodestocopy,2),k1:k2) &
                            *Coef_COPNOD(RNODE,nodestocopy,2) +    &
                  t_salt(t_f,Nghb_COPNOD(RNODE,nodestocopy,3),k1:k2) &
                            *Coef_COPNOD(RNODE,nodestocopy,3) +    &
                  t_salt(t_f,Nghb_COPNOD(RNODE,nodestocopy,4),k1:k2) &
                              *Coef_COPNOD(RNODE,nodestocopy,4) ) 
                nodestocopy = nodestocopy +1
                if(nodestocopy.gt.NUM_COPNOD(RNODE))exit
              enddo
            endif
            if(nodestocopy.gt.NUM_COPNOD(RNODE))exit
          enddo    
          if(nodestocopy.gt.NUM_COPNOD(RNODE))exit
        enddo
       endif ! (Zgrid)      

      !------------------------------------

      if(readTemp)then  
        call read_data_from_file(VAR_ID_temp,vi,uj,us,1,1,1,1,romTf,RNODE,stepf,1,do_not_interpolate)
      else
        romTf = constTemp
      endif

       !Reshape input to fit node numbers assigned to elements
       do j=t_ijruv(JMIN,RNODE),t_ijruv(JMAX,RNODE)
         do i=t_ijruv(IMIN,RNODE),t_ijruv(IMAX,RNODE)
           count = (j-1)*vi + i
           do k=1,us    
             kmask=min(k,us_tridim)
             t_temp(t_f,count,k) = romTf(i,j,k,1) * m_r(i,j,kmask)
          enddo
         enddo
       enddo

       ! COPY DATA OF FIRST CELL ABOVE THE BOTTOM TO ALL CELLS BELOW 
       if(Zgrid)then
        nodestocopy=1
        do j=t_ijruv(JMIN,RNODE),t_ijruv(JMAX,RNODE)
          do i=t_ijruv(IMIN,RNODE),t_ijruv(IMAX,RNODE)
            count = (j-1)*vi + i
            kbot=BottomK(i,j,1)
            if(kbot>1)then
              if(kbot<ws )then
               t_temp(t_f,count,kbot-1) =  t_temp(t_f,count,kbot)
              endif
              do searchnode=nodestocopy,NUM_COPNOD(RNODE)
                if(Node_COPNOD(RNODE,searchnode).le.count)nodestocopy=searchnode
                if(Node_COPNOD(RNODE,searchnode).ge.count)exit
              enddo
              do while(Node_COPNOD(RNODE,nodestocopy).eq.count)
                k1=klev_COPNOD(RNODE,nodestocopy,2)
                k2=klev_COPNOD(RNODE,nodestocopy,1)
                t_temp(t_f,count,k1:k2)=(                       &
                  t_temp(t_f,Nghb_COPNOD(RNODE,nodestocopy,1),k1:k2) &
                            *Coef_COPNOD(RNODE,nodestocopy,1) +    &
                  t_temp(t_f,Nghb_COPNOD(RNODE,nodestocopy,2),k1:k2) &
                            *Coef_COPNOD(RNODE,nodestocopy,2) +    &
                  t_temp(t_f,Nghb_COPNOD(RNODE,nodestocopy,3),k1:k2) &
                            *Coef_COPNOD(RNODE,nodestocopy,3) +    &
                  t_temp(t_f,Nghb_COPNOD(RNODE,nodestocopy,4),k1:k2) &
                              *Coef_COPNOD(RNODE,nodestocopy,4) ) 
                nodestocopy = nodestocopy +1
                if(nodestocopy.gt.NUM_COPNOD(RNODE))exit
              enddo
            endif
            if(nodestocopy.gt.NUM_COPNOD(RNODE))exit
          enddo    
          if(nodestocopy.gt.NUM_COPNOD(RNODE))exit
        enddo
       endif ! (Zgrid)      

      !------------------------------------

      if(readDens)then  
        call read_data_from_file(VAR_ID_den,vi,uj,us,1,1,1,1,romDf,RNODE,stepf,1,do_not_interpolate)
      else
        romDf = constDens
      endif

       !Reshape input to fit node numbers assigned to elements
       do j=t_ijruv(JMIN,RNODE),t_ijruv(JMAX,RNODE)
         do i=t_ijruv(IMIN,RNODE),t_ijruv(IMAX,RNODE)
           count = (j-1)*vi + i
           do k=1,us    
             kmask=min(k,us_tridim)
             t_den(t_f,count,k)  = (romDf(i,j,k,1) + DBLE(1000.0)) * m_r(i,j,kmask)
          enddo
         enddo
       enddo


       ! COPY DATA OF FIRST CELL ABOVE THE BOTTOM TO ALL CELLS BELOW 
       if(Zgrid)then
        nodestocopy=1
        do j=t_ijruv(JMIN,RNODE),t_ijruv(JMAX,RNODE)
          do i=t_ijruv(IMIN,RNODE),t_ijruv(IMAX,RNODE)
            count = (j-1)*vi + i
            kbot=BottomK(i,j,1)
            if(kbot>1)then
              if(kbot<ws )then
               t_den(t_f,count,kbot-1) =  t_den(t_f,count,kbot)
              endif
              do searchnode=nodestocopy,NUM_COPNOD(RNODE)
                if(Node_COPNOD(RNODE,searchnode).le.count)nodestocopy=searchnode
                if(Node_COPNOD(RNODE,searchnode).ge.count)exit
              enddo
              do while(Node_COPNOD(RNODE,nodestocopy).eq.count)
                k1=klev_COPNOD(RNODE,nodestocopy,2)
                k2=klev_COPNOD(RNODE,nodestocopy,1)
                t_den(t_f,count,k1:k2)=(                       &
                  t_den(t_f,Nghb_COPNOD(RNODE,nodestocopy,1),k1:k2) &
                           *Coef_COPNOD(RNODE,nodestocopy,1) +    &
                  t_den(t_f,Nghb_COPNOD(RNODE,nodestocopy,2),k1:k2) &
                           *Coef_COPNOD(RNODE,nodestocopy,2) +    &
                  t_den(t_f,Nghb_COPNOD(RNODE,nodestocopy,3),k1:k2) &
                           *Coef_COPNOD(RNODE,nodestocopy,3) +    &
                  t_den(t_f,Nghb_COPNOD(RNODE,nodestocopy,4),k1:k2) &
                           *Coef_COPNOD(RNODE,nodestocopy,4) ) 
                nodestocopy = nodestocopy +1
                if(nodestocopy.gt.NUM_COPNOD(RNODE))exit
              enddo
            endif
            if(nodestocopy.gt.NUM_COPNOD(RNODE))exit
          enddo    
          if(nodestocopy.gt.NUM_COPNOD(RNODE))exit
        enddo
       endif ! (Zgrid)      
 
      !------------------------------------

      if(readU)then  
        call read_data_from_file(VAR_ID_uvel,ui,uj,us,1,1,1,1,romUf,UNODE,stepf,1,do_not_interpolate)
      else
        romUf = constU
      endif

       do j=t_ijruv(JMIN,UNODE),t_ijruv(JMAX,UNODE)
         do i=t_ijruv(IMIN,UNODE),t_ijruv(IMAX,UNODE)
           count = (j-1)*ui + i
           do k=1,us    
             t_Uvel(t_f,count,k) = romUf(i,j,k,1) * m_u(i,j,min(k,us_tridim))
           enddo
         enddo
       enddo


       ! COPY DATA OF FIRST CELL ABOVE THE BOTTOM TO ALL CELLS BELOW 
       if(Zgrid)then
        nodestocopy=1
        do j=t_ijruv(JMIN,UNODE),t_ijruv(JMAX,UNODE)
          do i=t_ijruv(IMIN,UNODE),t_ijruv(IMAX,UNODE)
            count = (j-1)*ui + i
            kbot=BottomK(i,j,2)
            if(kbot>1)then
              if(kbot<ws )then
               t_Uvel(t_f,count,kbot-1) =  t_Uvel(t_f,count,kbot)
              endif
              do searchnode=nodestocopy,NUM_COPNOD(UNODE)
                if(Node_COPNOD(UNODE,searchnode).le.count)nodestocopy=searchnode
                if(Node_COPNOD(UNODE,searchnode).ge.count)exit
              enddo
              do while(Node_COPNOD(UNODE,nodestocopy).eq.count)
                k1=klev_COPNOD(UNODE,nodestocopy,2)
                k2=klev_COPNOD(UNODE,nodestocopy,1)
                  t_Uvel(t_f,count,k1:k2)=(                       &
                  t_Uvel(t_f,Nghb_COPNOD(UNODE,nodestocopy,1),k1:k2) &
                            *Coef_COPNOD(UNODE,nodestocopy,1) +    &
                  t_Uvel(t_f,Nghb_COPNOD(UNODE,nodestocopy,2),k1:k2) &
                            *Coef_COPNOD(UNODE,nodestocopy,2) +    &
                  t_Uvel(t_f,Nghb_COPNOD(UNODE,nodestocopy,3),k1:k2) &
                            *Coef_COPNOD(UNODE,nodestocopy,3) +    &
                  t_Uvel(t_f,Nghb_COPNOD(UNODE,nodestocopy,4),k1:k2) &
                            *Coef_COPNOD(UNODE,nodestocopy,4) ) 
                nodestocopy = nodestocopy +1
                if(nodestocopy.gt.NUM_COPNOD(UNODE))exit
              enddo
            endif
            if(nodestocopy.gt.NUM_COPNOD(UNODE))exit
          enddo    
          if(nodestocopy.gt.NUM_COPNOD(UNODE))exit
        enddo
       endif ! (Zgrid)      

      !------------------------------------
      if(readV)then  
          call read_data_from_file(VAR_ID_vvel,vi,vj,us,1,1,1,1,romVf,VNODE,stepf,1,do_not_interpolate)
      else
        romVf = constV
      endif

       do j=t_ijruv(JMIN,VNODE),t_ijruv(JMAX,VNODE)
         do i=t_ijruv(IMIN,VNODE),t_ijruv(IMAX,VNODE)
           count = (j-1)*vi + i
           do k=1,us    
             t_Vvel(t_f,count,k) = romVf(i,j,k,1) * m_v(i,j,min(k,us_tridim))
           enddo
         enddo    
       enddo
       ! COPY DATA OF FIRST CELL ABOVE THE BOTTOM TO ALL CELLS BELOW 
       if(Zgrid)then
        nodestocopy=1
        do j=t_ijruv(JMIN,VNODE),t_ijruv(JMAX,VNODE)
          do i=t_ijruv(IMIN,VNODE),t_ijruv(IMAX,VNODE)
            count = (j-1)*vi + i
            kbot=BottomK(i,j,3)
            if(kbot>1)then
              if(kbot<ws )then
               t_Vvel(t_f,count,kbot-1) =  t_Vvel(t_f,count,kbot)
              endif
              do searchnode=nodestocopy,NUM_COPNOD(VNODE)
                if(Node_COPNOD(VNODE,searchnode).le.count)nodestocopy=searchnode
                if(Node_COPNOD(VNODE,searchnode).ge.count)exit
              enddo
              do while(Node_COPNOD(VNODE,nodestocopy).eq.count)
                k1=klev_COPNOD(VNODE,nodestocopy,2)
                k2=klev_COPNOD(VNODE,nodestocopy,1)
                  t_Vvel(t_f,count,k1:k2)=(                       &
                  t_Vvel(t_f,Nghb_COPNOD(VNODE,nodestocopy,1),k1:k2) &
                            *Coef_COPNOD(VNODE,nodestocopy,1) +    &
                  t_Vvel(t_f,Nghb_COPNOD(VNODE,nodestocopy,2),k1:k2) &
                            *Coef_COPNOD(VNODE,nodestocopy,2) +    &
                  t_Vvel(t_f,Nghb_COPNOD(VNODE,nodestocopy,3),k1:k2) &
                            *Coef_COPNOD(VNODE,nodestocopy,3) +    &
                  t_Vvel(t_f,Nghb_COPNOD(VNODE,nodestocopy,4),k1:k2) &
                            *Coef_COPNOD(VNODE,nodestocopy,4) ) 
                nodestocopy = nodestocopy +1
                if(nodestocopy.gt.NUM_COPNOD(VNODE))exit
              enddo
            endif
            if(nodestocopy.gt.NUM_COPNOD(VNODE))exit
          enddo    
          if(nodestocopy.gt.NUM_COPNOD(VNODE))exit
        enddo

        endif ! (Zgrid)

      !------------------------------------

      if(readW)then  
        call read_data_from_file(VAR_ID_wvel,vi,uj,ws,1,1,1,1,romWf,RNODE,stepf,1,do_not_interpolate)
      else
        romWf = constW
      endif

       !Reshape input to fit node numbers assigned to elements
       do j=t_ijruv(JMIN,RNODE),t_ijruv(JMAX,RNODE)
         do i=t_ijruv(IMIN,RNODE),t_ijruv(IMAX,RNODE)
           count = (j-1)*vi + i
           do k=1,us    
             kmask=min(k,us_tridim)
             t_Wvel(t_f,count,k+1) = romWf(i,j,k+1,1) * m_r(i,j,kmask)
          enddo
           t_Wvel(t_f,count,1)  = romWf(i,j,1,1) * m_r(i,j,1)
         enddo
       enddo

       if (Zgrid)then
         nodestocopy=1
         do j=t_ijruv(JMIN,RNODE),t_ijruv(JMAX,RNODE)
          do i=t_ijruv(IMIN,RNODE),t_ijruv(IMAX,RNODE)
            count = (j-1)*vi + i
            kbot=BottomK(i,j,1)
            if(kbot>1)then
             if(kbot<ws)then
              t_Wvel(t_f,count,kbot-1) = t_Wvel(t_f,count,kbot)
             endif
             do searchnode=nodestocopy,NUM_COPNOD(RNODE)
               if(Node_COPNOD(RNODE,searchnode).le.count)nodestocopy=searchnode
               if(Node_COPNOD(RNODE,searchnode).ge.count)exit
             enddo
             do while(Node_COPNOD(RNODE,nodestocopy).eq.count)
               k1=klev_COPNOD(RNODE,nodestocopy,2)
               k2=klev_COPNOD(RNODE,nodestocopy,1)
                 t_Wvel(t_f,count,k1:k2)=(                        &
                 t_Wvel(t_f,Nghb_COPNOD(RNODE,nodestocopy,1),k1:k2) &
                           *Coef_COPNOD(RNODE,nodestocopy,1) +  &
                 t_Wvel(t_f,Nghb_COPNOD(RNODE,nodestocopy,2),k1:k2) &
                           *Coef_COPNOD(RNODE,nodestocopy,2) +  &
                 t_Wvel(t_f,Nghb_COPNOD(RNODE,nodestocopy,3),k1:k2) &
                           *Coef_COPNOD(RNODE,nodestocopy,3) +  &
                 t_Wvel(t_f,Nghb_COPNOD(RNODE,nodestocopy,4),k1:k2) &
                           *Coef_COPNOD(RNODE,nodestocopy,4) ) 
               nodestocopy = nodestocopy +1
               if(nodestocopy.gt.NUM_COPNOD(RNODE))exit
             enddo !while
            endif
            if(nodestocopy.gt.NUM_COPNOD(RNODE))exit
          enddo    
          if(nodestocopy.gt.NUM_COPNOD(RNODE))exit
         enddo      
       endif ! (Zgrid)

      !------------------------------------

      if(readAks)then  
        call read_data_from_file(VAR_ID_kh,vi,uj,us,1,1,1,1,romKHf,RNODE,stepf,1,do_not_interpolate)
      else
        romKHf = constAks
      endif

       !Reshape input to fit node numbers assigned to elements
       do j=t_ijruv(JMIN,RNODE),t_ijruv(JMAX,RNODE)
         do i=t_ijruv(IMIN,RNODE),t_ijruv(IMAX,RNODE)
           count = (j-1)*vi + i
           do k=1,us    
             kmask=min(k,us_tridim)
             t_KH(t_f,count,k+1)   = romKHf(i,j,k+1,1) * m_r(i,j,kmask)
          enddo
           t_KH(t_f,count,1)    = romKHf(i,j,1,1) * m_r(i,j,1)
         enddo
       enddo

       ! COPY DATA OF FIRST CELL ABOVE THE BOTTOM TO ALL CELLS BELOW 
       if(Zgrid)then
        nodestocopy=1
        do j=t_ijruv(JMIN,RNODE),t_ijruv(JMAX,RNODE)
          do i=t_ijruv(IMIN,RNODE),t_ijruv(IMAX,RNODE)
            count = (j-1)*vi + i
            kbot=BottomK(i,j,1)
            if(kbot>1)then
              if(kbot<ws )then
               t_KH(t_f,count,kbot-1) =  t_KH(t_f,count,kbot)
              endif
              do searchnode=nodestocopy,NUM_COPNOD(RNODE)
                if(Node_COPNOD(RNODE,searchnode).le.count)nodestocopy=searchnode
                if(Node_COPNOD(RNODE,searchnode).ge.count)exit
              enddo
              do while(Node_COPNOD(RNODE,nodestocopy).eq.count)
                k1=klev_COPNOD(RNODE,nodestocopy,2)
                k2=klev_COPNOD(RNODE,nodestocopy,1)
                t_KH(t_f,count,k1:k2)=(                       &
                  t_KH(t_f,Nghb_COPNOD(RNODE,nodestocopy,1),k1:k2) &
                          *Coef_COPNOD(RNODE,nodestocopy,1) +    &
                  t_KH(t_f,Nghb_COPNOD(RNODE,nodestocopy,2),k1:k2) &
                          *Coef_COPNOD(RNODE,nodestocopy,2) +    &
                  t_KH(t_f,Nghb_COPNOD(RNODE,nodestocopy,3),k1:k2) &
                          *Coef_COPNOD(RNODE,nodestocopy,3) +    &
                  t_KH(t_f,Nghb_COPNOD(RNODE,nodestocopy,4),k1:k2) &
                          *Coef_COPNOD(RNODE,nodestocopy,4) ) 
                nodestocopy = nodestocopy +1
                if(nodestocopy.gt.NUM_COPNOD(RNODE))exit
              enddo
            endif
            if(nodestocopy.gt.NUM_COPNOD(RNODE))exit
          enddo    
          if(nodestocopy.gt.NUM_COPNOD(RNODE))exit
        enddo
       endif ! (Zgrid)      

      !------------------------------------   

      if(Wind)then
       if(readUwind)then  
          call read_data_from_file(VAR_ID_uwind,ui,uj,1,1,1,1,1,modelUwindf,UNODE,stepf,1,interpol_from_cell_center_to_CArakawa)
       else
         modelUwindf = constUwind
       endif

       do j=t_ijruv(JMIN,UNODE),t_ijruv(JMAX,UNODE)
        do i=t_ijruv(IMIN,UNODE),t_ijruv(IMAX,UNODE)
          count = (j-1)*ui + i
          t_uwind(t_f,count) =    modelUwindf(i,j,1) *    m_u(i,j,us_tridim)  !---CL-OGS: replace m_r by m_u
        enddo
       enddo

      endif

      !------------------------------------

      if(Wind)then
       if(readVwind)then  
          call read_data_from_file(VAR_ID_vwind,vi,vj,1,1,1,1,1,modelVwindf,VNODE,stepf,1,interpol_from_cell_center_to_CArakawa)
       else
        modelVwindf = constVwind
       endif

       do j=t_ijruv(JMIN,VNODE),t_ijruv(JMAX,VNODE)
        do i=t_ijruv(IMIN,VNODE),t_ijruv(IMAX,VNODE)
          count = (j-1)*vi + i
          t_vwind(t_f,count) =    modelVwindf(i,j,1) *    m_v(i,j,us_tridim)  !---CL-OGS: replace m_r by m_v
        enddo
       enddo
      endif

      !------------------------------------

      if(readIwind)then  
       if(Zgrid)then
         call read_data_from_file(VAR_ID_iwind,vi,uj,1,1,1,1,1,modelIwindf,RNODE,stepf,1,do_not_interpolate)
       else
          write(*,*) ' ERROR Wind intensity not present in Roms files'
          write(*,*) ' setting modelIwindf = constIwind=',constIwind
          modelIwindf = constIwind
       endif
      else
       modelIwindf = constIwind
      endif
      if(WindIntensity .and. Zgrid)then
       do j=t_ijruv(JMIN,RNODE),t_ijruv(JMAX,RNODE)
        do i=t_ijruv(IMIN,RNODE),t_ijruv(IMAX,RNODE)
          count = (j-1)*vi + i
          t_iwind(t_f,count) =    modelIwindf(i,j,1) *    m_r(i,j,us_tridim)
        enddo
       enddo
      endif

      !------------------------------------

 
      ! Store the ranges of nodes that were updated
      updatenodesbuffer=0
      ! rho node range
      do j=t_ijruv(JMIN,RNODE),t_ijruv(JMAX,RNODE)
        updatenodesbuffer(1,j,1) = (j-1)*vi + t_ijruv(IMIN,RNODE) ! frst rnode at latitude j
        updatenodesbuffer(2,j,1) = (j-1)*vi + t_ijruv(IMAX,RNODE) ! last rnode at latitude j
      enddo
      ! u node range
      do j=t_ijruv(JMIN,UNODE),t_ijruv(JMAX,UNODE)
        updatenodesbuffer(1,j,2) = (j-1)*ui + t_ijruv(IMIN,UNODE) ! frst unode at latitude j
        updatenodesbuffer(2,j,2) = (j-1)*ui + t_ijruv(IMAX,UNODE) ! last unode at latitude j
      enddo
      ! v node range
      do j=t_ijruv(JMIN,VNODE),t_ijruv(JMAX,VNODE)
        updatenodesbuffer(1,j,3) = (j-1)*vi + t_ijruv(IMIN,VNODE)  ! frst vnode at latitude j
        updatenodesbuffer(2,j,3) = (j-1)*vi + t_ijruv(IMAX,VNODE) ! last vnode at latitude j
      enddo
   

    !  ***    IMIOM *****
    ! WIND WAVE MODEL DATA  ------------------------------------
    IF(OilOn)then
     write(*,*)'OilOn'
     if(WindWaveModel .and. (.not. Zgrid))THEN
          write(*,*)'WindWaveModel'
          nvarf=1
          scounter = iint + swan_filenum

          call set_filename(VAR_ID_swan,scounter,swannm)
         
          ! Read in data for first three external time steps
          STATUS = NF90_OPEN(TRIM(swannm), NF90_NOWRITE, NCID)
          if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem NF90_OPEN'
          if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)

              ! **** Hsig ****
              startz(1)=t_ijruv(IMIN,RNODE)
              startz(2)=t_ijruv(JMIN,RNODE)
              startz(3)=stepf

              countz(1)=t_ijruv(IMAX,RNODE)-t_ijruv(IMIN,RNODE)+1
              countz(2)=t_ijruv(JMAX,RNODE)-t_ijruv(JMIN,RNODE)+1
              countr(3)=1

              STATUS = NF90_INQ_VARID(NCID,'Hs',VID)
              if (STATUS .NE. NF90_NOERR) then
                write(*,*) 'Problem find Hs'
                write(*,*) NF90_STRERROR(STATUS)
                stop
              endif

              STATUS = NF90_GET_VAR(NCID,VID,swanHsf(t_ijruv(IMIN,RNODE):t_ijruv(IMAX,RNODE),    &
                                    t_ijruv(JMIN,RNODE):t_ijruv(JMAX,RNODE),:),STARTz,COUNTz)

              if (STATUS .NE. NF90_NOERR) then
                write(*,*) 'Problem read SwanHs array'
                write(*,*) NF90_STRERROR(STATUS)
                stop
              endif

              ! **** tm01 ****
              startz(1)=t_ijruv(IMIN,RNODE)
              startz(2)=t_ijruv(JMIN,RNODE)
              startz(3)=stepf

              countz(1)=t_ijruv(IMAX,RNODE)-t_ijruv(IMIN,RNODE)+1
              countz(2)=t_ijruv(JMAX,RNODE)-t_ijruv(JMIN,RNODE)+1
              countr(3)=1

              STATUS = NF90_INQ_VARID(NCID,'tm01',VID)
              if (STATUS .NE. NF90_NOERR) then
                write(*,*) 'Problem find tm01'
                write(*,*) NF90_STRERROR(STATUS)
                stop
              endif

              STATUS = NF90_GET_VAR(NCID,VID,swantm01f(t_ijruv(IMIN,RNODE):t_ijruv(IMAX,RNODE),  &
                                    t_ijruv(JMIN,RNODE):t_ijruv(JMAX,RNODE),:),STARTz,COUNTz)

              if (STATUS .NE. NF90_NOERR) then
                write(*,*) 'Problem read swantm01 array'
                write(*,*) NF90_STRERROR(STATUS)
                stop
              endif

              ! **** u10 ****
              startr(1)=t_ijruv(IMIN,UNODE)
              startr(2)=t_ijruv(JMIN,UNODE)
              startr(3)=stepf

              countr(1)=t_ijruv(IMAX,UNODE)-t_ijruv(IMIN,UNODE)+1
              countr(2)=t_ijruv(JMAX,UNODE)-t_ijruv(JMIN,UNODE)+1
              countr(3)=1

              STATUS = NF90_INQ_VARID(NCID,'u10',VID)
              if (STATUS .NE. NF90_NOERR) then
                write(*,*) 'Problem find u10'
                write(*,*) NF90_STRERROR(STATUS)
                stop
              endif

              STATUS = NF90_GET_VAR(NCID,VID,modelUwindf(t_ijruv(IMIN,UNODE):t_ijruv(IMAX,UNODE),&
                                    t_ijruv(JMIN,UNODE):t_ijruv(JMAX,UNODE),:),STARTz,COUNTz)

              if (STATUS .NE. NF90_NOERR) then
                write(*,*) 'Problem read swanU array'
                write(*,*) NF90_STRERROR(STATUS)
                stop
              endif

              ! **** V10 ****
              startr(1)=t_ijruv(IMIN,VNODE)
              startr(2)=t_ijruv(JMIN,VNODE)
              startz(3)=stepf

              countr(1)=t_ijruv(IMAX,VNODE)-t_ijruv(IMIN,VNODE)+1
              countr(2)=t_ijruv(JMAX,VNODE)-t_ijruv(JMIN,VNODE)+1
              countr(3)=1

              STATUS = NF90_INQ_VARID(NCID,'v10',VID)
              if (STATUS .NE. NF90_NOERR) then
                write(*,*) 'Problem find v10'
                write(*,*) NF90_STRERROR(STATUS)
                stop
              endif

              STATUS=NF90_GET_VAR(NCID,VID,modelVwindf(t_ijruv(IMIN,VNODE):t_ijruv(IMAX,VNODE), &
                                    t_ijruv(JMIN,VNODE):t_ijruv(JMAX,VNODE),:),STARTz,COUNTz)

              if (STATUS .NE. NF90_NOERR) then
                write(*,*) 'Problem read swanV array'
                write(*,*) NF90_STRERROR(STATUS)
                stop
              endif

              ! **** PeakDir ****
              startz(1)=t_ijruv(IMIN,RNODE)
              startz(2)=t_ijruv(JMIN,RNODE)
              startz(3)=stepf

              countz(1)=t_ijruv(IMAX,RNODE)-t_ijruv(IMIN,RNODE)+1
              countz(2)=t_ijruv(JMAX,RNODE)-t_ijruv(JMIN,RNODE)+1
              countr(3)=1

              STATUS = NF90_INQ_VARID(NCID,'Pd',VID)
              if (STATUS .NE. NF90_NOERR) then
                write(*,*) 'Problem find Pd'
                write(*,*) NF90_STRERROR(STATUS)
                stop
              endif

              STATUS = NF90_GET_VAR(NCID,VID,swanpdf(t_ijruv(IMIN,RNODE):t_ijruv(IMAX,RNODE),    &
                                    t_ijruv(JMIN,RNODE):t_ijruv(JMAX,RNODE),:),STARTz,COUNTz)

              if (STATUS .NE. NF90_NOERR) then
                write(*,*) 'Problem read swanpd array'
                write(*,*) NF90_STRERROR(STATUS)
                stop
              endif

              ! **** PeakWaveLength ****
              startz(1)=t_ijruv(IMIN,RNODE)
              startz(2)=t_ijruv(JMIN,RNODE)
              startz(3)=stepf

              countz(1)=t_ijruv(IMAX,RNODE)-t_ijruv(IMIN,RNODE)+1
              countz(2)=t_ijruv(JMAX,RNODE)-t_ijruv(JMIN,RNODE)+1
              countr(3)=1

              STATUS = NF90_INQ_VARID(NCID,'Pwl',VID)
              if (STATUS .NE. NF90_NOERR) then
                write(*,*) 'Problem find Pwl'
                write(*,*) NF90_STRERROR(STATUS)
                stop
              endif

              STATUS = NF90_GET_VAR(NCID,VID,swanwlf(t_ijruv(IMIN,RNODE):t_ijruv(IMAX,RNODE),    &
                                    t_ijruv(JMIN,RNODE):t_ijruv(JMAX,RNODE),:),STARTz,COUNTz)

              if (STATUS .NE. NF90_NOERR) then
                write(*,*) 'Problem read swanwl array'
                write(*,*) NF90_STRERROR(STATUS)
                stop
              endif

          !close the dataset and reassign the NCID
          STATUS = NF90_CLOSE(NCID)

      else   ! if WindWaveModel
            write(*,*)'not WindWaveModel or Zgrid',SigWaveHeight
            swanHsf   = SigWaveHeight
            write(*,*)'SigWaveHeight passed'
            swantm01f = MeanWavePeriod
            if(.not.readUwind)  modelUwindf = UWind_10
            if(.not.readVwind)  modelVwindf = VWind_10
            swanpdf   = PeakDirection
            swanwlf   = PeakWaveLength
      endif                ! if WindWaveModel

      !Reshape input to fit node numbers assigned to elements
       do j=t_ijruv(JMIN,RNODE),t_ijruv(JMAX,RNODE)
          do i=t_ijruv(IMIN,RNODE),t_ijruv(IMAX,RNODE)
            count = (j-1)*vi + i
            t_hsig(t_f,count) = swanHsf(i,j,1) * m_r(i,j,us_tridim)
            t_tm01(t_f,count) = swantm01f(i,j,1) * m_r(i,j,us_tridim)
            t_pdir(t_f,count) = swanpdf(i,j,1) * m_r(i,j,us_tridim)
            t_wlen(t_f,count) = swanwlf(i,j,1) * m_r(i,j,us_tridim)
          enddo
       enddo
      if(windwavemodel.or.((.not.readUwind) .and. (.not. Zgrid)))then     !only overwrite ROMS t_uwind (from sustr) if using windwavesmodel OR overwriting by constant U wind
            write(*,*)'overwrite ROMS U Wind'
              do j=t_ijruv(JMIN,UNODE),t_ijruv(JMAX,UNODE)
                do i=t_ijruv(IMIN,UNODE),t_ijruv(IMAX,UNODE)
                  count = (j-1)*ui + i
                  t_uwind(t_f,count) = modelUwindf(i,j,1) * m_u(i,j,us_tridim)
                enddo
              enddo
      end if
      if(windwavemodel.or.((.not.readVwind) .and. (.not. Zgrid)))then    !only overwrite ROMS t_vwind (from sustr) if using windwavesmodel OR overwriting by constant V wind
           write(*,*)'overwrite ROMS V Wind'
              do j=t_ijruv(JMIN,VNODE),t_ijruv(JMAX,VNODE)
                do i=t_ijruv(IMIN,VNODE),t_ijruv(IMAX,VNODE)
                  count = (j-1)*vi + i
              t_vwind(t_f,count) = modelVwindf(i,j,1) * m_v(i,j,us_tridim)
                enddo
              enddo
      end if
    END IF        ! if OilOn


     !DEALLOCATE SUBROUTINE VARIABLES
     if(OilOn)then
       DEALLOCATE(swanHsf,swantm01f,swanpdf,swanwlf)
     END IF        !if OilOn




    !DEALLOCATE SUBROUTINE VARIABLES
    DEALLOCATE(romZf,romSf,romTf,romUf,romVf,romWf,romKHf)
    DEALLOCATE(modelUwindf,modelVwindf,modelIwindf)

  END SUBROUTINE updateHydro


  SUBROUTINE setEle(Xpar,Ypar,Zpar,n,it,num,err,first)
    !This Subroutine determines which Rho, U, and V grid elements contain 
    !  the given particle
    USE PARAM_MOD, ONLY: numpar,ui,vi,us,ws,vj,uj,Zgrid 
    USE GRIDCELL_MOD, ONLY: gridcell
    USE CONVERT_MOD, ONLY: x2lon,y2lat
    !$ use OMP_LIB          
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: Xpar,Ypar,Zpar
    INTEGER, INTENT(IN) :: n,num,it
    INTEGER, INTENT(INOUT), OPTIONAL :: err
    LOGICAL, INTENT(IN), OPTIONAL :: first

    LOGICAL :: fst,same_vertical_level
    INTEGER :: i,triangle,checkele,P_r_ele,P_u_ele,P_v_ele,oP_ele,P_ele,error,k
    INTEGER :: updatednode ,j,err_in
    INTEGER :: rank
     rank=1
     !$ rank=OMP_GET_THREAD_NUM () +1
    err_in=0
    if(PRESENT(err)) then
       err_in= err
       err = 0
    endif

    error = 0
    if(Zgrid)then
      k=getKRlevel(Zpar)
    else
      k=1
    endif
    triangle=-1
    P_r_ele=-1
    P_u_ele=-1
    P_v_ele=-1
    !rite(*,'(a,i12,a,i1,a,i6,4(a,f13.4),2(a,f9.4),6(a,i7),4(a,f10.6),2(a,f9.4),4(a,i3),a)')&
    !'setEle warning it=',it,' at call num',num,' for part ',n, &
    !' of pos (X,Y,Z)=(',Xpar_at_setEle(n),&
    !' ->',Xpar,&
    !' ;',Ypar_at_setEle(n),&
    !' ->',Ypar,&
    !' ; ',Zpar_at_setEle(n),&
    !' ->',Zpar, &
    !' ) ; P_r,u,v=(',P_r_element(n),' ->',P_r_ele,' ; ',P_u_element(n),' ->',P_u_ele,' ; ',P_v_element(n),' ->',P_v_ele, &
    !' ) lon,lat,depth,k=( ',x2lon(Xpar_at_setEle(n),Ypar_at_setEle(n)), &
    !' ->',x2lon(Xpar,Ypar),  &
    !' ; ',y2lat(Ypar_at_setEle(n)),&
    !' ->',y2lat(Ypar),&
    !' ; ', Zpar_at_setEle(n),&
    !' ->',Zpar,&
    !' ; ',P_klev_old(n), &
    !' ->',k,' ) and triangle,error=(',triangle,' ; ',error, &
    !' ) ENTERING setEle ROUTINE'

    if(Xpar/=Xpar .or. Ypar/=Ypar .or. Zpar/=Zpar)then
     write(*,'(a,i12,a,i1,a,i6,4(a,f13.4),2(a,f9.4),6(a,i7),4(a,f10.6),2(a,f9.4),4(a,i3),a)')&
      'setEle warning it=',it,' at call num',num,' for part ',n, &
      ' of pos (X,Y,Z)=(',Xpar_at_setEle(n),&
      ' ->',Xpar,&
      ' ;',Ypar_at_setEle(n),&
      ' ->',Ypar,&
      ' ; ',Zpar_at_setEle(n),&
      ' ->',Zpar, &
      ' ) ; P_r,u,v=(',P_r_element(n),' ->',P_r_ele,' ; ',P_u_element(n),' ->',P_u_ele,' ; ',P_v_element(n),' ->',P_v_ele, &
      ' ) lon,lat,depth,k=( ',x2lon(Xpar_at_setEle(n),Ypar_at_setEle(n)), &
      ' ->',x2lon(Xpar,Ypar),  &
      ' ; ',y2lat(Ypar_at_setEle(n)),&
      ' ->',y2lat(Ypar),&
      ' ; ', Zpar_at_setEle(n),&
      ' ->',Zpar,&
      ' ; ',P_klev_old(n), &
      ' ->',k,' ) and triangle,error=(',triangle,' ; ',error, &
      ' ) ETYPE=invalid initial Xpar,Ypar,Zpar'
      error=-2
      return
    endif 

    P_klev(n)=k
    if( PRESENT(first) ) then
      fst = first
    else
      fst = .FALSE.
    endif
    if((P_r_element(n).eq.0.or.P_u_element(n).eq.0).or.P_v_element(n).eq.0)then
      fst=.True.
    endif

    
    if((.not.fst) .and. (P_klev(n).ne.P_klev_old(n)))then
      same_vertical_level=.False.
    else
      same_vertical_level=.True.
    endif 

     IF(P_r_element(n).eq.0.or.P_u_element(n).eq.0.or.P_v_element(n).eq.0)then
     write(*,*)'----------------------------------------------------------'
     write(*,'(a,i12,a,i1,a,i6,4(a,f13.4),2(a,f9.4),6(a,i7),4(a,f10.6),2(a,f9.4),4(a,i3),a)')&
      'setEle warning it=',it,' at call num',num,' for part ',n, &
      ' of pos (X,Y,Z)=(',Xpar_at_setEle(n),&
      ' ->',Xpar,&
      ' ;',Ypar_at_setEle(n),&
      ' ->',Ypar,&
      ' ; ',Zpar_at_setEle(n),&
      ' ->',Zpar, &
      ' ) ; P_r,u,v=(',P_r_element(n),' ->',P_r_ele,' ; ',P_u_element(n),' ->',P_u_ele,' ; ',P_v_element(n),' ->',P_v_ele, &
      ' ) lon,lat,depth,k=( ',x2lon(Xpar_at_setEle(n),Ypar_at_setEle(n)), &
      ' ->',x2lon(Xpar,Ypar),  &
      ' ; ',y2lat(Ypar_at_setEle(n)),&
      ' ->',y2lat(Ypar),&
      ' ; ', Zpar_at_setEle(n),&
      ' ->',Zpar,&
      ' ; ',P_klev_old(n), &
      ' ->',k,' ) and triangle,error=(',triangle,' ; ',error, &
      ' ) Invalid initial element, now searching among all elements '
      fst=.True. 
    ENDIF

    !--------------Find r element in which particle is located------------------
    IF(.not.fst)THEN !if not the first iteration 
      !Find rho element in which particle is located
      oP_ele = P_r_element(n)
      do i=1,10
        triangle = 0
        if(r_Adjacent(oP_ele,i,k).NE.0) then
           checkele = r_Adjacent(oP_ele,i,k)
           call gridcell(rho_kwele(k),r_kwele_y(1,1,k),           &
                      r_kwele_x(1,1,k), &
                    Xpar,Ypar,P_ele,triangle,checkele)
        endif
        if(triangle .NE. 0) then
          P_r_element(n) = P_ele
          exit
        endif
      enddo !r_singlecellloop
      if(triangle.EQ.0)error = 4

      if (error.eq.4 .and. same_vertical_level) then
     write(*,'(a,i12,a,i1,a,i6,4(a,f13.4),2(a,f9.4),6(a,i7),4(a,f10.6),2(a,f9.4),4(a,i3),a)')&
      'setEle warning it=',it,' at call num',num,' for part ',n, &
      ' of pos (X,Y,Z)=(',Xpar_at_setEle(n),&
      ' ->',Xpar,&
      ' ;',Ypar_at_setEle(n),&
      ' ->',Ypar,&
      ' ; ',Zpar_at_setEle(n),&
      ' ->',Zpar, &
      ' ) ; P_r,u,v=(',P_r_element(n),' ->',P_r_ele,' ; ',P_u_element(n),' ->',P_u_ele,' ; ',P_v_element(n),' ->',P_v_ele, &
      ' ) lon,lat,depth,k=( ',x2lon(Xpar_at_setEle(n),Ypar_at_setEle(n)), &
      ' ->',x2lon(Xpar,Ypar),  &
      ' ; ',y2lat(Ypar_at_setEle(n)),&
      ' ->',y2lat(Ypar),&
      ' ; ', Zpar_at_setEle(n),&
      ' ->',Zpar,&
      ' ; ',P_klev_old(n), &
      ' ->',k,' ) and triangle,error=(',triangle,' ; ',error, &
      ' ) Failed to find neighbor Relement, now searching among all elements'
      endif
    ENDIF

    IF(fst.or.(error.ne.0)) then !if the first iteration
      error=0
      !Find rho element in which particle is located
      P_r_ele=0
      triangle=0
      call gridcell(rho_kwele(k),r_kwele_y(1,1,k),                &
            r_kwele_x(1,1,k),Xpar,Ypar,P_r_ele,triangle)
      if(triangle .NE. 0) then
          P_r_element(n) = P_r_ele
      else
          error = 1
      endif
      if(error==1)write(*,*)'setEle failed to find containing Rele searching among all elements'
    ENDIF

    !--------------Find u element in which particle is located------------------
    IF(.not.fst)THEN !if not the first iteration 
      oP_ele = P_u_element(n)
      do i=1,10
        triangle = 0
        if(u_Adjacent(oP_ele,i,k).NE.0)then
           checkele = u_Adjacent(oP_ele,i,k)
           call gridcell(u_kwele(k),u_kwele_y(1,1,k),               &
                       u_kwele_x(1,1,k), &
                       Xpar,Ypar,P_ele,triangle,checkele)
        endif
        if(triangle .NE. 0) then
          P_u_element(n) = P_ele
          exit
        endif
      enddo
      if(triangle.EQ.0)error = 5


      if (error.eq.5 .and. same_vertical_level)then
     write(*,'(a,i12,a,i1,a,i6,4(a,f13.4),2(a,f9.4),6(a,i7),4(a,f10.6),2(a,f9.4),4(a,i3),a)')&
      'setEle warning it=',it,' at call num',num,' for part ',n, &
      ' of pos (X,Y,Z)=(',Xpar_at_setEle(n),&
      ' ->',Xpar,&
      ' ;',Ypar_at_setEle(n),&
      ' ->',Ypar,&
      ' ; ',Zpar_at_setEle(n),&
      ' ->',Zpar, &
      ' ) ; P_r,u,v=(',P_r_element(n),' ->',P_r_ele,' ; ',P_u_element(n),' ->',P_u_ele,' ; ',P_v_element(n),' ->',P_v_ele, &
      ' ) lon,lat,depth,k=( ',x2lon(Xpar_at_setEle(n),Ypar_at_setEle(n)), &
      ' ->',x2lon(Xpar,Ypar),  &
      ' ; ',y2lat(Ypar_at_setEle(n)),&
      ' ->',y2lat(Ypar),&
      ' ; ', Zpar_at_setEle(n),&
      ' ->',Zpar,&
      ' ; ',P_klev_old(n), &
      ' ->',k,' ) and triangle,error=(',triangle,' ; ',error, &
      ' ) Failed to find neighbor Uelement, now searching among all elements'
      endif
    ENDIF

    IF(fst.or.(error.ne.0)) then !if the first iteration
      error=0
      P_u_ele=0
      triangle=0
      call gridcell(u_kwele(k),u_kwele_y(1,1,k),                    &
                    u_kwele_x(1,1,k),                               &
                    Xpar,Ypar,P_u_ele,triangle)
      if(triangle .NE. 0) then
          P_u_element(n) = P_u_ele
      else
          error = 2
      endif
      if(error==2)write(*,*)'setEle failed to find containing Uele searching among all elements'
    ENDIF

    !--------------Find v element in which particle is located------------------
    IF(.not.fst)THEN !if not the first iteration 
      oP_ele = P_v_element(n)
      do i=1,10
        triangle = 0
        if(v_Adjacent(oP_ele,i,k).NE.0) then
           checkele = v_Adjacent(oP_ele,i,k)
           call gridcell(v_kwele(k),v_kwele_y(1,1,k),               &
                       v_kwele_x(1,1,k), &
                       Xpar,Ypar,P_ele,triangle,checkele)
        endif
        if(triangle .NE. 0) then
           P_v_element(n) = P_ele
           exit
        endif
      enddo
      if(triangle.EQ.0)error = 6

      if (error.eq.6 .and. same_vertical_level)then
     write(*,'(a,i12,a,i1,a,i6,4(a,f13.4),2(a,f9.4),6(a,i7),4(a,f10.6),2(a,f9.4),4(a,i3),a)')&
      'setEle warning it=',it,' at call num',num,' for part ',n, &
      ' of pos (X,Y,Z)=(',Xpar_at_setEle(n),&
      ' ->',Xpar,&
      ' ;',Ypar_at_setEle(n),&
      ' ->',Ypar,&
      ' ; ',Zpar_at_setEle(n),&
      ' ->',Zpar, &
      ' ) ; P_r,u,v=(',P_r_element(n),' ->',P_r_ele,' ; ',P_u_element(n),' ->',P_u_ele,' ; ',P_v_element(n),' ->',P_v_ele, &
      ' ) lon,lat,depth,k=( ',x2lon(Xpar_at_setEle(n),Ypar_at_setEle(n)), &
      ' ->',x2lon(Xpar,Ypar),  &
      ' ; ',y2lat(Ypar_at_setEle(n)),&
      ' ->',y2lat(Ypar),&
      ' ; ', Zpar_at_setEle(n),&
      ' ->',Zpar,&
      ' ; ',P_klev_old(n), &
      ' ->',k,' ) and triangle,error=(',triangle,' ; ',error, &
      ' ) Failed to find neighbor Velement, now searching among all elements'
      endif
    ENDIF

    IF(fst.or.(error.ne.0)) then !if the first iteration
      error=0
      P_v_ele=0
      triangle=0
      call gridcell(v_kwele(k),v_kwele_y(1,1,k),                    &
                     v_kwele_x(1,1,k),                              &
                     Xpar,Ypar,P_v_ele,triangle)
      if(triangle .NE. 0) then
          P_v_element(n) = P_v_ele
      else
          error = 3
      endif
      if(error==3)write(*,*)'setEle failed to find containing Vele searching among all elements'
    ENDIF

    ! Done searching in r u and v elements.

    IF(error>0)then
     write(*,*)'-------------------------------'
    !write(*,'(a,i12,a,i1,a,i6,4(a,f13.4),2(a,f9.4),6(a,i7),4(a,f10.6),2(a,f9.4),4(a,i3),a)')&
     write(*,'(a,i12,a,i1,a,i6,a)')&
      'setEle warning it=',it,' at call num',num,' for part ',n, &
    ! ' of pos (X,Y,Z)=(',Xpar_at_setEle(n),&
    ! ' ->',Xpar,&
    ! ' ;',Ypar_at_setEle(n),&
    ! ' ->',Ypar,&
    ! ' ; ',Zpar_at_setEle(n),&
    ! ' ->',Zpar, &
    ! ' ) ; P_r,u,v=(',P_r_element(n),' ->',P_r_ele,' ; ',P_u_element(n),' ->',P_u_ele,' ; ',P_v_element(n),' ->',P_v_ele, &
    ! ' ) lon,lat,depth,k=( ',x2lon(Xpar_at_setEle(n),Ypar_at_setEle(n)), &
    ! ' ->',x2lon(Xpar,Ypar),  &
    ! ' ; ',y2lat(Ypar_at_setEle(n)),&
    ! ' ->',y2lat(Ypar),&
    ! ' ; ', Zpar_at_setEle(n),&
    ! ' ->',Zpar,&
    ! ' ; ',P_klev_old(n), &
    ! ' ->',k,' ) and triangle,error=(',triangle,' ; ',error, &
      ') Searching containing elements among all elements failed'
      ENDIF

    if(num.ne.2 .and. &
    (P_r_element(n)<=0 .or. P_v_element(n)<=0 .or. P_u_element(n)<=0) )then
     if(error.eq.0) error=7
     write(*,*)'----------------------------------------------------------'
     write(*,'(a,i2)')'error while ending setEle call num ',num
     write(*,'(a,i4,3(a,i8))')'n=',n,', found P_r_element=',P_r_element(n),          &
             ', P_u_element=',P_u_element(n),', P_v_element=',P_v_element(n)
     write(*,*)'it may cause further segmentation fault!!!!!!!!!!!!'
     write(*,*)'triangle=',triangle,', checkele=',   &
                   checkele,', first=',first,', error =',error,', call ',num
     !write(*,'(a,i3,a,i12,a,l1,2(a,i2))')'triangle=',triangle,', checkele=',   &
     !              checkele,', first=',first,', error =',error,', call ',num
     write(*,'(a,3f15.9,a,i3)')'test pos is: lon,lat,Z=',x2lon(Xpar,Ypar),y2lat(Ypar),   &
             Zpar,' at level ',k
     write(*,'(2a,3f10.4,a,i3)')'last time that this part was found by setEle ',&
       'to be in an element it was at lon,lat=',&
      x2lon(Xpar_at_setEle(n),Ypar_at_setEle(n)),&
      y2lat(Ypar_at_setEle(n)),Zpar_at_setEle(n),' level ',P_klev_old(n)
     write(*,*)'----------------------------------------------------------'
    endif

    If(error==0)then
    ! Check that nodes are among the range of hydro nodes updated :
      do i=1,4
          updatednode=0
          do j=1,uj
                  if(updatenodesbuffer(1,j,1)<=RE(i,P_r_element(n),k) .and.    &
                     updatenodesbuffer(2,j,1)>=RE(i,P_r_element(n),k) )then
                          updatednode=1
                          exit
                  endif
          enddo
          if(updatednode==0)then
                  write(*,'(a,i6,a,3F10.5,a)')'particle',n,' of pos', &
                       x2lon(Xpar,Ypar),y2lat(Ypar),Zpar, &
                       ' is in a rho element but pb with nodes updates:'
                  write(*,'(a,i1,4(a,i8),a)')'rnode RE(',i,',',P_r_element(n), & 
                          ',',k,')=',RE(i,P_r_element(n),k), &
                          ' composing r_element ',P_r_element(n), &
                          ' is not among hydro-nodes-buffer updated'
                  write(*,'(a)')'You might need to increase ijbuff '
                  if(error==0) error=-3
          endif
          updatednode=0
          do j=1,uj
                  if(updatenodesbuffer(1,j,2)<=UE(i,P_u_element(n),k) .and.    &
                     updatenodesbuffer(2,j,2)>=UE(i,P_u_element(n),k) )then
                          updatednode=1
                          exit
                  endif
          enddo
          if(updatednode==0)then
                  write(*,'(a,i5,a,3f9.5,a)')'particle',n,' of pos',           &
                       x2lon(Xpar,Ypar),y2lat(Ypar),Zpar,                      &
                       ' is in a u element but pb with nodes updates'
                  write(*,'(a,i1,4(a,i8),a)')'rnode UE(',i,',',P_u_element(n), & 
                          ',',k,')=',UE(i,P_u_element(n),k), &
                          ' composing u_element ',P_u_element(n), &
                          ' is not among hydro-nodes-buffer updated'
                  write(*,'(a)')'You might need to increase ijbuff '
                  if(error==0) error=-2
          endif
          updatednode=0
          do j=1,vj
                  if(updatenodesbuffer(1,j,3)<=VE(i,P_v_element(n),k) .and.    &
                     updatenodesbuffer(2,j,3)>=VE(i,P_v_element(n),k) )then
                          updatednode=1
                          exit
                  endif
          enddo
          if(updatednode==0)then
            write(*,'(a,i5,a,3f9.5,a)')'particle',n,' of pos',&
                       x2lon(Xpar,Ypar),y2lat(Ypar),Zpar, & 
                       ' is in a v element but pb with nodes updates:'
            write(*,'(a,i1,2(a,i8),a)')'vnode',i,'=',VE(i,P_v_element(n),k),   &
           ' composing v_element ',P_v_element(n),                             &
           ' is not among hydro-nodes-buffer updated, INCREASE ijbuff'
            if(error==0) error=-3
          endif
      enddo

      !Assign node numbers for rho,u,v calculations
      call setRnode(n)
      Xpar_at_setEle(n)=Xpar
      Ypar_at_setEle(n)=Ypar
      Zpar_at_setEle(n)=Zpar

    EndIf

    if(PRESENT(err))then
            err = error
    endif
    P_klev_old(n)=k

  END SUBROUTINE setEle


  SUBROUTINE setRnode(n)
    !$ use OMP_LIB          
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    INTEGER :: klev
    INTEGER :: rank
     rank=1
     !$ rank=OMP_GET_THREAD_NUM () +1
      klev=P_klev(n)
      !Assign node numbers for rho,u,v calculations

      OMP_ruv(1,RNODE,rank) = RE(1,P_r_element(n),klev)
      OMP_ruv(2,RNODE,rank) = RE(2,P_r_element(n),klev)
      OMP_ruv(3,RNODE,rank) = RE(3,P_r_element(n),klev)
      OMP_ruv(4,RNODE,rank) = RE(4,P_r_element(n),klev)
     
      OMP_ruv(1,UNODE,rank) = UE(1,P_u_element(n),klev)
      OMP_ruv(2,UNODE,rank) = UE(2,P_u_element(n),klev)
      OMP_ruv(3,UNODE,rank) = UE(3,P_u_element(n),klev)
      OMP_ruv(4,UNODE,rank) = UE(4,P_u_element(n),klev)
     
      OMP_ruv(1,VNODE,rank) = VE(1,P_v_element(n),klev)
      OMP_ruv(2,VNODE,rank) = VE(2,P_v_element(n),klev)
      OMP_ruv(3,VNODE,rank) = VE(3,P_v_element(n),klev)
      OMP_ruv(4,VNODE,rank) = VE(4,P_v_element(n),klev)
  END SUBROUTINE setRnode

  SUBROUTINE setDeadOrOut(n)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    DeadOrOut(n)=.True. 
  END SUBROUTINE setDeadOrOut

  SUBROUTINE setEle_all(Xpar,Ypar,Zpar,err,par)
    !This Subroutine determines which Rho, U, and V grid elements contain
    !  each particle
    USE PARAM_MOD, ONLY: numpar,ui,vi,us,ws,Zgrid
    USE GRIDCELL_MOD, ONLY: gridcell
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN), DIMENSION(numpar) :: Xpar,Ypar,Zpar
    INTEGER, INTENT(OUT), OPTIONAL :: err,par ! err = error code
                                              ! par = problem particle
    INTEGER, DIMENSION(3) :: oP_ele_old
    INTEGER :: n,p,i,error,n1
    INTEGER :: triangle,checkele,P_r_ele,P_u_ele,P_v_ele,oP_ele,P_ele,k,kprevious

    error = 0
    p = 0
    n1=0
    do n=1,numpar
      if(DeadOrOut(n)) cycle   
      n1=n
      exit
    enddo

    if(Zgrid)then
      k=getKRlevel(Zpar(n1))
    else
      k=1
    endif
    P_klev(n1)=k

    !Find rho element in which first particle is located 
    P_r_ele=0
    triangle=0
    call gridcell(rho_kwele(k),r_kwele_y(1,1,k),                  &
                  r_kwele_x(1,1,k),                               &
                  Xpar(n1),Ypar(n1),P_r_ele,triangle)
    if (triangle.EQ.0) then
      error = 1
      p = 1
    endif
    P_r_element(n1)=P_r_ele

    !Find u element in which first particle is located
    P_u_ele=0
    triangle=0
    call gridcell(u_kwele(k),u_kwele_y(1,1,k),                      &
                  u_kwele_x(1,1,k),                                 &
                  Xpar(n1),Ypar(n1),P_u_ele,triangle)
    if (triangle.EQ.0) then
      error = 2
      p = 1
    endif
    P_u_element(n1)=P_u_ele

    !Find v element in which first particle is located
    P_v_ele=0
    triangle=0
    call gridcell(v_kwele(k),v_kwele_y(1,1,k),                      &
                       v_kwele_x(1,1,k),                            &
                       Xpar(n1),Ypar(n1),P_v_ele,triangle)
    if (triangle.EQ.0) then
      error = 3
      p = 1
    endif
    P_v_element(n1)=P_v_ele

    P_klev_old(n1)=k
    !Find rho, u, and v elements in which subsequent particles are located
    oP_ele_old(1)=P_r_element(n1)
    oP_ele_old(2)=P_u_element(n1)
    oP_ele_old(3)=P_v_element(n1)

    if(error.ne.0 .or.P_r_element(n1).le.0 .or.P_u_element(n1).le.0 .or.P_v_element(n1).le.0)then
       write(*,*)'part n1 error setEleall',error,P_r_element(n1),P_u_element(n1),P_v_element(n1), &
                          k,Xpar(n1),Ypar(n1)
       stop
    endif
    parloop: do n=n1+1,numpar
      if(DeadOrOut(n)) cycle   
      kprevious=k
      if(Zgrid)then
        k=getKRlevel(Zpar(n))
      else
        k=1
      endif
      P_klev(n)=k

      !Find rho element in which particle is located
      oP_ele = oP_ele_old(1) ! P_r_element(n-1)
      if (maxval(r_Adjacent(oP_ele,:,k)).EQ.0)then
        write(*,*)'ERROR in setEle_all : r_Adjacent is null for element oP_ele=',oP_ele
        stop
      endif
      do i=1,10
        if(r_Adjacent(oP_ele,i,k).EQ.0 .OR. kprevious.NE.k) then
          !If selective search based on previous particle location fails, 
          ! or previous particle had a different k-level
          !  search all
          P_r_ele=0
          triangle=0
          call gridcell(rho_kwele(k),r_kwele_y(1,1,k),            &
                  r_kwele_x(1,1,k),                               &
                  Xpar(n),Ypar(n),P_r_ele,triangle)
          if (triangle.EQ.0) then
            error = 1
            p = n
            write(*,*)'exit parloop 1'
            exit parloop
          endif
          P_r_element(n)=P_r_ele
        else
          triangle = 0
          checkele = r_Adjacent(oP_ele,i,k)
          call gridcell(rho_kwele(k),r_kwele_y(1,1,k),            &
                      r_kwele_x(1,1,k),                           &
                  Xpar(n),Ypar(n),P_ele,triangle,checkele)
          if(triangle .NE. 0) then
            P_r_element(n) = P_ele
            exit
          endif
        endif
      enddo


      !Find u element in which particle is located
      oP_ele = oP_ele_old(2) ! P_u_element(n-1)
      if (maxval(u_Adjacent(oP_ele,:,k)).EQ.0)then
        write(*,*)'ERROR in setEle_all : u_Adjacent is null for element oP_ele=',oP_ele
        stop
      endif
      do i=1,10
        if(u_Adjacent(oP_ele,i,k).EQ.0 .OR. kprevious.NE.k) then
          !If selective search based on previous particle location fails, 
          !  search all
          P_u_ele=0
          triangle=0
          call gridcell(u_kwele(k),u_kwele_y(1,1,k),                &
                  u_kwele_x(1,1,k),                                 &
                  Xpar(n),Ypar(n),P_u_ele,triangle)
          if (triangle.EQ.0) then
            error = 2
            p = n
            write(*,*)'exit parloop 2'
            exit parloop
          endif
          P_u_element(n)=P_u_ele
        else
          triangle = 0
          checkele = u_Adjacent(oP_ele,i,k)
          call gridcell(u_kwele(k),u_kwele_y(1,1,k),                &
                  u_kwele_x(1,1,k),                                 &
                  Xpar(n),Ypar(n),P_ele,triangle,checkele)
          if(triangle .NE. 0) then
            P_u_element(n) = P_ele
            exit
          endif
        endif
      enddo


      !Find v element in which particle is located
      oP_ele = oP_ele_old(3)  ! P_v_element(n-1)
      if (maxval(v_Adjacent(oP_ele,:,k)).EQ.0)then
        write(*,*)'ERROR in setEle_all : v_Adjacent is null for element oP_ele=',oP_ele
        stop
      endif
      do i=1,10
        if(v_Adjacent(oP_ele,i,k).EQ.0 .OR. kprevious.NE.k) then
          !If selective search based on previous particle location fails, 
          !  search all
          P_v_ele=0
          triangle=0
          call gridcell(v_kwele(k),v_kwele_y(1,1,k),                &
                  v_kwele_x(1,1,k),                                 &
                    Xpar(n),Ypar(n),P_v_ele,triangle)
          if (triangle.EQ.0) then
            error = 3
            p = n
            write(*,*)'exit parloop 3'
            exit parloop
          endif
          P_v_element(n)=P_v_ele
        else
          triangle = 0
          checkele = v_Adjacent(oP_ele,i,k)
          call gridcell(v_kwele(k),v_kwele_y(1,1,k),                &
                  v_kwele_x(1,1,k),                                 &
                    Xpar(n),Ypar(n),P_ele,triangle,checkele)
          if(triangle .NE. 0) then
            P_v_element(n) = P_ele
            exit
          endif
        endif
      enddo


     oP_ele_old(1)=P_r_element(n)
     oP_ele_old(2)=P_u_element(n)
     oP_ele_old(3)=P_v_element(n)
     
      P_klev_old(n)=k
    enddo parloop

    if(PRESENT(err)) err = error
    if(PRESENT(par)) par = p

  END SUBROUTINE setEle_all




  SUBROUTINE setInterp(xp,yp,n)
    !This subroutine calculates and stores the interpolation method and values 
    !  for the current particle

    USE PARAM_MOD, ONLY: FreeSlip
    !$ use OMP_LIB          

    IMPLICIT NONE

    INTEGER , INTENT(IN) :: n
    DOUBLE PRECISION, INTENT(IN) :: xp,yp
    double precision x1,x2,x3,x4,y1,y2,y3,y4
    double precision Dis1,Dis2,Dis3,Dis4,TDis
    !NTEGER :: rnode1,rnode2,rnode3,rnode4, &
    !          unode1,unode2,unode3,unode4, &
    !          vnode1,vnode2,vnode3,vnode4
    Integer :: RUV,i
    INTEGER :: tOK
    DOUBLE PRECISION :: t,u,Wgt1,Wgt2,Wgt3,Wgt4,masksum
    DOUBLE PRECISION :: ones(4,4)
    INTEGER :: elemask(4)
    INTEGER :: rank
    ones(:,:)=0.0
    ones(1,1)=1.0
    ones(2,2)=1.0
    ones(3,3)=1.0
    ones(4,4)=1.0
     rank=1
     !$ rank=OMP_GET_THREAD_NUM () +1

    call setRnode(n)
    i=P_klev(n)
    OMP_Xpar(rank) = xp
    OMP_Ypar(rank) = yp
    Do RUV=1,3
      OMP_Pint(:,:,RUV,rank)=ones(:,:)
      Wgt1 = 0
      Wgt2 = 0
      Wgt3 = 0
      Wgt4 = 0
      tOK = 0 !to store information on the interpolation method
      SELECT CASE(RUV)
       CASE(RNODE)
        x1 = rx(OMP_ruv(1,RNODE,rank))
        x2 = rx(OMP_ruv(2,RNODE,rank))
        x3 = rx(OMP_ruv(3,RNODE,rank))
        x4 = rx(OMP_ruv(4,RNODE,rank))
        y1 = ry(OMP_ruv(1,RNODE,rank))
        y2 = ry(OMP_ruv(2,RNODE,rank))
        y3 = ry(OMP_ruv(3,RNODE,rank))
        y4 = ry(OMP_ruv(4,RNODE,rank))
        elemask(1)=rho_mask(OMP_ruv(1,RNODE,rank),i)
        elemask(2)=rho_mask(OMP_ruv(2,RNODE,rank),i)
        elemask(3)=rho_mask(OMP_ruv(3,RNODE,rank),i)
        elemask(4)=rho_mask(OMP_ruv(4,RNODE,rank),i)
       CASE(UNODE)
        x1 = ux(OMP_ruv(1,UNODE,rank))
        x2 = ux(OMP_ruv(2,UNODE,rank))
        x3 = ux(OMP_ruv(3,UNODE,rank))
        x4 = ux(OMP_ruv(4,UNODE,rank))
        y1 = uy(OMP_ruv(1,UNODE,rank))
        y2 = uy(OMP_ruv(2,UNODE,rank))
        y3 = uy(OMP_ruv(3,UNODE,rank))
        y4 = uy(OMP_ruv(4,UNODE,rank))
        elemask(1)=u_mask(OMP_ruv(1,UNODE,rank),i)
        elemask(2)=u_mask(OMP_ruv(2,UNODE,rank),i)
        elemask(3)=u_mask(OMP_ruv(3,UNODE,rank),i)
        elemask(4)=u_mask(OMP_ruv(4,UNODE,rank),i)
       CASE(VNODE)
        x1 = vx(OMP_ruv(1,VNODE,rank))
        x2 = vx(OMP_ruv(2,VNODE,rank))
        x3 = vx(OMP_ruv(3,VNODE,rank))
        x4 = vx(OMP_ruv(4,VNODE,rank))
        y1 = vy(OMP_ruv(1,VNODE,rank))
        y2 = vy(OMP_ruv(2,VNODE,rank))
        y3 = vy(OMP_ruv(3,VNODE,rank))
        y4 = vy(OMP_ruv(4,VNODE,rank))
        elemask(1)=v_mask(OMP_ruv(1,VNODE,rank),i)
        elemask(2)=v_mask(OMP_ruv(2,VNODE,rank),i)
        elemask(3)=v_mask(OMP_ruv(3,VNODE,rank),i)
        elemask(4)=v_mask(OMP_ruv(4,VNODE,rank),i)
      END SELECT     

     !Ensure there is no friction near land (the free slip condition) !ewn.v.2
     !by setting values on land nodes equal to nearby water nodes
      IF (FreeSlip) THEN 
        ! determine if a land element is in present
        masksum = 0
        masksum = elemask(1)+elemask(2)+elemask(3)+elemask(4)
        ! masksum is an integer - the sum of the mask values of each node
        if (masksum .LT. 4) then  
          ! determine how many land nodes are present and reassign values 
          if (masksum .EQ. 3) then   !one land node
             if (elemask(1) .EQ. 0) then
                OMP_Pint(1,1,RUV,rank)=0.0
                OMP_Pint(2,1,RUV,rank)=0.5
                OMP_Pint(4,1,RUV,rank)=0.5   ! v1 = 0.5*(v2+v4)
             else if (elemask(2) .EQ. 0) then 
                OMP_Pint(2,2,RUV,rank)=0.0
                OMP_Pint(1,2,RUV,rank)=0.5
                OMP_Pint(3,2,RUV,rank)=0.5   ! v2 = 0.5*(v1+v3)
             else if (elemask(3) .EQ. 0) then 
                OMP_Pint(3,3,RUV,rank)=0.0
                OMP_Pint(2,3,RUV,rank)=0.5
                OMP_Pint(4,3,RUV,rank)=0.5   ! v3 = 0.5*(v2+v4)
             else if (elemask(4) .EQ. 0) then 
                OMP_Pint(4,4,RUV,rank)=0.0
                OMP_Pint(1,4,RUV,rank)=0.5
                OMP_Pint(3,4,RUV,rank)=0.5   ! v4 = 0.5*(v1+v3)
             end if         
          else if (masksum .EQ. 2) then   !two land nodes
             if (elemask(1).EQ.0 .AND. elemask(2).EQ.0) then
                OMP_Pint(1,1,RUV,rank)=0.0
                OMP_Pint(4,1,RUV,rank)=1.0   ! v1 = v4 
                OMP_Pint(2,2,RUV,rank)=0.0
                OMP_Pint(3,2,RUV,rank)=1.0  !  v2 = v3
             else if (elemask(2).EQ.0.AND.elemask(3).EQ.0) then
                OMP_Pint(2,2,RUV,rank)=0.0
                OMP_Pint(1,2,RUV,rank)=1.0  !  v2 = v1 
                OMP_Pint(3,3,RUV,rank)=0.0
                OMP_Pint(4,3,RUV,rank)=1.0  !  v3 = v4
             else if (elemask(3).EQ.0.AND.elemask(4).EQ.0) then
                OMP_Pint(3,3,RUV,rank)=0.0
                OMP_Pint(2,3,RUV,rank)=1.0  !  v3 = v2 
                OMP_Pint(4,4,RUV,rank)=0.0
                OMP_Pint(1,4,RUV,rank)=1.0  !  v4 = v1
             else if (elemask(4).EQ.0.AND.elemask(1).EQ.0) then
                OMP_Pint(4,4,RUV,rank)=0.0
                OMP_Pint(3,4,RUV,rank)=1.0  !  v4 = v3 
                OMP_Pint(1,1,RUV,rank)=0.0
                OMP_Pint(2,1,RUV,rank)=1.0  !  v1 = v2
             else if (elemask(1).EQ.0.AND.elemask(3).EQ.0) then
                OMP_Pint(1,1,RUV,rank)=0.0
                OMP_Pint(4,1,RUV,rank)=1.0  ! v1 = v4
                OMP_Pint(3,3,RUV,rank)=0.0
                OMP_Pint(2,3,RUV,rank)=1.0  ! v3 = v2
             else if (elemask(4).EQ.0.AND.elemask(2).EQ.0) then
                OMP_Pint(4,4,RUV,rank)=0.0
                OMP_Pint(1,4,RUV,rank)=1.0  ! v4 = v1
                OMP_Pint(2,2,RUV,rank)=0.0
                OMP_Pint(3,2,RUV,rank)=1.0  ! v2 = v3
             end if
          else if (masksum .EQ. 1) then   !three land nodes
             if (elemask(1) .EQ. 1) then
                OMP_Pint(2,2,RUV,rank)=0.0   ! v2 = v1
                OMP_Pint(1,2,RUV,rank)=1.0   ! v2 = v1
                OMP_Pint(3,3,RUV,rank)=0.0   ! v3 = v1  
                OMP_Pint(1,3,RUV,rank)=1.0   ! v3 = v1  
                OMP_Pint(4,4,RUV,rank)=0.0   ! v4 = v1
                OMP_Pint(1,4,RUV,rank)=1.0   ! v4 = v1
             else if (elemask(2) .EQ. 1) then
                OMP_Pint(1,1,RUV,rank)=0.0   ! v1 = v2
                OMP_Pint(2,1,RUV,rank)=1.0   ! v1 = v2
                OMP_Pint(3,3,RUV,rank)=0.0   ! v3 = v2  
                OMP_Pint(2,3,RUV,rank)=1.0   ! v3 = v2  
                OMP_Pint(4,4,RUV,rank)=0.0   ! v4 = v2
                OMP_Pint(2,4,RUV,rank)=1.0   ! v4 = v2
             else if (elemask(3) .EQ. 1) then
                OMP_Pint(1,1,RUV,rank)=0.0   ! v1 = v3
                OMP_Pint(3,1,RUV,rank)=1.0   ! v1 = v3
                OMP_Pint(2,2,RUV,rank)=0.0   ! v2 = v3  
                OMP_Pint(3,2,RUV,rank)=1.0   ! v2 = v3  
                OMP_Pint(4,4,RUV,rank)=0.0   ! v4 = v3
                OMP_Pint(3,4,RUV,rank)=1.0   ! v4 = v3
             else if (elemask(4) .EQ. 1) then
                OMP_Pint(1,1,RUV,rank)=0.0   ! v1 = v4
                OMP_Pint(4,1,RUV,rank)=1.0   ! v1 = v4
                OMP_Pint(2,2,RUV,rank)=0.0   ! v2 = v4  
                OMP_Pint(4,2,RUV,rank)=1.0   ! v2 = v4  
                OMP_Pint(3,3,RUV,rank)=0.0   ! v3 = v4
                OMP_Pint(4,3,RUV,rank)=1.0   ! v3 = v4
             end if
          end if
        end if
      END IF 


 
      ! bilinear interpolation of first triangle
      t = ((xp-x1)*(y3-y1)+(y1-yp)*(x3-x1)) / ((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))
      u = ((xp-x1)*(y2-y1)+(y1-yp)*(x2-x1)) / ((x3-x1)*(y2-y1)-(y3-y1)*(x2-x1))
      tOK = 1 !first triangle
    
      ! if outside triangle, then do bilinear interpolation of other triangle
      if( t.LT.0. .OR. u.LT.0. .OR. (t+u).GT.1.0 ) then
        t = ((xp-x3)*(y1-y3)+(y3-yp)*(x1-x3)) / ((x4-x3)*(y1-y3)-(y4-y3)*(x1-x3))
        u = ((xp-x3)*(y4-y3)+(y3-yp)*(x4-x3)) / ((x1-x3)*(y4-y3)-(y1-y3)*(x4-x3))
        tOK = 2 !second triangle
    
        !if bilinear techniques are undefined, then use inverse weighted distance
        if( t.LT.0. .OR. u.LT.0. .OR. (t+u).GT.1.0 ) then
          !if particle on node, then set equal to node value
          if (     abs(xp-x1)<1e-2 .AND. abs(yp-y1)<1e-2) then
            Wgt1 = 1.0
          elseif ( abs(xp-x2)<1e-2 .AND. abs(yp-y2)<1e-2) then 
            Wgt2 = 1.0
          elseif ( abs(xp-x3)<1e-2 .AND. abs(yp-y3)<1e-2) then 
            Wgt3 = 1.0
          elseif ( abs(xp-x4)<1e-2 .AND. abs(yp-y4)<1e-2) then 
            Wgt4 = 1.0
          else !use inverse weighted distance instead
            !write(*,*)'inverse weighted dist(',x1,'-',xp,')=',x1-xp,' (',y1,'-',yp,')=',y1-yp
            Dis1=1./( SQRT( (x1-xp)**2 + (y1-yp)**2 ) ) 
            Dis2=1./( SQRT( (x2-xp)**2 + (y2-yp)**2 ) ) 
            Dis3=1./( SQRT( (x3-xp)**2 + (y3-yp)**2 ) ) 
            Dis4=1./( SQRT( (x4-xp)**2 + (y4-yp)**2 ) ) 
            TDis = Dis1+Dis2+Dis3+Dis4
            Wgt1= Dis1/TDis
            Wgt2= Dis2/TDis
            Wgt3= Dis3/TDis
            Wgt4= Dis4/TDis
            tOK = 3 !no triangle - used inverse weighted distance
          endif
        endif
      endif     
      OMP_Wgt(1,RUV,rank) = Wgt1
      OMP_Wgt(2,RUV,rank) = Wgt2
      OMP_Wgt(3,RUV,rank) = Wgt3
      OMP_Wgt(4,RUV,rank) = Wgt4
    
      OMP_tOK(RUV,rank) = tOK 
      OMP_t(RUV,rank)   = t
      OMP_u(RUV,rank)   = u
    EndDo

  END SUBROUTINE setInterp


  DOUBLE PRECISION FUNCTION getInterp(xp,yp,var,i)
    !This Function returns the interpolated value at the particle's location 
    !  using the interpolation variables stored from function setInterp, and
    !  the hydrodynamic variables that have been read in

    !$ use OMP_LIB
#include "VAR_IDs.h"
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: xp,yp
    INTEGER, INTENT(IN) :: var
    INTEGER, INTENT(IN), OPTIONAL :: i

    INTEGER :: rmasksum   !ewn.v.2

    DOUBLE PRECISION :: v1,v2,v3,v4,vf1,vf2,vf3,vf4,Pint(4,4)
    INTEGER :: rnode1,rnode2,rnode3,rnode4, &
               unode1,unode2,unode3,unode4, &
               vnode1,vnode2,vnode3,vnode4
    

    INTEGER :: rank,RUV
     rank=1
     !$ rank=OMP_GET_THREAD_NUM () +1

     if( xp.ne.OMP_Xpar(rank) .or. &
         yp.ne.OMP_Ypar(rank) )then
         write(*,*)'ERROR running getInterp of var ',var
         write(*,*)'xp=',xp,' != OMP_Xpar(rank)=',OMP_Xpar(rank),' or yp=', &
         yp,' != OMP_Ypar(rank)=',OMP_Ypar(rank),' ; rank ',rank
         write(*,*)'could be that setInterp was not run previously'
         
         STOP 'END PROGRAM ERROR setInterp not done'
     endif 

 
     rnode1 = OMP_ruv(1,RNODE,rank)  
     rnode2 = OMP_ruv(2,RNODE,rank)  
     rnode3 = OMP_ruv(3,RNODE,rank)  
     rnode4 = OMP_ruv(4,RNODE,rank)  
     unode1 = OMP_ruv(1,UNODE,rank)  
     unode2 = OMP_ruv(2,UNODE,rank)  
     unode3 = OMP_ruv(3,UNODE,rank)  
     unode4 = OMP_ruv(4,UNODE,rank)  
     vnode1 = OMP_ruv(1,VNODE,rank)  
     vnode2 = OMP_ruv(2,VNODE,rank)  
     vnode3 = OMP_ruv(3,VNODE,rank)  
     vnode4 = OMP_ruv(4,VNODE,rank)  
  
     !Determine which data to interpolate from

     RUV = 1

     !determine which data to interpolate from
     SELECT CASE(var)
      !CASE("swdownb")
      !  v1 = t_Swdown(t_b,rnode1)
      !  v2 = t_Swdown(t_b,rnode2)
      !  v3 = t_Swdown(t_b,rnode3)
      !  v4 = t_Swdown(t_b,rnode4)
      !CASE("swdownc")
      !  v1 = t_Swdown(t_c,rnode1)
      !  v2 = t_Swdown(t_c,rnode2)
      !  v3 = t_Swdown(t_c,rnode3)
      !  v4 = t_Swdown(t_c,rnode4)
      !CASE("swdownf")
      !  v1 = t_Swdown(t_f,rnode1)
      !  v2 = t_Swdown(t_f,rnode2)
      !  v3 = t_Swdown(t_f,rnode3)
      !  v4 = t_Swdown(t_f,rnode4)
       CASE(VAR_ID_GrainSize)
         v1 = GrainSize(rnode1)
         v2 = GrainSize(rnode2)
         v3 = GrainSize(rnode3)
         v4 = GrainSize(rnode4)
       CASE(VAR_ID_depth)
         v1 = depth(rnode1)
         v2 = depth(rnode2)
         v3 = depth(rnode3)
         v4 = depth(rnode4)
       CASE(VAR_ID_angle)
         v1 = rho_angle(rnode1)
         v2 = rho_angle(rnode2)
         v3 = rho_angle(rnode3)
         v4 = rho_angle(rnode4)
       CASE(VAR_ID_zetab)
         v1 = t_zeta(t_b,rnode1)
         v2 = t_zeta(t_b,rnode2)
         v3 = t_zeta(t_b,rnode3)
         v4 = t_zeta(t_b,rnode4)
       CASE(VAR_ID_zetac)
         v1 = t_zeta(t_c,rnode1)
         v2 = t_zeta(t_c,rnode2)
         v3 = t_zeta(t_c,rnode3)
         v4 = t_zeta(t_c,rnode4)
       CASE(VAR_ID_zetaf)
         v1 = t_zeta(t_f,rnode1)
         v2 = t_zeta(t_f,rnode2)
         v3 = t_zeta(t_f,rnode3)
         v4 = t_zeta(t_f,rnode4)
       CASE(VAR_ID_saltb)
         if(.not. present(i))then
           write(*,*) 'Problem interpolating ',var
           write(*,*) 'Optional Argument (i) Required for this Variable'
           write(*,*) ' '
           write(*,*) 'The Program Cannot Continue and Will Terminate'
           stop
         endif
         v1 = t_salt(t_b,rnode1,i)
         v2 = t_salt(t_b,rnode2,i)
         v3 = t_salt(t_b,rnode3,i)
         v4 = t_salt(t_b,rnode4,i)
       CASE(VAR_ID_saltc)
         if(.not. present(i))then
           write(*,*) 'Problem interpolating ',var
           write(*,*) 'Optional Argument (i) Required for this Variable'
           write(*,*) ' '
           write(*,*) 'The Program Cannot Continue and Will Terminate'
           stop
         endif
         v1 = t_salt(t_c,rnode1,i)
         v2 = t_salt(t_c,rnode2,i)
         v3 = t_salt(t_c,rnode3,i)
         v4 = t_salt(t_c,rnode4,i)
       CASE(VAR_ID_saltf)
         if(.not. present(i))then
           write(*,*) 'Problem interpolating ',var
           write(*,*) 'Optional Argument (i) Required for this Variable'
           write(*,*) ' '
           write(*,*) 'The Program Cannot Continue and Will Terminate'
           stop
         endif
         v1 = t_salt(t_f,rnode1,i)
         v2 = t_salt(t_f,rnode2,i)
         v3 = t_salt(t_f,rnode3,i)
         v4 = t_salt(t_f,rnode4,i)
       CASE(VAR_ID_tempb)
         if(.not. present(i))then
           write(*,*) 'Problem interpolating ',var
           write(*,*) 'Optional Argument (i) Required for this Variable'
           write(*,*) ' '
           write(*,*) 'The Program Cannot Continue and Will Terminate'
           stop
         endif
         v1 = t_temp(t_b,rnode1,i)
         v2 = t_temp(t_b,rnode2,i)
         v3 = t_temp(t_b,rnode3,i)
         v4 = t_temp(t_b,rnode4,i)
       CASE(VAR_ID_tempc)
         if(.not. present(i))then
           write(*,*) 'Problem interpolating ',var
           write(*,*) 'Optional Argument (i) Required for this Variable'
           write(*,*) ' '
           write(*,*) 'The Program Cannot Continue and Will Terminate'
           stop
         endif
         v1 = t_temp(t_c,rnode1,i)
         v2 = t_temp(t_c,rnode2,i)
         v3 = t_temp(t_c,rnode3,i)
         v4 = t_temp(t_c,rnode4,i)
       CASE(VAR_ID_tempf)
         if(.not. present(i))then
           write(*,*) 'Problem interpolating ',var
           write(*,*) 'Optional Argument (i) Required for this Variable'
           write(*,*) ' '
           write(*,*) 'The Program Cannot Continue and Will Terminate'
           stop
         endif
         v1 = t_temp(t_f,rnode1,i)
         v2 = t_temp(t_f,rnode2,i)
         v3 = t_temp(t_f,rnode3,i)
         v4 = t_temp(t_f,rnode4,i)
       CASE(VAR_ID_denb)
         if(.not. present(i))then
           write(*,*) 'Problem interpolating ',var
           write(*,*) 'Optional Argument (i) Required for this Variable'
           write(*,*) ' '
           write(*,*) 'The Program Cannot Continue and Will Terminate'
           stop
         endif
         v1 = t_den(t_b,rnode1,i)
         v2 = t_den(t_b,rnode2,i)
         v3 = t_den(t_b,rnode3,i)
         v4 = t_den(t_b,rnode4,i)
       CASE(VAR_ID_denc)
         if(.not. present(i))then
           write(*,*) 'Problem interpolating ',var
           write(*,*) 'Optional Argument (i) Required for this Variable'
           write(*,*) ' '
           write(*,*) 'The Program Cannot Continue and Will Terminate'
           stop
         endif
         v1 = t_den(t_c,rnode1,i)
         v2 = t_den(t_c,rnode2,i)
         v3 = t_den(t_c,rnode3,i)
         v4 = t_den(t_c,rnode4,i)
       CASE(VAR_ID_denf)
         if(.not. present(i))then
           write(*,*) 'Problem interpolating ',var
           write(*,*) 'Optional Argument (i) Required for this Variable'
           write(*,*) ' '
           write(*,*) 'The Program Cannot Continue and Will Terminate'
           stop
         endif
         v1 = t_den(t_f,rnode1,i)
         v2 = t_den(t_f,rnode2,i)
         v3 = t_den(t_f,rnode3,i)
         v4 = t_den(t_f,rnode4,i)
       CASE(VAR_ID_uvelb)
         if(.not. present(i))then
           write(*,*) 'Problem interpolating ',var
           write(*,*) 'Optional Argument (i) Required for this Variable'
           write(*,*) ' '
           write(*,*) 'The Program Cannot Continue and Will Terminate'
           stop
         endif
         v1 = t_Uvel(t_b,unode1,i)
         v2 = t_Uvel(t_b,unode2,i)
         v3 = t_Uvel(t_b,unode3,i)
         v4 = t_Uvel(t_b,unode4,i)
         RUV = 2
       CASE(VAR_ID_uvelc)
         if(.not. present(i))then
           write(*,*) 'Problem interpolating ',var
           write(*,*) 'Optional Argument (i) Required for this Variable'
           write(*,*) ' '
           write(*,*) 'The Program Cannot Continue and Will Terminate'
           stop
         endif
         v1 = t_Uvel(t_c,unode1,i)
         v2 = t_Uvel(t_c,unode2,i)
         v3 = t_Uvel(t_c,unode3,i)
         v4 = t_Uvel(t_c,unode4,i)
         RUV = 2
       CASE(VAR_ID_uvelf)
         if(.not. present(i))then
           write(*,*) 'Problem interpolating ',var
           write(*,*) 'Optional Argument (i) Required for this Variable'
           write(*,*) ' '
           write(*,*) 'The Program Cannot Continue and Will Terminate'
           stop
         endif
         v1 = t_Uvel(t_f,unode1,i)
         v2 = t_Uvel(t_f,unode2,i)
         v3 = t_Uvel(t_f,unode3,i)
         v4 = t_Uvel(t_f,unode4,i)
         RUV = 2
       CASE(VAR_ID_vvelb)
         if(.not. present(i))then
           write(*,*) 'Problem interpolating ',var
           write(*,*) 'Optional Argument (i) Required for this Variable'
           write(*,*) ' '
           write(*,*) 'The Program Cannot Continue and Will Terminate'
           stop
         endif
         v1 = t_Vvel(t_b,vnode1,i)
         v2 = t_Vvel(t_b,vnode2,i)
         v3 = t_Vvel(t_b,vnode3,i)
         v4 = t_Vvel(t_b,vnode4,i)
         RUV = 3
       CASE(VAR_ID_vvelc)
         if(.not. present(i))then
           write(*,*) 'Problem interpolating ',var
           write(*,*) 'Optional Argument (i) Required for this Variable'
           write(*,*) ' '
           write(*,*) 'The Program Cannot Continue and Will Terminate'
           stop
         endif
         v1 = t_Vvel(t_c,vnode1,i)
         v2 = t_Vvel(t_c,vnode2,i)
         v3 = t_Vvel(t_c,vnode3,i)
         v4 = t_Vvel(t_c,vnode4,i)
         RUV = 3
       CASE(VAR_ID_vvelf)
         if(.not. present(i))then
           write(*,*) 'Problem interpolating ',var
           write(*,*) 'Optional Argument (i) Required for this Variable'
           write(*,*) ' '
           write(*,*) 'The Program Cannot Continue and Will Terminate'
           stop
         endif
         v1 = t_Vvel(t_f,vnode1,i)
         v2 = t_Vvel(t_f,vnode2,i)
         v3 = t_Vvel(t_f,vnode3,i)
         v4 = t_Vvel(t_f,vnode4,i)
         RUV = 3
       CASE(VAR_ID_wvelb)
         if(.not. present(i))then
           write(*,*) 'Problem interpolating ',var
           write(*,*) 'Optional Argument (i) Required for this Variable'
           write(*,*) ' '
           write(*,*) 'The Program Cannot Continue and Will Terminate'
           stop
         endif
         v1 = t_Wvel(t_b,rnode1,i)
         v2 = t_Wvel(t_b,rnode2,i)
         v3 = t_Wvel(t_b,rnode3,i)
         v4 = t_Wvel(t_b,rnode4,i)
       CASE(VAR_ID_wvelc)
         if(.not. present(i))then
           write(*,*) 'Problem interpolating ',var
           write(*,*) 'Optional Argument (i) Required for this Variable'
           write(*,*) ' '
           write(*,*) 'The Program Cannot Continue and Will Terminate'
           stop
         endif
         v1 = t_Wvel(t_c,rnode1,i)
         v2 = t_Wvel(t_c,rnode2,i)
         v3 = t_Wvel(t_c,rnode3,i)
         v4 = t_Wvel(t_c,rnode4,i)
       CASE(VAR_ID_wvelf)
         if(.not. present(i))then
           write(*,*) 'Problem interpolating ',var
           write(*,*) 'Optional Argument (i) Required for this Variable'
           write(*,*) ' '
           write(*,*) 'The Program Cannot Continue and Will Terminate'
           stop
         endif
         v1 = t_Wvel(t_f,rnode1,i)
         v2 = t_Wvel(t_f,rnode2,i)
         v3 = t_Wvel(t_f,rnode3,i)
         v4 = t_Wvel(t_f,rnode4,i)
       CASE(VAR_ID_khb)
         if(.not. present(i))then
           write(*,*) 'Problem interpolating ',var
           write(*,*) 'Optional Argument (i) Required for this Variable'
           write(*,*) ' '
           write(*,*) 'The Program Cannot Continue and Will Terminate'
           stop
         endif
         v1 = t_KH(t_b,rnode1,i)
         v2 = t_KH(t_b,rnode2,i)
         v3 = t_KH(t_b,rnode3,i)
         v4 = t_KH(t_b,rnode4,i)
       CASE(VAR_ID_khc)
         if(.not. present(i))then
           write(*,*) 'Problem interpolating ',var
           write(*,*) 'Optional Argument (i) Required for this Variable'
           write(*,*) ' '
           write(*,*) 'The Program Cannot Continue and Will Terminate'
           stop
         endif
         v1 = t_KH(t_c,rnode1,i)
         v2 = t_KH(t_c,rnode2,i)
         v3 = t_KH(t_c,rnode3,i)
         v4 = t_KH(t_c,rnode4,i)
       CASE(VAR_ID_khf)
         if(.not. present(i))then
           write(*,*) 'Problem interpolating ',var
           write(*,*) 'Optional Argument (i) Required for this Variable'
           write(*,*) ' '
           write(*,*) 'The Program Cannot Continue and Will Terminate'
           stop
         endif
         v1 = t_KH(t_f,rnode1,i)
         v2 = t_KH(t_f,rnode2,i)
         v3 = t_KH(t_f,rnode3,i)
         v4 = t_KH(t_f,rnode4,i)
       CASE(VAR_ID_uwindb)
         v1 = t_uwind(t_b,unode1)
         v2 = t_uwind(t_b,unode2)
         v3 = t_uwind(t_b,unode3)
         v4 = t_uwind(t_b,unode4)
       CASE(VAR_ID_uwindc)
         v1 = t_uwind(t_c,unode1)
         v2 = t_uwind(t_c,unode2)
         v3 = t_uwind(t_c,unode3)
         v4 = t_uwind(t_c,unode4)
       CASE(VAR_ID_uwindf)
         v1 = t_uwind(t_f,unode1)
         v2 = t_uwind(t_f,unode2)
         v3 = t_uwind(t_f,unode3)
         v4 = t_uwind(t_f,unode4)
       CASE(VAR_ID_vwindb)
         v1 = t_vwind(t_b,vnode1)
         v2 = t_vwind(t_b,vnode2)
         v3 = t_vwind(t_b,vnode3)
         v4 = t_vwind(t_b,vnode4)
       CASE(VAR_ID_vwindc)
         v1 = t_vwind(t_c,vnode1)
         v2 = t_vwind(t_c,vnode2)
         v3 = t_vwind(t_c,vnode3)
         v4 = t_vwind(t_c,vnode4)
       CASE(VAR_ID_vwindf)
         v1 = t_vwind(t_f,vnode1)
         v2 = t_vwind(t_f,vnode2)
         v3 = t_vwind(t_f,vnode3)
         v4 = t_vwind(t_f,vnode4)
       CASE(VAR_ID_iwindb)
         v1 = t_iwind(t_b,rnode1)
         v2 = t_iwind(t_b,rnode2)
         v3 = t_iwind(t_b,rnode3)
         v4 = t_iwind(t_b,rnode4)
       CASE(VAR_ID_iwindc)
         v1 = t_iwind(t_c,rnode1)
         v2 = t_iwind(t_c,rnode2)
         v3 = t_iwind(t_c,rnode3)
         v4 = t_iwind(t_c,rnode4)
       CASE(VAR_ID_iwindf)
         v1 = t_iwind(t_f,rnode1)
         v2 = t_iwind(t_f,rnode2)
         v3 = t_iwind(t_f,rnode3)
         v4 = t_iwind(t_f,rnode4)
       CASE DEFAULT
         write(*,*) 'Problem interpolating ',var
         write(*,*) ' '
         write(*,*) 'The Program Cannot Continue and Will Terminate'
         stop
     END SELECT

     Pint(:,:)=OMP_Pint(:,:,RUV,rank)
     vf1 = Pint(1,1)*v1 + Pint(2,1)*v2 + Pint(3,1)*v3 + Pint(4,1)*v4
     vf2 = Pint(1,2)*v1 + Pint(2,2)*v2 + Pint(3,2)*v3 + Pint(4,2)*v4
     vf3 = Pint(1,3)*v1 + Pint(2,3)*v2 + Pint(3,3)*v3 + Pint(4,3)*v4
     vf4 = Pint(1,4)*v1 + Pint(2,4)*v2 + Pint(3,4)*v3 + Pint(4,4)*v4

     !interpolate using the variables from setInterp
     if(OMP_tOK(RUV,rank) == 1) then 
       getInterp = vf1 + (vf2-vf1)*OMP_t(RUV,rank) + (vf3-vf1)*OMP_u(RUV,rank)
     elseif(OMP_tOK(RUV,rank) == 2) then
       getInterp = vf3 + (vf4-vf3)*OMP_t(RUV,rank)+ (vf1-vf3)*OMP_u(RUV,rank)
     else 
       getInterp = OMP_Wgt(1,RUV,rank)*vf1 + OMP_Wgt(2,RUV,rank)*vf2 + &
                   OMP_Wgt(3,RUV,rank)*vf3 + OMP_Wgt(4,RUV,rank)*vf4 
     endif

  END FUNCTION getInterp

  !This function creates a Water Column Tension Spline at back, center, and 
  !  forward hydrodynamic time then uses Polynomial Interpolation to determine
  !  Internal Time values to finally get the value of the particle in space and
  !  time.  The name is derived from: Water Column Tension Spline, Internal 
  !  Time Polynomial Interpolation.  The final variable v is for version
  !  (ie what is to be returned): 1-back, 2-center, 3-forward, 4-(b+4c+f)/6
  DOUBLE PRECISION FUNCTION WCTS_ITPI(var,Xpos,Ypos,deplvl,Pwc_zb,Pwc_zc,      &
                                Pwc_zf,slvls,P_zb,P_zc,P_zf,ex,ix,p,v,nP,nN)
    USE TENSION_MOD, ONLY: TSPSI,HVAL
    USE INT_MOD, ONLY: linint,polintd
    USE CONVERT_MOD, ONLY: x2lon,y2lat
    USE PARAM_MOD, ONLY: LinearVInterp 
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: deplvl,slvls,p,v,nP
    DOUBLE PRECISION, INTENT(IN) :: Xpos,Ypos,Pwc_zb(slvls),Pwc_zc(slvls),     &
                                    Pwc_zf(slvls),P_zb,P_zc,P_zf,ex(3),ix(3)
    INTEGER, INTENT(IN) :: var

    INTEGER,INTENT(INOUT), OPTIONAL :: nN 

    INTEGER :: varb,varc,varf
    INTEGER ::  i
    DOUBLE PRECISION ,ALLOCATABLE, DIMENSION(:) :: abb_zb,abb_zc,abb_zf,abb_vb,&
      abb_vc,abb_vf
    DOUBLE PRECISION P_vb,P_vc,P_vf,P_V,ey(3),vb,vc,vf,slope

    !TSPACK Variables
    INTEGER :: IER,SigErr
    DOUBLE PRECISION ,ALLOCATABLE, DIMENSION(:) :: YP,SIGM
    if(.not. present(nN)) nN=4

    !BEUG: not optimal to allocate and disalocate that way, would be better
    !to allocate just once at the beginning of the program and then use just
    !partially the array...
    ALLOCATE(abb_zb(nN))
    ALLOCATE(abb_zc(nN))
    ALLOCATE(abb_zf(nN))
    ALLOCATE(abb_vb(nN))
    ALLOCATE(abb_vc(nN))
    ALLOCATE(abb_vf(nN))
    ALLOCATE(SIGM(nN))
    ALLOCATE(YP(nN))
    varb = var+1
    varc = var+2
    varf = var+3
    do i=1,nN
      abb_zb(i) = Pwc_zb(i+deplvl-1)
      abb_zc(i) = Pwc_zc(i+deplvl-1)
      abb_zf(i) = Pwc_zf(i+deplvl-1)
      abb_vb(i) = getInterp(Xpos,Ypos,varb,i+deplvl-1)
      abb_vc(i) = getInterp(Xpos,Ypos,varc,i+deplvl-1)
      abb_vf(i) = getInterp(Xpos,Ypos,varf,i+deplvl-1)
    enddo
!       *********************************************************
!       *       5Aiic6b.  Fit Tension Spline to WC Profile      *
!       *********************************************************

    !ii. call TSPACK to fit a tension spline to water column profile 
    !  of U,V,W velocities at particle x-y location and find value at particle
    P_vb=0.0
    SigErr=0
    IF(nN>=3)then
      CALL TSPSI (nN,abb_zb,abb_vb,YP,SIGM,IER,SigErr)
      IF (SigErr.EQ.0) THEN
        P_vb = HVAL (P_zb,nN,abb_zb,abb_vb,YP,SIGM,IER)
      ENDIF
    ENDIF
    IF(((nN.eq.2..or.SigErr.ne.0).or.LinearVInterp))then
      CALL linint(abb_zb,abb_vb,nN,P_zb,P_vb,slope)
    ELSEIF(nN.lt.2)then
      P_vb=abb_vb(nN)
    ENDIF

    P_vc=0.0  
    SigErr=0
    IF(nN>=3 )then
        CALL TSPSI (nN,abb_zc,abb_vc,YP,SIGM,IER,SigErr)
        IF (SigErr.EQ.0) THEN
          P_vc = HVAL (P_zc,nN,abb_zc,abb_vc,YP,SIGM,IER)
        ENDIF
    ENDIF
    IF(((nN.eq.2..or.SigErr.ne.0).or.LinearVInterp))then
      CALL linint(abb_zc,abb_vc,nN,P_zc,P_vc,slope)
    ELSEIF(nN.lt.2)then
      P_vc=abb_vc(nN)
    ENDIF

    P_vf=0.0
    SigErr=0
    IF(nN>=3)then
        CALL TSPSI (nN,abb_zf,abb_vf,YP,SIGM,IER,SigErr)
        IF (SigErr.EQ.0) THEN
          P_vf = HVAL (P_zf,nN,abb_zf,abb_vf,YP,SIGM,IER)
        ENDIF

    ENDIF
    IF(((nN.eq.2..or.SigErr.ne.0).or.LinearVInterp))then
      CALL linint(abb_zf,abb_vf,nN,P_zf,P_vf,slope)
    ELSEIF(nN.lt.2)then
      P_vf=abb_vf(nN)
    ENDIF
!       *********************************************************
!       *               Find Internal b,c,f Values              *
!       *********************************************************

    !iii. fit polynomial to hydrodynamic model output and find 
    !  internal b,c,f values

    !    1. Prepare external time step values
    if (p .EQ. 1) then
      ey=0.0
      ey(1) = P_vb
      ey(2) = P_vb
      ey(3) = P_vc
    else
      ey=0.0
      ey(1) = P_vb
      ey(2) = P_vc
      ey(3) = P_vf
    endif

    !    2. Get value
    vb = polintd(ex,ey,3,ix(1))
    vc = polintd(ex,ey,3,ix(2))
    vf = polintd(ex,ey,3,ix(3))
    P_V = (vb + vc*4 + vf) / DBLE(6.0)


    SELECT CASE (v)
      CASE (1)
        WCTS_ITPI = vb
      CASE (2)
        WCTS_ITPI = vc
      CASE (3)
        WCTS_ITPI = vf
      CASE (4)
        WCTS_ITPI = P_V
      CASE DEFAULT
        write(*,*) 'ERROR: Illegal WCTS_ITPI version number'
        write(*,*) ' '
        write(*,*) 'The Program Cannot Continue and Will Terminate'
        stop
    END SELECT

  DEALLOCATE(abb_zb,abb_zc,abb_zf,abb_vb,abb_vc,abb_vf,SIGM,YP)
  END FUNCTION WCTS_ITPI

  SUBROUTINE getDepth(Xpos,Ypos,n,idt,Pdepth,kbottom,conflict)
    USE PARAM_MOD, ONLY: us,ws,ui,vi,uj,vj,Zgrid_depthinterp
    USE CONVERT_MOD, ONLY: x2lon,y2lat
    USE PIP_MOD, ONLY: inpoly
    !$ use OMP_LIB
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: Xpos,Ypos
    DOUBLE PRECISION, INTENT(OUT) :: Pdepth
    INTEGER, INTENT(IN) :: n,idt
    INTEGER, INTENT(INOUT) :: conflict
    INTEGER, INTENT(OUT) :: kbottom
    INTEGER :: kbA,kbB,kbC,conflictin
    INTEGER :: vnodeA,unodeB,neighborRnodeA,neighborRnodeB,neighborRnodeC
    INTEGER :: delta_i,delta_j,modrnode1,modrnode2,modrnode3,modrnode4
    INTEGER :: i,j,k,i_r,i_u,i_v,j_r,j_u,j_v,kbot,Unum,Vnum
    DOUBLE PRECISION :: alpha,beta,aA,bA,aB,bB,aC,bC,YpA,YpB,YpC
    DOUBLE PRECISION :: minx,miny,maxx,maxy,dx,dy
    DOUBLE PRECISION :: x1,x2,x3,y1,y2,y3,v1,v2,v3,tt,uu,vp,&
                        Dis1,Dis2,Dis3,TDis,Wt1,Wt2,Wt3
    DOUBLE PRECISION,DIMENSION(4,2)::trilatlon
    CHARACTER(LEN=200) :: string,title,marker
    LOGICAL :: tri_ext
    INTEGER :: closernode,DepthRefRnode !--- CL-OGS 
    INTEGER :: rnode1,rnode2,rnode3,rnode4, &
               unode1,unode2,unode3,unode4, &
               vnode1,vnode2,vnode3,vnode4

    INTEGER :: rank
     rank=1
     !$ rank=OMP_GET_THREAD_NUM () +1
    i=-1
    j=-1
    k=-1 
    rnode1 = OMP_ruv(1,RNODE,rank)  
    rnode2 = OMP_ruv(2,RNODE,rank)  
    rnode3 = OMP_ruv(3,RNODE,rank)  
    rnode4 = OMP_ruv(4,RNODE,rank)  
    unode1 = OMP_ruv(1,UNODE,rank)  
    unode2 = OMP_ruv(2,UNODE,rank)  
    unode3 = OMP_ruv(3,UNODE,rank)  
    unode4 = OMP_ruv(4,UNODE,rank)  
    vnode1 = OMP_ruv(1,VNODE,rank)  
    vnode2 = OMP_ruv(2,VNODE,rank)  
    vnode3 = OMP_ruv(3,VNODE,rank)  
    vnode4 = OMP_ruv(4,VNODE,rank)  
   conflictin=conflict
   if(conflictin<0)write(*,*)' HCONFLICT GetDepth at ',x2lon(Xpos,Ypos),y2lat(Ypos) 
   DepthRefRnode=-1
    conflict=1
! in the following, i and j are indexes of the elements made by 4 nodes:
!  do j=1,uj-1 ,  do i=1,vi-1 : rnode1 = i + (j-1)*vi -> i_r=mod(rnode1,vi) ; j_r = (rnode1-i_r)/vi+1
!  do j=1,uj-1 ,  do i=1,ui-1 : unode1 = i + (j-1)*ui -> i_u=mod(unode1,ui) ; j_u = (unode1-i_u)/ui+1
!  do j=1,vj-1,   do i=1,vi-1 : vnode1 = i + (j-1)*vi -> i_v=mod(vnode1,vi) ; j_v = (vnode1-i_v)/vi+1
   i_r=mod(rnode1,vi) 
   i_u=mod(unode1,ui)
   i_v=mod(vnode1,vi)
   j_r = (rnode1-i_r)/vi+1
   j_u = (unode1-i_u)/ui+1
   j_v = (vnode1-i_v)/vi+1
! make sure that rnode1 and vnode1 corrispond to the same index i :
    delta_i=mod(rnode1,vi)-mod(vnode1,vi)             ! i_r-i_v
! make sure that rnode1 and unode1 corrispond to the same index j :
    !delta_j=(((rnode1-mod(rnode1-1,vi))/vi +1)      &  
    !        -((unode1-mod(unode1-1,ui))/ui +1) )*ui        
    delta_j= j_r-j_u
    if((abs(delta_i)>1) .or. (abs(delta_j)>1))then
        write(*,*)'error part ',n,' ruv nodes deltaij'
        write(*,*)'di=',delta_i,' dj=',delta_j,' Ni_rv=',vi,' Ni_u=',ui, &
           'ir,iu,iv =(',i_r,i_u,i_v,'), jr,ju,jv=(',j_r,j_u,j_v,')'
        write(*,*)'error r,mr,u,v nodes1',rnode1,unode1,vnode1
        write(*,*)'error r,mr,u,v nodes2',rnode2,unode2,vnode2
        write(*,*)'error r,mr,u,v nodes3',rnode3,unode3,vnode3
        write(*,*)'error r,mr,u,v nodes4',rnode4,unode4,vnode4
        write(*,*)delta_i,delta_j,vi,ui,'(',i_r,i_u,i_v,'), (',j_r,j_u,j_v,')'
        write(*,*)'Xpos,Ypos,n,idt,Pdepth=',Xpos,Ypos,n,idt,Pdepth
        write(*,*)'x2lon(Xpos,Ypos),y2lat(Ypos)= ',x2lon(Xpos,Ypos),y2lat(Ypos)
        conflict=-2
        !stop     
    endif 
   
    modrnode1=rnode1 -delta_i -delta_j*vi
    modrnode2=rnode2 -delta_i -delta_j*vi
    modrnode3=rnode3 -delta_i -delta_j*vi
    modrnode4=rnode4 -delta_i -delta_j*vi
    
   i_r=mod(modrnode1,vi) 
   i_u=mod(unode1,ui)
   i_v=mod(vnode1,vi)
   j_r = (modrnode1-i_r)/vi+1
   j_u = (unode1-i_u)/ui+1
   j_v = (vnode1-i_v)/vi+1
   if( ((i_r.ne.i_v) .or. (j_r.ne.j_u)) .or. &
        (i_r.ne.i_u .and. i_r.ne.i_u+1) .or. &
        (j_r.ne.j_v .and. j_r.ne.j_v+1) .or. &
        (modrnode1.eq.0.or.unode1.eq.0.or.vnode1.eq.0) )then   
        write(*,*)'error r,mr,u,v nodes1',rnode1,modrnode1,unode1,vnode1
        write(*,*)'error r,mr,u,v nodes2',rnode2,modrnode2,unode2,vnode2
        write(*,*)'error r,mr,u,v nodes3',rnode3,modrnode3,unode3,vnode3
        write(*,*)'error r,mr,u,v nodes4',rnode4,modrnode4,unode4,vnode4
        write(*,*)delta_i,delta_j,vi,ui,'(',i_r,i_u,i_v,'), (',j_r,j_u,j_v,')'
        write(*,*)'Xpos,Ypos,n,idt,Pdepth=',Xpos,Ypos,n,idt,Pdepth
        write(*,*)'x2lon(Xpos,Ypos),y2lat(Ypos)= ',x2lon(Xpos,Ypos),y2lat(Ypos)
        conflict=-2
        !stop     
    endif 
! 00 ! 01 ! 02 ! 03 ! 04 ! 05 ! 06 ! 07 ! 08 ! 09 ! 10 ! 11 ! 12 ! 13 ! 14 ! 15 !  16  !  17  !  18  !  19        !
! .. ! x. ! .. ! x. ! .x ! xx !    ! xx ! .. !    ! .. ! x. ! .x ! xx ! .x ! xx ! \.\x ! ./x/ ! \x\. !  x/./      !
! .. ! .. ! x. ! x. ! .. ! .. !    ! x. ! .x !    ! xx ! xx ! .x ! .x ! xx ! xx ! x\.\ ! /x/. ! .\x\ !  /./x      !
    kbottom=-1
    tri_ext=.FALSE. 
     !--------------------------------------------------------------------------------------------------------- 
    IF(     (rx(modrnode1)>ux(unode1))        &                     !      . . . . . . . . . . . . . .
       .and.(ry(modrnode1)> vy(vnode1)))THEN                        !      .            .            .
       string='R1, V4, U2'                                          !      .            .            .
       closernode=modrnode1                                         !     U4     R4____U3_____R3     .    
       vnodeA=vnode4                                                !      .     |      .     |      .    
       unodeB=unode2                                                !      .     |      .     |      .    
       ! alpha=(yB-yA)/(xB-xA)                                      !      . . . V4 . . . . . V3 . . .    
       ! beta =yA-alpha*Xa                                          !      .     |  \ Ri.     |      .    
       alpha=( uy(unodeB)-vy(vnodeA) )/( ux(unodeB)-vx(vnodeA) )    !      .     |Re \  .     |      .    
       beta =vy(vnodeA)- alpha*vx(vnodeA)                           !     U1     R1___\U2_____R2     .    
       if(Ypos.lt.alpha*Xpos+beta)  tri_ext=.TRUE.                  !      .            .            . 
       Unum= unode2                                                 !      .            .            . 
       Vnum= vnode4                                                 !      . . . V1 . . . . . V2 . . .
       !write(*,'(a,i6,a)')'partonborder[',n-1,']=1   #PYTHON DEBUG'
     !-----------------------------------------------------------------------------------------------------------     
    ELSEIF( (rx(modrnode1)<ux(unode1))         &                    !      . . . . . . . . . . . . . . 
       .and.(ry(modrnode1)> vy(vnode1)))THEN                        !      .            .            . 
       string='R2, V3, U1'                                          !      .            .            . 
       closernode=modrnode2                                         !      .     R4____U4_____R3     U3 
       vnodeA=vnode3                                                !      .     |      .     |      .  
       unodeB=unode1                                                !      .     |      .     |      .  
       ! alpha=(yB-yA)/(xB-xA)                                      !      . . . V4 . . . . . V3 . . .  
       ! beta =yA-alpha*Xa                                          !      .     |      .Ri / |      .  
       alpha=( uy(unodeB)-vy(vnodeA) )/( ux(unodeB)-vx(vnodeA) )    !      .     |      .  /Re|      .  
       beta =vy(vnodeA)- alpha*vx(vnodeA)                           !      .     R1____U1_/___R2     U2
       if(Ypos.lt.alpha*Xpos+beta)    tri_ext=.TRUE.                !      .            .            . 
       Unum= unode1                                                 !      .            .            . 
       Vnum= vnode4                                                 !      . . . V1 . . . . . V2 . . . 
       !write(*,'(a,i6,a)')'partonborder[',n-1,']=3   #PYTHON DEBUG'
     !--------------------------------------------------------------------------------------------------------- 
     ELSEIF( (rx(modrnode1)<ux(unode1))        &                    !      . . . V4 . . . . . V3 . . . 
       .and.(ry(modrnode1)< vy(vnode1)))THEN                        !      .            .            . 
         string='R3, V2, U4'                                        !      .            .            . 
         closernode=modrnode3                                       !      .     R4____U4_____R3     U3    
         vnodeA=vnode2                                              !      .     |      . \   |      .    
         unodeB=unode4                                              !      .     |      .Ri\Re|      .    
         ! alpha=(yB-yA)/(xB-xA)                                    !      . . . V1 . . . . \ V2 . . .    
         ! beta =yA-alpha*Xa                                        !      .     |      .     |      .    
         alpha=( uy(unodeB)-vy(vnodeA) )/( ux(unodeB)-vx(vnodeA) )  !      .     |      .     |      .    
         beta =vy(vnodeA)- alpha*vx(vnodeA)                         !      .     R1____U1_____R2     U2   
         if(Ypos.gt.alpha*Xpos+beta)   tri_ext=.TRUE.               !      .            .            .    
         Unum= unode1                                               !      .            .            .    
         Vnum= vnode1                                               !      . . . .  . . . . . .  . . .    
       !write(*,'(a,i6,a)')'partonborder[',n-1,']=5   #PYTHON DEBUG'
     !--------------------------------------------------------------------------------------------------------- 
     ELSEIF( (rx(modrnode1)>ux(unode1))        &                    !      . . . V4 . . . . . V3 . . . 
       .and.(ry(modrnode1)< vy(vnode1)))THEN                        !      .            .            . 
         string='R4, V1, U3'                                        !      .            .            . 
         closernode=modrnode4                                       !     U4     R4____U3_____R3     .  
         vnodeA=vnode1                                              !      .     |Re /  .     |      .  
         unodeB=unode3                                              !      .     |  / Ri.     |      .  
         ! alpha=(yB-yA)/(xB-xA)                                    !      . . . V1 . . . . . V2 . . .  
         ! beta =yA-alpha*Xa                                        !      .     |      .     |      .  
         alpha=( uy(unodeB)-vy(vnodeA) )/( ux(unodeB)-vx(vnodeA) )  !      .     |      .     |      .    
         beta =vy(vnodeA)- alpha*vx(vnodeA)                         !     U1     R1____U2_____R2     .    
         if(Ypos.gt.alpha*Xpos+beta) tri_ext=.TRUE.                 !      .            .            .  
         Unum= unode2                                               !      .            .            .  
         Vnum= vnode1                                               !      . . . .  . . . . . .  . . .  
       !write(*,'(a,i6,a)')'partonborder[',n-1,']=7   #PYTHON DEBUG'
     !--------------------------------------------------------------------------------------------------------- 
    ELSE
        stop 'ERROR getclosernode'
    ENDIF
    !-------------------------------------------------------------
    if(conflictin<0)   write(*,*)' HCONFLICT str=', string
    IF (tri_ext)then
      kbottom=BottomKatRnode(closernode)
      !write(*,'(a,i6,a)')'partonborder[',n-1,']+=1  #PYTHON DEBUG'
    ELSE
      i=mod(modrnode1,vi)
      j=(modrnode1-i)/vi+1
      kbottom=1
      if(closernode.eq.modrnode1)then
        do k=us,1,-1
             select case(eleform(i,j,k))        ! 03 ! 07 ! 10 ! 11 ! 13 ! 14 ! 15 !  17  !  18      ! 
               case(3,7,10,11,13,14,15,17,18)   ! x. ! xx ! .. ! x. ! xx ! .x ! xx ! ./x/ ! \x\.     !
                kbottom=k+1                     ! x. ! x. ! xx ! xx ! .x ! xx ! xx ! /x/. ! .\x\     !    
                exit
               case default
                cycle
             end select 
          enddo
       elseif(closernode.eq.modrnode2)then
           do k=us,1,-1
              select case(eleform(i,j,k))           ! 07 ! 10 ! 11 ! 12 ! 13 ! 14 ! 15 !  17  !  18        !  
                case(7,10,11,12,13,14,15,17,18)     ! xx ! .. ! x. ! .x ! xx ! .x ! xx ! ./x/ ! \x\.       !
                 kbottom=k+1                        ! x. ! xx ! xx ! .x ! .x ! xx ! xx ! /x/. ! .\x\       !
                 exit
                case default
                 cycle
              end select 
           enddo
       elseif(closernode.eq.modrnode3)then
           do k=us,1,-1
              select case(eleform(i,j,k))         ! 05 ! 07 ! 11 ! 12 ! 13 ! 14 ! 15 !  17  !  18       ! 
                case(5,7,11,12,13,14,15,17,18)    ! xx ! xx ! x. ! .x ! xx ! .x ! xx ! ./x/ ! \x\.      !
                 kbottom=k+1                      ! .. ! x. ! xx ! .x ! .x ! xx ! xx ! /x/. ! .\x\      !
                 exit
                case default
                 cycle
              end select 
           enddo
       else
          do k=us,1,-1
             select case(eleform(i,j,k))         ! 03 ! 05 ! 07 ! 11 ! 13 ! 14 ! 15 !  17  !  18        ! 
               case(3,5,7,11,13,14,15,17,18)     ! x. ! xx ! xx ! x. ! xx ! .x ! xx ! ./x/ ! \x\.       !
                kbottom=k+1                      ! x. ! .. ! x. ! xx ! .x ! xx ! xx ! /x/. ! .\x\       !
                 exit
               case default
                cycle
             end select 
           enddo
       endif
    ENDIF
    !-------------------------------------------------------------

    IF (Zgrid_depthinterp)then

      x3=vx(vnodeA)                                      
      y3=vy(vnodeA)                                      
      v3=depthV(vnodeA)                                
      x2=ux(unodeB)                                    
      y2=uy(unodeB)                                       
      v2=depthU(unodeB)
      if(tri_ext)then                    !                
        x1=rx(closernode)                ! Re (external)    
        y1=ry(closernode)                !                  
        v1=depthR(closernode)
      else                               !                 
        x1=0.5*(ux(Unum)+ux(Unum+ui))   ! Ri (internal)   
        y1=0.5*(vy(Vnum)+vy(Vnum+1))   !
        v1=depthE(Unum)
      endif  

      ! bilinear interpolation of triangle
      tt = ((Xpos-x1)*(y3-y1)+(y1-Ypos)*(x3-x1)) &
               / ((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))
      uu = ((Xpos-x1)*(y2-y1)+(y1-Ypos)*(x2-x1)) &
               / ((x3-x1)*(y2-y1)-(y3-y1)*(x2-x1))
      vp = v1 + (v2-v1)*tt + (v3-v1)*uu

      !if bilinear techniques are undefined
      if( tt.LT.0. .OR. uu.LT.0. .OR. (tt+uu).GT.1.0 ) then
        !if particle on node, then set equal to node value
        if( (Xpos.EQ.x1 .AND. Ypos.EQ.y1).OR.(Xpos.EQ.x2 .AND. Ypos.EQ.y2) &
        .OR.(Xpos.EQ.x3 .AND. Ypos.EQ.y3))then
          if (Xpos.EQ.x1 .AND. Ypos.EQ.y1)  vp=v1
          if (Xpos.EQ.x2 .AND. Ypos.EQ.y2)  vp=v2
          if (Xpos.EQ.x3 .AND. Ypos.EQ.y3)  vp=v3
        else
        ! write(*,'(a,i6,a)')'partonborder[',n-1,']+=10  #PYTHON DEBUG'
          Dis1=1./( SQRT( (x1-Xpos)**2 + (y1-Ypos)**2 ) ) 
          Dis2=1./( SQRT( (x2-Xpos)**2 + (y2-Ypos)**2 ) ) 
          Dis3=1./( SQRT( (x3-Xpos)**2 + (y3-Ypos)**2 ) )
          !!trilatlon(1,1)=x1
          !!trilatlon(1,2)=y1
          !!trilatlon(2,1)=x2
          !!trilatlon(2,2)=y2
          !!trilatlon(3,1)=x3
          !!trilatlon(3,2)=y3
          !!trilatlon(4,1)=x1
          !!trilatlon(4,2)=y1
          !!if(inpoly(Xpos,Ypos,3,trilatlon)) then
          !! write(*,'(8(a,f14.6),2a)')'INPOLY=[[',Xpos,',',Ypos,&
          !!'],[',x1,',',y1,'],[',x2,',',y2,'],[',x3,',',y3,']] #',string 
          !aA =(y2-y1)    /(x2-x1)                                                    ! wrong out of the triangle 
          !bA =    y1-aA      *x1   
          !YpA=    aA*Xpos    +bA
          !aB =(y3-y2)    /(x3-x2)
          !bB =    y2-aB      *x2
          !YpB=    aB*Xpos    +bB
          !aC =(y1-y3)    /(x1-x3)
          !bC =    y3-aC      *x3
          !YpC=    aC*Xpos    +bC
          !if(abs(Ypos-YpA)<abs(Ypos-YpB).and.  abs(Ypos-YpA)<abs(Ypos-YpC))then
          !  Tdis=    1./( SQRT( (x1-Xpos)**2 + (y1-Ypos)**2 ) ) + &
          !           1./( SQRT( (x2-Xpos)**2 + (y2-Ypos)**2 ) )
          !  Dis3= 0.0
          !elseif(abs(Ypos-YpB)<abs(Ypos-YpA).and.  abs(Ypos-YpB)<abs(Ypos-YpC))then
          !  Tdis=    1./( SQRT( (x2-Xpos)**2 + (y2-Ypos)**2 ) ) + &
          !           1./( SQRT( (x3-Xpos)**2 + (y3-Ypos)**2 ) )
          !  Dis1= 0.0
          !elseif(abs(Ypos-YpC)<abs(Ypos-YpA).and.  abs(Ypos-YpC)<abs(Ypos-YpB))then
          !  Tdis=    1./( SQRT( (x1-Xpos)**2 + (y1-Ypos)**2 ) ) + &
          !           1./( SQRT( (x3-Xpos)**2 + (y3-Ypos)**2 ) )
          !  Dis2= 0.0
          !else !use inverse weighted distance instead  
          !  TDis = Dis1+Dis2+Dis3
          !endif
          !!else
          !! write(*,'(8(a,f14.6),2a)')'NOPOLY=[[',Xpos,',',Ypos,  &
          !! '],[',x1,',',y1,'],[',x2,',',y2,'],[',x3,',',y3,']] #',string
          !!endif
          TDis = Dis1+Dis2+Dis3
          Wt1= Dis1/TDis
          Wt2= Dis2/TDis
          Wt3= Dis3/TDis
          vp = Wt1*v1 + Wt2*v2 + Wt3*v3 
        endif
      endif
      Pdepth = -1.0*vp
       if(conflictin<0)then
         write(*,*)' HCONFLICT VnodeA ',vnodeA,x2lon(x1,y1),y2lat(y1),v1
         write(*,*)' HCONFLICT Unodeb ',unodeb,x2lon(x2,y2),y2lat(y2),v2
         if(tri_ext)then    
          write(*,*)' HCONFLICT Rnode ',closernode,x2lon(x3,y3),y2lat(y3),v3
         else                               !                 
          write(*,*)' HCONFLICT Enode ',Unum,x2lon(x3,y3),y2lat(y3),v3
         endif 
         write(*,*)' HCONFLICT tt,uu,tt+uu=',tt,uu,tt+uu
       endif
     ELSE

       if(tri_ext)then
        Pdepth=-1.0*depth(closernode)
       else
        !Pdepth=ZW(kbottom)  ! CHECK THAT THIS IS RIGHT.... !!
        Pdepth=-1.0*max(max(depthV(vnodeA),depthU(unodeB)),depthE(Unum))
        !if(Pdepth.ne.1.0*ZW(kbottom))then
        !   write(*,*)'check getDepth in hydromodule',kbottom,ZW(kbottom),&
        !           depthV(vnodeA),depthU(unodeB),depthE(Unum)
        !   write(*,*) closernode,rnode1,rnode3,rnode3,rnode4
        !   stop 'ERROR stop program'
        !endif
       endif
     ENDIF
    !-------------------------------------------------------------


    IF( kbottom>us      )then
     write(*,*)'-------------------------------------------------------------'
         If(k>us .or. j>uj .or. i>vi)then
             write(title,'(2(a,i10),3a,3(i5,a),i14)')&
              'n=',n,', idt=',idt,', ',&
               trim(string),&
             ', i=',i,', j=',j,', k=',k,', modrnode1=',modrnode1
         Else
         if(conflict==1)then
            write(title,'(2(a,i10),3a,3(i5,a),i14)')&
              'n=',n,', idt=',idt,', ',&
               trim(string),&
             
             ', i=',i,', j=',j,', k=',k,', modrnode1=',modrnode1
         else
            write(title,'(2(a,i10),2a,2(a,i4))')&
              'n=',n,', idt=',idt,', conflict ',&
             trim(string),& 
             ' delta_i=',delta_i,' delta_j=',delta_j
         endif
         Endif
         if(Pdepth==-999)then
             marker="*"
         else 
             marker="o"
         endif 
         conflict=-1
    ENDIF
    if(kbottom>us)then
     write(*,*)'kbottom=',kbottom,'  > us=',us
     write(*,*)'STOP PARTICLE'
     conflict=-1
    endif
    if(kbottom>us .or. conflictin<0)then
     write(*,*)' HCONFLICT string=',string
     write(*,*)' HCONFLICT kbottom,conflictin=',kbottom,conflictin
     write(*,*)' HCONFLICT Lon,Lat=',x2lon(Xpos,Ypos),y2lat(Ypos)
     write(*,*)' HCONFLICT conflict = ',conflict
     write(*,*)' HCONFLICT tri_ext=',tri_ext
     write(*,*)' HCONFLICT closernode,vnodeA,unodeB,',closernode,vnodeA,unodeB
     write(*,*)' HCONFLICT modrnode1,modrnode2,modrnode3,modrnode4',&
                modrnode1,modrnode2,modrnode3,modrnode4
     write(*,*)' HCONFLICT rnode1,rnode2,rnode3,rnode4',&
                rnode1,rnode2,rnode3,rnode4
     write(*,*)' HCONFLICT unode1,unode2,unode3,unode4',unode1,unode2,unode3,unode4
     write(*,*)' HCONFLICT vnode1,vnode2,vnode3,vnode4',vnode1,vnode2,vnode3,vnode4
     write(*,*)' HCONFLICT DepthRefRnode=',DepthRefRnode,' closernode=',closernode
     write(*,*)' HCONFLICT Ypos=',Ypos,' <=> alpha*Xpos+beta =',alpha*Xpos+beta
           i=mod(modrnode1,vi)                      
           j=(modrnode1-i)/vi+1                    
           write(*,*)'modrnode1=',modrnode1,' vi=',vi,'; i,j=',i,j 
           kbot=1
           do k=us,1,-1
              write(*,*)'k=',k,'eleform=',eleform(i,j,k)
              select case(eleform(i,j,k))         ! 03 ! 05 ! 07 ! 11 ! 13 ! 14 ! 15 !  17  !  18    !
                case(3,5,7,11,13,14,15,17,18)     ! .. ! .x ! .x ! x. ! xx ! xx ! xx ! ./x/ ! \x\.   !
                 kbot=k+1                         ! xx ! .x ! xx ! xx ! .x ! x. ! xx ! /x/. ! .\x\   !
                 exit
                case default
                 cycle
              end select 
           enddo
           write(*,*)'kbot=',kbot
   endif 
  END SUBROUTINE getDepth

  INTEGER FUNCTION getKRlevel(Zin)
    !This function returns the vertical RHO level of the Z value given in input
    USE PARAM_MOD, ONLY: us, ws, Zgrid
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: Zin
    INTEGER :: k

    if(Zgrid)then
      getKRlevel=0
      if ( Zin.LT.ZW(2) ) then
        getKRlevel = 1
      elseif( Zin.GE.ZW(ws-1) ) then
        getKRlevel = ws-1
      else
        do k=2,ws-2
           if ( ( Zin .GE. ZW(k) ) .and. ( Zin .LT. ZW(k+1) ) )then
             getKRlevel = k
             exit
           endif
        enddo
      endif
     
      if (getKRlevel==0) stop "ERROR computng KRlevel"
    else
      getKRlevel=1
    endif
  END FUNCTION getKRlevel

  DOUBLE PRECISION FUNCTION getSlevel(zeta,depth,i)
    !This function returns the depth of the current s-level
    USE PARAM_MOD, ONLY: hc,Vtransform,us
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    DOUBLE PRECISION, INTENT(IN) :: zeta,depth

    DOUBLE PRECISION :: S,h

    ! convert negative depth to positive depth
    h = DBLE(-1.0) * depth


    SELECT CASE(Vtransform)

      CASE(0)  ! MITgcm output contains data at z-coordinate nodes

        getSlevel = ZC(i)
        if (i==us) getSlevel= getSlevel+0.5*zeta

      CASE(1)  !Rutgers-ROMS formulation, eqn (1) of 
        !https://www.myroms.org/wiki/index.php/Vertical_S-coordinate

        S = hc*SC(i)+(h-hc)*CS(i)
        getSlevel = S+zeta*(DBLE(1.0)+S/h)

      CASE(2)  !UCLA-formulation, eqn(2) of 
        !https://www.myroms.org/wiki/index.php/Vertical_S-coordinate

        S = (hc*SC(i)+h*CS(i)) / (hc+h)
        getSlevel = zeta+(zeta+h)*S

      CASE(3)  !Song, Y. and D. B. Haidvogel, 1994: A semi-implicit
        !ocean circulation model using a generalized topography-following 
        !coordinate system, J. Comp. Phys., 115 (1), 228-244.

        getSlevel = zeta*(DBLE(1.0)+SC(i))+hc*SC(i)+(h-hc)*CS(i)

      CASE(4)                    ! 4: z*=H(z-Eta)/(H+Eta) Where H is bottom
                                 ! depth and Eta is the sea surface elevation
        getSlevel = h*(ZC(i)-zeta)/(h+zeta)

      CASE DEFAULT
        write(*,*) 'ERROR: Illegal Vtransform number'
        write(*,*) ' '
        write(*,*) 'The Program Cannot Continue and Will Terminate'
        stop

    END SELECT

  END FUNCTION getSlevel



  DOUBLE PRECISION FUNCTION getWlevel(zeta,depth,i)
    !This function returns the depth of the current w s-level
    USE PARAM_MOD, ONLY: hc,Vtransform,ws
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    DOUBLE PRECISION, INTENT(IN) :: zeta,depth

    DOUBLE PRECISION :: S,h

    ! convert negative depth to positive depth
    h = DBLE(-1.0) * depth


    SELECT CASE(Vtransform)

      CASE(0)  ! MITgcm output contains data at z-coordinate nodes

        getWlevel = ZW(i)
        if (i==ws) getWlevel= getWlevel+zeta

      CASE(1)  !Rutgers-ROMS formulation, eqn (1) of 
        !https://www.myroms.org/wiki/index.php/Vertical_S-coordinate

        S = hc*SCW(i)+(h-hc)*CSW(i)
        getWlevel = S+zeta*(DBLE(1.0)+S/h)

      CASE(2)  !UCLA-formulation, eqn(2) of 
        !https://www.myroms.org/wiki/index.php/Vertical_S-coordinate

        S = (hc*SCW(i)+h*CSW(i))/(hc+h)
        getWlevel = zeta+(zeta+h)*S

      CASE(3)  !Song, Y. and D. B. Haidvogel, 1994: A semi-implicit
        !ocean circulation model using a generalized topography-following 
        !coordinate system, J. Comp. Phys., 115 (1), 228-244.

        getWlevel = zeta*(DBLE(1.0)+SCW(i))+hc*SCW(i)+(h-hc)*CSW(i)

      CASE DEFAULT
        write(*,*) 'ERROR: Illegal Vtransform number'
        write(*,*) ' '
        write(*,*) 'The Program Cannot Continue and Will Terminate'
        stop

    END SELECT

  END FUNCTION getWlevel


  SUBROUTINE getMask_Rho(mask)
    !This subroutine returns the values in the variable mask_rho
    !This is used by createBounds() in the boundary module to make the 
    !  boundaries based on mask_rho
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(OUT) :: mask(:,:,:)

    if(GRD_SET)then
      mask = mask_rho
    else
      write(*,*) 'ERROR: Cannot create boundaries, mask_rho not yet read in'
      write(*,*) ' '
      write(*,*) 'The Program Cannot Continue and Will Terminate'
      stop
    endif

  END SUBROUTINE getMask_Rho

  SUBROUTINE getUVxy(ux,uy,vx,vy)
    !This subroutine returns the values in the variables x_u,y_u,x_v,y_v
    !This is used by createBounds() in the boundary module to make the 
    !  boundaries on the U & V node locations
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(OUT) :: ux(:,:),uy(:,:),vx(:,:),vy(:,:)

    if(GRD_SET)then
      ux = x_u
      uy = y_u
      vx = x_v
      vy = y_v
    else
      write(*,*) 'ERROR: Cannot create boundaries, x_u, y_u, x_v, or y_v is ', &
                 'not yet read in'
      write(*,*) ' '
      write(*,*) 'The Program Cannot Continue and Will Terminate'
      stop
    endif

  END SUBROUTINE getUVxy

  SUBROUTINE getR_ele(k,ele_x,ele_y)
    !This subroutine returns the values in the variables r_kwele_x, and r_kwele_y
    !This is used by createPolySpecs() in the settlement module to determine 
    !  which habitat polygons are in each element
    USE PARAM_MOD,    ONLY: us
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: k
    DOUBLE PRECISION, INTENT(OUT) :: ele_x(:,:),ele_y(:,:)

    if(GRD_SET)then
      ele_x = r_kwele_x(:,:,k)
      ele_y = r_kwele_y(:,:,k)
    else
      write(*,*)'ERROR: Cannot create Poly Specs, r_kwele_x or r_kwele_y ',    &
                 'is not yet created'
      write(*,*) ' '
      write(*,*) 'The Program Cannot Continue and Will Terminate'
      stop
    endif

  END SUBROUTINE getR_ele

  INTEGER FUNCTION getP_r_element(n)
    !This subroutine returns the id of the rho element the particle is 
    !  currently in
    !This is used by settlement() to determine which habitat polygons to check
    !  for settlement
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    getP_r_element = P_r_element(n)
  END FUNCTION getP_r_element


  INTEGER FUNCTION getP_klev(n)
    !This subroutine returns the k-level in which is particle n, determined in
    ! set_Ele or set_Ele_all
    !This is used by settlement() to determine which habitat polygons to check
    !  for settlement
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    getP_klev = P_klev(n)
  END FUNCTION getP_klev


!!SUBROUTINE GetCoastDist(n,coastdist,k)
!!  IMPLICIT NONE
!!  INTEGER, INTENT(IN) :: n,k
!!  INTEGER, INTENT(OUT):: coastdist
!!  coastdist=coast_ele(P_r_element(n),k)

!!END SUBROUTINE GetCoastDist


  SUBROUTINE setijruv()
    USE PARAM_MOD, ONLY: numpar,vi,ui,vj,uj,ijbuff
    USE CONVERT_MOD, ONLY: x2lon,y2lat
    IMPLICIT NONE

    INTEGER :: i,j,n,n1

    n1=1
    do n=1,numpar
      if(.not.DeadOrOut(n))then
        n1=n
        exit
      endif
    enddo

    !rho
    i = mod(RE(1,P_r_element(n1),P_klev(n1))-1,vi)+1
    j = (RE(1,P_r_element(n1),P_klev(n1))-1)/vi + 1
    t_ijruv(IMIN,RNODE) = max(i-ijbuff,1)
    t_ijruv(IMAX,RNODE) = min(i+ijbuff+1,vi)
    t_ijruv(JMIN,RNODE) = max(j-ijbuff,1)
    t_ijruv(JMAX,RNODE) = min(j+ijbuff+1,uj)

    do n=n1+1,numpar
      if(DeadOrOut(n)) cycle   
      i = mod(RE(1,P_r_element(n),P_klev(n))-1 ,vi)+ 1
      j = (RE(1,P_r_element(n),P_klev(n))-1)/vi + 1
      if((i-ijbuff  ) < t_ijruv(IMIN,RNODE)) t_ijruv(IMIN,RNODE) = max(i-ijbuff,1)
      if((i+ijbuff+1) > t_ijruv(IMAX,RNODE)) t_ijruv(IMAX,RNODE) = min(i+ijbuff+1,vi)
      if((j-ijbuff  ) < t_ijruv(JMIN,RNODE)) t_ijruv(JMIN,RNODE) = max(j-ijbuff,1)
      if((j+ijbuff+1) > t_ijruv(JMAX,RNODE)) t_ijruv(JMAX,RNODE) = min(j+ijbuff+1,uj)
    enddo

    !u
    i = mod(UE(1,P_u_element(n1),P_klev(n1))-1,ui)+1
    j = (UE(1,P_u_element(n1),P_klev(n1))-1)/ui + 1
    t_ijruv(IMIN,UNODE) = max(i-ijbuff,1)
    t_ijruv(IMAX,UNODE) = min(i+ijbuff+1,ui)
    t_ijruv(JMIN,UNODE) = max(j-ijbuff,1)
    t_ijruv(JMAX,UNODE) = min(j+ijbuff+1,uj)

    do n=n1+1,numpar
      if(DeadOrOut(n)) cycle   
      i = mod(UE(1,P_u_element(n),P_klev(n))-1,ui)+1
      j = (UE(1,P_u_element(n),P_klev(n))-1)/ui + 1
      if((i-ijbuff  ) < t_ijruv(IMIN,UNODE)) t_ijruv(IMIN,UNODE) = max(i-ijbuff,1)
      if((i+ijbuff+1) > t_ijruv(IMAX,UNODE)) t_ijruv(IMAX,UNODE) = min(i+ijbuff+1,ui)
      if((j-ijbuff  ) < t_ijruv(JMIN,UNODE)) t_ijruv(JMIN,UNODE) = max(j-ijbuff,1)
      if((j+ijbuff+1) > t_ijruv(JMAX,UNODE)) t_ijruv(JMAX,UNODE) = min(j+ijbuff+1,uj)
    enddo

    !v
    i = mod(VE(1,P_v_element(n1),P_klev(n1))-1,vi)+1
    j = (VE(1,P_v_element(n1),P_klev(n1))-1)/vi + 1
    t_ijruv(IMIN,VNODE) = max(i-ijbuff,1)
    t_ijruv(IMAX,VNODE) = min(i+ijbuff+1,vi)
    t_ijruv(JMIN,VNODE) = max(j-ijbuff,1)
    t_ijruv(JMAX,VNODE) = min(j+ijbuff+1,vj)

    do n=n1+1,numpar
      if(DeadOrOut(n)) cycle   
      i = mod(VE(1,P_v_element(n),P_klev(n))-1,vi)+1
      j = (VE(1,P_v_element(n),P_klev(n))-1)/vi + 1
      if((i-ijbuff  ) < t_ijruv(IMIN,VNODE)) t_ijruv(IMIN,VNODE) = max(i-ijbuff,1)
      if((i+ijbuff+1) > t_ijruv(IMAX,VNODE)) t_ijruv(IMAX,VNODE) = min(i+ijbuff+1,vi)
      if((j-ijbuff  ) < t_ijruv(JMIN,VNODE)) t_ijruv(JMIN,VNODE) = max(j-ijbuff,1)
      if((j+ijbuff+1) > t_ijruv(JMAX,VNODE)) t_ijruv(JMAX,VNODE) = min(j+ijbuff+1,vj)
    enddo


  END SUBROUTINE setijruv


  SUBROUTINE finHydro()
    USE PARAM_MOD, ONLY:Zgrid,read_GrainSize,OutDir,NCOutFile,&
           Zgrid_depthinterp,WindIntensity
    !This subroutine closes all the module's allocatable variables
    IMPLICIT NONE

    !ALLOCATE MODULE VARIABLES
    DEALLOCATE(r_Adjacent,u_Adjacent,v_Adjacent)
    DEALLOCATE(mask_rho,depth)
    if(Zgrid)then
      DEALLOCATE(ZC,ZW)
      DEALLOCATE(BottomK)
      DEALLOCATE(BottomKatRnode,BottomKatUnode,BottomKatVnode)
      DEALLOCATE(depthU)
      DEALLOCATE(depthV)
      DEALLOCATE(depthE)
      DEALLOCATE(depthR)
    else
      DEALLOCATE(SC,CS,SCW,CSW)
    endif
    DEALLOCATE(DeadOrOut)
    DEALLOCATE(rho_angle)
    DEALLOCATE(RE,UE,VE)
    DEALLOCATE(rx,ry,ux,uy,vx,vy)
    DEALLOCATE(rho_kwele,u_kwele,v_kwele)
    DEALLOCATE(r_kwele_x,r_kwele_y)
    DEALLOCATE(u_kwele_x,u_kwele_y)
    DEALLOCATE(v_kwele_x,v_kwele_y)
    DEALLOCATE(x_u,y_u,x_v,y_v)

    DEALLOCATE(P_r_element,P_u_element,P_v_element,P_klev,P_klev_old)
    DEALLOCATE(Xpar_at_setEle,Ypar_at_setEle,Zpar_at_setEle)             !--- CL-OGS
    DEALLOCATE(t_zeta,t_salt,t_temp,t_Wvel,t_Uvel,t_Vvel,t_den,t_KH)
    !DEALLOCATE(t_Swdown)
    DEALLOCATE(coast_eleij)  
    DEALLOCATE(coast_ele) 
    DEALLOCATE(eleform) 
    DEALLOCATE(nbounds,bnd_x,bnd_y)
    DEALLOCATE(updatenodesbuffer)
    DEALLOCATE(t_uwind)
    DEALLOCATE(t_vwind)
    if(WindIntensity .and. Zgrid)DEALLOCATE(t_iwind)
    DEALLOCATE(rho_mask,u_mask,v_mask)
    DEALLOCATE(m_r) 
    DEALLOCATE(m_u) 
    DEALLOCATE(m_v) 
    if(read_GrainSize)  DEALLOCATE(GrainSize)  !--- CL-OGS:  for behavior 8
    open (unit=fpy,file=TRIM(OutDir)//'/'//TRIM(NCOutFile)//'PartinEle.py',    &
                 POSITION='APPEND')
    write(fpy,'(a)')"plt.show()"
    close(fpy) 
    DEALLOCATE(lastplotpos)
    if(Zgrid)DEALLOCATE(Coef_COPNOD,Nghb_COPNOD,Node_COPNOD,klev_COPNOD)
    DEALLOCATE(OMP_ruv)
    DEALLOCATE(OMP_Xpar,OMP_Ypar)
    DEALLOCATE(OMP_tOK,OMP_t,OMP_u)
    DEALLOCATE(OMP_Wgt,OMP_Pint)

  END SUBROUTINE finHydro


  SUBROUTINE initNetCDF()
  
    !Initialize NetCDF Counters
    NCcount = 0
    NCstart = 0

  END SUBROUTINE initNetCDF

  SUBROUTINE createNodesDepthNetCDF()
    USE PARAM_MOD, ONLY: NCOutFile,outpath,outpathGiven,   &
        ExeDir,OutDir, ui,uj,vi,vj
   
   USE CONVERT_MOD, ONLY: x2lon,y2lat                             
    USE netcdf
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
    CHARACTER(LEN=300) :: ncFile
    INTEGER :: STATUS,NCID,LonRID,LatRID,LonUID,LatUID,LonVID,LatVID,  &
              LonEID,LatEID,DepthRID,DepthUId,DepthVID,DepthEID,       &
              niRID,njRID,niUID,njVID,count,countU,i,j
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DR,DU,DE,DV
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: LonR,LatR,LonU,LatU, &
                                               LonV,LatV,LonE,LatE
    ALLOCATE(LonR(vi,uj))
    ALLOCATE(LatR(vi,uj))
    ALLOCATE(LonU(ui,uj))
    ALLOCATE(LatU(ui,uj))
    ALLOCATE(LonV(vi,vj))
    ALLOCATE(LatV(vi,vj))
    ALLOCATE(LonE(ui,vj))
    ALLOCATE(LatE(ui,vj))
    ALLOCATE(DR(vi,uj))
    ALLOCATE(DU(ui,uj))
    ALLOCATE(DV(vi,vj))
    ALLOCATE(DE(ui,vj))


    IF(outpathGiven)THEN
        ncFile = TRIM(outpath) // '/' // TRIM(NCOutFile) // '_NodesDepth.nc'
    ELSE
        ncFile = TRIM(NCOutFile) // '_NodesDepth.nc'
    ENDIF

    write(*,*)'Creating NetCDF Output File: ',TRIM(ncFile)

    STATUS = NF90_CREATE(TRIM(ncFile), NF90_CLOBBER, NCID)
    IF(STATUS /= NF90_NOERR) THEN
      WRITE(*,*) 'Problem creating NetCDF output file'
      WRITE(*,*) NF_STRERROR(STATUS)
      STOP
    ENDIF

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !NF90_DEF_DIM

        STATUS = NF90_DEF_DIM(NCID,'niR',vi,niRID)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: niR dim'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
        STATUS = NF90_DEF_DIM(NCID,'njR',uj,njRID)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: njR dim'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
        STATUS = NF90_DEF_DIM(NCID,'niU',ui,niUID)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: niU dim'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
        STATUS = NF90_DEF_DIM(NCID,'njV',vj,njVID)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: njV dim'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !NF90_DEF_VAR

    STATUS = NF90_DEF_VAR(NCID,'LonR',NF_DOUBLE,(/niRID,njRID/),LonRID)
    IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: LonR var'
    IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

    STATUS = NF90_DEF_VAR(NCID,'LatR',NF_DOUBLE,(/niRID,njRID/),LatRID)
    IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: LatR var'
    IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

    STATUS = NF90_DEF_VAR(NCID,'LonU',NF_DOUBLE,(/niUID,njRID/),LonUID)
    IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: LonU var'
    IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

    STATUS = NF90_DEF_VAR(NCID,'LatU',NF_DOUBLE,(/niUID,njRID/),LatUID)
    IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: LatU var'
    IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

    STATUS = NF90_DEF_VAR(NCID,'LonV',NF_DOUBLE,(/niRID,njVID/),LonVID)
    IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: LonV var'
    IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

    STATUS = NF90_DEF_VAR(NCID,'LatV',NF_DOUBLE,(/niRID,njVID/),LatVID)
    IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: LatV var'
    IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

    STATUS = NF90_DEF_VAR(NCID,'LonE',NF_DOUBLE,(/niUID,njVID/),LonEID)
    IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: LonE var'
    IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

    STATUS = NF90_DEF_VAR(NCID,'LatE',NF_DOUBLE,(/niUID,njVID/),LatEID)
    IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: LatE var'
    IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


    STATUS = NF90_DEF_VAR(NCID,'DepthR',NF_DOUBLE,(/niRID,njRID/),DepthRID)
    IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: DepthR var'
    IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

    STATUS = NF90_DEF_VAR(NCID,'DepthE',NF_DOUBLE,(/niUID,njVID/),DepthEID)
    IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: DepthR var'
    IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

    STATUS = NF90_DEF_VAR(NCID,'DepthU',NF_DOUBLE,(/niUID,njRID/),DepthUID)
    IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: DepthR var'
    IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

    STATUS = NF90_DEF_VAR(NCID,'DepthV',NF_DOUBLE,(/niRID,njVID/),DepthVID)
    IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: DepthR var'
    IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !NF_ENDDEF

      STATUS = NF90_ENDDEF(NCID)
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: EndDef'
      IF(status /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !NF_PUT_VAR

     count=0
     countU=0
     do j=1,uj
       do i=1,vi
         count = count+1
         LonR(i,j)=x2lon(rx(count),ry(count)) 
         LatR(i,j)=y2lat(          ry(count))
         DR(i,j)=depthR(count) 
         if(i<vi .and. j<uj)then                     
            countU=countU+1                       
            DE(i,j)=depthE(countU)
            LonE(i,j)=x2lon(0.5*(ux(countU)+ux(countU+ui)), &
                            0.5*(vy(count )+vy(count +1 )) ) 
            LatE(i,j)=y2lat(0.5*(vy(count )+vy(count +1 )) )
         endif
       enddo
     enddo
      !Lon   R   
      STATUS = NF90_INQ_VARID(NCID, "LonR",            LonRID)
      STATUS = NF90_PUT_VAR(NCID,    LonRID,           LonR)
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put LonR '
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
      !Lat   R   
      STATUS = NF90_INQ_VARID(NCID, "LatR",            LatRID)
      STATUS = NF90_PUT_VAR(NCID,    LatRID,           LatR)
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put LatR '
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
      !Depth R   
      STATUS = NF90_INQ_VARID(NCID, "DepthR",          DepthRID)
      STATUS = NF90_PUT_VAR(NCID,    DepthRID,         DR)
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put DepthR '
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
      !Lon   E   
      STATUS = NF90_INQ_VARID(NCID, "LonE",            LonEID)
      STATUS = NF90_PUT_VAR(NCID,    LonEID,           LonE)
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put LonE '
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
      !Lat   E   
      STATUS = NF90_INQ_VARID(NCID, "LatE",            LatEID)
      STATUS = NF90_PUT_VAR(NCID,    LatEID,           LatE)
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put LatE '
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
      !Depth E   
      STATUS = NF90_INQ_VARID(NCID, "DepthE",          DepthEID)
      STATUS = NF90_PUT_VAR(NCID,    DepthEID,         DE)
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put DepthE '
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

     count=0
     do j=1,uj
       do i=1,ui
         count = count+1
         LonU(i,j)=x2lon(ux(count),uy(count)) 
         LatU(i,j)=y2lat(          uy(count))
         DU(i,j)=depthU(count) 
       enddo
     enddo
      !Lon   U   
      STATUS = NF90_INQ_VARID(NCID, "LonU",            LonUID)
      STATUS = NF90_PUT_VAR(NCID,    LonUID,           LonU)
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put LonU '
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
      !Lat   U   
      STATUS = NF90_INQ_VARID(NCID, "LatU",            LatUID)
      STATUS = NF90_PUT_VAR(NCID,    LatUID,           LatU)
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put LatU '
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
      !Depth U   
      STATUS = NF90_INQ_VARID(NCID, "DepthU",          DepthUID)
      STATUS = NF90_PUT_VAR(NCID,    DepthUID,         DU)
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put DepthU '
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

     count=0
     do j=1,vj
       do i=1,vi
         count = count+1
         LonV(i,j)=x2lon(vx(count),vy(count)) 
         LatV(i,j)=y2lat(          vy(count))
         DV(i,j)=depthV(count) 
       enddo
     enddo
      !Lon   V   
      STATUS = NF90_INQ_VARID(NCID, "LonV",            LonVID)
      STATUS = NF90_PUT_VAR(NCID,    LonVID,           LonV)
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put LonV '
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
      !Lat   V   
      STATUS = NF90_INQ_VARID(NCID, "LatV",            LatVID)
      STATUS = NF90_PUT_VAR(NCID,    LatVID,           LatV)
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put LatV '
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
      !Depth V   
      STATUS = NF90_INQ_VARID(NCID, "DepthV",          DepthVID)
      STATUS = NF90_PUT_VAR(NCID,    DepthVID,         DV)
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put DepthV '
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    !NF_CLOSE

    STATUS = NF_CLOSE(NCID)
    IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: Close'
    IF(status /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)


    DEALLOCATE( DR,DU,DE,DV,        &
               LonR,LatR,LonU,LatU, &
               LonV,LatV,LonE,LatE )

  END SUBROUTINE createNodesDepthNetCDF

  SUBROUTINE createNetCDF(dob)
    USE PARAM_MOD, ONLY: numpar,NCOutFile,outpath,outpathGiven,NCtime,         &
       SVN_Version,RunName,ExeDir,OutDir,RunBy,Institution,StartedOn,          &
       TrackCollisions,SaltTempOn,Behavior,Write_coastdist,Write_Poly_Presence
    USE netcdf
    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(numpar), OPTIONAL, INTENT(IN) :: dob

    INCLUDE 'netcdf.inc'

    CHARACTER(LEN=300) :: ncFile
    INTEGER :: STATUS,NCID,numparID,timeID,pageID,modtimeID,lonID,latID,       &
               depthID,colorID,hitBID,hitLID,dobID,saltID,tempID,              &
               grsizeID,psizeID,CoDiID,PolyID

    !   NF90_CREATE           ! create netCDF dataset: enter define mode
    !        ...
    !      NF90_DEF_DIM       ! define dimensions: from name and length
    !        ...
    !      NF90_DEF_VAR       ! define variables: from name, type, dims
    !        ...
    !      NF90_PUT_ATT       ! assign attribute values
    !        ...
    !   NF90_ENDDEF           ! end definitions: leave define mode
    !        ...
    !      NF90_PUT_VAR       ! provide values for variable
    !        ...
    !   NF90_CLOSE            ! close: save new netCDF dataset


    !NF90_CREATE

    !Reset Print Counter to 0
    prcount = 0

    IF(outpathGiven)THEN
      IF(NCtime == 0 ) THEN
        ncFile = TRIM(outpath) // '/' // TRIM(NCOutFile) // '.nc'
      ELSE
        NCcount = NCcount + 1
        write(ncFile,"(A,A,A,A,I3.3,A)")TRIM(outpath),'/',TRIM(NCOutFile),'_', &
                                      NCcount,'.nc'
      ENDIF
    ELSE
      IF(NCtime == 0 ) THEN
        ncFile = TRIM(NCOutFile) // '.nc'
      ELSE
        NCcount = NCcount + 1
        write(ncFile,"(A,A,I3.3,A)")TRIM(NCOutFile),'_',NCcount,'.nc'
      ENDIF
    ENDIF

    write(*,*)'Creating NetCDF Output File: ',TRIM(ncFile)

    STATUS = NF90_CREATE(TRIM(ncFile), NF90_CLOBBER, NCID)
    IF(STATUS /= NF90_NOERR) THEN
      WRITE(*,*) 'Problem creating NetCDF output file'
      WRITE(*,*) NF_STRERROR(STATUS)
      STOP
    ENDIF

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !NF90_DEF_DIM

        STATUS = NF90_DEF_DIM(NCID,'numpar',numpar,numparID)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: numpar dim'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_DEF_DIM(NCID,'time',NF90_UNLIMITED,timeID)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: time dim'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !NF90_DEF_VAR

        STATUS = NF90_DEF_VAR(NCID,'model_time',NF_DOUBLE,(/timeID/),modtimeID)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: time var'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        IF( PRESENT(dob) )THEN
          STATUS = NF90_DEF_VAR(NCID,'dob',NF_FLOAT,(/numparID/),dobID)
          IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: dob var'
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
        ENDIF

        STATUS = NF90_DEF_VAR(NCID,'age',NF_INT,(/numparID,timeID/),pageID)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: age var'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_DEF_VAR(NCID,'lon',NF_FLOAT,(/numparID,timeID/),lonID)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: lon var'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_DEF_VAR(NCID,'lat',NF_FLOAT,(/numparID,timeID/),latID)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: lat var'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_DEF_VAR(NCID,'depth',NF_FLOAT,(/numparID,timeID/),      &
                              depthID)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: depth var'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_DEF_VAR(NCID,'color',NF_INT,(/numparID,timeID/),      &
                              colorID)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: color var'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        IF(TrackCollisions)THEN
          STATUS =NF90_DEF_VAR(NCID,'hitBottom',NF_INT,(/numparID,timeID/), &
                               hitBID)
          IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: Bottom ', &
                                              'Collision var'
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_DEF_VAR(NCID,'hitLand',NF_INT,(/numparID,timeID/),  &
                                hitLID)
          IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: Land ',   &
                                              'Collision var'
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
        ENDIF

        IF(SaltTempOn)THEN
          STATUS = NF90_DEF_VAR(NCID,'salinity',NF_FLOAT,(/numparID,timeID/), &
                                saltID)
          IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
                                              'Salinity var'
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_DEF_VAR(NCID,'temperature',NF_FLOAT,                  &
                                (/numparID,timeID/),tempID)
          IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
                                              'Temperature var'
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
        ENDIF
        
        IF(Behavior.ge.8 .and. Behavior.le. 11)THEN
          STATUS = NF90_DEF_VAR(NCID,'GrainSize',NF_FLOAT,(/numparID,timeID/),&
                                grsizeID)
          IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
                                              'GrainSize var'
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_DEF_VAR(NCID,'Size',NF_FLOAT,                         &
                                (/numparID,timeID/),psizeID)
          IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
                                              'Size var'
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
        ENDIF        
        if(Write_coastdist)then
          STATUS = NF90_DEF_VAR(NCID,'CoastDist',NF_FLOAT,(/numparID,timeID/), &
                                CoDiID)
          IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
                                              'CoastDistance var'
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
        endif
        if(Write_Poly_Presence)then
          STATUS = NF90_DEF_VAR(NCID,'Polygon',NF_INT,(/numparID,timeID/), &
                                PolyID)
          IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: ',        &
                                              'Poly var'
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
        endif

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !NF90_PUT_ATT

        !Particle Time
        STATUS = NF90_PUT_ATT(NCID, modtimeID, "long_name",                    &
                              "time that has passed thus far in the model")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, modtimeID, "units", "seconds")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, modtimeID, "field",                        &
                              "model_time, scalar, series")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        IF( PRESENT(dob) )THEN
          !Particle Date of Birth
          STATUS = NF90_PUT_ATT(NCID, dobID, "long_name",                      &
                   "Date of Birth of particles in seconds from model start")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_PUT_ATT(NCID, dobID, "units", "seconds")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_PUT_ATT(NCID, dobID, "field", "age, scalar, series")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
        ENDIF

        !Particle Age
        STATUS = NF90_PUT_ATT(NCID, pageID, "long_name", "age of particles")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, pageID, "units", "seconds")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, pageID, "field", "age, scalar, series")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        !Longitude
        STATUS = NF90_PUT_ATT(NCID, lonID, "long_name",                        &
                              "longitude of particles")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, lonID, "units", "decimal degrees E")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, lonID, "field", "lon, scalar, series")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


        !Latitude
        STATUS =NF90_PUT_ATT(NCID, latID, "long_name", "latitude of particles")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, latID, "units", "decimal degrees N")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, latID, "field", "lat, scalar, series")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


        !Depth
        STATUS = NF90_PUT_ATT(NCID, depthID, "long_name", "depth of particles")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, depthID, "units", "meters below surface")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, depthID, "field", "depth, scalar, series")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


        !Color
        STATUS = NF90_PUT_ATT(NCID, colorID, "long_name",                      &
                 "identification number for particle behavior or status")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, colorID, "units",                          &
                              "nondimensional, see LTRANS User Guide")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, colorID, "field", "color, scalar, series")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


        IF(TrackCollisions)THEN
          write(*,*)'Output includes Hit bottom/land values'
          !hitBottom
          STATUS = NF90_PUT_ATT(NCID, hitBID, "long_name",                     &
                                "# of times Particle Collided with Bottom")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_PUT_ATT(NCID, hitBID, "units", "Number of Collisions")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_PUT_ATT(NCID, hitBID, "field",                         &
                                "hitBottom, scalar, series")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


          !hitLand
          STATUS = NF90_PUT_ATT(NCID, hitLID, "long_name",                     &
                                "# of times Particle Collided with Land")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_PUT_ATT(NCID, hitLID, "units", "Number of Collisions")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_PUT_ATT(NCID, hitLID, "field",                         &
                                "hitLand, scalar, series")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
        ENDIF

        IF(SaltTempOn)THEN
          !salt
          write(*,*)'Output includes Salt values'
          STATUS = NF90_PUT_ATT(NCID, saltID, "long_name",                     &
                                "Salinity at the particle's location")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS=NF90_PUT_ATT(NCID,saltID, "field", "salinity, scalar, series")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


          !temp
          write(*,*)'Output includes Temp values'
          STATUS = NF90_PUT_ATT(NCID, tempID, "long_name",                     &
                                "Temperature at the particle's location")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_PUT_ATT(NCID, tempID, "units", "° Celcius")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_PUT_ATT(NCID, tempID, "field",                         &
                                "temperature, scalar, series")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
        ENDIF
         
        IF(Behavior.ge.8 .and. Behavior.le.11)THEN
          !GrainSize
          write(*,*)'Output includes GrainSize values'
          STATUS = NF90_PUT_ATT(NCID, grsizeID, "long_name",                     &
                "Sediment Grain Size diameter at the particle's location")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_PUT_ATT(NCID, grsizeID, "units", "micrometer")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS=NF90_PUT_ATT(NCID,grsizeID, "field","GrainSize, scalar, series")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


          !size
          write(*,*)'Output includes Part Size values'
          STATUS = NF90_PUT_ATT(NCID, psizeID, "long_name",                     &
                                "Size of the larvae")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_PUT_ATT(NCID, psizeID, "units", "millimeter")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS = NF90_PUT_ATT(NCID, psizeID, "field",                         &
                                "Size, scalar, series")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
        ENDIF


        IF(Write_coastdist)THEN
          write(*,*)'Output includes coastdist values'
          STATUS = NF90_PUT_ATT(NCID, CoDiID, "long_name",                     &
                                "Closest Distance from coast at the particle's location")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS=NF90_PUT_ATT(NCID,CoDiID, "field", "CoastDist, scalar, series")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
        ENDIF
        IF(Write_Poly_Presence)THEN
          write(*,*)'Output includes Poly values'
          STATUS = NF90_PUT_ATT(NCID, PolyID, "long_name",                     &
           "ID number of the polygon where part n spent most of the time")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

          STATUS=NF90_PUT_ATT(NCID,PolyID, "field", "Polygon, scalar, series")
          IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
        ENDIF

        !Global
        STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "type",                       &
                              "Position and characteristics of particles")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID,NF90_GLOBAL,"title","LTRANS-Zlev v0(beta)")
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "svn", SVN_Version)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "run_name", RunName)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "executable", ExeDir)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "output_loc", OutDir)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "administrator", RunBy)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "institution", Institution)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

        STATUS = NF90_PUT_ATT(NCID, NF90_GLOBAL, "date", StartedOn)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !NF_ENDDEF

      STATUS = NF90_ENDDEF(NCID)
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: EndDef'
      IF(status /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !NF_PUT_VAR

      !Particle Date of Birth
      IF( PRESENT(dob) )THEN
        STATUS = NF90_INQ_VARID(NCID, "dob", dobID)
        STATUS = NF90_PUT_VAR(NCID, dobID, dob)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put dob'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)  
      ENDIF

    !NF_CLOSE

    STATUS = NF_CLOSE(NCID)
    IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: Close'
    IF(status /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

  END SUBROUTINE createNetCDF


  SUBROUTINE writeNetCDF(time,age,lon,lat,depth,colors,hitB,hitL,Salt,Temp,    &
                        GrSize,SizeP,CoDi,Poly)
    USE PARAM_MOD, ONLY: numpar,SaltTempOn,NCOutFile,outpath,outpathGiven,     &
        NCtime,TrackCollisions,Behavior,Write_coastdist
    USE netcdf
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: time
    DOUBLE PRECISION, INTENT(IN) :: lon(numpar),lat(numpar),       &
                                    depth(numpar)
    DOUBLE PRECISION,  INTENT(IN) :: age(numpar)
    DOUBLE PRECISION, INTENT(IN) :: colors(numpar)
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: Salt(numpar),Temp(numpar),       &
                                 GrSize(numpar),SizeP(numpar),CoDi(numpar)
   INTEGER, INTENT(IN), OPTIONAL ::  Poly(numpar)
    INTEGER, INTENT(IN), OPTIONAL :: hitB(numpar),hitL(numpar)
    REAL :: realarray(numpar)
    INTEGER :: integerarray(numpar)
    INCLUDE 'netcdf.inc'

    CHARACTER(LEN=300) :: ncFile
    INTEGER :: STATUS,NCID,modtimeID,pageID,lonID,latID,depthID,hitBID,hitLID, &
               colorID,saltID,tempID,grsizeID,psizeID,CoDiID,PolyID
    INTEGER :: NCelapsed

    !If only one NetCDF output file is being written to:
    IF(NCtime == 0) THEN

      IF(outpathGiven)THEN
        ncFile = TRIM(outpath) // '/' // TRIM(NCOutFile) // '.nc'
      ELSE
        ncFile = TRIM(NCOutFile) // '.nc'
      ENDIF

    !If sequentially numbered NetCDF output files are being written to:
    ELSE

      NCelapsed = time - NCstart

      !If specified time interval has been reached, create new NetCDF file
      IF(NCelapsed >= NCtime) THEN
        NCstart = time
        call createNetCDF()
      ENDIF

      IF(outpathGiven)THEN
        write(ncFile,"(A,A,A,A,I3.3,A)")TRIM(outpath),'/',TRIM(NCOutFile),'_', &
                                      NCcount,'.nc'
      ELSE
        write(ncFile,"(A,A,I3.3,A)")TRIM(NCOutFile),'_',NCcount,'.nc'
      ENDIF

    ENDIF

    prcount = prcount + 1

    STATUS = NF90_OPEN(TRIM(ncFile), NF90_WRITE, NCID)
    IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

      !Particle Time
      STATUS = NF90_INQ_VARID(NCID, "model_time", modtimeID)
      STATUS = NF90_PUT_VAR(NCID, modtimeID, DBLE(time), start = (/ prcount /))
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put model_time, time: ',time
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)  

      !Particle Age
      integerarray(:)=int(age(:))
      STATUS = NF90_INQ_VARID(NCID, "age", pageID)
      STATUS = NF90_PUT_VAR(NCID, pageID, integerarray,                       &
                            start = (/ 1, prcount /),                         &
                            count = (/ numpar,  1 /))
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put age, time: ',time
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

      !Longitude
      realarray(:)=real(lon(:))
      STATUS = NF90_INQ_VARID(NCID, "lon", lonID)
      STATUS = NF90_PUT_VAR(NCID, lonID, realarray,                            &
                            start = (/ 1, prcount /),                          &
                            count = (/ numpar,  1 /))
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put lon, time: ',time
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

      !Latitude
      realarray(:)=real(lat(:))
      STATUS = NF90_INQ_VARID(NCID, "lat", latID)
      STATUS = NF90_PUT_VAR(NCID, latID, realarray,                            &
                            start = (/ 1, prcount /),                          &
                            count = (/ numpar,  1 /))
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put lat, time: ',time
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

      !Depth
      realarray(:)=real(depth(:))
      STATUS = NF90_INQ_VARID(NCID, "depth", depthID)
      STATUS = NF90_PUT_VAR(NCID, depthID, realarray,                          &
                            start = (/ 1, prcount /),                          &
                            count = (/ numpar,  1 /))
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put depth, time: ',time
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

      !Color
      integerarray(:)=int(colors(:))
      STATUS = NF90_INQ_VARID(NCID, "color", colorID)
      STATUS = NF90_PUT_VAR(NCID, colorID, integerarray,                       &
                            start = (/ 1, prcount /),                          &
                            count = (/ numpar,  1 /))
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put color, time: ',time
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

      !hitBottom
      IF( PRESENT(hitB) )THEN
        STATUS = NF90_INQ_VARID(NCID, "hitBottom", hitBID)
        STATUS = NF90_PUT_VAR(NCID, hitBID, hitB,                              &
                              start = (/ 1, prcount /),                        &
                              count = (/ numpar,  1 /))
        IF(STATUS /= NF90_NOERR) WRITE(*,*)                                    &
          'Problem put # Bottom Collisions, time: ',time
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
      ENDIF

      !hitLand
      IF( PRESENT(hitL) )THEN
        STATUS = NF90_INQ_VARID(NCID, "hitLand", hitLID)
        STATUS = NF90_PUT_VAR(NCID, hitLID, hitL,                              &
                              start = (/ 1, prcount /),                        &
                              count = (/ numpar,  1 /))
        IF(STATUS /= NF90_NOERR) WRITE(*,*)                                    &
          'Problem put # Land Collisions, time: ',time
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
      ENDIF

      !hitBottom
      IF( PRESENT(salt) )THEN
        realarray(:)=real(salt(:))
        STATUS = NF90_INQ_VARID(NCID, "salinity", saltID)
        STATUS = NF90_PUT_VAR(NCID, saltID, realarray,                         &
                              start = (/ 1, prcount /),                        &
                              count = (/ numpar,  1 /))
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put salinity, time: ',time
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
      ENDIF

      !hitLand
      IF( PRESENT(temp) )THEN
        realarray(:)=real(temp(:))
        STATUS = NF90_INQ_VARID(NCID, "temperature", tempID)
        STATUS = NF90_PUT_VAR(NCID, tempID, realarray,                         &
                              start = (/ 1, prcount /),                        &
                              count = (/ numpar,  1 /))
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put temperature, time: ', &
                                            time
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
      ENDIF
  
      IF( PRESENT(GrSize) )THEN
        realarray(:)=real(GrSize(:))
        STATUS = NF90_INQ_VARID(NCID, "GrainSize", grsizeID)
        STATUS = NF90_PUT_VAR(NCID, grsizeID, realarray,                       &
                              start = (/ 1, prcount /),                        &
                              count = (/ numpar,  1 /))
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put GrainSize, time: ',   &
                                            time
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
      ENDIF

      IF( PRESENT(SizeP) )THEN
        realarray(:)=real(SizeP(:))
        STATUS = NF90_INQ_VARID(NCID, "Size", psizeID)
        STATUS = NF90_PUT_VAR(NCID, psizeID, realarray,                        &
                              start = (/ 1, prcount /),                        &
                              count = (/ numpar,  1 /))
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put Size, time: ',        &
                                            time
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
      ENDIF

      IF( PRESENT(CoDi) .and. Write_coastdist )THEN
        realarray(:)=real(CoDi(:))
        STATUS = NF90_INQ_VARID(NCID, "CoastDist", CoDiID)
        STATUS = NF90_PUT_VAR(NCID, CoDiID, realarray,                         &
                              start = (/ 1, prcount /),                        &
                              count = (/ numpar,  1 /))
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put CoastDist, time: ',time
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
      ENDIF
      IF( PRESENT(Poly) )THEN
        STATUS = NF90_INQ_VARID(NCID, "Polygon", PolyID)
        STATUS = NF90_PUT_VAR(NCID, PolyID, Poly,                              &
                              start = (/ 1, prcount /),                        &
                              count = (/ numpar,  1 /))
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put Polygon, time: ',time
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)
      ENDIF
  
    STATUS = NF_CLOSE(NCID)

  END SUBROUTINE writeNetCDF

   !*********************************************************
   !*                     seteleform                        *
   !*********************************************************
  SUBROUTINE seteleform(ni,nj,nk,eleform_tmp)  ! CL-OGS: created to get Depth of
                                               ! particles in hydro module
   USE PARAM_MOD, ONLY: ui,uj,vi,vj,us
   IMPLICIT NONE
   INTEGER,INTENT(IN):: ni,nj,nk
   INTEGER, INTENT(INOUT):: eleform_tmp(:,:,:)
   INTEGER :: i,j,k
      if((ni.ne.vi-1).or.(nj.ne.uj-1).or.(nk.ne.us_tridim))then
         write(*,*)'OR ni=',ni,' !=? vi-1=',vi-1
         write(*,*)'OR nj=',nj,' !=? uj-1=',uj-1
         write(*,*)'OR nk=',nk,' !=? us_tridim=',us_tridim
         write(*,*)'ERROR DIM eleform sent by boundary mod'
         stop
      endif
      allocate(eleform(vi-1,uj-1,us_tridim))
      eleform(:,:,:)=eleform_tmp
  END SUBROUTINE seteleform


  SUBROUTINE setbndeleform()  
   USE PARAM_MOD, ONLY: ui,uj,vi,vj,us
   IMPLICIT NONE
   INTEGER :: i,j,k

      ! Change eleform on open-sea borders faking land boundaries
      do k=1,us_tridim
        !---------------------
        i=1
        do j=1,uj-1
          select case(eleform(i,j,k))
           case(0,1,2)         !  x|o
             eleform(i,j,k)=3  !__x|o__ 
           case(4,5,6,16,17)   !  x x
             eleform(i,j,k)=7  !__x/o__
           case(8,9,10,18,19)  !  x\o
             eleform(i,j,k)=11 !__x x__
           case(12,13,14)      !  x x
             eleform(i,j,k)=15 !__x x__
          end select 
        enddo
        !---------------------
        i=vi-1
        do j=1,uj-1
          select case(eleform(i,j,k))
           case(0,4,8)          !  o|x     
             eleform(i,j,k)=12  !__o|x__
           case(1,5,9,18,19)    !  x x
             eleform(i,j,k)=13  !__o\x__
           case(2,6,10,16,17)   !  o/x
             eleform(i,j,k)=14  !__x x__
           case(3,7,11)         !  x x
             eleform(i,j,k)=15  !__x x__
          end select 
        enddo
        !---------------------
        j=1
        do i=1,vi-1
          select case(eleform(i,j,k))
           case(0,2,8)          !  o_o   
             eleform(i,j,k)=10  !__x x__
           case(1,3,9,18,19)    !  x\o
             eleform(i,j,k)=11  !__x x__
           case(4,6,12,16,17)   !  o/x
             eleform(i,j,k)=14  !__x x__
           case(5,7,13)         !  x x
             eleform(i,j,k)=15  !__x x__
          end select 
        enddo
        !---------------------
        j=uj-1
        do i=1,vi-1
          select case(eleform(i,j,k))
           case(0,1,4)          !  x_x   
             eleform(i,j,k)= 5  !__o o__
           case(2,3,6,16,17)    !  x x
             eleform(i,j,k)= 7  !__x/o__
           case(8,9,12,18,19)   !  x x
             eleform(i,j,k)=13  !__o\x__
           case(10,11,14)       !  x x
             eleform(i,j,k)=15  !__x x__
          end select 
        enddo
        !---------------------
      enddo
  END SUBROUTINE setbndeleform


   !*********************************************************
   !*                     setnodesdepth                     *
   !*********************************************************

  SUBROUTINE setnodesdepth()
   USE PARAM_MOD, ONLY: ui,uj,vi,vj,us,ws,Zgrid_depthinterp, &
               rho_nodes,u_nodes,v_nodes,BndOut
    USE CONVERT_MOD, ONLY: x2lon,y2lat  !--- CL-OGS
   IMPLICIT NONE
   INTEGER :: i,j,k,cR,cU!,occurences
   DOUBLE PRECISION :: dptR1,dptR2,dptR3,dptR4,dptUa,dptUb,dptVa,dptVb
   DOUBLE PRECISION :: factor(4)

   ALLOCATE(depthR(rho_nodes))    ! number of elements (both land and water)
  
   ALLOCATE(depthE(  u_nodes-ui)) ! number of elements (both land and water)
   ALLOCATE(depthU(  u_nodes))
   ALLOCATE(depthV(  v_nodes))
   write(*,*)'---------------------------------------------------------'
   write(*,*)'setnodesdepth in hydrodynamic module'
    depthR(:)=depth(:)

    cR = 0
    cU = 0                   
    !occurences=0                   
    !write(*,*)'lonlat=[]   #depth!=ZW(kbot)'
    do j=1,uj
      do i=1,vi
        cR = cR + 1   !move to next node number
       !if(depth(cR).ne.-ZW( BottomKatRnode(cR) ))then
       !     if(depth(cR)>-ZW(BottomKatRnode(cR)).or.&
       !        depth(cR)<-ZW(minBottomKatRnode(cR)+1 ))&
       !    !write(*,*)'lonlat=np.hstack((lonlat,(',x2lon(rx(cR),ry(cR)),',', &
       !    !y2lat(ry(cR)),"))) #depth!=ZW(kbot)", &
       !    !depth(cR)," != ",-ZW(BottomKatRnode(cR)),&
       !    ! -ZW(BottomKatRnode(cR)+1 )
       !    !!write(*,*)'partial bottom cell at i,j=',i,j,' depth=',depth(cR), &
       !    !!                 'instead of ',-ZW(BottomKatRnode(cR) )
       !     occurences=occurences+1
       !endif
        !if(j<uj)       depthV(cR)=-ZW( BottomKatVnode(cR) )
        if(j<uj)       depthV(cR)=min(depth(cR),depth(cR+vi))
        if(i<vi)then                          
                       cU=cU+1                       
                       !depthU(cU)=-ZW( BottomKatUnode(cU) )
                       depthU(cU)=min(depth(cR),depth(cR+1))
        endif
      enddo
    enddo
    !write(*,*)"ax6.scatter(lonlat[0::2],lonlat[1::2],c='r'",&
    !         ",s=10,edgecolor='none',alpha=0.5)   #depth!=ZW(kbot)"
    !write(*,*)'found ',occurences,' partial cells '

    !------------------------------------------------------

    IF(Zgrid_depthinterp)then
      write(*,*) 'ADJUSTING depthU and depthV before to ', &                    ! --U(cU+ui-1)--R(cR+vi)---U(cU+ui)  R(cR+vi+1)  
       'set depthE and reinterpolate depthR and depthE'                        !                 |                        |
    cR = 0                                                                      !                 |                        |
    cU = 0                                                                      !   E(cR-1)       V(cR)      E(cR)     V(cR+1)   
    do j=1,uj                                                                   !                 |                        |     
      do i=1,vi                                                                 !                 |                        |     
        cR = cR + 1   !move to next node number                                 ! --U(cU-1)-----  R(cR)------U(cU)-----R(cR+1)   
        if(j<uj)then                                                            !                 |                        |         
          if(depthV(cR)<-ZW( BottomKatVnode(cR) ))then                          !                 |                        |          
           depthV(cR)=min(-ZW( BottomKatVnode(cR) ), &                          !                 V(cR-vi)  E(cR-vi)  V(cR-vi+1)    
              max(depthV(cR),0.5*(depth(cR)+depth(cR+vi))))
          endif
        endif
        if(i<vi)then                          
          cU=cU+1                       
          if(depthU(cU)<-ZW( BottomKatUnode(cU) ))then
           depthU(cU)=min(-ZW( BottomKatUnode(cU) ) , &
              max(depthU(cU),0.5*(depth(cR)+depth(cR+1))))
          endif
        endif
        if((i>1.and.i<vi).and.j.eq.1)then
         if(BottomKatUnode(cU+ui-1)<ws)  depthV(cR)=min(depthV(cR),-ZW(BottomKatUnode(cU+ui-1))) 
         if(BottomKatUnode(cU+ui  )<ws)  depthV(cR)=min(depthV(cR),-ZW(BottomKatUnode(cU+ui  )))  
        elseif((i>1.and.i<vi).and.j.eq.vj)then
         if(BottomKatUnode(cU-1)<ws  )   depthV(cR)=min(depthV(cR),-ZW(BottomKatUnode(cU-1)))
         if(BottomKatUnode(cU)<ws    )   depthV(cR)=min(depthV(cR),-ZW(BottomKatUnode(cU)))
        endif
        if((j>1.and.j<uj) ) then
           if(i.eq.1.)then
                if(depthU(cU)>-ZW(BottomKatVnode(cR+1)) .or. &
                   depthU(cU)>-ZW(BottomKatVnode(cR-vi+1  )) ) & 
                      depthU(cU)=min(-ZW(BottomKatVnode(cR+1)),-ZW(BottomKatVnode(cR-vi+1)))
           elseif( i.eq.ui)then 
                if(  depthU(cU)>-ZW(BottomKatVnode(cR)) .or.   &
                   depthU(cU)>-ZW(BottomKatVnode(cR-vi  )))    &
                      depthU(cU)=min(-ZW(BottomKatVnode(cR)),-ZW(BottomKatVnode(cR-vi))) 
           endif
        endif
      enddo
    enddo
    ENDIF
! ------------------------------------------------------------------------------------------ 
    IF(Zgrid_depthinterp)then                                                                
    !write(*,*) 'set depthU,V along walls'                                                     ! --U(cU+ui-1)--R(cR+vi)---U(cU+ui) 
    cR = 0                                                                                    !                 |                
    cU = 0                                                                                    !                 |                
    do j=1,uj                                                                                 !   E(cR-1)      V(cR)      E(cR)  
      do i=1,vi                                                                               !                 |                
      cR = cR + 1                                                                             !                 |                
      cU=cU+1  ! add the end of the loop we substract 1 to cU when i.eq.vi                    ! --U(cU-1)-----R(cR)-------U(cU)--
      If(j<uj)then
      !write(*,*)' V node ',cR,' depthV=',depthV(cR),' kbot',BottomKatVnode(cR),ws           
      if(BottomKatVnode(cR).eq.ws .and. ( &                                        
       BottomKatRnode(cR).ne.ws .or. BottomKatRnode(cR+vi).ne.ws))then             
                                             depthV(cR)=-ZW(1)                                                                        
        if(i>1)then
         if(BottomKatUnode(cU-1).ne.ws)      depthV(cR)=min(depthU(cU-1   ),depthV(cR))  
        endif
        if(i<vi)then
           if( BottomKatUnode(cU   ).ne.ws)  depthV(cR)=min(depthU(cU   ),depthV(cR))    
        endif
        if(BottomKatRnode(cR   ).ne.ws)      depthV(cR)=min(depthR(cR   ),depthV(cR))    
        if(j<uj)then
          if( BottomKatUnode(cU+ui-1).ne.ws) depthV(cR)=min(depthU(cU+ui-1),depthV(cR))  
          if(BottomKatRnode(cR+vi).ne.ws)    depthV(cR)=min(depthR(cR+vi),depthV(cR))    
          if(i<vi)then
             if(BottomKatUnode(cU+ui).ne.ws) depthV(cR)=min(depthU(cU+ui),depthV(cR))    
          endif
       endif 

        do k = us_tridim-1,1,-1                                          
          if(eleform(i,j,k).ne.eleform(i,j,us_tridim))then 
            if(depthV(cR)>-ZW(k+1))then
             depthV(cR)= min(depthV(cR),-ZW(k+1))
            endif
            exit
          endif
        enddo
        if(i>1)then
        do k = us_tridim-1,1,-1                                          
          if(eleform(i-1,j,k).ne.eleform(i-1,j,us_tridim))then 
            if(depthV(cR)>-ZW(k+1))then
             depthV(cR)= min(depthV(cR),-ZW(k+1))
            endif
            exit
          endif
        enddo
        endif
        !!!!!!!!!!!!!!!!!!!
      endif
      Endif !j<vj
      if(i.eq.vi) then                                                                        
            cU=cU-1                                                                           
            cycle                                                                             
      endif                                                                                   
      if(BottomKatUnode(cU).eq.ws .and. ( &                                                   
       BottomKatRnode(cR).ne.ws .or. BottomKatRnode(cR+1).ne.ws))then                               ! V(cR)      E(cR)     V(cR+1)    
        depthU(cU)=-ZW(1)                                                                           ! |                        |
        if(j>1)then                                                                                 ! |                        |
         if( BottomKatVnode(cR-vi  ).ne.ws) depthU(cU)=min(depthV(cR-vi  ),depthU(cU))              ! R(cR)------U(cU)-----R(cR+1)
         if( BottomKatVnode(cR-vi+1).ne.ws) depthU(cU)=min(depthV(cR-vi+1),depthU(cU))              ! |                        |      
        endif                                                                                       ! |                        |      
        if(j<uj)then                                                                                ! V(cR-vi)  E(cR-vi)  V(cR-vi+1)  
           if( BottomKatVnode(cR  ).ne.ws)    depthU(cU)=min(depthV(cR   ),depthU(cU))       
           if( BottomKatVnode(cR+1).ne.ws)    depthU(cU)=min(depthV(cR+1 ),depthU(cU))       
        endif
        if(BottomKatRnode(cR  ).ne.ws)              depthU(cU)=min(depthR(cR   ),depthU(cU))      
        if(BottomKatRnode(cR+1).ne.ws)              depthU(cU)=min(depthR(cR+1 ),depthU(cU))       
        !write(*,*)' ele i,j=',i,j,' Udown  node brought to ',depthU(cU)                                      
        !!!!!!!!!!!!!!!!!!!                                                                                   
        if(i<vi .and. j<uj)then                                                                               
          do k = us_tridim-1,1,-1                                                                               
            if(eleform(i,j,k).ne.eleform(i,j,us_tridim))then                                                    
              if(depthU(cU)>-ZW(k+1))then                                                                       
              depthU(cU)=min(depthU(cU),-ZW(k+1))
             endif
             exit
            endif
          enddo
        endif
        if(i<vi .and. j>1)then
          do k = us_tridim-1,1,-1                                          
            if(eleform(i,j-1,k).ne.eleform(i,j-1,us_tridim))then 
              if(depthU(cU)>-ZW(k+1))then
              depthU(cU)=min(depthU(cU),-ZW(k+1))
             endif
             exit
            endif
          enddo
        endif
        !!!!!!!!!!!!!!!!!!!
      endif
     enddo
    enddo
    ENDIF
    !------------------------------------------------------

    !------------------------------------------------------
! 00 ! 01 ! 02 ! 03 ! 04 ! 05 ! 06 ! 07 ! 08 ! 09 ! 10 ! 11 ! 12 ! 13 ! 14 ! 15 !  16  !  17  !  18  !  19 
! .. ! x. ! .. ! x. ! .x ! xx !    ! xx ! .. !    ! .. ! x. ! .x ! xx ! .x ! xx ! \.\x ! ./x/ ! \x\. !  x/./ 
! .. ! .. ! x. ! x. ! .. ! .. !    ! x. ! .x !    ! xx ! xx ! .x ! .x ! xx ! xx ! x\.\ ! /x/. ! .\x\ !  /./x
    !write(*,*)'set depthE'                             ! R3---------Ub-----------R4
    cR = 0                                              ! |                       |
    cU = 0                                              ! |                       |
    depthE(:)=-999999.                                  ! Va         E(cR)     Vb
    do j=1,uj                                           ! |                       |
     do i=1,vi                                          ! |                       |
      cR = cR + 1   !move to next node number           ! R1---------Ua-----------R2
      if(i<vi .and. j<uj)then                          
       cU=cU+1                       
       dptR1=depth(cR)          
       dptR2=depth(cR+1)        
       dptR3=depth(cR+vi)       
       dptR4=depth(cR+vi+1)     
       dptUa=depthU(cU)         
       dptUb=depthU(cU+ui)      
       dptVa=depthV(cR)         
       dptVb=depthV(cR+1)
       do k = us_tridim,1,-1                                          
        if(  eleform(i,j,k).eq.3 .or. eleform(i,j,k).eq.12)then 
         depthE(cU)=max(dptUa,dptUb)
        elseif(eleform(i,j,k).eq.5 .or. eleform(i,j,k).eq.10)then  
         depthE(cU)=max(dptVa,dptVb)
        elseif(eleform(i,j,k).eq.7 )then                 
         depthE(cU)=max(dptUa,dptVb)
        elseif(eleform(i,j,k).eq.11)then                
         depthE(cU)=max(dptUb,dptVb)
        elseif(eleform(i,j,k).eq.13)then                      
         depthE(cU)=max(dptUa,dptVa)
        elseif(eleform(i,j,k).eq.14)then                       
         depthE(cU)=max(dptUb,dptVa)
        elseif(eleform(i,j,k).eq.15.or.eleform(i,j,k).eq.17   &
                            .or.eleform(i,j,k).eq.18)then
         depthE(cU)=max(max(dptUa,dptVa), max(dptUb,dptVb))
        elseif(eleform(i,j,k).eq.16)then                      
         depthE(cU)=min(dptR2,dptR3)  !0.5*(dptR2+dptR3)
        elseif(eleform(i,j,k).eq.19)then                    
         depthE(cU)=min(dptR1,dptR4)  !0.5*(dptR1+dptR4)
        endif
        if(depthE(cU)> -999998.) exit
       enddo
       if (k.eq.0)depthE(cU)=max(max(dptUa,dptVa), max(dptUb,dptVb))
       if(depthE(cU).eq.-999999.)then
        write(*,*)k,'error i,j=',i,j,' depthE=',depthE(cU),' elek='
        do k = us_tridim,1,-1
         write(*,*)eleform(i,j,k)
        enddo
        write(*,*)' ' 
       endif

      endif
     enddo
    enddo

    IF(Zgrid_depthinterp)then

    !write(*,*) 'set depthE reinterpolation, ni=',vi,ui,', nj=',uj,vj
    !write(*,'(a,61f7.1)')'Zw=',(Zw(i),i=1,ws)
    cR = 0
    cU = 0
    do j=1,uj-1
     do i=1,vi
      cR = cR + 1   !move to next node number
      if(( j.eq.uj) .or. (i.eq.vi)) cycle       
      cU=cU+1
      if(depthE(cU).eq.-999999.)then
        write(*,*)'error i,j=',i,j,' depthE=',depthE(cU)    
      endif
      dptUa=depthU(cU)                                                ! R3---------Ub-----------R4
      dptUb=depthU(cU+ui)                                             ! |                       |
      dptVa=depthV(cR)                                                ! |                       |
      dptVb=depthV(cR+1)                                              ! Va         E(cR)        Vb
      if ((dptUa<depthE(cU) .and. dptVb<depthE(cU)).or.  &            ! |                       |
          (dptVb<depthE(cU) .and. dptUb<depthE(cU)).or.  &            ! |                       |
          (dptUb<depthE(cU) .and. dptVa<depthE(cU)).or.  &            ! R1---------Ua-----------R2
          (dptVa<depthE(cU) .and. dptUa<depthE(cU))) then
          !depthE(cU)=min(depthE(cU),0.5*(0.5*(dptUa+dptUb)+0.5*(dptVa+dptVb)))
          depthE(cU)=min(depthE(cU),max(0.5*(dptUa+dptUb),0.5*(dptVa+dptVb)))
      endif
     !if( (dptUa.ge.depthE(cU).and.dptUb.ge.depthE(cU) &
     !                       .and.dptVb.ge.depthE(cU)).or. & !Va
     !    (dptVa.ge.depthE(cU).and.dptUb.ge.depthE(cU) &
     !                       .and.dptVb.ge.depthE(cU)).or. & !Ua
     !    (dptUa.ge.depthE(cU).and.dptVa.ge.depthE(cU) &
     !                       .and.dptVb.ge.depthE(cU)).or. & !Ub
     !    (dptUa.ge.depthE(cU).and.dptUb.ge.depthE(cU) &
     !                       .and.dptVa.ge.depthE(cU)) )then   
     !    factor(:)=1.0
     !    if(dptUa.ge.dptVa.and.dptUa.ge.dptVb.and.dptUa.ge.dptUb)then
     !         factor(1)=0.0
     !    elseif(dptUb.ge.dptVa.and.dptUb.ge.dptVb.and.dptUb.ge.dptUa)then
     !         factor(2)=0.0
     !    elseif(dptVa.ge.dptUa.and.dptVa.ge.dptVb.and.dptVa.ge.dptUb)then
     !         factor(3)=0.0
     !    else
     !         factor(4)=0.0 ! Vb
     !    endif
     !    if(dptVb.le.dptUa.and.dptVb.le.dptUb.and.dptVb.le.dptVa)then
     !         factor(4)=0.0
     !    elseif(dptVa.le.dptUa.and.dptVa.le.dptVb.and.dptVa.le.dptUb)then
     !         factor(3)=0.0
     !    elseif(dptUb.le.dptVa.and.dptUb.le.dptVb.and.dptUb.le.dptUa)then
     !         factor(2)=0.0
     !    else
     !         factor(1)=0.0 ! Ua
     !    endif
     !    if (factor(1)+factor(2)+factor(3)+factor(4).ne.2.0)then
     !         write(*,*)factor(1),factor(2),factor(3),factor(4)
     !         STOP 'ERROR factor in setnodesdepth'
     !    endif
     !    depthE(cU)=min(depthE(cU),&
     !    0.5*(factor(1)*dptUa+factor(2)*dptUb+factor(3)*dptVa+factor(4)*dptVb))
     !endif    
     enddo
    enddo

    !write(*,*) 'set depthE along walls'                                 ! R(cR+vi)---U(cU+ui)-----R(cR+vi+1)
    cR = 0                                                               ! |                       |
    cU = 0                                                               ! |                       |
    do j=1,uj-1                                                          ! V(cR)      E(cR)        V(cR+1)
     do i=1,vi                                                           ! |                       |
      cR = cR + 1   !move to next node number                            ! |                       |
      if(( j.eq.uj) .or. (i.eq.vi)) cycle                                ! R(cR)------U(cU)--------R(cR+1)
      cU=cU+1
      if( & !abs(depthE(cU))<0.001 .and.                   &
         (BottomKatUnode(cU   )+BottomKatVnode(cR)+    &
          BottomKatUnode(cU+ui)+BottomKatVnode(cR+1))<4*ws)then
          if(BottomKatVnode(cR   ).eq.ws.and.             &
             BottomKatUnode(cU+ui).eq.ws.and.             &
             BottomKatVnode(cR+1) .ne.ws.and.             &
             BottomKatUnode(cU   ).eq.ws) depthE(cU)=depthV(cR+1)
          if(BottomKatUnode(cU+ui).eq.ws.and.             &
             BottomKatVnode(cR+1) .eq.ws.and.             &
             BottomKatUnode(cU   ).ne.ws.and.             &
             BottomKatVnode(cR   ).eq.ws) depthE(cU)=depthU(cU)
          if(BottomKatVnode(cR+1) .eq.ws.and.             &
             BottomKatUnode(cU   ).eq.ws.and.             &
             BottomKatVnode(cR   ).ne.ws.and.             &
             BottomKatUnode(cU+ui).eq.ws) depthE(cU)=depthV(cR)
          if(BottomKatUnode(cU   ).eq.ws.and.             &
             BottomKatVnode(cR   ).eq.ws.and.             &
             BottomKatUnode(cU+ui).ne.ws.and.             &
             BottomKatVnode(cR+1).eq.ws)  depthE(cU)=depthU(cU+ui)
      endif
      if(j.eq.1.and.(i.ne.1 .and. j.ne.vi-1))then
        depthE(cU)=min(depthU(cU+ui),0.5*(depthV(cR)+depthV(cR+1)))
      elseif(j.eq.uj-1 .and.(i.ne.1 .and. j.ne.vi-1))then
        depthE(cU)=min(depthU(cU),0.5*(depthV(cR)+depthV(cR+1)))
      elseif(i.eq.1 .and.(j.ne.1 .and.j.ne.uj-1))then
        depthE(cU)=min(depthV(cR+1),0.5*(depthU(cU)+depthU(cU+ui)))
      elseif(i.eq.vi-1.and.(j.ne.1 .and.j.ne.uj-1))then
        depthE(cU)=min(depthV(cR),0.5*(depthU(cU)+depthU(cU+ui)))
      endif  
      if( (depth(cR).eq.-ZW(ws) .and.depth(cR+vi+1).eq.-ZW(ws))&
      .or.(depth(cR+1).eq.-ZW(ws) .and.depth(cR+vi).eq.-ZW(ws)))then
         depthE(cU)=-ZW(ws)
      endif
     enddo
    enddo


    !write(*,*) 'set depthR reinterpolation'
    cR = 0
    cU = 0
    do j=1,uj-1
     do i=1,vi
      cR = cR + 1   !move to next node number
      if(i<=ui)cU=cU+1
      if((j.eq.1 .or. j.eq.uj) .or. (i.eq.1 .or. i.eq.vi)) cycle      
      If(depth(cR).ne.-ZW(ws) )then
        dptUa=depthU(cU-1)                                              ! E3         Vb           E4
        dptUb=depthU(cU )                                               !            |             
        dptVa=depthV(cR-vi)                                             !            |             
        dptVb=depthV(cR   )                                             ! Ua---------R(cR)--------Ub
        if ((dptUa>depth(cR) .and. dptVb>depth(cR)).or. &               !            |             
            (dptVb>depth(cR) .and. dptUb>depth(cR)).or. &               !            |             
            (dptUb>depth(cR) .and. dptVa>depth(cR)).or. &               ! E1         Va           E2
            (dptVa>depth(cR) .and. dptUa>depth(cR)).or. &
            (dptUa<depth(cR) .and. dptVb<depth(cR)).or. &
            (dptVb<depth(cR) .and. dptUb<depth(cR)).or. &
            (dptUb<depth(cR) .and. dptVa<depth(cR)).or. &
            (dptVa<depth(cR) .and. dptUa<depth(cR)))then
            !depthR(cR)=min(depthR(cR),max(0.5*(dptUa+dptUb),0.5*(dptVa+dptVb)))
            depthR(cR)=max(0.5*(dptUa+dptUb),0.5*(dptVa+dptVb))
            !depthR(cR)=min(depthR(cR),0.5*(0.5*(dptUa+dptUb)+0.5*(dptVa+dptVb)))
            !write(*,*)i,j,'depthR',depthR(cR),depth(cR), &
            !dptUa,dptUb,dptVa,dptVb
            !if(depthR(cR)>depth(cR))stop 'ERROR depth R'
        endif
      Endif
      !if( (dptUa.gt.depth(cR).and.dptUb.gt.depth(cR) &
      !                          .and.dptVb.gt.depth(cR)).or. & !Va
      !    (dptVa.gt.depth(cR).and.dptUb.gt.depth(cR) &
      !                          .and.dptVb.gt.depth(cR)).or. & !Ua
      !    (dptUa.gt.depth(cR).and.dptVa.gt.depth(cR) &
      !                          .and.dptVb.gt.depth(cR)).or. & !Ub
      !    (dptUa.gt.depth(cR).and.dptUb.gt.depth(cR) &
      !                          .and.dptVa.gt.depth(cR)) )then   
      !    factor(:)=1.0
      !    if(dptUa.ge.dptVa.and.dptUa.ge.dptVb.and.dptUa.ge.dptUb)then
      !         factor(1)=0.0
      !    elseif(dptUb.ge.dptVa.and.dptUb.ge.dptVb.and.dptUb.ge.dptUa)then
      !         factor(2)=0.0
      !    elseif(dptVa.ge.dptUa.and.dptVa.ge.dptVb.and.dptVa.ge.dptUb)then
      !         factor(3)=0.0
      !    else
      !         factor(4)=0.0 ! Vb
      !    endif
      !    if(dptVb.le.dptUa.and.dptVb.le.dptUb.and.dptVb.le.dptVa)then
      !         factor(4)=0.0
      !    elseif(dptVa.le.dptUa.and.dptVa.le.dptVb.and.dptVa.le.dptUb)then
      !         factor(3)=0.0
      !    elseif(dptUb.le.dptVa.and.dptUb.le.dptVb.and.dptUb.le.dptUa)then
      !         factor(2)=0.0
      !    else
      !         factor(1)=0.0 ! Ua
      !    endif
      !    if (factor(1)+factor(2)+factor(3)+factor(4).ne.2.0)then
      !         STOP 'ERROR factor in setnodesdepth'
      !    endif
      !    depth(cR)=min(depth(cR),&
      !    0.5*(factor(1)*dptUa+factor(2)*dptUb+factor(3)*dptVa+factor(4)*dptVb))
      !endif    
     enddo
    enddo
    ENDIF

    call setbndeleform()
    if(BndOut)call createNodesDepthNetCDF()
   write(*,*)'---------------------------------------------------------'
  END SUBROUTINE setnodesdepth


   !*********************************************************
   !*                     setbounds                         *
   !*********************************************************
  SUBROUTINE setbounds(nbmax,us_tmp,nbounds_tmp,bndx1tmp,bndx2tmp,bndy1tmp,bndy2tmp)
   USE PARAM_MOD, ONLY: ui,uj,vi,vj,us,Zgrid,BndOut
   IMPLICIT NONE
   INTEGER,INTENT(IN):: nbmax,us_tmp
   INTEGER, INTENT(IN):: nbounds_tmp(us_tmp)
   DOUBLE PRECISION, INTENT(IN):: bndx1tmp(nbmax,us_tmp), &
                                  bndx2tmp(nbmax,us_tmp), &
                                  bndy1tmp(nbmax,us_tmp), &
                                  bndy2tmp(nbmax,us_tmp)
   INTEGER :: ff,i
     write(*,*)'setting bounds in hydrodynamic module'
    if(BndOut) write(*,*)'nbmax,us_tmp=',nbmax,us_tmp,nbounds_tmp
     if(Zgrid.and.(us_tmp.ne.us_tridim))then
         write(*,*)'us_tmp=',us_tmp,' != us_tridim=',us_tridim
         write(*,*)'ERROR DIM us_tmp sent by boundary mod'
         stop
     endif
     if(BndOut)write(*,*)'allocation nbounds(',us_tridim,')'
     allocate(nbounds(us_tridim))
 

      if(BndOut)write(*,*)'allocation bnd_x(2,',nbmax,',',us_tridim,')'
      allocate(bnd_x(2,nbmax,us_tridim))
     !write(*,*)'allocation bnd_x(2,2,1)'
     !allocate(bnd_x(2,2,1))

     if(BndOut)write(*,*)'allocation bnd_y(2,',nbmax,',',us_tridim,')'
     allocate(bnd_y(2,nbmax,us_tridim))
      if(BndOut)write(*,*)'allocation ok' 
      nbounds(:)=nbounds_tmp(:)
      if(BndOut)write(*,*)'nbounds ok'
      bnd_x(1,:,:)=bndx1tmp(:,:)
      bnd_x(2,:,:)=bndx2tmp(:,:)
      if(BndOut)write(*,*)'bnd_x ok'
      if(BndOut)write(*,*)'bnd_x ok'
      bnd_y(1,:,:)=bndy1tmp(:,:)
      bnd_y(2,:,:)=bndy2tmp(:,:)
      if(BndOut)write(*,*)'bnd_y ok'
     write(*,*)'bounds copied in hydrodynamic module'

  END SUBROUTINE setbounds



  SUBROUTINE printpython_adjele(n,k,Xpos,Ypos,eletype,flag)
   USE PARAM_MOD, ONLY: ui,uj,vi,vj,us,OutDir,NCOutFile
   USE CONVERT_MOD, ONLY: x2lon,y2lat                             
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: Xpos,Ypos
    INTEGER,INTENT(IN):: n,k,flag
    CHARACTER(LEN=1), INTENT(IN):: eletype
    INTEGER:: i,oP_ele,checkele,ff
    DOUBLE PRECISION :: Xpl,Ypl,Xpl2,Ypl2
    DOUBLE PRECISION,  DIMENSION(5)  :: Xplot,Yplot
     ff=fpy
     open (unit=fpy,file=TRIM(OutDir)//'/'//TRIM(NCOutFile)//'PartinEle.py',   &
                 POSITION='APPEND')

    write(ff,*)'#--------------------------------------'
    write(ff,'(a)')'fig,ax=plt.subplots()'
    write(ff,'(3a,i4,a)')'ax.set_title("Particle out of ',&
       trim(eletype),' elements, flag=',flag,'")'
    write(ff,'(2(a,f20.9),3a)')'im = ax.scatter(', &
          x2lon(Xpos,Ypos),',',y2lat(Ypos),        &
          ", c='r', s=65,marker='o',edgecolor='k')"

    oP_ele = P_r_element(n)
    do i=1,10
       if(r_Adjacent(oP_ele,i,k).EQ.0) exit
       checkele = r_Adjacent(oP_ele,i,k)
         Xplot(1)=x2lon(r_kwele_x(1,checkele,k),r_kwele_y(1,checkele,k))
         Xplot(2)=x2lon(r_kwele_x(2,checkele,k),r_kwele_y(2,checkele,k))
         Xplot(3)=x2lon(r_kwele_x(3,checkele,k),r_kwele_y(3,checkele,k))
         Xplot(4)=x2lon(r_kwele_x(4,checkele,k),r_kwele_y(4,checkele,k))
         Xplot(5)=x2lon(r_kwele_x(1,checkele,k),r_kwele_y(1,checkele,k))
         Yplot(1)=y2lat(                        r_kwele_y(1,checkele,k))
         Yplot(2)=y2lat(                        r_kwele_y(2,checkele,k))
         Yplot(3)=y2lat(                        r_kwele_y(3,checkele,k))
         Yplot(4)=y2lat(                        r_kwele_y(4,checkele,k))
         Yplot(5)=y2lat(                        r_kwele_y(1,checkele,k))
         write(*,'(5(a,f12.9),a)')'Xplot=[', &
              Xplot(1),',',Xplot(2),',',Xplot(3),',',Xplot(4),',',Xplot(5),']'
         write(*,'(5(a,f12.9),a)')'Yplot=[', &
              Yplot(1),',',Yplot(2),',',Yplot(3),',',Yplot(4),',',Yplot(5),']'
          write(ff,'(a)')&
        "im = ax.plot(Xplot,Yplot,linewidth=1,color='k')"
         write(ff,'(2a)')"im = ax.scatter(Xplot[:-1],Yplot[:-1],",&
            "c='w', s=50,marker='s',edgecolor='k')"
    enddo  

 
    oP_ele = P_u_element(n)
    do i=1,10
       if(u_Adjacent(oP_ele,i,k).EQ.0) exit
       checkele = u_Adjacent(oP_ele,i,k)
         Xplot(1)=x2lon(u_kwele_x(1,checkele,k),u_kwele_y(1,checkele,k))
         Xplot(2)=x2lon(u_kwele_x(2,checkele,k),u_kwele_y(2,checkele,k))
         Xplot(3)=x2lon(u_kwele_x(3,checkele,k),u_kwele_y(3,checkele,k))
         Xplot(4)=x2lon(u_kwele_x(4,checkele,k),u_kwele_y(4,checkele,k))
         Xplot(5)=x2lon(u_kwele_x(1,checkele,k),u_kwele_y(1,checkele,k))
         Yplot(1)=y2lat(                        u_kwele_y(1,checkele,k))
         Yplot(2)=y2lat(                        u_kwele_y(2,checkele,k))
         Yplot(3)=y2lat(                        u_kwele_y(3,checkele,k))
         Yplot(4)=y2lat(                        u_kwele_y(4,checkele,k))
         Yplot(5)=y2lat(                        u_kwele_y(1,checkele,k))
         write(*,'(5(a,f12.9),a)')'Xplot=[', &
              Xplot(1),',',Xplot(2),',',Xplot(3),',',Xplot(4),',',Xplot(5),']'
         write(*,'(5(a,f12.9),a)')'Yplot=[', &
              Yplot(1),',',Yplot(2),',',Yplot(3),',',Yplot(4),',',Yplot(5),']'
          write(ff,'(a)')&
      "im = ax.plot(Xplot,Yplot,linestyle='-.',linewidth=1,color='b',alpha=0.5)"
         write(ff,'(2a)')"im = ax.scatter(Xplot[:-1],Yplot[:-1],",&
            "c='w', s=50,marker='>',alpha=0.5,edgecolor='b')"
    enddo
 
    oP_ele = P_v_element(n)
    do i=1,10
       if(v_Adjacent(oP_ele,i,k).EQ.0) exit
       checkele = v_Adjacent(oP_ele,i,k)
         Xplot(1)=x2lon(v_kwele_x(1,checkele,k),v_kwele_y(1,checkele,k))
         Xplot(2)=x2lon(v_kwele_x(2,checkele,k),v_kwele_y(2,checkele,k))
         Xplot(3)=x2lon(v_kwele_x(3,checkele,k),v_kwele_y(3,checkele,k))
         Xplot(4)=x2lon(v_kwele_x(4,checkele,k),v_kwele_y(4,checkele,k))
         Xplot(5)=x2lon(v_kwele_x(1,checkele,k),v_kwele_y(1,checkele,k))
         Yplot(1)=y2lat(                        v_kwele_y(1,checkele,k))
         Yplot(2)=y2lat(                        v_kwele_y(2,checkele,k))
         Yplot(3)=y2lat(                        v_kwele_y(3,checkele,k))
         Yplot(4)=y2lat(                        v_kwele_y(4,checkele,k))
         Yplot(5)=y2lat(                        v_kwele_y(1,checkele,k))
         write(*,'(5(a,f12.9),a)')'Xplot=[', &
              Xplot(1),',',Xplot(2),',',Xplot(3),',',Xplot(4),',',Xplot(5),']'
         write(*,'(5(a,f12.9),a)')'Yplot=[', &
              Yplot(1),',',Yplot(2),',',Yplot(3),',',Yplot(4),',',Yplot(5),']'
          write(ff,'(a)')&
      "im = ax.plot(Xplot,Yplot,linestyle='-.',linewidth=1,color='g',alpha=0.5)"
         write(ff,'(2a)')"im = ax.scatter(Xplot[:-1],Yplot[:-1],",&
            "c='w', s=50,marker='^',alpha=0.5,edgecolor='g')"
    enddo
   call printpython_bounds(k,'r')
   if(k+1<us_tridim)   call printpython_bounds(k+1,'c')
   if(k>1)   call printpython_bounds(k-1,'y')

    write(ff,'(a)')"plt.show()"
    write(ff,*)'#--------------------------------------'
  close(fpy)
  END SUBROUTINE printpython_adjele

  SUBROUTINE printpython_bounds(k,col)
    IMPLICIT NONE
    INTEGER,INTENT(IN):: k
    CHARACTER(LEN=1), INTENT(IN):: col
    INTEGER:: i,ff
    ff=fpy
    do i=1,nbounds(k)
          write(ff,'(4(a,f20.9),3a)')&
               "im = ax.plot([",bnd_x(1,i,k),&
                ",",bnd_x(2,i,k), &
                "],[",bnd_y(1,i,k), &
                ",",bnd_y(2,i,k), &
             "],linewidth=2,color='",col,"',alpha=0.2)"
    enddo
  END SUBROUTINE printpython_bounds

  SUBROUTINE set_filename(var_id,counter,filename)
   USE PARAM_MOD, ONLY: numdigits,dirin, &
        prefix_Zeta,prefix_Salt,prefix_Temp,prefix_Uvel,prefix_Vvel,prefix_Wvel, & 
        prefix_Aks,prefix_Dens,prefix_Uwind,prefix_Vwind,prefix_Iwind,suffix
   IMPLICIT NONE
#include "VAR_IDs.h"
   integer, intent(in):: var_id,counter
   CHARACTER(len=200) :: prefix_var
   CHARACTER(len=200),intent(inout) :: filename

        SELECT CASE(var_id)
          CASE(VAR_ID_zeta ) 
                             prefix_var=prefix_Zeta 
          CASE(VAR_ID_salt ) 
                             prefix_var=prefix_Salt 
          CASE(VAR_ID_temp ) 
                             prefix_var=prefix_Temp 
          CASE(VAR_ID_den  ) 
                             prefix_var=prefix_Dens  
          CASE(VAR_ID_uvel ) 
                             prefix_var=prefix_Uvel 
          CASE(VAR_ID_vvel ) 
                             prefix_var=prefix_Vvel 
          CASE(VAR_ID_wvel ) 
                             prefix_var=prefix_Wvel 
          CASE(VAR_ID_kh   ) 
                             prefix_var=prefix_Aks   
          CASE(VAR_ID_uwind) 
                             prefix_var=prefix_Uwind
          CASE(VAR_ID_vwind) 
                             prefix_var=prefix_Vwind
          CASE(VAR_ID_iwind) 
                             prefix_var=prefix_Iwind
          CASE DEFAULT
           WRITE(*,*)'Model presently does not support var id ',var_id
           STOP
        END SELECT

        SELECT CASE(numdigits)
          CASE(0)
            WRITE(filename,'(A,A,A)')      TRIM(dirin),TRIM(prefix_var),      &
                                         TRIM(suffix)
          CASE(1)
            WRITE(filename,'(A,A,I1.1,A)') TRIM(dirin),TRIM(prefix_var),      &
                                         counter,TRIM(suffix)
          CASE(2)
            WRITE(filename,'(A,A,I2.2,A)') TRIM(dirin),TRIM(prefix_var),      &
                                         counter,TRIM(suffix)
          CASE(3)
            WRITE(filename,'(A,A,I3.3,A)') TRIM(dirin),TRIM(prefix_var),      &
                                         counter,TRIM(suffix)
          CASE(4)
            WRITE(filename,'(A,A,I4.4,A)') TRIM(dirin),TRIM(prefix_var),      &
                                         counter,TRIM(suffix)
          CASE(5)
            WRITE(filename,'(A,A,I5.5,A)') TRIM(dirin),TRIM(prefix_var),      &
                                         counter,TRIM(suffix)
          CASE(6)
            WRITE(filename,'(A,A,I6.6,A)') TRIM(dirin),TRIM(prefix_var),      &
                                         counter,TRIM(suffix)
          CASE(7)
            WRITE(filename,'(A,A,I7.7,A)') TRIM(dirin),TRIM(prefix_var),      &
                                         counter,TRIM(suffix)
          CASE(8)
            WRITE(filename,'(A,A,I8.8,A)') TRIM(dirin),TRIM(prefix_var),      &
                                         counter,TRIM(suffix)
          CASE(9)
            WRITE(filename,'(A,A,I9.9,A)') TRIM(dirin),TRIM(prefix_var),      &
                                         counter,TRIM(suffix)
          CASE(10)
            WRITE(filename,'(A,A,I10.10,A)') TRIM(dirin),TRIM(prefix_var),    &
                                         counter,TRIM(suffix)
          CASE DEFAULT
           WRITE(*,*)'Model presently does not support numdigits of ',numdigits
           WRITE(*,*)'Please use numdigit value from 1 to 10'
           WRITE(*,*)'  OR modify code in Hydrodynamic module'
           STOP
        END SELECT
  END SUBROUTINE

  CHARACTER(len=200) FUNCTION roms_netcdf_var_name(var_id)
#include "VAR_IDs.h"
   IMPLICIT NONE
   integer, intent(in):: var_id

        SELECT CASE(var_id)
          CASE(VAR_ID_zeta )
                             roms_netcdf_var_name='zeta' 
          CASE(VAR_ID_salt ) 
                             roms_netcdf_var_name='salt' 
          CASE(VAR_ID_temp ) 
                             roms_netcdf_var_name='temp' 
          CASE(VAR_ID_den  ) 
                             roms_netcdf_var_name='rho'  
          CASE(VAR_ID_uvel ) 
                             roms_netcdf_var_name='u' 
          CASE(VAR_ID_vvel ) 
                             roms_netcdf_var_name='v' 
          CASE(VAR_ID_wvel ) 
                             roms_netcdf_var_name='w' 
          CASE(VAR_ID_kh   ) 
                             roms_netcdf_var_name='AKs'   
          CASE(VAR_ID_uwind) 
                             roms_netcdf_var_name='sustr'
          CASE(VAR_ID_vwind) 
                             roms_netcdf_var_name='svstr'
          CASE(VAR_ID_iwind) 
                             roms_netcdf_var_name='wind_intensity'
          CASE DEFAULT
           WRITE(*,*)'Model presently does not support var id ',var_id
           STOP
        END SELECT
  END FUNCTION

  CHARACTER(len=200) FUNCTION explicit_var_name(var_id)
#include "VAR_IDs.h"
   IMPLICIT NONE
   integer, intent(in):: var_id

        SELECT CASE(var_id)
          CASE(VAR_ID_zeta ) 
                             explicit_var_name='zeta' 
          CASE(VAR_ID_salt ) 
                             explicit_var_name='salt' 
          CASE(VAR_ID_temp ) 
                             explicit_var_name='temperature' 
          CASE(VAR_ID_den  ) 
                             explicit_var_name='density'  
          CASE(VAR_ID_uvel ) 
                             explicit_var_name='u_velocity' 
          CASE(VAR_ID_vvel ) 
                             explicit_var_name='v_velocity' 
          CASE(VAR_ID_wvel ) 
                             explicit_var_name='w_velocity' 
          CASE(VAR_ID_kh   ) 
                             explicit_var_name='AKs'   
          CASE(VAR_ID_uwind) 
                             explicit_var_name='u_wind'
          CASE(VAR_ID_vwind) 
                             explicit_var_name='v_wind'
          CASE(VAR_ID_iwind) 
                             explicit_var_name='wind_intensity'
          CASE DEFAULT
           WRITE(*,*)'Model presently does not support var id ',var_id
           STOP
        END SELECT
  END FUNCTION

  SUBROUTINE read_data_from_file(var_id,ni,nj,nk,nt,tarray,tf1,tff,  &
                                 field,RUVnod,recordnum,incrstepf,interpolate)
   USE PARAM_MOD, ONLY: ui,uj,vi,vj,us,ws,Zgrid,hydrobytes,Zgrid,filenum
   USE RANDOM_MOD, ONLY: genrand_real1
   USE netcdf
   IMPLICIT NONE
   INCLUDE 'netcdf.inc'
   integer,intent(in):: var_id
   integer, intent(in):: ni,nj,nk,nt,tarray,tf1,tff
   double precision, intent(inout):: field(ni,nj,nk,nt)
   integer,intent(in):: RUVnod,recordnum,incrstepf
   integer,intent(in):: interpolate
   character(200) :: filenm
   integer, allocatable, dimension(:):: start_index,count_index
   integer:: nk_MITfile,ios,waiting
   integer:: interpol_uv,one_if_interpol_v,rand15,t,k,j,i,ktlev
   REAL, ALLOCATABLE, DIMENSION(:) :: vec_prev
   REAL :: real_vec_read(vi)
   DOUBLE PRECISION :: dbl_vec_read(vi)
   INTEGER:: STATUS,NCID,VID
   CHARACTER(len=200) :: varname
  
   call set_filename(var_id,iint+filenum,filenm)

   if(Zgrid .or. nk>1)then
     allocate(start_index(4))
     allocate(count_index(4))
   else
     allocate(start_index(3))
     allocate(count_index(3))
   endif
   start_index(1)=t_ijruv(IMIN,RUVnod)
   start_index(2)=t_ijruv(JMIN,RUVnod)
   count_index(1)=t_ijruv(IMAX,RUVnod)-t_ijruv(IMIN,RUVnod)+1
   count_index(2)=t_ijruv(JMAX,RUVnod)-t_ijruv(JMIN,RUVnod)+1
   if(Zgrid)then
     nk_MITfile=min(us,nk)! here using us even when nk=uw=us+1 as W output has dim us=uw-1 instead of uw for MITGCM (no bottom value)
     start_index(3)=1
     start_index(4)=recordnum
     count_index(3)=nk_MITfile
     count_index(4)=incrstepf
     if(interpolate>0)then
        interpol_uv=RUVnod
        write(*,'(4a,3(a,i8))')'read in MITgcm file var ',trim(explicit_var_name(var_id)),'  at cell center, interpolating it on the cell borders from file',TRIM(filenm), &
          ' for time record num=',start_index(4),':',start_index(4)+count_index(4)-1,' nk=',nk

     else
        interpol_uv=0
        write(*,'(4a,3(a,i8))')'read in MITgcm file var ',trim(explicit_var_name(var_id)),' from file',TRIM(filenm), &
          ' for time record num=',start_index(4),':',start_index(4)+count_index(4)-1,' nk=',nk
     endif
   else
     interpol_uv=0
     if(nk==1)then
       start_index(3)=recordnum
       count_index(3)=incrstepf
     else
       start_index(3)=1
       start_index(4)=recordnum
       count_index(3)=nk
       count_index(4)=incrstepf
     endif
     varname=trim(roms_netcdf_var_name(var_id))
     write(*,'(4a,3(a,i8))')'read in ROMs file var ',trim(varname),' from file',TRIM(filenm),' for time record num=',recordnum,':',recordnum+incrstepf-1,' nk=',nk
   endif
     if(interpol_uv==VNODE)then
       ALLOCATE(vec_prev(vi))
       one_if_interpol_v=1 
     else
       one_if_interpol_v=0 
     endif

     if(Zgrid)then
       !write(*,*)'time record num=',start_index(4)
       open (unit=110,file=TRIM(filenm),form='unformatted',status='old',   & 
             action='read',access='direct', recl=hydrobytes*vi, iostat=ios,convert='little_endian')          !--- CL-OGS: all var are read with dim vi !
        if ( ios /= 0 ) then
            do waiting=1,10
                    rand15=int( 15.0*genrand_real1() )
                    call sleep(rand15)
                    write(*,*)'waiting as opening of file ',trim(filenm),  &
                              ' failed . New trial after sleep ',rand15
                    open(unit=110,file=TRIM(filenm),form='unformatted',    &
                      status='old', action='read',access='direct',         &
                      recl=hydrobytes*vi, iostat=ios,convert='little_endian')
                    if ( ios == 0 ) exit
            enddo
        endif
       if ( ios /= 0 )  then
           write(*,*) " ERROR OPENING ",TRIM(filenm)
           stop
       endif

        !write(*,*) ' i=',start_index(1),':',start_index(1)+count_index(1)-1,' and (1:',ni,')=vec(',vi-ni+1,':',vi,')'
        !write(*,*) ' j=',start_index(2),':',start_index(2)+count_index(2)-1+one_if_interpol_v
        !write(*,*) ' k=',1,':',nk_MITfile
        !write(*,*) ' t=',start_index(3),':',start_index(3)+count_index(3)-1

       do t=start_index(4),start_index(4)+count_index(4)-1
         real_vec_read=0.0
         dbl_vec_read=0.0
         field(:,:,:,tarray+t-start_index(4))=0.0
         do k=start_index(3),start_index(3)+count_index(3)-1
           ktlev=(t-1)*(nk_MITfile*uj)+(k-1)*uj
           do j=start_index(2),start_index(2)+count_index(2)-1+one_if_interpol_v
             if(hydrobytes.eq.4)then 
               read(110,rec=(ktlev+j+(uj-nj-one_if_interpol_v)),IOSTAT=ios)real_vec_read(1:vi) ! adding +(uj-nj)=1 for v nodes as mitgcm output is on uj nodes = vj+1
               if ( ios == 0 ) then
                 if(interpol_uv==UNODE) real_vec_read(2:vi)=0.5*(real_vec_read(1:vi-1)+real_vec_read(2:vi)) 
                 ! adding +(nk-nk_MITfile)=+1 to skip bottom values for Wtype nodes that have nk-us=1
                 if(interpol_uv==VNODE.and.j.gt.start_index(2))then
                    field(1:ni,j-1,(nk-k+1),tarray+t-start_index(4))=             &
                               0.5*(real_vec_read(vi-ni+1:vi)+vec_prev(vi-ni+1:vi))
                 else
                    field(1:ni,j,(nk-k+1),tarray+t-start_index(4))=real_vec_read(vi-ni+1:vi)
                 endif
                 if(interpol_uv==VNODE) vec_prev=real_vec_read 
               endif
             else
               read(110,rec=(ktlev+j+(uj-nj-one_if_interpol_v)),IOSTAT=ios)dbl_vec_read(1:vi) ! adding +(uj-nj)=1 for v nodes as mitgcm output is on uj nodes = vj+1
               if ( ios == 0 ) then
                 if(interpol_uv==UNODE)dbl_vec_read(2:vi)=0.5*(dbl_vec_read(1:vi-1)+dbl_vec_read(2:vi))
                 ! for nk=ws, nk-k+1=us-k+2 : skipping bottom values for Wtype nodes that have nk-us=ws-us=1
                 if(interpol_uv==VNODE.and.j.gt.start_index(2))then
                    field(1:ni,j-1,(nk-k+1),tarray+t-start_index(4))=             &
                               0.5*(dbl_vec_read(vi-ni+1:vi)+vec_prev(vi-ni+1:vi))
                 else
                    field(1:ni,j,(nk-k+1),tarray+t-start_index(4))=dbl_vec_read(vi-ni+1:vi)
                 endif
                 if(interpol_uv==VNODE) vec_prev=dbl_vec_read 
               endif
             endif
             if ( ios /= 0 ) then
                write(*,*) 'Problem reading  ',varname
                write(*,'(4(a,i4),2(a,i8),6(a,i4))')'i=',vi-ni+1,':',vi,' t=',t,' k=',k, &
                             ' ktlev=',ktlev,' rec=',ktlev+j+(uj-nj-one_if_interpol_v), &
                ' j=',j,' -> i=',1,':',ni,' j=',j,' k=',(nk-k+1),' t=',tarray+t-start_index(4)
                write(*,*) ' i=',start_index(1),':',start_index(1)+count_index(1)-1
                write(*,*) ' j=',start_index(2),':',start_index(2)+count_index(2)-1
                write(*,*) ' k=',start_index(3),':',start_index(3)+count_index(3)-1
                write(*,*) ' t=',start_index(4),':',start_index(4)+count_index(4)-1
                stop " ERROR READING FIELD "//varname
             endif
           enddo
         enddo
       enddo
       if(interpol_uv==VNODE)then
         DEALLOCATE(vec_prev)
       endif

       CLOSE(110)

     else ! Roms NETcdf outputs
       STATUS = NF90_OPEN(TRIM(filenm), NF90_NOWRITE, NCID)
       if (STATUS .NE. NF90_NOERR) write(*,*) 'Problem NF90_OPEN'
       if (STATUS .NE. NF90_NOERR) write(*,*) NF90_STRERROR(STATUS)
      
       STATUS = NF90_INQ_VARID(NCID,varname,VID)
       if (STATUS .NE. NF90_NOERR) then
         write(*,*) 'Problem finding ',varname
         write(*,*) NF90_STRERROR(STATUS)
         stop
       endif
         STATUS = NF90_GET_VAR(NCID,VID,field(t_ijruv(IMIN,RUVnod):t_ijruv(IMAX,RUVnod), &
                         t_ijruv(JMIN,RUVnod):t_ijruv(JMAX,RUVnod),1:nk,tf1:tff),     &
                         start_index,count_index )
       if (STATUS .NE. NF90_NOERR) then
         write(*,*) 'Problem reading ',varname
         write(*,*) ' i=',start_index(1),':',start_index(1)+count_index(1)-1
         write(*,*) ' j=',start_index(2),':',start_index(2)+count_index(2)-1
         write(*,*) ' k=',1,':',nk
         write(*,*) ' t=',recordnum,':',recordnum+incrstepf-1
         write(*,*) NF90_STRERROR(STATUS)
         stop
       endif
       if(var_id==VAR_ID_uwind)then
          do i=t_ijruv(IMIN,RUVnod),t_ijruv(IMAX,RUVnod)
           do j=t_ijruv(JMIN,RUVnod),t_ijruv(JMAX,RUVnod)
             DO t=tf1,tff
             if(field(i,j,1,t).lt.0.0)then
             field(i,j,1,t) = (-20.659 * (abs(field(i,j,1,t))**0.4278))
             else
             field(i,j,1,t)= (20.659 * (field(i,j,1,t)**0.4278)) 
             end if
             ENDDO
           enddo
          enddo
       elseif(var_id==VAR_ID_vwind)then
          do i=t_ijruv(IMIN,RUVnod),t_ijruv(IMAX,RUVnod)
           do j=t_ijruv(JMIN,RUVnod),t_ijruv(JMAX,RUVnod)
             DO t=tf1,tff
             if(field(i,j,1,t).lt.0.0)then
               field(i,j,1,t)= (-20.659 * (abs(field(i,j,1,t))**0.4278)) 
             else
               field(i,j,1,t)= (20.659 * (field(i,j,1,t)**0.4278))
             end if
             ENDDO
           enddo
          enddo
       endif
       !close the dataset and reassign the NCID
       STATUS = NF90_CLOSE(NCID)
     endif
     deallocate(start_index,count_index)

  END SUBROUTINE
END MODULE HYDRO_MOD
