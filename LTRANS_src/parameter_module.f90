MODULE PARAM_MOD 
 
!  The Parameter Module reads in the include file, LTRANS.data, making the 
!  parameters declared within available to all the other modules. It also 
!  reads in information from the NetCDF grid file and calculates values of
!  grid specific parameters, making them available to all the other modules.
! 
!  Created by:            Zachary Schlag 
!  Created on:            28 Jul 2008 
!  Last Modified on:         Feb 2013
!--- CL-OGS:   Modified on: June 2016 by Celia Laurent - OGS
!--- CL-OGS:   main modifications are:
!--- CL-OGS:   - coupling of LTRANS v2b to the Z-grid and fields of the MITgcm
!--- CL-OGS:   - allow running LTRANS using a data file other than LTRANS.data
!--- CL-OGS:   - reading of 3d masks instead of 2d
!--- CL-OGS:   - not computing wetR[U,V] neither rho[u,v]_elements 
!--- CL-OGS:     as this will be made later in hydrodynamic_module computing 
!--- CL-OGS:     rho_kwele(k), u_kwele(k) and v_kwele(k) at each vertical level k
!--- CL-OGS:   - implementation of the Irish Marine Institute Oil Model Module in LTRANS v2bZ-MITgcm
 
IMPLICIT NONE 
PUBLIC 
SAVE 

  include 'LTRANS.h'

CONTAINS


  SUBROUTINE getParams(inputdatafile)
  !Subroutine to read all input parameters from LTRANS.data 

    character(len=120) :: header
    integer :: istat,err,count
    !--- CL-OGS: modified to allow running LTRANS using the data file which name is given in argument
    !--- CL-OGS: use  "$ LTRANS.exe LTRANS_inputdatafile.data"
    CHARACTER(len=100), INTENT(in) :: inputdatafile     
    write(*,*)'Reading parameters from file ',trim(inputdatafile)
    err = 0
    prefix_Salt  = '' 
    prefix_Temp  = '' 
    prefix_Uvel  = '' 
    prefix_Vvel  = '' 
    prefix_Wvel  = '' 
    prefix_Aks   = '' 
    prefix_Dens  = '' 
    prefix_Uwind = '' 
    prefix_Vwind = '' 
    prefix_Iwind = '' 
    Adjele_fname='Adjacentelements.data'
    ADJele_file= .FALSE.
    habitatfile='NONE'
    holefile='NONE'
    swdown_ASCIIfname='NONE'
    SeabedRelease=.FALSE.
    SeabedRelease_meters=0.0
    Stokes=.False.
    OPEN(1,file=trim(inputdatafile))                  !--- read control variables:
      IF(err == 0) THEN
        READ(1,nml=numparticles ,IOSTAT=istat)  !--- number of particles
        IF(istat/=0)err = 10
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=timeparam    ,IOSTAT=istat)  !--- time info
        !write(*,timeparam)
        IF(istat/=0)err = 20
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=hydroparam   ,IOSTAT=istat)  !--- hydrodynamics info
        !write(*,hydroparam)
        IF(istat/=0)err = 30
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=turbparam    ,IOSTAT=istat)  !--- turbulence info
        !write(*,turbparam)
        IF(istat/=0)err = 40
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=behavparam   ,IOSTAT=istat)  !--- behavior info
        !write(*,behavparam)
        IF(istat/=0)err = 50
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=dvmparam     ,IOSTAT=istat)  !--- diurnal vertical migration 
        !write(*,dvmparam)
        IF(istat/=0)err = 60
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=settleparam  ,IOSTAT=istat)  !--- settlement info
        !write(*,settleparam)
        IF(istat/=0)err = 70
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=convparam    ,IOSTAT=istat)  !--- unit conversion
        !write(*,convparam)
        IF(istat/=0)err = 80
      ENDIF
      !--- CL-OGS: name of model "roms"->"hydromodel" modified as this version of LTRANS runs as well with the MITgcm
      IF(err == 0) THEN
        READ(1,nml=hydromodelgrid     ,IOSTAT=istat)  !--- hydromodel grid
        !write(*,hydromodelgrid)
        IF(istat/=0)err = 90
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=hydromodeloutput   ,IOSTAT=istat)  !--- hydromodel history output file
        !write(*,hydromodeloutput)
        IF(istat/=0)err = 100
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=parloc       ,IOSTAT=istat)  !--- particle locations
        !write(*,parloc)
        IF(istat/=0)err = 110
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=HabPolyLoc   ,IOSTAT=istat)  !--- habitat polygon info
        !write(*,HabPolyLoc)
        IF(istat/=0)err = 120
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=output       ,IOSTAT=istat)  !--- output related info
        !write(*,output)
        IF(istat/=0)err = 130
      ENDIF
      IF(err == 0) THEN
        READ(1,nml=other        ,IOSTAT=istat)  !--- other misc 
        !write(*,other)
        IF(istat/=0)err = 140
      ENDIF

!      *****   IMIOM      *****
        IF(OilOn)THEN
        IF(err == 0)THEN
          READ(1,nml=oilparams ,IOSTAT=istat)   !--- Oil Model Properties
          !write(*,oilparams)
          IF(istat/=0)err = 170
        ENDIF
        IF(err == 0)THEN
          READ(1,nml=oilprocs   ,IOSTAT=istat)  !--- Oil Processes
          !write(*,oilprocs)
           IF(istat/=0)err = 180
        ENDIF
        IF(err == 0)THEN
          READ(1,nml=windswaves ,IOSTAT=istat)  !--- SWAN Model Properties
          !write(*,windswaves)
          IF(istat/=0)err = 190
        ENDIF
      END IF
!      ***** END IMIOM *****

    CLOSE(1)

    IF(err == 0) THEN
      call gridData(IOSTAT=istat)
      IF(istat/=0)err = 150
    ENDIF

    SELECT CASE(err)
      CASE(0)
        header='No Errors'
      CASE(10)
        header='Error when reading numparticles, pls check LTRANS.data'
      CASE(20)
        header='Error when reading timeparam, pls check LTRANS.data'
      CASE(30)
        header='Error when reading hydroparam, pls check LTRANS.data'
      CASE(40)
        header='Error when reading turbparam, pls check LTRANS.data'
      CASE(50)
        header='Error when reading behavparam, pls check LTRANS.data'
      CASE(60)
        header='Error when reading behavdvm, pls check LTRANS.data'
      CASE(70)
        header='Error when reading settleparam, pls check LTRANS.data'
      CASE(80)
        header='Error when reading convparam, pls check LTRANS.data'
      CASE(90)
        header='Error when reading hydromodelgrid, pls check LTRANS.data'
      CASE(100)
        header='Error when reading hydromodeloutput, pls check LTRANS.data'
      CASE(110)
        header='Error when reading parloc, pls check LTRANS.data'
      CASE(120)
        header='Error when reading HabPolLoc, pls check LTRANS.data'
      CASE(130)
        header='Error when reading output, pls check LTRANS.data'
      CASE(140)
        header='Error when reading other, pls check LTRANS.data'
      CASE(150)
        header='Error when reading gridinfo, pls check GRID.data'
      CASE(180)
        header='Error when reading oil process, pls check LTRANS.data'
      CASE(190)
        header='Error when reading wave properties, pls check LTRANS.data'
      CASE DEFAULT
        header='Error: unexpected err number'
    END SELECT

    IF(err/=0) then
        if(oilOn)then
                write(*,*)header,err
                err = (err/10 + 1000)
!                  open(999,file='temporary',status='rewind')
                  open(999,file='temporary',POSITION='rewind')
                  write(999,*)err !was write(999,*)err       on 13/01/2012
                  close(999)
                CALL errorHandler(Header,-1)        !print the error message and stop
        else
                CALL errorHandler(Header,-1)        !print the error message and stop
        end if
    ENDIF

    !Subract 1 from lonmin and latmin to eliminate roundoff error
    lonmin = lonmin - 1
    latmin = latmin - 1

  END SUBROUTINE getParams


  SUBROUTINE errorHandler(header, flag)
    IMPLICIT NONE
    CHARACTER(LEN=120), INTENT(IN) :: header
    INTEGER, INTENT(IN) :: flag
   
    IF (flag .eq. -1) THEN
      WRITE(*,"(A120)")header               !print error message in report.txt
      STOP
    ELSE
      WRITE(*,"('***** WARNING *****')")    !print warning message to screen
      WRITE(*,"(A120)")header
    ENDIF
   
  END SUBROUTINE errorHandler


  SUBROUTINE gridData(IOSTAT)
    ! This subroutine was originally a separate program: Grid_Generator.f90
    ! Created by:           Zachary Schlag
    ! Program Created:      28 Aug 2008
    ! Subroutine Created:   12 Nov 2010
    ! Last Modified on:        Feb 2011
    USE netcdf
    IMPLICIT NONE

    INTEGER, INTENT(OUT), OPTIONAL :: IOSTAT

    INCLUDE 'netcdf.inc'

    !NetCDF Variables
    INTEGER :: STATUS,GF_ID,VID,dimid,dimcount
    INTEGER :: xi_rho,xi_u,xi_v,eta_rho,eta_u,eta_v
    INTEGER :: s_rho,s_w    !--- CL-OGS: number of third-dimension vertical levels of rho and w grids
    !--- CL-OGS: extending masks to third dimension and adding mask_w
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: mask_rho,mask_u,mask_v,mask_w

    !Grid File Output Variables
    INTEGER :: nR,nU,nV,maxR,maxU,maxV    !--- CL-OGS: deleted wetR,wetU,wetV 

    !Iteration Variables
    INTEGER :: err

    err = 0

    ! *********************** GET GRID INFO ***********************

    ! OPEN NETCDF FILE - GET GF_ID VALUE

    STATUS = NF90_OPEN(NCgridfile,NF90_NOWRITE,GF_ID)
    if (STATUS .NE. NF90_NOERR) then
      write(*,*) 'Problem NF90_OPEN Error opening',NCgridfile
      err = 10
    endif

    ! GET VALUES FOR xi_rho,xi_u,xi_v,eta_rho,eta_u,eta_v,s_rho

      STATUS = NF90_INQ_DIMID(GF_ID,'xi_rho',dimid)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem dimid xi_rho in INQ_DIMID'
      endif
      STATUS = NF90_INQUIRE_DIMENSION(GF_ID,dimid,len=dimcount)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem dimid xi_rho'
        err = 20 
      endif
      xi_rho = dimcount

      STATUS = NF90_INQ_DIMID(GF_ID,'eta_rho',dimid)
      STATUS = NF90_INQUIRE_DIMENSION(GF_ID,dimid,len=dimcount)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem dimid eta_rho'
        err = 20 
      endif
      eta_rho = dimcount

      STATUS = NF90_INQ_DIMID(GF_ID,'xi_u',dimid)
      STATUS = NF90_INQUIRE_DIMENSION(GF_ID,dimid,len=dimcount)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem dimid xi_u'
        err = 20 
      endif
      xi_u = dimcount

      STATUS = NF90_INQ_DIMID(GF_ID,'eta_u',dimid)
      STATUS = NF90_INQUIRE_DIMENSION(GF_ID,dimid,len=dimcount)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem dimid eta_u'
        err = 20 
      endif
      eta_u = dimcount

      STATUS = NF90_INQ_DIMID(GF_ID,'xi_v',dimid)
      STATUS = NF90_INQUIRE_DIMENSION(GF_ID,dimid,len=dimcount)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem dimid xi_v'
        err = 20 
      endif
      xi_v = dimcount

      STATUS = NF90_INQ_DIMID(GF_ID,'eta_v',dimid)
      STATUS = NF90_INQUIRE_DIMENSION(GF_ID,dimid,len=dimcount)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem dimid eta_v'
        err = 20 
      endif
      eta_v = dimcount

     !--- CL-OGS: the rho,u,v grids have the third dimension s_rho.
     !--- CL-OGS: the w grid has dimensions xi_rho, eta_rho, and s_w
     !--- CL-OGS: so there is no need to read xi_w and eta_w
     if(Zgrid)then
       STATUS = NF90_INQ_DIMID(GF_ID,'Z',dimid)
       STATUS = NF90_INQUIRE_DIMENSION(GF_ID,dimid,len=dimcount)
       if (STATUS .NE. NF90_NOERR) then
         write(*,*) 'Problem dimid Z'
         err = 20 
       endif
       s_rho = dimcount

       STATUS = NF90_INQ_DIMID(GF_ID,'Zi',dimid)
       STATUS = NF90_INQUIRE_DIMENSION(GF_ID,dimid,len=dimcount)
       if (STATUS .NE. NF90_NOERR) then
         write(*,*) 'Problem dimid Zi'
         err = 20 
       endif
       s_w = dimcount
     else
       s_rho=1
       s_w=1 
     endif
 
    ! ALLOCATE VARIABLE ARRAY DIMENSIONS

      ALLOCATE (mask_rho(xi_rho,eta_rho,s_rho),STAT=STATUS) !--- CL-OGS: added third dim s_rho
      if(STATUS /= 0) then
        write(*,*) 'Problem allocating mask_rho'
        err = 30 
      endif
      ALLOCATE (mask_u(xi_u,eta_u,s_rho),STAT=STATUS)       !--- CL-OGS: added third dim s_rho
      if(STATUS /= 0) then
        write(*,*) 'Problem allocating mask_u'
        err = 30 
      endif
      ALLOCATE (mask_v(xi_v,eta_v,s_rho),STAT=STATUS)       !--- CL-OGS: added third dim s_rho
      if(STATUS /= 0) then
        write(*,*) 'Problem allocating mask_v'
        err = 30 
      endif
      !--- CL-OGS: added allocation of mask_w
      ALLOCATE (mask_w(xi_rho,eta_rho,s_w),STAT=STATUS)
      if(STATUS /= 0) then
        write(*,*) 'Problem allocating mask_w'
        err = 30 
      endif
      write(*,*)' xi_rho =', xi_rho
      write(*,*)'eta_rho =',eta_rho
      write(*,*)' xi_u   =', xi_u
      write(*,*)'eta_u   =',eta_u
      write(*,*)' xi_v   =', xi_v
      write(*,*)'eta_v   =',eta_v
      write(*,*)'s_rho   =',  s_rho
      write(*,*)'s_w     =',  s_w
      if (Zgrid .and. (us .ne. s_rho)) then ! rho-grid dimension in z direction
        write(*,*)'Problem s_rho != us'
        us=s_rho
      endif      
      if (Zgrid .and. (ws .ne. s_w)) then ! rho-grid dimension in z direction
        write(*,*)'Problem s_w != ws'
        ws=s_w
      endif    
        

      ! rho grid mask
      STATUS = NF90_INQ_VARID(GF_ID,'mask_rho',VID)
      STATUS = NF90_GET_VAR(GF_ID,VID,mask_rho)
      if (STATUS .NE. NF90_NOERR) then
        write(*,*) 'Problem read mask_rho'
        err = 40 
      endif

      ! u grid mask
      STATUS = NF90_INQ_VARID(GF_ID,'mask_u',VID)
      STATUS = NF90_GET_VAR(GF_ID,VID,mask_u)
      if(STATUS .NE. NF90_NOERR) then
        write(*,*)'Problem read mask_u'
        err = 40 
      endif

      ! v grid mask
      STATUS = NF90_INQ_VARID(GF_ID,'mask_v',VID)
      STATUS = NF90_GET_VAR(GF_ID,VID,mask_v)
      if(STATUS .NE. NF90_NOERR) then
        write(*,*)'Problem read mask_v'
        err = 40 
      endif

     !--- CL-OGS: added read mask_w
     !! w grid mask
     !STATUS = NF90_INQ_VARID(GF_ID,'mask_w',VID)
     !STATUS = NF90_GET_VAR(GF_ID,VID,mask_w)
     !if (STATUS .NE. NF90_NOERR) then
     !  write(*,*) 'Problem read mask_w'
     !  err = 40 
     !endif


    STATUS = NF90_CLOSE(GF_ID)
    if(STATUS /= NF90_NOERR) then
      write(*,*)'Problem closing GF_ID'
      err = 50
    endif


  ! ********************** MAKE GRID FILE **********************
    !Calculate number of nodes in each grid
    nR = xi_rho * eta_rho
    nU = xi_u * eta_u
    nV = xi_v * eta_v

    !Calculate Maximum number of elements in each grid
    maxR = (xi_rho-1)*(eta_rho-1)
    maxU = (xi_u-1)*(eta_u-1)
    maxV = (xi_v-1)*(eta_v-1)

    !--- CL-OGS:  deleted the computation of the wet elements

    !Set Parameter Values
    ui = xi_u               ! u-grid dimension in x direction
    uj = eta_u              ! u-grid dimension in y direction
    vi = xi_v               ! v-grid dimension in x direction
    vj = eta_v              ! v-grid dimension in y direction
    rho_nodes = nR          ! number of rho points (vi*uj)
    u_nodes = nU            ! number of u points   (ui*uj)
    v_nodes = nV            ! number of v points   (vi*vj)
    max_rho_elements = maxR ! Number of elements using four rho points as corners
    max_u_elements = maxU   ! Number of elements using four  u  points as corners
    max_v_elements = maxV   ! Number of elements using four  v  points as corners
    !--- CL-OGS: deleted assignation of wetR,wetU and wetV to rho_elements, u_elements and v_elements
    write(*,*)'ui=',ui
    write(*,*)'uj=',uj
    write(*,*)'vi=',vi
    write(*,*)'vj=',vj
    write(*,*)'rho_nodes=',rho_nodes
    write(*,*)'u_nodes=',u_nodes
    write(*,*)'v_nodes=',v_nodes
    write(*,*)'max_rho_elements =',max_rho_elements 
    write(*,*)'max_u_elements   =',max_u_elements   
    write(*,*)'max_v_elements   =',max_v_elements   
    !--- CL-OGS: added deallocate(mask_rho,mask_u,mask_v,mask_w)
    DEALLOCATE(mask_rho,mask_u,mask_v,mask_w) 


    !If IOSTAT is present, set return value to error code
    IF(PRESENT(IOSTAT)) IOSTAT = err
    !  0=No Errors                 30=Error allocating arrays
    ! 10=Error Opening NCgridfile  40=Error getting variables
    ! 20=Error getting dimensions  50=Error Closing NCgridfile

  END SUBROUTINE gridData

END MODULE PARAM_MOD 
