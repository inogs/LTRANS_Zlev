!
! LTRANS-Zlev version 0(beta)
!
! **********************************************************************
! **********************************************************************
! **      Original LTRANS v.2b source code : Copyright (c) 2013       **
! **   The University of Maryland Center for Environmental Science .  **
! **********************************************************************
! **         OILTRANS oil module : Copyright (c) 2016                 **
! **               The Marine Institute, Ireland                      **
! **********************************************************************
! **    Other code modified and added in the LTRANS-Zlev version :    **
! **                    Copyright (c) 2019 OGS                        **
! **   Istituto Nazionale di Oceanografia e Geofisica Sperimentale    ** 
! **(National Institute of Oceanography and Applied Geophysics, Italy)**
! **********************************************************************
! **                                                                  **
! ** This Software is open-source and licensed under the following    **
! ** conditions as stated by MIT/X License:                           **
! **                                                                  **
! **  (See http://www.opensource.org/licenses/mit-license.php ).      **
! **                                                                  **
! ** Permission is hereby granted, free of charge, to any person      **
! ** obtaining a copy of this Software and associated documentation   **
! ** files (the "Software"), to deal in the Software without          **
! ** restriction, including without limitation the rights to use,     **
! ** copy, modify, merge, publish, distribute, sublicense,            **
! ** and/or sell copies of the Software, and to permit persons        **
! ** to whom the Software is furnished to do so, subject to the       **
! ** following conditions:                                            **
! **                                                                  **
! ** The above copyright notice and this permission notice shall      **
! ** be included in all copies or substantial portions of the         **
! ** Software.                                                        **
! **                                                                  **
! ** THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,  **
! ** EXPRESSED OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE           **
! ** WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE  **
! ** AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT  **
! ** HOLDERS BE LIABLE FOR ANY CLAIMS, DAMAGES OR OTHER LIABILITIES,  **
! ** WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     **
! ** FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR    **
! ** OTHER DEALINGS IN THE SOFTWARE.                                  **
! **                                                                  **
! **********************************************************************
! ** The most current official versions of the LTRANS-Zlev software   **
! ** and associated tools and documentation are available at:         **
! **                                                                  **
! **         http://github.com/inogs/LTRANS_Zlev                      **
! **                                                                  **
! **********************************************************************
! ** The most current official versions of the original v.2b software **
! ** and associated tools and documentation are available at:         **
! **                                                                  **
! **  http://northweb.hpl.umces.edu/LTRANS.htm                        **
! **                                                                  **
! **********************************************************************
! ** The most current official versions of the OILTRANS oil module    **
! ** and associated tools and documentation are available from the    ** 
! ** authors by e-mail:                                               **
! **                                                                  **
! **      ocean.modelling@marine.ie                                   **
! **                                                                  **
! **********************************************************************
! ** We ask that users make appropriate acknowledgement of            **
! ** The University of Maryland Center for Environmental Science,     **
! ** individual developers, participating agencies and institutions,  **
! ** and funding agencies. One way to do this is to cite one or       **
! ** more of the relevant publications listed at:                     **
! **                                                                  **
! **  http://northweb.hpl.umces.edu/LTRANS.htm#Description            **
! **********************************************************************
! ** We ask that users make appropriate acknowledgement of            **
! ** The Marine Institute, Ireland,                                   **
! ** individual developers, participating agencies and institutions,  **
! ** and funding agencies. One way to do this is to cite              **
! ** the following publication:                                       **
! **                                                                  **
! ** Berry, A., Dabrowski, T., Lyons, K., 2012. The oil spill model   **
! ** OILTRANS and its application to the Celtic Sea. Marine Pollution **
! ** Bulletin, 64(11) : 2489-2501.                                    **
! **********************************************************************
! ** We ask that users make appropriate acknowledgement of the OGS,   **
! ** Istituto Nazionale di Oceanografia e Geofisica Sperimentale      **
! **(National Institute of Oceanography and Applied Geophysics, Italy)**
! ** individual developers, participating agencies and institutions,  **
! ** and funding agencies. One way to do this is to cite one or       **
! ** more of the relevant publications  :                             **
! **                                                                  **
! **  Laurent, C. , Querin, S. Solidoro, C., Melaku Canu, D.,         **
! **" Modelling marine particle dynamics with LTRANS-Zlev:            **
! **  implementation and validation" (in preparation)                 **
! **********************************************************************
! ********************************************************************** 
! As the original v.2b version, the LTRANS-Zlev includes the Mersenne 
! Twister random number generator and the tension spline
!  curve-fitting package (TSPACK) which have their own
! licenses and are freely available for non-commercial uses. 


PROGRAM main

! LTRANS.f90 contains the main structure of the particle-tracking program. 
! It executes the external time step, internal time step, and particle loops, 
! advects particles, and writes output. It calls modules that read in 
! hydrodynamic model information, move particles due to turbulence and 
! behavior, test if particles are in habitat polygons, and apply boundary 
! conditions to keep particles in the model domain. 
!
! Developers:
!   Elizabeth North: enorth@umces.edu
!   Zachary Schlag: zschlag@umces.edu
!   Ian Mitchell: imitchell@umces.edu
!   Celia Laurent: claurent@inogs.it - OGS (Zlev version)

IMPLICIT NONE
#define ID_TEMP 1
#define ID_DEPTH 2
#define ID_U_WIND 3
#define ID_V_WIND 4
#define NUM_ID_VALUES 4
#define ID_STRANDED 5
#define ID_ALIVE 6
#define NUM_ID_NUMPARTS 6
!   *************************************************************************
!   *                                                                       *
!   *                       Variable Declarations                           *
!   *                                                                       *
!   *************************************************************************

  INTEGER, PARAMETER :: nAttrib   = 13

  INTEGER, PARAMETER :: pX        =  1  ! Particle X-coordinate
  INTEGER, PARAMETER :: pY        =  2  ! Particle Y-coordinate
  INTEGER, PARAMETER :: pZ        =  3  ! Particle Z-coordinate
  INTEGER, PARAMETER :: pnX       =  4  ! Particle new X-coordinate
  INTEGER, PARAMETER :: pnY       =  5  ! Particle new Y-coordinate
  INTEGER, PARAMETER :: pnZ       =  6  ! Particle new Z-coordinate
  INTEGER, PARAMETER :: ppX       =  7  ! Particle previous X-coordinate
  INTEGER, PARAMETER :: ppY       =  8  ! Particle previous Y-coordinate
  INTEGER, PARAMETER :: ppZ       =  9  ! Particle previous Z-coordinate
  INTEGER, PARAMETER :: pStatus   = 10  ! Status of particle (previously Color)
  INTEGER, PARAMETER :: pDOB      = 11  ! Particle Date Of Birth
  INTEGER, PARAMETER :: pAge      = 12  ! Particle Age (s)
  INTEGER, PARAMETER :: pLifespan = 13  ! Age at which particle settled or died

  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: par
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: P_Salt,P_Temp,P_coastdist
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: P_GrainSize,P_Size
  INTEGER, ALLOCATABLE, DIMENSION( : ) :: P_MainPoly

  INTEGER, ALLOCATABLE, DIMENSION( : ) :: P_oldLev
  INTEGER, ALLOCATABLE, DIMENSION(:) :: startpoly,endpoly,hitBottom,hitLand, &
                                        ID_Poly
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: Time_in_Poly
  INTEGER :: npoly,poly0,polyN
  double precision, allocatable, dimension(:) :: parIniDepth  !--- CL-OGS: to keep particles at constant depth under the moving sea surface

  DOUBLE PRECISION :: ex(3),ix(3)
  INTEGER :: prcount,printdt,p,it,fpy
  REAL :: timeCounts(8),times(9)

  DOUBLE PRECISION :: pTS   !pTS = Particle release Time of Spill (seconds) [for temporal/spatial spills]
  INTEGER:: Average_Numpart(NUM_ID_NUMPARTS)     
  DOUBLE PRECISION :: Average_Value(NUM_ID_VALUES)
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: PartAtSurf
!   *************************************************************************
!   *                                                                       *
!   *                             Execution                                 *
!   *                                                                       *
!   *************************************************************************

  call run_LTRANS()

contains



! MODEL IRF SUBROUTINES
  subroutine run_LTRANS()
    ! *************************************************************************
    ! *                                                                       *
    ! *                              Run Model                                *
    ! *                                                                       *
    ! *************************************************************************
    use param_mod, only: days,dt,  &
                         outpathGiven,outpath,NCOutFile, &      !--- CL-OGS: for file name management
                         OilOn,WriteCurrents,numpar            !--- CL-OGS

    integer :: seconds,stepT,n
    character(10) :: realtime                                   !--- CL-OGS: computational time tracking
    character(250)::OilPropfilename                             !--- CL-OGS: for OIL transport module
    character(LEN=40) :: fstring                                     !--- CL:OGS
    realtime='-'
    !call date_and_time(TIME=realtime)                           !--- CL-OGS: computational time tracking
    !write(*,*)' before ini_LTRANS time is ',realtime(1:2),'h',realtime(3:4),'m',realtime(5:6),'s'
  

    call ini_LTRANS()
      call date_and_time(TIME=realtime)                       !--- CL-OGS: computational time tracking
      write(*,*)' after ini_LTRANS time is ',realtime(1:2),'h',realtime(3:4),'m',realtime(5:6),'s'
      if(WriteCurrents)then
       DO n=1,numpar
        write(fstring,'(a,i5.5,a)')'Currents',n,''//TRIM(NCOutFile)//'.csv'
        OPEN(110,FILE=trim(fstring),POSITION='APPEND',status='replace')
        CLOSE(110)
       ENDDO
      endif 

      write(*,'(/,A)') '****** BEGIN ITERATIONS *******'

      ! days*24*60*60 = total number of seconds to run the model
      ! divide that by dt to get the number of external time steps
      seconds = int(days*86400.0) !Total seconds to run model
      !--- CL-OGS: added abs() on next line to allow backward-in-time simulation
      stepT   = abs(seconds/dt)        !number of external time steps
      do p=1,stepT
        call date_and_time(TIME=realtime)                  !--- CL-OGS: computational time tracking
        write(*,'(a,i6,7a)')'External timestep num',p,' time is ', &
                             realtime(1:2),'h',realtime(3:4),'m',realtime(5:6),'s'
        call run_External_Timestep()
      enddo
    call fin_LTRANS()
    call date_and_time(TIME=realtime)                      !--- CL-OGS: computational time tracking
    write(*,*)' after  fin_LTRANS time is',realtime(1:2),'h',realtime(3:4),'m',realtime(5:6),'s'
  end subroutine run_LTRANS



  subroutine ini_LTRANS()
    ! *************************************************************************
    ! *                                                                       *
    ! *                           Initialize Model                            *
    ! *                                                                       *
    ! *************************************************************************
    use behavior_mod, only: initBehave,setOut,die
    use boundary_mod, only: createBounds,mbounds,ibounds,Get_coastdist
    use convert_mod,  only: lon2x,lat2y,x2lon,y2lat                                  !--- CL-OGS   

    use random_mod,   only: init_genrand
    use hydro_mod,    only: initGrid,initHydro,setEle_all,initNetCDF,       &
                      createNetCDF,writeNetCDF,                             &
                      getKRlevel,setnodesdepth, &                 !--- CL-OGS: for coupling with MITgcm'Z-grid bathymetry and fields
                      setInterp,getInterp,getDepth
    use param_mod,    only: numpar,days,dt,idt,seed,parfile,settlementon,   &
                      Behavior,TrackCollisions,SaltTempOn,writeNC,          &
                      WriteHeaders,WriteModelTiming,ErrorFlag,getParams,    & 
                      WriteCurrents,NCOutFile,OutDir,                       & !--- CL-OGS: for file name management
                      Zgrid,                                          & !--- CL-OGS: for coupling with MITgcm'Z-grid bathymetry and fields 
                      Ext0,iprinto,mortality,                         & !--- CL-OGS: for backward-in-time simulation, restarts, ...
                      read_GrainSize,                                 & !--- CL-OGS: for behavior type 8
                      constTemp,constUwind,constVwind,                &
                      SeabedRelease,SeabedRelease_meters,             &
                      Write_coastdist,StrandingDist,us,Write_Poly_Presence, &
!        *****   IMIOM        *****
                      OilOn,WindWeatherFac
    use oil_mod, only: InitOilModel,OilModel
    use tension_mod, only : initTensionModule
    use settlement_mod, only :  get_NumPoly,get_IDPoly
#include "VAR_IDs.h"
    IMPLICIT NONE

    double precision :: pDep,P_depth
    integer:: m,STATUS
!        ***** END IMIOM *****

    integer :: n,ele_err
    double precision, allocatable, dimension(:) :: pLon,pLat

    ! Initial Boundary check
    integer :: in_island,inbounds,m_nestdeg,i_nestdeg,counterrors
    integer:: island
  !--- CL-OGS: for coupling with MITgcm'Z-grid bathymetry and fields:
  !--- CL-OGS: as does the variable "island" for the id number of the island,
  !--- CL-OGS: mainbound identifies the id number of the main water boundary basin where a particle 
  !--- CL-OGS: is located, as using Z-grid bathymetry induced the existence of multiples water basins
    integer:: mainbound
    CHARACTER(len=100) :: inputdatafile  !--- CL-OGS: store name of ".data" input file
    integer :: klev,Fstlev,conflict               !--- CL-OGS: store the vertical level containing a particle
    INTEGER :: Phase1Time
    conflict=1
 
  ! ***************************************************************************
  ! *                          Get Parameter Values                           *
  ! ***************************************************************************
  !--- CL-OGS: get name of LTRANS.data file as first argument after the name of the executable
    CALL getarg(1, inputdatafile) 
    CALL getParams(inputdatafile)
    CALL writeModelInfo()
      ex=0.0
      ex(1) = (1-2)*dt  +Ext0
      ex(2) = (1-1)*dt  +Ext0
      ex(3) = 1*dt      +Ext0
      ix(1) = dt+Ext0 + DBLE((1-2)*idt)
      ix(2) = dt+Ext0 + DBLE((1-1)*idt)
      ix(3) = dt+Ext0 + DBLE(1*idt)

    write(*,*) ' '
    write(*,*) ' *************** LTRANS INITIALIZATION ************** '

  ! ***************************************************************************
  ! *                       Allocate Dynamic Variables                        *
  ! ***************************************************************************
    ! CL-OGS -------------
    if(Behavior.eq.997)  then
      ALLOCATE(PartAtSurf(numpar,4) )
      PartAtSurf(:,:)=0.0
      mortality=.TRUE.
      write(*,*)'Mortality set to true for Behavior 997'
    endif

    ALLOCATE(par(numpar,nAttrib))

!        *****   IMIOM        *****
         IF(SettlementOn)THEN
            ALLOCATE(startpoly(numpar))
            ALLOCATE(endpoly(numpar))
            endpoly = 0               !initialize end polygon location to zero
         ENDIF
         If(.not.OilOn)then
          !Local variables for read-in of Latitude and Longitude
           ALLOCATE(pLon(numpar))
           ALLOCATE(pLat(numpar))
         Endif
!        ***** END IMIOM *****

    !--- CL-OGS: took allocations of P_Salt,P_Temp,hitBottom and hitLand out of if(.not.OilOn) paragraph
    IF(SaltTempOn)THEN
    ALLOCATE(P_Salt(numpar))
    ALLOCATE(P_Temp(numpar))
      P_Temp = 0.0
      P_Salt = 0.0
      IF(Behavior.ge.8.and.Behavior.le.10)THEN
        P_Temp=-999.0
      ELSE
        P_Temp = 0.0
      ENDIF
    ENDIF
    
    IF( StrandingDist>=0 .and. (.not.settlementon ))then
      write(*,*)'error StrandingDist>=0 .and. (.not.settlementon )'
      stop 'set settlementon to True or StrandingDist<0'
    ENDIF
    IF( StrandingDist>=0 .or. Write_coastdist )  then
      ALLOCATE(P_coastdist(numpar))
      P_coastdist= 0.0
    ENDIF
    
    if(read_GrainSize)then
      ALLOCATE(P_GrainSize(numpar))
      IF(Behavior.ge.8.and.Behavior.le.11)THEN
         P_GrainSize = -999.0
      ELSE
         P_GrainSize = 0.0
      ENDIF
    ENDIF

    IF(Behavior.GE.8 .and. Behavior.LE.11) THEN
      ALLOCATE(P_Size(numpar))
      P_Size = 0.0
    ENDIF 

    IF((.not.read_GrainSize) .and. &
        (Behavior.ge.8.and.Behavior.le.11) )THEN
     write(*,*) "WARNING, Missing GrainSize file for behavior 8"
     write(*,*) "Constant value -999.0 is taken"
    ENDIF
    ALLOCATE(P_oldLev(numpar))     
    P_oldLev(:)=-1 

    IF(TrackCollisions)THEN
      ALLOCATE(hitBottom(numpar))
      ALLOCATE(hitLand(numpar))
      hitBottom = 0
      hitLand = 0
    ENDIF

    ALLOCATE(parIniDepth(numpar))  !--- CL-OGS: store the initial depth of each particle
    ! *************************************************************************
    ! *         Initialize print counters and random number generator         *
    ! *************************************************************************

    ! THE FOLLOWING VARIABLE INITIALIZATIONS SHOULD NOT BE CHANGED:
    prcount=0                  !print counter; number of external time steps
    printdt=-abs(iprinto)      ! first instant of printing
    CALL init_genrand(seed)    !set random number generator Seed Value

    ! *************************************************************************
    ! *         Initialize the tension spline module                          *
    ! *************************************************************************

    CALL initTensionModule() !CL-OGS : added to allocate OMP_IL array for OMP simulations

  
    ! *************************************************************************
    ! *                   Initialize Particle Attributes                      *
    ! *************************************************************************

    ! Read-in lat/long of particles. If settlement module is on, read in    
    ! the habitat polygon on which the particle start                       
    write(*,*) 'read in particle locations for numpar=', numpar

    OPEN(1,FILE=TRIM(parfile))
!        write(*,*)TRIM(parfile)
!        *****   IMIOM        *****
        if(OilOn)then
!                read(1,*)nts                                                                                                                        !numbers of times for spill events
!                ALLOCATE(pLon(nts),pLat(nts),pDep(nts),pTS(nts))                                                !allocate appropriate array sizes
!                ALLOCATE(pLon(1),pLat(1))                                                                !allocate appropriate array sizes
                ALLOCATE(pLon(numpar),pLat(numpar),STAT=STATUS)                                                                !allocate appropriate array sizes

                if(STATUS /= 0) then
                  write(*,*) 'Problem allocating pLon,pLat'
                  stop
                else
                  
                  write(*,*) 'allocating pLon,pLat went right'
                endif

!                do n = 1, nts                                                                                                                        !for each spill time (DOB)
!                        read(1,*)pLon(n),pLat(n),pDep(n),pTS(n)                                                                !Lon,Lat,Dep,Time
!                        do m = ((numpar/nts)*(n-1))+1, ((numpar/nts)*n)
                        read(1,*)pLon(1),pLat(1),pDep,pTS                                                                !Lon,Lat,Dep,Time
                        write(*,*)'Oil Slick at pos ',pLon(1),pLat(1),pDep,pTS
                        do m = 1, numpar
!                        par(m,pX)  = lon2x(pLon(n),pLat(n))
!                        par(m,pY)  = lat2y(pLat(n))
!                                par(m,pZ)  = pDep(n)
                        pLon(m)=pLon(1)
                        pLat(m)=pLat(1) 
                        parIniDepth(m)= pDep
                        par(m,pX)  = lon2x( pLon(1),  pLat(1))
                        par(m,pY)  = lat2y(pLat(1))
                                par(m,pZ)  = pDep
                        par(m,pnX) = par(m,pX)
                        par(m,pnY) = par(m,pY)
                        par(m,pnZ) = par(m,pZ)
                        par(m,ppX) = par(m,pX)  !0.0                                                                                            !initialize to 0.0 to indicate no previous location
                        par(m,ppY) = par(m,pY)  !0.0                                                                                            !initialize to 0.0 to indicate no previous location
                        par(m,ppZ) = par(m,pZ)  !0.0                                                                                            !initialize to 0.0 to indicate no previous location
                        par(m,pStatus)   = 0                                                                                        !floating oil
                        par(m,pAge)      = 0.0
!                        par(m,pDOB)      = pTS(n)                                                                                !set particle DoB to Time of Spill (pTS)
                        par(m,pDOB)      = pTS                                                                                        !set particle DoB to Time of Spill (pTS)
                        par(m,pLifespan) = 0.0
                    end do
!                end do
        else
!        ***** END IMIOM *****

            do n=1,numpar
              if(settlementon .and.  StrandingDist<0)then
                read (1,*) pLon(n),pLat(n),par(n,pZ),par(n,pDOB),startpoly(n) 
              else
                read (1,*) pLon(n),pLat(n),par(n,pZ),par(n,pDOB)
              endif
              !--- CL-OGS: store in parIniDepth the initial depth for simulations with constant depth
              parIniDepth(n)=par(n,pZ)
              par(n,pX)  = lon2x(pLon(n),pLat(n))
              par(n,pY)  = lat2y(pLat(n))
              par(n,pnX) = par(n,pX)
              par(n,pnY) = par(n,pY)
              par(n,pnZ) = par(n,pZ)
              par(n,ppX) = 0.0    !initialize to 0.0 to indicate no previous location
              par(n,ppY) = 0.0    !initialize to 0.0 to indicate no previous location
              par(n,ppZ) = 0.0    !initialize to 0.0 to indicate no previous location
              par(n,pStatus)   = Behavior
              par(n,pAge)      = 0.0
              par(n,pLifespan) = 0.0
            enddo
        endif

    CLOSE(1)

   write(*,*) '  Particle n=1 Latitude=',pLat(1),'Longitude=',pLon(1)
   write(*,*) '  Particle n=1 Depth=',par(1,pZ)
   write(*,*) '  Particle n=1 X=',par(1,pX),'Y=',par(1,pY)
   write(*,*) '  Particle n=1 X2lon=',x2lon(par(1,pX),par(1,pY)),'Y2lat=',y2lat(par(1,pY))
!   if(settlementon) write(*,*) '  Particle n=5 Start Polygon=',startpoly(5)
 
    ! *************************************************************************
    ! *                                                                       *
    ! *           Initial Read-In of Hydrodynamic Model Information           *
    ! *                                                                       *
    ! *************************************************************************

    !Initialize Grid / Create Elements
    CALL initGrid()


    !Initialize Behavior
    CALL initBehave()

     IF(SettlementOn.and.Write_Poly_Presence) then
                 write(*,*)'Initialize polygon related arrays'
                 ALLOCATE(P_MainPoly(numpar))
                 do n=1,numpar
                  P_MainPoly(n)=int(startpoly(n))
                 enddo
                 call get_NumPoly(poly0,polyN,npoly)
                 ALLOCATE(Time_in_Poly(polyN-poly0+1,numpar))
                 ALLOCATE(ID_Poly(npoly))
                 call get_IDPoly(ID_Poly) 
                 write(*,*)'done'
      ENDIF

  
    ! *******************************************************************
    ! *                    Initialize NetCDF Output                     *
    ! *******************************************************************

    !Create NetCDF Output File if Needed
    IF(writeNC) then
      CALL initNetCDF()
      CALL createNetCDF(par(:,pDOB))
    ENDIF

    prcount = 1
    if(writeNC)then !Write to NetCDF Output File
      CALL writeNetCDF(Ext0,par(:,pAge),pLon,pLat,par(:,pZ),par(:,pStatus),   &
                      SALT=P_Salt,TEMP=P_Temp,GrSize=P_GrainSize,SizeP=P_Size,   &
                 HITB=hitBottom,HITL=hitLand,CoDi=P_coastdist,POLY=P_MainPoly)
     !if (SaltTempOn) then
     !  if(TrackCollisions)then
     !    !--- CL-OGS: Ext0 contains time in seconds at the beginning of the simulation
     !    !--- CL-OGS: (at the first external timestep). To be given in LTRANSinputfile.data
     !    if((Behavior.ge.8.and.Behavior.le.10))THEN
     !    CALL writeNetCDF(Ext0,par(:,pAge),pLon,pLat,par(:,pZ),par(:,pStatus),   &
     !                    SALT=P_Salt,TEMP=P_Temp,GrSize=P_GrainSize,PSize=P_Size,   &
     !                    HITB=hitBottom,HITL=hitLand)
     !    else
     !    CALL writeNetCDF(Ext0,par(:,pAge),pLon,pLat,par(:,pZ),par(:,pStatus),   &
     !                    SALT=P_Salt,TEMP=P_Temp,HITB=hitBottom,HITL=hitLand)
     !    endif
     !  else
     !    if((Behavior.ge.8.and.Behavior.le.10))THEN
     !    CALL writeNetCDF(Ext0,par(:,pAge),pLon,pLat,par(:,pZ),par(:,pStatus),   &
     !                    SALT=P_Salt,TEMP=P_Temp,GrSize=P_GrainSize,PSize=P_Size)
     !    else
     !    CALL writeNetCDF(Ext0,par(:,pAge),pLon,pLat,par(:,pZ),par(:,pStatus),   &
     !                    SALT=P_Salt,TEMP=P_Temp)
     !    endif
     !  endif
     !else
     !  if(TrackCollisions)then
     !    if((Behavior.ge.8.and.Behavior.le.10))THEN
     !    CALL writeNetCDF(Ext0,par(:,pAge),pLon,pLat,par(:,pZ),par(:,pStatus),   &
     !                    GrSize=P_GrainSize,PSize=P_Size,HITB=hitBottom,HITL=hitLand)
     !    else
     !    CALL writeNetCDF(Ext0,par(:,pAge),pLon,pLat,par(:,pZ),par(:,pStatus),   &
     !                    HITB=hitBottom,HITL=hitLand)
     !    endif
     !  else
     !    if((Behavior.ge.8.and.Behavior.le.10))THEN
     !    CALL writeNetCDF(Ext0,par(:,pAge),pLon,pLat,par(:,pZ),par(:,pStatus),   &
     !                    GrSize=P_GrainSize,PSize=P_Size)
     !    else
     !    CALL writeNetCDF(Ext0,par(:,pAge),pLon,pLat,par(:,pZ),par(:,pStatus))
     !    endif
     !  endif
     !endif
    endif

   
    ! *************************************************************************
    ! *                                                                       *
    ! *                      Prepare for Particle Tracking                    *
    ! *                                                                       *
    ! *************************************************************************

    write(*,*) 'prepare boundary arrays'

    !Create Boundaries
    CALL createBounds()
    
    if(Zgrid) CALL setnodesdepth()

    !Create file to output information if a problem is encountered
    !--- CL-OGS:  OutDir and NCOutFile strings are used to create the name of the ErrorLog file
    SELECT CASE(ErrorFlag)
    CASE(1)
      OPEN(210,FILE=TRIM(OutDir)//'/ErrorLog'//TRIM(NCOutFile)//'.txt',STATUS='REPLACE')
        write(210,*) &
          'The following particles were returned to their previous locations:'
        write(210,*) ' '
      CLOSE(210)
    CASE(2)
      OPEN(210,FILE=TRIM(OutDir)//'/ErrorLog'//TRIM(NCOutFile)//'.txt',STATUS='REPLACE')
        write(210,*) 'The following particles were killed:'
        write(210,*) ' '
      CLOSE(210)
    CASE(3)
      OPEN(210,FILE=TRIM(OutDir)//'/ErrorLog'//TRIM(NCOutFile)//'.txt',STATUS='REPLACE')
        write(210,*) 'The following particles were set out of bounds:'
        write(210,*) ' '
      CLOSE(210)
    END SELECT
    !Get Start Elements
    write(*,*) "finding each particle's initial element"
    counterrors=0
    do n=1,numpar

      inbounds = 0
      !--- CL-OGS:  for Z-Grid model get the vertical level of the particle, that will be 
      !--- CL-OGS:  needed in mbounds and ibounds routines given that for the Z-grid the
      !--- CL-OGS:  boundaries are level-dependant
      if(Zgrid)then
        klev=getKRlevel(par(n,Pz))
      else
        klev=1
      endif
      !Determine if particle is within model boundaries
      call mbounds(par(n,pY),par(n,pX),inbounds,mainbound,klev,m_nestdeg)
      if (inbounds.EQ.0) then 
        counterrors=counterrors+1
        if(ErrorFlag < 1 .OR. ErrorFlag > 3)then
          write(*,*) 'Particle initial location outside main bounds, n=',n
          write(*,*) 'x:   ',par(n,pX),' y:   ',par(n,pY)
          write(*,*) 'lon: ',pLon(n)  ,' lat: ',pLat(n)
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and is Terminating'
          stop
        else
          write(*,*) 'Particle initial location outside main bounds, n=',n
          if(ErrorFlag == 2)then
            call die(n)
          else
            call setOut(n)
          endif
          OPEN(210,FILE=TRIM(OutDir)//'/ErrorLog'//TRIM(NCOutFile)//'.txt',POSITION='APPEND')
            write(210,"('Particle ',I10,' initially outside main bounds')") n
            write(210,*) 'lon: ',pLon(n)  ,' lat: ',pLat(n),' k=',klev
            write(210,*)' '
          CLOSE(210)
          cycle
        endif
      endif

      in_island = 0
      call ibounds(in_island,par(n,pY),par(n,pX),island,klev,i_nestdeg)
      if (in_island.EQ.1) then
        counterrors=counterrors+1
          write(*,*) 'Particle ',n,' initial location is within an island ',&
              ' at lev ',klev
         ! write(*,*) 'x:   ',par(n,pX),' y:   ',par(n,pY)
         ! write(*,*) 'lon: ',pLon(n)  ,' lat: ',pLat(n)
        if(ErrorFlag < 1 .OR. ErrorFlag > 3)then
          write(*,*) ' '
          write(*,*) 'The Program Cannot Continue and is Terminating'
          stop
        else
          if(ErrorFlag == 2)then
            call die(n)
          else
            call setOut(n)
          endif
          OPEN(210,FILE=TRIM(OutDir)//'/ErrorLog'//TRIM(NCOutFile)//'.txt',POSITION='APPEND')
            write(210,"('Particle ',I10,' initially inside island bounds')") n
            write(210,*) 'lon: ',pLon(n)  ,' lat: ',pLat(n),' k=',klev
            write(210,*)' '
          CLOSE(210)
          cycle
        endif
      endif
      IF( StrandingDist>=0 .or. Write_coastdist ) then
      call Get_coastdist(par(n,pY),par(n,pX),us,P_coastdist(n))
      ENDIF

    enddo
    if(counterrors>0)then
     write(*,*)'Some particles have an initial location out of the water basin'
     write(*,*)'The Program Cannot Continue and is Terminating'
     stop
    endif
    !Determine which Rho, U, & V elements the particles are in

    !--- CL-OGS: par(:,pZ) had to be given to setEle_all as the element number became dependant 
    !--- CL-OGS: on the vertical level of the particle, recomputed in setEle_all based on par(:,pZ)
    write(*,*)'setEle_all'
    CALL setEle_all(par(:,pX),par(:,pY),par(:,pZ),ele_err,n)

    !If the particle was not found to be within an element,
    ! write a message to the screen and discontinue the program
    IF(ele_err .ne. 0)THEN
      if(ErrorFlag < 1 .OR. ErrorFlag > 3)then

        write(*,*) " "

        SELECT CASE (ele_err)
          CASE(1)
            write(*,*) n,' start - particle not inside any rho element'
          CASE(2)
            write(*,*) n,' start - particle not inside any u element'
          CASE(3)
            write(*,*) n,' start - particle not inside any v element'
        END SELECT

        write(*,*) ' '
        write(*,*) ' - Now Stopped - '
        write(*,*) ' '
        write(*,*) 'Current Location:',  '      k:   ',klev
        write(*,*) '  x:   ',par(n,pX),  ' y:   ',par(n,pY)
        write(*,*) '  lon: ',pLon(n),' lat: ',pLat(n)
        write(*,*) ' '
        write(*,*) 'The Program Cannot Continue and Will Terminate'
        stop
      else
        if(ErrorFlag == 2)then
          call die(n)
        else
          call setOut(n)
        endif
        OPEN(210,FILE=TRIM(OutDir)//'/ErrorLog'//TRIM(NCOutFile)//'.txt',POSITION='APPEND')
          SELECT CASE (ele_err)
            CASE(1)
              write(210,"('Particle ',I10,' initially not in rho element')") n
            CASE(2)
              write(210,"('Particle ',I10,' initially not in u element')") n
            CASE(3)
              write(210,"('Particle ',I10,' initially not in v element')") n
          END SELECT
            write(210,*)' '
        CLOSE(210)
      endif
    ENDIF

   if(SeabedRelease)then
    do m = 1, numpar
      if(Zgrid)then 
        klev=getKRlevel(par(m,pZ))
      else
        klev=1
      endif     
      !Set Interpolation Values for the current particle
      CALL setInterp(par(m,pX),par(m,pY),m)
      !Find depth, angle, and sea surface height at particle location
      conflict=1
      if(.not. Zgrid)then ! (ROMS sigma-level-grid :)
       P_depth = DBLE(-1.0)* getInterp(par(m,pX),par(m,pY),VAR_ID_depth,klev)
      else !             (Zgrid :)      
        call getDepth(par(m,pX),par(m,pY),m,it,P_depth,Fstlev,conflict)  
      endif
      parIniDepth(m)=P_depth+SeabedRelease_meters
      par(m,pZ) =parIniDepth(m)
      par(m,pnZ) = par(m,pZ)
      if(m.eq.1)then
      write(*,*)' Releasing particles at ',SeabedRelease_meters,' above seabed'
      write(*,*)' Example part 1: depth=',P_depth,' Zpart=',par(m,pZ)
      endif
     enddo
    endif !if(SeabedRelease)


    !Read in initial hydrodynamic model data
    CALL initHydro()

    !Create files to output 'land hits' and 'bottom hits'
    IF(TrackCollisions) then
      OPEN(100,FILE=TRIM(OutDir)//'/LandHits'//TRIM(NCOutFile)//'.csv',STATUS='REPLACE')
        write(100,*)'numpar,lon,lat,depth,age,time,hitLand'
      CLOSE(100)
      OPEN(101,FILE=TRIM(OutDir)//'/BottomHits'//TRIM(NCOutFile)//'.csv',STATUS='REPLACE')
        write(101,*)'numpar,lon,lat,depth,age,time,hitBottom'
      CLOSE(101)
    ENDIF

  
    !Create file to track model timing
    IF(WriteModelTiming)then
      OPEN(300,FILE=TRIM(OutDir)//'/Timing'//TRIM(NCOutFile)//'.csv',STATUS='REPLACE')
        write(300,"(2(a14),7(a11,a8,','))")'Daytime ,','Elapsed ,', &
        'Hydro',' Hyd%','Set Up',     &
        '  SetUp%','Advection','   Adv%','HTurb','  HTurb% ,','VTurb',     &
        ' VTurb%','Behavior',' Behave%,','Update',' Update%'
      CLOSE(300)

    ENDIF

    !Create files with Header Information
    IF(WriteHeaders)then
      OPEN(100,FILE=TRIM(OutDir)//'/para_Headers.txt',STATUS='REPLACE')
        write(100,*)'column 01: depth   -Depth of particle at end of time ',   &
                    'step (meters)'
        write(100,*)'column 02: color   -integer value to indicate the ',      &
                    'status of the particle'
        write(100,*)'column 03: lon     -Longitude of particle at end of ',    &
                    'time step (decimal °)'
        write(100,*)'column 04: lat     -Latitude  of particle at end of ',    &
                    'time step (decimal °)'
        IF(SaltTempOn)then
          write(100,*)'column 05: salt    -Salinity at particle location at ', &
                    'end of time step'
          write(100,*)'column 06: temp    -Temperature at particle location ', &
                    'at end of time step (°C)'
        ENDIF
      CLOSE(100)

      IF(TrackCollisions)then
        OPEN(100,FILE=TRIM(OutDir)//'/LandHits_Headers.txt',STATUS='REPLACE')
          write(100,*)'column 01: numpar  -Particle identification number ',   &
                      '(dimensionless)'
          write(100,*)'column 02: lon     -Longitude of particle at end of ',  &
                      'time step (decimal °)'
          write(100,*)'column 03: lat     -Latitude  of particle at end of ',  &
                      'time step (decimal °)'
          write(100,*)'column 04: depth   -Depth of particle at end of time ', &
                      'step (meters)'
          write(100,*)'column 05: age     -Age of particle (in days since ',   &
                      'released)'
          write(100,*)'column 06: time    -Model time (in days since the ',    &
                      'start of the model)'
          write(100,*)'column 07: hitLand -Number of times the particle ',     &
                      'struck land in the last print interval time step'
        CLOSE(100)

        OPEN(101,FILE=TRIM(OutDir)//'/BottomHits_Headers.txt',STATUS='REPLACE')
          write(101,*)'column 01: numpar    -Particle identification number ', &
                      '(dimensionless)'
          write(101,*)'column 02: lon       -Longitude of particle at end ',   &
                      'of time step (decimal °)'
          write(101,*)'column 03: lat       -Latitude  of particle at end ',   &
                      'of time step (decimal °)'
          write(101,*)'column 04: depth     -Depth of particle at end of ',    &
                      'time step (meters)'
          write(101,*)'column 05: age       -Age of particle (in days since ', &
                      'released)'
          write(101,*)'column 06: time      -Model time (in days since the ',  &
                      'start of the model)'
          write(101,*)'column 07: hitBottom -Number of times the particle ',   &
                      'struck bottom in the last print interval time step'
        CLOSE(101)
      ENDIF

      if(WriteCurrents)then
        OPEN(102,FILE=TRIM(OutDir)//'/Currents_Headers.txt',STATUS='REPLACE')
          write(102,*)'column 01: daytime   = t/86400' 
          write(102,*)'column 02: Xposition'
          write(102,*)'column 03: Yposition'
          write(102,*)'column 04: Zposition'
          write(102,*)'column 05: Advection currents direction X (m/s)'
          write(102,*)'column 06: Advection currents direction Y (m/s) '
          write(102,*)'column 07: Advection currents direction Z (m/s) '
          write(102,*)'column 08: Turbulent motion   direction X (m/s) '
          write(102,*)'column 09: Turbulent motion   direction Y (m/s) '
          write(102,*)'column 10: Turbulent motion   direction Z (m/s) '
          write(102,*)'column 11: Number of reflections  '
        CLOSE(102)
      endif
    ENDIF

    !Deallocate local variables
    DEALLOCATE(pLon,pLat)

    fpy=1234567
    open (unit = fpy,file = TRIM(OutDir)//'/'//TRIM(NCOutFile)//               &
                         'PartErrors'//'.py',STATUS='REPLACE')
    write(fpy,'(a)')'import matplotlib.pyplot as plt'
    write(fpy,'(a)')'import numpy as np'
    write(fpy,'(a)')'#fig,ax=plt.subplots()'
    close(fpy)
    Average_Numpart(:)=0                                             
    Average_Value(:)=0.0
    !write(*,*)' Set Initial Wind and Temperature '
    call set_Average_Wind()
    call set_Average_Temperature()
    !write(*,*)' Set Initial Wind and Temperature DONE'
    !write(*,*)'Average_Value(ID_U_WIND)=',Average_Value(ID_U_WIND)
    !write(*,*)'Average_Value(ID_V_WIND)=',Average_Value(ID_V_WIND)
    !write(*,*)'Average_Value(ID_TEMP)=',Average_Value(ID_TEMP)



!        *****   IMIOM        *****
        IF(OilOn)then                                                                                                                                !Initialise oil properties and open output files to write headers
          Average_Numpart(ID_STRANDED) = 0
          Average_Numpart(ID_ALIVE) = numpar
         
          CALL InitOilModel(Average_Value(ID_TEMP),pTS)
          !call OilModel(2,par(:,:2),Average_Numpart(ID_ALIVE),Average_Numpart(ID_STRANDED),sqrt(Average_Value(ID_U_WIND)**2+Average_Value(ID_V_WIND)**2)*WindWeatherFac,&
          !      Average_Value(ID_TEMP),Angle_wrtEast(Average_Value(ID_U_WIND),Average_Value(ID_V_WIND)))
              
               write(*,*)'!--- CL-OGS: Oil Model initiated'
               write(*,*)'!--- CL-OGS: AvWind=',sqrt(Average_Value(ID_U_WIND)**2+Average_Value(ID_V_WIND)**2),'*',WindWeatherFac, &
                                     ' Average_Value(ID_TEMP)=',Average_Value(ID_TEMP)
               write(*,*)'!--- CL-OGS: OILTRANS initialization made at the first iteration, t=',ix(2)
        End If
!        ***** END IMIOM *****


    !Output time spent initializing model
    call CPU_TIME(times(1))
    write(*,'("Time to initialize model = ",f20.10," seconds.")') times(1)
    timeCounts = 0.0      !Initialize time counters to 0.0

  end subroutine ini_LTRANS

  

  subroutine run_External_Timestep()
    use param_mod, only: dt,idt,WriteModelTiming,                     & 
                        numdigits,suffix,numpar,                     & !--- CL-OGS 
                        OpenOceanBoundary, mortality,                & !--- CL-OGS 
                        settlementon,Ext0,WriteParfile,VTurbOn,StrandingDist          !--- CL-OGS 
    use convert_mod, only: x2lon,y2lat                                  !--- CL-OGS   
    use hydro_mod, only: updateHydro,                                 &
                        filenm                                          !--- CL-OGS 
    USE BEHAVIOR_MOD,   ONLY: isOut,isDead                              !--- CL-OGS  
    USE SETTLEMENT_MOD, ONLY: isSettled,isStranded
    integer :: stepIT,ios  
    CHARACTER(len=200) :: namefile                                      !--- CL-OGS 
    real :: before,after
    integer :: n,nunborn,npsetl,npdead,npout,npstrd                  !--- CL-OGS
      
        npsetl=0
        npstrd=0
        npdead=0
        npout=0
        DO n=1,numpar
           !If particle settled or dead
            if(settlementon)then
              if(isSettled(n)) npsetl=npsetl+1
              if(isStranded(n)) npstrd=npstrd+1
            endif
            if(mortality)then
              if ( isDead(n) ) npdead=npdead+1
            endif
            if(OpenOceanBoundary)then
              if(isOut(n)) npout=npout+1
            endif
        ENDDO
        nunborn=numpar-npsetl-npstrd-npdead-npout-Average_Numpart(ID_ALIVE)
        write(*,'(6(a,i8))')'Number of parts alive=',Average_Numpart(ID_ALIVE),' settled=',npsetl,&
                    ' stranded=',npstrd,' dead=',npdead,' out=',npout,' unborn=',nunborn 


      stepIT  = int(dt/idt)                     !number of internal time steps

      IF(WriteModelTiming) call CPU_TIME(before)

      !Read in hydrodynamic model data 
      IF(p > 2) then
        CALL updateHydro()   !do not start updating until 3rd iteration
      ENDIF

      IF(p>2 .and. WriteParfile)then
        write(namefile,'(A,A,A)')'Parfile_', &
        filenm(len(trim(filenm))-len(trim(suffix))-numdigits+1:len(trim(filenm))-len(trim(suffix))),'.csv'
        OPEN(1,FILE=TRIM(namefile),form='formatted',status='unknown', &
             action='write', iostat=ios)
        if ( ios /= 0 ) stop " ERROR OPENING Parfile"
      
        do n=1,numpar
          if( (.not.isOut(n)) .and. (.not.isDead(n)))then
            if(settlementon .and.  StrandingDist<0)then
              write(1,"(3(F18.7,','),2(i13,','),i20)") x2lon(par(n,pX),par(n,pY)),y2lat(par(n,pY)),par(n,pZ), &
                         int(max(par(n,pDOB)-ix(2),0.0)),startpoly(n),n
            else
              write(1,"(3(F18.7,','),2(i13,','),i20)") x2lon(par(n,pX),par(n,pY)),y2lat(par(n,pY)),par(n,pZ), &
                         int(max(par(n,pDOB)-ix(2),0.0)),0,n
            endif
          endif
        enddo
        CLOSE(1)
      ENDIF

      IF(WriteModelTiming) then
        call CPU_TIME(after)
        timeCounts(1) = timeCounts(1) + (after-before)
      ENDIF

      !Prepare external time step values to be used for 
      !  calculating Advection and Turbulence
      ex=0.0
      ex(1) = (p-2)*dt  +Ext0
      ex(2) = (p-1)*dt  +Ext0
      ex(3) = p*dt      +Ext0
      write(*,*)'###################################################'
      do it=1,stepIT

        call run_Internal_Timestep()
      
      enddo !ITloop
     
      write(*,*)'###################################################'

  end subroutine run_External_Timestep



  subroutine run_Internal_Timestep()
    use param_mod, only: idt,iprint,                 &
                         Ext0,SaltTempOn,Wind,       &       !--- CL-OGS 
                         OilOn,WindWeatherFac,numpar !        *****   IMIOM        *****

    use convert_mod, only: x2lon,y2lat                                  !--- CL-OGS   

    use oil_mod, only: OilModel

    INTEGER :: Phase1Time,n,m,ElapsedTime,SprdCase
    !Prepare internal time step values to be used for 
    !  calculating Advection and Turbulence
    ix(1) = ex(2) + DBLE((it-2)*idt)
    ix(2) = ex(2) + DBLE((it-1)*idt)
    ix(3) = ex(2) + DBLE(it*idt)
          
    !********************************************************
    !*                    Particle Loop                     *
    !********************************************************

    call update_particles()

       if(Average_Numpart(ID_ALIVE)>0)then
          if( Wind .and. Average_Numpart(ID_U_WIND).eq.0) then
                call set_Average_Wind()
                !Average_Numpart(ID_U_WIND)=1
          endif
          if( SaltTempOn .and.Average_Numpart(ID_TEMP).eq.0) then
                call set_Average_Temperature()
                !Average_Numpart(ID_TEMP)=1
          endif
          if( Average_Numpart(ID_DEPTH).eq.0 ) then
                call set_Average_Water_Depth()
                !Average_Numpart(ID_DEPTH)=1
          endif
       endif

      !**********************************************************************
      !*     Irish Marine Institute Oil Model                               *
      !**********************************************************************
      if(OilOn)then
       !write(*,*)'call Oil, AvWind,Angle,Temp=',sqrt(Average_Value(ID_U_WIND)**2+Average_Value(ID_V_WIND)**2), &
       !    Angle_wrtEast(Average_Value(ID_U_WIND),Average_Value(ID_V_WIND)),Average_Value(ID_TEMP),Average_Numpart(ID_U_WIND),Average_Numpart(ID_TEMP)
       Average_Numpart(ID_ALIVE)       = Average_Numpart(ID_ALIVE) - Average_Numpart(ID_STRANDED)            !no of particles left in simulation

       call OilModel(2,par(:,:2),Average_Numpart(ID_ALIVE),Average_Numpart(ID_STRANDED), &
            sqrt(Average_Value(ID_U_WIND)**2+Average_Value(ID_V_WIND)**2)*WindWeatherFac,Average_Value(ID_TEMP), &
            Angle_wrtEast(Average_Value(ID_U_WIND),Average_Value(ID_V_WIND)),abs(Average_Value(ID_DEPTH)))

       Average_Numpart(ID_STRANDED)    = 0           !reset no of particles stranded for this timestep
      end if                                         ! OilOn

      !********************************************************

    !********************************************************
    !*                 PRINT OUTPUT TO FILE                 *
    !********************************************************

    printdt=printdt+abs(idt) !--- CL-OGS: added abs() on next line to allow backward-in-time simulation 

    if(printdt.eq.0 .or. printdt.GE.iprint) then
      !write(*,*) 'write output to file, day = ',(DBLE(ix(3))/DBLE(86400)),' ix(3)=',ix(3)
      !ix(3)/86400 = (current model time in seconds) /
      !              (# of seconds in a day)
                call printOutput()

      printdt=0  !reset print counter
    endif
!    close(111)


  end subroutine run_Internal_Timestep  



  subroutine fin_LTRANS()
    use param_mod, only: numpar,writeCSV,outpathGiven,outpath,settlementon, & 
                         NCOutFile,Behavior,OutDir,StrandingDist                           !--- CL-OGS
    use behavior_mod, only: finBehave,getStatus
    use convert_mod, only: x2lon,y2lat
    use hydro_mod, only: finHydro
    use boundary_mod, only: finBoundary
    use random_mod,   only: fin_genrand
    use tension_mod, only : finTensionModule

    !OUTPUT ENDFILE NAME CONSTRUCTION VARIABLE
    CHARACTER(LEN=200) :: efile

    integer :: n,d,h,m,ios
    real :: fintime,s
    double precision :: pLon,pLat
    double precision :: default_stat
    open (unit = fpy,file = TRIM(OutDir)//'/'//TRIM(NCOutFile)//               &
                         'PartErrors'//'.py',POSITION='APPEND')
    write(fpy,'(a)')'#plt.show()'
    close(fpy)
    !Write final positions and status to endfile.csv
    IF(writeCSV)THEN

      if(outpathGiven)then
        efile = TRIM(outpath)//'/'//TRIM(NCOutFile)//'-endfile.csv'
      else
        efile = TRIM(NCOutFile)//'-endfile.csv'
      endif 

      write(*,*) 'write ',efile
      !--- CL-OGS: OILTRANS : moved iunit 3 -> 333 of as 3 already used as well by OILTRANS developments
      open(333,FILE=efile,STATUS='REPLACE')
        3 format(I7,',',I7,',',I7,',', F13.8,",",F13.8,",",F13.8,",",I7)
        4 format(I7,',',F13.8,',',F13.8,',',F13.8,',',I7)
        do n=1,numpar
          pLon = x2lon(par(n,pX),par(n,pY))
          pLat = y2lat(par(n,pY))
          default_stat=par(n,pStatus)
          par(n,pStatus) = getStatus(n,default_stat)
          if(settlementon .and.  StrandingDist<0)then
            write(333,3) startpoly(n),endpoly(n),int(par(n,pStatus)),pLat,pLon,  &
                       par(n,pZ),int(par(n,pLifespan))
          else
            write(333,4) int(par(n,pStatus)),pLat,pLon,par(n,pZ),int(par(n,pLifespan))
          endif
        enddo
      close(333)

    ENDIF
    IF(Behavior.eq.997)then
        OPEN(1,FILE=TRIM(outpath)//'/PartAtSurf_'//TRIM(NCOutFile)//'.csv',&
           form='formatted',status='unknown', action='write', iostat=ios)
        if ( ios /= 0 ) stop " ERROR OPENING Parfile" 
        do n=1,numpar
              write(1,"(3(F18.7,','),2(i13,','),i20)") &
                         PartAtSurf(n,1),PartAtSurf(n,2),PartAtSurf(n,3), &
                         int(PartAtSurf(n,4)), 0,n
        enddo
        CLOSE(1)
    ENDIF



    !DEALLOCATE LOCAL VARIABLES
    DEALLOCATE(par)
    IF(ALLOCATED(hitBottom)) DEALLOCATE(hitBottom)
    IF(ALLOCATED(startpoly)) DEALLOCATE(startpoly)
    IF(ALLOCATED(endpoly  )) DEALLOCATE(endpoly)
    IF(ALLOCATED(hitLand  )) DEALLOCATE(hitLand)
    IF(ALLOCATED(P_Salt   )) DEALLOCATE(P_Salt)
    IF(ALLOCATED(P_Temp   )) DEALLOCATE(P_Temp)
    IF(ALLOCATED(P_coastdist)) DEALLOCATE(P_coastdist)
    IF(ALLOCATED(parIniDepth )) DEALLOCATE(parIniDepth)
    IF(ALLOCATED(PartAtSurf )) DEALLOCATE(PartAtSurf)
    IF(ALLOCATED(P_oldLev)) DEALLOCATE(P_oldLev)     

    !DEALLOCATE MODULE VARIABLES
    call fin_genrand()
    call finTensionModule()
    call finBehave()
    call finHydro()
    call finBoundary()
    !Calculate model run time and output to screen before exiting
    call CPU_TIME(fintime)
    d = int(fintime/86400.0)             !# of full days that the model ran
    h = int((fintime - real(d*86400))/3600.0)       !# of hours   (minus days)
    m = int((fintime - real(d*86400 - h*3600))/60.0)!# of minutes (minus days and hours)
    s =  fintime - REAL(d*86400) - REAL(h*3600) - REAL(m*60) !# of seconds (- days, hrs and mins)

    11 format('Time to run model = ',i3,' days ',i5,' hours ',i5,              &
              ' minutes and ',f10.3,' seconds.')
    12 format('Time to run model = ',i5,' hours ',i5,' minutes and ',f10.3,     &
              ' seconds.')
    13 format('Time to run model = ',i5,' minutes and ',f10.3,' seconds.')
    14 format('Time to run model = ',f10.3,' seconds.')

    if(fintime > 86400.0)then
      write(*,11) d,h,m,s
    elseif(fintime > 3600.0)then
      write(*,12) h,m,s
    elseif(fintime > 60.0)then
      write(*,13) m,s
    else
      write(*,14) fintime
    endif

    write(*,'(/,A)') '****** END LTRANS *******'

  end subroutine fin_LTRANS


  subroutine update_particles()
    USE PARAM_MOD,      ONLY: numpar,us,ws,OilOn
   !$ use OMP_LIB

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(us) :: Pwc_zb,Pwc_zc,Pwc_zf
    DOUBLE PRECISION, DIMENSION(ws) :: Pwc_wzb,Pwc_wzc,Pwc_wzf
    DOUBLE PRECISION frstpriv_Average_Value(NUM_ID_VALUES)
    INTEGER frstpriv_Average_Numpart(NUM_ID_NUMPARTS)
    INTEGER n
    frstpriv_Average_Numpart(:)=0
    frstpriv_Average_Value(:)=0.0
    Average_Numpart(:)=0       
    Average_Value(:) = 0.0     

   !$OMP PARALLEL DEFAULT(NONE) &
   !$OMP SHARED(numpar,par) &
   !$OMP PRIVATE(Pwc_zb,Pwc_zc,Pwc_zf,Pwc_wzb,Pwc_wzc,Pwc_wzf) &
   !$OMP FIRSTPRIVATE (frstpriv_Average_Value,frstpriv_Average_Numpart)  &
   !$OMP REDUCTION(+:Average_Value) &    
   !$OMP REDUCTION(+:Average_Numpart) 

   !$OMP DO PRIVATE(n)  
    DO n=1,numpar
      call update_single_particle(n,frstpriv_Average_Value,frstpriv_Average_Numpart, &
                                  Pwc_zb,Pwc_zc,Pwc_zf,Pwc_wzb,Pwc_wzc,Pwc_wzf)
      par(n,ppX) = par(n,pX)
      par(n,ppY) = par(n,pY)
      par(n,ppZ) = par(n,pZ)
      par(n,pX)  = par(n,pnX) 
      par(n,pY)  = par(n,pnY)
      par(n,pZ)  = par(n,pnZ)
       if ((isnan(par(n,pX)) .or. isnan(par(n,pY)) )          &
        .or.  (  isnan(par(n,pnX)).or.isnan(par(n,pnY)) ) )then
         write(*,*)'n=',n
         write(*,*)'Xpar=',par(n,pX),par(n,pnX)
         write(*,*)'Ypar=',par(n,pY),par(n,pnY)
         write(*,*)'Zpar=',par(n,pZ),par(n,pnZ)
         stop 'found NaN in updateparloc'
       endif
    ENDDO !end loop for each particle
   !$OMP END DO          

   !DEALLOCATE(Pwc_zb,Pwc_zc,Pwc_zf)
   !DEALLOCATE(Pwc_wzb,Pwc_wzc,Pwc_wzf)

    Average_Value(:)=Average_Value(:)+frstpriv_Average_Value(:)
    Average_Numpart(:)=Average_Numpart(:)+frstpriv_Average_Numpart(:)

   !$OMP END PARALLEL 

      IF(Average_Numpart(ID_U_WIND).ge.1)then
         Average_Value(ID_U_WIND)=Average_Value(ID_U_WIND)/float(Average_Numpart(ID_U_WIND))
         Average_Value(ID_V_WIND)=Average_Value(ID_V_WIND)/float(Average_Numpart(ID_V_WIND))
         !write(*,*)'***** AvUwd=',Average_Value(ID_U_WIND),' AvVwd=',Average_Value(ID_V_WIND)
      ENDIF
      IF(Average_Numpart(ID_TEMP).ge.1)then
         Average_Value(ID_TEMP)=Average_Value(ID_TEMP)/float(Average_Numpart(ID_TEMP))
      ENDIF
      IF(Average_Numpart(ID_DEPTH).ge.1)then
         Average_Value(ID_DEPTH)=Average_Value(ID_DEPTH)/float(Average_Numpart(ID_DEPTH))
      ENDIF

    !IMIOM
    if(OilOn)then  !set oil particle vertical location to surface
      do n=1,numpar
        par(n,ppZ) = 0.0
        par(n,pZ) = 0.0
      end do
    end if
    !END IMIOM


  end subroutine update_particles

  !-----------------------------------------------------------------------

  subroutine update_single_particle(n,private_Average_Value,private_Average_Numpart, &
                                  Pwc_zb,Pwc_zc,Pwc_zf,Pwc_wzb,Pwc_wzc,Pwc_wzf)

    USE PARAM_MOD,      ONLY: numpar,us,ws,idt,HTurbOn,VTurbOn,settlementon,   &
                              Behavior,SaltTempOn,OpenOceanBoundary,Swimdepth, &
                              TrackCollisions,WriteModelTiming,mortality,      &
                              ErrorFlag,                                       &
                              WriteCurrents,NCOutFile,Outdir,                  &!--- CL-OGS: output Currents if requested
                              StrandingDist,storedincolor,                     &!--- CL-OGS: Stranding/settlement optionals
                              Zgrid,Zgrid_depthinterp,                         &!--- CL-OGS: coupling with MITgcm'Z-grid bathymetry and fields
                              iprint,WindDriftFac,WindDriftDev,StokDriftFac,   &!--- CL-OGS 
                              BndOut,dt,read_GrainSize,WindIntensity,z0,  &
                              !readNetcdfSwdown, &
                              Write_Temp_min_max_ins,Write_Salt_min_max_ins, &
                              Write_Poly_Presence, &                            !--- CL-OGS
!        ***** IMIOM *****
                              OilOn,Wind,SigWaveHeight,MeanWavePeriod,UWind_10,&
                              VWind_10,PeakDirection,PeakWaveLength,pi,Stokes, &
                              constDens,Write_coastdist

!        ***** END IMIOM *****
    USE SETTLEMENT_MOD, ONLY: isSettled,testSettlement,isStranded,p_Stranding
    USE BEHAVIOR_MOD,   ONLY: updateStatus,behave,setOut,isOut,isDead,die
    USE BOUNDARY_MOD,   ONLY: mbounds,ibounds,intersect_reflect,Get_coastdist
    USE CONVERT_MOD,    ONLY: x2lon,y2lat
    USE HTURB_MOD,      ONLY: HTurb
    USE VTURB_MOD,      ONLY: VTurb
    USE INT_MOD,        ONLY: linint,polintd
    USE HYDRO_MOD,      ONLY: setEle,setInterp,getInterp,getSlevel,getWlevel,  &
                              WCTS_ITPI, &
                              getKRlevel,getDepth,                             &!--- CL:OGS
                              !outputdetails_closernode, &
                              getP_r_element          !--- CL:OGS
    use oil_mod,  only: STOKESDRIFT

    !$ use OMP_LIB

#include "VAR_IDs.h"

    IMPLICIT NONE

    DOUBLE PRECISION, DIMENSION(us) :: Pwc_zb,Pwc_zc,Pwc_zf
    DOUBLE PRECISION, DIMENSION(ws) :: Pwc_wzb,Pwc_wzc,Pwc_wzf
    DOUBLE PRECISION private_Average_Value(NUM_ID_VALUES)
    INTEGER private_Average_Numpart(NUM_ID_NUMPARTS)
    INTEGER n

    ! Iteration Variables
    INTEGER :: i,deplvl
    INTEGER :: klev,Fstlev,NumInterpLvl,ixnum,nklev,k                     !--- CL:OGS

    ! Particle tracking
    DOUBLE PRECISION :: Xpar,Ypar,Zpar,newXpos,newYpos,newZpos,P_zb,P_zc,P_zf, &
      nP_depth,P_depth,P_angle,P_zeta,P_zetab,P_zetac,P_zetaf,P_surfdist,ey(3),        &
      AdvecUwind,AdvecVwind
    
    ! Behavior and Turbulence
    DOUBLE PRECISION :: TurbHx,TurbHy,TurbV,Behav,XBehav,YBehav,ZBehav
    LOGICAL :: bott   ! for Behavior 7 along with XBehav,YBehav,ZBehav
    !DOUBLE PRECISION ::P_swdown

    ! Boundaries
    INTEGER :: intersectf,skipbound,in_island,inbounds,reflects,inpoly
    INTEGER :: saveintersectf
    DOUBLE PRECISION :: fintersectX,fintersectY,freflectX,freflectY,coastdist, &
      Xpos,Ypos,nXpos,nYpos,pushedup,reflectsup,reflectinf
    INTEGER :: island,mainbound
    LOGICAL :: isWater,waterFlag

    ! Advection
    DOUBLE PRECISION :: AdvectX,AdvectY,AdvectZ,maxpartdepth,minpartdepth,     &
      kn1_u,kn1_v,kn1_w,kn2_u,kn2_v,kn2_w,kn3_u,kn3_v,kn3_w,kn4_u,kn4_v,kn4_w, &
      P_V,P_U,P_W,UAD,VAD,WAD,x1,x2,x3,y1,y2,y3,z1,z2,z3,slope,length,& 
      Ttemp,Tsalt,LarvSize
    ! Error Handling Variables
    INTEGER :: ele_err
    INTEGER :: conflict,conflict_tmp,m_nestdeg,i_nestdeg    !--- CL:OGS
    character(LEN=40) :: fstring                                     !--- CL:OGS
    logical:: dbg                                            !--- CL:OGS 
    DOUBLE PRECISION :: bndintersect(12)                      !--- CL:OGS
    INTEGER foundsetEleproblem
    CHARACTER(LEN=1) :: col(4)
    CHARACTER(LEN=10) :: ErrorName,annotate
    LOGICAL :: docycle
    !IMIOM
        DOUBLE PRECISION::      P_hsig,P_tm01,P_Uwind,P_Vwind,P_pdir,P_wlen,UWindDrift,&
                                VWindDrift,alpha,PWind,P_Uw,P_Vw,Uadw,Vadw,            &
                                UStokesDrift,VStokesDrift
        DOUBLE PRECISION::  kn1_uw,kn1_vw,kn2_uw,kn2_vw,kn3_uw,kn3_vw,kn4_uw,kn4_vw
        INTEGER:: ElapsedTime
    !END IMIOM
    INTEGER:: posfactor
    DOUBLE PRECISION :: PTemptmp
    !P_swdown = 0.0

    dbg=.FALSE.
    col(1)='b'
    col(2)='g'
    col(3)='r'
    col(4)='k'
    coastdist=9e12

      !$OMP MASTER
      IF(WriteModelTiming) call CPU_TIME(times(2))
      !$OMP END MASTER

      ! *********************************************************
      ! *                                                       *
      ! *        Update Particle Age and Characteristics        *
      ! *                                                       *
      ! *********************************************************

      !If the particle is not yet released, set new location to 
      !  current location, and cycle to next particle

!--- CL-OGS:   modified to allow backward-in-time simulation :
      if(idt>0)then !--- CL-OGS: for forward run
        if(ix(3) < par(n,pDOB))then 
          par(n,pnX) = par(n,pX)
          par(n,pnY) = par(n,pY)
          par(n,pnZ) = par(n,pZ)
          return
        endif
      else !--- CL-OGS: for backward run
        if(ix(3) > par(n,pDOB))then 
          par(n,pnX) = par(n,pX)
          par(n,pnY) = par(n,pY)
          par(n,pnZ) = par(n,pZ)
          return
        endif
      endif
      !Update particle age
      par(n,pAge) = par(n,pAge) + float(idt)

      !Update particle status if settled or dead
      CALL updateStatus(par(n,pAge),n)

      !If particle settled or dead, skip tracking
      if(settlementon)then
        if ( isSettled(n) ) return
      endif
 
!--- CL-OGS:   !If particle stranded, skip tracking
      if(settlementon)then
        if(isStranded(n)) return
      endif

      if(mortality)then
        if ( isDead(n) ) return
      endif

      !If there are open ocean boundaries and the current
      !  particle has exited the model domain via them, skip it
      if(OpenOceanBoundary)then
        if(isOut(n)) return
      endif

      !IMIOM
      if(OilOn)then
          if(par(n,pStatus) == -2) then
              return                               !if stranded
         else
            par(n,pStatus)   = 1          !floating oil
        end if
      end if
      !END IMIOM
      ! *********************************************************
      ! *                                                       *
      ! *          Find Element that Contains Particle          *
      ! *                                                       *
      ! *********************************************************

      !Find rho/u/v elements in which the particle is located.

      Xpar = par(n,pX)
      Ypar = par(n,pY)

!--- CL-OGS: computing Zpar and klev needed for setEle getInterp, testSettlement, 
!--- CL-OGS:                               intersect_reflect, mbounds and ibounds 
      Zpar = par(n,pZ)
      if(Zgrid)then 
        klev=getKRlevel(par(n,pZ))
        if(P_oldLev(n)>=0) P_oldLev(n)=klev
      else
        klev=1
      endif     
      !Determine which Rho, U, & V elements the particle is in
      ele_err=-666
      CALL setEle(Xpar,Ypar,Zpar,n,it,1,ele_err)
      !If the particle was not found to be within an element,
      !  write a message to the screen and discontinue the program
      IF(ele_err .ne.  0)THEN
        !write(ErrorName,'(a)')'setEle'
              write(*,'(a,i12,a,i2,a,3f7.2,6(a,f7.2))') &
                  'ERROR it',it,'on 1st setEle call, error number is ',ele_err, &
                 'local old, present and next depth is ',&
                    par(n,pZ),P_depth,nP_depth, &
                 'at old Z=',par(n,pZ),' ZW(',getKRlevel(par(n,pZ)),')=', &
                            Pwc_wzc(getKRlevel(par(n,pZ))),  &
                  'at read Z=',Zpar  ,' ZW(',getKRlevel(Zpar),')=', &
                            Pwc_wzc(getKRlevel(Zpar)) 
        call handleERROR('setEle    ',1,ele_err,n,     &
                     par(n,pX),par(n,pY),par(n,pZ), &
                     par(n,pX),par(n,pY),par(n,pZ),docycle  )
        if(docycle) return
      ENDIF

      
      !Set Interpolation Values for the current particle
      CALL setInterp(Xpar,Ypar,n)

      ! *********************************************************
      ! *                                                       *
      ! *       Ensure Particle is Within Verticle Bounds       *
      ! *                                                       *
      ! *********************************************************

      !Find depth, angle, and sea surface height at particle location
      conflict=1
      if(.not. Zgrid)then ! (ROMS sigma-level-grid :)
        P_depth = DBLE(-1.0)* getInterp(Xpar,Ypar,VAR_ID_depth,klev)
        P_angle = getInterp(Xpar,Ypar,VAR_ID_angle,klev)
        Fstlev=1
      else !             (Zgrid :)                                    !--- CL-OGS:  for Z-grid model P_depth is computed by the 
        call getDepth(Xpar,Ypar,n,it,P_depth,Fstlev,conflict)         !--- CL-OGS:   new routine getDepth returning as well the
        if(Fstlev>us)then                                    
          call handleERROR('Fstlev    ',1,conflict,n,     &
                     par(n,pX),par(n,pY),par(n,pZ), &
                     Xpar,Ypar,par(n,pZ),docycle)
           call die(n)
          if(docycle) return 
        elseif(conflict.ne.1)then
          call handleERROR('dDepth    ',1,conflict,n,     &
                     par(n,pX),par(n,pY),par(n,pZ), &
                     Xpar,Ypar,par(n,pZ),docycle)
           call die(n)
           if(docycle) return 
        endif
        P_angle=0
      endif
      
      !**********   Sea level zeta  *****************************
      if(p.eq.1)then                        !--- CL-OGS: added in v.Zlev, but is it right ?? 
      P_zetab = getInterp(Xpar,Ypar,VAR_ID_zetab,klev)     !--- CL-OGS: added in v.Zlev, but is it right ?? 
      P_zetac = getInterp(Xpar,Ypar,VAR_ID_zetab,klev)     !--- CL-OGS: added in v.Zlev, but is it right ?? 
      P_zetaf = getInterp(Xpar,Ypar,VAR_ID_zetac,klev)     !--- CL-OGS: added in v.Zlev, but is it right ?? 
      else                                  !--- CL-OGS: added in v.Zlev, but is it right ?? 
      P_zetab = getInterp(Xpar,Ypar,VAR_ID_zetab,klev)
      P_zetac = getInterp(Xpar,Ypar,VAR_ID_zetac,klev)
      P_zetaf = getInterp(Xpar,Ypar,VAR_ID_zetaf,klev)
      endif              
      P_surfdist=Zpar-P_zetac
      if(P_surfdist>0)P_surfdist=Zpar-P_zetab
      if(P_surfdist>0)P_surfdist=Zpar-P_zetaf
      P_surfdist=min(-0.000001,P_surfdist)

!--- CL-OGS : shouldn'be moved here the SETTLEMENT SECTION ?            

      ! *********************************************************
      ! *                                                       *
      ! *             Create Matrix of Z-coordinates            *
      ! *                                                       *
      ! *********************************************************

      !Create matrix of z-coordinates at particle and at each node for
      !  back, center, forward times
      do i=1,us

        !Rho-coordinate depths at particle location
        Pwc_zb(i)=getSlevel(P_zetab,P_depth,i)
        Pwc_zc(i)=getSlevel(P_zetac,P_depth,i)
        Pwc_zf(i)=getSlevel(P_zetaf,P_depth,i)

        !W-coordinate depths at particle location
        Pwc_wzb(i)= getWlevel(P_zetab,P_depth,i)
        Pwc_wzc(i)= getWlevel(P_zetac,P_depth,i)
        Pwc_wzf(i)= getWlevel(P_zetaf,P_depth,i)

      enddo

      !W-coordinate depths at particle location (cont.)
      Pwc_wzb(ws)= getWlevel(P_zetab,P_depth,ws)
      Pwc_wzc(ws)= getWlevel(P_zetac,P_depth,ws)
      Pwc_wzf(ws)= getWlevel(P_zetaf,P_depth,ws)
      if(Zgrid .and. .not.Zgrid_depthinterp)then !--- CL-OGS
        if(Pwc_wzb(Fstlev).ne. P_depth .or. &
        Pwc_wzc(Fstlev).ne. P_depth .or. &
        Pwc_wzf(Fstlev).ne. P_depth ) then
         write(*,*)'ERROR OR NOT?'
         stop
        endif
      endif      


      ! *********************************************************
      ! *                                                       *
      ! *             Prepare for Particle Movement             *
      ! *                                                       *
      ! *********************************************************

      AdvectX = 0.0
      AdvectY = 0.0
      AdvectZ = 0.0
      TurbHx = 0.0
      TurbHy = 0.0
      TurbV = 0.0
      Behav = 0.0 !BEUG? is not used...
      XBehav=0.0 ! BEUG? wasn't initialised...
      YBehav=0.0 ! BEUG? wasn't initialised...
      ZBehav=0.0 ! BEUG? wasn't initialised...
      P_Uwind=0.0
      P_Vwind=0.0
      AdvecUwind=0.0
      AdvecVwind=0.0
 
      !IMIOM
      !If(OilOn)then
              if(Wind)then
                      UWindDrift = 0.0
                      VWindDrift = 0.0
                      UStokesDrift = 0.0
                      VStokesDrift = 0.0
              end if
      !end if     
 
      ! *********************************************************
      ! *                                                       *
      ! *                       ADVECTION                       *
      ! *                                                       *
      ! *********************************************************


      !$OMP MASTER
      IF(WriteModelTiming) call CPU_TIME(times(3))
      !$OMP END MASTER
      

      !-------------------------------
      ! Check if particle location above or below boundary  
      maxpartdepth = Pwc_wzb(Fstlev)
      if(Zgrid .and.Zgrid_depthinterp) maxpartdepth = max(maxpartdepth,P_depth)
      minpartdepth = Pwc_wzb(ws)
      if(Pwc_wzc(ws) .LT. minpartdepth) minpartdepth = Pwc_wzc(ws)
      if(Pwc_wzf(ws) .LT. minpartdepth) minpartdepth = Pwc_wzf(ws) 

      if( Behavior.eq.999 .or.Behavior.eq.998 )then 
        if(Behavior.eq.999)then                                   !--- CL-OGS: keep constant depth below the sea surface height zeta
           if( minpartdepth+parIniDepth(n) >=maxpartdepth) then 
             Zpar = minpartdepth+parIniDepth(n) 
           else !--- CL-OGS: unless sea bottom higher than zeta+parIniDepth, then:
             Zpar = maxpartdepth + DBLE(0.00001)      
             IF(TrackCollisions) hitBottom(n) = hitBottom(n) + 1
           endif
        elseif(Behavior.eq.998 )then                            !--- CL-OGS: keep constant depth below zero
           Zpar=parIniDepth(n)
           if( minpartdepth<Zpar ) then
             Zpar = minpartdepth - DBLE(0.00001)
           endif 
           if(maxpartdepth>Zpar)then
             Zpar = maxpartdepth + DBLE(0.00001)      
             IF(TrackCollisions) hitBottom(n) = hitBottom(n) + 1
           endif
        endif
        par(n,pZ)  = Zpar
        par(n,pnZ) = Zpar
        P_zb = Zpar
        P_zc = Zpar
        P_zf = Zpar 
      else 
        ! If particle location above or below boundary,
        ! place just within boundary (1 mm)
        if( par(n,pZ) .LT. maxpartdepth ) then
          par(n,pZ) = maxpartdepth + DBLE(0.00001)
        endif             
        if (par(n,pZ).LT.P_depth) then
          par(n,pZ) = P_depth + DBLE(0.00001)
          IF(TrackCollisions) hitBottom(n) = hitBottom(n) + 1
        endif

        P_zb = par(n,pZ)
        P_zc = par(n,pZ)
        P_zf = par(n,pZ)
        if (par(n,pZ).GT.P_zetab) P_zb = P_zetab - DBLE(0.00001)
        if (par(n,pZ).GT.P_zetac) P_zc = P_zetac - DBLE(0.00001)
        if (par(n,pZ).GT.P_zetaf) P_zf = P_zetaf - DBLE(0.00001)

        ey(1) = P_zb
        ey(2) = P_zc
        ey(3) = P_zf
        Zpar = polintd(ex,ey,3,ix(2))

        !update particle
        par(n,pZ) = Zpar

      endif

      do i=1,us 
        if ((Zpar .LT. Pwc_zb(i)) .OR. &   
            (Zpar .LT. Pwc_zc(i)) .OR. &
            (Zpar .LT. Pwc_zf(i))      ) exit
      enddo
        
      !-------------------------------
! it appears that, in the original v2b version, 
! the last three calls of find_currents when particle is not in the log layer 
! compute the interpolation of the currents by calling WCTS_ITPI at z position
! P_zc/b/f BASED ON Zpar INSTEAD OF z1 (resp z2,z3) !
! (while correctly interpolating at z1 (resp z2,z3) if particle is in the log
! layer) a suggestion of modification to solve this problem is made below:

      !Find advection currents at original coordinates
      ixnum=1
      CALL find_currents(Xpar,Ypar,Zpar,Pwc_zb,Pwc_zc,Pwc_zf,Pwc_wzb,Pwc_wzc, &
        Pwc_wzf,P_zb,P_zc,P_zf,ex,ix,p,ixnum,Uad,Vad,Wad,Fstlev,P_depth,n)
      !Store advection currents at original coordinates
      kn1_u = Uad
      kn1_v = Vad
      kn1_w = Wad
      !-------------------------------

      !Estimate new coordinates for next RK position
      x1 = Xpar + (Uad*cos(P_angle) - Vad*sin(P_angle)) * DBLE(idt)/DBLE(2) 
      y1 = Ypar + (Uad*sin(P_angle) + Vad*cos(P_angle)) * DBLE(idt)/DBLE(2) 
      if( Behavior.ne.999 .and.Behavior.ne.998 )then
        z1 = Zpar +  Wad * DBLE(idt)/DBLE(2) 
        if(z1 .GT. minpartdepth) z1 = minpartdepth - DBLE(0.000001)
        if(z1 .LT. maxpartdepth) z1 = maxpartdepth + DBLE(0.000001)
        P_zb = max(z1,P_depth+0.00001)                         !SUGGESTION TO SOLVE PROBLEM
        P_zc = max(z1,P_depth+0.00001)                         !SUGGESTION TO SOLVE PROBLEM
        P_zf = max(z1,P_depth+0.00001)                         !SUGGESTION TO SOLVE PROBLEM
        if (P_zb.GT.P_zetab) P_zb = P_zetab - DBLE(0.00001)    !SUGGESTION TO SOLVE PROBLEM
        if (P_zc.GT.P_zetac) P_zc = P_zetac - DBLE(0.00001)    !SUGGESTION TO SOLVE PROBLEM
        if (P_zf.GT.P_zetaf) P_zf = P_zetaf - DBLE(0.00001)    !SUGGESTION TO SOLVE PROBLEM
      else
        z1=Zpar
      endif 

      !Find advection currents at estimated next RK position
      ixnum=2
      CALL find_currents(x1,y1,z1,Pwc_zb,Pwc_zc,Pwc_zf,Pwc_wzb,Pwc_wzc,        &
           Pwc_wzf,P_zb,P_zc,P_zf,ex,ix,p,ixnum,Uad,Vad,Wad,Fstlev,P_depth,n)
      !Store advection currents at 2nd RK position
      kn2_u = Uad
      kn2_v = Vad
      kn2_w = Wad
      !-------------------------------

      !Estimate new coordinates for next RK position
      x2 = Xpar + (Uad*cos(P_angle) - Vad*sin(P_angle)) * DBLE(idt)/DBLE(2) 
      y2 = Ypar + (Uad*sin(P_angle) + Vad*cos(P_angle)) * DBLE(idt)/DBLE(2) 
      if( Behavior.ne.999  .and.Behavior.ne.998 )then
        z2 = Zpar +  Wad * DBLE(idt)/DBLE(2)
        if(z2 .GT. minpartdepth) z2 = minpartdepth - DBLE(0.000001)
        if(z2 .LT. maxpartdepth) z2 = maxpartdepth + DBLE(0.000001)
        P_zb = max(z2,P_depth+0.00001)                       !SUGGESTION TO SOLVE PROBLEM
        P_zc = max(z2,P_depth+0.00001)                       !SUGGESTION TO SOLVE PROBLEM
        P_zf = max(z2,P_depth+0.00001)                       !SUGGESTION TO SOLVE PROBLEM
        if (P_zb.GT.P_zetab) P_zb = P_zetab - DBLE(0.00001)    !SUGGESTION TO SOLVE PROBLEM
        if (P_zc.GT.P_zetac) P_zc = P_zetac - DBLE(0.00001)    !SUGGESTION TO SOLVE PROBLEM
        if (P_zf.GT.P_zetaf) P_zf = P_zetaf - DBLE(0.00001)    !SUGGESTION TO SOLVE PROBLEM
      else
        z2=Zpar
      endif

      !Find advection currents at estimated next RK position
      ixnum=2
      CALL find_currents(x2,y2,z2,Pwc_zb,Pwc_zc,Pwc_zf,Pwc_wzb,Pwc_wzc,        &
           Pwc_wzf,P_zb,P_zc,P_zf,ex,ix,p,ixnum,Uad,Vad,Wad,Fstlev,P_depth,n)
      kn3_u = Uad
      kn3_v = Vad
      kn3_w = Wad
      !-------------------------------

      !Calculate the coordinates at the final position
      x3 = Xpar + (Uad*cos(P_angle) - Vad*sin(P_angle)) * DBLE(idt) 
      y3 = Ypar + (Uad*sin(P_angle) + Vad*cos(P_angle)) * DBLE(idt) 
      if( Behavior.ne.999  .and.Behavior.ne.998 )then
        z3 = Zpar + Wad * DBLE(idt)
        if(z3 .GT. minpartdepth) z3 = minpartdepth - DBLE(0.000001)
        if(z3 .LT. maxpartdepth) z3 = maxpartdepth + DBLE(0.000001)
        P_zb = max(z3,P_depth+0.00001)                         !SUGGESTION TO SOLVE PROBLEM
        P_zc = max(z3,P_depth+0.00001)                         !SUGGESTION TO SOLVE PROBLEM
        P_zf = max(z3,P_depth+0.00001)                         !SUGGESTION TO SOLVE PROBLEM
        if (P_zb.GT.P_zetab) P_zb = P_zetab - DBLE(0.00001)    !SUGGESTION TO SOLVE PROBLEM
        if (P_zc.GT.P_zetac) P_zc = P_zetac - DBLE(0.00001)    !SUGGESTION TO SOLVE PROBLEM
        if (P_zf.GT.P_zetaf) P_zf = P_zetaf - DBLE(0.00001)    !SUGGESTION TO SOLVE PROBLEM
      else
        z3=Zpar
      endif

      !Find advection currents at the final position
      ixnum=3

      CALL find_currents(x3,y3,z3,Pwc_zb,Pwc_zc,Pwc_zf,Pwc_wzb,Pwc_wzc,        &
           Pwc_wzf,P_zb,P_zc,P_zf,ex,ix,p,ixnum,Uad,Vad,Wad,Fstlev,P_depth,n) 

      !Store advection currents at final position
      kn4_u = Uad
      kn4_v = Vad
      kn4_w = Wad
      !-------------------------------

      !Use the RK formula to get the final Advection values
      P_U = (kn1_u + DBLE(2.0)*kn2_u + DBLE(2.0)*kn3_u + kn4_u)/DBLE(6.0)
      P_V = (kn1_v + DBLE(2.0)*kn2_v + DBLE(2.0)*kn3_v + kn4_v)/DBLE(6.0)
      if( Behavior.ne.999 .and. Behavior.ne.998)then
        P_W = (kn1_w + DBLE(2.0)*kn2_w + DBLE(2.0)*kn3_w + kn4_w)/DBLE(6.0)
      else
        P_W=0.0
      endif
      AdvectX = idt*(P_U*cos(P_angle) - P_V*sin(P_angle))
      AdvectY = idt*(P_U*sin(P_angle) + P_V*cos(P_angle))
      AdvectZ = idt*P_W


      ! *********************************************************
      ! *                 IMIOM                                 *
      ! *                WIND DRIFT                                                *
      ! *                                                       *
      ! *********************************************************
          if(Wind)then
            CALL find_winds(Xpar,Ypar,ex,ix,p,1,Uadw,Vadw,n)
         
            !Store advection currents at original coordinates
            kn1_uw = Uadw
            kn1_vw = Vadw
         
            !Estimate new coordinates for next RK position
            x1 = Xpar + (Uadw*cos(P_angle) - Vadw*sin(P_angle)) * DBLE(idt)/DBLE(2)
            y1 = Ypar + (Uadw*sin(P_angle) + Vadw*cos(P_angle)) * DBLE(idt)/DBLE(2)
         
            !Find advection currents at estimated next RK position
            CALL find_winds(x1,y1,ex,ix,p,2,Uadw,Vadw,n)
         
            !Store advection currents at 2nd RK position
            kn2_uw = Uadw
            kn2_vw = Vadw
         
            !Estimate new coordinates for next RK position
            x2 = Xpar + (Uadw*cos(P_angle) - Vadw*sin(P_angle)) * DBLE(idt)/DBLE(2)
            y2 = Ypar + (Uadw*sin(P_angle) + Vadw*cos(P_angle)) * DBLE(idt)/DBLE(2)
         
            !Find advection currents at estimated next RK position
            CALL find_winds(x2,y2,ex,ix,p,2,Uadw,Vadw,n)
         
            !Store advection currents at 3rd RK position
            kn3_uw = Uadw
            kn3_vw = Vadw
         
            !Calculate the coordinates at the final position
            x3 = Xpar + (Uadw*cos(P_angle) - Vadw*sin(P_angle)) * DBLE(idt)
            y3 = Ypar + (Uadw*sin(P_angle) + Vadw*cos(P_angle)) * DBLE(idt)
         
         
            !Find advection currents at the final position
            CALL find_winds(x3,y3,ex,ix,p,3,Uadw,Vadw,n)
         
            !Store advection currents at final position
            kn4_uw = Uadw
            kn4_vw = Vadw
         
            !Use the RK formula to get the final Advection values
            P_Uw = (kn1_uw + DBLE(2.0)*kn2_uw + DBLE(2.0)*kn3_uw + kn4_uw)/DBLE(6.0)
            P_Vw = (kn1_vw + DBLE(2.0)*kn2_vw + DBLE(2.0)*kn3_vw + kn4_vw)/DBLE(6.0)
         
!          if(P_Uw .eq. 0.0)P_Uw = 0.0001
!          if(P_Vw .eq. 0.0)P_Vw = 0.0001
!        
!          Pwind = sqrt(P_Uw**2.0 + P_Vw**2.0)
!          alpha = atan(P_Vw/P_Uw)
!        
!          if(P_Uw .lt. 0)then
!              if(P_Vw .lt. 0)then
!                  alpha = ((3.0/2.0)*pi) - alpha
!              elseif(P_Vw .gt. 0)then
!                  alpha = ((1.0/2.0)*pi) - alpha
!              end if
!          elseif(P_Uw .gt. 0)then
!              if(P_Vw .lt. 0)then
!                  alpha = ((3.0/2.0)*pi) - alpha
!              elseif(P_Vw .gt. 0)then
!                  alpha = ((1.0/2.0)*pi) - alpha
!              end if
!          end if
!        
!          alpha = alpha + P_angle
 
         
                  !calculate wind vector magnitude
                  PWind = sqrt((P_Uw*cos(P_angle) - P_Vw*sin(P_angle))**2.0    & !u rectified to E-W orientation
                              + (P_Uw*sin(P_angle) + P_Vw*cos(P_angle))**2.0)    !v rectified to N-S orientation
!        
                  !calculate vector angle alpha
                  alpha = Angle_wrtEast( (P_Uw*cos(P_angle) - P_Vw*sin(P_angle)), &        !true 0angle of wind vector in cartesian (X/Y) plane
                                 (P_Uw*sin(P_angle) + P_Vw*cos(P_angle))   )        !not ROMS U/V plane.

              CALL setInterp(Xpar,Ypar,n)

                if ((isnan(alpha) ))then 
                  write(*,*)'n=',n,'it=',it,' interpolated null wind speed at pos (',&
                           x2lon(par(n,pX),par(n,pY)),',', &
                          y2lat(par(n,pY)),')  >> alpha set to zero'
                  alpha=0.0
                else 
                  !write(*,*)n,' Uwind=',(P_Uw*cos(P_angle) - P_Vw*sin(P_angle)),&
                  !            ' Vwind=',(P_Uw*sin(P_angle) + P_Vw*cos(P_angle))
                  private_Average_Value(ID_U_WIND)=private_Average_Value(ID_U_WIND)+ &
                                                (P_Uw*cos(P_angle) - P_Vw*sin(P_angle)) 
                  private_Average_Value(ID_V_WIND)=private_Average_Value(ID_V_WIND)+ &
                                                (P_Uw*sin(P_angle) + P_Vw*cos(P_angle)) 
                  private_Average_Numpart(ID_U_WIND)=private_Average_Numpart(ID_U_WIND)+1   
                  private_Average_Numpart(ID_V_WIND)=private_Average_Numpart(ID_V_WIND)+1   
                endif  


         
              !15deg offset to RHS of wind vector (Coriolis to the
              !right in northern hemisphere)
                  UWindDrift = idt * WindDriftFac * PWind * &
                        cos(alpha - (WindDriftDev *(pi/180.0)))
                  VWindDrift = idt * WindDriftFac * PWind * &
                        sin(alpha - (WindDriftDev *(pi/180.0)))

              if(Stokes)then
                  CALL STOKESDRIFT(P_Uw,P_Vw,P_angle,P_surfdist,  &
                                   UStokesDrift,VStokesDrift,alpha )
                  UStokesDrift = idt * (StokDriftFac/0.016)*UStokesDrift
                  VStokesDrift = idt * (StokDriftFac/0.016)*VStokesDrift
              end if


             
             if ((isnan(UWindDrift) .or. isnan(UStokesDrift) ))then 
               write(*,*)'n=',n,'it=',it
               write(*,*)'Xpar=',par(n,pX),par(n,pnX)
               write(*,*)'Ypar=',par(n,pY),par(n,pnY)
               write(*,*)'Zpar=',par(n,pZ),par(n,pnZ)
               write(*,*)'P_Uw=',P_Uw,' P_Vw=',P_Vw
               write(*,*)'UWindDrift,UStokesDrift=',UWindDrift,UStokesDrift
               write(*,*)'Pwind=',Pwind,' P_depth=',P_depth, &
                         ' alpha=',alpha,' P_angle=',P_angle
               write(*,*)'alpha computed is ',atan( (P_Uw*sin(P_angle) + P_Vw*cos(P_angle))        &
                                    / (P_Uw*cos(P_angle) - P_Vw*sin(P_angle)) )
               write(*,*)'kn1_uw=',kn1_uw,' kn2_uw=',kn2_uw,       &
                        ' kn3_uw=',kn3_uw,' kn4_uw=',kn4_uw
               stop 'found NaN after NewPos'
             endif


          else      !if wind
              CALL setInterp(Xpar,Ypar,n)
          end if    !if wind


          !----------------------------------------------------------
      IF (Behavior.ne.999 .and. Behavior.ne.998 ) THEN
  ! RESET P_zb,P_zc,P_zf based on old Zpar position for further  !SUGGESTION TO SOLVE PROBLEM
  ! use in computing salinity,temperature,behavior,Vturb         !SUGGESTION TO SOLVE PROBLEM
        P_zb = Zpar                                              !SUGGESTION TO SOLVE PROBLEM
        P_zc = Zpar                                              !SUGGESTION TO SOLVE PROBLEM
        P_zf = Zpar                                              !SUGGESTION TO SOLVE PROBLEM
        if (P_zb.GT.P_zetab) P_zb = P_zetab - DBLE(0.00001)      !SUGGESTION TO SOLVE PROBLEM
        if (P_zc.GT.P_zetac) P_zc = P_zetac - DBLE(0.00001)      !SUGGESTION TO SOLVE PROBLEM
        if (P_zf.GT.P_zetaf) P_zf = P_zetaf - DBLE(0.00001)      !SUGGESTION TO SOLVE PROBLEM
      ENDIF

      ! *********************************************************
      ! *                                                       *
      ! *                Salinity and Temperature               *
      ! *                                                       *
      ! *********************************************************


      IF (SaltTempOn) THEN

        !Calculate salinity and temperature at the particle location       
        do i=Fstlev+2,us-2
          if ((Zpar .LT. Pwc_zb(i)) .OR.    &
              (Zpar .LT. Pwc_zc(i)) .OR.    &
              (Zpar .LT. Pwc_zf(i))         ) exit
        enddo
        deplvl = i-2
        if(Zgrid .and. deplvl+3>us)then
            NumInterpLvl=us-Fstlev+1
            deplvl=Fstlev
        else
            NumInterpLvl=4
        endif  
        !P_Salt(n) = WCTS_ITPI("salt",Xpar,Ypar,deplvl,Pwc_zb,Pwc_zc,Pwc_zf,    &
        !                      us,P_zb,P_zc,P_zf,ex,ix,p,4,n,NumInterpLvl)
        Tsalt=WCTS_ITPI(VAR_ID_salt,Xpar,Ypar,deplvl,Pwc_zb,Pwc_zc,Pwc_zf,    &
                              us,P_zb,P_zc,P_zf,ex,ix,p,4,n,NumInterpLvl)        
        Ttemp=WCTS_ITPI(VAR_ID_temp,Xpar,Ypar,deplvl,Pwc_zb,Pwc_zc,Pwc_zf,          &
                              us,P_zb,P_zc,P_zf,ex,ix,p,4,n,NumInterpLvl)
        IF(Behavior.ge.8.and.Behavior.le.10)THEN
          if(Write_Temp_min_max_ins.eq.'max')then
            P_Temp(n) = max(Ttemp,P_Temp(n))
          elseif(Write_Temp_min_max_ins.eq.'min')then
            P_Temp(n) = min(Ttemp,P_Temp(n))
          elseif(Write_Temp_min_max_ins.eq.'ins')then
            P_Temp(n) = Ttemp 
          else
            write(*,*)'Write_Temp_min_max_ins not specified'
            stop
          endif

          if(Write_Salt_min_max_ins.eq.'max')then
            P_Salt(n) = max(Tsalt,P_Salt(n))
          elseif(Write_Salt_min_max_ins.eq.'min')then
            P_Salt(n) = min(Tsalt,P_Salt(n))
          elseif(Write_Salt_min_max_ins.eq.'ins')then
            P_Salt(n) = Tsalt 
          else
            write(*,*)'Write_Salt_min_max_ins not specified'
            stop
          endif

        ELSE
          P_Salt(n) = Tsalt
          P_Temp(n) = Ttemp 
        ENDIF
        private_Average_Value(ID_TEMP)=private_Average_Value(ID_TEMP)+Ttemp
        private_Average_Numpart(ID_TEMP)=private_Average_Numpart(ID_TEMP)+1
      ENDIF

      ! *********************************************************
      ! *                                                       *
      ! *                  Horizontal Turbulence                *
      ! *                                                       *
      ! *********************************************************

      !$OMP MASTER
      IF(WriteModelTiming) call CPU_TIME(times(4))
      !$OMP END MASTER

        IF (HTurbOn) CALL HTurb(TurbHx,TurbHy)


      ! *********************************************************
      ! *                                                       *
      ! *                   Verticle Turbulence                 *
      ! *                                                       *
      ! ********************************************************* 


      !$OMP MASTER
      IF(WriteModelTiming) call CPU_TIME(times(5))
      !$OMP END MASTER

      IF (VTurbOn) CALL VTurb(P_zc,P_depth,P_zetac,p,ex,ix,Pwc_wzb,Pwc_wzc,    &
                              Pwc_wzf,Xpar,Ypar,TurbV)

  
      ! *********************************************************
      ! *                                                       *
      ! *                       Behavior                        *
      ! *                                                       *
      ! *********************************************************


      !$OMP MASTER
      IF(WriteModelTiming) call CPU_TIME(times(6))
      !$OMP END MASTER
      
      IF((Behavior.ge.8.and.Behavior.le.11).and.read_GrainSize)then
       P_GrainSize(n)=max(P_GrainSize(n),getInterp(Xpar,Ypar,VAR_ID_GrainSize,klev))
       !if(readNetcdfSwdown)P_swdown = getInterp(Xpar,Ypar,VAR_ID_swdown,klev)
      ENDIF
      PTemptmp=0.0
      if(SaltTempOn)PTemptmp=P_Temp(n)
     
      IF (Behavior.NE.0 .AND. Behavior.LE.1000) CALL behave(Xpar,Ypar,Zpar,Pwc_zb,Pwc_zc,Pwc_zf,      &
           P_zb,P_zc,P_zf,P_zetac,par(n,pAge),P_depth,P_U,P_V,P_angle,P_Temp(n),   &
           n,it,ex,ix,ix(3)/DBLE(86400),p,bott,XBehav,YBehav,ZBehav,LarvSize)!,P_swdown)
           IF(Behavior.ge.8.and.Behavior.le.11) P_Size(n)=LarvSize

      !--- CL-OGS: implemented then cancelled by OILTRANS RK scheme
!     IF(Wind)then
!       !**********   UWind  ***************************** 
!       ey(1) = getInterp("uwindb",klev)
!       ey(2) = getInterp("uwindc",klev)
!       ey(3) = getInterp("uwindf",klev)
!       P_Uwind = polintd(ex,ey,3,ix(2))
!       !**********   VWind  ***************************** 
!       ey(1) = getInterp("vwindb",klev)
!       ey(2) = getInterp("vwindc",klev)
!       ey(3) = getInterp("vwindf",klev)
!       P_Vwind = polintd(ex,ey,3,ix(2))
!
!       if(Behavior.EQ.1000)then
!         !AdvecUwind=idt*(P_Uwind*cos(P_angle) - P_Vwind*sin(P_angle))*(10.0-min(P_zetac-par(n,pZ),10.0))/10.0
!         !AdvecVwind=idt*(P_Uwind*sin(P_angle) + P_Vwind*cos(P_angle))*(10.0-min(P_zetac-par(n,pZ),10.0))/10.0
!         AdvecUwind=idt*(P_Uwind*cos(P_angle) - P_Vwind*sin(P_angle)) * 0.01
!         AdvecVwind=idt*(P_Uwind*sin(P_angle) + P_Vwind*cos(P_angle)) * 0.01
!         !if(n==11) write(*,*)'depth below zeta=',P_zetac-par(n,pZ),' percent of wind=',(10.0-min(P_zetac-par(n,pZ),10.0))/10.0
!       endif
!     ENDIF


  
      ! *********************************************************
      ! *                                                       *
      ! *     Update Particle Locations and Check Boundaries    *
      ! *                                                       *
      ! *********************************************************

      !$OMP MASTER
      IF(WriteModelTiming) call CPU_TIME(times(7))
      !$OMP END MASTER

      newXpos = 0.0
      newYpos = 0.0
      newZpos = 0.0

      !--- CL-OGS: implemented then cancelled by OILTRANS RK scheme
!     !Update due to Advection Turbulence and wind
!     newXpos = par(n,pX) + AdvectX + TurbHx +AdvecUwind
!     newYpos = par(n,pY) + AdvectY + TurbHy +AdvecVwind
!     newZpos = par(n,pZ) + AdvectZ + TurbV

      !Update due to Advection and Turbulence
      newXpos = par(n,pX) + AdvectX + TurbHx
      newYpos = par(n,pY) + AdvectY + TurbHy
      newZpos = par(n,pZ) + AdvectZ + TurbV
      !write(*,*)it,n,AdvectZ,kn1_u,kn2_u,kn3_u,kn4_u,kn1_v,kn2_v,kn3_v,kn4_v,kn1_w,kn2_w,kn3_w,kn4_w
           
  !IF (Behavior.eq.8 .and.( P_Size(n)>14 &
      !     .and.abs(par(n,pZ)-P_depth)<2)) then
      !   newXpos = par(n,pX) + AdvectX*0.1  + TurbHx
      !   newYpos = par(n,pY) + AdvectY*0.1  + TurbHy
      !ENDIF

         if(Wind)then
          newXpos = newXpos + UWindDrift
          newYpos = newYpos + VWindDrift
          if(Stokes)then
            newXpos = newXpos + UStokesDrift
            newYpos = newYpos + VStokesDrift
          end if
         end if !if(Wind)


      !Horizontal boundary tests. Ensure particle still within domain
      !If not, reflect particle off main boundary or island walls (at old klev
      !for Zgrid)
      Xpos = par(n,pX)
      Ypos = par(n,pY)
      nXpos = newXpos
      nYpos = newYpos

      skipbound = -1
      reflects = 0
      waterFlag = .FALSE.
      saveintersectf=0
      do
        call intersect_reflect(Xpos,Ypos,nXpos,nYpos,fintersectX,fintersectY,  &
         freflectX,freflectY,intersectf,skipbound,isWater,klev,bndintersect) 
        saveintersectf=max(saveintersectf,intersectf)
        if(intersectf == 0)exit

       IF(TrackCollisions) hitLand(n) = hitLand(n) + 1
       
       coastdist=0

       IF(OilOn .and. StrandingDist>=0)then
        P_coastdist(n)=0
        nXpos = Xpos+(fintersectX-Xpos)*0.9 
        nYpos = Ypos+(fintersectY-Ypos)*0.9 
        par(n,pnZ) = newZpos
        par(n,pStatus) = 2
        call p_Stranding(n) ! Apply stranding in "is_Stranded"
        private_Average_Numpart(ID_STRANDED) =private_Average_Numpart(ID_STRANDED) + 1
        exit
       ELSE
        if(OpenOceanBoundary .AND. isWater)then
          par(n,pnX) = Xpos+(fintersectX-Xpos)*0.9  !fintersectX
          par(n,pnY) = Ypos+(fintersectY-Ypos)*0.9  !fintersectY
          par(n,pnZ) = newZpos
          call setOut(n)
          waterFlag = .TRUE.
          exit
        endif

        reflects = reflects + 1
        if(BndOut)then
          write(*,*)'p',p,' it',it,' part ',n,                     &
        ' got reflected at lower level',      &
        ' while going from pos ',x2lon(Xpos,Ypos),',',y2lat(Ypos),             &
        ' to pos ',x2lon(nXpos,nYpos),',',y2lat(nYpos)
          length=sqrt((x2lon(nXpos,nYpos)-x2lon(Xpos,Ypos))**2+                &
                      (y2lat(nYpos)-y2lat(Ypos))**2)
         write(annotate,'(a,i6)')'I',it+p*int(dt/idt)
         call writeErrortopython(annotate,'c',Xpos,Ypos,nXpos,nYpos) 
        endif
        if(reflects > 3) then
          call handleERROR('intersect ',1,reflects,n,     &
                      Xpos,Ypos,par(n,pZ), &
                      nXpos,nYpos,par(n,pZ),docycle)
          waterFlag = .TRUE.
          exit
          if(docycle) return
        endif

        Xpos = fintersectX
        Ypos = fintersectY
        nXpos = freflectX
        nYpos = freflectY

       ENDIF
      enddo
      
      if(waterFlag) return

      newXpos = nXpos
      newYpos = nYpos

      !Check vertical boundaries and reflect
      !  if particle above surface, particle reflects off surface
      reflectsup=0.0
      if (newZpos.GT.P_zetac) then
        reflectsup = P_zetac - newZpos
        NewZpos = P_zetac + reflectsup
      endif

!      newZpos = par(n,pZ) + AdvectZ + TurbV


      IF(Zgrid)then
       CALL setEle(newXpos,newYpos,max(par(n,pZ),newZpos),n,it,2,ele_err)
       nklev=getKRlevel(newZpos)
       if(ele_err.ne.0)then
        IF(OilOn)then
          do posfactor=8,0,-1
            newXpos = Xpos+(fintersectX-Xpos)*posfactor/10.0 
            newYpos = Ypos+(fintersectY-Ypos)*posfactor/10.0
            CALL setEle(newXpos,newYpos,max(par(n,pZ),newZpos),n,it,2,ele_err)
            if(ele_err.eq.0)then
               write(*,*)'Stranded particle ',n,' put at ',posfactor*10., &
               '% of the distance from the boundary for setEle reasons'
               exit
            endif
          enddo
        ENDIF
        if(ele_err.ne.0)then
         !  if(nklev<klev)then
         !   CALL setEle(newXpos,newYpos,par(n,pZ),n,3,ele_err)
         !   if(ele_err.ne.0)then
         !   endif
         !  else
         !    write(*,*)'setEle error at newZpos=',newZpos,&
         !     ' for n, it= ',n,it
         !    write(*,*)'but same level for oldZpos=',par(n,pZ), &
         !                  ', killing particle ',n,it
         !  endif
         !endif
         !If(ele_err .ne.  0)then 
         ! write(*,*) "pb setEle UNSOLVED !!! "
              write(*,'(a,i8,a,i2,a,3f8.2,2(a,f7.2,a,i3,a,f7.2),(a,i5),10(a,f8.2),a)') &
                  'ERROR it ',it,' on 2nd setEle call, error number is ',ele_err, &
                 '. Local old, present and next depth are ',&
                    par(n,pZ),P_depth,nP_depth, &
                 ' at old Z=',par(n,pZ),' ZW(',getKRlevel(par(n,pZ)),')=', &
                            Pwc_wzc(getKRlevel(par(n,pZ))),  &
                  ' at new Z=',newZpos  ,' ZW(',getKRlevel(newZpos),')=', &
                            Pwc_wzc(getKRlevel(newZpos)) ,  &
                  '; PART n=',n,' WENT DOWN of ',-(newZpos-par(n,pZ)), &
                  'meters : newZ=',newZpos,              & 
                  ' =min(depth[',P_depth, ' ], oldZ[',par(n,pZ),            &
                   ' ] + Adv[', AdvectZ,' ] + Tu[',TurbV,' ] + Bhv[',ZBehav,       &
                   ' ]) + reflectsup[',reflectsup,' ] + reflectinf[',reflectinf,&
                   ' ] + pushedup[',pushedup,']'
           call handleERROR('setEle    ',2,ele_err,n,     &
                        par(n,pX),par(n,pY),par(n,pZ), &
                        newXpos,newYpos,par(n,pZ),docycle  )
           if(docycle) return
        endif
      Endif !ele_err

       conflict_tmp=1
       call getDepth(newXpos,newYpos,n,it,nP_depth,Fstlev,conflict_tmp)
        If(Fstlev>us)then                                    
         if(OilOn)then
          do posfactor=8,0,-1
            newXpos = Xpos+(fintersectX-Xpos)*posfactor/10.0 
            newYpos = Ypos+(fintersectY-Ypos)*posfactor/10.0
            CALL setEle(newXpos,newYpos,max(par(n,pZ),newZpos),n,it,3,ele_err)
            If(ele_err .ne.  0)then 
              write(*,'(a,i12,a,i2,a,3f7.2,6(a,f7.2),(a,i5),10(a,f6.2),a)') &
                  'ERROR it',it,'on 3rd setEle call, error number is ',ele_err, &
                 'local old, present and next depth is ',&
                    par(n,pZ),P_depth,nP_depth, &
                 'at old Z=',par(n,pZ),' ZW(',getKRlevel(par(n,pZ)),')=', &
                            Pwc_wzc(getKRlevel(par(n,pZ))),  &
                  'at new Z=',newZpos  ,' ZW(',getKRlevel(newZpos),')=', &
                            Pwc_wzc(getKRlevel(newZpos)) ,  &
                  'PART n=',n,' WENT DOWN of ',-(newZpos-par(n,pZ)), &
                  'm : newZ=',newZpos,              & 
                  '=min(depth[',P_depth, '], oldZ[',par(n,pZ),            &
                   ']+Adv[', AdvectZ,']+Tu[',TurbV,']+Bhv[',ZBehav,       &
                   ']) +reflectsup[',reflectsup,']+ reflectinf[',reflectinf,&
                   ']+ pushedup[',pushedup,']'
           
              call handleERROR('setEle    ',3,ele_err,n,     &
                           par(n,pX),par(n,pY),par(n,pZ), &
                           newXpos,newYpos,newZpos,docycle  )
              write(*,*)'-------------------------------'
              if(docycle) return
            Endif !ele_err
            call getDepth(newXpos,newYpos,n,it,nP_depth,Fstlev,conflict_tmp)
            if(Fstlev<=us)then
               write(*,*)'Stranded particle ',n,' put at ',posfactor*10, &
               '% of the distance from the boundary for depth reasons'
               exit
            endif
          enddo
         endif
         if(Fstlev>us)then                                    
            call handleERROR('Fstlev    ',2,conflict_tmp,n,     &
                     par(n,pX),par(n,pY),par(n,pZ), &
                     newXpos,newYpos,newZpos,docycle)
            if(docycle) return 
         endif
        Endif
        If(conflict_tmp.ne.1)then
          call handleERROR('dDepth    ',2,conflict_tmp,n,     &
                     par(n,pX),par(n,pY),par(n,pZ), &
                     newXpos,newYpos,newZpos,docycle)
          call die(n)
           if(docycle) return 
        Else
              P_depth=nP_depth
        Endif
      ENDIF !Zgrid
      private_Average_Value(ID_DEPTH)=private_Average_Value(ID_DEPTH)+P_depth
      private_Average_Numpart(ID_DEPTH)=private_Average_Numpart(ID_DEPTH)+1
      if (newZpos.GT.P_zetac) then
        reflectsup = P_zetac - newZpos
        NewZpos = P_zetac + reflectsup
      endif

      pushedup=0.0
      if (newZpos.LT.P_depth+z0) then
        pushedup= P_depth - newZpos + z0+1e-5
        NewZpos = P_depth + z0 +1e-5
        IF(TrackCollisions) hitBottom(n) = hitBottom(n) + 1
      endif 


      !Update due to Behavior
      newZpos = NewZpos + ZBehav


      !  if particle deeper than bottom, particle reflects off bottom
      reflectinf=0.0
      if (newZpos.LT.P_depth) then
        reflectinf = P_depth - newZpos
        NewZpos = P_depth + reflectinf
        IF(TrackCollisions) hitBottom(n) = hitBottom(n) + 1
      endif 

      if(Behavior == 7)then

        if(bott)then
          !Particle on bottom and doesn't move
          newXpos = par(n,pX)
          newYpos = par(n,pY)
          newZpos = P_depth
        else
          !Update due to Behavior 7
          !Set Z position to swimdepth
          newXpos = newXpos + XBehav      !ewn, FO
          newYpos = newYpos + YBehav      !ewn, FO
          newZpos = P_depth + Swimdepth
        endif

      endif
  
      !Check vertical boundaries and move back in domain
      !  if particle above surface, particle moves just below surface
      if (newZpos.GT.P_zetac) NewZpos = P_zetac - DBLE(0.000001)

      !  if particle deeper than bottom, particle moves just off bottom
      if (newZpos.LT.P_depth) then
        NewZpos = P_depth + DBLE(0.000001)
        IF(TrackCollisions) hitBottom(n) = hitBottom(n) + 1
      endif

      if(Behavior.eq.997 .and.NewZpos .GT. minpartdepth)then
        PartAtSurf(n,1)= x2lon(par(n,pX),par(n,pY))
        PartAtSurf(n,2)= y2lat(par(n,pY))
        PartAtSurf(n,3)= par(n,pZ)
        PartAtSurf(n,4)= ix(3)
        call die(n)
       return
      endif
      !if(NewZpos .GT. minpartdepth) NewZpos = minpartdepth - DBLE(0.000001)
      !if(NewZpos .LT. maxpartdepth) NewZpos = maxpartdepth + DBLE(0.000001)

   !  if (n==5 .and. abs(AdvectX)>0.000001) &
   ! !    write(*,*)'p',p,' it',it,' part ',n,          &
   ! !    '; old pos = (',             &
   ! !    x2lon(par(n,pX),par(n,pY)),',',     &
   ! !    y2lat(par(n,pY)),par(n,pZ),         &
   ! !  ') ; new pos =(',x2lon(newXpos,newYpos),',',y2lat(newYpos),&
   ! ! newZpos, &
   ! !  ') AdvectX, Y =',AdvectX,AdvectY,' P_depth=',P_depth
   !       write(*,'(3(a,i5),17(a,f11.6))')       &
   !       'p',p,' it',it,' part n=',n,' old pos = (',             &
   !      x2lon(par(n,pX),par(n,pY)),',',     &
   !      y2lat(par(n,pY)),',',par(n,pZ),         &
   !    ') ; new pos =(',x2lon(newXpos,newYpos),',',y2lat(newYpos),&
   !   ',',newZpos,' WENT UP of ',newZpos-par(n,pZ), &
   !       'm : newZ=',newZpos,              & 
   !       'm =min(depth[',P_depth, '], oldZ[',par(n,pZ),            &
   !        ']+Adv[', AdvectZ,']+Tu[',TurbV,']+Bhv[',ZBehav,       &
   !        ']) +reflectsup[',reflectsup,']+ reflectinf[',reflectinf,&
   !        ']+ pushedup[',pushedup,'], m above bott=',newZpos-P_depth

      IF(Zgrid)then
         nklev=getKRlevel(newZpos)
         P_oldLev(n)=nklev
        !if(newZpos-par(n,pZ)>3.0)then
        !   write(*,'((a,i5),10(a,f8.2),a)')       &
        !   'WARNING n=',n,' WENT UP of ',newZpos-par(n,pZ), &
        !   'm : newZ=',newZpos,              & 
        !   '=min(depth[',P_depth, '], oldZ[',par(n,pZ),            &
        !    ']+Adv[', AdvectZ,']+Tu[',TurbV,']+Bhv[',ZBehav,       &
        !    ']) +reflectsup[',reflectsup,']+ reflectinf[',reflectinf,&
        !    ']+ pushedup[',pushedup,']'
        !endif
      ELSE
        nklev=1
      ENDIF
      
      ! Check to make sure new position is within a rho, u and v element
      CALL setEle(newXpos,newYpos,newZpos,n,it,4,ele_err)
      If(ele_err .ne.  0)then 
        write(*,'(a,i12,a,i2,a,3f7.2,6(a,f7.2),(a,i5),10(a,f6.2),a)') &
            'ERROR it',it,'on 4th setEle call, error number is ',ele_err, &
           'local old, present and next depth is ',&
              par(n,pZ),P_depth,nP_depth, &
           'at old Z=',par(n,pZ),' ZW(',getKRlevel(par(n,pZ)),')=', &
                      Pwc_wzc(getKRlevel(par(n,pZ))),  &
            'at new Z=',newZpos  ,' ZW(',getKRlevel(newZpos),')=', &
                      Pwc_wzc(getKRlevel(newZpos)) ,  &
            'PART n=',n,' WENT DOWN of ',-(newZpos-par(n,pZ)), &
            'm : newZ=',newZpos,              & 
            '=min(depth[',P_depth, '], oldZ[',par(n,pZ),            &
             ']+Adv[', AdvectZ,']+Tu[',TurbV,']+Bhv[',ZBehav,       &
             ']) +reflectsup[',reflectsup,']+ reflectinf[',reflectinf,&
             ']+ pushedup[',pushedup,']'

        call handleERROR('setEle    ',4,ele_err,n,     &
                     par(n,pX),par(n,pY),par(n,pZ), &
                     newXpos,newYpos,newZpos,docycle  )
        write(*,*)'-------------------------------'
        if(docycle) return
      Endif !ele_err

      if(ele_err.ne.0)then
      endif

      !Check to make sure new position is within model boundaries
      call mbounds(newYpos,newXpos,inbounds,mainbound,nklev,m_nestdeg)
      if(inbounds /= 1) then
        write(*,*) 'ERROR: Particle',n,' Outside Main Boundaries After ',         &
                     'intersect_reflect ',ErrorFlag
        call handleERROR('mbounds   ',1,ele_err,n,     &
                     par(n,pX),par(n,pY),par(n,pZ), &
                     newXpos,newYpos,newZpos,docycle  )
        if(docycle) return
      endif

       if ((isnan(par(n,pX)) .or. isnan(par(n,pY)) )          &
        .or.  (  isnan(par(n,pnX)).or.isnan(par(n,pnY)) ) )then
         write(*,*)'n=',n,'it=',it
         write(*,*)'Xpar=',par(n,pX),par(n,pnX)
         write(*,*)'Ypar=',par(n,pY),par(n,pnY)
         write(*,*)'Zpar=',par(n,pZ),par(n,pnZ)
         stop 'found NaN after mbounds'
       endif

      !Check to make sure new position is not within an island
      call ibounds(in_island,newYpos,newXpos,island,nklev,i_nestdeg)
      if(in_island == 1) then
       if(i_nestdeg<m_nestdeg)then
        write(*,*) 'WARNING: Particle',n,' Inside Island Boundaries nestdeg=',&
                i_nestdeg ,' P_depth=',P_depth,FstLev
        write(*,*)'p',p,' it',it,' part ',n,                     &
        ' got reflected at lower level',      &
        ' while going from pos ',             &
          x2lon(par(n,pX),par(n,pY)),',',     &
          y2lat(par(n,pY)),                   &
        ' to pos ',x2lon(newXpos,newYpos),',',y2lat(newYpos)
        write(*,*)'local old, present and next depth is ',&
              par(n,pZ),P_depth,nP_depth
        write(*,*)'at old Z=',par(n,pZ),' ZW(',getKRlevel(par(n,pZ)),')=', &
                      Pwc_wzc(getKRlevel(par(n,pZ)))
        write(*,*)'at new Z=',newZpos  ,' ZW(',getKRlevel(newZpos),')=', &
                      Pwc_wzc(getKRlevel(newZpos))
       else        
        write(*,*) 'ERROR: Particle',n,' Inside Island Boundaries After ',     &
                     'intersect_reflect ',ErrorFlag,' P_depth=',P_depth,FstLev
          write(*,*)'p',p,' it',it,' part ',n,                     &
        ' got reflected at lower level',      &
        ' while going from pos ',             &
          x2lon(par(n,pX),par(n,pY)),',',     &
          y2lat(par(n,pY)),                   &
        ' to pos ',x2lon(newXpos,newYpos),',',y2lat(newYpos)
        write(*,*)'local old, present and next depth is ',&
              par(n,pZ),P_depth,nP_depth
        write(*,*)'at old Z=',par(n,pZ),' ZW(',getKRlevel(par(n,pZ)),')=', &
                      Pwc_wzc(getKRlevel(par(n,pZ)))
        write(*,*)'at new Z=',newZpos  ,' ZW(',getKRlevel(newZpos),')=', &
                      Pwc_wzc(getKRlevel(newZpos))
        call handleERROR('ibounds   ',1,ele_err,n,     &
                     par(n,pX),par(n,pY),par(n,pZ), &
                     newXpos,newYpos,newZpos,docycle  )
        if(docycle) return
       endif
      endif


       if ((isnan(par(n,pX)) .or. isnan(par(n,pY)) )          &
        .or.  (  isnan(par(n,pnX)).or.isnan(par(n,pnY)) ) )then
         write(*,*)'n=',n,'it=',it
         write(*,*)'Xpar=',par(n,pX),par(n,pnX)
         write(*,*)'Ypar=',par(n,pY),par(n,pnY)
         write(*,*)'Zpar=',par(n,pZ),par(n,pnZ)
         stop 'found NaN at end BC tests'
       endif

      ! End boundary condition tests ******************* 

      !Assign new particle positions
      par(n,pnX) = newXpos
      par(n,pnY) = newYpos
      par(n,pnZ) = newZpos
       if ((isnan(par(n,pX)) .or. isnan(par(n,pY)) )          &
        .or.  (  isnan(par(n,pnX)).or.isnan(par(n,pnY)) ) )then
         write(*,*)'n=',n,'it=',it
         write(*,*)'Xpar=',par(n,pX),par(n,pnX)
         write(*,*)'Ypar=',par(n,pY),par(n,pnY)
         write(*,*)'Zpar=',par(n,pZ),par(n,pnZ)
         stop 'found NaN before writecurrents'
       endif

      if(WriteCurrents)then
         write(fstring,'(a,i5.5,a)')'Currents',n,''//TRIM(NCOutFile)//'.csv'
         OPEN(110,FILE=trim(fstring),POSITION='APPEND')
         write(110,"(F8.3,'  ,   ',3(F10.5,','),'   ',6(F10.6,','),'    ',i5)")(DBLE(ix(3))/DBLE(86400)), & 
                  x2lon(par(n,pnX),par(n,pnY)),y2lat(par(n,pnY)),par(n,pnZ),  &
                  AdvectX/float(idt), &
                  AdvectY/float(idt), &
                  AdvectZ/float(idt) , & 
                  TurbHx/float(idt),  &
                  TurbHy/float(idt),  &
                  TurbV/float(idt),reflects
         CLOSE(110)


      endif
     ! *                      Settlement                       *
     ! BEUG:moved up to avoid computing newpos if n already settling at currentpos



!--- CL-OGS: SETTLEMENT SECTION !------------------------------------------------------
!--- CL-OGS: depending on value of storedincolor (defined in LTRANS.data)
!--- CL-OGS: the status (and in output the color variable) can contain the water element number
!--- CL-OGS:       or the First level number


      ! *********************************************************
      ! *                                                       *
      ! *                      Settlement                       *
      ! * !--- CL-OGS:BEUG?: moved here as check settlement at current pos *
      ! *********************************************************


     ! If stranding is active or Write_coastdist requested, compute coast dist
     ! at upper level 
      if( (StrandingDist>0) .or. Write_coastdist )  then
          call Get_coastdist(newYpos,newXpos,us,P_coastdist(n))
          coastdist=P_coastdist(n)
      endif


      if(settlementon) then
        CALL testSettlement(par(n,pAge),n,par(n,pX),par(n,pY),par(n,pZ),inpoly, &
                            P_depth,klev,saveintersectf,coastdist)
        if(Write_Poly_Presence .and. inpoly.gt.0)then
          Time_in_Poly(inpoly-poly0+1,n)=Time_in_Poly(inpoly-poly0+1,n)+int(idt)
        endif
      endif

      if(settlementon) then
        !if(inpoly>0.and.isStranded(n))write(*,*)'it=',it,' n=',n,' is in poly',inpoly
        if (inpoly .GT. 0 .and. (isSettled(n) .or. isStranded(n)) ) then
          if(isStranded(n)) then
            par(n,pnZ) = par(n,pZ) 
          else
            par(n,pnZ) = P_depth
          endif
          endpoly(n) = inpoly
          par(n,pLifespan) = par(n,pAge)
          ! When n is settled or stranded, color(n)= -99999:
          ! W hen n is settled or stranded, color(n)= -1*polynumber:
          if(.not.OilOn ) then 
            par(n,pStatus) = - inpoly
          else
            par(n,pStatus) = 2
          endif
          return
        endif
      endif
!--- CL-OGS: END OF SETTLEMENT SECTION !-----------------------------------------------
  
      private_Average_Numpart(ID_ALIVE)=private_Average_Numpart(ID_ALIVE)+1

      ! *****************************************************************
      ! *                      End of Particle Loop                     *
      ! *****************************************************************

      !$OMP MASTER
      IF(WriteModelTiming) then
        call CPU_TIME(times(8))

        timeCounts(2) = timeCounts(2) + (times(3)-times(2))
        timeCounts(3) = timeCounts(3) + (times(4)-times(3))
        timeCounts(4) = timeCounts(4) + (times(5)-times(4))
        timeCounts(5) = timeCounts(5) + (times(6)-times(5))
        timeCounts(6) = timeCounts(6) + (times(7)-times(6))
        timeCounts(7) = timeCounts(7) + (times(8)-times(7))
      ENDIF
      !$OMP END MASTER
  end subroutine update_single_particle

  SUBROUTINE find_currents(Xpar,Ypar,Zpar,Pwc_zb,Pwc_zc,Pwc_zf,Pwc_wzb,   &
    Pwc_wzc,Pwc_wzf,P_zb,P_zc,P_zf,ex,ix,p,version,Uad,Vad,Wad,&
    kwBotLvl,P_depth,n)
    !This Subroutine calculates advection currents at the particle's 
    !  location in space and time

    USE PARAM_MOD,  ONLY: us,ws,z0,idt,VInterpUVinSurfWater, &
     Zgrid_depthinterp,BottomLayerThickness,PercentVelinBottomLayer, &
     PercentVel_under_z0
    USE HYDRO_MOD,  ONLY: WCTS_ITPI,setInterp,getInterp
    USE TENSION_MOD, ONLY: TSPSI,HVAL
    USE INT_MOD,    ONLY: linint,polintd
#include "VAR_IDs.h"
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: p,version,kwBotLvl,n
    DOUBLE PRECISION, INTENT(IN) :: Xpar,Ypar,Zpar,P_zb,P_zc,P_zf,ex(3),ix(3), &
      Pwc_zb(:),Pwc_zc(:),Pwc_zf(:),Pwc_wzb(:),Pwc_wzc(:),Pwc_wzf(:),P_depth
    DOUBLE PRECISION, INTENT(OUT) :: Uad,Vad,Wad

    !Number of Depth Levels to create tension spline with
    INTEGER, PARAMETER :: nN = 4

    INTEGER :: i,ii,iii,NumInterpLvlii,NumInterpLvliii
    DOUBLE PRECISION :: P_Ub,P_Uc,P_Uf,P_Vb,P_Vc,P_Vf,P_Wb,P_Wc,P_Wf,ey(3),    &
      Pwc_ub,Pwc_uc,Pwc_uf,Pwc_vb,Pwc_vc,Pwc_vf,Pwc_wb,Pwc_wc,Pwc_wf,slope, &
      BottDepth, BottLayerHeight 

        !version: 1 = return b, 2 = return c, 3 = return f

    !Set Interpolation Values for the current particle
    CALL setInterp(Xpar,Ypar,n)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !       Determine the Lowest Numbered US-Level of the Closest Four
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    do i = kwBotLvl+2 , us-2
      if ( (Zpar .LT. Pwc_zb(i)) .OR.  &
           (Zpar .LT. Pwc_zc(i)) .OR.  &
           (Zpar .LT. Pwc_zf(i))       )then
                exit
      endif
    enddo
    ii = i-2

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !       Determine the Lowest Numbered WS-Level of the Closest Four
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    do i=kwBotLvl+2,ws-2
      if ( (Zpar .LT. Pwc_wzb(i)) .OR. &
           (Zpar .LT. Pwc_wzc(i)) .OR. &
           (Zpar .LT. Pwc_wzf(i))      )then
                exit
      endif 
    enddo
    iii = i - 2

    if(ii+3>us)then
        ii=max(KwBotLvl,us-3)
        NumInterpLvlii=us-ii+1
    else
        NumInterpLvlii=nN
    endif
    if(iii+3>ws)then
        iii=max(KwBotLvl,ws-3)
        NumInterpLvliii=ws-iii+1
    else
        NumInterpLvliii=nN
    endif

    !           *********************************************************
    !           *                                                       *
    !           *       Calculate U,V,W in Water Column Profile         *
    !           *                                                       *
    !           *********************************************************

    if(Zgrid_depthinterp)then
      BottDepth=P_depth
      if(BottomLayerThickness>0)then
        BottLayerHeight=P_depth+BottomLayerThickness
      else
        BottLayerHeight=P_depth+abs(Pwc_zc(us)-Pwc_wzc(us))
      endif
    else
      BottDepth=   &
      max(max(Pwc_wzb(kwBotLvl),Pwc_wzc(kwBotLvl)),Pwc_wzf(kwBotLvl))
      BottLayerHeight=   &
      max(max(Pwc_zb(kwBotLvl),Pwc_zc(kwBotLvl)),Pwc_zf(kwBotLvl))
    endif

    !i. Determine if particle is deep enough that velocities are affected by 
    !  the bottom.
    ! If so, apply log layer between deepest current velocity predicitons 
    ! (deepest rho s-level for u,v and deepest w s-level for w) and bottom.
    ! OR, if below z0, set advection velocities to 0.0
   !if( (Zpar .LT. BottDepth+z0)     ) then

   !  Uad = 0.0
   !  Vad = 0.0
   !  Wad = 0.0

   !elseif (Zpar .LT. BottLayerHeight  ) then ! apply log layer
    if (Zpar .LT. BottLayerHeight  ) then ! apply log layer

      Pwc_Ub = getInterp(Xpar,Ypar,VAR_ID_uvelb,kwBotLvl  )
      Pwc_Uc = getInterp(Xpar,Ypar,VAR_ID_uvelc,kwBotLvl  )
      Pwc_Uf = getInterp(Xpar,Ypar,VAR_ID_uvelf,kwBotLvl  )
      Pwc_Vb = getInterp(Xpar,Ypar,VAR_ID_vvelb,kwBotLvl  )
      Pwc_Vc = getInterp(Xpar,Ypar,VAR_ID_vvelc,kwBotLvl  )
      Pwc_Vf = getInterp(Xpar,Ypar,VAR_ID_vvelf,kwBotLvl  )
      Pwc_Wb = getInterp(Xpar,Ypar,VAR_ID_wvelb,kwBotLvl+1)
      Pwc_Wc = getInterp(Xpar,Ypar,VAR_ID_wvelc,kwBotLvl+1)
      Pwc_Wf = getInterp(Xpar,Ypar,VAR_ID_wvelf,kwBotLvl+1)

      !  u(z)  = [ u(zB) / (log(zB/zo) ] * (log (z/zo) 
      !where:
      !  u is current velocity
      !  zB is height of first rho-water level above bottom
      !  z0 is roughness height of model
      !  z is height of desired velocity
      !
      !  Note that Pwc_wzb(kwBotLvl) = P_depth = Depth at particle location
      if(PercentVelinBottomLayer .le.0)then
        P_Ub=Pwc_Ub*log10((Zpar-BottDepth)/z0)/log10((BottLayerHeight-BottDepth)/z0)
        P_Uc=Pwc_Uc*log10((Zpar-BottDepth)/z0)/log10((BottLayerHeight-BottDepth)/z0)
        P_Uf=Pwc_Uf*log10((Zpar-BottDepth)/z0)/log10((BottLayerHeight-BottDepth)/z0)
        P_Vb=Pwc_Vb*log10((Zpar-BottDepth)/z0)/log10((BottLayerHeight-BottDepth)/z0)
        P_Vc=Pwc_Vc*log10((Zpar-BottDepth)/z0)/log10((BottLayerHeight-BottDepth)/z0)
        P_Vf=Pwc_Vf*log10((Zpar-BottDepth)/z0)/log10((BottLayerHeight-BottDepth)/z0)
        P_Wb=Pwc_Wb*log10((Zpar-BottDepth)/z0)/log10((BottLayerHeight-BottDepth)/z0)
        P_Wc=Pwc_Wc*log10((Zpar-BottDepth)/z0)/log10((BottLayerHeight-BottDepth)/z0)
        P_Wf=Pwc_Wf*log10((Zpar-BottDepth)/z0)/log10((BottLayerHeight-BottDepth)/z0)
      else
        P_Ub=Pwc_Ub*PercentVelinBottomLayer
        P_Uc=Pwc_Uc*PercentVelinBottomLayer
        P_Uf=Pwc_Uf*PercentVelinBottomLayer
        P_Vb=Pwc_Vb*PercentVelinBottomLayer
        P_Vc=Pwc_Vc*PercentVelinBottomLayer
        P_Vf=Pwc_Vf*PercentVelinBottomLayer
        P_Wb=Pwc_Wb*PercentVelinBottomLayer
        P_Wc=Pwc_Wc*PercentVelinBottomLayer
        P_Wf=Pwc_Wf*PercentVelinBottomLayer
      endif

      if(Zpar .LT. BottDepth+z0 .or. P_Ub/Pwc_Ub< PercentVel_under_z0) P_Ub=PercentVel_under_z0*Pwc_Ub
      if(Zpar .LT. BottDepth+z0 .or. P_Uc/Pwc_Uc< PercentVel_under_z0) P_Uc=PercentVel_under_z0*Pwc_Uc
      if(Zpar .LT. BottDepth+z0 .or. P_Uf/Pwc_Uf< PercentVel_under_z0) P_Uf=PercentVel_under_z0*Pwc_Uf
      if(Zpar .LT. BottDepth+z0 .or. P_Vb/Pwc_Vb< PercentVel_under_z0) P_Vb=PercentVel_under_z0*Pwc_Vb
      if(Zpar .LT. BottDepth+z0 .or. P_Vc/Pwc_Vc< PercentVel_under_z0) P_Vc=PercentVel_under_z0*Pwc_Vc
      if(Zpar .LT. BottDepth+z0 .or. P_Vf/Pwc_Vf< PercentVel_under_z0) P_Vf=PercentVel_under_z0*Pwc_Vf
      if(Zpar .LT. BottDepth+z0 .or. P_Wb/Pwc_Wb< PercentVel_under_z0) P_Wb=PercentVel_under_z0*Pwc_Wb
      if(Zpar .LT. BottDepth+z0 .or. P_Wc/Pwc_Wc< PercentVel_under_z0) P_Wc=PercentVel_under_z0*Pwc_Wc
      if(Zpar .LT. BottDepth+z0 .or. P_Wf/Pwc_Wf< PercentVel_under_z0) P_Wf=PercentVel_under_z0*Pwc_Wf
    
      !     *********************************************************
      !     *        Find Internal b,c,f and Advection Values       *
      !     *********************************************************
      !
      ! ii. fit polynomial to hydrodynamic model output and find internal 
      !     b,c,f values

      !a. U velocity
      ! 1. Prepare external time step values
      if (p .EQ. 1) then
        ey=0.0
        ey(1) = P_Ub
        ey(2) = P_Ub
        ey(3) = P_Uc
      else
        ey=0.0
        ey(1) = P_Ub
        ey(2) = P_Uc
        ey(3) = P_Uf
      endif

      ! 2. Get Advection value
      Uad = polintd(ex,ey,3,ix(version))

      !b. V velocity
      ! 1. Prepare external time step values
      if (p .EQ. 1) then
        ey=0.0
        ey(1) = P_Vb
        ey(2) = P_Vb
        ey(3) = P_Vc
      else
        ey=0.0
        ey(1) = P_Vb
        ey(2) = P_Vc  ! WAS CANCELLED, WHY ???
        ey(3) = P_Vf
      endif

      ! 2. Get Advection value
      Vad = polintd(ex,ey,3,ix(version))    ! WAS CANCELLED, WHY ??? 

      !c. W velocity
      ! 1. Prepare external time step values
      if (p .EQ. 1) then
        ey=0.0
        ey(1) = P_Wb
        ey(2) = P_Wb
        ey(3) = P_Wc
      else
        ey=0.0
        ey(1) = P_Wb
        ey(2) = P_Wc
        ey(3) = P_Wf
      endif

      ! 2. Get Advection value
      Wad = polintd(ex,ey,3,ix(version))

    elseif (NuminterpLvlii>0 .and. NumInterpLvliii>0) then
      if (( (Zpar .GT. Pwc_zb(us)) .or. (Zpar .GT. Pwc_zc(us)) &
          .or. (Zpar .GT. Pwc_zf(us) .and. p.ne.1) ).and. &
           (.not.VInterpUVinSurfWater) )then
 
          P_Ub = getInterp(Xpar,Ypar,VAR_ID_uvelb,us)
          P_Uc = getInterp(Xpar,Ypar,VAR_ID_uvelc,us)
          P_Uf = getInterp(Xpar,Ypar,VAR_ID_uvelf,us)
          P_Vb = getInterp(Xpar,Ypar,VAR_ID_vvelb,us)
          P_Vc = getInterp(Xpar,Ypar,VAR_ID_vvelc,us)
          P_Vf = getInterp(Xpar,Ypar,VAR_ID_vvelf,us)
                 
          !a. U velocity
          ! 1. PrgetIepare external time step values
          if (p .EQ. 1) then
            ey=0.0
            ey(1) = P_Ub
            ey(2) = P_Ub
            ey(3) = P_Uc
          else
            ey=0.0
            ey(1) = P_Ub
            ey(2) = P_Uc
            ey(3) = P_Uf
          endif

          ! 2. Get Advection value
          Uad = polintd(ex,ey,3,ix(version))

          !b. V velocity
          ! 1. Prepare external time step values
          if (p .EQ. 1) then
            ey=0.0
            ey(1) = P_Vb
            ey(2) = P_Vb
            ey(3) = P_Vc
          else
            ey=0.0
            ey(1) = P_Vb
            ey(2) = P_Vc  ! WAS CANCELLED, WHY ???
            ey(3) = P_Vf
          endif
          ! 2. Get Advection value
          Vad = polintd(ex,ey,3,ix(version))    ! WAS CANCELLED, WHY ??? 
      
      else
          Uad = WCTS_ITPI(VAR_ID_uvel,Xpar,Ypar,ii ,Pwc_zb ,Pwc_zc ,Pwc_zf ,us,P_zb,    &
                      P_zc,P_zf,ex,ix,p,version,n,NumInterpLvlii)
          Vad = WCTS_ITPI(VAR_ID_vvel,Xpar,Ypar,ii ,Pwc_zb ,Pwc_zc ,Pwc_zf ,us,P_zb,    &
                      P_zc,P_zf,ex,ix,p,version,n,NumInterpLvlii)
      endif
      Wad = WCTS_ITPI(VAR_ID_wvel,Xpar,Ypar,iii,Pwc_wzb,Pwc_wzc,Pwc_wzf,ws,P_zb,    &
                      P_zc,P_zf,ex,ix,p,version,n,NumInterpLvliii)
    else
        write(*,*)'ERROR NuminterpLvl ',NuminterpLvlii,NuminterpLvliii
        stop 'ERROR NuminterpLvl '
    endif


    RETURN
  END SUBROUTINE find_currents

SUBROUTINE find_winds(Xpar,Ypar,ex,ix,p,version,Uadw,Vadw,n)
    !This Subroutine calculates wind induced surface drift currents at the particle's
    !  location in space and time
    USE HYDRO_MOD,  ONLY: setInterp,getInterp
    USE INT_MOD,    ONLY: linint,polintd
#include "VAR_IDs.h"
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: p,version,n
    DOUBLE PRECISION, INTENT(IN) :: Xpar,Ypar,ex(3),ix(3)
    DOUBLE PRECISION, INTENT(OUT) :: Uadw,Vadw

    DOUBLE PRECISION :: P_Uwindb,P_Uwindc,P_Uwindf,P_Vwindb,P_Vwindc,P_Vwindf, &
                        ey(3),slope

    !Set Interpolation Values for the current particle
    CALL setInterp(Xpar,Ypar,n)

        P_Uwindb = getInterp(Xpar,Ypar,VAR_ID_uwindb,1)
        P_Uwindc = getInterp(Xpar,Ypar,VAR_ID_uwindc,1)
        P_Uwindf = getInterp(Xpar,Ypar,VAR_ID_uwindf,1)
        P_Vwindb = getInterp(Xpar,Ypar,VAR_ID_vwindb,1)
        P_Vwindc = getInterp(Xpar,Ypar,VAR_ID_vwindc,1)
        P_Vwindf = getInterp(Xpar,Ypar,VAR_ID_vwindf,1)

      !     *********************************************************
      !     *        Find Internal b,c,f and Advection Values       *
      !     *********************************************************
      !
      ! ii. fit polynomial to hydrodynamic model output and find internal
      !     b,c,f values

      !a. Uwind velocity
      ! 1. Prepare external time step values
      if (p .EQ. 1) then
        ey=0.0
        ey(1) = P_Uwindb
        ey(2) = P_Uwindb
        ey(3) = P_Uwindc
      else
        ey=0.0
        ey(1) = P_Uwindb
        ey(2) = P_Uwindc
        ey(3) = P_Uwindf
      endif

      ! 2. Get Advection value
      Uadw = polintd(ex,ey,3,ix(version))

      !b. Vwind velocity
      ! 1. Prepare external time step values
      if (p .EQ. 1) then
        ey=0.0
        ey(1) = P_Vwindb
        ey(2) = P_Vwindb
        ey(3) = P_Vwindc
      else
        ey=0.0
        ey(1) = P_Vwindb
        ey(2) = P_Vwindc
        ey(3) = P_Vwindf
      endif

      ! 2. Get Advection value
      Vadw = polintd(ex,ey,3,ix(version))
      !if(Uadw==0.0 .or.Vadw ==0.0) then
      !     write(*,*)'U or V current null!',p,Uadw,Vadw
      !endif

    RETURN
  END SUBROUTINE find_winds




  subroutine printOutput()
    use param_mod,   only: numpar,SaltTempOn,TrackCollisions,WriteModelTiming, &
                           NCOutFile,OutDir,Behavior,Write_coastdist, &
                          Write_Poly_Presence, read_GrainSize
    use convert_mod, only: x2lon,y2lat
    integer :: n,mainpoly,poly
    double precision :: pLon,pLat,timeinpoly

    ! increment file number
    prcount = prcount + 1
    if(Write_Poly_Presence)then
     do n=1,numpar
       mainpoly=0
       timeinpoly=0.0
       do poly=poly0,polyN
         if(Time_in_Poly(poly-poly0+1,n)>timeinpoly)then
           timeinpoly=Time_in_Poly(poly-poly0+1,n)
           mainpoly=poly
         endif
       enddo
       P_MainPoly(n)=mainpoly
     enddo  
     Time_in_Poly(:,:)=0
    endif
    ! NON-ALLOCATED arrays set optionals in writeOutput will be perceive as NOT present by writeOutput
    CALL writeOutput(par(:,pX),par(:,pY),par(:,pZ),par(:,pAge),par(:,pStatus),    &
             prcount,HITBOTTOM=hitBottom,HITLAND=hitLand,PSALT=P_Salt,PTEMP=P_Temp,  &
             PCoDi=P_coastdist,PPoly=P_MainPoly,PGrainSize=P_GrainSize,PSize=P_Size)
   !!Based on user options, write specified data to output
   !IF(SaltTempOn)THEN
   !  IF(TrackCollisions)THEN
   !    if(Write_coastdist)then
   !    CALL writeOutput(par(:,pX),par(:,pY),par(:,pZ),par(:,pAge),par(:,pStatus),    &
   !         prcount,HITBOTTOM=hitBottom,HITLAND=hitLand,P_SALT=P_Salt,P_TEMP=P_Temp,  &
   !         P_CoDi=P_coastdist )
   !    else
   !    CALL writeOutput(par(:,pX),par(:,pY),par(:,pZ),par(:,pAge),par(:,pStatus),    &
   !         prcount,HITBOTTOM=hitBottom,HITLAND=hitLand,P_SALT=P_Salt,P_TEMP=P_Temp  )
   !    endif
   !  ELSE
   !    if(Write_coastdist)then
   !    CALL writeOutput(par(:,pX),par(:,pY),par(:,pZ),par(:,pAge),par(:,pStatus),   &
   !         prcount,P_SALT=P_Salt,P_TEMP=P_Temp,P_CoDi=P_coastdist)
   !    else
   !    CALL writeOutput(par(:,pX),par(:,pY),par(:,pZ),par(:,pAge),par(:,pStatus),   &
   !         prcount,P_SALT=P_Salt,P_TEMP=P_Temp)
   !    endif
   !  ENDIF
   !ELSE
   !  IF(TrackCollisions)THEN
   !    if(Write_coastdist)then
   !    CALL writeOutput(par(:,pX),par(:,pY),par(:,pZ),par(:,pAge),par(:,pStatus),    &
   !         prcount,HITBOTTOM=hitBottom,HITLAND=hitLand,P_CoDi=P_coastdist)
   !    else
   !    CALL writeOutput(par(:,pX),par(:,pY),par(:,pZ),par(:,pAge),par(:,pStatus),    &
   !         prcount,HITBOTTOM=hitBottom,HITLAND=hitLand)
   !    endif
   !  ELSE
   !    if(Write_coastdist)then
   !    CALL writeOutput(par(:,pX),par(:,pY),par(:,pZ),par(:,pAge),par(:,pStatus), &
   !                     prcount,P_CoDi=P_coastdist)
   !    else
   !    CALL writeOutput(par(:,pX),par(:,pY),par(:,pZ),par(:,pAge),par(:,pStatus), &
   !                     prcount)
   !    endif
   !  ENDIF
   !ENDIF

    if((Behavior.ge.8.and.Behavior.le.11).and.read_GrainSize)then
        P_GrainSize = -999.0
    endif
    if (SaltTempOn) P_Temp = -999.0

    !If Tracking Collisions, write update to .csv files
    IF(TrackCollisions)then

      OPEN(100,FILE=TRIM(OutDir)//'/LandHits'//TRIM(NCOutFile)//'.csv'  ,POSITION='APPEND')
      OPEN(101,FILE=TRIM(OutDir)//'/BottomHits'//TRIM(NCOutFile)//'.csv',POSITION='APPEND')
      101 format(I7,2(',',F9.4),',',F10.3,2(',',F10.5),',',I7)

        do n=1,numpar
          pLon = x2lon(par(n,pX),par(n,pY))
          pLat = y2lat(par(n,pY))
          !write(*,*)par(n,pAge)
          !write(*,*)hitLand
          !write(*,*)hitBottom
          if(hitLand(n)   > 0)write(100,101) n,pLon,pLat,par(n,pZ),            &
             par(n,pAge)/DBLE(3600*24),(DBLE(ix(3))/DBLE(86400)),hitLand(n)
          if(hitBottom(n) > 0)write(101,101) n,pLon,pLat,par(n,pZ),            &
             par(n,pAge)/DBLE(3600*24),(DBLE(ix(3))/DBLE(86400)),hitBottom(n)
        enddo

      CLOSE(100)
      CLOSE(101)

      !Reset Collision counters
      hitBottom = 0
      hitLand = 0
    ENDIF

    !If Tracking Model Timing, write Time data to file
    IF(WriteModelTiming)then
      call CPU_TIME(times(9))

      timeCounts(8) = times(9)-times(1)

      OPEN(300,FILE=TRIM(OutDir)//'/Timing'//TRIM(NCOutFile)//'.csv',POSITION='APPEND')

        write(300,"(2(F13.2,','),7(F13.2,',',F5.1,'%'))")   &
          (DBLE(ix(3))/DBLE(86400)),timeCounts(8),          &
          timeCounts(1),(timeCounts(1)/timeCounts(8))*DBLE(100.00),            &
          timeCounts(2),(timeCounts(2)/timeCounts(8))*DBLE(100.00),            &
          timeCounts(3),(timeCounts(3)/timeCounts(8))*DBLE(100.00),            &
          timeCounts(4),(timeCounts(4)/timeCounts(8))*DBLE(100.00),            &
          timeCounts(5),(timeCounts(5)/timeCounts(8))*DBLE(100.00),            &
          timeCounts(6),(timeCounts(6)/timeCounts(8))*DBLE(100.00),            &
          timeCounts(7),(timeCounts(7)/timeCounts(8))*DBLE(100.00)

      CLOSE(300)

      timeCounts = 0
      call CPU_TIME(times(1))
    ENDIF

  
  end subroutine printOutput


  SUBROUTINE writeOutput(x,y,z,Page,PStat,prcount,hitBottom,hitLand, &
                                      PSalt,PTemp,PCoDi,PPoly,PGrainSize,PSize)
    USE PARAM_MOD, ONLY: numpar,SaltTempOn,outpathGiven,outpath,writePARA,      &
                  writeNC,TrackCollisions,Behavior,Write_coastdist,storedincolor
    USE BEHAVIOR_MOD, ONLY: getStatus
    USE CONVERT_MOD, ONLY: x2lon,y2lat
    USE HYDRO_MOD, ONLY: writeNetCDF

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: x(numpar),y(numpar),z(numpar),Page(numpar)
    DOUBLE PRECISION, INTENT(INOUT) :: PStat(numpar)
    INTEGER         , INTENT(IN) :: prcount
    INTEGER, DIMENSION(numpar), INTENT(IN), OPTIONAL :: hitBottom,hitLand
    DOUBLE PRECISION, DIMENSION(numpar), INTENT(IN), OPTIONAL ::   &
                      PSalt,PTemp,PCoDi,PGrainSize,PSize
    INTEGER, DIMENSION(numpar), INTENT(IN), OPTIONAL ::  PPoly

    INTEGER :: n
    DOUBLE PRECISION :: statuses(numpar)
    double precision, dimension(numpar) :: pLon,pLat
    double precision :: default_stat

    !INPUT/OUTPUT FILE NAME CONSTRUCTION VARIABLES
    CHARACTER(LEN=100) :: filenm2
    CHARACTER(LEN=4  ) :: prefix2,suffix2
    INTEGER :: counter2

    !Convert particle position (in meters) to latitude and longitude &
    !Find identification number that describes a particle's behavior 
    !  type or status for use in visualization routines
    do n=1,numpar
      pLon(n) = x2lon(x(n),y(n))
      pLat(n) = y2lat(y(n))
      default_stat= par(n,pStatus)
      statuses(n) = getStatus(n,default_stat)
      !if(statuses(n)>=0) statuses(n)=PStat(n) 
    enddo

    !if writing to .csv files:
    if(writePARA)then
      !Create a filename and unit number for each iteration    
      counter2=prcount+10000000
      prefix2='para'
      suffix2='.csv'
      if(outpathGiven)then
        write(filenm2,'(A,A,I8,A)') TRIM(outpath),prefix2,counter2,suffix2
      else
        write(filenm2,'(A,I8,A)') prefix2,counter2,suffix2
      endif
      open(2,FILE=TRIM(filenm2),STATUS='REPLACE')
      5 format(F10.3,',',I7,2(',',F9.4),2(',',F8.4)) !SaltTemp

      !Based on user options, Write specified data to file
      do n = 1,numpar

        if (SaltTempOn) then
          write(2,5) z(n),int(statuses(n)),pLon(n),pLat(n),PSalt(n),PTemp(n)
        else
          write(2,5) z(n),int(statuses(n)),pLon(n),pLat(n)
        endif 

      enddo

      close(2)

    endif


    !if writing to NetCDF file(s):
    if(writeNC)then
          ! Non-allocated arrays will be seen as NOT present by writeNetCDF()
          CALL writeNetCDF(int(ix(3)),Page(:),pLon,pLat,z,statuses,           &
               SALT=PSalt,TEMP=PTemp,GrSize=PGrainSize,SizeP=PSize,           &
               HITB=hitBottom,HITL=hitLand,CoDi=PCoDi,POLY=PPoly)

     !!Based on user options, Write specified data to file
     !if (SaltTempOn) then
     !  if(TrackCollisions)then
     !   if(Behavior.ge.8.and.Behavior.le.10)then
     !    if(Write_coastdist)then
     !    CALL writeNetCDF(int(ix(3)),P_age(:),pLon,pLat,z,statuses,           &
     !         SALT=P_Salt,TEMP=P_Temp,GrSize=P_GrainSize,PSize=P_Size,           &
     !         HITB=hitBottom,HITL=hitLand,CoDi=P_CoDi)
     !    else
     !    CALL writeNetCDF(int(ix(3)),P_age(:),pLon,pLat,z,statuses,           &
     !         SALT=P_Salt,TEMP=P_Temp,GrSize=P_GrainSize,PSize=P_Size,           &
     !         HITB=hitBottom,HITL=hitLand)
     !    endif
     !   else
     !    if(Write_coastdist)then
     !    CALL writeNetCDF(int(ix(3)),P_age(:),pLon,pLat,z,statuses,           &
     !         SALT=P_Salt,TEMP=P_Temp,HITB=hitBottom,HITL=hitLand,CoDi=P_CoDi)
     !    else
     !    CALL writeNetCDF(int(ix(3)),P_age(:),pLon,pLat,z,statuses,           &
     !         SALT=P_Salt,TEMP=P_Temp,HITB=hitBottom,HITL=hitLand)
     !    endif
     !   endif
     !  else
     !   if(Behavior.ge.8.and.Behavior.le.10)then
     !    if(Write_coastdist)then
     !    CALL writeNetCDF(int(ix(3)),P_age(:),pLon,pLat,z,statuses,           &
     !         SALT=P_Salt,TEMP=P_Temp,GrSize=P_GrainSize,PSize=P_Size,CoDi=P_CoDi)
     !    else
     !    CALL writeNetCDF(int(ix(3)),P_age(:),pLon,pLat,z,statuses,           &
     !         SALT=P_Salt,TEMP=P_Temp,GrSize=P_GrainSize,PSize=P_Size)
     !    endif
     !   else
     !    if(Write_coastdist)then
     !    CALL writeNetCDF(int(ix(3)),P_age(:),pLon,pLat,z,statuses,           &
     !         SALT=P_Salt,TEMP=P_Temp,CoDi=P_CoDi)
     !    else
     !    CALL writeNetCDF(int(ix(3)),P_age(:),pLon,pLat,z,statuses,           &
     !         SALT=P_Salt,TEMP=P_Temp)
     !    endif
     !   endif
     !  endif
     !else
     !  if(TrackCollisions)then
     !   if(Behavior.ge.8.and.Behavior.le.10)then
     !    if(Write_coastdist)then
     !    CALL writeNetCDF(int(ix(3)),P_age(:),pLon,pLat,z,statuses,           &
     !         GrSize=P_GrainSize,PSize=P_Size,HITB=hitBottom,HITL=hitLand,       &
     !         CoDi=P_CoDi)
     !    else
     !    CALL writeNetCDF(int(ix(3)),P_age(:),pLon,pLat,z,statuses,           &
     !         GrSize=P_GrainSize,PSize=P_Size,HITB=hitBottom,HITL=hitLand)
     !    endif
     !   else
     !    if(Write_coastdist)then
     !    CALL writeNetCDF(int(ix(3)),P_age(:),pLon,pLat,z,statuses,           &
     !         HITB=hitBottom,HITL=hitLand,CoDi=P_CoDi)
     !    else
     !    CALL writeNetCDF(int(ix(3)),P_age(:),pLon,pLat,z,statuses,           &
     !         HITB=hitBottom,HITL=hitLand)
     !    endif
     !   endif
     !  else
     !   if(Behavior.ge.8.and.Behavior.le.10)then
     !    if(Write_coastdist)then
     !    CALL writeNetCDF(int(ix(3)),P_age(:),pLon,pLat,z,statuses,           &
     !         GrSize=P_GrainSize,PSize=P_Size,CoDi=P_CoDi)
     !    else
     !    CALL writeNetCDF(int(ix(3)),P_age(:),pLon,pLat,z,statuses,           &
     !         GrSize=P_GrainSize,PSize=P_Size)
     !    endif
     !   else
     !    if(Write_coastdist)then
     !    CALL writeNetCDF(int(ix(3)),P_age(:),pLon,pLat,z,statuses,CoDi=P_CoDi)
     !    else
     !    CALL writeNetCDF(int(ix(3)),P_age(:),pLon,pLat,z,statuses)
     !    endif
     !   endif
     !  endif
     !endif
    endif

  END SUBROUTINE writeOutput

  SUBROUTINE writeModelInfo()
    !This subroutine simply writes model information to standard output
    USE PARAM_MOD
    IMPLICIT NONE

    CHARACTER(len=10) :: tmp !For Converting Integers to Characters
    CHARACTER(len=200) :: filenm
    INTEGER :: nf,nfilesin

    write(*,*) ' ******************** Model Info ******************** '
    write(*,*) ' '
    write(*,*) ' LTRANS-Zlev version 0(beta)'
    write(*,*) ' '

    write(*,*) ' Run Name:              = ',TRIM(RunName)
    write(*,*) ' Executable Directory:  = ',TRIM(ExeDir)
    write(*,*) ' Output Directory:      = ',TRIM(OutDir)
    write(*,*) ' Run By:                = ',TRIM(RunBy)
    write(*,*) ' Institution:           = ',TRIM(Institution)
    write(*,*) ' Started On:            = ',TRIM(StartedOn)
    write(*,*) ' '

    write(tmp,'(F10.3)') days
    tmp = ADJUSTL(tmp)
    write(*,*) ' Days:                  = ',TRIM(tmp)
    write(tmp,'(I10)') numpar
    tmp = ADJUSTL(tmp)
    write(*,*) ' Particles:             = ',TRIM(tmp)
    write(*,*) ' Particle File:         = ',TRIM(parfile)
    write(*,*) ' '

    SELECT CASE(Behavior)
      CASE(0)
      write(*,*) ' Behavior:        = Passive'
      CASE(1)
      write(*,*) ' Behavior:        = Near-Surface'
      CASE(2)
      write(*,*) ' Behavior:        = Near-Bottom'
      CASE(3)
      write(*,*) ' Behavior:        = Diurnal Vertical Migration'
      CASE(4)
      write(*,*) ' Behavior:        = C.virginica oyster larvae'
      CASE(5)
      write(*,*) ' Behavior:        = C.ariakensis oyster larvae'
      CASE(6)
      write(*,*) ' Behavior:        = Constant sink/float'
      CASE(7)
      write(*,*) ' Behavior:        = Tidal Stream Transport'
      CASE(8)
      write(*,*) ' Behavior:        = Nephrops larvae'
      CASE(9)
      write(*,*) ' Behavior:        = Solea Solea larvae'
      CASE(10)
      write(*,*) ' Behavior:        = Mullus Barbatus larvae'
      CASE(11)
      write(*,*) ' Behavior:        = Parameterizable larvae behavior'
      CASE(997)
      write(*,*) ' Behavior:        = Die when reaching surface'
      CASE(998)
      write(*,*) ' Behavior:        = Keep constant absolute depth under zero'
      CASE(999)
      write(*,*) ' Behavior:        = Keep constant depth under zeta'
      CASE(1000)
      write(*,*) ' Behavior:        = Oil spill'
      CASE DEFAULT
      write(*,*) ' Behavior:        = ',Behavior,' UNKNOWN !'
      write(*,*) 'The Program Cannot Continue and is Terminating'
      stop
    END SELECT

    if(mortality)then
      write(*,*) ' Particle Mortality:    = On'
    else
      write(*,*) ' Particle Mortality:    = Off'
    endif

    if(settlementon)then
      write(*,*) ' Settlement:            = On'
      write(*,*) ' Habitat File:          = ',TRIM(habitatfile)
      if(holesExist)write(*,*) ' Hole File:             = ',TRIM(holefile)
    else
      write(*,*) ' Settlement:            = Off'
    endif
    write(*,*) ' '

    if(HTurbOn)then
      write(*,*) ' Horizontal Turbulence: = On'
    else
      write(*,*) ' Horizontal Turbulence: = Off'
    endif
    if(VTurbOn)then
      write(*,*) ' Vertical Turbulence:   = On'
    else
      write(*,*) ' Vertical Turbulence:   = Off'
    endif
    if(SphericalProjection)then
      write(*,*) ' Projection:            = Spherical'
    else
      write(*,*) ' Projection:            = Mercator'
    endif
    if(OpenOceanBoundary)then
      write(*,*) ' Ocean Boundary:        = Open'
    else
      write(*,*) ' Ocean Boundary:        = Closed'
    endif
    if(SaltTempOn)then
      write(*,*) ' Salt & Temp Output:    = On'
    else
      write(*,*) ' Salt & Temp Output:    = Off'
    endif
    if(TrackCollisions)then
      write(*,*) ' Track Collisions:      = Yes'
    else
      write(*,*) ' Track Collisions:      = No'
    endif
    if(WriteModelTiming)then
      write(*,*) ' Track Model Timing:    = Yes'
    else
      write(*,*) ' Track Model Timing:    = No'
    endif

    if(Zgrid)then
        nfilesin=8
       if(Wind) nfilesin=nfilesin+2 
       if(WindIntensity) nfilesin=nfilesin+1
       !if(readNetcdfSwdown) nfilesin=nfilesin+1
    else 
      nfilesin=1
    endif
    DO nf=1,nfilesin
    SELECT CASE(numdigits)
       CASE(0)
         WRITE(filenm,'(A,A,A)')     TRIM(dirin),TRIM(prefix(nf)),TRIM(suffix)
      CASE(1)
        WRITE(filenm,'(A,A,I1.1,A)') TRIM(dirin),TRIM(prefix(nf)),filenum,TRIM(suffix)
      CASE(2)
        WRITE(filenm,'(A,A,I2.2,A)') TRIM(dirin),TRIM(prefix(nf)),filenum,TRIM(suffix)
      CASE(3)
        WRITE(filenm,'(A,A,I3.3,A)') TRIM(dirin),TRIM(prefix(nf)),filenum,TRIM(suffix)
      CASE(4)
        WRITE(filenm,'(A,A,I4.4,A)') TRIM(dirin),TRIM(prefix(nf)),filenum,TRIM(suffix)
      CASE(5)
        WRITE(filenm,'(A,A,I5.5,A)') TRIM(dirin),TRIM(prefix(nf)),filenum,TRIM(suffix)
      CASE(6)
        WRITE(filenm,'(A,A,I6.6,A)') TRIM(dirin),TRIM(prefix(nf)),filenum,TRIM(suffix)
      CASE(7)
        WRITE(filenm,'(A,A,I7.7,A)') TRIM(dirin),TRIM(prefix(nf)),filenum,TRIM(suffix)
      CASE(8)
        WRITE(filenm,'(A,A,I8.8,A)') TRIM(dirin),TRIM(prefix(nf)),filenum,TRIM(suffix)
      CASE(9)
        WRITE(filenm,'(A,A,I9.9,A)') TRIM(dirin),TRIM(prefix(nf)),filenum,TRIM(suffix)
      CASE(10)
        WRITE(filenm,'(A,A,I10.10,A)') TRIM(dirin),TRIM(prefix(nf)),filenum,TRIM(suffix)
      CASE DEFAULT
        WRITE(*,*) 'Model presently does not support numdigits of ',numdigits
        WRITE(*,*) 'Please use numdigit value from 1 to 10'
        WRITE(*,*) '  OR modify code in Hydrodynamic module'
        STOP
    END SELECT
    !write(*,*)'file to open is ',TRIM(filenm)
    ENDDO

    write(*,*) ' '
    write(*,*) ' Grid File:             = ',NCgridfile
    write(*,*) ' First Hydro File:      = ',TRIM(filenm)

    write(*,*) ' '
    write(tmp,'(I10)') seed
    tmp = ADJUSTL(tmp)
    write(*,*) ' Seed:                  = ',TRIM(tmp)
    write(*,*) ' '

  END SUBROUTINE writeModelInfo

  SUBROUTINE handleERROR(ErrorName,callnum,errornum,n, & 
             Xpos,Ypos,Zpos,nXpos,nYpos,nZpos,docycle)
    USE PARAM_MOD,      ONLY: ErrorFlag,dt,idt,OutDir,NCOutFile
    use convert_mod,  only: x2lon,y2lat
    USE BEHAVIOR_MOD,   ONLY: setOut,die
    USE HYDRO_MOD,      ONLY: getKRlevel
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: callnum,errornum,n
    DOUBLE PRECISION, INTENT(IN) ::Xpos,Ypos,Zpos,nXpos,nYpos,nZpos 
    CHARACTER(LEN=10), INTENT(IN):: ErrorName
    CHARACTER(LEN=10) :: annotate
    LOGICAL, INTENT(OUT) :: docycle
    INTEGER:: klev,nklev
    !Error Handling Formats
    21 FORMAT ('Particle ',I10,' not in rho element after ',F15.2,' seconds')
    22 FORMAT ('Particle ',I10,' not in u element after '  ,F15.2,' seconds')
    23 FORMAT ('Particle ',I10,' not in v element after '  ,F15.2,' seconds')
    24 FORMAT ('Particle ',I10,' out after 3rd reflection after ',F15.2,         &
               ' seconds')
    25 FORMAT ('Particle ',I10,                                                &
               ' outside main bounds after intersect_reflect after ',F15.2,      &
               ' seconds')
    26 FORMAT ('Particle ',I10,                                                &
               ' inside island bounds after intersect_reflect after ',F15.2,     &
               ' seconds')
    27 FORMAT ('Particle ',I10,' jumped over rho element after ',F15.2,          &
               ' seconds')
    28 FORMAT ('Particle ',I10,' jumped over u element after ',F15.2,' seconds')
    29 FORMAT ('Particle ',I10,' jumped over v element after ',F15.2,' seconds')
    30 FORMAT ('Prev Loc: x=',F11.3,' y=',F11.3,' lon=',F8.3,' lat=',F8.3, &
               ' Z=',F8.3,' k=',i3,/, &
               'New  Loc: x=',F11.3,' y=',F11.3,' lon=',F8.3,' lat=',F8.3, & 
               ' Z=',F8.3,' k=',i3)
    31 FORMAT ('AdvX=',F8.3,' TurbHx=',F8.3,' XBehav=',F8.3,/, & 
               'AdvY=',F8.3,' TurbHy=',F8.3,' YBehav=',F8.3,/, & 
               'AdvZ=',F8.3,' TurbV =',F8.3,' ZBehav=',F8.3, &
               ' reflect=',F8.3,' zeta=',F8.3,' Depth=',F8.3,/ )
    32 FORMAT ('Particle ',I10,' got out on land  after ',F15.2,' seconds')
    33 FORMAT ('Particle ',I10,' not found in a rho element after ',F15.2,          &
               ' seconds')
    34 FORMAT ('Particle ',I10,' not found in a u element after ',F15.2,' seconds')
    35 FORMAT ('Particle ',I10,' not found in a v element after ',F15.2,' seconds')




    klev=getKRlevel(Zpos)
    nklev=getKRlevel(nZpos)
    docycle=.False.
    IF(Trim(ErrorName).eq.'setEle')then
        write(annotate,'(i1,a,i6)')callnum,'sEl',it+p*int(dt/idt)
        call writeErrortopython(annotate,'r',Xpos,Ypos)   
        if(ErrorFlag < 1 .OR. ErrorFlag > 3)then
          write(*,*) " "
          write(*,*) "part ",n," not in element"
          SELECT CASE (errornum)
            CASE(1)
              write(*,*) 'part ',n,"not found in a rho element"
            CASE(2)
              write(*,*) 'part ',n,"not found in a u element"
            CASE(3)
              write(*,*) 'part ',n,"not found in a v element"
            CASE(4)
              write(*,*) 'part ',n,' Jumped over a rho element'
            CASE(5)
              write(*,*) 'part ',n,' Jumped over a u element'
            CASE(6)
              write(*,*) 'part ',n,' Jumped over a v element'
            CASE(7)
              write(*,*) 'part ',n," setEle error 7 ."
          END SELECT
          write(*,*) ' '
          write(*,30) Xpos,Ypos,x2lon(Xpos,Ypos),     &
                     y2lat(Ypos),Zpos,klev,     &
                      nXpos,nYpos,x2lon(nXpos,nYpos),     &
                     y2lat(nYpos),nZpos,nklev    
          write(*,*) 'The Program Cannot Continue and Will Terminate'
          stop "error setEle"
        else
          OPEN(210,FILE=TRIM(OutDir)//'/ErrorLog'//TRIM(NCOutFile)//'.txt',POSITION='APPEND')
            SELECT CASE (errornum)
              CASE(1)
                write(210,33) n,ix(3)
              CASE(2)
                write(210,34) n,ix(3)
              CASE(3)
                write(210,35) n,ix(3)
              CASE(4)
                write(210,21) n,ix(3)
              CASE(5)
                write(210,22) n,ix(3)
              CASE(6)
                write(210,23) n,ix(3)
              CASE(7)
              write(210,*) 'part ',n," setEle error 7 ."
            END SELECT
            write(210,30) Xpos,Ypos,x2lon(Xpos,Ypos),     &
                     y2lat(Ypos),Zpos,klev,     &
                      nXpos,nYpos,x2lon(nXpos,nYpos),     &
                     y2lat(nYpos),nZpos,nklev    
          write(210,*)' after set_Ele 1 '
          CLOSE(210)
          docycle=.True.
        endif
    ELSEIF(trim(ErrorName).eq.'Fstlev' .or. trim(ErrorName).eq.'gDepth')then
           write(*,*)'getDepth conflict1 part ',n,it
           write(*,'(3(a,f14.8),3a)') &
             'mlab.points3d([',x2lon(Xpos,Ypos),'],[',y2lat(Ypos), &
             '],[(np.max(depth)-depthfactor*',abs(Zpos),")],",&
              "scale_mode='none',scale_factor=scalefac*3,",        &
              "color=(0.0,1.0,1.0))# old pos CONFLICTplot"
             write(*,'(3(a,f14.8),3a)') &
             'mlab.points3d([',x2lon(nXpos,nYpos),'],[',y2lat(nYpos),&
             '],[(np.max(depth)-depthfactor*',abs(nZpos),")],",&
              "scale_mode='none',scale_factor=scalefac*3,",        &
              "color=(1.0,0.0,0.0))# new pos CONFLICTplot"
             write(*,'(6(a,f14.8),3a)') &
             'mlab.plot3d([',x2lon(Xpos,Ypos),',',x2lon(nXpos,nYpos),&
              '],[',y2lat(Ypos),',',y2lat(nYpos),                    &
             '],[np.max(depth)-depthfactor*',                       &
              abs(Zpos),',np.max(depth)-depthfactor*',abs(nZpos), &
               "],","representation='wireframe',tube_radius=None,",     &
              "color=(0.0,0.0,0.0),line_width=4.0) # CONFLICTplot"
             write(*,'(3(a,f14.8),a,i4,a,i8,a)') &
              'mlab.text3d(',0.5*(x2lon(Xpos,Ypos)+x2lon(nXpos,nYpos)),&
              ',',0.5*(y2lat(Ypos)+y2lat(nYpos)),                      &
             ',(np.max(depth)-depthfactor*',                           &
              0.5*(abs(Zpos)+abs(nZpos)),"),'GD1-",n,'-',it,     &
              "') # CONFLICTplot"
        write(annotate,'(i1,a,i1,i6)')callnum,'gDp',it+p*int(dt/idt)
        call writeErrortopython(annotate,'r',Xpos,Ypos)   
           if(ErrorFlag < 1 .OR. ErrorFlag > 3)then
            write(*,*)'Particle ',n,' error on Pdepth computation',errornum
            stop 'Part out before new pos computation'
           endif
           OPEN(210,FILE=TRIM(OutDir)//'/ErrorLog'//TRIM(NCOutFile)//'.txt',POSITION='APPEND')
                 write(210,32) n,ix(3)
                 write(210,'(2(a,f35.30))')'at  LONpos=',x2lon(nXpos,nYpos), &
                                      '    LATpos=',y2lat(nYpos)
                 !write(210,*)'At this last pos Fstlev=',Fstlev,'>us=',us,' at the beginning of the loop'
                 write(210,*)' '
           CLOSE(210)
          docycle=.True.

    ELSEIF(trim(ErrorName).eq.'intersect')then
      write(annotate,'(a,i6)')'>3Rx',it+p*int(dt/idt)
      !call writeErrortopython(annotate,'r',fintersectX,fintersectY,        &
       !                                 freflectX,freflectY)
        !                                       freflectX,freflectY)
          if(ErrorFlag < 1 .OR. ErrorFlag > 3)then
            write(*,*) n,'still out after 3rd reflection.'
            write(*,*) 'Particle set Out of Bounds and no longer tracked'
            write(*,*) ' '
            write(*,*) 'The Program Cannot Continue and Will Terminate'
            stop 'Particle set Out of Bounds and no longer tracked'
          else
            OPEN(210,FILE=TRIM(OutDir)//'/ErrorLog'//TRIM(NCOutFile)//'.txt',POSITION='APPEND')
              write(210,24) n,ix(3)
              write(210,30) Xpos,Ypos,x2lon(Xpos,Ypos),     &
                     y2lat(Ypos),Zpos,klev,     &
                     nXpos,nYpos,x2lon(nXpos,nYpos),              &
                     y2lat(nYpos),nZpos,nklev
              !write(210,31)   AdvectX,TurbHx,XBehav, &
              !                AdvectY,TurbHy,YBehav, &
              !                AdvectZ,TurbV,ZBehav, reflect,P_zetac,P_depth
            write(210,*)' after intersect_reflect'
            CLOSE(210)
          endif
         docycle=.True.
    ELSEIF(ErrorName.eq.'mbounds')THEN
        write(annotate,'(a,i6)')'mbnd',it+p*int(dt/idt)
        call writeErrortopython(annotate,'r',Xpos,Ypos,Xpos,Ypos,        &
                                                           nXpos,nYpos)       

        if(ErrorFlag < 1 .OR. ErrorFlag > 3)then
          write(*,*) 'Model Run Cannot Continue'
          write(*,*) ' '
          write(*,30) Xpos,Ypos,x2lon(Xpos,Ypos),     &
                     y2lat(Ypos),Zpos,klev,     &
                     nXpos,nYpos,x2lon(nXpos,nYpos),              &
                     y2lat(nYpos),nZpos,nklev
          stop 'Model Run Cannot Continue after mbounds'
        else
          write(*,*) 'ERROR: Particle Outside Main Boundaries ',ErrorFlag
          OPEN(210,FILE=TRIM(OutDir)//'/ErrorLog'//TRIM(NCOutFile)//'.txt',POSITION='APPEND')
            write(210,25) n,ix(3)
            write(210,30) Xpos,Ypos,x2lon(Xpos,Ypos),     &
                     y2lat(Ypos),Zpos,klev,     &
                     nXpos,nYpos,x2lon(nXpos,nYpos),              &
                     y2lat(nYpos),nZpos,nklev
            write(210,*)' after mbounds'
          CLOSE(210)
          docycle=.True.
        endif
    ELSEIF(ErrorName.eq.'ibounds')THEN
        write(annotate,'(a,i6)')'ibnd',it+p*int(dt/idt)
        call writeErrortopython(annotate,'r',Xpos,Ypos,Xpos,Ypos,    &
                                                           nXpos,nYpos)       

        if(ErrorFlag < 1 .OR. ErrorFlag > 3)then
          write(*,*) 'Model Run Cannot Continue'
          write(*,*) ' '
          write(*,30) Xpos,Ypos,x2lon(Xpos,Ypos),     &
                     y2lat(Ypos),Zpos,klev,     &
                     nXpos,nYpos,x2lon(nXpos,nYpos),              &
                     y2lat(nYpos),nZpos,nklev
          stop 'Model Run Cannot Continue after ibounds'
        else
          write(*,*) 'ERROR: Particle Inside Island Boundaries ',ErrorFlag
          OPEN(210,FILE=TRIM(OutDir)//'/ErrorLog'//TRIM(NCOutFile)//'.txt',POSITION='APPEND')
            write(210,26) n,ix(3)
            write(210,30) par(n,pX),Ypos,x2lon(par(n,pX),Ypos),     &
                     y2lat(Ypos),Zpos,klev,     &
                     nXpos,nYpos,x2lon(nXpos,nYpos),              &
                     y2lat(nYpos),nZpos,nklev
            write(210,*)' after ibounds'
          CLOSE(210)
          docycle=.True.
        endif

    ELSEIF(ErrorName.eq.'KRlevjump')THEN
         !f(P_oldLev(n)>=0 .and. nklev>P_oldLev(n)+1 &
         !  ! .and. (.not.Zgrid_depthinterp)    &
         !    )then
         !  write(*,*)'#################################' 
         !  write(*,'(3(a,i5),7(a,f7.2))')                   &
         !  'ERROR particle ',n,' Jumped of vlev from',    &
         !  P_oldLev(n),' to ',nklev,' nZ[',nZpos,         &
         !  ']=min(depth[',P_depth, '], Z[',Zpos,       &
         !   ']+Adv[', AdvectZ,']+Tu[',TurbV,']+Bhv[',ZBehav,&
         !   ']) +reflect[',reflect,']'
         !   write(*,'(3(a,f14.8),3a)') &
         !   'mlab.points3d([',x2lon(Xpos,Ypos),'],[',y2lat(Ypos), &
         !   '],[(np.max(depth)-depthfactor*',abs(Zpos),")],",&
         !    "scale_mode='none',scale_factor=scalefac*3,",        &
         !    "color=(0.0,1.0,1.0),mode='cube')# old pos JUMPLEVplot"
         !   write(*,'(3(a,f14.8),3a)') &
         !   'mlab.points3d([',x2lon(nXpos,nYpos),'],[',y2lat(nYpos),&
         !   '],[(np.max(depth)-depthfactor*',abs(nZpos),")],",&
         !    "scale_mode='none',scale_factor=scalefac*3,",        &
         !    "color=(1.0,0.0,0.0),mode='cube')# new pos JUMPLEVplot"
         !   write(*,'(6(a,f14.8),3a)') &
         !   'mlab.plot3d([',x2lon(Xpos,Ypos),',',x2lon(nXpos,nYpos),&
         !    '],[',y2lat(Ypos),',',y2lat(nYpos),                    &
         !   '],[np.max(depth)-depthfactor*',                       &
         !    abs(Zpos),',np.max(depth)-depthfactor*',abs(nZpos),&
         !    "],","representation='wireframe',tube_radius=None,",     &
         !    "color=(0.0,0.0,0.0),line_width=4.0) # JUMPLEVplot"
         !   write(*,'(3(a,f14.8),a,i4,a,i8,a)') &
         !    'mlab.text3d(',0.5*(x2lon(Xpos,Ypos)+x2lon(nXpos,nYpos)),&
         !    ',',0.5*(y2lat(Ypos)+y2lat(nYpos)),                      &
         !   ',(np.max(depth)-depthfactor*',                           &
         !    0.5*(abs(Zpos)+abs(nZpos)),"),'",n,'-',it,        &
         !    "') # JUMPLEVplot"

         !   write(*,'(6(a,f12.7))') &
         !     ' while going from pos ',x2lon(Xpos,Ypos),',',y2lat(Ypos),&
         !     ',',Zpos,   ' to pos ',x2lon(nXpos,nYpos),',',y2lat(nYpos),&
         !     ',',nZpos
         !  If(  Zpos<P_depth)then
         !  conflict_tmp=-1 
         !  write(*,*)'---- call getDepth on old pos ----' 
         !  call getDepth(Xpos,Ypos,n,it,P_depth,Fstlev,conflict_tmp)
         !  if(conflict_tmp<0) write(*,*)'getDepth conflict2 part ',n,it
         !  write(*,*)' old depth was ',P_depth,Fstlev,conflict_tmp
         !  conflict_tmp=-1  
         !  write(*,*)'---- call getDepth on new pos ----' 
         !  call getDepth(nXpos,nYpos,n,it,P_depth,Fstlev,conflict_tmp)
         !  if(conflict_tmp<0) write(*,*)'getDepth conflict3 part ',n,it
         !  write(*,*)' new depth is  ',P_depth,Fstlev,conflict_tmp
         !  write(*,*)'---- call getDepth DONE       ----' 
         !  if(ErrorFlag < 1 .OR. ErrorFlag > 3)then
         !    write(*,*) " "
         !    write(*,*) "part ",n," Jumped of vertical level"
         !    write(*,*) 'The Program Cannot Continue and Will Terminate'
         !    stop "error klev  "
         !  else
         !    write(*,*) "part ",n," Jumped of vertical level. ",&
         !     " End of tracking for this particle."
         !    write(*,*)'#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-'
         !    OPEN(210,FILE=TRIM(OutDir)//'/ErrorLog'//TRIM(NCOutFile)//'.txt',POSITION='APPEND')
         !      write(210,* ) 'particle ',n,' jumped vertlev at time ',ix(3)
         !      write(210,30) par(n,pX),Ypos,x2lon(par(n,pX),Ypos),     &
         !               y2lat(Ypos),Zpos,P_oldLev(n),     &
         !               nXpos,nYpos,x2lon(nXpos,nYpos),     &
         !               y2lat(nYpos),nZpos,nklev
         !    CLOSE(210)
         !    cycle
         !  endif
         !  Endif 
         !ndif

    ENDIF


    IF(ErrorFlag == 1)then
       par(n,pnX) = par(n,pX)
       par(n,pnY) = par(n,pY)
       par(n,pnZ) = par(n,pZ)
    ELSEIF(ErrorFlag == 2)then
       call die(n)
    ELSEIF(ErrorFlag ==3) then
       call setOut(n)
    ENDIF

 
  END SUBROUTINE  handleERROR

  !-----------------------------------------------------------------------------

  SUBROUTINE writeErrortopython(annotate,col,X1,Y1,X2,Y2,X3,Y3,X4,Y4,X5,Y5)
    use param_mod, only :    NCOutFile,OutDir
    use convert_mod,  only: x2lon,y2lat  
    IMPLICIT NONE
    CHARACTER(LEN=10), INTENT(IN):: annotate
    CHARACTER(LEN=1), INTENT(IN):: col
    DOUBLE PRECISION, INTENT(IN) :: X1,Y1
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: X2,Y2,X3,Y3,X4,Y4,X5,Y5
    DOUBLE PRECISION, DIMENSION(5) :: X,Y
    INTEGER :: i,num
    DOUBLE PRECISION :: length
    X(1)=x2lon(X1,Y1)
    Y(1)=y2lat(   Y1)
    num=1
    if(present(X2).and.present(Y2))then
     X(2)=x2lon(X2,Y2)
     Y(2)=y2lat(   Y2)
     num=2
    endif
    if(present(X3).and.present(Y3))then
     num=3
     X(3)=x2lon(X3,Y3)
     Y(3)=y2lat(   Y3)
    endif
    if(present(X4).and.present(Y4))then
     num=4
     X(4)=x2lon(X4,Y4)
     Y(4)=y2lat(   Y4)
    endif
    if(present(X5).and.present(Y5))then
     num=5
     X(5)=x2lon(X5,Y5)
     Y(5)=y2lat(   Y5)
    endif
    open (unit = fpy,file = TRIM(OutDir)//'/'//TRIM(NCOutFile)//               &
                         'PartErrors'//'.py',POSITION='APPEND')
    IF(num>1)then
      DO i=1,num-1
         length=sqrt((X(i+1)-X(i))**2+(Y(i+1)-Y(i))**2)
         write(fpy,'(2a,4(a,f12.9),3a)')                                       &
            "plt.annotate('",TRIM(annotate),                                   &
            "', xy = (",   X(i)+0.5*(X(i+1)-X(i)),",",Y(i)+0.5*(Y(i+1)-Y(i)),  &
            "), xytext= (",X(1)+0.6*(X(num)-X(1)),",",Y(1)+0.6*(Y(num)-Y(1)),  &
        "), bbox = dict(boxstyle = 'round,pad=0.05', fc ='",col,               &
        "', alpha = 0.3),arrowprops=dict(facecolor='black', shrink=0.01))"
          if(length>1e-9)then
          write(fpy,'(7(a,f12.9),3a)')                                         &
                   'ax.arrow(',X(i),',',Y(i),',',                              &
              0.9*(X(i+1)-X(i)),',',0.9*(Y(i+1)-Y(i)),                         &
              ", head_width=",0.025*length,                                    &
              ", head_length=",0.1*length,                                     &
              ", width=",0.005*length,                                         &
              ", fc='none',ec='",col,"')"
          endif
          write(fpy,'(2(a,f12.9),3a)')                                         &
              "ax.scatter(",X(i),",",Y(i),",c='none',edgecolor='",             &
               col,"',marker='o',s=60)"
          write(fpy,'(2(a,f12.9),3a)')                                         &
               "ax.scatter(",X(i+1),",",Y(i+1),",c='none',edgecolor='",        &
               col,"',marker='x',s=60)"
    ENDDO
  ELSE
         i=1
         write(fpy,'(2(a,f12.9),3a)')                                          &
              "ax.scatter(",X(i),",",Y(i),",c='none',edgecolor='",             &
               col,"',marker='*',s=60)"
         write(fpy,'(2a,4(a,f12.9),3a)')                                       &
            "plt.annotate('",TRIM(annotate),                                   &
            "', xy = (",   X(i),",",Y(i),                                      &
            "), xytext= (",X(i)+0.0001           ,",",Y(i)+0.0001           ,  &
        "), bbox = dict(boxstyle = 'round,pad=0.05', fc ='",col,               &
        "', alpha = 0.3),arrowprops=dict(facecolor='black', shrink=0.01))"
  ENDIF
  close(fpy)
  END  SUBROUTINE writeErrortopython

  !-----------------------------------------------------------------------------


  SUBROUTINE   set_Average_Water_Depth()
    USE PARAM_MOD, ONLY: ui,vi,uj,vj,us,ws,constTemp,constUwind,constVwind, &
                         numpar,idt,Zgrid,Wind,SaltTempOn,WindIntensity,pi, &
                         OilOn,settlementon,mortality,OpenOceanBoundary
    USE HYDRO_MOD, ONLY: WCTS_ITPI,getKRlevel,getDepth, &
                         getSlevel,getWlevel,setInterp,getInterp
    USE INT_MOD,    ONLY: polintd
    use behavior_mod, only: die,isOut,isDead
    USE SETTLEMENT_MOD, ONLY: isSettled,testSettlement,isStranded

#include "VAR_IDs.h"
    IMPLICIT NONE

    INTEGER :: i,deplvl,n
    INTEGER :: klev,Fstlev,NumInterpLvl,ixnum,nklev,k,conflict                     !--- CL:OGS

    ! Particle tracking
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: Pwc_zb,Pwc_zc,Pwc_zf
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: Pwc_wzb,Pwc_wzc,Pwc_wzf
    DOUBLE PRECISION :: Xpar,Ypar,Zpar,newXpos,newYpos,newZpos,P_zb,P_zc,P_zf, &
      P_depth,nP_depth,P_angle,P_zeta,P_zetab,P_zetac,P_zetaf,ey(3),  &
      AdvecUwind,AdvecVwind

    DOUBLE PRECISION :: x1,x2,x3,y1,y2,y3,z1,z2,z3,slope,length,Ttemp

        DOUBLE PRECISION::      P_hsig,P_tm01,P_Uwind,P_Vwind,P_pdir,P_wlen,UWindDrift,&
                                VWindDrift,alpha,PWind,P_Uw,P_Vw,Uadw,Vadw,            &
                                UStokesDrift,VStokesDrift
        DOUBLE PRECISION::  kn1_uw,kn1_vw,kn2_uw,kn2_vw,kn3_uw,kn3_vw,kn4_uw,kn4_vw
    LOGICAL :: docycle
    
    !Allocate Dynamic Variables
    IF(Average_Numpart(ID_DEPTH).eq.0)THEN
      Average_Value(ID_DEPTH)=0.0
      DO n=1,numpar

        if(settlementon)then
          if ( isSettled(n) ) cycle
        endif
        if(settlementon)then
          if(isStranded(n)) cycle
        endif
        if(mortality)then
          if ( isDead(n) ) cycle
        endif
        if(OpenOceanBoundary)then
          if(isOut(n)) cycle
        endif
        if(OilOn)then
            if(par(n,pStatus) == -2) cycle !stranded
        end if
 
        Xpar = par(n,pX)
        Ypar = par(n,pY)
        Zpar = par(n,pZ)
        CALL setInterp(Xpar,Ypar,n)
        conflict=1
        if(Zgrid)then 
            call getDepth(Xpar,Ypar,n,it,P_depth,Fstlev,conflict)         !--- CL-OGS:   new routine getDepth returning as well the
            if(conflict.ne.1)then
              call handleERROR('gDepth    ',999,conflict,n,     &
                         par(n,pX),par(n,pY),par(n,pZ), &
                         Xpar,Ypar,par(n,pZ),docycle)
               call die(n)
               if(docycle) cycle 

            endif
        else
            klev=1
            P_depth = DBLE(-1.0)* getInterp(Xpar,Ypar,VAR_ID_depth,klev)
            Fstlev=1
        endif     
        Average_Value(ID_DEPTH)=Average_Value(ID_DEPTH)+P_depth
        Average_Numpart(ID_DEPTH)= Average_Numpart(ID_DEPTH)+1                           !--- CL-OGS
      ENDDO
       if((Average_Numpart(ID_DEPTH))>0) then
       Average_Value(ID_DEPTH)=Average_Value(ID_DEPTH)/float(Average_Numpart(ID_DEPTH))
       else
       write(*,*)"In set_Average_Water_Depth Average_Numpart(ID_DEPTH)=0"
       endif
    ENDIF

  END SUBROUTINE  set_Average_Water_Depth


  SUBROUTINE   set_Average_Wind()
    USE PARAM_MOD, ONLY: ui,vi,uj,vj,us,ws,constTemp,constUwind,constVwind, &
                         numpar,idt,Zgrid,Wind,SaltTempOn,WindIntensity,pi, &
                         OilOn,settlementon,mortality,OpenOceanBoundary
    USE HYDRO_MOD, ONLY: WCTS_ITPI,getKRlevel,getDepth, &
                         getSlevel,getWlevel,setInterp,getInterp
    USE INT_MOD,    ONLY: polintd
    use behavior_mod, only: die,isOut,isDead
    USE SETTLEMENT_MOD, ONLY: isSettled,testSettlement,isStranded

#include "VAR_IDs.h"
    IMPLICIT NONE

    INTEGER :: i,deplvl,n
    INTEGER :: klev,Fstlev,NumInterpLvl,ixnum,nklev,k,conflict                     !--- CL:OGS

    ! Particle tracking
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: Pwc_zb,Pwc_zc,Pwc_zf
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: Pwc_wzb,Pwc_wzc,Pwc_wzf
    DOUBLE PRECISION :: Xpar,Ypar,Zpar,newXpos,newYpos,newZpos,P_zb,P_zc,P_zf, &
      P_depth,nP_depth,P_angle,P_zeta,P_zetab,P_zetac,P_zetaf,ey(3),  &
      AdvecUwind,AdvecVwind

    DOUBLE PRECISION :: x1,x2,x3,y1,y2,y3,z1,z2,z3,slope,length,Ttemp

        DOUBLE PRECISION::      P_hsig,P_tm01,P_Uwind,P_Vwind,P_pdir,P_wlen,UWindDrift,&
                                VWindDrift,alpha,PWind,P_Uw,P_Vw,Uadw,Vadw,            &
                                UStokesDrift,VStokesDrift
        DOUBLE PRECISION::  kn1_uw,kn1_vw,kn2_uw,kn2_vw,kn3_uw,kn3_vw,kn4_uw,kn4_vw
    LOGICAL :: docycle
    
    IF(Average_Numpart(ID_U_WIND).eq.0 .and. Wind)THEN
      Average_Value(ID_U_WIND)=0.0
      Average_Value(ID_V_WIND)=0.0
      DO n=1,numpar

        if(settlementon)then
          if ( isSettled(n) ) cycle
        endif
        if(settlementon)then
          if(isStranded(n)) cycle
        endif
        if(mortality)then
          if ( isDead(n) ) cycle
        endif
        if(OpenOceanBoundary)then
          if(isOut(n)) cycle
        endif
        if(OilOn)then
            if(par(n,pStatus) == -2) cycle !stranded
        end if
 
        Xpar = par(n,pX)
        Ypar = par(n,pY)
        Zpar = par(n,pZ)
        CALL setInterp(Xpar,Ypar,n)
        conflict=1
        if(Zgrid)then 
            klev=getKRlevel(par(n,pZ))
            if(P_oldLev(n)>=0) P_oldLev(n)=klev
            call getDepth(Xpar,Ypar,n,it,P_depth,Fstlev,conflict)         !--- CL-OGS:   new routine getDepth returning as well the
            if(conflict.ne.1)then
              call handleERROR('gDepth    ',999,conflict,n,     &
                         par(n,pX),par(n,pY),par(n,pZ), &
                         Xpar,Ypar,par(n,pZ),docycle)
               call die(n)
               if(docycle) cycle 

            endif
            P_angle=0
        else
            klev=1
            P_depth = DBLE(-1.0)* getInterp(Xpar,Ypar,VAR_ID_depth,klev)
            P_angle = getInterp(Xpar,Ypar,VAR_ID_angle,klev)
            Fstlev=1
        endif     

          CALL find_winds(Xpar,Ypar,ex,ix,p,1,Uadw,Vadw,n)
          
          !Store advection currents at original coordinates
          kn1_uw = Uadw
          kn1_vw = Vadw
          
          !Estimate new coordinates for next RK position
          x1 = Xpar + (Uadw*cos(P_angle) - Vadw*sin(P_angle)) * DBLE(idt)/DBLE(2)
          y1 = Ypar + (Uadw*sin(P_angle) + Vadw*cos(P_angle)) * DBLE(idt)/DBLE(2)
          
          !Find advection currents at estimated next RK position
          CALL find_winds(x1,y1,ex,ix,p,2,Uadw,Vadw,n)
          
          !Store advection currents at 2nd RK position
          kn2_uw = Uadw
          kn2_vw = Vadw
          
          !Estimate new coordinates for next RK position
          x2 = Xpar + (Uadw*cos(P_angle) - Vadw*sin(P_angle)) * DBLE(idt)/DBLE(2)
          y2 = Ypar + (Uadw*sin(P_angle) + Vadw*cos(P_angle)) * DBLE(idt)/DBLE(2)
          
          !Find advection currents at estimated next RK position
          CALL find_winds(x2,y2,ex,ix,p,2,Uadw,Vadw,n)
          
          !Store advection currents at 3rd RK position
          kn3_uw = Uadw
          kn3_vw = Vadw
          
          !Calculate the coordinates at the final position
          x3 = Xpar + (Uadw*cos(P_angle) - Vadw*sin(P_angle)) * DBLE(idt)
          y3 = Ypar + (Uadw*sin(P_angle) + Vadw*cos(P_angle)) * DBLE(idt)
          
          
          !Find advection currents at the final position
          CALL find_winds(x3,y3,ex,ix,p,3,Uadw,Vadw,n)
          
          !Store advection currents at final position
          kn4_uw = Uadw
          kn4_vw = Vadw
          
          !Use the RK formula to get the final Advection values
          P_Uw = (kn1_uw + DBLE(2.0)*kn2_uw + DBLE(2.0)*kn3_uw + kn4_uw)/DBLE(6.0)
          P_Vw = (kn1_vw + DBLE(2.0)*kn2_vw + DBLE(2.0)*kn3_vw + kn4_vw)/DBLE(6.0)
         
          !calculate wind vector magnitude
          PWind = sqrt((P_Uw*cos(P_angle) - P_Vw*sin(P_angle))**2.0    & !u rectified to E-W orientation
                      + (P_Uw*sin(P_angle) + P_Vw*cos(P_angle))**2.0)    !v rectified to N-S orientation
!        
       
          Average_Value(ID_U_WIND)=Average_Value(ID_U_WIND)+(P_Uw*cos(P_angle) - P_Vw*sin(P_angle)) 
          Average_Value(ID_V_WIND)=Average_Value(ID_V_WIND)+(P_Uw*sin(P_angle) + P_Vw*cos(P_angle)) 
          Average_Numpart(ID_U_WIND)= Average_Numpart(ID_U_WIND)+1                           !--- CL-OGS

      ENDDO 

      IF(Average_Numpart(ID_U_WIND)>0 .and. abs(Average_Value(ID_U_WIND))+abs(Average_Value(ID_V_WIND))>0.0001)then
         Average_Value(ID_U_WIND)=Average_Value(ID_U_WIND)/float(Average_Numpart(ID_U_WIND))
         Average_Value(ID_V_WIND)=Average_Value(ID_V_WIND)/float(Average_Numpart(ID_U_WIND))
      ELSE
         write(*,'(2(a,f9.4),2a,i4)')'In set_Average_Wind, USING constUWind=',constUWind, &
        ' constVWind=',ConstVWind,' instead of computed AvWind because',&
        ' Average_Numpart(ID_U_WIND)=',Average_Numpart(ID_U_WIND)
         Average_Value(ID_U_WIND) = constUwind
         Average_Value(ID_V_WIND) = constVwind    !v rectified to N-S orientation
      ENDIF


    ENDIF !Average_Numpart(ID_U_WIND)==0 .and. Wind)

  END SUBROUTINE  set_Average_Wind

  SUBROUTINE   set_Average_Temperature()
    USE PARAM_MOD, ONLY: ui,vi,uj,vj,us,ws,constTemp,constUwind,constVwind, &
                         numpar,idt,Zgrid,Wind,SaltTempOn,WindIntensity,pi, &
                         OilOn,settlementon,mortality,OpenOceanBoundary
    USE HYDRO_MOD, ONLY: WCTS_ITPI,getKRlevel,getDepth, &
                         getSlevel,getWlevel,setInterp,getInterp
    USE INT_MOD,    ONLY: polintd
    use behavior_mod, only: die,isOut,isDead
    USE SETTLEMENT_MOD, ONLY: isSettled,testSettlement,isStranded

#include "VAR_IDs.h"
    IMPLICIT NONE

    INTEGER :: i,deplvl,n
    INTEGER :: klev,Fstlev,NumInterpLvl,ixnum,nklev,k,conflict                     !--- CL:OGS

    ! Particle tracking
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: Pwc_zb,Pwc_zc,Pwc_zf
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION( : ) :: Pwc_wzb,Pwc_wzc,Pwc_wzf
    DOUBLE PRECISION :: Xpar,Ypar,Zpar,newXpos,newYpos,newZpos,P_zb,P_zc,P_zf, &
      P_depth,nP_depth,P_angle,P_zeta,P_zetab,P_zetac,P_zetaf,ey(3),  &
      AdvecUwind,AdvecVwind

    DOUBLE PRECISION :: x1,x2,x3,y1,y2,y3,z1,z2,z3,slope,length,Ttemp

        DOUBLE PRECISION::      P_hsig,P_tm01,P_Uwind,P_Vwind,P_pdir,P_wlen,UWindDrift,&
                                VWindDrift,alpha,PWind,P_Uw,P_Vw,Uadw,Vadw,            &
                                UStokesDrift,VStokesDrift
        DOUBLE PRECISION::  kn1_uw,kn1_vw,kn2_uw,kn2_vw,kn3_uw,kn3_vw,kn4_uw,kn4_vw
    LOGICAL :: docycle
    
 
    IF(Average_Numpart(ID_TEMP).eq.0 .and. SaltTempOn)THEN
      Average_Value(ID_TEMP)=0.0
      ALLOCATE(Pwc_zb(us))
      ALLOCATE(Pwc_zc(us))
      ALLOCATE(Pwc_zf(us))
      ALLOCATE(Pwc_wzb(ws))
      ALLOCATE(Pwc_wzc(ws))
      ALLOCATE(Pwc_wzf(ws))
      DO n=1,numpar ! loop for each particle

        if(settlementon)then
          if ( isSettled(n) ) cycle
        endif
        if(settlementon)then
          if(isStranded(n)) cycle
        endif
        if(mortality)then
          if ( isDead(n) ) cycle
        endif
        if(OpenOceanBoundary)then
          if(isOut(n)) cycle
        endif
        if(OilOn)then
            if(par(n,pStatus) == -2) cycle !stranded
        end if
 
          Xpar = par(n,pX)
          Ypar = par(n,pY)
          Zpar = par(n,pZ)
          CALL setInterp(Xpar,Ypar,n)
          conflict=1
          if(Zgrid)then 
            klev=getKRlevel(par(n,pZ))
            if(P_oldLev(n)>=0) P_oldLev(n)=klev
            call getDepth(Xpar,Ypar,n,it,P_depth,Fstlev,conflict)         !--- CL-OGS:   new routine getDepth returning as well the
            P_angle=0
          else
            klev=1
            P_depth = DBLE(-1.0)* getInterp(Xpar,Ypar,VAR_ID_depth,klev)
            P_angle = getInterp(Xpar,Ypar,VAR_ID_angle,klev)
            Fstlev=1
          endif    
          if(Zgrid.and.klev.eq.us)then
              ey(1) =getInterp(Xpar,Ypar,VAR_ID_tempb,klev)
              ey(2) =getInterp(Xpar,Ypar,VAR_ID_tempc,klev)
              ey(3) =getInterp(Xpar,Ypar,VAR_ID_tempf,klev)
              !write(*,*)'Interp Temp',ex,ey,3,ix(2)
              Ttemp= polintd(ex,ey,3,ix(2))
              !write(*,*)'Interp Temp result is ',Ttemp
          else 
             !**********   Sea level zeta  *****************************
             if(p.eq.1)then                        !--- CL-OGS: added in v.Zlev, but is it right ?? 
             P_zetab = getInterp(Xpar,Ypar,VAR_ID_zetab,klev)     !--- CL-OGS: added in v.Zlev, but is it right ?? 
             P_zetac = getInterp(Xpar,Ypar,VAR_ID_zetab,klev)     !--- CL-OGS: added in v.Zlev, but is it right ?? 
             P_zetaf = getInterp(Xpar,Ypar,VAR_ID_zetac,klev)     !--- CL-OGS: added in v.Zlev, but is it right ?? 
             else                                !--- CL-OGS: added in v.Zlev, but is it right ?? 
             P_zetab = getInterp(Xpar,Ypar,VAR_ID_zetab,klev)
             P_zetac = getInterp(Xpar,Ypar,VAR_ID_zetac,klev)
             P_zetaf = getInterp(Xpar,Ypar,VAR_ID_zetaf,klev)
             endif                                 !--- CL-OGS: added in v.Zlev, but is it right ??
            
             !Create matrix of z-coordinates at particle and at each node for
             !  back, center, forward times
             do i=1,us
               !Rho-coordinate depths at particle location
               Pwc_zb(i)=getSlevel(P_zetab,P_depth,i)
               Pwc_zc(i)=getSlevel(P_zetac,P_depth,i)
               Pwc_zf(i)=getSlevel(P_zetaf,P_depth,i)
               !W-coordinate depths at particle location
               Pwc_wzb(i)= getWlevel(P_zetab,P_depth,i)
               Pwc_wzc(i)= getWlevel(P_zetac,P_depth,i)
               Pwc_wzf(i)= getWlevel(P_zetaf,P_depth,i)
             enddo
             !W-coordinate depths at particle location (cont.)
             Pwc_wzb(ws)= getWlevel(P_zetab,P_depth,ws)
             Pwc_wzc(ws)= getWlevel(P_zetac,P_depth,ws)
             Pwc_wzf(ws)= getWlevel(P_zetaf,P_depth,ws)
             do i=Fstlev+2,us-2
               if ((Zpar .LT. Pwc_zb(i)) .OR.    &
                   (Zpar .LT. Pwc_zc(i)) .OR.    &
                   (Zpar .LT. Pwc_zf(i))         ) exit
             enddo
             deplvl = i-2
             if(Zgrid .and. deplvl+3>us)then
                 NumInterpLvl=us-Fstlev+1
                 deplvl=Fstlev
             else
                 NumInterpLvl=4
             endif 
             Ttemp=WCTS_ITPI(VAR_ID_temp,Xpar,Ypar,deplvl,Pwc_zb,Pwc_zc,Pwc_zf,     &
                              us,P_zb,P_zc,P_zf,ex,ix,p,4,n,NumInterpLvl)
           endif
           Average_Value(ID_TEMP)=Average_Value(ID_TEMP)+Ttemp                                            !--- CL-OGS 
           Average_Numpart(ID_TEMP)= Average_Numpart(ID_TEMP)+1                           !--- CL-OGS
      ENDDO !n=1,numpar ! loop for each particle
      
      IF(Average_Numpart(ID_TEMP)>0 .and. abs(Average_Value(ID_TEMP))>0.0001 .and.    &
       abs(Average_Value(ID_TEMP)/float(Average_Numpart(ID_TEMP)))<50.)then
         Average_Value(ID_TEMP)=Average_Value(ID_TEMP)/float(Average_Numpart(ID_TEMP))
      ELSE
         Average_Value(ID_TEMP)=constTemp
         if(SaltTempOn) &
         write(*,'(a,f9.4,2a,i4)')'In set_Average_Temperature, USING constTemp=',constTemp, &
        ' instead of computed Average_Value(ID_TEMP) because',&
        ' Average_Numpart(ID_TEMP)=',Average_Numpart(ID_TEMP)
      ENDIF

      DEALLOCATE(Pwc_zb)
      DEALLOCATE(Pwc_zc)
      DEALLOCATE(Pwc_zf)
      DEALLOCATE(Pwc_wzb)
      DEALLOCATE(Pwc_wzc)
      DEALLOCATE(Pwc_wzf)
    ENDIF
    !------

  END SUBROUTINE  set_Average_Temperature

  !-----------------------------------------------------------------------------

  DOUBLE PRECISION FUNCTION Angle_wrtEast(VecX,VecY)  ! Angle of the direction
      USE PARAM_MOD,  ONLY: pi                        ! of a vector with respect
      implicit none                                   ! to the Eastward direction

      double precision, INTENT(IN) :: VecX,VecY
      double precision :: alpha
      IF(abs(VecX)<10e-10)THEN
            if(VecY.gt.0)then
                alpha = pi/2.0
            elseif(VecY.lt.0)then
                alpha = pi * (3.0/2.0)
            else
                alpha = 0.0
            end if
      ELSE
        alpha=atan( (VecY / VecX ) )
        if(Vecx .gt. 0)then
            if(VecY.gt.0)then
                alpha = alpha
            elseif(VecY.lt.0)then
                alpha = 2.0*pi + alpha                    !was alpha = 2.0*pi + alpha and pi * alpha
            else
                alpha = 0.0
            end if
        else if(Vecx .lt. 0)then
            if(VecY.gt.0)then
                alpha = pi + alpha        !                    was alpha = pi + alpha and 2pi * alpha
            elseif(VecY.lt.0)then
                alpha = pi + alpha
            else
                alpha = pi
            end if
        end if
      ENDIF 
!      alpha = alpha + (17.5*(pi/180.0))
      Angle_wrtEast=alpha

  END FUNCTION Angle_wrtEast

  !-----------------------------------------------------------------------------

end program
