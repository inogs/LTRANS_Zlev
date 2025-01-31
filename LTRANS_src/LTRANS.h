!
!This file defines many variables that are read from the GRID.data and LTRANS.data files
!And also groups the input parameters/variables into namelists
!
!NOTE: variables in a namelist can NOT be dynamic variables!!!!! 
!      dynamic namelists are NOT yet supported in FORTRAN90/95 standard
!
!Basically, this file replaces GRID.inc and LTRANS.inc in previous code
!

!--these are the grid dimension definitions, they should read from GRID.data
!
!--- CL-OGS:   Last Modified on: June 2016 by Celia Laurent - OGS
!--- CL-OGS:   main modifications are:
!--- CL-OGS:   - coupling with MITgcm'Z-grid bathymetry and fields
!--- CL-OGS:     including reading of merged-in-time MITgcm files
!--- CL-OGS:     including 
!--- CL-OGS:   - allow backward-in-time simulation and addition of several parameters to 
!--- CL-OGS:     control the time variable of the simulation and the Runge-Kutta scheme
!--- CL-OGS:   - allow stranding along coasts depending on coast distance
!--- CL-OGS:   - allow running LTRANS using a data file other than LTRANS.data
!--- CL-OGS:   - allow keeping particles at constant depth below the sea surface
!--- CL-OGS:   - writing and reading Adjacent elements binary files to save computational time
!--- CL-OGS:   - integration of OILTRANS developments for oil transport and weathering 
!--- CL-OGS:     implemented in LTRANS by Alan Berry, Marine Institute, Ireland
!--- CL-OGS:     (details of development: search for IMIOM for modifications by Alan Berry)


  INTEGER :: ui                     ! number u points in x direction
  INTEGER :: uj                     !        u           y 
  INTEGER :: vi                     !        v           x
  INTEGER :: vj                     !        v           y

  INTEGER :: rho_nodes              ! number rho nodes   !--- CL-OGS: by k-plane
  INTEGER :: u_nodes                ! number of u nodes  !--- CL-OGS: by k-plane
  INTEGER :: v_nodes                ! number of v nodes  !--- CL-OGS: by k-plane

  INTEGER :: max_rho_elements       ! total number of rho elements in any k-plane
  INTEGER :: max_u_elements         ! total number of u elements in any k-plane
  INTEGER :: max_v_elements         ! total           v

!--- CL-OGS: comment declarations:   !--- CL-OGS:   cancelled rho_elements,u_elements,v_elements 
!  INTEGER :: rho_elements           !--- CL-OGS:   as the number of wet elements became dependant on the vertical level
!  INTEGER :: u_elements             !--- CL-OGS:   substitution by rho_kwele(k), u_kwele(k) and v_kwele(k)             
!  INTEGER :: v_elements             !--- CL-OGS:   declared and computed at each vertical level in hydrodynamic_module 
!group the grid info section in a namelist:

  namelist/gridinfo/ ui, uj, vi, vj,rho_nodes, u_nodes, v_nodes,     &
     &               max_rho_elements,max_u_elements,max_v_elements !, &
!    &               rho_elements, u_elements, v_elements           !--- CL-OGS: commented


!The following used to be in LTRANS.inc:


!*** NUMBER OF PARTICLES ***

  INTEGER :: numpar

  namelist/numparticles/numpar      ! Number of particles



!*** TIME PARAMETERS ***

  REAL             :: days          ! Number of days to run the model
  INTEGER          :: iprint        ! Print interval for LTRANS output (s); 3600 = every hour
  INTEGER          :: dt            ! External time step (duration between hydro model predictions) (s) 
  INTEGER          :: idt           ! Internal (particle tracking) time step (s)
!--- CL-OGS: additional parameters to handle time as written in netcdf output file 
  INTEGER          :: iprinto       ! Initial delay in printing
  INTEGER          :: Ext0          ! Initial time step (s)
!--- CL-OGS: added to give the possibility to centre in time the Runge-Kutta scheme

  namelist/timeparam/ days,iprint,dt,idt,       & 
        Ext0,iprinto                                 !--- CL-OGS: added 



!!*** HYDRODYNAMIC MODULE PARAMETERS ***
  INTEGER          :: us            ! Number of Rho grid s-levels in hydro model
  INTEGER          :: ws            ! Number of W grid s-levels in hydro model
  INTEGER          :: tdim          ! Number of time steps per hydro predictions file (=0 if all timesteps are in the same hydro file)
!--- CL-OGS: recordnum added to decide the first record to be read in the case of MITgcm input files merged in time
  INTEGER          :: recordnum       ! first record number to be read (when no previous record: recordnum=1) 
  REAL             :: hc            ! Min Depth - used in ROMS S-level transformations
  DOUBLE PRECISION :: z0            ! ROMS roughness parameter
  INTEGER          :: Vtransform    ! 1-WikiRoms Eq. 1 ; 2-WikiRoms Eq. 2 ; 3-Song/Haidvogel 1994 Eq. ; 0-Zcoordinate system (MITgcm)
  LOGICAL          :: VInterpUVinSurfWater  ! If False keep U,V velocities of last level in surface water instead of running the vertical interpolation 
  DOUBLE PRECISION :: BottomLayerThickness           ! If > 0.0 this value will be used instead of the half height of the bottom grid cell
  DOUBLE PRECISION :: PercentVelinBottomLayer        ! If > 0.0 this value will be used instead of the logarithmic law (must be in range ]0.0,1.0[) between z0 and the bottom layer thickness
  DOUBLE PRECISION :: PercentVel_under_z0            ! where :  0.0 >= PercentVel_under_z0 >= 1.0 : defines the percentage of the velocity used between the bottom and z0
  LOGICAL          :: readZeta      ! If .TRUE. read in sea-surface height   (zeta) from NetCDF file, else use constZeta
  DOUBLE PRECISION :: constZeta     ! Constant value for Zeta if readZeta is .FALSE.
  LOGICAL          :: readSalt      ! If .TRUE. read in salinity             (salt) from NetCDF file, else use constSalt
  DOUBLE PRECISION :: constSalt     ! Constant value for Salt if readSalt is .FALSE.
  LOGICAL          :: readTemp      ! If .TRUE. read in temperature          (temp) from NetCDF file, else use constTemp
  DOUBLE PRECISION :: constTemp     ! Constant value for Temp if readTemp is .FALSE.
  LOGICAL          :: readU         ! If .TRUE. read in u-momentum component (U   ) from NetCDF file, else use constU
  DOUBLE PRECISION :: constU        ! Constant value for U if readU is .FALSE.
  LOGICAL          :: readV         ! If .TRUE. read in v-momentum component (V   ) from NetCDF file, else use constV
  DOUBLE PRECISION :: constV        ! Constant value for V if readV is .FALSE.
  LOGICAL          :: readW         ! If .TRUE. read in w-momentum component (W   ) from NetCDF file, else use constW
  DOUBLE PRECISION :: constW        ! Constant value for W if readW is .FALSE.
  LOGICAL          :: readAks       ! If .TRUE. read in salinity vertical diffusion coefficient (Aks ) from NetCDF file, else use constAks
  DOUBLE PRECISION :: constAks      ! Constant value for Aks if readAks is .FALSE.
  LOGICAL          :: readDens
  DOUBLE PRECISION :: constDens
  LOGICAL          :: readChl         ! If .TRUE. read in Chlorophyl from NetCDF file, else use constChl
  DOUBLE PRECISION :: constChl        ! Constant value for Chl if readU is .FALSE.
!--- CL-OGS: added reading of wind in hydrodynamic_module
  LOGICAL          :: readUwind     ! If .TRUE. read in Uwind from NetCDF file, else use constUwind
  DOUBLE PRECISION :: constUwind    ! Constant value for Uwind if readUwind is .FALSE.
  LOGICAL          :: readVwind     ! If .TRUE. read in Vwind from NetCDF file, else use constVwind
  DOUBLE PRECISION :: constVwind    ! Constant value for Vwind if readVwind is .FALSE.
  LOGICAL          :: readIwind      ! If .TRUE. read in Iwind from NetCDF file, else use constIwind
  DOUBLE PRECISION :: constIwind     ! Constant value for Iwind if readIwind is .FALSE.
  LOGICAL :: Wind            ! Include wind drift effects (T/F)
  LOGICAL :: WindIntensity   ! For Oil weathering use input Iwind intensity instead of calculating it from Uwind,Vwind
  DOUBLE PRECISION :: WindWeatherFac ! Wind Reduction Factor (between 0.0 and 1.0) applied to Oil Spill WEATHERING only 
  DOUBLE PRECISION :: WindDriftFac ! Wind Drift factor, initially 0.035 in OILTRANS 
  DOUBLE PRECISION :: WindDriftDev ! Wind Drift deviation in degrees (offset to RHS of wind vector), initially 5.0 in OILTRANS
  DOUBLE PRECISION :: Wind_hc      ! Wind_hc and Wind_ke used to create a vertical profile of wind drift with 100% wind drift in depth range [-hc, 0] m,
  DOUBLE PRECISION :: Wind_ke      ! while bellow hc the exponential decay is ruled by Wind_ke following WindDrift = WindDrift0 * exp( |Wind_ke| * min( 0 , |Wind_hc|-Z) ) 
  LOGICAL          :: LinearVInterp  ! if True uses linear interpolation scheme on the vertical direction instead of the tension spline fitting.
  DOUBLE PRECISION :: StokDriftFac ! Stokes Drift factor, initially 0.016 in OILTRANS 
  LOGICAL :: Stokes             ! Include Stokes drift effects (T/F)
!
  namelist/hydroparam/us,ws,tdim,hc,z0,Vtransform,readZeta,constZeta,readSalt,   &
                    & constSalt,readTemp,constTemp,readU,constU,readV,           & !--- CL-OGS: cancelled readU which was mentionned twice
                    & constV,readW,constW,readAks,constAks,readDens,constDens,   &
                    & readUwind,constUwind,readVwind,constVwind,readChl,ConstChl, & !--- CL-OGS: additional variables
                    & recordnum,Wind,WindDriftFac,WindDriftDev,StokDriftFac,     &  !--- CL-OGS: moved here from oilprocs
                    & VInterpUVinSurfWater,LinearVInterp,             &  !--- CL-OGS
                    & BottomLayerThickness,PercentVelinBottomLayer,              &   !--- CL-OGS
                    & PercentVel_under_z0,  &  !--- CL-OGS
                    & readIwind,constIwind,WindIntensity,WindWeatherFac,         &  !--- CL-OGS
                    & Wind_hc, Wind_ke, &
                    & Stokes                                                        !--- CL-OGS

!*** TURBULENCE MODULE PARAMETERS ***
  LOGICAL          :: HTurbOn       ! Horizontal Turbulence on (.TRUE.) or off (.FALSE.)
  LOGICAL          :: VTurbOn       ! Vertical   Turbulence on (.TRUE.) or off (.FALSE.)
  DOUBLE PRECISION :: ConstantHTurb ! Constant value of horizontal turbulence (m2/s)

  namelist/turbparam/HTurbOn,VTurbOn,ConstantHTurb

!*** BEHAVIOR MODULE PARAMETERS ***
  INTEGER :: Behavior               ! Behavior type (specify a number)
                                    !   Note: The behavior types numbers are: 
                                    !     0 Passive, 1 near-surface, 2 near-bottom, 3 DVM, 
                                    !     4 C.virginica oyster larvae, 5 C.ariakensis oyster larvae, 
                                    !     6 constant, 7 Tidal Stream Transport !)
                                    !--- CL-OGS:           999 = CstDepth below the (varying) sea surface,
                                                           !particles are not advected by vertical currents 
                                                           ! and are sensible only to vertical turbulence,
                                                           ! vertical displacement due to the variations of the sea surface height
                                                           ! or displacements due to the raising of the bottom sea depth 
                                    ! *****  IMIOM *****  1000 = oil particles
  LOGICAL :: BottomRelease          ! Set to .true. to release particles 1m above the bottom, .False. by default 
  LOGICAL :: OpenOceanBoundary      ! Note: If you want to allow particles to "escape" via open ocean 
                                    !   boundaries, set this to TRUE; Escape means that the particle 
                                    !   will stick to the boundary and stop moving
  LOGICAL :: mortality              ! TRUE if particles can die; else FALSE
  DOUBLE PRECISION :: deadage       ! Age at which a particle stops moving (i.e., dies) (s)
                                    !   Note: deadage stops particle motion for all behavior types (0-6)
  DOUBLE PRECISION :: pediage       ! Age when particle reaches max swim speed and can settle (s)
                                    !   Note: for oyster larvae behavior (types 4 & 5):
                                    !     pediage = age at which a particle becomes a pediveliger
                                    !   Note: pediage does not cause particles to settle if the Settlement module is not on
  DOUBLE PRECISION :: swimstart     ! Age that swimming or sinking begins (s) 1 day = 1.*24.*3600.
  DOUBLE PRECISION :: swimslow      ! Swimming speed when particle begins to swim (m/s)
  DOUBLE PRECISION :: swimfast      ! Maximum swimming speed (m/s)  0.05 m/s for 5 mm/s
                                    !   Note: for constant swimming speed for behavior types 1,2 & 3: 
                                    !     set swimslow = swimfast = constant speed
  DOUBLE PRECISION :: Sgradient     ! Salinity gradient threshold that cues larval behavior (psu/m)
                                    !   Note: This parameter is only used if Behavior = 4 or 5. 
  DOUBLE PRECISION :: rise          ! Rising velocity for behavior type 8
  DOUBLE PRECISION :: sink          ! Sinking velocity for behavior type 6 and 8
                                    !   Note: This parameter is only used if Behavior = 6. or 8.
! Tidal Stream Transport behavior type:
  DOUBLE PRECISION :: Hswimspeed    ! Horizontal swimming speed (m/s)
  DOUBLE PRECISION :: Swimdepth     ! Depth at which fish swims during flood time 
                                    ! in meters above bottom (this should be a positive value)

  CHARACTER(LEN=200) :: GrainSize_fname ! for Behavior 8
  LOGICAL            :: read_GrainSize  ! for Behavior 8
  DOUBLE PRECISION :: surflayer_upperdepth  ! minimal depth of the surface layer for Behavior 9 particles
  DOUBLE PRECISION :: surflayer_lowerdepth  ! Depth of the surface layer for Behavior 8 and 9 particles
  DOUBLE PRECISION :: surflayer_upperdepth_night  ! minimal depth of the surface layer for Behavior 9 particles
  DOUBLE PRECISION :: surflayer_lowerdepth_night  ! Depth of the surface layer for Behavior 8 and 9 particles
  DOUBLE PRECISION :: DVMtime
  DOUBLE PRECISION :: SettlementSize
  LOGICAL          :: SeabedRelease
  DOUBLE PRECISION :: SeabedRelease_meters
  DOUBLE PRECISION :: Seabed_layerheight
  CHARACTER(LEN=3) :: Write_Temp_min_max_ins ! values 'min', 'max', or 'ins' per "instantaneous"
  CHARACTER(LEN=3) :: Write_Salt_min_max_ins ! values 'min', 'max', or 'ins' per "instantaneous"
  LOGICAL          :: Write_coastdist        ! If True write closest distance to coast in netcdf output file
  namelist/behavparam/Behavior,OpenOceanBoundary, &
           mortality,deadage,pediage,swimstart,   & 
           swimslow,swimfast,Sgradient,sink,      &
           Hswimspeed,Swimdepth,                  &
           GrainSize_fname,read_GrainSize,rise,   & !--- CL-OGS additional parameters for behavior 8
           surflayer_upperdepth,surflayer_lowerdepth,  &  !--- CL-OGS additional parameters for behavior 8 and 9 
           surflayer_lowerdepth_night,surflayer_upperdepth_night,     &
           DVMtime,SettlementSize,  &
           SeabedRelease,SeabedRelease_meters,Seabed_layerheight, &
           Write_Temp_min_max_ins,Write_Salt_min_max_ins,Write_coastdist, &
          BottomRelease


!*** DVM. The following are parameters for the Diurnal Vertical Migration (DVM) behavior type:
  DOUBLE PRECISION :: twistart      ! Time of twilight start (hr) **
  DOUBLE PRECISION :: twiend        ! Time of twilight end (hr) **
  DOUBLE PRECISION :: daylength     ! Length of day (hr) **
  DOUBLE PRECISION :: Em            ! Irradiance at solar noon (microE m^-2 s^-1) **
  DOUBLE PRECISION :: Kd            ! Vertical attenuation coefficient
  DOUBLE PRECISION :: thresh        ! Light threshold that cues behavior (microE m^-2 s^-1)
  !  Note: These values were calculated for September 1 at the latitude of 37.0 (Chesapeake Bay mouth)
  !  Note: Variables marked with ** were calculated with light_v2BlueCrab.f (not included in LTRANS yet)
  !  Note: These parameters are only used if Behavior = 3 
  ! --- CL-OGS: following parameters for dvm from MITgcm-STDOUT-swdown averages used for behavior 8 and 10
  LOGICAL          :: swdown_ASCII      ! FLAG to activate/disactivate dial vertical migration from swdown
  CHARACTER(LEN=200) :: swdown_ASCIIfname  ! swdown file name (txt file, one value per line corrisponding to the average value over the whole basin) 
  INTEGER          :: swdown_dt       ! swdown external time step (duration between swdown prints/lines    ) (s) 
  INTEGER          :: swdown_rec      ! swdown first line/record to be read (min 1)                              
  DOUBLE PRECISION :: swdown_thresh   ! swdown threshold for dial vertical migration
  DOUBLE PRECISION :: swdown_t0       ! time in seconds (same referential as Ext0) of the first record indicated as swdown_rec 

  namelist/dvmparam/twistart,twiend,daylength,Em,Kd,thresh,&
     swdown_ASCII,swdown_ASCIIfname,swdown_dt,swdown_rec,swdown_thresh,swdown_t0



!*** SETTLEMENT MODULE PARAMETERS ***
  LOGICAL :: settlementon           ! settlement module on (.TRUE.) or off (.FALSE.)
  !  Note: If settlement is off: set minholeid, maxholeid, minpolyid, maxpolyid, pedges, & hedges to 1
  !        to avoid both wasted variable space and errors due to arrays of size 0.
  !        If settlement is on and there are no holes: set minholeid, maxholeid, & hedges to 1
  LOGICAL :: holesExist             ! Are there holes in habitat? yes(TRUE) no(FALSE)
  INTEGER :: minpolyid              ! Lowest habitat polygon id number
  INTEGER :: maxpolyid              ! Highest habitat polygon id number
  INTEGER :: minholeid              ! Lowest hole id number
  INTEGER :: maxholeid              ! Highest hole id number
  INTEGER :: pedges                 ! Number of habitat polygon edge points (# of rows in habitat polygon file)
  INTEGER :: hedges                 ! Number of hole edge points (number of rows in holes file)
  LOGICAL :: Write_Poly_Presence  
  INTEGER :: storedincolor  !--- CL-OGS : if ==1: If particle is in a polygon, store the poly-number in color array
                            !--- CL-OGS : if ==2: store the element number in which is the particle in color array
                            !--- CL-OGS : if ==0: store Status in color array
  namelist/settleparam/settlementon,holesExist,minpolyid,maxpolyid,minholeid,maxholeid,pedges,hedges, &
          storedincolor,Write_Poly_Presence   !--- CL-OGS : additional parameters

!*** STRANDING MODULE PARAMETERS ***
  LOGICAL :: stranding_on           ! stranding module on (.TRUE.) or off (.FALSE.)
!--- CL-OGS : allow stranding along coast when particle enters the water squared elements bordering the coast
  DOUBLE PRECISION :: StrandingDist          !--- CL-OGS : maximal distance from land at which stranding can occur (in meters)
  DOUBLE PRECISION :: strandingMaxDistFromSurf          !--- CL-OGS : maximal depth (distance from surface) at which stranding can occur                  
  DOUBLE PRECISION :: StrandingMaxDistFromBott  !--- CL-OGS : maximal distance from bottom depth at which stranding can occur                 
  namelist/strandingparam/stranding_on,StrandingDist,strandingMaxDistFromSurf,StrandingMaxDistFromBott  !--- CL-OGS : additional parameters

!*** CONVERSION MODULE PARAMETERS ***
  DOUBLE PRECISION :: PI            ! Pi
  DOUBLE PRECISION :: Earth_Radius  ! Equatorial radius
  LOGICAL :: SphericalProjection    ! Spherical Projection from ROMS (T/F)
  DOUBLE PRECISION :: lonmin        ! minimum longitude value, only used if SphericalProjection is .TRUE.
  DOUBLE PRECISION :: latmin        ! minimum  latitude value, only used if SphericalProjection is .TRUE.

  namelist/convparam/PI,Earth_Radius,SphericalProjection,lonmin,latmin



!*** INPUT FILE NAME AND LOCATION PARAMETERS ***; 

!  ** NetCDF Hydro Model Grid file **
  CHARACTER(LEN=200) :: NCgridfile
    !Note: the path to the file is necessary if the file is not in the same folder as the code
    !Note: if .nc file in separate folder in Linux, then include path. For example:
    !      CHARACTER(LEN=29), PARAMETER :: NCgridfile = '/share/enorth/CPB_GRID_wUV.nc' 
    !Note: if .nc file in separate folder in Windows, then include path. For example:
    !      CHARACTER(LEN=23), PARAMETER :: NCgridfile = 'D:\ROMS\CPB_GRID_wUV.nc'
!--- CL-OGS : main keyword to allow the model to read MITgcm Z-grid files and create 3d grid water element numbering
  LOGICAL            :: Zgrid    ! If .TRUE. read in MITgcm z-coordinate grids+fields, else use ROMS s-coordinate grids+fields
  LOGICAL            :: Zgrid_depthinterp    ! If .TRUE. interpolates depth to make the Z grid bathymetry smoother (for Zgrid bathymetry only)
!--- CL-OGS : option implemented to save computational time when running multiple simulations on the same grid:
!--- CL-OGS : at the first run give a name to ADJele_fname and set ADJele_file=False 
!--- CL-OGS :            -> creates a binary file containing the Adjacent elements list computed by the boundary module
!--- CL-OGS : for all successive runs give the name to ADJele_fname and set ADJele_file=True 
!--- CL-OGS :            -> reads the previously created binary file containing the Adjacent elements list 
!--- CL-OGS :               instead of re-computing them in the boundary module
  CHARACTER(LEN=200) :: ADJele_fname
  LOGICAL            :: ADJele_file  



  namelist/hydromodelgrid/NCgridfile, &
           Zgrid,Zgrid_depthinterp,ADJele_fname,ADJele_file              !--- CL-OGS additional parameters




!  ** Hydro Model Predictions NetCDF Input File **
!  Filename = dirin + prefix + filenum + suffix
  CHARACTER(LEN=200) :: dirin       !--- CL-OGS : for MITgcm files : ! NetCDF Input Directory
!--- CL-OGS : for MITgcm prefix contains the file names of the 10 NetCDF field variables without directory neither suffix
  CHARACTER(LEN=200) :: prefix_Zeta       ! NetCDF Input Filename prefix Zeta (Sea level height) 
  CHARACTER(LEN=200) :: prefix_Salt       ! NetCDF Input Filename prefix Salt  
  CHARACTER(LEN=200) :: prefix_Temp       ! NetCDF Input Filename prefix Temp  
  CHARACTER(LEN=200) :: prefix_Uvel       ! NetCDF Input Filename prefix Uvel (Current velocity in X direction)
  CHARACTER(LEN=200) :: prefix_Vvel       ! NetCDF Input Filename prefix Vvel (Current velocity in X direction) 
  CHARACTER(LEN=200) :: prefix_Wvel       ! NetCDF Input Filename prefix Wvel (Current velocity in X direction)
  CHARACTER(LEN=200) :: prefix_Aks        ! NetCDF Input Filename prefix Aks  (salinity vertical diffusion coefficient) 
  CHARACTER(LEN=200) :: prefix_Dens       ! NetCDF Input Filename prefix Dens (Sea water density)
  CHARACTER(LEN=200) :: prefix_Chl        ! NetCDF Input Filename prefix Chl  (chlorophyl) 
  CHARACTER(LEN=200) :: prefix_Uwind      ! NetCDF Input Filename prefix Uwind (Wind in X direction)
  CHARACTER(LEN=200) :: prefix_Vwind      ! NetCDF Input Filename prefix Vwind (Wind in y direction)
  CHARACTER(LEN=200) :: prefix_Iwind      ! NetCDF Input Filename prefix Iwind (Wind intensity)
  CHARACTER(LEN=200) :: namevar_Zeta       ! NetCDF Input Filename prefix Zeta (Sea level height) 
  CHARACTER(LEN=200) :: namevar_Salt       ! NetCDF Input Filename prefix Salt  
  CHARACTER(LEN=200) :: namevar_Temp       ! name of var Temp  
  CHARACTER(LEN=200) :: namevar_Uvel       ! name of var Uvel  in NetCDF input file (Current velocity in X direction)
  CHARACTER(LEN=200) :: namevar_Vvel       ! name of var Vvel  in NetCDF input file (Current velocity in X direction) 
  CHARACTER(LEN=200) :: namevar_Wvel       ! name of var Wvel  in NetCDF input file (Current velocity in X direction)
  CHARACTER(LEN=200) :: namevar_Aks        ! name of var Aks   in NetCDF input file (salinity vertical diffusion coefficient) 
  CHARACTER(LEN=200) :: namevar_Dens       ! name of var Dens  in NetCDF input file (Sea water density)
  CHARACTER(LEN=200) :: namevar_Chl        ! name of var Chl   in NetCDF input file (Chlorophyl)
  CHARACTER(LEN=200) :: namevar_Uwind      ! name of var Uwind in NetCDF input file (Wind in X direction)
  CHARACTER(LEN=200) :: namevar_Vwind      ! name of var Vwind in NetCDF input file (Wind in y direction)
  CHARACTER(LEN=200) :: namevar_Iwind      ! name of var Iwind in NetCDF input file (Wind intensity)
  CHARACTER(LEN=200) :: suffix      ! NetCDF Input Filename suffix
  INTEGER :: filenum                ! Number in First NetCDF Input Filename
  INTEGER :: filestep               !--- CL-OGS : Number between successive NetCDF input filename needed for MITgcm files 
!--- CL-OGS : with MITgcm merged-in-time files : numdigits=0 
  INTEGER :: numdigits              ! Number of digits in number portion of file name (with leading zeros) 
  INTEGER :: hydrobytes             ! Number of bytes of float variables: single =4 bytes, double = 8 bytes
  LOGICAL :: startfile              ! .TRUE. means the first file has an additional time step
  LOGICAL :: Hydro_NetCDF           ! .TRUE. means the Hydro model prediction Input files are in NetCDF format
  LOGICAL :: First_vertical_layer_is_surface ! .TRUE.  means that in the netcdfs hydro files the first vertical layer is surface (for MITgcm netcdfs)
                                             ! .FALSE. means that in the netcdfs hydro files the first vertical layer is bottom (for ROMS netcdfs)  
  !Note: the path to the file is necessary if the file is not in the same folder as the code
  !Note: if .nc file in separate folder in Windows, then include path in prefix. For example:
  !      CHARACTER(LEN=15), PARAMETER :: prefix='D:\ROMS\y95hdr_'   
  !      if .nc file in separate folder in Linux, then include path in prefix. For example:
  !      CHARACTER(LEN=26), PARAMETER :: prefix='/share/lzhong/1995/y95hdr_'   

  namelist/hydromodeloutput/prefix_Zeta,prefix_Salt,prefix_Temp,prefix_Uvel,prefix_Vvel, &
                            prefix_Wvel,prefix_Aks,prefix_Dens,prefix_Uwind,prefix_Vwind, &
                            prefix_Iwind,prefix_Chl,suffix,filenum,numdigits,startfile,              &
                            namevar_Zeta,namevar_Salt,namevar_Temp,namevar_Uvel,namevar_Vvel,  &
                            namevar_Wvel,namevar_Aks,namevar_Dens,namevar_Uwind,namevar_Vwind, &
                            namevar_Iwind,namevar_Chl,Hydro_NetCDF,First_vertical_layer_is_surface,       &
                            filestep,dirin,hydrobytes  !--- CL-OGS additional parameters


!  ** Particle Location Input File **
  CHARACTER(LEN=200) :: parfile     ! Particle locations file
  !Note: the path to the file is necessary if the file is not in the same folder as the code

  namelist/parloc/parfile



!  ** Habitat Polygon Location Input Files **
  CHARACTER(LEN=200) :: habitatfile ! Habitat polygon file
  CHARACTER(LEN=200) :: holefile    ! Holes in habitat file
  !Note: the path to the file is necessary if the file is not in the same folder as the code

  namelist/HabPolyLoc/habitatfile,holefile



!  ** Output Related Variables **
  CHARACTER(LEN=200) :: outpath     ! Location to write output files
  CHARACTER(LEN=200) :: NCOutFile   ! Name of the NetCDF output file if outputting to NetCDF
  LOGICAL :: outpathGiven           ! If TRUE files are written to the path given in outpath
  LOGICAL :: writeCSV               ! If TRUE write CSV output files
  LOGICAL :: writePARA              !--- CL-OGS : ! If TRUE write PARA.CSV output files
  LOGICAL :: writeNC                ! If TRUE write .NC output files
  INTEGER :: NCtime                 ! Time interval between creation of new NetCDF output files

  !NetCDF Model Metadata:
  CHARACTER(LEN=200) :: SVN_Version ! SVN Repository and Version #
  CHARACTER(LEN=100) :: RunName     ! Unique Identifier for this particular model run
  CHARACTER(LEN=200) :: ExeDir      ! Location of the model run executable
  CHARACTER(LEN=200) :: OutDir      ! Location of the model run output files
  CHARACTER(LEN=100) :: RunBy       ! Name of person who setup/run the model
  CHARACTER(LEN=100) :: Institution ! Place the model is run
  CHARACTER(LEN=200) :: StartedOn   ! Date the model run began

  namelist/output/outpath,NCOutFile,outpathGiven,writeCSV,writeNC,NCtime,        &
                  SVN_Version,RunName,ExeDir,OutDir,RunBy,Institution,StartedOn, &
                  writePARA             !--- CL-OGS added 



!*** OTHER PARAMETERS *** 
  INTEGER :: seed                   ! Seed value for random number generator (Mersenne Twister)
  INTEGER :: ErrorFlag              ! What to do if an error is encountered: 0=stop, 1=return particle to previous location
                                    ! 2=kill particle & stop tracking, 3=set particle out of bounds & stop tracking
                                                      ! Note: Options 1-3 will output information to ErrorLog.txt
  LOGICAL :: BndOut                 ! Write Bounds.out file 
  LOGICAL :: BoundaryBLNs           ! Create Surfer Blanking Files of boundaries? .TRUE.=yes, .FALSE.=no
  LOGICAL :: SaltTempOn             ! Calculate salinity and temperature at particle 
                                    ! location: yes (.TRUE.) or no (.FALSE.)
  LOGICAL :: TrackCollisions        ! Write Bottom and Land Hit Files? .TRUE.=yes, .FALSE.=no
  LOGICAL :: WriteHeaders           ! Write .txt files with column headers? .TRUE.=yes, .FALSE.=no
  LOGICAL :: WriteModelTiming       ! Write .csv file with model timing data? .TRUE.=yes, .FALSE.=no
  LOGICAL :: WriteProblemFile       ! Write a file with problem particles for debugging? .TRUE.=yes, .FALSE.=no
  LOGICAL :: WriteCurrents          ! Write .csv file with particle locations for each external timestep ? .TRUE.=yes, .FALSE.=no                     
  LOGICAL :: WriteParfile           ! Write .csv file with particle advection currents and locations at each internal timestep? .TRUE.=yes, .FALSE.=no

  INTEGER :: ijbuff                 ! number of extra elements to read in on every side of the particles

  LOGICAL :: FreeSlip               ! use free slip condition?
  LOGICAL :: OilOn                ! *****   IMIOM ***** : Irish Marine Institute Oil Model Selector

  namelist/other/seed,BoundaryBLNs,SaltTempOn,TrackCollisions,WriteHeaders, &
                 WriteModelTiming,WriteProblemFile,ijbuff,ErrorFlag,FreeSlip, &
                 WriteCurrents,WriteParfile,BndOut, &       !--- CL-OGS
                 OilOn                                      ! *****   IMIOM *****
                
!*******************************************************
!*** Marine Institute Oil Model Additional Variables ***
!*******************************************************

!*** OIL MODEL SIMULATION PARAMETERS
DOUBLE PRECISION :: VolumeSpill
INTEGER :: SecSpill
DOUBLE PRECISION :: API
DOUBLE PRECISION :: Oil_Dens
DOUBLE PRECISION :: Oil_Dens_RefT
DOUBLE PRECISION :: WaterTemp
DOUBLE PRECISION :: Dyn_Visc
DOUBLE PRECISION :: Kin_Visc
DOUBLE PRECISION :: Dyn_Visc_RefT
DOUBLE PRECISION :: Kin_Visc_RefT
CHARACTER(LEN=6) :: Cut_Unit
DOUBLE PRECISION :: Oil_Asph
DOUBLE PRECISION :: Oil_Resin
DOUBLE PRECISION :: Oil_Sat
DOUBLE PRECISION :: Fingas_B
DOUBLE PRECISION :: Fingas_T
INTEGER :: Fingas_TYP
DOUBLE PRECISION, DIMENSION(15) :: Cut_Temp
DOUBLE PRECISION, DIMENSION(15) :: Cut_Frac

namelist/oilparams/VolumeSpill,SecSpill,API,Oil_Dens,Oil_Dens_RefT,WaterTemp,  &
                              Dyn_Visc,Kin_Visc,Dyn_Visc_RefT,Kin_Visc_RefT,&
                              Cut_Unit,Oil_Asph, Oil_Resin,Oil_Sat,Fingas_B, Fingas_T, &
                              Fingas_TYP,Cut_Temp,Cut_Frac



!*** OIL PROCESS OPTIONS
LOGICAL :: Spreading                    ! Include gravity-viscous spreading (T/F)
CHARACTER(LEN=6) :: AreaOption        ! Choose from ADIOS2 / MOHID2 / CONCAW / OILPOL
CHARACTER(LEN=6) :: SprdOption    ! Choose from ADIOS2 / MOHID2 / CONCAW / OILPOL
LOGICAL :: Evaporation                    ! Include evaporation (T/F)
CHARACTER(LEN=6) :: EvapOption    ! Choose from: FINGAS / MACKAY / PSEUDO
LOGICAL :: Emulsification              ! Include emulsification (T/F)
LOGICAL :: Dispersion                    ! Include dispersion due to wave breaking (T/F)
LOGICAL :: Remove_Stranded_Oil            ! Stranded oil is removed from weathering oil               
LOGICAL :: Langmuir                            ! Include Langmuir circulation effects (T/F)

namelist/oilprocs/Spreading,AreaOption,SprdOption,Evaporation,EvapOption, &
                  Emulsification,Dispersion,Remove_Stranded_Oil,Langmuir


!*** WAVE MODEL DATA
LOGICAL :: WindWaveModel
CHARACTER(LEN=200) :: swan_prefix
INTEGER :: swan_filenum
CHARACTER(LEN=3) :: swan_suffix
DOUBLE PRECISION :: SigWaveHeight
DOUBLE PRECISION :: SigWavePeriod
DOUBLE PRECISION :: SigWaveLength
DOUBLE PRECISION :: MeanWavePeriod
DOUBLE PRECISION :: PeakDirection
DOUBLE PRECISION :: PeakWaveLength
DOUBLE PRECISION :: MixingDepth
DOUBLE PRECISION :: Cd
DOUBLE PRECISION :: Disper

namelist/windswaves/WindWaveModel,swan_prefix,swan_filenum,swan_suffix,SigWaveHeight, &
                  SigWavePeriod,SigWaveLength,MeanWavePeriod,     &
                  PeakDirection,PeakWaveLength,MixingDepth,Cd,Disper

