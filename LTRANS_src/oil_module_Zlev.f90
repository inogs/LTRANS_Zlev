MODULE OIL_MOD
! ********************************************************************
! *                                                                  *
! *             Irish Marine Institute Oil Model Module              *
! *                                                                                                          *
! * Developer: Alan Berry                                             *
! * Release  : v0.1                                                   *
! * Date     : July 2011                                             *
! *                                                                  *
! ********************************************************************
!  This module handles all the Marine Institute Oil Model formulations
!
!  Created by:            Alan Berry
!  Created on:            01 July 2011
!  Last Modified on:      01 July 2011MODULE OIL_MOD
!
! For referencing please use:
!
! Berry, A., Dabrowski, T., Lyons, K., 2012. The oil spill model
! OILTRANS and its application to the Celtic Sea. Marine Pollution 
! Bulletin, 64(11) : 2489-2501.
!
! The OILTRANS model was developed as part of ARCOPOL – The Atlantic
! Regions’ Coastal Pollution Response Project (contract nr.2008-1/061)
! funded by the Atlantic Area Trans-National Programme (Priority 2: 
! Marine Environment and Renewable Energy) with support by the European
! Regional Development Fund (ERDF).
!
!
!
! **********************************************************************
! **********************************************************************
! **                      Copyright (c) 2016                          **
! **               The Marine Institute, Ireland                      **
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
! ** The most current official versions of this Software and          **
! ** associated tools and documentation are available from the        ** 
! ** authors by e-mail:                                               **
! **                                                                  **
! **      ocean.modelling@marine.ie                                   **
! **                                                                  **
! ** We ask that users make appropriate acknowledgement of            **
! ** The Marine Institute, Ireland,                                   **
! ** individual developers, participating agencies and institutions,  **
! ** and funding agencies. One way to do this is to cite              **
! ** the following publication:                                       **
! **                                                                  **
! ** Berry, A., Dabrowski, T., Lyons, K., 2012. The oil spill model   **
! ** OILTRANS and its application to the Celtic Sea. Marine Pollution **
! ** Bulletin, 64(11) : 2489-2501.                                    **
! **                                                                  **
! **********************************************************************
! **********************************************************************



      IMPLICIT NONE
      SAVE


    !General
    ! OGS for LTRANS-Zlev: uncommented the next five lines
    DOUBLE PRECISION :: MassOil, MassSpill, MassEvap, MassDisp, RhoOil,        &
                        DeltaRho, AreaOil,ViscOil, AsphOil, ResinOil, SatOil,  &
                        ThickLimit, lenR, OilThickness,SlickThickness,         &
                        VolumeOil, VolumeSlick, VolumeEvap, WaterContent,      &
                        xviscemul, VolumeDisp, OilDensity, OilDensity_RefT,    &
                        AvWindAngle,VolumeBeached,MassBeached, & !,VolumeDiss   ! Added by OGS
                        MassEvapPrev,MassDispPrev,MassBeachedPrev
    !Timing variables 
    INTEGER :: Phase1Time
    LOGICAL :: FIRSTAP
    INTEGER :: iT                                            !internal timestep counters

    !*** OIL CONSTANTS
    DOUBLE PRECISION, PARAMETER :: CDensT = 8.0E-04          !NOAA
    DOUBLE PRECISION, PARAMETER :: CDensE = 0.18            !Mackay
    DOUBLE PRECISION, PARAMETER :: Emul_C4 = 2.0E-06        !Reed 1998 Oil & Chem Poll, after Mackay 1982
    DOUBLE PRECISION, PARAMETER :: Emul_C3 = 0.70            !Reed 1998 Oil & Chem Poll, after Mackay 1982
    DOUBLE PRECISION, PARAMETER :: Emul_C1 = 0.65            !Mackay 1982, after Mooney 1951
    DOUBLE PRECISION, PARAMETER :: Evap_C4 = 10.0            !Reed 1998 Oil & Chem Poll, after Mackay 1982
    DOUBLE PRECISION, PARAMETER :: Disp_Cb = 0.032            !Delvigne 1988 Oil & Chem Poll
    DOUBLE PRECISION, PARAMETER :: eps = 1.0E-06            !smallest number > 0.0
    DOUBLE PRECISION, PARAMETER :: Gravity = 9.8056        !gravity (m2/s)
    DOUBLE PRECISION, PARAMETER :: KinViscWater = 1.15E-06    !kinematic viscosity of seawater (m2/s)
    DOUBLE PRECISION, PARAMETER :: Pa = 101325                !Atmospheric Pressure (Pascals)
    DOUBLE PRECISION, PARAMETER :: R = 8.314472            !Gas constant
    DOUBLE PRECISION, PARAMETER :: RhoWater = 1026.00        !density of water at 15degC & 35psu (kg/m3)
    DOUBLE PRECISION, PARAMETER :: sec2min = 0.016666         !convert seconds to minutes
    DOUBLE PRECISION, PARAMETER :: SpreadCoeff = 25.00        !CONCAWE spreading coefficient (dynes/cm)
    DOUBLE PRECISION, PARAMETER :: ViscCt = 5000.00        !Viscosity Ct proportionality constant, see ADIOS (degK)
    DOUBLE PRECISION, PARAMETER :: m32bbl = 6.2933            !convert cumecs to barrels
    DOUBLE PRECISION, PARAMETER :: ms2kts = 1.9438            !convert m/s to knots
    DOUBLE PRECISION, PARAMETER :: WaterTemp15 = 15.0        !reference water temperature (degC)

    DOUBLE PRECISION:: WindSpeed  !Zlev-OGS
    DOUBLE PRECISION:: WaterTempInst   !Zlev-OGS
    DOUBLE PRECISION :: TIMESPILLSTARTS
    LOGICAL :: SPILLSTARTED
    INTEGER :: oil_time
      CONTAINS

      !*****************************
      !*     Subroutine OilModel   *
    !*****************************
    SUBROUTINE OilModel(o_attrib,par,nParWater,nParBeaching,WindInst,TempInst,WAngleInst,WaterDepth)
        USE PARAM_MOD,    ONLY:    numpar,idt,SecSpill,spreading,emulsification,  &
               evaporation,dispersion,beaching,pi,Wind,SaltTempOn,Uwind_10,Vwind_10,  &
               WaterTemp,Ext0,VolumeSpill,iprint


        IMPLICIT NONE
        INTEGER, INTENT(IN) :: o_attrib
        INTEGER, INTENT(IN) :: nParWater ! Added by OGS: Number of UNBEACHED Particles = in circulation in water at current timestep
        INTEGER, INTENT(IN) :: nParBeaching ! Added by OGS: number of beaching particles at current timestep
        DOUBLE PRECISION, INTENT(INOUT) :: par(numpar,o_attrib)
        DOUBLE PRECISION, INTENT(IN) :: WindInst,TempInst,WAngleInst,WaterDepth

        !copied from LTRANS main program
        INTEGER, PARAMETER :: pX        =  1  ! Particle X-coordinate
        INTEGER, PARAMETER :: pY        =  2  ! Particle Y-coordinate
!        INTEGER, PARAMETER :: pZ        =  3  ! Particle Z-coordinate
!        INTEGER, PARAMETER :: pnX       =  4  ! Particle new X-coordinate
!        INTEGER, PARAMETER :: pnY       =  5  ! Particle new Y-coordinate
!        INTEGER, PARAMETER :: pnZ       =  6  ! Particle new Z-coordinate
!        INTEGER, PARAMETER :: ppX       =  7  ! Particle previous X-coordinate
!        INTEGER, PARAMETER :: ppY       =  8  ! Particle previous Y-coordinate
!        INTEGER, PARAMETER :: ppZ       =  9  ! Particle previous Z-coordinate
!        INTEGER, PARAMETER :: pStatus   = 10  ! Status of particle (previously Color)
!        INTEGER, PARAMETER :: pDOB      = 11  ! Particle Date Of Birth
!        INTEGER, PARAMETER :: pAge      = 12  ! Particle Age (s)
!        INTEGER, PARAMETER :: pLifespan = 13  ! Age at which particle settled or died

        !Local subroutine variables
        DOUBLE PRECISION :: param_a,param_b
        DOUBLE PRECISION :: x_diff,y_diff,Ud,Vd
        DOUBLE PRECISION :: ran1,ran2
        INTEGER :: SprdCase
        INTEGER :: n                            !counter for particles
        INTEGER :: ElapsedTime        !time in seconds
        DOUBLE PRECISION :: radius
        DOUBLE PRECISION :: DM_Evap,DM_Disp,DM_Beach,WeatherPercent
        AvWindAngle=WAngleInst

        oil_time = Ext0+iT * idt

        ! Check if time to spill oil
        ! Currently only modelling a single instantaneous spill
        IF(Wind)then
           WindSpeed=WindInst
        else
           WindSpeed = sqrt((Uwind_10**2.0) + (Vwind_10**2.0))
        endif
        !WRITE(*,*)'OilTime & SecSpill = ',Oil_Time,'  ',SecSpill,'Wind=',WindSpeed
        IF(.not. SaltTempOn)then
           WaterTempInst=WaterTemp
        else
           WaterTempInst =TempInst 
        endif

        IF(SPILLSTARTED.eqv..False. .and. oil_time >= TIMESPILLSTARTS) THEN
           
           ! write(*,*)'CALLING EMISSION'
            CALL Emission(Phase1Time,radius)

            DO n=1,numpar
                CALL Distribute(radius, x_diff,y_diff)    !Assign 'spilled' particle locations
                CALL update_oil_particles(0,par(n,pX),par(n,pY),x_diff,y_diff)
            END DO
            SPILLSTARTED=.TRUE.

        END IF

        ! If simulation time is greater than Phase One time then
        ! begin spreading and weathering the oil
        IF(VolumeOil>0 .and. oil_time >= TIMESPILLSTARTS+Phase1Time) THEN

            ElapsedTime = oil_time - Phase1Time                    !Time since Phase1Time reached (s)

            !Mechanical spreading of oil (Fay, Lehr etc)
            IF (Spreading) THEN

                CALL SpreadOptions(ElapsedTime, SprdCase, param_a, param_b)
               !!Added by OGS : updating SlickTickness and OilThickness
               !! after Beaching with AreaOil updated by StreadOptions
               !SlickThickness = VolumeSlick / AreaOil 
               !OilThickness = VolumeOil / AreaOil

                IF(SprdCase > 1)THEN        !using CoefR, CoefQ multiplier values

                    DO n = 1, numpar
                        CALL update_oil_particles(SprdCase,par(n,pX),par(n,pY),param_a,param_b)
                    END DO

                ELSE                        !using MOHID dispersion relationships

                    DO n = 1, numpar
                        CALL random_number(ran1)
                        CALL random_number(ran2)

                        Ud = ran1 * cos(2.0 * pi * ran2) * param_a    !param_a = DiffVelocity
                        Vd = ran1 * sin(2.0 * pi * ran2) * param_a

                        !dummy param_a and param_b as Ud,Vd (not tidy but will do!)

                        CALL update_oil_particles(SprdCase,par(n,pX),par(n,pY),Ud,Vd)
                    END DO

                END IF    ! SprdCase


            END IF    ! Spreading

            IF (Emulsification) THEN
                CALL Emulsify(ElapsedTime)
                VolumeSlick = VolumeOil / (1 - WaterContent)
                SlickThickness = VolumeSlick / AreaOil
                OilThickness = VolumeOil / AreaOil
            ELSE
                WaterContent = 0.0
            END IF    ! Emulsification

            IF (Evaporation) THEN
                CALL Evaporate(ElapsedTime, FirstAP)
            END IF    ! Evaporation

            CALL Density

            CALL Viscosity

            IF (Dispersion) THEN
                call Disperse(ElapsedTime, FirstAP,WaterDepth)


            END IF

            IF (Beaching .and. MassSpill-MassEvap-MassDisp>0) THEN
              !Removing Beached Oil  (Added by OGS)
              CALL Remove_Beached_Oil_from_Weathering_Oil(nParWater,nParBeaching)
            ENDIF
             !update finally to reflect losses due oil weathering         ! OGS:added from old oiltrans version 

               if(MassSpill-MassEvap-MassDisp-MassBeached<0)then
                  DM_Evap=MassEvap-MassEvapPrev
                  DM_Disp=MassDisp-MassDispPrev
                  DM_Beach=MassBeached-MassBeachedPrev
                  WeatherPercent=(MassSpill-MassEvapPrev-MassDispPrev-MassBeachedPrev) / (DM_Evap+DM_Disp+DM_Beach)
                  MassEvap= MassEvapPrev + WeatherPercent*DM_Evap
                  MassDisp= MassDispPrev + WeatherPercent*DM_Disp
                  MassBeached= MassBeachedPrev + WeatherPercent*DM_Beach
               endif
               MassEvapPrev=MassEvap
               MassDispPrev=MassDisp
               MassBeachedPrev=MassBeached


!             VolumeOil       =  VolumeSpill - VolumeBeached - VolumeEvap & ! OGS:added from old oiltrans version 
!                              - VolumeDisp                              ! OGS:added from old oiltrans version  
!             MassOil       = VolumeOil * RhoOil                           ! OGS:added from old oiltrans version 

             MassOil = MassSpill - MassEvap - MassDisp - MassBeached   ! changed by OGS as Mass is conserved, not volume
             VolumeOil = MassOil / RhoOil                              ! and all variables depend on mass variations  

             OilThickness   = VolumeOil / AreaOil                         ! OGS:added from old oiltrans version 
             VolumeSlick    = VolumeOil / (1 - WaterContent)              ! OGS:added from old oiltrans version   
             if(VolumeOil < 0.0)VolumeOil = 0.0                           ! OGS:added from old oiltrans version 
             FirstAP         = .false.                                    ! OGS:added from old oiltrans version 

        END IF
        if(mod(oil_time,iprint)==0) call WriteOilOutputs() ! Added by OGS for LTRANS-Zlev

        iT = iT + 1                                             !increment timestep counter
    END SUBROUTINE OilModel

      !*****************************************
      !*     Subroutine update_oil_particles   *
    !*****************************************
    SUBROUTINE update_oil_particles(n, loc_x, loc_y, param_1, param_2)
    USE PARAM_MOD, ONLY: idt

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION, INTENT(INOUT):: loc_x, loc_y
    DOUBLE PRECISION, INTENT(IN):: param_1, param_2

    SELECT CASE (n)

        !Initial distribution after release
        CASE (0)
            loc_x  = loc_x + param_1
            loc_y  = loc_y + param_2

        !MOHID Diffusion
        CASE(1)
            loc_x = loc_x + (param_1 * idt)
            loc_y = loc_y + (param_2 * idt)

        !ADIOS/CONCAW circular spreading radius (only using CoefR: param_1)
        CASE(2:3)
            loc_x = loc_x + (loc_x * param_1)
            loc_y = loc_y + (loc_y * param_1)

        !OILPOL ellipsical spreading radius
        CASE(4)
            loc_x = loc_x + (loc_x * param_1)    !CoefR
            loc_y = loc_y + (loc_y * param_2)    !coefQ

        CASE DEFAULT
            loc_x = Loc_x    !don't move particle
            loc_y = loc_y    !don't move particle

    END SELECT

    END SUBROUTINE update_oil_particles

    !*****************************
    !*     Subroutine Density    *
    !*****************************
    SUBROUTINE Density
    !Buchanan, I. 1988
    !"Methods for predicting the physical changes of oil spilled at sea"
    !Oil and Chemical pollution vol 4(4) pp311-328
        USE PARAM_MOD,    ONLY:     Evaporation,  &
                                Emulsification
        IMPLICIT NONE

        !First re-calculate RhoOil based on temperature correction.
        RhoOil = OilDensity *  (1.0 - CDensT * ((WaterTempInst + 273.15) - OilDensity_RefT))            !actual density of spilled oil at ocean temperature (kg/m3)

        !Check if Evaporation enabled - if so include Evap effect on RhoOil already calculated
        IF (Evaporation) THEN
            RhoOil = RhoOil * (1.0 + CDensE * (MassEvap / MassSpill))
        END IF

        !Check if Emulsification enabled - if so include Emuls effect on RhoOil already calculated
        IF (Emulsification) THEN
            RhoOil = (WaterContent * RhoWater) + (RhoOil * (1.0 - WaterContent))
        END IF

        DeltaRho  = (RhoWater - RhoOil) / RhoWater

        RETURN

    END SUBROUTINE Density

    !*******************************
    !*     Subroutine Distribute   *
    !*******************************
    SUBROUTINE Distribute(lenR, XDiff, YDiff)
    !Initially randomly distribute particles using a normal distribution throughout the theoretical Fay area of spill
        USE PARAM_MOD,    ONLY: PI
        IMPLICIT NONE

        !I/O variables
        DOUBLE PRECISION, INTENT(IN):: lenR
        DOUBLE PRECISION, INTENT(OUT):: XDiff
        DOUBLE PRECISION, INTENT(OUT):: YDiff

        !Local variables
        DOUBLE PRECISION:: ran        !holder for random number
        DOUBLE PRECISION:: ang        !randomly generated angle (0 < ang < 2pi)
        DOUBLE PRECISION:: dist        !randomly generated distance (0 < dist < Radius)

        CALL random_number(ran)
        ang = ran * pi * 2.0

        CALL random_number(ran)
        dist = ran * lenR

        XDiff = dist * cos(ang)
        YDiff = dist * sin(ang)

    END SUBROUTINE Distribute


    !******************************
    !*     Subroutine Emission    *
    !******************************
    SUBROUTINE Emission(Phase1Time,lenR)
        USE PARAM_MOD, ONLY: VolumeSpill,PI,SecSpill,Oil_Resin,Oil_Asph

        IMPLICIT NONE

        !I/O variables
        INTEGER, INTENT(OUT):: Phase1Time
        DOUBLE PRECISION, INTENT(OUT):: lenR

        !Local variables
!        DOUBLE PRECISION:: XDiff
!        DOUBLE PRECISION:: YDiff
        WRITE(*,*)'IN EMISSION'
        WRITE(*,*)'VOLUMESPILL = ',VOLUMESPILL
        WRITE(*,*)'RHOOIL = ',RHOOIL
        MassSpill   = VolumeSpill * RhoOil                    !(kg)

        CALL InitialArea(Phase1Time)                        !Calculate initial area at end of gravity-inertia phase (CONCAWE)

        Phase1Time   = (SecSpill * 3600) + Phase1Time        !(sec) Re-calculate Phase1Time in terms of simulation time (sec)
        ThickLimit   = F_APIThickness()                        !(m) Thickness limit to end viscous spreading
        lenR         = sqrt(AreaOil / pi)                    !(m) Initial Radius of spilled oil at end of inertial spreading
        OilThickness = VolumeSpill / AreaOil                !(m) Initial thickness of oil at end of inertial spreading
        SlickThickness = OilThickness                        !(m) Initial thickness of oil slick at end of inertial spreading
        VolumeOil      = VolumeSpill                            !(m3)
        MassOil         = MassSpill                            !(kg)
       ! write(*,*) 'Phase1Time=',Phase1Time
        ResinOil = Oil_Resin
        AsphOil = Oil_Asph 
   END SUBROUTINE Emission

    !******************************
    !*     Subroutine InitialArea *
    !******************************
    SUBROUTINE InitialArea(Phase1Time)
    !Calculates the time to the end of the gravity-inertial spreading phase (which we don't model)
    !and calculates the area of the spreaded oil at the end of that phase.
        USE PARAM_MOD,    ONLY: AreaOption, pi, VolumeSpill

        IMPLICIT NONE

        !I/O variables
        INTEGER, INTENT(OUT):: Phase1Time        !Computed time to gravity-inertial spreading phase (sec)

        !Local variables
        DOUBLE PRECISION:: Area1
        DOUBLE PRECISION:: Area2
        DOUBLE PRECISION:: Time

        !When oil is denser than surrounding water, oil will sink => no spreading
        IF (RhoOil > RhoWater) THEN
            WRITE(*,*)'Oil is denser than water => no surface spreading'
            !INCLUDE CODE FOR VERTICAL MOVEMENT OF OIL PARTICLES
        END IF

        WRITE(*,*)'IN INITIAL AREA'
        WRITE(*,*)'AREA_OPTION = ',F_UpCase(trim(adjustl(AreaOption)))
        SELECT CASE (F_UpCase(trim(adjustl(AreaOption))))

        CASE ("MOHID2")
            !MOHID Description, 2003
            !downloaded from http://maretec.mohid.com/PublicData/Products/Manuals/Mohid_Description.pdf
            AreaOil = pi * (1.45**4.0 / 1.14**2.0) * 0.25 * (VolumeSpill**5.0 * Gravity * DeltaRho / (KinViscWater**2.0))**(1.0/6.0)
            Phase1Time = nint((0.725/0.570)**4.0 *(VolumeSpill / (KinViscWater * Gravity * DeltaRho))**(1.0/3.0))

        CASE ("ADIOS2")
            !NOAA, 2000
            !"ADIOS (Automated Data Inquiry for Oil Spills) version 2.0.1 online help manual"
            !Hazardous Materials Response and Assessment Division,NOAA.
            !Prepared for the U.S. Coast Guard Research and Development Center, Groton Connecticut
            AreaOil = pi * (1.21**4.0 / 1.53**2.0) * (VolumeSpill**5.0 * Gravity * DeltaRho / (KinViscWater**2.0))**(1.0/6.0)
            Phase1Time = nint((1.21/1.53)**4.0 * (VolumeSpill / (KinViscWater * Gravity * DeltaRho))**(1.0/3.0))

        CASE ("CONCAW")
            !van Oudenhoven, J., Draper, V., Ebbon, G., Holmes, P., Nooyen, J. 1983
            !"Characteristics of petroleum and its behavior at sea"
            !CONCAWE Report No.8/83. Den Haag, November 1983.
            Time = 0.10            !seconds (ie: almost instantaneous)
            Area1 = 0.0            !gravity-inertia
            Area2 = 1.0            !gravity-viscous
            DO WHILE (Area1 < Area2)    !until gravity-viscous regime reached
                Area1 = pi * (1.14**2.0) * (DeltaRho * Gravity * VolumeSpill)**0.5 * Time
                Area2 = pi * (0.98**2.0) * ((DeltaRho * Gravity * (VolumeSpill**2.0))/(KinViscWater**0.5))**(1.0/3.0) * (Time**0.5)
                Time = Time + 10.0
            END DO
            AreaOil = Area2        !Area at start of gravity-viscous regime
            Phase1Time = nint(Time)    !Time since spill to start of gravity-viscous regime

        CASE ("OILPOL")    !GULFSPILL
            !Rabeh, A.H., Lardner, R.W., Gunay, N. 2000
            !"GulfSpill Version 2.0: a software package for oil spills in the Arabian Gulf"
            !Environmental Modelling and Software 15 (2000) 425-442
            AreaOil = pi * (2.81 * sqrt(VolumeSpill))**2.0
            Phase1Time = nint(0.00)

        CASE DEFAULT
            WRITE(*,*)'       Case not encoded'
            WRITE(*,*)'No INITIAL SPILL AREA calculated'
            WRITE(*,*)'****** PROGRAM TERMINATING ******'
            STOP    !Terminate program as no value for oil area computed to run further weathering processes with

        END SELECT

    END SUBROUTINE InitialArea

    !**************************************
    !*     Subroutine InitOil             *
    !**************************************
    ! OGS for LTRANS-Zlev: commented the next 2 lines and replaced them by
    ! following 2 lines
!    SUBROUTINE InitOilModel(RhoOil,OilDensity,OilDensity_RefT,DeltaRho,     &
!                                ViscOil)
    SUBROUTINE InitOilModel(TempInst,pTS)

        USE PARAM_MOD,     ONLY:     Cut_Unit,Cut_Temp,Oil_Dens,       &
                                Oil_Dens_RefT,API,Dyn_Visc,Dyn_Visc_RefT,   &
                                Kin_Visc,Kin_Visc_RefT,Oil_Asph,Oil_Resin,  &
                                Oil_Sat,WaterTemp,SaltTempOn

        IMPLICIT NONE
        !I/O variables
        DOUBLE PRECISION, INTENT(IN):: TempInst    ! Added by OGS for LTRANS-Zlev
        DOUBLE PRECISION, INTENT(IN) :: pTS        ! time of spill starts as read in iniparloc (Added by OGS)
    ! OGS for LTRANS-Zlev: commented the next 2 lines
!        DOUBLE PRECISION, INTENT(OUT):: RhoOil,OilDensity,OilDensity_RefT,  &
!                                        DeltaRho,ViscOil

        !Local variables
        DOUBLE PRECISION:: SGOil15                            !Specific gravity of oil at 15degC
        DOUBLE PRECISION:: RhoOil15                            !Density of oil at 15degC (kg/m3)
        INTEGER:: i, n, NumCuts                                !counters
        LOGICAL:: novalue                                    !logic controllers
        DOUBLE PRECISION:: LowPresTemps(5)                    !array for holding converted 40mmHG pressure temperatures
        DOUBLE PRECISION, PARAMETER :: RhoFWater15 = 1000.0 !density of water at 15degC & 0psu (kg/m3)
        TIMESPILLSTARTS=pTS
        SPILLSTARTED=.FALSE.
        FirstAP         = .true.
        WRITE(*,*)''
        WRITE(*,*)'IN INITOILMODEL'
        !***************************
        !*   INTIALISE VARIABLES   *
        !***************************
        ! Equates to     200.0, 225.0, 250.0, 275.0, 300.0 degC at 40mmHg
        LowPresTemps = (/307.8, 337.7, 365.8, 394.9, 424.0/) !(K)
        ! see BPO Crude Oil Analysis Data Bank User's Guide Methods
        ! Taken from standard pressure-temperature nomograph

         MassOil = 0.0
         MassSpill = 0.0
         MassEvap = 0.0
         MassDisp = 0.0
         VolumeOil = 0.0
         VolumeSlick = 0.0
         VolumeEvap = 0.0
         VolumeDisp = 0.0
         AreaOil = 0.0
         ViscOil = 0.0
         RhoOil = 0.0
         OilThickness = 0.0
         SlickThickness = 0.0
         DeltaRho = 0.0
         ThickLimit = 0.0
         lenR = 0.0
         WaterContent = 0.0
         xviscemul = 0.0
         OilDensity = 0.0
         OilDensity_RefT = 0.0
        VolumeBeached = 0.0 ! Added by OGS

        novalue = .FALSE.
        FIRSTAP = .TRUE.                                    ! first time for oil model processes
        iT = 0                                                ! zero internal timestep counter

        !*******************************
        !*    INTIALISE OIL PROPERTIES   *
        !*******************************

       !Determine no. of cuts for volume distillation only
       !All other options use the ADIOS correlation eqn for pseudocomponent evaporation
        SELECT CASE (trim(adjustl(Cut_Unit)))
            CASE ("-----")
                NumCuts = 0

            CASE ("weight")
                NumCuts = 0

            CASE ("volume")
                n = 1
                NumCuts = 0
                !determine number of cuts
                   DO WHILE (.NOT.novalue .AND. n < 16)    !n = 15 [max number of cuts (10@760mmHg & 5@40mmHg)]
                    IF (Cut_Temp(n) == 0.0) THEN
                        novalue = .TRUE.
                    ELSE
                        NumCuts = n
                        n = n + 1
                    END IF
                END DO
                !determine if cuts have been made at reduced pressure
                i = 1
                DO n = 1, NumCuts-1
                    IF (Cut_Temp(n) > Cut_Temp(n+1)) THEN            !if temperature steps down => reduced pressure distillation
                        Cut_Temp(n+1) = LowPresTemps(i)    + 273.15    !assign temp correction based on P-T nomograph
                        i = i + 1
                    END IF
                END DO

            CASE DEFAULT
                WRITE(*,*)'Error: Cut_Unit not defined'
                WRITE(*,*)'*** FATAL ERROR - STOP ***'
                STOP

        END SELECT

        IF(.not.SaltTempOn)then
           WaterTempInst=WaterTemp
        else
           WaterTempInst =TempInst 
        endif

        !Determine initial oil density
        IF (Oil_Dens > 0.0) THEN
            RhoOil      = Oil_Dens * (1.0 - CDensT * ((WaterTempInst+273.15) - Oil_Dens_RefT))    !actual density of spilled oil at ocean temperature (kg/m3)
            OilDensity = Oil_Dens
            OilDensity_RefT = Oil_Dens_RefT
        ELSEIF (API > 0.0) THEN
            SGOil15  = 141.5 / (131.5 + API)                                    !specific gravity of spilled oil at reference temp (15.5degC) - from API standards
            RhoOil15 = SGOil15 * RhoFWater15                                    !density of spilled oil at reference temp (15.5degC)
            RhoOil   = RhoOil15 * (1.0 - CDensT * (WaterTempInst - WaterTemp15))     !actual density of spilled oil at ocean temperature (kg/m3)
            OilDensity = RhoOil15                            !Reference oil density for subroutine DENSITY
            OilDensity_RefT = WaterTemp15 + 273.15        !reference oil temperature for subroutine DENSITY
        ELSE
            WRITE(*,*)'No value DENSITY associated with this oil'
            WRITE(*,*)'No modelling can be done'
            WRITE(*,*)'*** FATAL ERROR - STOP ***'
            STOP
        END IF

        DeltaRho  = (RhoWater - RhoOil) / RhoWater        !relative density difference between water and oil densities

        !Determine initial oil dynamic viscosity
        IF (Dyn_Visc > 0.0) THEN
            ViscOil = 1000.0 * Dyn_Visc * exp(ViscCt * ((1.0/(WaterTempInst+273.15)) - (1.0/Dyn_Visc_RefT)))!Dynamic (kg/ms -> cP)
        ELSE
            ViscOil = Kin_Visc * exp(ViscCt * ((1.0/(WaterTempInst+273.15)) - (1.0/Kin_Visc_RefT)))            !Kinematic (m2/s -> cSt)
            ViscOil = 1000.0 * ViscOil * RhoOil        !Dynamic (cP)
        END IF

    call InitOilOutputs() ! Added by OGS for LTRANS-Zlev

        RETURN

    END SUBROUTINE InitOilModel

    !***********************************
    !*     Subroutine SpreadOptions    *
    !***********************************
    SUBROUTINE SpreadOptions(ElapsedTime,SprdCase,CoefR,CoefQ)
    USE PARAM_MOD,     ONLY: PI,SprdOption,idt,numpar,Uwind_10,Vwind_10

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ElapsedTime
!    DOUBLE PRECISION, INTENT(IN) :: OilThickness  ! Commented by OGS
    INTEGER, INTENT(OUT) :: SprdCase
    DOUBLE PRECISION, INTENT(OUT) :: CoefR, CoefQ

    DOUBLE PRECISION:: DiffCoef        !MOHID
    DOUBLE PRECISION:: DiffVelocity !MOHID
    DOUBLE PRECISION:: Ud            !MOHID
    DOUBLE PRECISION:: Vd            !MOHID
    DOUBLE PRECISION:: beta            !MOHID
    DOUBLE PRECISION:: beta_old        !MOHID
    DOUBLE PRECISION:: aux1            !MOHID
    DOUBLE PRECISION:: aux2            !MOHID
    DOUBLE PRECISION:: ran1            !MOHID
    DOUBLE PRECISION:: ran2            !MOHID
    DOUBLE PRECISION:: Area            !CONCAWE & ADIOS
    DOUBLE PRECISION:: Area2        !CONCAWE & ADIOS
    DOUBLE PRECISION:: Area3        !CONCAWE & ADIOS
    DOUBLE PRECISION:: MaxArea        !CONCAWE & ADIOS
    DOUBLE PRECISION:: DeltaR        !CONCAWE & ADIOS
!    DOUBLE PRECISION:: CoefR        !CONCAWE & ADIOS
    DOUBLE PRECISION:: DeltaQ        !OILPOL
!    DOUBLE PRECISION:: CoefQ        !OILPOL
    DOUBLE PRECISION:: lenQ            !OILPOL
    INTEGER:: n

    !save values on exit
    SAVE:: beta_old

!    IF (OilThickness > ThickLimit) THEN

        SELECT CASE (F_UpCase(SprdOption))

        !FAY based method
        CASE ("MOHID2")
            !MOHID Description, 2003
            !downloaded from http://maretec.mohid.com/PublicData/Products/Manuals/Mohid_Description.pdf
            !& with reference to source code.
            SprdCase = 1
               beta = pi * ((1.45**2)/4.0) * (DeltaRho * Gravity * (VolumeOil**2) &
                        / sqrt(KinViscWater))**(1.0/3.0)

               IF (FirstAP) THEN
                   beta_old = beta
               END IF

               aux1       = beta_old
              aux2       = AreaOil/aux1
            beta_old = beta
            AreaOil  = beta * sqrt((aux2)**2 + IDT)
            OilThickness = VolumeOil / AreaOil                !(m) Updated average thickness of oil

            DiffCoef = (pi * (0.725**2.0))/16.0 &
                        * ((DeltaRho * Gravity * (VolumeOil**2.0)) / sqrt(KinViscWater))**(1.0/3.0) &
                        * (1.0/ sqrt (REAL(ElapsedTime,kind(1))))

            DiffVelocity = sqrt ((2.0 * DiffCoef) / IDT)

            CoefR = DiffVelocity
            CoefQ = 0.0

        !Fay based method
        CASE ("ADIOS2")
            !NOAA, 2000
            !"ADIOS (Automated Data Inquiry for Oil Spills) version 2.0.1 online help manual"
            !Hazardous Materials Response and Assessment Division,NOAA.
            SprdCase = 2
            Area          = pi * (1.21**2.0) * ((VolumeOil**2.0) * Gravity * DeltaRho * (ElapsedTime**(3.0/2.0)) &
                                / sqrt(KinViscWater))**(1.0/3.0)
            AreaOil      = Area                                            !Updated area of oil
            OilThickness = VolumeOil / AreaOil                            !Updated average thickness of oil

            lenR          = sqrt(Area/pi)
            DeltaR          = lenR - sqrt(AreaOil/pi)                        !NewRadius - OldRadius

            CoefR          = DeltaR / lenR                                !DeltaR / NewRadius
            CoefQ = 0.0

        !FAY based method
        CASE ("CONCAW")
            !van Oudenhoven, J., Draper, V., Ebbon, G., Holmes, P., Nooyen, J. 1983
            !"Characteristics of petroleum and its behavior at sea"
            !CONCAWE Report No.8/83. Den Haag, November 1983.
            SprdCase = 3
            Area2 = pi * (0.98**2.0) * ((DeltaRho * Gravity * (VolumeOil**2.0))/(KinViscWater**0.5))**(1.0/3.0) * (ElapsedTime**0.5)
            Area3 = pi * (1.60**2.0) * (((SpreadCoeff * 10**(-3.0))**2.0) &
                /((RhoWater**2.0) * KinViscWater))**0.5 * (ElapsedTime**1.5)
            MaxArea = 10**5.0 * (VolumeOil**0.75)
            Area = Area2
            IF (Area2 < Area3) THEN
                Area = Area3
                IF (Area3 > MaxArea) THEN
                    Area = MaxArea
                END IF
            END IF
            AreaOil      = Area                                            !Updated area of oil
            OilThickness = VolumeOil / AreaOil                            !Updated average thickness of oil

            lenR          = sqrt(Area/pi)
            DeltaR          = lenR - sqrt(AreaOil/pi)                        !NewRadius - OldRadius

            CoefR          = DeltaR / lenR                                !DeltaR / NewRadius
            CoefQ = 0.0

        !Lehr based method
        CASE ("OILPOL")     !GULFSPILL
            !OILPOL_2 solution - see OILPOL2.xls
            !Chao, X., Shankar, J., Wang, S. 2003
            !"Development and application of oil spill model for Singapore coastal waters"
            !Journal of Hydraulic Engineering 129:7 (2003) 495-503
            SprdCase = 4
            !WindSpeed = sqrt((Uwind_10**2.0) + (Vwind_10**2.0))
            AreaOil = 2270.0 * (DeltaRho * (RhoWater/RhoOil))**(2.0/3.0) &
                * (VolumeOil * m32bbl)**(2.0/3.0) * (ElapsedTime * sec2min)**(1.0/2.0) &
                + (40.0 * (DeltaRho * (RhoWater/RhoOil))**(1.0/3.0) * (ElapsedTime &
                 * sec2min) * (WindSpeed * ms2kts)**(4.0/3.0))
            OilThickness   = VolumeOil / AreaOil                        !Updated average thickness of oil

            lenQ = 53.76 * (DeltaRho * (RhoWater/RhoOil))**(1.0/3.0) *  &
                (VolumeOil * m32bbl)**(1.0/3.0) * (ElapsedTime * sec2min)**(1.0/4.0)
            lenR = lenQ + 0.95 * (WindSpeed * ms2kts)**(4.0/3.0) *  &
                (ElapsedTime * sec2min)**(3.0/4.0)

            DeltaQ = 53.76 * (DeltaRho * (RhoWater/RhoOil))**(1.0/3.0) &
                 * (VolumeOil * m32bbl)**(1.0/3.0) * 0.25 * (ElapsedTime * sec2min)**(-3.0/4.0) * IDT
            DeltaR = DeltaQ + 0.95 * (WindSpeed * ms2kts)**(4.0/3.0) * &
                 0.75 * (ElapsedTime * sec2min)**(-1.0/4.0) * IDT

            CoefR = DeltaR / lenR
            CoefQ = DeltaQ / lenQ

        CASE DEFAULT
            WRITE(*,*)'Case not encoded'
            WRITE(*,*)'No SPREADING processes modelled'
        END SELECT

!    ELSE
!        AreaOil = AreaOil
!        WRITE(*,*)'******* Thickness limit reached *******'
!        write(*,*)'No further viscous spreading of oil slick'
!    END IF

    RETURN

    END SUBROUTINE SpreadOptions

    !*******************************
    !*     Subroutine Viscosity    *
    !*******************************
    SUBROUTINE Viscosity
        USE PARAM_MOD, ONLY: Dyn_Visc, Kin_Visc, Dyn_Visc_RefT, Kin_Visc_RefT, &
                            Evaporation,  &
                            Emulsification
        IMPLICIT NONE
        !NewViscosity     = RefViscosity * exp(dViscTemp + dViscEvap) * xViscEmul
        !                = RefViscosity * exp(dViscTemp) * exp(dViscEvap) * xViscEmul
       ! write(*,*)'Visc t=',ViscOil
        IF (Dyn_Visc > 0.0) THEN
            ViscOil = 1000 * Dyn_Visc * exp(ViscCt * ((1.0/(WaterTempInst+273.15)) - (1.0/Dyn_Visc_RefT)))    !Dynamic (kg/ms -> cP)
         ! write(*,*)'Visc using Dyn_Visc',ViscOil,Dyn_Visc,ViscCt,WaterTempInst,Dyn_Visc_RefT
        ELSE
            ViscOil = 1000 * Kin_Visc * exp(ViscCt * ((1.0/(WaterTempInst+273.15)) - (1.0/Kin_Visc_RefT)))    !Kinematic (m2/s -> cSt)
            ViscOil = ViscOil * RhoOil !Dynamic (cP)
         ! write(*,*)'Visc using Kin_Visc',ViscOil,Kin_Visc,ViscCt,WaterTempInst,Kin_Visc_RefT,RhoOil
          END IF

        IF (Evaporation) THEN !(Mackay,1980)
           ViscOil = ViscOil * exp((Evap_C4 * (MassEvap / MassSpill)))
         ! write(*,*)'Visc post Evap',ViscOil,Evap_C4,MassEvap,MassSpill    
         END IF

        IF (Emulsification) THEN !(Fingas 2011)
           ViscOil = ViscOil * xViscEmul
         ! write(*,*)'Visc post Emul',ViscOil,xViscEmul

         END IF
       ! write(*,*)'Visc t+1 =',ViscOil


        RETURN

    END SUBROUTINE Viscosity

    !*********************************
    !*     Function APIThickness   *
    !*********************************
    DOUBLE PRECISION FUNCTION F_APIThickness() RESULT(APIThickness)
        USE PARAM_MOD, ONLY: API

        IMPLICIT NONE

        IF (API < 17.5 ) THEN
            APIThickness = 1.0E-04    !(m)
        ELSEIF ((API > 17.5 ) .AND. (API < 45.0 )) THEN
            APIThickness = 1.5727E-04 - 3.2727E-06 * API
        ELSEIF (API > 45.0 ) THEN
            APIThickness = 1.0E-05    !(m)
        END IF

    END FUNCTION F_APIThickness

    !*********************************
    !*     Function WindAngle      *
    !*********************************
    DOUBLE PRECISION FUNCTION F_WindAngle(wangle) RESULT(WindAngle)
    ! Function to give the Cartesian direction in radians of a wind or wave direction specified as 'direction coming from in degrees'
        USE PARAM_MOD, ONLY: PI

        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: wangle
        DOUBLE PRECISION :: Ang,Theta

        Ang = wangle*PI/180.0

        write(*,*)"WARNING until formulations used by F_WindAngle "
        write(*,*)"won't be verified, this subroutine must not be used."
        write(*,*)"program stops"
        stop
        IF (Ang.LT.PI/2.0) THEN          ! Check formulation  
           Theta = -1.0*PI/2.0 - Ang     ! Check formulation 
        ELSE                             ! Check formulation 
           Theta = 3.0*PI/2.0 - Ang      ! Check formulation 
        END IF                           ! Check formulation 

        WindAngle = Theta

    END FUNCTION F_WindAngle

    !****************************
    !*     Function F_UpCase    *
    !****************************
    FUNCTION F_UpCase(string) RESULT(upper)
    IMPLICIT NONE

    !I/O variables
    CHARACTER(LEN=*), INTENT(IN) :: string
    CHARACTER(LEN=len(string)):: upper

    !Local variables
    INTEGER:: j

    DO j = 1,len(string)
        IF(string(j:j) >= "a" .AND. string(j:j) <= "z") THEN
            upper(j:j) = achar(iachar(string(j:j)) - 32)
        ELSE
            upper(j:j) = string(j:j)
        END IF
    END DO

    END FUNCTION F_UpCase

    !*****************************
    !*     Subroutine Emulsify   *
    !*****************************
    SUBROUTINE Emulsify(ElapsedTime)
    !Fingas,M.,2011, "Models for Water-in-Oil Emulsion Formation" in
    !Chpt. 10 of Oil Spill Science and Technology, 2011. Gulf Professional Publishing, UK. ISBN:978-1-85617-943-0
    USE PARAM_MOD,    ONLY: Oil_Sat,idt,SigWaveHeight,windwavemodel

    IMPLICIT NONE

    !I/O variables
    INTEGER, INTENT(IN):: ElapsedTime

    !Local variables
    DOUBLE PRECISION:: RelRhoOil            !Relative density (decimal)
    DOUBLE PRECISION:: D_t                    !transformed density
    DOUBLE PRECISION:: V_t                    !transformed viscosity
    DOUBLE PRECISION:: S_t                    !transformed saturates
    DOUBLE PRECISION:: R_t                    !transformed resins
    DOUBLE PRECISION:: A_t                    !transformed asphaltenes
    DOUBLE PRECISION:: AR                    !asphaltene/resin ratio
    DOUBLE PRECISION:: AR_t                    !transformed asphaltene/resin ratio
    DOUBLE PRECISION:: StabilityC            !StabilityC index
    DOUBLE PRECISION:: Fingas_Class(4,8)            !Array to hold Fingas empricial values
    DOUBLE PRECISION:: FormTime                !Time for emulsion state to form (min)
    DOUBLE PRECISION:: StartTime            !Start time from which to begin timing emulsion formation, etc (sec)
    DOUBLE PRECISION:: FingasDay            !Seconds from StartTime to end of one day (including formation time)
    DOUBLE PRECISION:: FingasWeek            !as above to end of one week
    DOUBLE PRECISION:: FingasYear            !as above to end of one year
    DOUBLE PRECISION:: waveheight
    INTEGER:: ClassIndex        !Class Index (1 - 4)
    INTEGER:: StartClass        !class index and integer
    INTEGER::n        !class index and integer
    INTEGER :: Ctemp
    !save values on exit
    SAVE:: ClassIndex, FormTime, StartClass, StartTime, FingasDay, FingasWeek, FingasYear,Fingas_Class

    IF (FirstAP) THEN
    !Fingas,M.,2011, "Models for Water-in-Oil Emulsion Formation" in
    !Chpt. 10 of Oil Spill Science and Technology, 2011. Gulf Professional Publishing, UK. ISBN:978-1-85617-943-0
        ! write(*,*)'FirstAp',FirstAP,' initialize emulsification'
                        !  | U     E      M      S   |
         Fingas_Class = reshape((/ 0.06,  0.42,  0.64,  0.76, &     ! DayWaterContent <- Figure 10.4
                            0.06,  0.37,  0.32,  0.76, &     ! WeekWaterContent <- Figure 10.4
                            0.06,  0.37,  0.20,  0.68, &     ! YearWaterContent <- Figure 10.4
                            0.99,   1.9,   7.2, 405.0, &     ! DayViscosityIncrease <-Table 10.4
                             1.0,   1.9,  11.0,1054.0, &     ! WeekViscosityIncrease <-Table 10.4
                             1.0,   2.1,  32.0, 991.0, &     ! YearViscosityIncrease <-Table 10.4
                             0.0,  30.8,  47.0,  27.1, &     ! FormTimeParamA <- Table 10.5
                             0.0, 18300., 49100.,  7520. /), &    ! FormTimeParamB <- Table 10.5
                             shape(Fingas_Class))

        !do Ctemp=1,4
        !   write(*,*)'Fingas_Class(',Ctemp,',:)=',Fingas_Class(Ctemp,:)
        !enddo
        WaterContent = 0.06
        xViscEmul = 1.0
        VolumeSlick = 0.0
        ClassIndex = 0
        StartClass = 0
        StartTime = 0.0
    END IF

    IF(ClassIndex < 4 ) THEN        !check until stable emulsion forms
       ! write(*,*)'EMUL ClassIndex < 4'
        !Fingas,M.,2011, "Models for Water-in-Oil Emulsion Formation" in
        !Chpt. 10 of Oil Spill Science and Technology, 2011. Gulf Professional Publishing, UK. ISBN:978-1-85617-943-0

        !Density transform (ranges checked)
        RelRhoOil = RhoOil / 1000.0
        IF (exp(RelRhoOil) < 2.5) THEN
            D_t = max(0.01d0, (2.5d0 - exp(RelRhoOil)))
        ELSE
            D_t = max(0.01d0, (exp(RelRhoOil) - 2.5d0))
        END IF

        !Viscosity transform (ranges checked)
        IF (log(ViscOil) < 5.8) THEN
            V_t = 5.8 - log(ViscOil)
        ELSE
            V_t = log(ViscOil) - 5.8    !changed from 8.7 in text
        END IF

        !Saturate transform (ranges checked)
        IF (Oil_Sat < 45.0) THEN
            S_t = 45.0 - Oil_Sat
        ELSE
            S_t = max(eps, Oil_Sat - 45.0)
        END IF

        !Resin transform (ranges checked)
        IF (ResinOil < 10.0) THEN
            R_t = max(0.1d0, 10.0 - ResinOil)
        ELSE
            R_t = max(0.1d0, ResinOil - 10.0)
        END IF

        !Asphaltene transform (ranges checked)
        IF (AsphOil < 4.0) THEN
            A_t = 4.0 - AsphOil
        ELSE
            A_t = max(eps, AsphOil - 4.0)
        END IF

        !A/R ratio transform (ranges checked)
        AR = AsphOil / ResinOil
        IF (AR < 0.6) THEN
            AR_t = max(0.025d0, 0.6 - AR)
        ELSE
            AR_t = max(0.025d0, AR - 0.6)
        END IF

        !calculate stabilityC (Eqn 15)
        StabilityC = 12.3 + (0.259 * S_t) - (1.601 * R_t) - (17.2 * AR_t) &
                    - (0.5 * (V_t**3.0)) + (0.002 * (R_t**3.0)) + (0.001 * (A_t**3.0)) + (8.51 * (AR_t**3.0)) &
                    - (1.12 * log(V_t)) + (0.7 * log(R_t)) + (2.97 * log(AR_t)) &
                    + (6.0E-08 * (exp(V_t)**2.0)) - (1.96 * (exp(AR_t)**2.0)) &
                    - (4.0E-06 * (log10(D_t)/(D_t**2.0))) - (1.5E-04 * (log10(AR_t)/(AR_t**2.0)))

        !determine emulsion state based on stabilityC (Table 10.3)
        IF(StabilityC >= 2.2 .AND. StabilityC <= 15) THEN
            ClassIndex = 4        !Stable
        ELSEIF (StabilityC >= -12.0 .AND. StabilityC <= -0.7) THEN
            ClassIndex = 3        !Mesostable
        ELSEIF (StabilityC >= -18.3 .AND. StabilityC <= -9.1) THEN
            IF(RelRhoOil > 0.96 .AND. ViscOil > 6000.0) THEN
                ClassIndex = 2    !Entrained
            END IF
        ELSEIF (StabilityC >= -39.1 .AND. StabilityC <= -7.1) THEN
            IF ( (RelRhoOil < 0.85 .OR. RelRhoOil > 1.0) .AND.  &
                (ViscOil < 100 .OR. ViscOil > 800000) .AND.  &
                (AsphOil < 1.0 .OR. ResinOil < 1.0) )THEN
                ClassIndex = 1    !Unstable
            END IF
        ELSE
!           ! write(*,*)'Unresolved Emulsion Class - assume UNSTABLE'
            ClassIndex = 1
        END IF

        !Only update emulsion state if progressively more stable emulsions form
        IF (ClassIndex > StartClass) THEN

            StartClass = ClassIndex
            StartTime = ElapsedTime

            !determine formation time for emulsion state
            IF (ClassIndex >= 2)THEN
                                    if(windwavemodel)then
                                       waveheight = SigWaveHeight                                                                                                !SWAN model output
                                    else
                                       waveheight = 0.243 * (0.71*WindSpeed**1.23)**2.0 / Gravity                                                        !ADIOS formulation (NOAA 1994)
                                    end if         
                                    FormTime = Fingas_Class(ClassIndex, 7) + &
                                      (Fingas_Class(ClassIndex, 8) / ((waveheight*100.0)**1.5)) !(minutes)
                                   ELSE
                                       FormTime = 0.0
                                   END IF
                       
            !calculate times to end of day, week and year for Fingas values (from Table 10.4)
            FingasDay = StartTime + (FormTime * 60.0) + (1.0 * 86400.0)
            FingasWeek = StartTime + (FormTime * 60.0) + (7.0 * 86400.0)
            FingasYear = StartTime + (FormTime * 60.0) + (365.0 * 86400.0)

        END IF

    END IF

    !Calculate incremental increase in viscosity and water content per DT based on times above.
    IF (StartClass == 0) THEN
       ! write(*,*)'EMUL ClassIndex ==0'
        WaterContent = WaterContent
        xViscEmul = 1.0
    ELSE
       ! write(*,*)'StartTime=',StartTime,' ElapsedTime=',ElapsedTime,' FingasDay=', FingasDay,'FingasWeek=',FingasWeek
        !do Ctemp=1,4
        !  ! write(*,*)'Fingas_Class(',Ctemp,',:)=',Fingas_Class(Ctemp,:)
        !enddo
        IF (ElapsedTime <= FingasDay) THEN !emulsion time <= FormTime + 1 Day
           ! write(*,*)'CP watCont using ',StartClass,WaterContent, &
            ! Fingas_Class(StartClass, 1),IDT,(FingasDay - StartTime)
            WaterContent = min(WaterContent + (( Fingas_Class(StartClass, 1)  &
                / (FingasDay - StartTime) ) * IDT), Fingas_Class(StartClass,1))
           ! write(*,*)'CP xViscEmul using ',StartClass,xViscEmul, &
           !  Fingas_Class(StartClass, 4),IDT,(FingasDay - StartTime)
            xViscEmul = min(xViscEmul + (( (Fingas_Class(StartClass, 4) - 1.0) &
                 / (FingasDay - StartTime) ) * IDT), Fingas_Class(StartClass,4))
        ELSEIF (ElapsedTime <= FingasWeek) THEN !Form Time + 1 Day <= emulsion time <= 1 Week
            WaterContent = max(WaterContent + (( (Fingas_Class(StartClass, 2) &
                 - Fingas_Class(StartClass, 1)) / (FingasWeek - FingasDay) ) * IDT), &
                 Fingas_Class(StartClass,2))
            xViscEmul = min(xViscEmul + (( (Fingas_Class(StartClass, 5) -  &
                Fingas_Class(StartClass, 4)) / (FingasWeek - FingasDay) )  &
                * IDT), Fingas_Class(StartClass,5))
        ELSEIF (ElapsedTime <= FingasYear) THEN !emulsion time <= 1 Year (+ FormTime)
            WaterContent = max(WaterContent + (( (Fingas_Class(StartClass, 3)  &
                - Fingas_Class(StartClass, 2)) / (FingasYear - FingasWeek) ) * &
                 IDT), Fingas_Class(StartClass,3))
            xViscEmul = min(xViscEmul + (( (Fingas_Class(StartClass, 6) -  &
                Fingas_Class(StartClass, 5)) / (FingasYear - FingasWeek) ) * &
                 IDT), Fingas_Class(StartClass,6))
        ELSE
            WaterContent = Fingas_Class(StartClass,3)
            xViscEmul = Fingas_Class(StartClass,6)
        END IF
       ! write(*,*)StartClass,'EMUL ClassIndex else WatCont=',WaterContent,' xViscEmul=',xViscEmul
    END IF

    RETURN

    END SUBROUTINE Emulsify
    !*****************************
    !*     Subroutine Evaporate  *
    !*****************************
    SUBROUTINE Evaporate(ElapsedTime,FirstAP)
    USE PARAM_MOD, ONLY: Uwind_10,Vwind_10,EvapOption,idt,VolumeSpill,Oil_Resin,Oil_Asph, &
                         Fingas_B,Fingas_T,Fingas_TYP  ! c.laurent-OGS:  Fingas_B,Fingas_T,Fingas_TYP 
                         ! reintroduced to avoid using fingas with fixed coefficients 
    IMPLICIT NONE

    !I/O variables
    INTEGER, INTENT(IN):: ElapsedTime
    LOGICAL, INTENT(IN):: FirstAP

    !Local variables
    DOUBLE PRECISION:: PercentEvap            !percent of oil evaporated in timestep
    DOUBLE PRECISION:: PrevPercentEvap            !percent of oil evaporated at previous timestep ! added by OGS as PercentEvap was decreasing with time
    DOUBLE PRECISION:: CumPercentEvap        !cumulative percentage of oil evaporated to date
    DOUBLE PRECISION:: PercentDist            !percent of oil mass distilled at 180degC (FINGAS only)
    INTEGER:: n                    !counter
    DOUBLE PRECISION:: Mol(16)                !molecular weight of pseudocomponent
    DOUBLE PRECISION:: MolFrac(16)            !molar fraction of pseudocomponent
    DOUBLE PRECISION:: VolFrac(16)            !volume fraction of pc
    DOUBLE PRECISION:: AvgMW                !average molecular weight
    DOUBLE PRECISION:: Ke                    !mass transfer coefficient (m/s)
    DOUBLE PRECISION:: VEvap(16)            !volume of each pc evaporated
    DOUBLE PRECISION:: dTdFe                !rate of change of temperature versus fraction evaporated
    DOUBLE PRECISION:: InitBP                !initial boiling point
    INTEGER:: nPC                !number of pseudocomponent
    !DOUBLE PRECISION :: WindSpeed

    !save values on exit
    SAVE:: CumPercentEvap, dTdFe, InitBP, nPC, PercentDist,PrevPercentEvap

    IF (FirstAP) THEN
        CALL InitialEvap(dTdFe, InitBP, nPC)
        PercentEvap = 0.0
        CumPercentEvap = 0.0
        MassEvap = 0.0
        VolumeEvap = 0.0
        PercentDist = 16.0        !Need to either read in value or calculate a propoer value
    END IF


    SELECT CASE (F_UpCase(EvapOption))

        CASE ("FINGAS")
            !Fingas, M. 1997
            !"The Evaporation of Oil Spills: Prediction of equations using distillation data"
            !Arctic and Marine OilSpill Program Technical Seminar, Environment Canada. 1997 Vol1:20 pp1-20

            !Assume Logarithmic (vast majority of oil types
!            PercentEvap = ((0.165 * PercentDist) + (0.045 * (WaterTempInst - 15.0))) * log(REAL(ElapsedTime,KIND(1)) / 60.0) &
!                            * (1.0 - WaterContent)
            SELECT CASE (Fingas_TYP)
                CASE (1)
                    !Logarithmic
                    PercentEvap = max(PrevPercentEvap,  &  ! Added by OGS to prevent PercentEvap to decrease in time
                            (Fingas_B + (Fingas_T * (WaterTempInst-15.0)))   &
                            * log(REAL(ElapsedTime,KIND(1)) / 60.0) &
                            * (1.0-WaterContent) )
                CASE (2)
                    !Sqrt
                    PercentEvap =  max(PrevPercentEvap,  &  ! Added by OGS to prevent PercentEvap to decrease in time
                              (Fingas_B + (Fingas_T * (WaterTempInst-15.0)))   &
                            * sqrt(ElapsedTime / 60.0)          &
                            * (1.0-WaterContent)   )
                CASE DEFAULT
                    write(*,*)'Case not encoded'
                    stop

            END SELECT
            PrevPercentEvap=PercentEvap
            !SquareRoot only applies to a few refined products - eqn below (may include later)
            !CumPercentEvap = ((0.0254 * PercentDist) + (0.01 * (WaterTempInst - 15.0))) * sqrt(ElapsedTime / 60.0)
            MassEvap = MassSpill * (PercentEvap / 100.0)
            VolumeEvap = MassEvap / RhoOil

           ! write(*,*)'ResinOil<',Oil_Resin,MassEvap,MassSpill
            ResinOil = Oil_Resin / (1.0 - (MassEvap/MassSpill))    !Increase in resin % (Assuming no evaporation)
            AsphOil = Oil_Asph / (1.0 - (MassEvap/MassSpill))    !Increase in asphaltene % (Assuming no evaporation)

!        CASE ("PSEUDO")
!            !ADIOS2
!            DO n = 1, nPC
!                Mol(n) = Vol(n) / Vbar(n)        !moles in volume of pseudo component
!            end do
!
!            do n = 1, nPC
!                IF (Vol(n) == 0.0)THEN
!                    MolFrac(n) = 0.0
!                    VolFrac(n) = 0.0
!                ELSE
!                    MolFrac(n) = Mol(n) / sum(Mol)
!                    VolFrac(n) = Vol(n) / sum(Vol)
!                    AvgMW = AvgMW + MolFrac(n)*MW(n)
!                END IF
!            END DO
!            Ke = 0.0048 * WindSpeed**(7.0/9.0) * (1.3676 * (sqrt(0.018/AvgMW))**(2.0/3.0)) * (lenR*2.0)**(-1.0/9.0)
!            DO n = 1, nPC
        !                VEvap(n) = DT * (Ke * VolumeOil * VapP(n) * Vbar(n) * VolFrac(n)) / &
!                                (R * Slickthickness * (WaterTempInst+273.15))
!                IF (Vol(n) - VEvap(n) > 0) THEN
!                    Vol(n) = Vol(n) - VEvap(n)
!                ELSE
!                    Vol(n) = 0.0
!                END IF
!            END DO
!            VolumeEvap = VolumeEvap + sum(VEvap)
!            MassEvap = VolumeEvap * RhoOil

        CASE ("MACKAY")
            !Stiver W., Mackay, D. 1984
            !"Evaporation rate of spills of hydrocarbons and petroleum mixtures"
            !Environmental Science and Technology, 1984. vol 18, pp 834-480
            !as modified by ADIOS2 (NOAA 2000)
            !WindSpeed = sqrt((Uwind_10**2.0) + (Vwind_10**2.0))
            Ke = 1.5E-3 * WindSpeed**0.78
            PercentEvap = ((Ke * AreaOil * idt) / VolumeSpill) * &
                            exp(6.3 - ((10.3 * (InitBP + (dTdFe * CumPercentEvap))) / (WaterTempInst + 273.15))) &
                            * (1.0 - WaterContent)

            CumPercentEvap = CumPercentEvap + PercentEvap        !Cumulative percentage evaporated

!            VolumeEvap  = VolumeSpill * CumPercentEvap            !Total volume evaporated to date
!            MassEvap = VolumeEvap * RhoOil                        !Total mass evaporated to date
            MassEvap = MassSpill * CumPercentEvap  ! changed by OGS as Mass is conserved, not volume 
            VolumeEvap  = MassEvap / RhoOil        ! and all variables depend on mass variations  

           ! write(*,*)'ResinOil<',Oil_Resin,MassEvap,MassSpill
            ResinOil = Oil_Resin / (1.0 - (MassEvap/MassSpill))    !Increase in resin % (Assuming no evaporation of resins)
            AsphOil = Oil_Asph / (1.0 - (MassEvap/MassSpill))    !Increase in asphaltene % (Assuming no evaporation of asphaltenes)

        CASE DEFAULT
            VolumeEvap = 0.0
            MassEvap = 0.0
            WRITE(*,*)'Case not encoded'
            WRITE(*,*)'No EVAPORATION processes modelled'

    END SELECT

    END SUBROUTINE Evaporate

    !*******************************
    !*     Subroutine InitialEvap  *
    !*******************************
    SUBROUTINE InitialEvap(dTdFe, InitBP, nPC)
    !*** CRUDE OILS ONLY ***!
    !if refined oil is chosen then need to write code to automatically select a crude oil
    !based on common density and viscosity values (at the least) and echo selection to user
    USE PARAM_MOD, ONLY: API,EvapOption

    IMPLICIT NONE
    !I/O variables
    INTEGER, OPTIONAL, INTENT(IN):: nPC
    DOUBLE PRECISION, INTENT(OUT):: dTdFe
    DOUBLE PRECISION, INTENT(OUT):: InitBP

    !Local variables
    INTEGER:: n                    !counter
    INTEGER:: nPts                !number of interpolation points
    DOUBLE PRECISION:: Cut_T(17)            !temperature of cut (degC)
    DOUBLE PRECISION:: Cut_F(17)            !volume fraction of cut
    DOUBLE PRECISION:: deltazb                !deltaZb value in Antoines eqn (ADIOS, 2001)
    DOUBLE PRECISION:: C2(16)                 !C2 value array in Antoines eqn (ADIOS, 2001)
    DOUBLE PRECISION:: deltaS(16)            !deltaS value array in Antoines eqn (ADIOS, 2001)

    deltazb = 0.97

    SELECT CASE (F_UpCase(EvapOption))

    CASE ("MACKAY")
        !ADIOS1 (1994)
        InitBP = 532.98 - 3.1295 * API
        dTdFe = 985.62 - 13.597 * API

!        CASE ("PSEUDO")
!        *****************************************
!        **** THERE IS SOMETHING WRONG HERE   ****
!        **** EVAPORATION X10 TIMES TOO FAST  ****
!        **** NEEDS FURTHER INVESTIGATION     ****
!        *****************************************

!            !ADIOS2 (2000)
!             InitBP = 456.17 - 3.3447 * API
!            dTdFe = 1356.7 - 247.36 * log(API)
!            !CONSTRUCT PSEUDO-COMPONENTS
!            !if 2 or more distillation cuts available in seleceted oil then use distillation data
!            !otherwise use simple correlation scheme below based on constructing 5 pseudo components:
!            IF (NumCuts < 2) THEN
!                nPC = 5
!                DO n = 1,nPC
!                    BP(n) = InitBP + (dTdFe * (n-0.5)/ nPC)
!                    VolInit(n) = VolumeSpill / nPC
!                END DO
!               ELSE
!                   nPts = NumCuts + 2
!                   Cut_T(1) = Cut_Temp(1) - Cut_Frac(1) * ((Cut_Temp(2) - Cut_Temp(1)) / (Cut_Frac(2) - Cut_Frac(1)))
!                   Cut_F(1) = 0.0
!
!                   DO n = 2, nPts-1
!                    Cut_T(n) = Cut_Temp(n-1)
!                    Cut_F(n) = Cut_Frac(n-1)
!                END DO
!
!                Cut_T(nPts) = Cut_Temp(NumCuts) + (1 - Cut_Frac(NumCuts)) * &
!                                ((Cut_Temp(NumCuts) - Cut_Temp(1)) / (Cut_Frac(NumCuts) - Cut_Frac(1)))
!                Cut_F(nPts) = 1.0
!                nPC = nPts - 1
!
!                DO n = 1, nPC
!                    BP(n) = Cut_T(n+1)
!                    VolInit(n) = (Cut_F(n+1) - Cut_F(n)) * VolumeSpill
!                END DO
!            END IF
!
!            DO n = 1, nPC
!                C2(n) = (0.19 * BP(n)) - 18.0
!                deltaS(n) = 8.75 + 1.987 * log10(BP(n))
!                Vbar(n) = 7.0E-05 - (2.102E-07 * BP(n)) + (1.0E-09 * (BP(n)**2.0))
!                MW(n) = 0.04132 - (1.985E-04 * BP(n)) + (9.494E-07 * (BP(n)**2.0))
!                VapP(n) = Pa * exp( ( (deltaS(n) * (BP(n) - C2(n))**2.0) / (deltaZb * R * BP(n)) ) * &
!                                    ( (1.0/(BP(n) - C2(n))) - (1.0/((WaterTempInst+273.15)-C2(n)))   )   )
!                Vol(n) = VolInit(n)
!            END DO

    CASE ("FINGAS")
        !do nothing
        !no pre-processing required (as of yet)
        !may look at calculating PercentDistCut at some stage

    CASE DEFAULT
        WRITE(*,*)'Case not encoded'
        WRITE(*,*)'No EVAPORATION processes modelled'

    END SELECT

    END SUBROUTINE InitialEvap

    !*****************************
    !*     Subroutine Disperse   *
    !*****************************
    SUBROUTINE Disperse(ElapsedTime,FirstAP,WaterDepth)
    !French-McKay, 2004
    !"Oil Spill Impact Modelling: Development and Validation"
    !Environmental Toxicology and Chemistry, Vol 23, No. 10, pp 2441-2456. SETAC, USA
    !after
    !Delvigne, G., & Sweeney, C. 1988
    !"Natural Dispersion of Oil"
    !Oil and Chemical Pollution, Vol 4, pp 281-310

    USE PARAM_MOD,    ONLY:    SigWaveHeight,SigWaveLength,SigWavePeriod,idt,     &
                            UWind_10,VWind_10,Wind,windwavemodel
    IMPLICIT NONE

    !I/O variables
    INTEGER, INTENT(IN):: ElapsedTime
    LOGICAL, INTENT(IN):: FirstAP
    DOUBLE PRECISION, INTENT(IN):: WaterDepth

    !Local variables
    DOUBLE PRECISION:: Dbwe            !dissipated breaking wave energy per unit area (J/m2)
    DOUBLE PRECISION:: Zintrude        !intrusion depth (m)
    DOUBLE PRECISION:: RelDepth        !relative depth used to determine deep / shallow water
    DOUBLE PRECISION:: Hbreak        !breaking wave height (m)
    DOUBLE PRECISION:: Cstar        !empirical entrainment constant
    DOUBLE PRECISION:: FracWave        !fraction of sea surface hit by breaking waves
    DOUBLE PRECISION:: Oil_d50        !mean oil droplet diameter (um)
    DOUBLE PRECISION:: DropDiam(10)    !Do droplet diameter per size class (m)
    DOUBLE PRECISION:: Qd(10)        !Entrainment rate per size class (kg/m2s)
    DOUBLE PRECISION:: Zmix(10)        !mixing layer depth per size class(m)
    DOUBLE PRECISION:: Qdtotal        !Total entrainment rate for all size classes
    DOUBLE PRECISION:: DeltaDiam    !oil droplet interval diameter (m)
    DOUBLE PRECISION:: Dv            !vertical disperson (m2/s) (Thorpe)
    DOUBLE PRECISION:: RiseVel(10)    !droplet terminal rise velocity (m/s)
    DOUBLE PRECISION:: Rmin            !minimum and maximum droplet radius (m)
    DOUBLE PRECISION:: Rmax            !minimum and maximum droplet radius (m)
    !DOUBLE PRECISION:: WindSpeed
    INTEGER:: i            !droplet distribution array counter
    DOUBLE PRECISION:: kw1            !kw for Re < 50
    DOUBLE PRECISION:: kw2            !kw for Re > 50
    DOUBLE PRECISION, PARAMETER:: Uth = 6.0        !Threshold wind speed for the onset of breaking waves (m/s)
    DOUBLE PRECISION, PARAMETER:: Ewave = 5000.0     !mean energy dissipation rate per unit volume (J/m3-s) ranges between 1,000 and 10,000 (J/m3.s) - See Delvigne (1988)
    DOUBLE PRECISION, PARAMETER:: FracOil = 1.0    !fraction of sea surface covered by oil.
    DOUBLE PRECISION:: waveheight,waveperiod
    IF (FirstAP) THEN
        MassDisp = 0.0        !initialise mass dispersed
        VolumeDisp = 0.0    !initialise volume dispersed
        !Waterdepth = 10.0
    END IF
    if(windwavemodel)then
            waveheight = SigWaveHeight        
            waveperiod = SigWavePeriod
            !SWAN model output
    else
            waveheight = 0.243 * (0.71*WindSpeed**1.23)**2.0 / Gravity                                                        !ADIOS formulation (NOAA 1994)
            waveperiod = 8.13 * WindSpeed / gravity
    end if       

    !Calculate breaking wave height based on deep / shallow water wave conditions
    !Actual HBREAK value will be read from SWAN model and the shallow/deep eqns removed
    RelDepth = WaterDepth / SigWaveLength
    IF (RelDepth < 0.05 ) THEN                                !Deep Water
        Hbreak = 0.2185 * (waveperiod**2.0)
    ELSE IF (RelDepth > 0.50) THEN                            !Shallow Water
        Hbreak = 0.78 * WaterDepth
    ELSEIF (RelDepth <= 0.50 .AND. RelDepth >= 0.05)THEN    !Intermediate waters
        Hbreak = 1.5 * waveheight                            !Default values used in SeatrackWeb technical documentation
    END IF

    !Calculate dissipated breaking wave energy (J/m2)
    Dbwe = 0.0034 * RhoWater * Gravity * (Hbreak**2.0)

    !Calculate fraction of sea surface hit by breaking waves
    !WindSpeed = sqrt((Uwind_10**2.0) + (Vwind_10**2.0))
    IF (WindSpeed <= Uth ) THEN
        FracWave = 3E-06 * (WindSpeed**3.5 / waveperiod)
    ELSE
        FracWave = 0.032 * ( (WindSpeed - Uth) / waveperiod)
    END IF

    !Calculate intrusion depth
    Zintrude = 1.5 * Hbreak

    !Calculate vertical dispersion coefficient (m2/s)
    Dv = 0.0015 * WindSpeed

    !Calculate empirical entrainment constant Cstar
    IF (ViscOil < 132.0 ) THEN
        Cstar = exp( (-0.1023 * log(ViscOil)) + 7.572)
    ELSE
        Cstar = exp( (-1.8927 * log(ViscOil)) + 16.313)
    END IF

    !Calculate mean oil droplet diameter (um), minimum radius (m) and maximum radius (m)
    Oil_d50 = 1818.0 * (Ewave**(-0.5)) * (ViscOil**0.34) !based on viscosity
    Rmin = 0.1 * Oil_d50 * 1E-06 !(convert from micrometers to meters)
    Rmax = Oil_d50 * 1E-06 !(convert from micrometers to meters)

    !Construct droplet size distribution
    !Adopt 10No. size classes between Rmin and Rmax, equally spaced on diameter
    DeltaDiam = ((Rmax * 2.0) - (Rmin * 2.0)) / 10.0

    !initialise total entrainment rate every timestep
    Qdtotal = 0.0

    DO i = 1, 10
        !for each droplet interval, calculate centred droplet diameter, Do
        DropDiam(i) = ((2.0 * Rmin) + (0.5 * DeltaDiam)) + (DeltaDiam * (i-1))
        !calculate entrainment rate for each droplet size class
    !write(*,*)'Disp',ViscOil,Dbwe,FracWave,WindSpeed,'|', &
    !        waveperiod,RhoWater,Gravity,Hbreak
         Qd(i) = Cstar * Dbwe * FracOil * FracWave * DropDiam(i)**0.7 * DeltaDiam !(kg/ms)
        !sum over all droplet classes for total entrainment rate
        Qdtotal = Qdtotal + Qd(i)
    END DO
    !calculate mass dispered (kg) and volume dispersed (m3)
    !write(*,*)' Dispersion :',MassDisp,Qdtotal,AreaOil,RhoOil
    MassDisp = MassDisp + (Qdtotal * AreaOil * idt)
    VolumeDisp = MassDisp / RhoOil

    !Calculate critical droplet diameter Proctor et al (1994)


!    !Check Reynolds criteria according to Tkalich & Chan (2002)
!    kw1 = 2.0 * Gravity * (1 - (RhoOil/RhoWater)) / (9.0 * KinViscWater)
!    kw2 = sqrt(16.0 * Gravity * (1 - (RhoOil/RhoWater)) / 3.0)
!
!    do i = 1, 10
!        if (2.0 * kw1 * ((DropDiam(i)/2.0)**3.0) / KinViscWater <= 50.0) then        !Low Reynolds => use Stokes Law
!            RiseVel(i) = kw1 * (DropDiam(i)/2.0)**2.0
!            Zmix(i) = max(Dv/RiseVel(i), Zintrude)
!        end if
!        if (2.0 * kw2 * ((DropDiam(i)/2.0)**(3.0/2.0)) / KinViscWater >= 50.0)then    !High Reynolds => use Reynolds Law
!            RiseVel(i) = kw2 * (DropDiam(i)/2.0)**0.5
!            Zmix(i) = max(Dv/RiseVel(i), Zintrude)
!        end if
!    end do

    !STILL TO FINISH


    END SUBROUTINE Disperse

    !********************************
    !*     Subroutine StokesDrift   *
    !********************************
    SUBROUTINE STOKESDRIFT(uwind,vwind,pdir,depth,ustoke,vstoke,WVecAngle)
    ! Give the Stokes drift at depth of a particle assuming a Pierson Moskowitz spectrum and Langmuir circulation
    ! Author Marcel Cure (www.numericswarehouse.com) Apr. 2010
    !
    DOUBLE PRECISION, INTENT(IN) :: uwind,vwind,pdir,depth
    DOUBLE PRECISION, INTENT(OUT) :: ustoke,vstoke                       ! velocity components of particle at depth
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: WVecAngle
    !
    DOUBLE PRECISION :: V10,knum,phi,VS0,VS
    DOUBLE PRECISION, PARAMETER :: grav=9.81
    !
    !porint *,'Calculating Stokes Drift'
    !
    V10 = SQRT(uwind**2.0 + vwind**2.0)
    !
    ! We use the modified Stokes profile according to eqn 9 Carniel, S. Sclavo, M., Kantha, L.H. and C.A. Clayson 2005
    ! ' Langmuir cells and mixing in the upper ocean', Il Nuevo Cimento 28 33-54
    ! KLC: mag of Stokes drift vel |VS| = VS0*EXP(2kz) (top of p.35)
    ! KLC: eqns for VS0 and knum are shown in eqn(9) p.42
    !
    ! *** Direct use of WVecAngle added to Zlev version by OGS ***    
    if(present(WVecAngle))then
      phi = WVecAngle       
    else
    ! *** End of OGS modifications for direct use of WVecAngle ***    
      phi = F_WindAngle(pdir)   ! *** WARNING formulations used by F_WindAngle must be verified *** 
    endif
    !
    knum = 1.25*(grav/(V10**2.0))
    VS0 = 0.016*V10
    VS = VS0*EXP(2.0*knum*depth)
    ustoke = VS*COS(phi)
    vstoke = VS*SIN(phi)
    !
    END SUBROUTINE STOKESDRIFT

!*********************************************************************
!        *** BEGINNING OF supplements added to Zlev version by OGS ***    
!*********************************************************************
    ! Remove_Beached_Oil_from_Weathering_Oil routine is not yet validated
    SUBROUTINE Remove_Beached_Oil_from_Weathering_Oil(nParWater,nParBeaching)

        USE PARAM_MOD, ONLY: VolumeSpill

        implicit none

        INTEGER, INTENT(IN) :: nParWater 
        INTEGER, INTENT(IN) :: nParBeaching
        DOUBLE PRECISION :: MassOil_LeftAtSurf,MassOil_Beaching

        MassOil_LeftAtSurf=MassSpill-MassEvap-MassDisp

!       MassOil_Beaching=MassOil           * &
        MassOil_Beaching=MassOil_LeftAtSurf* & 
                           (DBLE(nParBeaching)/DBLE(nParWater+nParBeaching)) 

        MassBeached    = MassBeached + MassOil_Beaching 
        VolumeBeached  = MassBeached / RhoOil


    END SUBROUTINE Remove_Beached_Oil_from_Weathering_Oil

    SUBROUTINE InitOilOutputs()

        USE PARAM_MOD,     ONLY:  outpathGiven,outpath,NCOutFile
        IMPLICIT NONE

        character(250)::Ratesfilename,OilPropfilename

           ! -------------------------------
           if(outpathGiven)then
             Ratesfilename  = TRIM(outpath) //'/'//TRIM(NCOutFile)//'-OilRates.csv'
           else
             Ratesfilename  = TRIM(NCOutFile)//'-OilRates.csv'
           endif 
           open(997,file=TRIM(Ratesfilename),status='replace')
           write(997,"(a12,4(',',a15))")'s','AvWindSpeed', &
                  ' Oil Viscosity','Oil Density'   
           close(997)   

           ! -------------------------------
           if(outpathGiven)then                                                   
             OilPropfilename = TRIM(outpath) //'/'//TRIM(NCOutFile)//'-OilPropsOut.csv'     
           else                                                                   
             OilPropfilename = TRIM(NCOutFile)//'-OilPropsOut.csv'                
           endif                                                                  
           open(998,file=TRIM(OilPropfilename),status='replace')                  
           write(998,"(a10,14(',',a15))")'time','VolumeOil','VolumeEvap','VolumeDisp', &
                         'VolumeBeached','WaterContent', &                              
                         'RhoOil','ViscOil','AreaOil','OilThickness',   &                           
                         'VolumeSlick','MassOil','AverageWind','WaterTemp(oC)',&                         
                          'AvWindAngle'                                
           close(998)      
           ! -------------------------------

    END SUBROUTINE InitOilOutputs

    SUBROUTINE WriteOilOutputs()
        USE PARAM_MOD,     ONLY:  outpathGiven,outpath,NCOutFile,idt,numpar
        use convert_mod, only: x2lon,y2lat

        implicit none


        character(250)::Ratesfilename,OilPropfilename
        INTEGER :: ElapsedTime        !time in seconds

        ElapsedTime = iT * idt - Phase1Time                    !Time since Phase1Time reached (s)
       ! write(*,*)'-----------------------------------------------------'
       ! write(*,*)'*** output Wind =',WindSpeed,' Temp=',WaterTempInst, &
       !           ' Angle=',AvWindAngle,'  ***'
       ! write(*,*)'-----------------------------------------------------'
           ! -------------------------------
           if(outpathGiven)then
             Ratesfilename  = TRIM(outpath) //'/'//TRIM(NCOutFile)//'-OilRates.csv'
           else
             Ratesfilename  = TRIM(NCOutFile)//'-OilRates.csv'
           endif 
           open(997,file=TRIM(Ratesfilename),position='append')
           write(997,"(I12,',',4(f15.9,','))")oil_time,WindSpeed, &
                     ViscOil,RhoOil   
           close(997)    
           ! -------------------------------
           if(outpathGiven)then                                                   
             OilPropfilename = TRIM(outpath) //'/'//TRIM(NCOutFile)//'-OilPropsOut.csv'     
           else                                                                   
             OilPropfilename = TRIM(NCOutFile)//'-OilPropsOut.csv'                
           endif                                                                  
           open(998,file=TRIM(OilPropfilename),position='append')                  

          write(998,99)oil_time,VolumeOil,VolumeEvap,VolumeDisp, &
                       VolumeBeached,WaterContent, &
                       RhoOil,ViscOil,AreaOil,OilThickness,&
                       VolumeSlick,MassOil,WindSpeed,WaterTempInst,AvWindAngle
          close(998)
           ! -------------------------------
        99 format(I10,',',14(f15.6,','),f15.6)

    END SUBROUTINE WriteOilOutputs

!        *** END of supplements of the Zlev version by OGS ***

END MODULE OIL_MOD
