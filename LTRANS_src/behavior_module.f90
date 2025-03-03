MODULE BEHAVIOR_MOD

! The behavior module is used to assign biological or physical characteristics
!   to particles. 
! Currently particle movement is in the vertical direction. 
!
! Particle characteristics can include a swimming/sinking speed component and 
! a behavioral cue component that can depend upon particle age. The swimming/
! sinking speed component controls the speed of particle motion and can be 
! constant or set with a function. The behavioral cue component regulates the 
! direction of particle movement. For biological behaviors, a random component 
! is added to the swimming speed and direction to simulate random variation in
! the movements of individuals (in behavior types 1 - 5, see list below). 
! Physical characteristics can also be assigned to particles, like constant 
! sinking velocity, without the additional random movements (behavior type 6). 
! The following behavior types are currently available in LTRANS and are 
! specified using the Behavior parameter in the LTRANS.inc file:
!
!
! Passive (no behavior): Behavior = 0. In this case, the behavior module is not
!   executed. Particle motion is based on advection, and, if turned on, 
!   horizontal and vertical turbulence.
!
! Near-surface orientation: Behavior = 1. Particles swim up if they are deeper 
!   than 1 m from the surface.  
!
! Near-bottom orientation: Behavior = 2. Particles swim down if they are 
!   shallower than 1 m from the bottom.  
!
! Diurnal vertical migration: Behavior = 3. Particles swim down if light levels 
!   at the particle location exceed a predefined threshold value.  
!
! Crassostrea virginica oyster larvae: Behavior = 4. Swimming speeds and 
!   direction of motion vary depending upon age (stage) according to field and 
!   laboratory observations (see North et al. 2008). 
!
! C. ariakensis oyster larvae: Behavior = 5. Swimming speeds and direction of 
!   motion vary depending upon age (stage) according to field and laboratory 
!   observations (see North et al. 2008).
!
! Sinking velocity: Behavior = 6. Particles move up or down with constant
!   sinking (or floating) speeds without individual random motion. Code that 
!   calculates salinity and temperature at the particle location is included 
!   (but commented out) as a basis for calculating density-dependent sinking 
!   velocities.  
!
! Tidal Stream Transport: Behavior = 7.
!
! Nephrops Norvegicus : Behavior = 8.
! Solea Solea         : Behavior = 9.
! Mullus Barbatus     : Behavior = 10.
! Parameterizable larvae  : Behavior = 11.
! Oyster Ostrea Edulis larvae  : Behavior = 12.
!
! Behavior algorithms and code created by: Elizabeth North
! Module structure created by:             Zachary Schlag
! Created on:                              2004
! Last Modified on:                        22 March 2011
!

  IMPLICIT NONE
  PRIVATE
  SAVE

  !Timer for C. ariakensis downward swimming behavior
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: timer

  !Behavior of each particle
  INTEGER, ALLOCATABLE, DIMENSION(:) :: P_behave


  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:  ) ::    &
    P_pediage, & !Age at which the particle will settle (become a pediveliger)
    P_deadage, & !Age at which the particle will stop moving (die)
      !The following are for calculating salt gradient:
    P_Sprev,   & !Salinity at particle's previous location 
    P_zprev      !Particle's previous depth 

  DOUBLE PRECISION ::   maxtimeatsurf,maxsizeatsurf 


  !Swimming speed (age-dependent, linear increase unless constant)   
  !(n,1)slope, (n,2)intercept, (n,3) speed at current age
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: P_swim

  !For behavior 8, Larval Size increasing according to growth function
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: LarvSize

  !For behavior 8 , 10 and 11 swdown radiation for dvm 
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: swdown
  INTEGER :: swdown_numrec
  !For behavior 7, tracks if particle is on the bottom
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: bottom

  !Tracks if the particle is dead (TRUE) or alive (FALSE)
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: dead

  !Tracks if particles are Out Of Bounds (ie cross open ocean bound)
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: oob

  !The following procedures have been made public:
  PUBLIC :: initBehave,updateStatus,behave,getStatus,finBehave, &
            setOut,isOut,die,isDead

CONTAINS

  SUBROUTINE initBehave()    !Initialize the behavior module
    USE PARAM_MOD, ONLY: numpar,Behavior,swimfast,swimslow,swimstart,     &
                      pediage,deadage,Sgradient,settlementon, &
                      swdown_ASCIIfname,swdown_dt,swdown_rec,swdown_ASCII,days,&
                      DVMtime,SettlementSize,stranding_on!,readNetCdfSwdown
    USE SETTLEMENT_MOD, ONLY: initSettlement
    USE STRANDING_MOD, ONLY: initStranding
    USE NORM_MOD,   ONLY: norm
    IMPLICIT NONE
    INTEGER :: n,count
    REAL :: swdown_read

    write(*,*) 'initialize behavior'    

    !Allocate Behavior Variables
    ALLOCATE(timer(numpar))
    ALLOCATE(P_behave(numpar))
    ALLOCATE(P_pediage(numpar))
    ALLOCATE(P_deadage(numpar))
    ALLOCATE(P_Sprev(numpar))
    ALLOCATE(P_zprev(numpar))
    ALLOCATE(P_swim(numpar,3))
    ALLOCATE(dead(numpar))
    ALLOCATE(oob(numpar))
    ALLOCATE(LarvSize(numpar))
    IF(Behavior == 7) THEN
      ALLOCATE(bottom(numpar))
      bottom = .TRUE.
    ENDIF
    !IF(swdown_ASCII .and. (.not. readNetCdfSwdown)) THEN
    IF(swdown_ASCII) THEN
      swdown_numrec=ceiling(days*86400./float(swdown_dt)+1)
      write(*,*)swdown_numrec,' records to read from file ',trim(swdown_ASCIIfname)
      write(*,*)' start reading from record ',swdown_rec
      ALLOCATE(swdown(swdown_numrec))
      count=1
      OPEN(unit=88,FILE=trim(swdown_ASCIIfname), status='old',    &
             access='sequential', form='formatted', action='read' )
      do
        if(count.gt.swdown_numrec+swdown_rec-1)exit
        READ(88,*) swdown_read
        !write(*,*)count,swdown_read,count.ge.swdown_rec
        if(count.ge.swdown_rec)then
          swdown(count-swdown_rec+1)=swdown_read
        endif
        count=count+1
      enddo
      !write(*,*)swdown
      CLOSE(88)
    ENDIF  

    maxtimeatsurf=DVMtime
    maxsizeatsurf=SettlementSize
    do n=1,numpar
      !Set behavior to the one specified in LTRANS.inc
      P_behave(n) = Behavior  !Behavior
      P_pediage(n) = pediage  !age at which particle reaches maximum swimming
                              !speed and can settle (becomes a pediveliger) (s)
      P_deadage(n) = deadage  !age at which particle stops moving (dies) (s)
      !Note: the following code assigns different veliger and pediveliger
      !  stage durations
      !P_pediage(n) = (14. + norm()*0.5)*24.*3600.
      !P_deadage(n) = P_pediage(n) + (7. + norm()*0.5)*24.*3600.

      !Calculate slope and intercept for age-dependent linear swimming speed
      if(abs(P_pediage(n) - swimstart)<1e-7)then
         write(*,*)'error pediage must be greater than swimstart'
         stop 'error increate pediage or decrease swimstart'
      endif
      P_swim(n,1) = (swimfast - swimslow)/(P_pediage(n) - swimstart) !slope
      P_swim(n,2) = swimfast - P_swim(n,1)*P_pediage(n)              !intercept
      P_swim(n,3) = 0.0                                  !swimming speed (m/s)
      IF(Behavior.eq.8) then  ! Nephrops Norvegicus
         LarvSize(n) = 6.0
         P_behave(n) = 0
         if(DVMtime==0) then
           maxtimeatsurf=100*86400
           write(*,*)'DVMtime not given by the user, taking ',maxtimeatsurf
         else 
           maxtimeatsurf=DVMtime    
         endif
         if(settlementSize==0) then
            maxsizeatsurf=14
         else 
            maxsizeatsurf=settlementSize
         endif
      ELSEIF(Behavior.eq.9) then ! Solea Solea
         LarvSize(n) = 1.0
         P_behave(n) = 0
         if(DVMtime==0) then
           maxtimeatsurf=30*86400
           write(*,*)'DVMtime not given by the user, taking ',maxtimeatsurf
        else 
           maxtimeatsurf=DVMtime    
         endif
         if(settlementSize==0) then
           maxsizeatsurf=8 
         else 
           maxsizeatsurf=SettlementSize
         endif
      ELSEIF(Behavior.eq.10) then ! Mullus Barbatus
           P_behave(n) = 0
           maxtimeatsurf=35*86400
           maxsizeatsurf=9999
      ELSEIF(Behavior.eq.11) then ! Parameterizable larvae
           P_behave(n) = 0
      ENDIF
      if(n==0)&
        write(*,*)'Behavior ',Behavior,' intial larvae size=',LarvSize(n),&
        ' maxtimeatsurf=',maxtimeatsurf,' maxsizeatsurf=',maxsizeatsurf
      !Note: P_swim(n,3) is updated at each time step in Subroutine behave
    enddo

    IF(Behavior.eq.12) then  ! Ostrea Edulis Oyster
     call Ostrea_Edulis_initialize()
    ENDIF

    if(maxtimeatsurf==0 )maxtimeatsurf=999999999
    !The following variables are used by the C. virginica and C. ariakensis 
    !  behavior routines
    timer = DBLE(0.0)         !to count how long C. arikensis particles swim down

    ! Initialize salt storage matrices 
    P_Sprev = 0.0       !Initialized to 0.0
    P_zprev = 0.0       !Initialized to 0.0

    ! Initialize dead to .FALSE. i.e. all particles are initially alive
    dead = .FALSE.

    ! Initialize out of bounds tracker to .FALSE.
    !   (i.e. all particles start in bounds)
    oob = .FALSE.

    !if Settlement is turned on then inform Settlement module of the age at 
    !  which particle can settle (i.e., become pediveligers)
    if(settlementon)then
      CALL initSettlement(P_pediage)
    endif
    if(stranding_on)then
      CALL initStranding()
    endif

  END SUBROUTINE initBehave

  SUBROUTINE Ostrea_Edulis_initialize()
    USE PARAM_MOD, ONLY: numpar,DVMtime,SettlementSize
    IMPLICIT NONE
    INTEGER :: n
    do n=1,numpar

         LarvSize(n) = 0.17 ! mm

         P_behave(n) = 0

         if(DVMtime==0) then ! maximal time the larvae can spend in pre-settlement phase doing dial vertical migration
           maxtimeatsurf=30*86400
           write(*,*)'DVMtime not given by the user, taking ',maxtimeatsurf
         else 
           maxtimeatsurf=DVMtime    
         endif

         if(settlementSize==0) then
            maxsizeatsurf=0.30 ! mm
         else 
            maxsizeatsurf=settlementSize
         endif

      enddo
  END SUBROUTINE Ostrea_Edulis_initialize

  SUBROUTINE Ostrea_Edulis_behave(Xpar,Ypar,Zpar,Pwc_zb,Pwc_zc,Pwc_zf,P_zb,P_zc,P_zf,  &
                    P_zetac,P_age,P_depth,P_U,P_V,P_angle,P_T,                         &
                    n,it,ex,ix,                                                        &
                    daytime,p,bott,XBehav,YBehav,ZBehav,P_size,Fstlev,klev            )
    USE PARAM_MOD, ONLY: us,dt,idt,twistart,twiend,Em,pi,daylength,Kd,thresh,  &
                    Sgradient,swimfast,swimstart,sink,Hswimspeed,              &
                    Swimdepth,rise,surflayer_upperdepth,surflayer_lowerdepth,  &
                    surflayer_lowerdepth_night,surflayer_upperdepth_night,     &
                    Behavior,Zgrid_depthinterp,BottomLayerThickness,  &
                    swdown_dt,swdown_thresh,Seabed_layerheight, &
                    swdown_ASCII,swdown_t0,Ext0,Zgrid!,readNetCdfSwdown

    USE HYDRO_MOD, ONLY: WCTS_ITPI,getInterp
    USE RANDOM_MOD, ONLY: genrand_real1
#include "VAR_IDs.h"
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: daytime
    DOUBLE PRECISION, INTENT(IN) :: Xpar,Ypar,Zpar,Pwc_zb(:),Pwc_zc(:),        &
                                    Pwc_zf(:),P_zb,P_zc,P_zf,P_zetac,P_age,    &
                                    P_depth,P_U,P_V,P_angle,ex(3),ix(3) !,       &
                                    !P_grainsize ,P_swdown
    INTEGER, INTENT(IN) :: n,it,p,Fstlev,klev
    LOGICAL, INTENT(OUT) :: bott
    DOUBLE PRECISION, INTENT(OUT) :: XBehav,YBehav,ZBehav,P_size
    
    INTEGER :: btest,i,deplvl_minus2,NumInterpLvl,kdown,kup
    INTEGER :: swdown_r1,swdown_r2                                  ! for swdown dvm
    DOUBLE PRECISION :: swdown_alpha,swdown_interp                  ! for swdown dvm
    DOUBLE PRECISION :: targetlayerUPPERdepth,targetlayerLOWERdepth ! for Behavior 8 and 9
    DOUBLE PRECISION :: negpos,dev1,devB,switch,switchslope
    DOUBLE PRECISION :: P_chl,parBehav,Sslope,deltaS,deltaz
    DOUBLE PRECISION :: P_T !not needed unless temperature code below is enabled
    DOUBLE PRECISION :: dtime,tst,E0,P_light
    DOUBLE PRECISION :: currentspeed,Hdistance,theta,X,Y,Growth,P_chl_gradient
    DOUBLE PRECISION :: ws_speed, Chl_up,Chl_down,dZ
     P_size=0.0
     XBehav = 0.0
     YBehav = 0.0
     ZBehav = 0.0
     parBehav = 0.0
     do i=Fstlev+2,us-2
       if ((Zpar .LT. Pwc_zb(i)) .OR.    &
           (Zpar .LT. Pwc_zc(i)) .OR.    &
           (Zpar .LT. Pwc_zf(i))         ) exit
     enddo
     deplvl_minus2 = i-2
     if(Zgrid .and. deplvl_minus2+3>us)then
         NumInterpLvl=us-Fstlev+1
         deplvl_minus2=Fstlev
     else
         NumInterpLvl=4
     endif  
    
    ! Get Chlorophyl at particle position

     P_chl = WCTS_ITPI(VAR_ID_chl,Xpar,Ypar,deplvl_minus2,Pwc_zb,Pwc_zc,Pwc_zf,us,P_zb,    &
                     P_zc,P_zf,ex,ix,p,4,n,NumInterpLvl)
     if(Zpar<Pwc_zc(klev))then
      kup=min(us,klev)
      kdown=max(Fstlev,klev-1)
      if(kdown.eq.kup) kup=kdown+1
     else
      kup=min(us,klev+1)
      kdown=max(Fstlev,klev)
      if(kdown.eq.kup) kdown=kup-1
     endif
     
     Chl_down=getInterp(Xpar,Ypar,VAR_ID_chlc,kdown)
     i=0
     do i=0,us-kup,1
       Chl_up  = getInterp(Xpar,Ypar,VAR_ID_chlc,kup+i)
       if(abs(Chl_up-Chl_down) > 1e-8) exit
     enddo
     if(kup+i>us)then
       Chl_up  = getInterp(Xpar,Ypar,VAR_ID_chlc,kup)
       do i=0,kdown-Fstlev
         Chl_down  = getInterp(Xpar,Ypar,VAR_ID_chlc,kdown-i)
         if(abs(Chl_up-Chl_down) > 1e-8) exit
       enddo
       kdown=max(Fstlev,kdown-i)
     else 
       kup=kup+i
     endif
     dZ=(Pwc_zc(kup)-Pwc_zc(kdown))
     P_chl_gradient = ( Chl_up-Chl_down ) / dZ ! positive means chlorphyl greater upward respect to downward
 
     !Temperature at particle location 
     P_T = WCTS_ITPI(VAR_ID_temp,Xpar,Ypar,deplvl_minus2,Pwc_zb,Pwc_zc,Pwc_zf,us,     &
                     P_zb,P_zc,P_zf,ex,ix,p,4,n,NumInterpLvl)


     ! Compute Growth and increment larval size
     Growth = exp( -8354.9/(P_T+273.15) + 30.7) / 1000.0 ! mm/day
     LarvSize(n)=LarvSize(n)+ Growth * DBLE(idt)/DBLE(86400.0)
     ws_speed=(0.07 * LarvSize(n) + 0.00009 * (P_T+273.15) + 0.0006 * LarvSize(n)* (P_T+273.15) + 0.0017) /1000.0  ! m/s  ! CHECK THAT FORMULATION IS IN mm/s
     P_size=LarvSize(n)
     targetlayerLOWERdepth=(P_zetac-abs(surflayer_lowerdepth))
!    targetlayerUPPERdepth=(P_zetac-abs(surflayer_upperdepth))

!    if(P_behave(n).lt.3)then ! Determine targetlayerLOWERdepth and targetlayerUPPERdepth
!              if(swdown_ASCII )then  !.or. readNetCdfSwdown
!                 swdown_r1=int((ix(2)-float(Ext0)-swdown_t0)/float(swdown_dt))
!                 swdown_r2=swdown_r1+1
!                 if(swdown_r2>swdown_numrec)&
!                      stop 'ERROR missing records in file swdown'
!                 swdown_alpha =&
!                    (float(swdown_r2)-(ix(2)- &
!                    float(Ext0)-swdown_t0)/float(swdown_dt)) &
!                    / float(swdown_r2-swdown_r1)
!                 swdown_interp=swdown_alpha*swdown(swdown_r1+1) &
!                        +(1.0-swdown_alpha)*swdown(swdown_r2+1)
!               if(swdown_interp>swdown_thresh)then
!                 targetlayerLOWERdepth=(P_zetac-abs(surflayer_lowerdepth))
!                 targetlayerUPPERdepth=(P_zetac-abs(surflayer_upperdepth))
!               else
!                 targetlayerLOWERdepth=&
!                                 (P_zetac-abs(surflayer_lowerdepth_night))
!                 targetlayerUPPERdepth=&
!                                 (P_zetac-abs(surflayer_upperdepth_night))
!               endif
!              else
!                  targetlayerLOWERdepth=(P_zetac-abs(surflayer_lowerdepth))
!                  targetlayerUPPERdepth=(P_zetac-abs(surflayer_upperdepth))
!              endif
!    endif

     IF(P_behave(n).eq.0)THEN !TYPE 0. First instant of life
         parBehav=(P_depth-P_zc+1.0)/float(idt)    ! Release at 1m from bottom
         P_behave(n)=2   

   ! ELSEIF(P_behave(n).eq.1)THEN !TYPE 1. Surface oriented. Particle swims up
   !     parBehav=ws_speed
   !     if(P_zc .GT. targetlayerLOWERdepth) then
   !        P_behave(n)=2  
   !     endif

      ELSEIF(P_behave(n).eq.2)THEN !TYPE 2. Follow chlorophyl gradient
         timer(n) =  timer(n)+DBLE(idt) 
         if((timer(n)>maxtimeatsurf) .or.    &   !IF((Behavior.EQ.9.and.timer(n)>30*86400) .or. &   
            (LarvSize(n)>maxsizeatsurf) )then   !   (Behavior.EQ.10.and.LarvSize(n)>8.0) )THEN    ! P_size>14.0 for nephrops
               parBehav=-abs(ws_speed)
               P_behave(n)=3
         else
!              if (P_zc .LT. targetlayerLOWERdepth) then
!               parBehav=rise  !-(P_zc+abs(surflayer_upperdepth))
!              elseif (P_zc .GT. targetlayerUPPERdepth)then
!               parBehav=-abs(sink)  !-(P_zc+abs(surflayer_upperdepth))
!              else
                   !swim up or down according to Chl gradient 
                   negpos = 1.0
                   if(abs(P_chl_gradient)>1e-8) then
                     if(P_chl_gradient<0) negpos = -1.0
                   !P_Weight= 4.0 * PI * LarvSize(n)**3 / 3.0
                   !L = 0.0541 * log(P_Weight) + 0.6154
                     parBehav= negpos * ws_speed ! m/s 
                   else
                     parBehav=0.0
                   endif
!              endif
         endif

       if(P_zc+parBehav*float(idt)<(P_depth+BottomLayerThickness))&
                 parBehav=(P_depth-P_zc+BottomLayerThickness+0.1)/float(idt)

      ELSEIF(P_behave(n).eq.3)THEN !TYPE 3. Sink to bottom 
         parBehav=-abs(ws_speed)
         if(P_zc .LT. (P_depth+1.0)) then
            write(*,'(i2)',advance='no')P_behave(n)
            P_behave(n) = 4
         endif

      ELSEIF(P_behave(n).eq.4)THEN !TYPE 4. Near-bottom 
         btest = 0   !switch to control behavior
         
         !particle has 80% change of swimming down if greater than "Seabed_layerheight" m 
         !  from bottom
         if ((Zgrid_depthinterp .and. P_zc .GT.P_depth+Seabed_layerheight).or.&
             ((.not.Zgrid_depthinterp) .and.                   &
             P_zc .GT. ( P_depth+1.35*sqrt(abs(P_depth))) ) )then
            negpos = 1.0
            dev1=genrand_real1()
            switch = 0.20 
            if (dev1.GT.switch) negpos = -1.0
            devB=genrand_real1()
            parBehav=negpos*devB*abs(ws_speed)
            btest = 1
         end if
         
         !if within "Seabed_layerheight" m of bottom, just swim randomly
         if (btest.EQ.0) then    
            negpos = 1.0
            dev1=genrand_real1()
            switch = 0.5 
            if (dev1.GT.switch) negpos = -1.0
            devB=genrand_real1()
            parBehav=negpos*devB*abs(ws_speed)   
         end if
       if(P_zc+parBehav*float(idt)>(P_depth+Seabed_layerheight))&
                 parBehav=(P_depth-P_zc+0.5)/float(idt)
      ENDIF
    ZBehav = parBehav * idt

  END SUBROUTINE Ostrea_Edulis_behave

  SUBROUTINE updateStatus(P_age,n)  !Update particle status
    USE SETTLEMENT_MOD, ONLY: isSettled
    USE PARAM_MOD, ONLY: settlementon,mortality
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION, INTENT(IN) :: P_age

    !Determine if particle dies from old age, if so kill it
    if ((P_age .GE. P_deadage(n)) .AND. mortality) then
      if(settlementon)then
        if(.NOT. isSettled(n)) call die(n)
      else
        call die(n)
      endif
    endif

  END SUBROUTINE updateStatus

  SUBROUTINE behave(Xpar,Ypar,Zpar,Pwc_zb,Pwc_zc,Pwc_zf,P_zb,P_zc,P_zf,        &
                    P_zetac,P_age,P_depth,P_U,P_V,P_angle,P_T,     &
                    n,it,ex,ix,          &
                    daytime,p,bott,XBehav,YBehav,ZBehav,P_size,Fstlev,klev)!,P_swdown)
    USE PARAM_MOD, ONLY: us,dt,idt,twistart,twiend,Em,pi,daylength,Kd,thresh,  &
                    Sgradient,swimfast,swimstart,sink,Hswimspeed,              &
                    Swimdepth,rise,surflayer_upperdepth,surflayer_lowerdepth,  &
                    surflayer_lowerdepth_night,surflayer_upperdepth_night,     &
                    Behavior,Zgrid_depthinterp,BottomLayerThickness,  &
                    swdown_dt,swdown_thresh,Seabed_layerheight, &
                    swdown_ASCII,swdown_t0,Ext0!,readNetCdfSwdown

    USE HYDRO_MOD, ONLY: WCTS_ITPI
    USE RANDOM_MOD, ONLY: genrand_real1
#include "VAR_IDs.h"
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: daytime
    DOUBLE PRECISION, INTENT(IN) :: Xpar,Ypar,Zpar,Pwc_zb(:),Pwc_zc(:),        &
                                    Pwc_zf(:),P_zb,P_zc,P_zf,P_zetac,P_age,    &
                                    P_depth,P_U,P_V,P_angle,ex(3),ix(3) !,       &
                                    !P_grainsize ,P_swdown
    INTEGER, INTENT(IN) :: n,it,p,Fstlev,klev
    LOGICAL, INTENT(OUT) :: bott
    DOUBLE PRECISION, INTENT(OUT) :: XBehav,YBehav,ZBehav,P_size
    
    INTEGER :: btest,i,deplvl
    INTEGER :: swdown_r1,swdown_r2                                  ! for swdown dvm
    DOUBLE PRECISION :: swdown_alpha,swdown_interp                  ! for swdown dvm
    DOUBLE PRECISION :: targetlayerUPPERdepth,targetlayerLOWERdepth ! for Behavior 8 and 9
    DOUBLE PRECISION :: negpos,dev1,devB,switch,switchslope
    DOUBLE PRECISION :: P_S,parBehav,Sslope,deltaS,deltaz
    DOUBLE PRECISION :: P_T !not needed unless temperature code below is enabled
    DOUBLE PRECISION :: dtime,tst,E0,P_light
    DOUBLE PRECISION :: currentspeed,Hdistance,theta,X,Y

    P_size=0.0
    !   ***************** Initialize Return Values
    XBehav = 0.0
    YBehav = 0.0
    ZBehav = 0.0
    P_S= -9999
        if(P_age .GE. swimstart) P_swim(n,3) = P_swim(n,1)*P_age+P_swim(n,2)
        if(P_age .GE. P_pediage(n)) P_swim(n,3) = swimfast

    !   ***************** Update vertical swimming speeds based on particle age
    IF(Behavior.lt.8)THEN
      if(Behavior.le.4)then
        P_swim(n,3) = swimfast
      else
        if(P_age .GE. swimstart) P_swim(n,3) = P_swim(n,1)*P_age+P_swim(n,2)
        if(P_age .GE. P_pediage(n)) P_swim(n,3) = swimfast
      endif
    
     !   ***************** Prepare for TYPE 4 & 5 (Oyster Larvae) Behaviors
      !Update pediveliger behavior/status and timer
      IF(P_behave(n) .EQ. 4 .OR. P_behave(n) .EQ. 5) THEN
         
        !Set behavior code for pediveligers
        if (P_age .GE. P_pediage(n) .AND. P_age .LT. P_deadage(n)) then
          P_behave(n) = 2
        endif
         
        !decrement timer
        timer(n) = max(DBLE(0.0), timer(n)-DBLE(dt))
      ENDIF
      
      !obtain salinity at particle location (P_S) to cue oyster larvae or tidal
      !  stream transport behavior
      IF ((P_behave(n).EQ.4) .OR. (P_behave(n).EQ.5 .AND. timer(n).EQ.0.0) .OR.  &
           P_behave(n).EQ.7) THEN 
         
        do i=3,us-2
          if ((Zpar .LT. Pwc_zb(i)) .OR. (Zpar .LT. Pwc_zc(i)) .OR.              &
              (Zpar .LT. Pwc_zf(i))) exit
        enddo
        deplvl = i-2   !depth level
         
        !Salinity at particle location
        P_S = WCTS_ITPI(VAR_ID_salt,Xpar,Ypar,deplvl,Pwc_zb,Pwc_zc,Pwc_zf,us,P_zb,    &
                        P_zc,P_zf,ex,ix,p,4,n)
         
      ENDIF
    ENDIF

    
    !           *********************************************************
    !           *                                                       *
    !           *                   Behaviors                           *
    !           *                                                       *
    !           *********************************************************
    
    parBehav = 0.0

    IF(Behavior.GE.8 .and. Behavior.LE.11) THEN
     IF(Behavior.EQ.10 .AND. P_behave(n).lt.2) THEN
         timer(n) =  timer(n)+DBLE(idt)
     ENDIF 
     ! GROWTH
     !--------------- TYPE 8. Nephrops-----------------------------------
     If(Behavior.EQ.8)then
       LarvSize(n)=LarvSize(n)+   &
                (0.02*max(DBLE(4.0),P_T)+0.04)*DBLE(idt)/DBLE(86400.0)
     !--------------- TYPE 9. Solea Solea -------------------------------
     ElseIf(Behavior.EQ.9)then
       if(LarvSize(n)<4.0)then
         LarvSize(n)=LarvSize(n)+  DBLE(idt)/DBLE(86400.0)* &
            max(0.,       0.1023*EXP(0.099*P_T)    )
!        LarvSize(n)=LarvSize(n)+  DBLE(idt)/DBLE(86400.0)* &
!               (0.0436*max(DBLE(8.0),P_T)-0.178)
       else
         LarvSize(n)=LarvSize(n)+  LarvSize(n)*DBLE(idt)/DBLE(86400.0)* &
            max(0.     ,0.0459*LOG(max(9.0,P_T))-0.0882      )
       endif
     Endif
     P_size=LarvSize(n)

     !--------------- TYPE 8. Nephrops, 9. Solea Solea and 10. Mullus Barbatus and 11. OTHERS ----------
     targetlayerLOWERdepth=(P_zetac-abs(surflayer_lowerdepth))
     targetlayerUPPERdepth=(P_zetac-abs(surflayer_upperdepth))
     if(P_behave(n).lt.3)then ! Determine targetlayerLOWERdepth and targetlayerUPPERdepth
               if(swdown_ASCII )then  !.or. readNetCdfSwdown
                !if(readNetCdfSwdown)then
                !  swdown_interp=P_swdown
                !else
                  swdown_r1=int((ix(2)-float(Ext0)-swdown_t0)/float(swdown_dt))
                  swdown_r2=swdown_r1+1
                  if(swdown_r2>swdown_numrec)&
                       stop 'ERROR missing records in file swdown'
                  swdown_alpha =&
                     (float(swdown_r2)-(ix(2)- &
                     float(Ext0)-swdown_t0)/float(swdown_dt)) &
                     / float(swdown_r2-swdown_r1)
                  swdown_interp=swdown_alpha*swdown(swdown_r1+1) &
                         +(1.0-swdown_alpha)*swdown(swdown_r2+1)
                !endif
                if(swdown_interp>swdown_thresh)then
                  targetlayerLOWERdepth=(P_zetac-abs(surflayer_lowerdepth))
                  targetlayerUPPERdepth=(P_zetac-abs(surflayer_upperdepth))
                else
                  targetlayerLOWERdepth=&
                                  (P_zetac-abs(surflayer_lowerdepth_night))
                  targetlayerUPPERdepth=&
                                  (P_zetac-abs(surflayer_upperdepth_night))
                endif
               else
                   targetlayerLOWERdepth=(P_zetac-abs(surflayer_lowerdepth))
                   targetlayerUPPERdepth=(P_zetac-abs(surflayer_upperdepth))
               endif
     endif

     if(P_behave(n).eq.0)then !TYPE 0. First instant of life
         parBehav=(P_depth-P_zc+1.0)/float(idt)    ! Release at 1m from bottom
         P_behave(n)=1   
         !write(*,*)'part ',n,' on bottom , surface oriented ',idt,timer(n),P_depth
     elseif(P_behave(n).eq.1)then !TYPE 1. Surface oriented. Particle swims up
         parBehav=rise
         if(P_zc .GT. targetlayerLOWERdepth) then
            P_behave(n)=2  
            !write(*,*)'part ',n,' reached surface target layer',timer(n)
         endif
     elseif(P_behave(n).eq.2)then !TYPE 2. Stay within the Surface Layer
         timer(n) =  timer(n)+DBLE(idt) 
         IF((timer(n)>maxtimeatsurf) .or.    &   !IF((Behavior.EQ.9.and.timer(n)>30*86400) .or. &   
            (LarvSize(n)>maxsizeatsurf) )THEN    !   (Behavior.EQ.10.and.LarvSize(n)>8.0) )THEN    ! P_size>14.0 for nephrops
               parBehav=-abs(sink)
               P_behave(n)=3
          !write(*,*)'part ',n,' sinking',idt,timer(n),LarvSize(n),Behavior
         ELSE
               if (P_zc .LT. targetlayerLOWERdepth) then
                parBehav=rise  !-(P_zc+abs(surflayer_upperdepth))
               elseif (P_zc .GT. targetlayerUPPERdepth)then
                parBehav=-abs(sink)  !-(P_zc+abs(surflayer_upperdepth))
               else
                  !just swim randomly (50% chance of swimming up)
                   negpos = 1.0
                   dev1=genrand_real1()
                   switch = 0.5 
                   if (dev1.GT.switch) negpos = -1.0
                   devB=genrand_real1()
                   parBehav=negpos*devB*0.5*(abs(sink)+abs(rise))   
               endif
         ENDIF
       if(P_zc+parBehav*float(idt)<(P_depth+BottomLayerThickness))&
                 parBehav=(P_depth-P_zc+BottomLayerThickness+0.1)/float(idt)
     elseif(P_behave(n).eq.3)then !TYPE 3. Sink to bottom 
         parBehav=-abs(sink)
         if(P_zc .LT. (P_depth+1.0)) then
            write(*,'(i2)',advance='no')P_behave(n)
            P_behave(n) = 4
            !write(*,'(2(a,i5),3(a,f8.2))') &
            !'it',it,'Part',n,' reached bottom arriving at P_zc=',P_zc,&
            !'< P_depth=',P_depth,'+1.0 =',(P_depth+1.0)
         endif
     elseif(P_behave(n).eq.4)then !TYPE 4. Near-bottom 
         btest = 0   !switch to control behavior
         
         !particle has 80% change of swimming down if greater than "Seabed_layerheight" m 
         !  from bottom
         if ((Zgrid_depthinterp .and. P_zc .GT.P_depth+Seabed_layerheight).or.&
             ((.not.Zgrid_depthinterp) .and.                   &
             P_zc .GT. ( P_depth+1.35*sqrt(abs(P_depth))) ) )then
            negpos = 1.0
            dev1=genrand_real1()
            switch = 0.20 
            if (dev1.GT.switch) negpos = -1.0
            devB=genrand_real1()
            parBehav=negpos*devB*abs(sink)
            btest = 1
         end if
         
         !if within "Seabed_layerheight" m of bottom, just swim randomly
         if (btest.EQ.0) then    
            negpos = 1.0
            dev1=genrand_real1()
            switch = 0.5 
            if (dev1.GT.switch) negpos = -1.0
            devB=genrand_real1()
            parBehav=negpos*devB*abs(sink)   
         end if
       if(P_zc+parBehav*float(idt)>(P_depth+Seabed_layerheight))&
                 parBehav=(P_depth-P_zc+0.5)/float(idt)
     endif

    ELSE

     !-----------Original LTRANS v.2b behaviors -------------------------------
     !TYPE 1. Surface oriented. Particle swims up if deeper than 1 m.
     IF (P_behave(n).EQ.1) THEN 
        btest = 0   !switch to control behavior
        
        !particle has 80% chance of swimming up if deeper than 1.0 m of bottom
        if (P_zc .LT. (P_zetac-1.0)) then
           negpos = 1.0
           dev1=genrand_real1()
           switch = 0.80 
           if (dev1.GT.switch) negpos = -1.0
           devB=genrand_real1()
           parBehav=negpos*devB*P_swim(n,3) 
           btest = 1
        end if
        
        !if within 1 m of surface, swim randomly (50% chance of swimming up)
        if (btest.EQ.0) then    
           negpos = 1.0
           dev1=genrand_real1()
           switch = 0.5 
           if (dev1.GT.switch) negpos = -1.0
           devB=genrand_real1()
           parBehav=negpos*devB*P_swim(n,3)  
        end if
        
     END IF
     
     
     !TYPE 2. Near-bottom. Particle swim down if not within 1 m of bottom.
     IF (P_behave(n).EQ.2 .OR. (P_behave(n).EQ.5 .AND. timer(n).GT.0.0)) THEN
        btest = 0   !switch to control behavior
        
        !particle has 80% change of swimming down if greater than 1.0 m 
        !  from bottom
        if (P_zc .GT. (P_depth+1.0)) then
           negpos = 1.0
           dev1=genrand_real1()
           switch = 0.20 
           if (dev1.GT.switch) negpos = -1.0
           devB=genrand_real1()
           parBehav=negpos*devB*P_swim(n,3)
           btest = 1
        end if
        
        !if within 1 m of bottom, just swim randomly
        if (btest.EQ.0) then    
           negpos = 1.0
           dev1=genrand_real1()
           switch = 0.5 
           if (dev1.GT.switch) negpos = -1.0
           devB=genrand_real1()
           parBehav=negpos*devB*P_swim(n,3)   
        end if
        
     END IF
     
     !TYPE 3: Diurnal Vertical Migration
     IF (P_behave(n).EQ.3) THEN
        
        !A. Find daytime in hrs since midnight (dtime)
        dtime = (daytime - aint(daytime))*DBLE(24.0)  !time of day 
        !This assumes that model simulations start at midnight
        
        !B. Calcluate irradiance at the water's surface (E0)
        tst = 0.0  !seconds since twilight start
        E0 = 0.0   !irradiance at the water's surface
        if (dtime.GT.twiStart .AND. dtime.LT.twiEnd) then
           tst=(dtime-twiStart)*DBLE(3600.0)
           E0= Em*SIN(PI*tst/(daylength*DBLE(3600.0)))*     &
                  SIN(PI*tst/(daylength*DBLE(3600.0)))
        else 
           E0 = 0.0
        end if
        
        !C. Calcluate irradiance at depth of the particle
        P_light = E0 * exp(Kd*P_zc)
        
        !If light at particle location is less than threshold, random swimming
        if (P_light.LT.thresh ) then
           negpos = 1.0
           dev1=genrand_real1()
           switch = 0.5 
           if (dev1.GT.switch) negpos = -1.0
           devB=genrand_real1()
           parBehav=negpos*devB*P_swim(n,3)   
        end if
        !If light at particle > threshold, then have 80% chance of swimming down
        if (P_light.GT.thresh ) then
           negpos = 1.0
           dev1=genrand_real1()
           switch = 0.20 
           if (dev1.GT.switch) negpos = -1.0
           devB=genrand_real1()
           parBehav=negpos*devB*P_swim(n,3)  
        end if
        
     END IF
     
     
     !TYPE 4. Crassostrea virginica -- above the halocline
     IF (P_behave(n).EQ.4) THEN 
        if (it.EQ.1) then
           P_Sprev(n) = P_S  !for first iteration
           P_zprev(n) = P_zc
        endif
        btest = 0   !switch to control behavior
        Sslope = 0.0  !salinity gradient that larvae swam through
        
        !determine if larva swam through salinity gradient large enough to 
        !  cue behavior;  if so, then 80% chance of swimming up                                            
        deltaS = P_Sprev(n) - P_S
        deltaz = P_zprev(n) - P_zc
        if (it.GT.1) Sslope = deltaS/deltaz
        if (abs(Sslope).GT.Sgradient) then
           negpos = 1.0
           dev1=genrand_real1()
           switch = 0.80 
           if (dev1.GT.switch) negpos = -1.0
           parBehav=negpos*P_swim(n,3)
           btest = 1
        endif
        
        !if no directed swimming, swim randomly with probabilities that result
        !  in particles moving up initially, then slowly moving toward bottom 
        !  with increasing age
        if (btest.EQ.0) then    
           negpos = 1.0
           dev1=genrand_real1()
           if (P_age .LT. 1.5*24.*3600.) then     !if Age < 1.5 Days
             switch = 0.1
           elseif (P_age .LT. 5.*24.*3600.) then  !if 1.5 Days <= Age < 5.0 Days
             switch = 0.49
           elseif (P_age .LT. 8.*24.*3600.) then  !if 5.0 Days <= Age < 8.0 Days
             switch = 0.50
           else                                   !if Age >= 8.0 Days
              switchslope = (DBLE(0.50)-DBLE(0.517)) /                          &
                            (DBLE(8.0)*DBLE(24.0)*DBLE(3600.0) - P_pediage(n))
              switch = switchslope*P_age + DBLE(0.50) -                         &
                       switchslope*DBLE(8.0)*DBLE(24.0)*DBLE(3600.0) 
              if (P_zc .LT. P_depth+1.) switch = 0.5
           endif
           if (dev1.GT.(1-switch)) negpos = -1.0
           devB=genrand_real1()
           parBehav=negpos*devB*P_swim(n,3)  
        endif
        
        !update previous salt and depth matrix for next iteration
        P_Sprev(n) = P_S
        P_zprev(n) = P_zc
     ENDIF
     
     
     !TYPE 5. Crassostrea ariakensis -- below the halocline
     IF (P_behave(n).EQ.5 .AND. timer(n).EQ.0.0) THEN 
        if (it.EQ.1) then
           P_Sprev(n) = P_S  !for first iteration
           P_zprev(n) = P_zc
        endif
        btest = 0   !switch to control behavior
        Sslope = 0.0  !salinity gradient that larvae swam through
        
        !determine if larva swam through salinity gradient large enough to 
        !  cue behavior.  If so, then 80% chance of swimming down. Set timer
        !  to keep particle near bottom for 2 hrs
        deltaS = P_Sprev(n) - P_S
        deltaz = P_zprev(n) - P_zc
        if (it.GT.1) Sslope = deltaS/deltaz
        if (abs(Sslope).GT.Sgradient) then
           negpos = 1.0
           dev1=genrand_real1()
           switch = 0.20 
           btest = 1
           timer(n) = DBLE(2.0)*DBLE(3600.0)  !2 hr times 3600 s
           if (dev1.GT.switch) negpos = -1.0
           parBehav=negpos*P_swim(n,3) 
           !keep bottom oriented behavior from starting until after particle
           !  is 3.5 days old  
           if (P_age .LT. 3.5*24.*3600.) then  
              btest = 0
              timer(n) = 0.
           endif
        endif
        
        !if no directed swimming, just swim randomly with probabilities that 
        !  result in particles moving up initially, then moving toward bottom 
        !  with increasing age
        if (btest.EQ.0) then    
           negpos = 1.0
           dev1=genrand_real1()
           switch = 0.495 
           if (P_age .LT. 1.5*24.*3600.) switch = 0.9
           if (P_age .GT. 2.0*24.*3600. .AND. P_age .LT. 3.5*24.*3600.) then
              switchslope = (DBLE(0.3)-DBLE(0.495)) /                 &
                            (DBLE(2.0)*DBLE(24.0)*DBLE(3600.0) -      &
                             DBLE(3.5)*DBLE(24.0)*DBLE(3600.0))
              switch = switchslope*P_age+DBLE(0.3) -                  &
                       switchslope*DBLE(2.0)*DBLE(24.0)*DBLE(3600.0) 
           endif
           if (dev1.GT.switch) negpos = -1.0
           devB=genrand_real1()
           parBehav=negpos*devB*P_swim(n,3)  
        endif
        
        !update previous salt and depth matrix for next iteration        
        P_Sprev(n) = P_S
        P_zprev(n) = P_zc
     ENDIF
     
     !TYPE 6. Constant -- no random motion to vertical movement
     IF ((P_behave(n).EQ.6)) THEN               
        if(P_age .GE. swimstart) then
           parBehav = sink 
        else
           parBehav = P_swim(n,3)
        endif
        
        !Note: the code below is included if someone wants to calculate density
        ! ! To calculate salinity (P_S) and temperature (P_T) at particle location
        !  do i=3,us-2
        !    if ((Zpar .LT. Pwc_zb(i)) .OR. (Zpar .LT. Pwc_zc(i)) .OR.     & 
        !        (Zpar .LT. Pwc_zf(i))) exit
        !  enddo
        !  deplvl = i-2   !depth level
        !
        !  !Salinity at particle location
        !  P_S = WCTS_ITPI("salt",Xpar,Ypar,deplvl,Pwc_zb,Pwc_zc,Pwc_zf,us,     &
        !                  P_zb,P_zc,P_zf,ex,ix,p,4)
        ! 
        !  !Temperature at particle location 
        !  P_T = WCTS_ITPI("temp",Xpar,Ypar,deplvl,Pwc_zb,Pwc_zc,Pwc_zf,us,     &
        !                  P_zb,P_zc,P_zf,ex,ix,p,4)
     ENDIF
    ENDIF !-----------(Nephrops or Original v.2b behaviors)--------------------

    !Calculate movement due to behavior for all behaviors other than 7
    ZBehav = parBehav * idt
    
    !TYPE 7. Tidal Stream Transport: if flooding, then swim in direction of
    !  currents, else sit on bottom
    IF ((P_behave(n).EQ.7)) THEN 
       ! Set initial values for the first iteration
       if (it.EQ.1) then         
          P_Sprev(n) = P_S
          currentspeed = 0.0    
       endif
       !Find current speed at the particle location ( c = sqrt(a**2 + b**2) )
       currentspeed = sqrt( (P_U*cos(P_angle) - P_V*sin(P_angle))**2 +         &
                            (P_U*sin(P_angle) + P_V*cos(P_angle))**2 )
       if (bottom(n) .EQV. .TRUE.) then          !CRS
          !if particle is on bottom, test if salinity is increasing
          if (P_Sprev(n).LT.P_S) then       !if salinity is increasing:
             bottom(n) = .FALSE.            !  come off bottom 
             ZBehav = P_depth + Swimdepth   !  and swim to the swimming depth
          else
             ZBehav = -9999  !if salinity is not increasing, stay on bottom
          end if
       else        
          !if particle is not on bottom, test if currents are not slack 
          !  (defined as 0.05 m/s)
          if (currentspeed.GT.0.05) then
             !if the current speed is greater than 0.05 m/s, then swim in the 
             !  direction of the current
             Hdistance = Hswimspeed*idt
             !find theta of currents
             theta = atan( (P_U*sin(P_angle) + P_V*cos(P_angle)) /             &
                           (P_U*cos(P_angle) - P_V*sin(P_angle)) )
             X = (P_U*cos(P_angle) - P_V*sin(P_angle))
             Y = (P_U*sin(P_angle) + P_V*cos(P_angle))
             if(X.GT.0.0) then
                XBehav =  Hdistance*cos(theta)     
                YBehav =  Hdistance*sin(theta)     
             end if
             if(X.LT.0.0) then
                XBehav =  DBLE(-1.0)*Hdistance*cos(theta)     
                YBehav =  DBLE(-1.0)*Hdistance*sin(theta)     
             end if
             if(X.EQ.0 .AND. Y.GE.0.0) then
                XBehav =  0.0     
                YBehav =  Hdistance
             end if
             if(X.EQ.0 .AND. Y.LE.0.0) then
                XBehav =  0.0     
                YBehav =  DBLE(-1.0)*Hdistance
             end if
             !keep vertical position of particle at swim depth
             ZBehav = P_depth + Swimdepth
          else
             ZBehav = -9999       !if the current speed is less than 0.05 m/s,
             bottom(n) = .TRUE.   !  then swim to bottom
          end if
       end if
       bott = bottom(n)
    ENDIF

    IF(Behavior.eq.12)THEN
      call Ostrea_Edulis_behave(Xpar,Ypar,Zpar,Pwc_zb,Pwc_zc,Pwc_zf,P_zb,P_zc,P_zf,        &
                    P_zetac,P_age,P_depth,P_U,P_V,P_angle,P_T,     &
                    n,it,ex,ix,          &
                    daytime,p,bott,XBehav,YBehav,ZBehav,P_size,Fstlev,klev)
    ENDIF

! ******************* End Particle Behavior ******************************
  END SUBROUTINE behave


  INTEGER FUNCTION getStatus(n,default_status)
    !This function returns an identification number that describes a particle's  
    !behavior type or status for use in visualization routines. It was
    !initially developed to contain the color code for plotting in Surfer.)                
    USE PARAM_MOD, ONLY: SETTLEMENTON,OPENOCEANBOUNDARY,OilOn,stranding_on
    USE SETTLEMENT_MOD, ONLY: isSettled
    USE STRANDING_MOD, ONLY: isStranded
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION, INTENT(IN) :: default_status
    if(OilOn)then
    getStatus = default_status
    else
    getStatus = P_behave(n)          ! Set Status to behavior ID
    endif
                                     ! Change if Dead, Settled, or OutOfBounds


    if(dead(n)) getStatus = -1         ! -1 = Dead
    if(settlementon)then
      if(isSettled(n)) getStatus = -2  ! -2 = Settled
    endif
    if(stranding_on)then
      if(isStranded(n)) getStatus = -4  ! -2 = Stranded
    endif
    if(OpenOceanBoundary)then
      if(oob(n)) getStatus = -3        ! -3 = Out of Bounds
    endif

  END FUNCTION getStatus


  LOGICAL FUNCTION isDead(n)
  !This function returns .TRUE. if the particle is "dead", and FALSE if not
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n

    isDead = dead(n)

  END FUNCTION isDead


  SUBROUTINE die(n)
  !This subroutine sets the value of dead(n) to TRUE, indicating 
  !  the particle is "dead"
    USE HYDRO_MOD, ONLY: setDeadOrOut
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n

    dead(n) = .TRUE.
    call setDeadOrOut(n)

  END SUBROUTINE die


  SUBROUTINE setOut(n)
    !This subroutine changes particle n's status to out of bounds
    USE HYDRO_MOD, ONLY: setDeadOrOut
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n

    oob(n) = .TRUE.
    call setDeadOrOut(n)

  END SUBROUTINE setOut

  LOGICAL FUNCTION isOut(n)
    !This function returns the value of oob for particle n
    ! i.e. Returns True if the particle is out of bounds, False if in bounds
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n

    isOut = oob(n)

  END FUNCTION isOut


  SUBROUTINE finBehave()    !Finish the behavior module
  USE PARAM_MOD     , ONLY: settlementon,swdown_ASCII
  USE SETTLEMENT_MOD, ONLY: finSettlement
  IMPLICIT NONE

    !Deallocate Behavior Variables
    DEALLOCATE(timer)
    DEALLOCATE(P_behave)
    DEALLOCATE(P_pediage)
    DEALLOCATE(P_deadage)
    DEALLOCATE(P_Sprev)
    DEALLOCATE(P_zprev)
    DEALLOCATE(P_swim)
    DEALLOCATE(dead)
    DEALLOCATE(oob)

    if(ALLOCATED(bottom))DEALLOCATE(bottom)
    if(ALLOCATED(LarvSize))DEALLOCATE(LarvSize)
    IF(swdown_ASCII)DEALLOCATE(swdown)

    !If Settlement is on, Deallocate Settlement Variables
    if(settlementon) call finSettlement()

  END SUBROUTINE finBehave


END MODULE BEHAVIOR_MOD
