MODULE HTURB_MOD

!   A random walk model is used to simulate turbulent particle motion in the 
!   horizontal direction (x- and y- directions).
!
! Created by:             Elizabeth North
! Modified by:            Zachary Schlag
! Created on:             2003
! Last Modified on:       18 Aug 2008

! ********************************************************************
! *                                                                  *
! *             Irish Marine Institute Oil Model Module              *
! *                         										 *
! * Langmuir circulation code										 *
! *                                                                  *
! * Created by: 	Marcel Cure (www.numericswarehouse.com)          *
! * Created on:		April 2010                                       *
! * Modified by:    Alan Berry (Marine Institute)                    *
! * Modified on:    01 July 2011                                     *
! *                                                                  *
! ********************************************************************
!--- CL-OGS: Modified on:    June 2016 by Celia Laurent - OGS
!--- CL-OGS: coupling of LTRANS v2b to the Z-grid and fields of the MITgcm
!--- CL-OGS: implementation of the Irish Marine Institute Oil Model Module in LTRANS v2bZ-MITgcm
IMPLICIT NONE
PUBLIC

CONTAINS

!       ***************** Horizontal Turbulence (RWM) *********************
!     ***********************************************************************
!    ** Random Walk Model (Visser 1997 MEPS 158:275-281) for simulating     **
!    ** turbulent diffusion, applied in the horizontal direction            **
!    **    z(t+1)= z + R[2/r K(z)dt]**0.5                                   **
!    **    where z = particle vertical location at time t                   **
!    **      K  = horizontal diffusivity (KM from QUODDY)                   **
!    **      dt = time step of RWM (deltat)                                 **
!    **      R  = random process with mean = 0 and standard deviation = r.  **
!    **                                                                     **
!    ** Programmed by EW North February 2003 UMCES HPL enorth@hpl.umces.edu **
!    *************************************************************************

  SUBROUTINE HTurb(TurbHx,TurbHy, &
!	*****   IMIOM   *****
                   uwind, vwind, pdir,WVecAngle)
!	***** END IMIOM *****

    USE PARAM_MOD, ONLY: ConstantHTurb,idt, &
!	*****   IMIOM   *****
                         Langmuir,Cd,OilOn,pi
    USE STOKES_DRIFT_MOD, ONLY: F_WindAngle
!	***** END IMIOM *****
    USE NORM_MOD,  ONLY: norm
    IMPLICIT NONE

!	*****   IMIOM   *****
    DOUBLE PRECISION, OPTIONAL, INTENT(IN) :: uwind,vwind,pdir,WVecAngle
!	***** END IMIOM *****
    DOUBLE PRECISION, INTENT(OUT) :: TurbHx,TurbHy

    DOUBLE PRECISION :: devX,devY,r
    DOUBLE PRECISION :: KM
!	*****   IMIOM   *****
    DOUBLE PRECISION :: Theta,WaveAng,WindWaveAng,V10,ustar,ustoke,VCL,dv,KH_Dw,KH_Xw,Turb_Dw,Turb_Xw
    DOUBLE PRECISION, PARAMETER :: RhoAir = 1.23
    DOUBLE PRECISION, PARAMETER :: RefRho = 1000.0
!	***** END IMIOM *****


    r=1.0             ! standard deviation of the random deviate

!   *****   IMIOM   *****
    IF(OilOn .AND. Langmuir)THEN

      Theta = ATAN2(vwind,uwind)
      if(present(WVecAngle))then
        WaveAng=WVecAngle
      else
        WaveAng = F_WindAngle(pdir)        ! convert wave direction to cartesian (radians)
      endif
      WindWaveAng = Theta - WaveAng    ! angle between wind and wave directions

      IF (WindWaveAng .GT. pi) WindWaveAng = (2.0 * pi) - WindWaveAng

      V10 = SQRT(uwind**2.0 + vwind**2.0)
      ustar = SQRT(Cd * RhoAir / RefRho ) * V10   ! Assume a neutral atmospheric boundary layer

      ! Use Pierson Moskowitz spectrum to calculate the Stokes drift for a fully developed sea at this wind speed.
      ! See Thorpe, S.A. 2004 'Langmuir Circulation', Ann. Rev. Fluid Mech. 36 55-79 for discussion
      ustoke = 0.014 * V10

      ! We use the component of the stokes drift in the wind direction.
      ! If wind and waves oppose, then set to zero
      IF (ABS(WindWaveAng) .LT. pi/2.0) THEN
        ustoke = ustoke * COS(WindWaveAng)
      ELSE
        ustoke = 0.0
      END IF

      ! The Leibovich Langmuir strength scale is VCL
      VCL = SQRT(ustar * ustoke)

      ! Use the sonar derived measurements of Plueddemann et. al. 1996 'Structure and variability of Langmuir Circulation
      ! during the Surface Waves Processes Program', J. Geograph. Res. C 101 3525-3543
      ! dv is the LC inflow speed
      dv = 2.0 * VCL   ! based on a fit to Fig. 12 of Plueddemann et. al.

      ! Li, M. 2000 'Estimating horizontal dispersion of floating particles in wind-driven upper ocean', Spill Science and Tech. Bull. 6 255-261
      KH_Dw = 20.0 * (dv**2.0)          ! based on a fit to Fig. 7.
      KH_Xw = 315.0 * (dv**(5.0/2.0))

      devX = norm()  ! the random deviate in the cross wind direction
      devY = norm()  ! the random deviate in the down wind  direction

      ! Apply random walk model to calculate horizontal turbulent particle displacement
      Turb_Xw = devX*(DBLE(2.0)/r * KH_Xw *idt)**0.5
      Turb_Dw = devY*(DBLE(2.0)/r * KH_Dw *idt)**0.5

      ! Reorient to X and Y directions
      TurbHx = Turb_Xw*COS(theta) - Turb_Dw*SIN(theta)
      TurbHy = Turb_Xw*SIN(theta) + Turb_Dw*COS(theta)

    ELSE
!   ***** END IMIOM *****

      KM=ConstantHTurb  ! constant horizontal diffusivity

      devX=norm()       ! the random deviate in the X direction
      devY=norm()       ! the random deviate in the Y direction

      !Apply random walk model to calculate horizontal turbulent 
      !  particle displacement
      TurbHx= devX*(DBLE(2.0)/r * KM *idt)**0.5 
      TurbHy= devY*(DBLE(2.0)/r * KM *idt)**0.5

!   *****   IMIOM   *****
    END IF
!   ***** END IMIOM *****

  END SUBROUTINE HTurb

END MODULE HTURB_MOD
