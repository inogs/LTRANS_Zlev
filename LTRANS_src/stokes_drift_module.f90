  !-----------------------------------------------------------------------------
  ! * Stokes drift and F_WindAngle routines taken from                         *
  ! * the Irish Marine Institute Oil Model Module developed for LTRANS         *
  ! * by Alan Berry, July 2011                                                 *
  !-----------------------------------------------------------------------------
MODULE STOKES_DRIFT_MOD
  IMPLICIT NONE
  PUBLIC

  CONTAINS

  !***************************************************
  !*     Subroutine StokesDrift_Estimate_from_Wind   *
  !***************************************************
  SUBROUTINE StokesDrift_Estimate_from_Wind(uwind,vwind,pdir,depth,ustoke,vstoke,StokesDriftFac,WVecAngle)
    ! Give the Stokes drift at depth of a particle assuming a Pierson Moskowitz spectrum and Langmuir circulation
    ! Author Marcel Cure (www.numericswarehouse.com) Apr. 2010
    !
    DOUBLE PRECISION, INTENT(IN) :: uwind,vwind,pdir,depth
    DOUBLE PRECISION, INTENT(OUT) :: ustoke,vstoke                       ! velocity components of particle at depth
    DOUBLE PRECISION, INTENT(IN) :: StokesDriftFac
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
    knum = 1.25*(grav/(max(V10,1e-5)**2.0))
    VS0 =StokesDriftFac*V10
    VS = VS0*EXP(2.0*knum*depth)
    ustoke = VS*COS(phi)
    vstoke = VS*SIN(phi)
    !
  END SUBROUTINE StokesDrift_Estimate_from_Wind

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

END MODULE STOKES_DRIFT_MOD

