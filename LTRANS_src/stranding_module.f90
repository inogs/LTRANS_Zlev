MODULE STRANDING_MOD
  IMPLICIT NONE
  SAVE
  PRIVATE

  LOGICAL, ALLOCATABLE, DIMENSION(:) :: stranded
  PUBLIC :: initStranding,testStranding,finStranding,isStranded, &
             p_Stranding,testRefloating
CONTAINS

  SUBROUTINE initStranding()
    USE PARAM_MOD, ONLY: numpar
    ALLOCATE(stranded(numpar))
    stranded = .FALSE.
  END SUBROUTINE initStranding

  SUBROUTINE finStranding()
    IMPLICIT NONE

    DEALLOCATE(stranded)

  END SUBROUTINE finStranding

  LOGICAL FUNCTION isStranded(n)
  !This function returns .TRUE. if the particle is "stranded", and FALSE if not
    USE PARAM_MOD, ONLY:stranding_on
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n

    if(stranding_on)then
     isStranded = stranded(n)
    else
     isStranded = .FALSE.
    endif

  END FUNCTION isStranded

  SUBROUTINE testStranding(n,Px,Py,Pz, &
                            P_depth,k,intersectf,coastdist)
    USE PARAM_MOD, ONLY: StrandingDist,StrandingMaxDistFromSurf,    &
                         StrandingMaxDistFromBott,OilOn,stranding_on
    USE HYDRO_MOD, ONLY: getP_r_element,getP_klev
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n,k
    DOUBLE PRECISION, INTENT(IN) :: Px,Py,Pz,P_depth,coastdist
    INTEGER, INTENT(IN) :: intersectf                                    ! CL-OGS
    INTEGER :: klev

    INTEGER :: polyin,R_ele

    R_ele = getP_r_element(n)
    klev= getP_klev(n)

    if(stranding_on)then
    if((P_depth+abs(StrandingMaxDistFromBott)>=Pz) &
                 .and.( Pz>=-abs(StrandingMaxDistFromSurf) )          ) then 
      if(coastdist<StrandingDist .or.  intersectf>0)then! If (StrandingDist>0) n can strand
            stranded(n)=.TRUE.
      elseif(abs(StrandingDist)<1e-5 .and. intersectf>0)then! If (StrandingDist=0 and bound intersected) n can strand
            stranded(n)=.TRUE.
      endif
    endif
    endif
  END SUBROUTINE testStranding

  SUBROUTINE p_Stranding(n)
    USE PARAM_MOD, ONLY:stranding_on
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n

    if(stranding_on)then
     stranded(n) = .TRUE.
    else
     stranded(n) = .FALSE.
    endif

  END SUBROUTINE p_Stranding


  SUBROUTINE testRefloating(n,tS)
    ! With`refloat=.True.` the stranded particles may return into the water and be further advected.
    ! The probability P that a given particle has to refloat decreases exponentially with the time tS that this particle spent stranded:
    !   P_refloat = 1.0 - 0.5 exp(−tS/refloat_Th)
    ! Where Th is the half-life time.
    ! At each time step, for each stranded particle a random number generator, Rrefloat, is called up and the particle is released
    ! back into the water if Rrefloat < Prefloat
    USE PARAM_MOD, ONLY:stranding_on,refloat,refloat_Pc,refloat_Po,refloat_Tc
    USE RANDOM_MOD, ONLY: genrand_real1
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION, INTENT(IN) :: tS ! time spent stranded
    DOUBLE PRECISION :: P_refloat, Rand_refloat

    if(stranding_on .and. refloat)then
      P_refloat = refloat_Pc + refloat_Po * exp(-tS/refloat_Tc)
      Rand_refloat=genrand_real1()
      if(Rand_refloat<P_refloat)then
        stranded(n) = .FALSE.
      endif
    endif

  END SUBROUTINE testRefloating


END MODULE STRANDING_MOD
