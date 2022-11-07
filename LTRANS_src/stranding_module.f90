MODULE STRANDING_MOD
  IMPLICIT NONE
  SAVE
  PRIVATE

  LOGICAL, ALLOCATABLE, DIMENSION(:) :: stranded
  PUBLIC :: initStranding,testStranding,finStranding,isStranded, &
             p_Stranding
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
      if(coastdist<StrandingDist)then! If (StrandingDist>0) n can strand
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



END MODULE STRANDING_MOD
