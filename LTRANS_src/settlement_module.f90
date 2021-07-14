MODULE SETTLEMENT_MOD

!  The Settlement Module handles all code related to the settlement routine.  
!  This includes reading in the habitat polygons and holes, creating variables
!  containing the specifications of the habitat polygons and holes, keeping
!  track of the settlement status of every particle, and checking if the 
!  particle is within a habitat polygon and can settle.
!
!  Original concepts and code by:             Elizabeth North
!  Module creation and code modification by:  Zachary Schlag
!  Created on:                                2005
!  Last modified on:                          8 Mar 2011

  IMPLICIT NONE
  SAVE
  PRIVATE

  TYPE :: polyPerEle
    
    !number of polygons in each element OR number of holes in each polygon
    INTEGER :: numpoly

    !id of each polygon in the element OR id of each hole in the polygon
    INTEGER, ALLOCATABLE, DIMENSION(:) :: poly

  END TYPE polyPerEle


  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: polys,holes
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:)   :: maxbdis,maxhdis,settletime
  INTEGER, ALLOCATABLE, DIMENSION(:,:)          :: polyspecs,holespecs
  INTEGER, ALLOCATABLE, DIMENSION(:)          :: polygonidlist
  LOGICAL, ALLOCATABLE, DIMENSION(:)            :: settle
  TYPE (polyPerEle), ALLOCATABLE, DIMENSION(:,:):: elepolys
  TYPE (polyPerEle), ALLOCATABLE, DIMENSION(:):: polyholes 
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: stranded
  ! vertical level numbers for MITgcm hydro using varying water cell at each level
  INTEGER :: us_tridim,npoly

  !The following procedures have been made public:
  PUBLIC :: initSettlement,testSettlement,isSettled,finSettlement,isStranded, &
             get_NumPoly,get_IDPoly,p_Stranding


CONTAINS

  SUBROUTINE initSettlement(P_pediage)
    USE PARAM_MOD, ONLY: numpar,max_rho_elements,minholeid,maxholeid,minpolyid,    &
                         maxpolyid,pedges,hedges,us,Zgrid
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: P_pediage(:)
    INTEGER :: n

    !ALLOCATE VARIABLES

    ALLOCATE(polys(pedges,5))
      ! polys(i,1)=habitat polygon number
      ! polys(i,2)=center x
      ! polys(i,3)=center y
      ! polys(i,4)=edge x
      ! polys(i,5)=edge y

    ALLOCATE(holes(hedges,6))
      ! Holes(i,1)=hole number
      ! Holes(i,2)=center longitude
      ! Holes(i,3)=center latitude
      ! Holes(i,4)=edge longitude
      ! Holes(i,5)=edge latitude
      ! Holes(i,6)=habitat polygon number

      !Maximum distance from each habitat polygon's center to its 
      !  farthest edge point
    ALLOCATE(maxbdis(minpolyid:maxpolyid))

      !Maximum distance from each hole's center to its farthest edge point
    ALLOCATE(maxhdis(minholeid:maxholeid))

      !FALSE=not settled, TRUE=successful settlement
    ALLOCATE(settle(numpar))
    ALLOCATE(stranded(numpar))

      !age particles are competent to settle
    ALLOCATE(settletime(numpar))

      !polyspecs(n,1) = location in polys of first point of habitat polygon n
      !polyspecs(n,2) = # of edge points that make up habitat polygon n
    ALLOCATE(polyspecs(minpolyid:maxpolyid,2))

    if(Zgrid)then
      us_tridim=us
    else
      us_tridim=1 
    endif
      !id # of every polygon in each element
    ALLOCATE(elepolys(us_tridim,max_rho_elements))

      !holespecs(n,1) = location in holes of first point of hole n
      !holespecs(n,2) = # of edge points that make up hole n
    ALLOCATE(holespecs(minholeid:maxholeid,2))

      !id # of every hole in each polygon
    ALLOCATE(polyholes(minpolyid:maxpolyid))


    settle = .FALSE.  !initialize to FALSE=not settled
    stranded = .FALSE.

    do n=1,numpar
      settletime(n) = P_pediage(n)
    enddo

    call getHabitat()      ! Read In Habitat Polygon Information

    call createPolySpecs() ! Organize Habitat data to speed up
                           !   testing for particle settlement

  END SUBROUTINE initSettlement

  SUBROUTINE finSettlement()
    IMPLICIT NONE

    !DEALLOCATE VARIABLES
    DEALLOCATE(settletime)
    DEALLOCATE(polyspecs)
    DEALLOCATE(holespecs)
    DEALLOCATE(polyholes)
    DEALLOCATE(elepolys)
    DEALLOCATE(maxbdis)
    DEALLOCATE(maxhdis)
    DEALLOCATE(settle)
    DEALLOCATE(stranded)
    DEALLOCATE(polys)
    DEALLOCATE(holes)

  END SUBROUTINE finSettlement


  SUBROUTINE getHabitat()
    USE PARAM_MOD, ONLY: pedges,hedges,habitatfile,holefile,holesexist
    USE CONVERT_MOD, ONLY: lon2x,lat2y
    IMPLICIT NONE

    INTEGER :: i,curpoly,curpoly2
    DOUBLE PRECISION :: dise
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: P_lonlat,H_lonlat

    ALLOCATE(P_lonlat(pedges,5))
    ALLOCATE(H_lonlat(hedges,6))

    ! ***********************************************************************
    ! *                     Initialize settlement model                     *
    ! ***********************************************************************

    write(*,*) 'read in habitat polygon locations'

    ! need to import center and edge coordinate data
    ! User specified file name
    OPEN(1,FILE=TRIM(habitatfile))

      do i=1,pedges
        read (1,*) curpoly,P_lonlat(i,2),P_lonlat(i,3),    &
                           P_lonlat(i,4),P_lonlat(i,5)
        P_lonlat(i,1) = DBLE(curpoly)
      enddo
      ! P_lonlat(i,1)=habitat polygon number,
      ! P_lonlat(i,2)=center longitude, P_lonlat(i,3)=center latitude,
      ! P_lonlat(i,4)=edge longitude,   P_lonlat(i,5)=edge latitude
    CLOSE(1)

    write(*,*) '  Edge i=5 Polygon ID=',P_lonlat(5,1)
    write(*,*) '  Edge i=5 Center Lat=',P_lonlat(5,3),'Long=',P_lonlat(5,2)
    write(*,*) '  Edge i=5   Edge Lat=',P_lonlat(5,5),'Long=',P_lonlat(5,4)


    ! Convert particle latitude and longitude to meters using equations from  
    ! sg_mercator.m and seagrid2roms.m in Seagrid. 
    ! polys(i,1)=habitat polygon number
    ! polys(i,2)=center x
    ! polys(i,3)=center y
    ! polys(i,4)=edge x
    ! polys(i,5)=edge y
    do i=1,pedges
      polys(i,1) = P_lonlat(i,1)

      polys(i,2) = lon2x(P_lonlat(i,2),P_lonlat(i,3))
      polys(i,3) = lat2y(P_lonlat(i,3))

      polys(i,4) = lon2x(P_lonlat(i,4),P_lonlat(i,5))
      polys(i,5) = lat2y(P_lonlat(i,5))
    enddo

    ! Determine maximum distance between center and edge coordinates for each
    !   habitat polygon; this is used to restrict settlement model search
    maxbdis = 1.0
    do i=1,pedges
      curpoly = nint(polys(i,1))
      dise=sqrt((polys(i,2)-polys(i,4))**2+(polys(i,3)-polys(i,5))**2)
      if (dise.GT.maxbdis(curpoly)) maxbdis(curpoly)=dise
    enddo


    if(holesExist)then

      ! User specified file name
      OPEN(2,FILE=TRIM(holefile))

        do i=1,hedges
          read (2,*) curpoly,H_lonlat(i,2),H_lonlat(i,3),H_lonlat(i,4),   &
                     H_lonlat(i,5),curpoly2
          H_lonlat(i,1) = curpoly
          H_lonlat(i,6) = curpoly2
        enddo
        ! H_lonlat(i,1)=hole number      | H_lonlat(i,4)=edge longitude
        ! H_lonlat(i,2)=center longitude | H_lonlat(i,5)=edge latitude
        ! H_lonlat(i,3)=center latitude  | H_lonlat(i,6)=habitat polygon number

      CLOSE(2)

      write(*,*) '  Hole i=5 Center Lat=',H_lonlat(1,3),'Long=',H_lonlat(1,2)
      write(*,*) '  Hole i=5   Edge Lat=',H_lonlat(1,5),'Long=',H_lonlat(1,4)

      ! Convert particle latitude and longitude to meters using equations from  
      ! sg_mercator.m and seagrid2roms.m in Seagrid. 
      ! Holes(i,1)=hole number
      ! Holes(i,2)=center x
      ! Holes(i,3)=center y
      ! Holes(i,4)=edge x
      ! Holes(i,5)=edge y
      ! Holes(i,6)=habitat polygon number
      do i=1,hedges
        holes(i,1) = H_lonlat(i,1)

        holes(i,2) = lon2x(H_lonlat(i,2),H_lonlat(i,3))
        holes(i,3) = lat2y(H_lonlat(i,3))

        holes(i,4) = lon2x(H_lonlat(i,4),H_lonlat(i,5))
        holes(i,5) = lat2y(H_lonlat(i,5))

        holes(i,6) = H_lonlat(i,6)
      enddo 


      ! Determine maximum distance between center and edge coordinates for
      !   each hole; this is used to restrict settlement model hole search
      maxhdis = 1.0
      do i=1,hedges
        curpoly = nint(holes(i,1))
        dise=sqrt((holes(i,2)-holes(i,4))**2+(holes(i,3)-holes(i,5))**2)
        if (dise.GT.maxhdis(curpoly)) maxhdis(curpoly)=dise
      enddo

    endif

    DEALLOCATE(P_lonlat)
    DEALLOCATE(H_lonlat)

  END SUBROUTINE getHabitat

  SUBROUTINE identify_polygons()
    USE PARAM_MOD, ONLY: pedges
    IMPLICIT NONE

    INTEGER :: i,countpoly 

    !store list of polynom identification numbers
    countpoly=1
    do i=2,pedges
      if(polys(i,1).ne.polys(i-1,1)) countpoly=countpoly+1
    enddo
    npoly=countpoly
    allocate(polygonidlist(npoly))
    polygonidlist(1)=int(polys(1,1))
    countpoly=1
    do i=2,pedges
      if(polys(i,1).ne.polys(i-1,1)) then
       countpoly=countpoly+1
       polygonidlist(countpoly)=int(polys(i,1))
      endif
    enddo

  END SUBROUTINE identify_polygons



  SUBROUTINE  get_NumPoly(poly0,polyN,npoly_out)
    IMPLICIT NONE
    integer, intent(out) :: poly0,polyN,npoly_out
    call identify_polygons()
    npoly_out=npoly
    poly0=minval(polygonidlist)
    polyN=maxval(polygonidlist)

  END SUBROUTINE  get_NumPoly


  SUBROUTINE get_IDPoly(polyidlist)
    IMPLICIT NONE

    INTEGER, INTENT(OUT) :: polyidlist(npoly)
    
    polyidlist=polygonidlist

  END SUBROUTINE get_IDPoly

 

  SUBROUTINE createPolySpecs()
    USE PARAM_MOD,    ONLY: max_rho_elements,holesExist,minpolyid,maxpolyid,  &
                            pedges,hedges
    USE HYDRO_MOD,    ONLY: getR_ele
    USE GRIDCELL_MOD, ONLY: gridcell
    USE PIP_MOD,      ONLY: inpoly
    IMPLICIT NONE


    TYPE :: polynum
      INTEGER :: num
      TYPE (polynum), POINTER :: p => NULL()
    ENDTYPE

    TYPE (polynum), POINTER :: polyhead => NULL(),polytail => NULL(),ptr => NULL()
    INTEGER :: ISTAT

    LOGICAL :: check
    INTEGER :: i,j,k,count,triangle,checkele,P_ele,klev
    DOUBLE PRECISION :: dis
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: poly,r_ele_x,r_ele_y

    ALLOCATE(r_ele_x(4,max_rho_elements))
    ALLOCATE(r_ele_y(4,max_rho_elements))


    write(*,*) 'find polygons in elements'
   !------------------
   DO klev=1,us_tridim
    CALL getR_ele(klev,r_ele_x,r_ele_y) 
    polyspecs = 0

    !iterate through each habitat polygon and store information:
    !  polyspecs(n,1) = location in polys of first point of habitat polygon n
    !  polyspecs(n,2) = # of edge points that make up habitat polygon n
    count = 1
    polyspecs(int(polys(1,1)),1) = 1
    do i=2,pedges
      if(polys(i,1) /= polys(i-1,1))then
        polyspecs(nint(polys(i,1)),1) = i
        polyspecs(nint(polys(i-1,1)),2) = count
        count = 1
      elseif(i==pedges)then
        polyspecs(nint(polys(i,1)),2) = count + 1
      else
        count = count + 1
      endif
    enddo
    !iterate through each element determining which habitat polygons are
    !  inside them 
    do i=1,max_rho_elements
      count = 0
      if(r_ele_x(1,i)==0.0 .and. r_ele_x(2,i)==0.0 .and. r_ele_x(3,i)==0.0)then
         elepolys(klev,i)%numpoly = count
         exit
      endif
      do j=1,pedges
        !if the current polygon edge point has the same id as the last one
        !  added to the list: skip this iteration (so same polygon isnt added
        !  multiple times)
        if(count.NE.0)then
          if(polys(j,1).EQ.polytail%num) cycle
        endif

        !check if each habitat polygon edge point is in the element
        triangle = 0
        checkele = i
        CALL gridcell(max_rho_elements,r_ele_y,r_ele_x,polys(j,4),polys(j,5),      &
                      P_ele,triangle,checkele)
        if(triangle /= 0) then
          count = count + 1

          IF (.NOT. ASSOCIATED(polyhead)) THEN
            ALLOCATE(polyhead,STAT=ISTAT)
            polytail => polyhead
            NULLIFY(polytail%p)
            polytail%num = NINT(polys(j,1))
          ELSE
            ALLOCATE(polytail%p,STAT=ISTAT)
            polytail => polytail%p
            NULLIFY(polytail%p)
            polytail%num = NINT(polys(j,1))
          ENDIF

          cycle
        endif

        !if none of the current habitat polygon edge points are in the element
        !  check to make sure none of the element edge points are in the
        !  habitat polygon
        check = .FALSE.
        if(j==pedges)then
          check = .TRUE.
        else
          if(polys(j,1)/=polys(j+1,1)) check = .TRUE.
        endif
        if(check)then
          check = .FALSE.
          !check if any of the element edge points are in range of the habitat
          !  polygon
          do k=1,4
            dis = sqrt( (r_ele_x(k,i)-polys(j,2))**2 +     &
                        (r_ele_y(k,i)-polys(j,3))**2 )
            if(dis<maxbdis(nint(polys(j,1)))) check = .TRUE.
          enddo

          !if the element edge points are in range:
          if(check)then
            !allocate poly to contain the current habitat polygon edge points
            ALLOCATE(poly(2,polyspecs(nint(polys(j,1)),2)))
            do k=1,polyspecs(nint(polys(j,1)),2)
              poly(1,k) = polys(polyspecs(nint(polys(j,1)),1)+k-1,4)
              poly(2,k) = polys(polyspecs(nint(polys(j,1)),1)+k-1,5)
            enddo
            !check if any of the four element edge points are in the polygon
            do k=1,4
              if(inpoly(r_ele_x(k,i),r_ele_y(k,i),    &
                 polyspecs(nint(polys(j,1)),2),poly))then
                !if one is inside the polygon, add the polygons id to polynums
                count = count + 1

                IF (.NOT. ASSOCIATED(polyhead)) THEN
                  ALLOCATE(polyhead,STAT=ISTAT)
                  polytail => polyhead
                  NULLIFY(polytail%p)
                  polytail%num = NINT(polys(j,1))
                ELSE
                  ALLOCATE(polytail%p,STAT=ISTAT)
                  polytail => polytail%p
                  NULLIFY(polytail%p)
                  polytail%num = NINT(polys(j,1))
                ENDIF

                exit
              endif
            enddo
            !make sure poly is deallocated
            DEALLOCATE(poly)
          endif
        endif          

      enddo
      
      ! write(*,*) 'about to transfer info, element',i,' of',rho_elements !CRS
      !if any polygons were inside this element, transfer that information
      !  from polynums to elepolys
      elepolys(klev,i)%numpoly = count
      if(count>0)then
        ALLOCATE(elepolys(klev,i)%poly(count))
        j=0
        do
          j=j+1
          if(.NOT. ASSOCIATED(polyhead) ) exit
          ptr => polyhead
          polyhead => polyhead%p
          elepolys(klev,i)%poly(j) = ptr%num
          DEALLOCATE(ptr)
        enddo
        NULLIFY(polyhead,polytail)
      endif
    enddo
    NULLIFY(ptr,polyhead,polytail)

   ENDDO ! klev=1,us_tridim
   !------------------

    write(*,*) 'checking if holes exist' !CRS

    if(holesExist)then

      !iterate through each hole and store information:
      !  holespecs(n,1) = location in holes of first point of hole n
      !  holespecs(n,2) = # of edge points that make up hole n
      count = 1
      holespecs(nint(Holes(1,1)),1) = 1
      do i=2,hedges
        if(Holes(i,1) /= Holes(i-1,1))then
          holespecs(nint(Holes(i,1)),1) = i
          holespecs(nint(Holes(i-1,1)),2) = count
          count = 1
        elseif(i==hedges)then
          holespecs(nint(Holes(i,1)),2) = count + 1
        else
          count = count + 1
        endif
      enddo


      !iterate through every habitat polygon determining which holes are 
      !  inside them
      do i=minpolyid,maxpolyid
        if(polyspecs(i,1)==0)cycle

        count = 0

        !iterate through all the hole edge points
        do j=1,hedges
          if(j == 1 .OR. holes(j,1) /= holes(j-1,1)) then
            if(holes(j,6)==polys(polyspecs(i,1),1))then
              !if the hole is in current habitat polygon, add it to polynums
              count = count + 1

              IF (.NOT. ASSOCIATED(polyhead)) THEN
                ALLOCATE(polyhead,STAT=ISTAT)
                polytail => polyhead
                NULLIFY(polytail%p)
                polytail%num = NINT(holes(j,1))
              ELSE
                ALLOCATE(polytail%p,STAT=ISTAT)
                polytail => polytail%p
                NULLIFY(polytail%p)
                polytail%num = NINT(holes(j,1))
              ENDIF

              cycle
            endif
          endif
        enddo

        !if there were any holes in the current habitat polygon, transfer that
        !  information from polynums to polyholes
        polyholes(i)%numpoly = count
        if(count>0)then
          ALLOCATE(polyholes(i)%poly(count))
          j=0
          do
            j=j+1
            if(.NOT. ASSOCIATED(polyhead) ) exit
            ptr => polyhead
            polyhead => polyhead%p
            polyholes(i)%poly(j) = ptr%num
            DEALLOCATE(ptr)
          enddo
          NULLIFY(polyhead,polytail)
        endif
      enddo

    endif

    DEALLOCATE(r_ele_x)
    DEALLOCATE(r_ele_y)
    write(*,*) 'about to return from createPolySpecs' !CRS

  END SUBROUTINE createPolySpecs

  LOGICAL FUNCTION isStranded(n)
  !This function returns .TRUE. if the particle is "stranded", and FALSE if not
    USE PARAM_MOD, ONLY:settlementon
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n

    if(settlementon)then
     isStranded = stranded(n)
    else
     isStranded = .FALSE.
    endif

  END FUNCTION isStranded

  !Subroutine to determine if the particle is on any oyster polys
  !  (including holes)
  SUBROUTINE testSettlement(P_age,n,Px,Py,Pz,inpoly, &
                            P_depth,k,intersectf,coastdist)
    USE PARAM_MOD, ONLY: holesExist,StrandingDist,strandingmaxdepth,    &
                         strandingmaxdistfromdepth,OilOn
    USE HYDRO_MOD, ONLY: getP_r_element,getP_klev
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n,k
    INTEGER, INTENT(OUT) :: inpoly
    DOUBLE PRECISION, INTENT(IN) :: P_age,Px,Py,Pz,P_depth,coastdist
    INTEGER, INTENT(IN) :: intersectf                                    ! CL-OGS
    INTEGER :: klev

    INTEGER :: polyin,R_ele

    R_ele = getP_r_element(n)
    klev= getP_klev(n)
    polyin = 0
    inpoly = 0

      !check if the particle is within the boundaries of any habitat polygon
      CALL psettle(Px,Py,R_ele,klev,polyin)
      !if within a habitat polygon, check if particle is within the boundaries
      ! of any hole that is in that particular habitat polygon
      if (polyin .GT. 0) then          !if it is within a habitat polygon
        inpoly = polyin                !  set inpoly to the polygon id #
        if(holesExist)then
          CALL hsettle(Px,Py,polyin)
          if (polyin /= 0) inpoly = 0  !if its within a hole, reset inpoly to 0
        endif
      endif

    if( (inpoly .GT. 0) .and. &
      ( (abs(P_age) >= settletime(n)) .and.(P_depth+abs(strandingmaxdistfromdepth)>=Pz )  ) )then 
      if (.not.OilOn .and. StrandingDist<0) then
            settle(n) = .TRUE.       ! If (StrandingDist=0) n can settle   
      elseif(StrandingDist>0)then! If (StrandingDist>0) n can strand
            if(coastdist<=StrandingDist)   stranded(n)=.TRUE.
      elseif(StrandingDist==0 .and. abs(coastdist)<0.5)then! If (StrandingDist>0) n can strand
            stranded(n)=.TRUE.
      endif

    endif

  END SUBROUTINE testSettlement

  SUBROUTINE p_Stranding(n)
    IMPLICIT NONE
    INTEGER :: n
    stranded(n)=.TRUE.
  END SUBROUTINE p_Stranding


  SUBROUTINE psettle(Px,Py,R_ele,klev,polyin)
    USE PIP_MOD, ONLY: inpoly
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: R_ele,klev
    INTEGER, INTENT(OUT) :: polyin
    DOUBLE PRECISION, INTENT(IN) :: Px,Py

    INTEGER :: i,j,start,size
    DOUBLE PRECISION :: dis
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: polybnds

    polyin = 0               !initialize polyin to 0

    !if there are any habitat polygons in the element that the particle is in:
    if(elepolys(klev,R_ele)%numpoly > 0)then

      !iterate through all the habitat polygons in that element
      do i=1,elepolys(klev,R_ele)%numpoly
        start = polyspecs(elepolys(klev,R_ele)%poly(i),1)
        size =  polyspecs(elepolys(klev,R_ele)%poly(i),2)

        !if the particle is not within range of the habitat polygon,
        !  skip this habitat polygon
        dis = sqrt( (Px-polys(start,2))**2 + (Py-polys(start,3))**2 )
        if(dis>maxbdis(NINT(polys(start,1))))cycle

        !allocate polybnds and fill it with the boundary point locations of 
        !  the current habitat polygon
        ALLOCATE(polybnds(2,size))
        do j=1,size
          polybnds(1,j) = polys(start + j - 1, 4)
          polybnds(2,j) = polys(start + j - 1, 5)
        enddo
        ! call inpoly to see if the point is in the polygon
        if(inpoly(Px,Py,size,polybnds))then ! if it is in the polygon:
          polyin = NINT(polys(start,1))           !   set polyin to the polygon id
          DEALLOCATE(polybnds)              !   deallocate polybnds
          exit                              !   and exit the subroutine
        endif

        !make sure polybnds is deallocated
        DEALLOCATE(polybnds)
      enddo
    endif

  END SUBROUTINE psettle


  SUBROUTINE hsettle(Px,Py,holein)
    USE PIP_MOD, ONLY: inpoly
    IMPLICIT NONE

    INTEGER, INTENT(INOUT) :: holein
    DOUBLE PRECISION, INTENT(IN) :: Px,Py

    INTEGER :: i,j,start,size,polyin
    DOUBLE PRECISION :: dis
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: polybnds

    polyin = holein          !initialize polyin to holein
    holein = 0               !initialize holein to 0

    !if there are any holes in the habitat polygon that the particle is inside:
    if(polyholes(polyin)%numpoly > 0)then

      !iterate through each hole
      do i=1,polyholes(polyin)%numpoly
        start = holespecs(polyholes(polyin)%poly(i),1)
        size =  holespecs(polyholes(polyin)%poly(i),2)

        !if the particle is not within range of the hole, skip this hole
        dis = sqrt( (Px-Holes(start,2))**2 + (Py-Holes(start,3))**2 )
        if(dis>maxhdis(NINT(Holes(start,1))))cycle

        !allocate polybnds and fill it with the boundary point locations of 
        !  the current hole
        ALLOCATE(polybnds(2,size))
        do j=1,size
          polybnds(1,j) = Holes(start + j - 1, 4)
          polybnds(2,j) = Holes(start + j - 1, 5)
        enddo

        ! call inpoly to see if the point is in the hole
        if(inpoly(Px,Py,size,polybnds,.FALSE.))then
        !  NOTICE: onin is set .FALSE. meaning a 
        !  particle on the edge of a hole in a 
        !  habitat polygon is not considered to be 
        !  in the hole and thus will still settle
          holein = NINT(Holes(start,1))      ! if in: set holein to the hole id
          DEALLOCATE(polybnds)         !        deallocate polybnds
          exit                         !        and exit the subroutine
        endif

        !make sure polybnds is deallocated
        DEALLOCATE(polybnds)
      enddo
    endif

  END SUBROUTINE hsettle


  LOGICAL FUNCTION isSettled(n)
  !This function returns .TRUE. if the particle has "settled", and FALSE if not
    USE PARAM_MOD, ONLY:settlementon
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n

    if(settlementon)then
      isSettled = settle(n)
    else
      isSettled = .FALSE.
    endif


  END FUNCTION isSettled


END MODULE SETTLEMENT_MOD
