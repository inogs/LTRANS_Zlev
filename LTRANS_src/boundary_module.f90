MODULE BOUNDARY_MOD

! This module contains variables and subroutines associated with the 
! creation of the land/sea boundaries.  The main purpose of this module 
! is to create the land/sea boundaries from a given masked rho grid.  
! The main subroutine in the module, createBounds, determines the number
! of boundary points, allocates an array of that size, and fills it with 
! the boundary points in order.

! Module and most subroutines created by:  Zachary Schlag
! Created on:           03 Apr 2008
! Last Modified on:        Feb 2011
!
! Subroutine intersect_reflect created by: Elizabeth North
! Created on:                  2005
! Last Modified on:     08 Aug 2008

IMPLICIT NONE
PRIVATE
SAVE

!*****************************************************************
!*                          VARIABLES                            *
!*****************************************************************

  TYPE :: boundary  !boundary point variable
    LOGICAL :: onU      !TRUE if point is on U grid, FALSE if on V Grid
    INTEGER :: ii       !i grid position of point
    INTEGER :: jj       !j grid position of point
    INTEGER :: poly     !polygon of point (1 = main body of water)
                        !                 (>1 = island)
  END TYPE boundary

  !initial holder of all boundary points, 
  !  deallocated after reformatted into hid,hx,hy,bx,by
  TYPE (boundary), ALLOCATABLE,  DIMENSION (:,:) :: bnds

  !final boundary variables, after reformatting from bnds
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: hx,hy,bx,by
  INTEGER, ALLOCATABLE, DIMENSION (:,:) :: mid,hid
  INTEGER, ALLOCATABLE, DIMENSION (:,:) :: mdeg,hdeg
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION (:,:) :: bnd_x1,bnd_x2,bnd_y1,bnd_y2

  !TRUE if the boundary is land, FALSE if it is open ocean
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: land

  INTEGER, ALLOCATABLE, DIMENSION (:) :: bnum !number of boundary points in k-level
  INTEGER, ALLOCATABLE, DIMENSION (:) :: num_mbnd ! number of main boundary polynoms found in k-level
  INTEGER , ALLOCATABLE, DIMENSION (:) :: tot_mbnd_pts    ! number of main boundary elements in k-level

  INTEGER :: bnum_max       ! max number of boundary points among k-levels

  INTEGER, ALLOCATABLE, DIMENSION(:) ::  & 
             tot_ibnd_pts,   & !number of boundary points around islands
             num_ibnd,  & !number of islands
             nbounds        !total number of boundary segments

  INTEGER::  & 
             max_mbnd_pts,    & !number of boundary points around water
             max_ibnd_pts,   & !number of boundary points around islands
             nbounds_max        !total number of boundary segments

  LOGICAL, ALLOCATABLE, DIMENSION(:) :: BND_SET  !Boundaries Set?
 
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: polynestdeg

  ! CL-OGS: following declarations moved here to make ele avaible to add routine
  TYPE :: element       !element variable
    INTEGER :: form     !0-19 depends on masking of 4 nodes & situation
    LOGICAL :: used     !TRUE if element is already used, FALSE if not
                        !Note that in a 'cross' element, used will not
                        !  change to TRUE until it has been used twice
    LOGICAL :: unused   !TRUE if a 'cross' element has never been used,
                        !  FALSE if it has been used at least once
  END TYPE element

  TYPE (element), ALLOCATABLE, DIMENSION(:,:,:) :: ele        !elements

  LOGICAL :: wf = .false.   !short for waterfall, indicates path is going
                            !  through an element on the edge of rho grid
  INTEGER :: dir = 5        !direction boundary left the previous element
                            !  based on numeric keypad (8-up,6-right,etc)  
  ! CL-OGS:  -------------end of moved section
  INTEGER :: pyout=1000               ! CL-OGS: python output file 
  INTEGER :: us_tridim,ws_tridim      ! CL-OGS: =us,ws for Zgrid, =1,1 for ROMS

  !The following procedures have been made public:
  PUBLIC :: isBndSet,createBounds,mbounds,ibounds,&
            intersect_reflect,finBoundary,Get_coastdist

CONTAINS

!*****************************************************************
!*                    FUNCTIONS & SUBROUTINES                    *
!*****************************************************************



    !*********************************************************
    !*                      isBndSet                         *
    !*********************************************************

LOGICAL FUNCTION isBndSet()
  USE PARAM_MOD, ONLY:us
  INTEGER :: k
!This function simply returns TRUE if the boundaries have been created
!  and FALSE if they have not yet been created
  isBndSet=.true.
  do k=1,us
    if (.NOT.BND_SET(k))  isBndSet = BND_SET(k)
  enddo
  RETURN
END FUNCTION isBndSet


    !*********************************************************
    !*                     createBounds                      *
    !*********************************************************

SUBROUTINE createBounds()
  USE PARAM_MOD, ONLY:ui,uj,vi,vj,us,ws,BoundaryBLNs,NCOutFile,BndOut, &
                      Zgrid
  USE HYDRO_MOD, ONLY:getMask_Rho,getUVxy,seteleform,setbounds
!This subroutine creates boundaries based on the masking of the rho grid
!  The boundary points are U & V grid points directly between the rho points
!  This subroutine assumes there is only one body of water in the grid 
!    but there may exist multiple islands within the water
!  It creates the boundary points to go clockwise around the body of water
!    It then will create boundary points clockwise around any islands
!
!The subroutine does the following things:
!  1) Determines the number of boundary points in the given rho grid
!  2) Determines the form of each rho element
!    a) Elements consist of the rho nodes:  (i+1,j) --- (i+1,j+1)
!                                              |            |
!                                            (i,j)  ---  (i,j+1)
!
!    b) Element form is based on how the four nodes in the element are masked
!
!    c) If any elements exist where their form contains water and land crossing
!         diagonally, which are referred to as 'crosses' then these forms must
!         be solved as to which direction the boundarys are going through them
!
!  3) Starting at the first element it encounters with a water masked node, 
!       it goes around the body of water clockwise adding each boundary point
!       until it reaches the boundary point that it started on
!
!  4) If the number of boundary points used is not equal to the total number of
!       boundary points originally found in the grid, then find an unused
!       boundary element to begin on and go clockwise around the island it is 
!       part of, this is continued until all the boundary points are used


  LOGICAL :: U = .true.     !sent to subroutine ADD to indicate on U grid
  LOGICAL :: V = .false.    !sent to subroutine ADD to indicate on V grid
  !LOGICAL :: wf = .false.   !short for waterfall, indicates path is going
  !                          !  through an element on the edge of rho grid
  !INTEGER :: dir = 5        !direction boundary left the previous element
  !                          !  based on numeric keypad (8-up,6-right,etc)

  INTEGER :: ipos,jpos      !i and j position of current element when 
                            ! making bounds

  !The following variables are the same as the similarly named variables
  INTEGER :: ipos1,jpos1,ipos2,jpos2,dir1,dir2  ! above, only these are
  LOGICAL :: wf1, wf2                           ! used for solving crosses
  INTEGER :: itmp,jtmp

  INTEGER, ALLOCATABLE, DIMENSION(:) :: crossnum
  INTEGER :: i,j,STATUS,m,oldcrossnum,count,forms0,forms15,otherforms
  LOGICAL :: found, deadend1, deadend2

  LOGICAL :: polydone       !TRUE if current polygon is closed, else FALSE
  LOGICAL :: done = .false. !TRUE if all boundaries are used, else FALSE
  LOGICAL :: continentalborder   ! to look for continental borders
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: is_continent
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: nestingdegree
  
  !Used to reformat from bnds to bx,by,hx,hy,hid
  INTEGER :: k,klevel,c,pstart,pend
  INTEGER :: bnd,isl
  INTEGER :: numpoly_max
  INTEGER, ALLOCATABLE, DIMENSION(:) :: numpoly
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: polysizes
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: x_u,y_u,x_v,y_v

  !Rho Mask Used to create boundaries
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: mask_rho

  !Keep track of boundary points on open ocean
  LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: oo
  CHARACTER(LEN=250) :: formatoutput,namefile
  INTEGER :: io,nmp,counter,count_land,count_water,a,b,d,degree,dmin,dmax
  
  if(Zgrid)then
    us_tridim=us
    ws_tridim=ws
  else
    us_tridim=1
    ws_tridim=1
  endif

  ALLOCATE(x_u(ui,uj))
  ALLOCATE(y_u(ui,uj))
  ALLOCATE(x_v(vi,vj))
  ALLOCATE(y_v(vi,vj))
  ALLOCATE(mask_rho(vi,uj,us_tridim))
  ALLOCATE(crossnum(us_tridim))
  ALLOCATE(BND_SET(us_tridim))
  if(BndOut)then
  DO k=1,us_tridim 
    write(namefile,'(a,a,i3.3,a)')TRIM(NCOutFile),"_boundary_module",k,".py"
    open (unit = pyout+k, file = trim(namefile))
  ENDDO
  endif
  BND_SET=.false.

  CALL getMask_Rho(mask_rho)

  !Edit Rho Mask - remove nodes that don't have at least 2 neighbors
  !  This removes a situation that may occur where an area that is 
  !  within the created boundaries, is not covered by either the U 
  !  or V grid.
  write(*,*)'remove nodes that dont have at least 2 neighbors' 
  do
   m=0
   DO k=1,us_tridim
    do j=1,uj
      do i=1,vi
        if(mask_rho(i,j,k)== 1)then
          count = 0
          if(i> 1)then
            if(mask_rho(i-1,j,k)==1) count = count + 1
          endif
          if(i<vi)then
            if(mask_rho(i+1,j,k)==1) count = count + 1
          endif
          if(j> 1)then
            if(mask_rho(i,j-1,k)==1) count = count + 1
          endif
          if(j<uj)then
            if(mask_rho(i,j+1,k)==1) count = count + 1
          endif
          if(count < 2) then
            mask_rho(i,j,k) = 0
            m = m + 1
            write(*,*)'Deleting mask_rho at node i,i,k=',i,j,k
            if(BndOut)then
              write(pyout+k,'(a,i3,a)')"if (thirddimindex==",k-1,"):"
              write(pyout+k,'(a,i4,a,i4,a,i4,a,i4,a)')                         &
                                     " if (v=='mask_rho'):datafield[",j-1,",", &
                                                   i-1,"]= 0.2 # PYTHON"
              if(i> 1 )then
                      if(  mask_rho(i-1,j,k)==0 )                             &
                                     write(pyout+k,'(a,i4,a,i4,a)')            &
                                     "	if (v=='mask_u'):datafield[",j-1,",",  &
                                     i-2,"]= 0.2 # PYTHON"                   
              endif 
              if(i<vi )then
                      if( mask_rho(i+1,j,k)==0 )                              &
                                     write(pyout+k,'(a,i4,a,i4,a)')            &
                                      "	if (v=='mask_u'):datafield[",j-1,",",  &
                                      i-1,"]= 0.2 # PYTHON"                   
              endif 
              if(j> 1)then
                      if( mask_rho(i,j-1,k)==0 )                               &
                                     write(pyout+k,'(a,i4,a,i4,a)')            &
                                      "	if (v=='mask_v'):datafield[",j-2,",",  &
                                      i-1,"]= 0.2 # PYTHON"                   
              endif 
              if(j<uj)then
                      if( mask_rho(i,j+1,k)==0 )                               &
                                     write(pyout+k,'(a,i4,a,i4,a)')            &
                                      "	if (v=='mask_v'):datafield[",j-1,",",  &
                                      i-1,"]= 0.2 # PYTHON"
              endif 
            endif!BndOut
          endif
        endif
      enddo
    enddo
   ENDDO
   if(m==0)exit
  enddo

    !*****************************************
    !*              SET BOUNDS               *
    !*****************************************

  !This section determines the number of boundary points in the given rho grid

!  write(*,*) 'Calculating Number of Boundary Points'
!  write(*,*) ' '

  !Allocate the variable bnum to the total number of k-levels
  ALLOCATE (bnum(us_tridim),STAT=STATUS)
  if(STATUS /= 0) write(*,*) 'Problem allocating bnum'
  DO k=1,us_tridim
    bnum(k) = 0  !initialize bnum to 0
  ENDDO


  !Allocate the variable num_mbnd to the total number of k-levels
  ALLOCATE (num_mbnd(us_tridim),STAT=STATUS)
  if(STATUS /= 0) write(*,*) 'Problem allocating num_mbnd'
  DO k=1,us_tridim
    num_mbnd(k) = 0  !initialize to 0
  ENDDO 


  DO k=1,us_tridim
  do i=2,vi-1
    do j=2,uj-1

      !Excluding the edges of the rho grid, increment bnum every time adjacent
      !  nodes are found that are of differing mask values
      if(i<vi-1)then
        if((mask_rho(i,j,k).EQ.0 .AND. mask_rho(i+1,j,k).EQ.1) .OR.            &
       (mask_rho(i,j,k).EQ.1 .AND. mask_rho(i+1,j,k).EQ.0)) bnum(k) = bnum(k)+1
      endif

      if(j<uj-1)then
        if((mask_rho(i,j,k).EQ.0 .AND. mask_rho(i,j+1,k).EQ.1) .OR.            &
       (mask_rho(i,j,k).EQ.1 .AND. mask_rho(i,j+1,k).EQ.0)) bnum(k) = bnum(k)+1
      endif

    enddo
  enddo
  ENDDO 

  !bnum must also be incremented where the second node from the edge of the 
  !  rho grid is masked as water because the boundary will need to go between
  !  this point and the edge point
  DO k=1,us_tridim
  do i=2,vi-1
    if(mask_rho(i,2,k).EQ.1) bnum(k) = bnum(k)+1
    if(mask_rho(i,uj-1,k).EQ.1) bnum(k) = bnum(k)+1
  enddo
  
  do j=2,uj-1
    if(mask_rho(2,j,k).EQ.1) bnum(k) = bnum(k)+1
    if(mask_rho(vi-1,j,k).EQ.1) bnum(k) = bnum(k)+1
  enddo
  ENDDO

  ! compute max(bnum(k)) for bnds(bnum_max,us_tridim) allocation
  bnum_max=0
  DO k=1,us_tridim
    if ( bnum(k)>bnum_max )  bnum_max=bnum(k)
  ENDDO

  !Allocate the variable bnds to the total number of boundaries (bnum)
  ALLOCATE (bnds(bnum_max,us_tridim),STAT=STATUS)
  if(STATUS /= 0) write(*,*) 'Problem allocating bnds'
  if(BndOut)write(*,*)'allocated bnds(',bnum_max,',',us_tridim,')'
  !Initialize all the values in newly allocated bnds
  DO k=1,us_tridim
  do i=1,bnum(k)
    bnds(i,k)%onU  = .false.
    bnds(i,k)%ii   = 0
    bnds(i,k)%jj   = 0
    bnds(i,k)%poly = 0
  enddo
  ENDDO
  
  if(BndOut)write(*,*) 'GET FORMS'

    !*****************************************
    !*              GET FORMS                *
    !*****************************************

  !This section determines the form of each rho element
  !
  !Remember that:
  !  a) Elements consist of the rho nodes:  (i+1,j,k) --- (i+1,j+1,k)
  !                                            |            |
  !                                          (i,j,k)  ---  (i,j+1,k)
  !
  !  b) Element form is based on how the four nodes in the element are masked
  !
  !  c) If any elements exist where their form contains water and land crossing
  !       diagonally, which are referred to as 'crosses' then these forms must
  !       be solved as to which direction the boundarys are going through them

  !Allocate the variable ele to the total number of elements
  ALLOCATE (ele(vi-1,uj-1,us_tridim),STAT=STATUS)
  if(BndOut) write(*,*) 'ele(',vi-1,uj-1,us_tridim,') has been allocated'
  if(STATUS /= 0) write(*,*) 'Problem allocating ele'

  ALLOCATE (is_continent(vi,uj),STAT=STATUS)
  if(STATUS /= 0) write(*,*) 'Problem allocating is_continent'
  ALLOCATE (nestingdegree(vi,uj,us_tridim),STAT=STATUS)
  if(STATUS /= 0) write(*,*) 'Problem allocating nestingdegree'

  if(BndOut)write(*,*) 'Getting Element Forms'
  if(BndOut)write(*,*) ' k=1,',us_tridim

  crossnum = 0  !Initialize 'cross' element counter to 0

  if(BndOut)write(*,*) 'loop k=1,',us_tridim
  if(BndOut)write(*,*)''

  DO k=1,us_tridim
  if(BndOut)write(*,*)' '
  if(BndOut)write(*,*)' '
  if(BndOut)write(*,*)' '
  if(BndOut)write(*,*)'----------------------------------------------------'
  if(BndOut)write(*,'(a,i3,a)')'-----------LEVEL ',k,' ----------------------'
  if(BndOut)write(*,*)'----------------------------------------------------'
  !f(BndOut)then 
  !do jtmp=1,uj
  !   if(jtmp==1)then
  !     !
  !     write(*,'(a)',advance="no")'----|'
  !     do itmp=1,vi
  !        write(*,'(a)',advance="no")'---'
  !     enddo 
  !     write(*,*)' |'
  !     !
  !     write(*,'(a)',advance="no")'    |'
  !     do itmp=1,vi
  !        write(*,'(i3)',advance="no")itmp
  !     enddo
  !     write(*,*)' |'
  !     !
  !     write(*,'(a)',advance="no")'----|'
  !     do itmp=1,vi
  !        write(*,'(a)',advance="no")'---'
  !     enddo 
  !     write(*,*)' |'
  !   endif
  !   !
  !   write(*,'(i3,a)',advance="no")jtmp,' |'
  !   do itmp=1,vi
  !      if(int(mask_rho(itmp,jtmp,k))==0) then
  !            write(*,'(a3)',advance="no")'X'
  !      else  
  !            write(*,'(a3)',advance="no")'.'
  !      endif
  !   enddo
  !   write(*,*)' |'
  !   !
  !enddo
  !!
  !write(*,'(a)',advance="no")'----|'
  !do itmp=1,vi
  !   write(*,'(a)',advance="no")'----'
  !enddo 
  !write(*,*)' |'
  !!
  !ndif
 
  if(BndOut)write(pyout+k,'(a,i3,a)')"if (thirddimindex==",k-1,") : #PYTHON"

  forms0=0
  forms15=0
  otherforms=0
  do i=1,vi-1
    do j=1,uj-1
      ele(i,j,k)%used = .false.   !Initialize used to FALSE
      ele(i,j,k)%unused = .true.  !Initialize unused to TRUE (for 'crosses')
      !Determine the form of the current element
      ele(i,j,k)%form = 0
      if(mask_rho(i  ,j+1,k) == 0) ele(i,j,k)%form = ele(i,j,k)%form + 1
      if(mask_rho(i  ,j  ,k) == 0) ele(i,j,k)%form = ele(i,j,k)%form + 2
      if(mask_rho(i+1,j+1,k) == 0) ele(i,j,k)%form = ele(i,j,k)%form + 4
      if(mask_rho(i+1,j  ,k) == 0) ele(i,j,k)%form = ele(i,j,k)%form + 8
  
      SELECT CASE(ele(i,j,k)%form)
  
        CASE(6,9) !ud
        !Elements with forms 6 or 9 are 'crosses', increment the cross counter
          crossnum(k) = crossnum(k) + 1
        CASE(0) !0
        !Elements with forms 0 or 15 are all water or all land
        !Initialized used to TRUE because there are no boundaries through them
          ele(i,j,k)%used = .true.
          forms0=forms0+1
        CASE(15) !0
        !Elements with forms 0 or 15 are all water or all land
        !Initialized used to TRUE because there are no boundaries through them
          ele(i,j,k)%used = .true. ! BEUG ERROR: in i=1 or j=1, boundary pass through water element! 
          forms15=forms15+1
        CASE DEFAULT
          otherforms=otherforms+1
      END SELECT
    enddo
  enddo


  if(BndOut)write(*,'(a,i7)') '       crosses at this level    =',crossnum(k)
  if(BndOut)write(*,'(a,i7)') '       forms=0  at this level   =',forms0 
  if(BndOut)write(*,'(a,i7)') '       forms=15 at this level   =',forms15
  if(BndOut)write(*,'(a,i7)') '       other forms at this level=',otherforms
  if(BndOut)&
  write(*,'(6(a,i7))')'k=',k,',bnum(k)=',bnum(k),', crossnum(k)=',crossnum(k), &
      ' , forms0=',forms0,' , form15=',forms15,' , otherforms=',otherforms

  if(BndOut)write(*,*) ''
  if(BndOut)write(*,*) 'SOLVE CROSSES'
      !*********************************
      !*         SOLVE CROSSES         *
      !*********************************

  oldcrossnum = 0   !Initialize oldcrossnum to 0
                  !This variable is to know if an endless loop has occurred
  if(crossnum(k) > 0) then
    if(BndOut)write(*,*) ' '
  endif
  counter=0
  do 
    if(crossnum(k) == 0) exit          !If no crosses exist unsolved, exit
    if(BndOut)                                                                 &
      write(*,*)' at this level :',crossnum(k),' crosses remain to be solved'
    if(oldcrossnum == crossnum(k))then !If crosses still exist, but none were
      if(BndOut)write(*,*)'  Cross loop - using defaults'    !  solved in prior loop
      if(BndOut)write(*,*)' '                                !  use defaults of 16 & 19
      do i=1,vi-1
        do j=1,uj-1
          if(ele(i,j,k)%form == 6) then
              ele(i,j,k)%form = 16
              crossnum(k)=crossnum(k)-1 ! BEUG? this was added
          endif
          if(ele(i,j,k)%form == 9) then
              ele(i,j,k)%form = 19
              crossnum(k)=crossnum(k)-1  ! BEUG? this was added 
          endif
        enddo
      enddo
    endif    
    oldcrossnum = crossnum(k)

    if(crossnum(k) == 0)exit         !NOT BEUG but time gain: If no crosses exist unsolved, exit
    do i=1,vi-1
      if(crossnum(k) == 0)exit         !If no crosses exist unsolved, exit

      do j=1,uj-1
        if(crossnum(k) == 0)exit       !If no crosses exist unsolved, exit


        if(ele(i,j,k)%form == 6) then
          !Form 6 has water nodes in the top left and bottom right.
          !  Boundary edges travelling clockwise around water therefore
          !  enter the element from top or bottom and exit left or right


          if(i==1 .OR. j==uj-1 .OR. i==vi-1 .OR. j==1)then ! BEUG: CL: why uj-1
                                                           ! and not uj? why vi-1 and not vi?
            !if the cross is in an element on the edge of the rho grid,
            !  then the boundary will either enter the top and exit left,
            !  or enter the bottom and exit right
            ele(i,j,k)%form = 17 
            crossnum(k) = crossnum(k) - 1
          else


            !if its not on the edge of the grid it will have to be solved
            !  the hard way.  This is done by exiting the two possible
            !  exits (left & right) and following the boundaries until
            !  either, one returns to this cross or, both hit another
            !  unsolved cross.  If one returns to this cross then this 
            !  cross is solved.  If neither direction returns to this 
            !  cross, skip this cross.
            !  As other crosses are solved this one becomes more likely
            !  to find a path back to itself.

            ipos1 = i           !initialize two paths to the cross location
            jpos1 = j
            ipos2 = i
            jpos2 = j
            wf1 = .false.       !initialize wf & deadend variables to false
            wf2 = .false.
            deadend1 = .false.
            deadend2 = .false.

            ! Initialized ipos1,jpos1 to the element Left of the cross
            if(jpos1 == 2)then         !If element to left is on edge
              wf1 = .true.             !  switch wf to TRUE and move up
              jpos1 = jpos1 - 1
              ipos1 = ipos1 + 1
              dir1 = 8
              if(ipos1 == vi-1)then    !If element to left and up from
                jpos1 = jpos1 + 1      !  cross is the corner, move right
                dir1 = 6               !  from corner
              endif
            else
              jpos1 = jpos1 - 1        !Else just move left from cross
              dir1 = 4
            endif
          
            ! Initialized ipos2,jpos2 to the element Right of the cross
            if(jpos2 == uj-2)then      !If element to right is on edge
              wf2 = .true.             !  switch wf to TRUE and move down
              jpos2 = jpos + 1
              ipos2 = ipos - 1
              dir2 = 2
              if(ipos2 == 1)then       !If element to right and down from
                jpos2 = jpos2 - 1      !  cross is corner, move left from
                dir2 = 4               !  corner
              endif
            else
              jpos2 = jpos2 + 1        !Else just move right from cross
              dir2 = 6
            endif
            counter=0
            do
              counter=counter+1
              if(ipos1 == i .AND. jpos1 == j)then   !If left path returned to
                if(dir1 == 2)then                   !  cross, set new form
                  ele(i,j,k)%form = 16                !  based on the direction
                elseif(dir1 == 8)then               !  the path returns from
                  ele(i,j,k)%form = 17
                else
                  write(*,*)'Problem Form 6 Left Solution'
                endif
                crossnum(k) = crossnum(k) - 1             !  and decrement crossnum
                exit

              elseif(ipos2 == i .AND. jpos2 == j)then  !If right path
                if(dir2 == 2)then                   !  returned to cross, set
                  ele(i,j,k)%form = 17                !  new form based on the
                elseif(dir2 == 8)then               !  direction the path
                  ele(i,j,k)%form = 16                !  returns from
                else
                  write(*,*)'Problem Form 6 Right Solution'
                endif
                crossnum(k) = crossnum(k) - 1             !  and decrement crossnum
                exit
              endif

              !if the left path has hit a dead end, switch deadend1 to TRUE
              if(ele(ipos1,jpos1,k)%form==6 .OR.             &
                 ele(ipos1,jpos1,k)%form==9) deadend1 = .true.
              !if the right path has hit a dead end, switch deadend2 to TRUE
              if(ele(ipos2,jpos2,k)%form==6 .OR.             &
                 ele(ipos2,jpos2,k)%form==9) deadend2 = .true.

              !if both paths have hit a dead end, this cross cannot be
              !  solved yet, so move on and come back to it later
              if(deadend1 .AND. deadend2) exit
            
              !if left path has not hit dead end, move to next point on path
              if(.NOT. deadend1) then
                if(BndOut) write(*,'(2i5,a,l2,2i5)',advance='no')               &
                           ipos1,jpos1,'   ',wf1,dir1,ele(ipos1,jpos1,k)%form
                CALL getNext(ipos1,jpos1,wf1,dir1,ele(ipos1,jpos1,k)%form)
                if(BndOut)write(*,'(4(a,i3),a,l2)')                             &
                          '   deadend2 TRUE ; i=',ipos1,' , j=',jpos1,         & 
               ', k=',k,' >> element form is ',ele(ipos1,jpos1,k)%form,'  ',wf1
              endif
              !if right path has not hit dead end, move to next point on path
              if(.NOT. deadend2) then
                CALL getNext(ipos2,jpos2,wf2,dir2,ele(ipos2,jpos2,k)%form)
                if(BndOut)write(*,'(4(a,i3))')                                 &
                        'deadend1 TRUE ; i=',ipos1,' , j=',jpos1, &
                        ', k=',k,' >> element form is ',ele(ipos1,jpos1,k)%form
              endif

              if(counter>10000)then
                write(*,*)'program end due to element conflict in boundary'
                write(*,*)'writing mask_rho to help understand the situation'
                do jtmp=max(1,jpos1-20),min(uj,jpos1+20)
                   if(jtmp==max(1,jpos1-20))then
                     !
                     write(*,'(a)',advance="no")'----|'
                     do itmp=max(1,ipos1-20),min(vi,ipos1+20)
                        write(*,'(a)',advance="no")'----'
                     enddo 
                     write(*,*)' |'
                     !
                     write(*,'(a)',advance="no")'    |'
                     do itmp=max(1,ipos1-20),min(vi,ipos1+20)
                        write(*,'(i4)',advance="no")itmp
                     enddo
                     write(*,*)' |'
                     !
                     write(*,'(a)',advance="no")'----|'
                     do itmp=max(1,ipos1-20),min(vi,ipos1+20)
                        write(*,'(a)',advance="no")'----'
                     enddo 
                     write(*,*)' |'
                   endif
                   !
                   write(*,'(i3,a)',advance="no")jtmp,' |'
                   do itmp=max(1,ipos1-20),min(vi,ipos1+20)
                      write(*,'(i4)',advance="no")int(mask_rho(itmp,jtmp,k))
                   enddo
                   write(*,*)' |'
                   !
                enddo
                !
                write(*,'(a)',advance="no")'----|'
                do itmp=max(1,ipos1-20),min(vi,ipos1+20)
                   write(*,'(a)',advance="no")'----'
                enddo 
                write(*,*)' |'
                !
                stop
              endif
            enddo

          endif


        elseif(ele(i,j,k)%form == 9)then
          !Form 9 has water nodes in the bottom left and top right.
          !  Boundary edges travelling clockwise around water therefore
          !  enter the element from left or right and exit top or bottom


          if(i==1 .OR. j==1 .OR. i==vi-1 .OR. j==uj-1)then
            !if the cross is in an element on the edge of the rho grid,
            !  then the boundary will either enter right and exit the top,
            !  or enter left and exit the bottom
            ele(i,j,k)%form = 18
            crossnum(k) = crossnum(k) - 1


            !if its not on the edge of the grid it will have to be solved
            !  the hard way.  This is done by exiting the two possible
            !  exits (up & down) and following the boundaries until
            !  either, one returns to this cross or, both hit another
            !  unsolved cross.  If one returns to this cross then this 
            !  cross is solved.  If neither direction returns to this 
            !  cross, skip this cross.
            !  As other crosses are solved this one becomes more likely
            !  to find a path back to itself.
          else
            ipos1 = i           !initialize two paths to the cross location
            jpos1 = j
            ipos2 = i
            jpos2 = j
            wf1 = .false.       !initialize wf & deadend variables to false
            wf2 = .false.
            deadend1 = .false.
            deadend2 = .false.


            ! Initialized ipos1,jpos1 to the element Above the cross
            if(ipos1 == vi-2)then      !If element above is on edge
              wf1 = .true.             !  switch wf to TRUE & move right
              ipos1 = ipos1 + 1
              jpos1 = jpos1 + 1
              dir1 = 6
              if(jpos1 == uj-1)then    !If element above and to right
                ipos1 = ipos1 - 1      !  of cross is the corner, move
                dir1 = 2               !  down from corner
              endif
            else
              ipos1 = ipos1 + 1        !Else just move above the cross
              dir1 = 8
            endif


            ! Initialized ipos1,jpos1 to the element Below the cross
            if(ipos2 == 2)then         !If element below is on edge
              wf2 = .true.             !  switch wf to TRUE & move left
              ipos2 = ipos2 - 1
              jpos2 = jpos2 - 1
              if(jpos2 == 1)then       !If element below and left of
                ipos2 = ipos2 + 1      !  cross is the corner, move
                dir2 = 8               !  up from corner
              else
                dir2 = 4
              endif
            else
              ipos2 = ipos2 - 1        !Else just move below the cross
              dir2 = 2
            endif


            do
              if(ipos1 == i .AND. jpos1 == j)then   !If up path returned to
                if(dir1 == 4)then                   !  cross, set new form
                  ele(i,j,k)%form = 19                !  based on the direction
                elseif(dir1 == 6)then               !  the path returns from
                  ele(i,j,k)%form = 18
                else
                  write(*,*)'Problem Form 9 Up Solution'
                endif
                crossnum(k) = crossnum(k) - 1             !  and decrement crossnum
                exit

              elseif(ipos2 == i .AND. jpos2 == j)then  !If down path returned
                if(dir2 == 4)then                   !  to cross, set new form
                  ele(i,j,k)%form = 18                !  based on the direction
                elseif(dir2 == 6)then               !  the path returns from
                  ele(i,j,k)%form = 19
                else
                  write(*,*)'Problem Form 9 Down Solution'
                endif
                crossnum(k) = crossnum(k) - 1             !  and decrement crossnum
                exit
              endif

              !if the up path has hit a dead end, switch deadend1 to TRUE
              if(ele(ipos1,jpos1,k)%form==6 .OR.             &
                 ele(ipos1,jpos1,k)%form==9) deadend1 = .true.

              !if the down path has hit a dead end, switch deadend2 to TRUE
              if(ele(ipos2,jpos2,k)%form==6 .OR.             &
                 ele(ipos2,jpos2,k)%form==9) deadend2 = .true.

              !if both paths have hit a dead end, this cross cannot be
              !  solved yet, so move on and come back to it later
              if(deadend1 .AND. deadend2) exit

              !if up path has not hit dead end, move to next point on path
              if(.NOT. deadend1) then
                CALL getNext(ipos1,jpos1,wf1,dir1,ele(ipos1,jpos1,k)%form)
              endif

              !if down path has not hit dead end, move to next point on path
              if(.NOT. deadend2) then
                CALL getNext(ipos2,jpos2,wf2,dir2,ele(ipos2,jpos2,k)%form)
              endif

            enddo

          endif
        endif
      enddo ! i loop
    enddo ! j loop

  enddo ! while (crossnum(k) != 0)

  if(BndOut)write(*,*) 'Crosses Solved'
  if(BndOut)write(*,*) ' '
  if(BndOut)write(*,*) 'Crosses Solved'
        
        

  if(BndOut)write(*,*) '.........'
  if(BndOut)write(*,*) ' '
  
  if(BndOut)write(*,*) 'MAIN BOUNDARY'


  !---------------IDENTIFY CONTINENTAL NODES----------------------  ! CL-OGS
  is_continent(:,:)=.false.
  count_land=0
  nestingdegree(:,:,k)=-1
  count_water=0
  j=1 !LOWER EDGE 
  do i=1,vi
   if(mask_rho(i  ,j  ,k)==0) then
     is_continent(i,j)=.true.
     count_land =  count_land + 1
     nestingdegree(i,j,k)=0           ! main continental land
    else
     count_water = count_water+1
     nestingdegree(i,j,k)=1           ! main water basin
   endif
  enddo
  j=uj !UPPER EDGE 
  do i=1,vi
   if(mask_rho(i  ,j  ,k)==0) then
      is_continent(i,j)=.true.
     count_land =  count_land + 1
     nestingdegree(i,j,k)=0           ! main continental land
    else
     count_water = count_water+1
     nestingdegree(i,j,k)=1           ! main water basin
   endif
  enddo
  i=1  !LEFT EDGE
  do j=1,uj
    if(nestingdegree(i,j,k)<0)then
    if(mask_rho(i  ,j  ,k)==0) then
      is_continent(i,j)=.true.
     count_land =  count_land + 1
     nestingdegree(i,j,k)=0           ! main continental land
    else
     count_water = count_water+1
     nestingdegree(i,j,k)=1           ! main water basin
   endif
   endif
  enddo
  i=vi  !RIGHT EDGE
  do j=1,uj
    if(nestingdegree(i,j,k)<0)then
    if(mask_rho(i  ,j  ,k)==0) then
      is_continent(i,j)=.true.
     count_land =  count_land + 1
     nestingdegree(i,j,k)=0           ! main continental land
    else
     count_water = count_water+1
     nestingdegree(i,j,k)=1           ! main water basin
    endif
    endif
  enddo
  counter=-1
   degree=0
   a=degree-1
   b=degree  
   d=degree+1
   if(BndOut)write(*,*)'after searching nested land/water nodes ',&
                  ' along the boundary set to nesting degree ',&
                     b,d,' max nesting value is', maxval(nestingdegree(:,:,k))
   if(BndOut)write(*,*)'count_land=',count_land,& 
                       'count_water=',count_water,'/',vi*uj
   if(BndOut)write(*,*)'Now propagate g nested land/water nodes ',&
                  ' along the boundary for nesting degree ',&
                     b,d,' max nesting value is', maxval(nestingdegree(:,:,k))
  DO    ! propagate nesting degree on neighboor nodes
   if(counter.eq.0)exit
   counter=0
   do i=2,vi-1
    do j=2,uj-1
     if(nestingdegree(i,j,k).ge.0)cycle
     if(mask_rho(i ,j,k).eq.0)then ! propagate main_continental_land property on neighboor nodes
      if(nestingdegree(i-1,j,k).eq.b .or. nestingdegree(i+1,j,k).eq.b .or.  &
         nestingdegree(i,j-1,k) .eq.b   .or.nestingdegree(i,j+1,k).eq. b )then
         is_continent(i,j)=.true.   ! "no cross" ele form 3,5,7,10,11,12,13,14,15
         counter=counter+1
         count_land =  count_land + 1
         nestingdegree(i,j,k)=b
      elseif((nestingdegree(i-1,j-1,k).eq.b.and.ele(i-1,j-1,k)%form.eq.17 ).or.&
             (nestingdegree(i+1,j+1,k).eq.b.and.ele(i  ,j  ,k)%form.eq.17 ).or.&
             (nestingdegree(i-1,j+1,k).eq.b.and.ele(i-1,j  ,k)%form.eq.18 ).or.&
             (nestingdegree(i+1,j-1,k).eq.b.and.ele(i  ,j-1,k)%form.eq.18 )) then   
         is_continent(i,j)=.true.   ! crossing ele form 17 or 18
         counter=counter+1
         count_land =  count_land + 1
         nestingdegree(i,j,k)=b
      endif
     else   ! propagate is_main_water_basin property on neighboor nodes
      if((nestingdegree(i-1,j,k).eq.b .or.nestingdegree(i+1,j,k).eq.b .or.  &
          nestingdegree(i,j-1,k) .eq.b .or.nestingdegree(i,j+1,k).eq.b ).or. &
         (nestingdegree(i-1,j,k).eq.d .or.nestingdegree(i+1,j,k).eq.d .or.  &
          nestingdegree(i,j-1,k) .eq.d .or.nestingdegree(i,j+1,k).eq.d ))then
         counter=counter+1
         count_water =  count_water + 1
         nestingdegree(i,j,k)=d
      elseif((nestingdegree(i-1,j-1,k).eq.d.and.ele(i-1,j-1,k)%form.eq.19 ).or.&
             (nestingdegree(i+1,j+1,k).eq.d.and.ele(i  ,j  ,k)%form.eq.19 ).or.&
             (nestingdegree(i-1,j+1,k).eq.d.and.ele(i-1,j  ,k)%form.eq.16 ).or.&
             (nestingdegree(i+1,j-1,k).eq.d.and.ele(i  ,j-1,k)%form.eq.16 )) then   
         counter=counter+1
         count_water =  count_water + 1
         nestingdegree(i,j,k)=d
      endif
     endif
    enddo
   enddo
  ENDDO
   if(BndOut)write(*,*)'after propagating boundary land/water nodes ',&
                  ' along the boundary set to nesting degree ',&
                     b,d,' max nesting value is', maxval(nestingdegree(:,:,k))
   if(BndOut)write(*,*)'count_land=',count_land,& 
                       'count_water=',count_water,'/',vi*uj
  counter=-1
  DO
   degree=degree+2
   if((count_land+count_water).eq.(vi*uj))exit  ! exit if all nodes were
                                                       ! identified as either continent, island or water nodes
   counter=-1
   a=degree-1
   b=degree  
   d=degree+1
   if(maxval(nestingdegree(:,:,k))<a)then
         count_land=0
         count_water=0
         do i=1,vi
          do j=1,uj
           if(mask_rho(i ,j,k).eq.0)then 
             count_land=count_land+1
           else
             count_water=count_water+1
           endif
           if(nestingdegree(i,j,k).eq.-1)write(*,'(2(a,i4),9(a,i2),a)')&
              'for node i=',i,'j=',j,' no nesting found, its mask is ',&
              mask_rho(i ,j,k),&
             '. Its neighbours nestings are ', & 
             nestingdegree(i-1 ,j,k),'(i-1) ; ',nestingdegree(i+1,j,k),'(i+1) ; ',&
             nestingdegree(i ,j-1,k),'(j-1) ; ',nestingdegree(i,j+1,k),  &
             '(j+1). Its neighbour elements are :', & 
              ele(i-1,j-1,k)%form,'(i-1,j-1) ; ',ele(i ,j ,k)%form,'(i,j) ; ', &
              ele(i-1,j  ,k)%form,'(i-1,j) ; ',  ele(i ,j-1,k)%form,'(i,j-1) .'
          enddo
         enddo
         write(*,*)'The program was expecting to find land/water nodes =',& 
              count_land,count_water
         call createNestingDegreeNetCDF(nestingdegree)
         write(*,*)'ERROR while searching for nested nodes in boundary module'
         write(*,*)' created NestingDegree file search for error level ',&
                      k,'/',us_tridim
         stop 'ERROR while searching for nested nodes in boundary module'
   endif
   if(BndOut)write(*,*)'Start searching nested land/water regions of degree',&
                    a,b,d,' max nesting value is', maxval(nestingdegree(:,:,k))
   if(BndOut)write(*,*)'count_land=',count_land,& 
                       'count_water=',count_water,'/',vi*uj
   Do
   if(counter.eq.0)exit
   counter=0
   do i=2,vi-1
    do j=2,uj-1
     if(nestingdegree(i,j,k).ge.0)cycle
     if(mask_rho(i ,j,k).eq.0)then ! propagate or increment land property on neighboor nodes
      if(nestingdegree(i-1,j,k).eq.a .or.nestingdegree(i+1,j,k).eq.a .or.  &
         nestingdegree(i,j-1,k).eq.a .or.nestingdegree(i,j+1,k).eq.a .or.  &
         nestingdegree(i-1,j,k).eq.b .or.nestingdegree(i+1,j,k).eq.b .or.  &
         nestingdegree(i,j-1,k).eq.b .or.nestingdegree(i,j+1,k).eq.b)then
         is_continent(i,j)=.true.   ! "no cross" ele form 3,5,7,10,11,12,13,14,15
         counter=counter+1
         count_land =  count_land + 1
         nestingdegree(i,j,k)=b
      elseif((nestingdegree(i-1,j-1,k).eq.b.and.ele(i-1,j-1,k)%form.eq.17 ).or.&
             (nestingdegree(i+1,j+1,k).eq.b.and.ele(i  ,j  ,k)%form.eq.17 ).or.&
             (nestingdegree(i-1,j+1,k).eq.b.and.ele(i-1,j  ,k)%form.eq.18 ).or.&
             (nestingdegree(i+1,j-1,k).eq.b.and.ele(i  ,j-1,k)%form.eq.18 )) then   
         is_continent(i,j)=.true.   ! crossing ele form 17 or 18
         counter=counter+1
         count_land =  count_land + 1
         nestingdegree(i,j,k)=b
      endif
     else   ! propagate or increment water_basin property on neighboor nodes
      if(nestingdegree(i-1,j,k).eq.b.or.nestingdegree(i+1,j,k).eq.b.or.  &
         nestingdegree(i,j-1,k).eq.b.or.nestingdegree(i,j+1,k).eq.b.or.  &
         nestingdegree(i-1,j,k).eq.d.or.nestingdegree(i+1,j,k).eq.d.or.  &
         nestingdegree(i,j-1,k).eq.d.or.nestingdegree(i,j+1,k).eq.d )then
         counter=counter+1
         count_water =  count_water + 1
         nestingdegree(i,j,k)=d
      elseif((nestingdegree(i-1,j-1,k).eq.d.and.ele(i-1,j-1,k)%form.eq.19).or.&
             (nestingdegree(i+1,j+1,k).eq.d.and.ele(i  ,j  ,k)%form.eq.19).or.&
             (nestingdegree(i-1,j+1,k).eq.d.and.ele(i-1,j  ,k)%form.eq.16).or.&
             (nestingdegree(i+1,j-1,k).eq.d.and.ele(i  ,j-1,k)%form.eq.16)) then   
         counter=counter+1
         count_water =  count_water + 1
         nestingdegree(i,j,k)=d
      endif
     endif
    enddo
   enddo
   EndDo
   if(BndOut)write(*,*)'IDENTIFIED ',count_land+count_water,' nodes among ',vi*uj
   if(BndOut)write(*,*)'Maximal nesting level is ',maxval(nestingdegree(:,:,k))
   if(BndOut)write(*,*)'-------------------------------------------------'
  ENDDO
  !---------------------------------------------------------------




      !*********************************
      !*  FIND ELEMENT TO START FROM   *
      !*********************************
 
  ipos = 0
  jpos = 0
  i=0
  j=0
  found = .false.
  done=.false.

  if(mask_rho(2,2,k) == 1) then    !check to see if the boundary will have to
    if(BndOut)write(*,*)'mask_rho(2,2,k) == 1'
    ipos = 1                    !  pass through the bottom left element
    jpos = 1
    found = .true.
  endif
 

  do i=2,vi-2
    if(found)exit               !if found, exit loop
    do j=2,uj-2
      if( & 
       !((j==2.or.j==uj-2).or.(i==2.or.i==vi-2)) &    ! LOWER,UPPER,LEFT,RIGHT EDGES (BEUG!?)
       !.AND. &
         (ele(i,j,k)%form > 0 .AND. ele(i,j,k)%form < 15) )then
        ipos = i                ! if an element that isn't all water or all
        jpos = j                !  land has been found, store its location
        found = .true.          !  in ipos & jpos and set found to TRUE
        exit                    !  to exit the loops
      endif
    enddo
  enddo

        !*********************************
        !*       CREATE BOUNDARIES       *
        !*********************************
  polydone = .false.                   !Initialize polydone and wf to FALSE
  wf = .false.
  if(found)then
    polydone=.false.
    if(BndOut)write(*,'(a,i3,a,i3,a,i3)')'Create boundary starting in (i,j)=(',&
      ipos,',',jpos,'), ele%form=',ele(ipos,jpos,k)%form
    num_mbnd(k)=num_mbnd(k)+1
  elseif(bnum(k)==0) then
    polydone=.true.
    if(BndOut)write(*,*) 'No boundary exist at this level'
  else
    polydone=.true.
    if(BndOut)write(*,*) 'ERROR, bnum(k)>0 but no boundary found'
  endif

!  write(*,*) 'Creating Boundaries'
!  write(*,*) ' '

  if(ipos == 1 .AND. jpos == 1)then    !If it starts in bottom left corner:
    CALL add(V,ipos+1,jpos,polydone,done,k)  !  Add V(2,1)
    wf = .true.                        !  Initialize wf to TRUE
    dir = 8                            !  Initialize dir to up (8)
    ipos = ipos + 1                    !  And move to element(2,1)
  endif

  do
   !if(polydone)exit         !If water boundaries have been completed, exit
    if( polydone )then
       if(bnum(k)==0.or.BND_SET(k))then
        if(BndOut)write(*,*)'Main bounds Set, ',                               &
           'no looking for another piece of water',bnum(k),BND_SET(k)
        exit         !If water boundaries have been completed, exit
       endif
       ! SEARCH FOR A NEW CONTINENTAL BOUNDARY
        found=.false.
        wf=.false.
        dir=5
        do i=2,vi-2
         if(found)exit               !if found, exit loop
         do j=2,uj-2
           if((ele(i,j,k)%form>0 .and. ele(i,j,k)%form < 15)               &
                            .AND. (ele(i,j,k)%used .EQV. .false.))then
                dmin=minval(nestingdegree(i:i+1,j:j+1,k))
                dmax=maxval(nestingdegree(i:i+1,j:j+1,k))
                if(dmax-dmin.ne.1)then
                   write(*,*)'eleform(',i,',',j,')= ',ele(i,j,k)%form
                   write(*,*)'nestingdegree(i:i+1,j+1,k)=',nestingdegree(i:i+1,j+1,k)
                   write(*,*)'nestingdegree(i:i+1,j,k  )=',nestingdegree(i:i+1,j,k)
                   stop 'ERROR in nestingdegree'
                endif
                ! the following if exclude cases (dmin,dmax)=(1,2):(3,4);(5,6);... corrisponding to island external boundary polygon
                if(mod(dmin,2).eq.0) then ! keep case (dmin,dmax)=(0,1);(2,3);(4,5); ... corrisponding to water basin exteral polygon
                     found=.true.
                     ipos = i
                     jpos = j
                     exit
                endif
            endif
         enddo
        enddo
        if (found) then
          polydone=.false.
          if(BndOut)write(*,*)'FOUND another main piece of water',             &
                             ' starting near element (',ipos,',',jpos,')'
          if(BndOut.and.dmin.eq.2)write(*,*)'nested water basin at i,j=',i,j,&
                                 ' dmin,dmax=',dmin,dmax
          num_mbnd(k)=num_mbnd(k)+1
        else
          if(BndOut)write(*,*)"Bounds NOT Set, but no other is_continental bound"
          exit         ! main water boundaries have been completed, exit
        endif 
    endif
    !write(*,'(a,a,i3,a,i3,a,i3,a,i3)')'                                    ', &
    !    'treat ele(',ipos,',',jpos,')%form=',ele(ipos,jpos,k)%form,' dir=',dir
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11------------------------------------------------------------
    if(ele(ipos,jpos,k)%form > 15 .AND. ele(ipos,jpos,k)%unused)then
      !if the element is a cross and has not been used, set unused to FALSE
      !  this indicates that it has been used once
      ele(ipos,jpos,k)%unused = .false.
    else
      !if the element is NOT a cross, or if it is & has been used once prior
      !  set used to TRUE to indicate that the element will never need to be
      !  used again
      ele(ipos,jpos,k)%used = .true.
    endif


    !NOTE: for more comprehensive explaination of following code see the
    !  subroutine getNext.  The following code is the same as the code
    !  found there, except this code contains calls to the subroutine
    !  add, which adds the boundaries to bnds.  getNext does not add
    !  boundaries so it does not contain calls to add.

    if(wf)then                              !if element on an edge
      if(dir == 2) then                     !right edge
        SELECT CASE(ele(ipos,jpos,k)%form)
          CASE(0,1,4,5)
            CALL add(V,ipos,jpos,polydone,done,k)
            ipos = ipos - 1
            if(ipos == 1)then
              if(.NOT. polydone) CALL add(U,ipos,jpos,polydone,done,k)
              ele(ipos,jpos,k)%used = .true.
              dir = 4
              jpos = jpos - 1
            endif
          CASE DEFAULT
            wf = .false.
        END SELECT
      elseif(dir == 4) then                 !bottom edge
        SELECT CASE(ele(ipos,jpos,k)%form)
          CASE(0,1,2,3)
            CALL add(U,ipos,jpos,polydone,done,k)
            jpos = jpos - 1
            if(jpos == 1)then
              if(.NOT. polydone)CALL add(V,ipos+1,jpos,polydone,done,k)
              ele(ipos,jpos,k)%used = .true.
              dir = 8
              ipos = ipos + 1
            endif
          CASE DEFAULT
          wf = .false.
        END SELECT
      elseif(dir == 6) then                 !top edge
        SELECT CASE(ele(ipos,jpos,k)%form)
          CASE(0,4,8,12)
            CALL add(U,ipos,jpos+1,polydone,done,k)
            jpos = jpos + 1
            if(jpos == uj-1)then
              if(.NOT. polydone) CALL add(V,ipos,jpos,polydone,done,k)
              ele(ipos,jpos,k)%used = .true.
              dir = 2
              ipos = ipos - 1
            endif
          CASE DEFAULT
            wf = .false.
        END SELECT
      elseif(dir == 8) then                 !left edge
        SELECT CASE(ele(ipos,jpos,k)%form)
          CASE(0,2,8,10) 
            CALL add(V,ipos+1,jpos,polydone,done,k)
            ipos = ipos + 1
            if(ipos == vi-1)then
              if(.NOT. polydone)CALL add(U,ipos,jpos+1,polydone,done,k)
              ele(ipos,jpos,k)%used = .true.
              dir = 6
              jpos = jpos + 1
            endif
          CASE DEFAULT
            wf = .false.
        END SELECT
      else
        write(*,*) 'Error: wf direction not one of 2,4,6,8'
      endif


    else                                    !if element NOT on edge

      SELECT CASE(ele(ipos,jpos,k)%form)


        CASE(1,5,13)                        !These forms exit element bottom
          CALL add(V,ipos,jpos,polydone,done,k)
          if(ipos == 2)then
            wf = .true.
            ipos = ipos - 1
            if(.NOT. polydone) CALL add(U,ipos,jpos,polydone,done,k)
            ele(ipos,jpos,k)%used = .true.
            jpos = jpos - 1
            if(jpos == 1)then
              if(.NOT. polydone)CALL add(V,ipos+1,jpos,polydone,done,k)
              ele(ipos,jpos,k)%used = .true.
              ipos = ipos + 1
              dir = 8
            else
              dir = 4
            endif
          else
            ipos = ipos - 1
            dir = 2
          endif


        CASE(2,3,7)                         !These forms exit element left
          CALL add(U,ipos,jpos,polydone,done,k)
          if(jpos == 2)then
            wf = .true.
            jpos = jpos - 1
            if(.NOT. polydone) CALL add(V,ipos+1,jpos,polydone,done,k)
            ele(ipos,jpos,k)%used = .true.
            ipos = ipos + 1
            if(ipos == vi-1)then
              if(.NOT. polydone)CALL add(U,ipos,jpos+1,polydone,done,k)
              ele(ipos,jpos,k)%used = .true.
              jpos = jpos + 1
              dir = 6
            else
              dir = 8
            endif
          else
            jpos = jpos - 1
            dir = 4
          endif


        CASE(4,12,14)                       !These forms exit element right
          CALL add(U,ipos,jpos+1,polydone,done,k)
          if(jpos == uj-2)then
            wf = .true.
            jpos = jpos + 1
            if(.NOT. polydone) CALL add(V,ipos,jpos,polydone,done,k)  ! BEUG FIXED U-> V !
            if(BndOut)write(*,*)'   (DEBEUG added V instead of U',             &
                                    ' at ipos=',ipos,' jpos=',jpos,' )'
            ele(ipos,jpos,k)%used = .true.
            ipos = ipos - 1
            if(ipos == 1)then
              if(.NOT. polydone) CALL add(U,ipos,jpos,polydone,done,k) ! BEUG GUESSED, modif done: V-> U !
              if(BndOut)write(*,*)'   (DEBEUG added U instead of V',           &
                                    ' at ipos=',ipos,' jpos=',jpos,' )'
              ele(ipos,jpos,k)%used = .true.
              jpos = jpos - 1
              dir = 4
            else
              dir = 2
            endif
          else
            jpos = jpos + 1
            dir = 6
          endif


        CASE(8,10,11)                       !These forms exit element top
          CALL add(V,ipos+1,jpos,polydone,done,k)
          if(ipos == vi-2)then
            wf = .true.
            ipos = ipos + 1
            if(.NOT. polydone) CALL add(U,ipos,jpos+1,polydone,done,k)
            ele(ipos,jpos,k)%used = .true.
            jpos = jpos + 1
            if(jpos == uj-1)then
              if(.NOT. polydone) CALL add(V,ipos,jpos,polydone,done,k)
              ele(ipos,jpos,k)%used = .true.
              ipos = ipos - 1
              dir = 2
            else
              dir = 6
            endif
          else
            ipos = ipos + 1
            dir = 8
          endif


        CASE(16,17)                         !These forms exit right or left
                                            !16: u->r & d->l  17: d->r & u->l

          if((ele(ipos,jpos,k)%form == 16 .AND. dir == 2).OR.   & !if exit right
             (ele(ipos,jpos,k)%form == 17 .AND. dir == 8)) then

            CALL add(U,ipos,jpos+1,polydone,done,k)
            if(jpos == uj-2)then
              wf = .true.
              jpos = jpos + 1
              if(.NOT. polydone) CALL add(U,ipos,jpos,polydone,done,k)
              ele(ipos,jpos,k)%used = .true.
              ipos = ipos - 1
              if(ipos == 1)then
                if(.NOT. polydone)CALL add(V,ipos,jpos,polydone,done,k)
                ele(ipos,jpos,k)%used = .true.
                jpos = jpos - 1
                dir = 4
              else
                dir = 2
              endif
            else
              jpos = jpos + 1
              dir = 6
            endif
          
          else                                                  !if exit left

            CALL add(U,ipos,jpos,polydone,done,k)
            if(jpos == 2)then
              wf = .true.
              jpos = jpos - 1
              if(.NOT.polydone)CALL add(V,ipos+1,jpos,polydone,done,k)
              ele(ipos,jpos,k)%used = .true.
              ipos = ipos + 1
              if(ipos == vi-1)then
                if(.NOT. polydone) CALL add(U,ipos,jpos+1,polydone,done,k)
                ele(ipos,jpos,k)%used = .true.
                jpos = jpos + 1
                dir = 6
              else
                dir = 8
              endif
            else
              jpos = jpos - 1
              dir = 4
            endif
          
          endif


        CASE(18,19)                         !These forms exit up or down
                                            !18: r->u & l->d  19: l->u & r->d

          if((ele(ipos,jpos,k)%form == 18 .AND. dir == 4).OR.   & !if exit top
             (ele(ipos,jpos,k)%form == 19 .AND. dir == 6)) then

            CALL add(V,ipos+1,jpos,polydone,done,k)
            if(ipos == vi-2)then
              wf = .true.
              ipos = ipos + 1
              if(.NOT.polydone)CALL add(U,ipos,jpos+1,polydone,done,k)
              ele(ipos,jpos,k)%used = .true.
              jpos = jpos + 1
              if(jpos == uj-1)then
                if(.NOT. polydone)CALL add(V,ipos,jpos,polydone,done,k)
                ele(ipos,jpos,k)%used = .true.
                ipos = ipos - 1
                dir = 2
              else
                dir = 6
              endif
            else
              ipos = ipos + 1
              dir = 8
            endif

          else                                                  !if exit down

            CALL add(V,ipos,jpos,polydone,done,k)
            if(ipos == 2)then
              wf = .true.
              ipos = ipos - 1
              if(.NOT. polydone) CALL add(U,ipos,jpos,polydone,done,k)
              ele(ipos,jpos,k)%used = .true.
              jpos = jpos - 1
              if(jpos == 1)then
                if(.NOT. polydone) CALL add(V,ipos+1,jpos,polydone,done,k)
                ele(ipos,jpos,k)%used = .true.
                ipos = ipos + 1
                dir = 8
              else
                dir = 4
              endif
            else
              ipos = ipos - 1
              dir = 2
            endif

          endif

      END SELECT
    endif

  enddo


   if(BndOut)write(*,*) ' '
   if(BndOut)write(*,*) '.........'
   if(BndOut)write(*,*) ''
   if(BndOut)write(*,*) 'ISLANDS '
 
        !*********************************
        !*         ISLAND SECTION        *
        !*********************************
  if(BND_SET(k).or.bnum(k)==0) then
    done=.true.
    if(BndOut)write(*,*) 'No boundary points left so no islands can exist',    &
                          ' at this level'
  else
    if(BndOut)write(*,*) 'looking for islands'
  endif

  do
    if(done)exit                  !if not all boundary points have been used
                                  !  Islands must exist...


        !*********************************
        !*  FIND ELEMENT TO START ISLAND *
        !*********************************

    found = .false.               !initialize found to false

    do i=2,vi-2
      do j=2,uj-2
        if(ele(i,j,k)%form < 15 .AND. (ele(i,j,k)%used .EQV. .false.))then
          ipos = i                !find an element that hasn't been used yet
          jpos = j                !store its location in ipos,jpos
          found = .true.
          exit
        endif
      enddo
      if(found)exit               !if one has been found, exit
    enddo

   if(BndOut)write(*,'(a,i3,a,i3,a,i3)')                                       &
        'ISLAND : found in (',ipos,',',jpos,')%form=',ele(ipos,jpos,k)%form
        !*********************************
        !*    CREATE ISLAND BOUNDARIES   *
        !*********************************


    !initialize polydone to FALSE
    if(found)then
      polydone=.false.
    else
      polydone=.true.
      if(BndOut)write(*,*) 'No more Island boundary found at this level'
    endif

    do

      if(polydone)exit            !if the current island is finished, exit


      if(ele(ipos,jpos,k)%form > 15 .AND. ele(ipos,jpos,k)%unused)then
        !if the element is a cross and has not been used, set unused to FALSE
        !  this indicates that it has been used once
        ele(ipos,jpos,k)%unused = .false.
      else
        !if the element is NOT a cross, or if it is & has been used once prior
        !  set used to TRUE to indicate that the element will never need to be
        !  used again
        ele(ipos,jpos,k)%used = .true.
      endif

      !NOTE: there is no wf section this time because islands cannot be on
      !  the edges of the grid... land on the edges of the grid would be
      !  already used in the previous section

      SELECT CASE(ele(ipos,jpos,k)%form)

        CASE(1,3,11)              !These forms exit element right
          CALL add(U,ipos,jpos+1,polydone,done,k)
          jpos = jpos + 1
          dir = 6

        CASE(2,10,14)             !These forms exit element bottom
          CALL add(V,ipos,jpos,polydone,done,k)
          ipos = ipos - 1
          dir = 2

        CASE(4,5,7)               !These forms exit element top
          CALL add(V,ipos+1,jpos,polydone,done,k)
          ipos = ipos + 1
          dir = 8

        CASE(8,12,13)             !These forms exit element left
          CALL add(U,ipos,jpos,polydone,done,k)
          jpos = jpos - 1
          dir = 4

        CASE(16,17)               !These forms exit up or down
                                  !16: r->u & l->d    17: l->u & r->d

          if((ele(ipos,jpos,k)%form == 16 .AND. dir == 4) .OR.  & !if exit up
             (ele(ipos,jpos,k)%form == 17 .AND. dir == 6)) then
            CALL add(V,ipos+1,jpos,polydone,done,k)
            ipos = ipos + 1
            dir = 8
          else                                                  !if exit down
            CALL add(V,ipos,jpos,polydone,done,k)
            ipos = ipos - 1
            dir = 2
          endif

        CASE(18,19)               !These forms exit right or left
                                  !18: u->r & d->l    19: d->r & u->l

          if((ele(ipos,jpos,k)%form == 18 .AND. dir == 2) .OR.  & !if exit right
             (ele(ipos,jpos,k)%form == 19 .AND. dir == 8)) then
            CALL add(U,ipos,jpos+1,polydone,done,k)
            jpos = jpos + 1
            dir = 6           
          else                                                  !if exit left
            CALL add(U,ipos,jpos,polydone,done,k)
            jpos = jpos - 1
            dir = 4           
          endif

      END SELECT

    enddo

  enddo ! Island detection
  if(BndOut)write(pyout+k,'(a,i3,a)')"	print (str(",k,")+' done')    #PYTHON"
  if(BndOut)write(*,'(a,i3,a)')"------------------------- DONE with level",k, &
            "--------------------------------"
  ENDDO ! k=1,us_tridim

  if(BndOut)call createNestingDegreeNetCDF(nestingdegree)


!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   ~              Communicate eleform to HYDRO module                      ~
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   call seteleform(vi-1,uj-1,us_tridim,ele(:,:,:)%form)

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   ~               Prepare Boundary and Island Coordinates                 ~
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Set up ROMS boundaries. Use location of u and v grid points to demarcate the 
! land boundary between rho points where mask_rho = 0 and mask_rho = 1.

  !get the x/y locations of the U and V grid nodes
  CALL getUVxy(x_u,y_u,x_v,y_v)

  if(BoundaryBLNs) then
    OPEN(10,FILE='OpenOceanBoundaryMidpoints.csv',STATUS='REPLACE')
    CLOSE(10)
    OPEN(11,FILE='LandBoundaryMidpoints.csv',STATUS='REPLACE')
    CLOSE(11)
  endif

  ALLOCATE(tot_mbnd_pts(us_tridim),STAT=STATUS)
  IF(STATUS /= 0) write(*,*) 'Problem allocate tot_mbnd_pts'
  ALLOCATE(tot_ibnd_pts(us_tridim),STAT=STATUS)
  IF(STATUS /= 0) write(*,*) 'Problem allocate tot_ibnd_pts'
  ALLOCATE(num_ibnd(us_tridim),STAT=STATUS)
  IF(STATUS /= 0) write(*,*) 'Problem allocate num_ibnd'
  ALLOCATE(nbounds(us_tridim),STAT=STATUS)
  IF(STATUS /= 0) write(*,*) 'Problem allocate nbounds'
  ALLOCATE(numpoly(us_tridim),STAT=STATUS)
  IF(STATUS /= 0) write(*,*) 'Problem allocate numpoly'

  max_mbnd_pts=0
  max_ibnd_pts=0
  nbounds_max=0
  numpoly_max=0

  ! get numpoly_max to allocate polysizes
  do klevel=1,us_tridim
    if(bnum(klevel)>0)then
      numpoly(klevel) = bnds(bnum(klevel),klevel)%poly       ! # of polygons
      if (numpoly(klevel)>numpoly_max) numpoly_max=numpoly(klevel)
    else
      numpoly(klevel) = 0
    endif
    
    if(BndOut)write(*,'(a,i2,a,i5,a,i2,a,i5,a,i2,a,i5)')'bnum(',klevel,')=',   &
                bnum(klevel),', numpoly(',klevel,')=',numpoly(klevel),         &
                ', num_mbnd(',klevel,')=',num_mbnd(klevel)
  enddo

  !allocate polysizes to contain the # of boundary points in each polygon for each k-level
  if(BndOut)write(*,*)'allocate polysizes (',numpoly_max,us_tridim,')'
  ALLOCATE(polysizes(numpoly_max,us_tridim),STAT=STATUS)
  IF(STATUS /= 0) write(*,*) 'Problem allocate polysizes'
  polysizes = 0    !initialize to 0

  !allocate polynestdeg to contain the degree of nesting of each polygon for every k-level
  if(BndOut)write(*,*)'allocate polynestdeg (',numpoly_max,us_tridim,')'
  ALLOCATE(polynestdeg(numpoly_max,us_tridim),STAT=STATUS)
  IF(STATUS /= 0) write(*,*) 'Problem allocate polynestdeg'
  polynestdeg = -1
  DO klevel=1,us_tridim
   if(numpoly(klevel)>0)then
    do bnd=1,bnum(klevel)
      i = bnds(bnd,klevel)%ii
      j = bnds(bnd,klevel)%jj
      if( bnds(bnd,klevel)%onU ) then
        polynestdeg(bnds(bnd,klevel)%poly,klevel)= &
         max(nestingdegree(i,j,klevel),nestingdegree(i+1,j,klevel))
      else
        polynestdeg(bnds(bnd,klevel)%poly,klevel)= &
         max(nestingdegree(i,j,klevel),nestingdegree(i,j+1,klevel))
      endif
    enddo
   endif
  ENDDO

!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   ~              Create Boundary Blanking Files for Surfer                ~
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  if(BoundaryBLNs) then
   CALL output_llBounds()
   CALL output_xyBounds()
  endif
!   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


 
  tot_mbnd_pts=0
  ! get max_mbnd_pts to allocate bx,by and oo
  DO klevel=1,us_tridim
  !iterate through boundary points, counting the # in each polygon
  if(numpoly(klevel)>0)then

    do bnd=1,bnum(klevel)
      polysizes(bnds(bnd,klevel)%poly,klevel) =                                &
                 polysizes(bnds(bnd,klevel)%poly,klevel) + 1
    enddo
    !Sum up the total # of main boundary points
    do m=1,num_mbnd(klevel)
     tot_mbnd_pts(klevel)=tot_mbnd_pts(klevel)+polysizes(m,klevel) + 1
     if( tot_mbnd_pts(klevel)>max_mbnd_pts) max_mbnd_pts=tot_mbnd_pts(klevel)
    enddo
    if(BndOut)then
      if(numpoly(klevel)<10)write(formatoutput,'(a,i1,a)')                     &
                                 '(a,i2,a,i3,a,i2,a,',numpoly(klevel),'(i6))'
      if(numpoly(klevel)>9 .and. numpoly(klevel)<100 )                         &
     write(formatoutput,'(a,i2,a)')'(a,i2,a,i3,a,i2,a,',numpoly(klevel),'(i6))'
      if(numpoly(klevel)>99 .and. numpoly(klevel)<1000 )                       &
     write(formatoutput,'(a,i3,a)')'(a,i2,a,i3,a,i2,a,',numpoly(klevel),'(i6))'
      if(numpoly(klevel)>999 .and. numpoly(klevel)<10000 )                     &
     write(formatoutput,'(a,i4,a)')'(a,i2,a,i3,a,i2,a,',numpoly(klevel),'(i6))'
      !write(*,*)formatoutput
      write(*,formatoutput)'klevel=',klevel,                                 &
         ' polysizes(1:',numpoly(klevel),',',klevel,')=',                        &
         polysizes(1:numpoly(klevel),klevel)
    endif!BndOut
  endif
  ENDDO ! klevel=1,us_tridim

  if(BndOut)write(*,*)'max_mbnd_pts=',max_mbnd_pts
  ALLOCATE(mid(max_mbnd_pts,us_tridim),STAT=STATUS)         !   #
  IF(STATUS /= 0) write(*,*) 'Problem allocate mid'
  ALLOCATE(bx(max_mbnd_pts,us_tridim),STAT=STATUS)         !  and allocate bx and by to this #
  IF(STATUS /= 0) write(*,*) 'Problem allocate bx'
  ALLOCATE(by(max_mbnd_pts,us_tridim),STAT=STATUS)         !the +1 is to close the polygon
  IF(STATUS /= 0) write(*,*) 'Problem allocate by'
  ALLOCATE(oo(max_mbnd_pts,us_tridim),STAT=STATUS)         !Track open ocean boundary points
  IF(STATUS /= 0) write(*,*) 'Problem allocate oo'
  !initialize oo
  oo = .FALSE.

  ! get numisland_max to allocate hid, hx, hy
  DO klevel=1,us_tridim
    num_ibnd(klevel) = numpoly(klevel) - num_mbnd(klevel)        ! # of islands
    !initialize # of island boundaries to the # of islands
    !  this takes care of +1 for each island (see comment above about +1)
    tot_ibnd_pts(klevel) = num_ibnd(klevel)         
    if(num_ibnd(klevel)>0)then
      !Sum up the total # of island boundary points
      do isl=1,num_ibnd(klevel)
        tot_ibnd_pts(klevel) = tot_ibnd_pts(klevel) +                          &
                               polysizes(num_mbnd(klevel)+isl,klevel)
      enddo
      if (tot_ibnd_pts(klevel)>max_ibnd_pts) max_ibnd_pts=tot_ibnd_pts(klevel)
    endif
  ENDDO 
 
  if(BndOut)write(*,*)'max_ibnd_pts=',max_ibnd_pts

  !  and allocate hid,hx,hy to the # of island boundary points
  ALLOCATE(hid(max_ibnd_pts,us_tridim),STAT=STATUS)  
  IF(STATUS /= 0) write(*,*) 'Problem allocating hid'
  ALLOCATE(hx(max_ibnd_pts,us_tridim),STAT=STATUS)  
  IF(STATUS /= 0) write(*,*) 'Problem allocating hx'
  ALLOCATE(hy(max_ibnd_pts,us_tridim),STAT=STATUS)  
  IF(STATUS /= 0) write(*,*) 'Problem allocating hy'
  !
  DO klevel=1,us_tridim
  if(BndOut)write(*,'(3(a,i3))')'k=',klevel,', tot_mbnd_pts=',                 &
            tot_mbnd_pts(klevel),', tot_ibnd_pts(klevel)=',tot_ibnd_pts(klevel)

  pstart = 1

  !iterate through the polygons adding the boundary points to bx,by,hid,hx,hy
  if(BndOut)write(*,*)'iterate through the ',numpoly(klevel),                  & 
              ' polygons of level',klevel
  do m=1,numpoly(klevel)
    pend = pstart + polysizes(m,klevel) - 1
    if( (m==1) .or. (m==(num_mbnd(klevel)+1)) ) c = 0
    if(BndOut)write(*,'(5(a,i4))')'       poly number',m,', c=[',c+1,' :',     &
                                c+1+pend-pstart+1,                             & 
                                '], pstart=',pstart,', pend=',pend
    do bnd=pstart,pend
      c = c + 1
      i = bnds(bnd,klevel)%ii
      j = bnds(bnd,klevel)%jj
      if(m .le. num_mbnd(klevel)) then             !if it is the main body of water boundary
        mid(c,klevel) = m           !  add main boundary to mid

        if( bnds(bnd,klevel)%onU ) then    !  if its on the U grid
          bx(c,klevel) = x_u(i,j)        !  add the U grid node location to bx and by
          by(c,klevel) = y_u(i,j)
          !If this boundary point is open ocean, set oo for this point to TRUE
          if(    ( (i==1   ) .AND. mask_rho( 1,j,klevel)==1) &
          .OR.((i==vi-1).AND.mask_rho(vi,j,klevel)==1) ) oo(c,klevel) = .TRUE.
        else                      !  if its on the V grid
          bx(c,klevel) = x_v(i,j)        !  add the V grid node location to bx and by
          by(c,klevel) = y_v(i,j)
          !If this boundary point is open ocean, set oo for this point to TRUE
          if(    ( (j==1   ) .AND. mask_rho(i, 1,klevel)==1) &
            .OR. ((j==uj-1).AND.mask_rho(i,uj,klevel)==1))oo(c,klevel) = .TRUE.
        endif
      else                        !if it is an island boundary
        hid(c,klevel) = m-num_mbnd(klevel)+1+999           !  add island # to hid
                                  !  note: island IDs start at 1001
        if( bnds(bnd,klevel)%onU ) then    !  if its on the U grid
          hx(c,klevel) = x_u(i,j)        !  add the U grid node location to hx and hy
          hy(c,klevel) = y_u(i,j)
        else                      !  if its on the V grid
          hx(c,klevel) = x_v(i,j)        !  add the V grid node location to hx and hy
          hy(c,klevel) = y_v(i,j)
        endif
      endif
    enddo ! bnd=pstart,pend

    !At the end of each polygon, put the first boundary point on the end to 
    !  close the polygon
    c = c + 1
    bnd = pstart
    i = bnds(bnd,klevel)%ii
    j = bnds(bnd,klevel)%jj
    if(m.le.num_mbnd(klevel))then                !if its mainbay bound
      mid(c,klevel) = m           !  add main boundary to mid
      if( bnds(bnd,klevel)%onU ) then
          bx(c,klevel) = x_u(i,j)
          by(c,klevel) = y_u(i,j)
          !If this boundary point is open ocean, set oo for this point to TRUE
          if(    ( (i==1   ) .AND. mask_rho( 1,j,klevel)==1) &
            .OR. ( (i==vi-1) .AND. mask_rho(vi,j,klevel)==1) )                 &
                                     oo(c,klevel) = .TRUE.
      else
          bx(c,klevel) = x_v(i,j)
          by(c,klevel) = y_v(i,j)
          !If this boundary point is open ocean, set oo for this point to TRUE
          if(    ( (j==1   ) .AND. mask_rho(i, 1,klevel)==1) &
            .OR. ( (j==uj-1) .AND. mask_rho(i,uj,klevel)==1) )                 &
                                     oo(c,klevel) = .TRUE.
      endif
    else                          !if its island bound
      if( bnds(bnd,klevel)%onU ) then
          hid(c,klevel) = m-num_mbnd(klevel)+1+999
          hx(c,klevel) = x_u(i,j)
          hy(c,klevel) = y_u(i,j)
      else
          hid(c,klevel) = m-num_mbnd(klevel)+1+999
          hx(c,klevel) = x_v(i,j)
          hy(c,klevel) = y_v(i,j)
      endif
    endif

    pstart = pstart + polysizes(m,klevel)

  enddo ! m=1,numpoly(klevel)


  ! Construct bnd_x and bnd_y matrices that contain boundary line segments.
  ! bnd_x1(i),bnd_y1(i) = coordinates of first point in line segment.
  ! bnd_x2(i),bnd_y2(i) = coordinates of second point in line segment.
  ! The bnd_x/bnd_y matrices are used to 
  !   1) determine if a particle intersects a boundary segment and 
  !   2) reflect the particle off the boundary segment.
  !write(*,*)'tot_ibnd_pts(klevel)=',tot_ibnd_pts(klevel)
  !nbounds(klevel)= tot_mbnd_pts(klevel) - num_mbnd(klevel) -1 & ! beug nbound too small
  !              + tot_ibnd_pts(klevel) - num_ibnd(klevel) - 1
  nbounds(klevel)= tot_mbnd_pts(klevel) - num_mbnd(klevel)     & !CL-OGS beug solved
                 + tot_ibnd_pts(klevel) - num_ibnd(klevel)    
  if(nbounds(klevel)> nbounds_max) nbounds_max=nbounds(klevel)
  ENDDO  ! klevel=1,us_tridim

  ALLOCATE(bnd_x1(nbounds_max,us_tridim),STAT=STATUS)  
  IF(STATUS /= 0) write(*,*) 'Problem allocating bnd_x1'
  ALLOCATE(bnd_x2(nbounds_max,us_tridim),STAT=STATUS)  
  IF(STATUS /= 0) write(*,*) 'Problem allocating bnd_x2'
  ALLOCATE(bnd_y1(nbounds_max,us_tridim),STAT=STATUS)  
  IF(STATUS /= 0) write(*,*) 'Problem allocating bnd_y1'
  ALLOCATE(bnd_y2(nbounds_max,us_tridim),STAT=STATUS)  
  IF(STATUS /= 0) write(*,*) 'Problem allocating bnd_y2'

  !Allocate land to track land/ocean boundaries
  ALLOCATE (land(nbounds_max,us_tridim),STAT=STATUS)
  if(STATUS /= 0) write(*,*) 'Problem allocating land'

  !Initialize newly allocated land
  !TRUE if the boundary is land, FALSE if it is open ocean
  land = .TRUE.

  if(BoundaryBLNs) then
    OPEN(11,FILE='LandBoundaryMidpoints.csv')
    OPEN(10,FILE='OpenOceanBoundaryMidpoints.csv')
  endif

  ! First create the line segments for the land and open ocean boundaries by
  !   simply writing each x/y location and the next x/y location. The fact 
  !   that the first bx/by coordinate equals the last ensures that all line 
  !   segments surrounding the model are included.
  do klevel=1,us_tridim
    if(BndOut)write(*,*)' ----------------------------------------'
    if(BndOut)write(*,'(6(a,i4))')'k=',klevel,', tot_mbnd_pts=',               &
                           tot_mbnd_pts(klevel),' nbounds_max=',               &
                           nbounds_max,' num_mbnd=',num_mbnd(klevel), &
   ', tot_ibnd_pts(klevel)=',tot_ibnd_pts(klevel),' num_ibnd=',num_ibnd(klevel)
 
    if(BndOut)write(*,*)' '

    j=0  !NOT BEUG : j=0 should be right and should not be j=1 
    k=mid(1,klevel)
    nmp=1
    io=1
    do i=1,tot_mbnd_pts(klevel)-1
        if(mid(i,klevel) .NE. k)then       !if moving on to a new main bound
          if(BndOut)write(*,'(10(a,i5))')'k=',klevel,' polisize=',             &
                        polysizes(nmp,klevel),' mid=',mid(i-1,klevel),' j=',j, &
                                      '      segments i-j=',io-j,' :',(i-1)-j, &
                                         '      madebynodes i=',io,' :',(i-1), &
                                        '      and i+1=',io+1,' :',(i-1)+1
          nmp=nmp+1                    
          io=i                         
          k = mid(i,klevel)                !update k to new main bound id
          j = j + 1                 !and increment j
        endif
      bnd_x1(i-j,klevel) = bx(i  ,klevel)   
      bnd_y1(i-j,klevel) = by(i  ,klevel)
      bnd_x2(i-j,klevel) = bx(i+1,klevel)
      bnd_y2(i-j,klevel) = by(i+1,klevel)
      if(oo(i,klevel).AND.oo(i+1,klevel))then
        land(i,klevel) = .FALSE. !Not a land boundary
        if(BoundaryBLNs) WRITE(10,*) bx(i,klevel),',',bx(i+1,klevel),',',      &
                            by(i,klevel),',',by(i+1,klevel),',',klevel,',',i-j
      endif
      if(BoundaryBLNs) WRITE(11,*) bx(i,klevel),',',bx(i+1,klevel),',',        &
                      by(i,klevel),',',by(i+1,klevel),',',klevel,',',i-j,',',0
     
    enddo
   !if(BndOut)write(*,*)'line segment stored'
   !if(BndOut)write(*,*)'klevel=',klevel,'nmp=',nmp
   !if(BndOut)write(*,*)' polysize=',polysizes
   !if(BndOut)write(*,*)' mid=',mid     
   !     !(i-1,klevel),' at',i-1,klevel
   !if(BndOut)write(*,*)'line segment stored'
   !if(BndOut)write(*,*)'line segment stored'
   !if(BndOut)write(*,*)'line segment stored'
   !if(BndOut)write(*,*)'line segment stored'
   !if(BndOut)write(*,'(a,i5)')' j=',j
   !if(BndOut)write(*,'(a,i5)')' segments i-j=',io-j
   !if(BndOut)write(*,'(a,i5)')' :',(i-1)-j
   !if(BndOut)write(*,'(a,i5)')' madebynodes i=',io
   !if(BndOut)write(*,'(a,i5)')' :',(i-1)
   !if(BndOut)write(*,'(a,i5)')' and i+1=',io+1
   !if(BndOut)write(*,'(a,i5)')' :',(i-1)+1
    if(num_ibnd(klevel)>0)then         
      ! Now create the line segments for the island boundaries. Same principle as
      !   above, but need to account for the fact that there are islands in the
      !   hx/hy matrices. Use hid matrix of island id numbers to guide process. 
      j=0                         
      k=hid(1,klevel)
      io=1
         
      do i=1,tot_ibnd_pts(klevel)-1
        if(hid(i,klevel) .NE. k)then       !if moving on to a new island
         ! if(BndOut)write(*,'(10(a,i5))')'k=',klevel,' polisize=',polysizes(nmp,klevel),' hid=',hid(i-1,klevel),' j=',j, &
         !                       '      segments i-j=', tot_mbnd_pts(klevel)-num_mbnd(klevel)+io-j,' :', & 
         !                                              tot_mbnd_pts(klevel)-num_mbnd(klevel)+(i-1)-j, &
         !                       '      madebynodes i=',io,' :', & 
         !                                              (i-1), &
         !                       '      and i+1=',      io+1,' :', & 
         !                                         (i-1)+1
          io=i
          nmp=nmp+1
          k = hid(i,klevel)                !update k to new island id
          j = j + 1                 !and increment j
          if(BndOut)write(*,*)'imax=',tot_ibnd_pts(klevel)-1, &
         ', tot_mbnd_pts=',tot_mbnd_pts(klevel),' num_mbnd=',num_mbnd(klevel), &
         ' j=',j
        endif
        bnd_x1(tot_mbnd_pts(klevel)-num_mbnd(klevel)+i-j,klevel)=hx(i  ,klevel)
        bnd_y1(tot_mbnd_pts(klevel)-num_mbnd(klevel)+i-j,klevel)=hy(i  ,klevel)
        bnd_x2(tot_mbnd_pts(klevel)-num_mbnd(klevel)+i-j,klevel)=hx(i+1,klevel)
        bnd_y2(tot_mbnd_pts(klevel)-num_mbnd(klevel)+i-j,klevel)=hy(i+1,klevel)
        if(BoundaryBLNs) WRITE(11,*) hx(i,klevel),',',hx(i+1,klevel),',',      &
                                     hy(i,klevel),',',hy(i+1,klevel),',',      &
              klevel,',',tot_mbnd_pts(klevel)-num_mbnd(klevel)+i-j,',',1
      enddo
      if(BndOut)write(*,'(10(a,i5))')'k=',klevel,' polisize=',                 &
                        polysizes(nmp,klevel),' hid=',hid(i-1,klevel),' j=',j, &
            ' segments i-j=', tot_mbnd_pts(klevel)-num_mbnd(klevel)+io-j,' :', & 
                               tot_mbnd_pts(klevel)-num_mbnd(klevel)+(i-1)-j,  &
             '      madebynodes i=',io,' :',                                   & 
                               (i-1),                                          &
             '      and i+1=',      io+1,' :',                                 & 
                               (i-1)+1
    else
     if(BndOut)write(*,*)'No island boundaries'
    endif
  enddo
  if(BndOut)write(*,*)'now set boundaries calling setbounds'

  call setbounds(nbounds_max,us_tridim,nbounds,bnd_x1,bnd_x2,bnd_y1,bnd_y2)

  if(BoundaryBLNs) CLOSE(10)
  if(BoundaryBLNs) CLOSE(11)
  
  if(BndOut)write(*,*)'End Subroutine createBounds'

  DEALLOCATE(bnds,polysizes)
  DEALLOCATE(x_u)
  DEALLOCATE(y_u)
  DEALLOCATE(x_v)
  DEALLOCATE(y_v)
  DEALLOCATE(mask_rho)
  DEALLOCATE(oo)
  DEALLOCATE(is_continent)
  DEALLOCATE(nestingdegree)
  DEALLOCATE(numpoly)
  DEALLOCATE(crossnum)
  
  DO k=1,us_tridim 
    close(pyout+k)
  enddo

END SUBROUTINE createBounds



    !*********************************************************
    !*                         add                           *
    !*********************************************************

SUBROUTINE add(isU,iii,jjj,polydone,done,klev)
!This subroutine is for adding boundary points to bnds
! INPUT:
!   isU - TRUE if point is on U grid, FALSE if it is on V grid
!   iii - i position on the grid
!   jjj - j position on the grid

  USE PARAM_MOD, ONLY:BndOut
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iii,jjj,klev
  LOGICAL, INTENT(IN) :: isU
  LOGICAL, INTENT(INOUT) :: polydone,done
  

  INTEGER :: polystart = 1  !bnd location at start of current polygon
  INTEGER :: polynum   = 1  !# of the current polygon
  INTEGER :: b = 0          !# of boundary points used so far
  !Note that any local variable that is initialized in a type declaration
  !  statement, as done above, is automatically given the SAVE attribute,
  !  this means that the value in these variables will be retained in 
  !  subsequent calls to the subroutine

  !Check if the current polygon has returned to its starting location
! CRS commented out next two lines
!  if( (polystart /= (b+1)) .AND. (bnds(polystart)%onU == isU) .AND.           &  
!      (bnds(polystart)%ii == iii).AND.(bnds(polystart)%jj == jjj) )then
! CRS replaced them with:
  if( (polystart .NE. (b+1)) .AND. (bnds(polystart,klev)%onU .EQV. isU) .AND.  &
      (bnds(polystart,klev)%ii==iii).AND.(bnds(polystart,klev)%jj==jjj))then

    !switch polydone to TRUE to indicate the current polygon is complete
    polydone = .true.

    !since the current polygon is complete, increment polynum 
    polynum = polynum + 1
    !and set polystart to the first position of the subsequent polygon
    polystart = b + 1

    if(BndOut)write(*,'(a,i6,a,i7,a,i7)')' Polynum ',polynum-1,                &
    ' returned to its started position after point ',b,'/',bnum(klev)
    if(BndOut)write(*,*)' '

    if(b == bnum(klev))then
      done = .true.
      BND_SET(klev) = .true.
      polystart = 1
      polynum = 1
      b = 0
      if(BndOut)write(*,*)'Boundary parameters reinitialised.'
    endif

  else
    if(mod(b-1,bnum_max)==0 .and.b>bnum_max) write(*,*)'bndserror after:: ',   &
                    bnum_max,us_tridim,b,klev
    !If the current polygon isn't complete
    b = b+1             !increment b to the next boundary point and
    bnds(b,klev)%onU  = isU  !store values for onU, ii, jj, and polynum
    bnds(b,klev)%ii   = iii
    bnds(b,klev)%jj   = jjj
    bnds(b,klev)%poly = polynum
    if (isU) then
      if(BndOut)write(pyout+klev,'(a,a,i4,a,i4,a,i2,a,l2,i4,i4,i10,a,i10)')     &
                           "	if (v=='mask_u'",                              &
                           "): datafield[",jjj-1,",",iii-1,                    &
                           "]= 0.7 # PYTHON ele=",ele(iii,jjj,klev)%form,'  ', &
                           wf,dir,polynum,b,'/',bnum(klev)
    else
      if(BndOut)write(pyout+klev,'(a,a,i4,a,i4,a,i2,a,l2,i4,i4,i10,a,i10)')     &
                                   "	if (v=='mask_v'",                      &
                                   "): datafield[",jjj-1,",",iii-1,            &
                           "]= 0.7 # PYTHON ele=",ele(iii,jjj,klev)%form,'  ', &
                                   wf,dir,polynum,b,'/',bnum(klev)
    endif
    if (b>2) then
      if ( ((bnds(b,klev)%onU  .eqv. bnds(b-1,klev)%onU)   .AND.               &
            (bnds(b,klev)%ii     ==  bnds(b-1,klev)%ii))   .AND.               &
           ((bnds(b,klev)%jj     ==  bnds(b-1,klev)%jj)    .AND.               &
            (bnds(b,klev)%poly   ==  bnds(b-1,klev)%poly)) ) then 
          WRITE(*,*) 'BEUG ERROR? 2 identical bnds points ',                   &
                         'identified in add() routine!!!'
      endif
    endif

  endif

END SUBROUTINE add



    !*********************************************************
    !*                       getNext                         *
    !*********************************************************

SUBROUTINE getNext(i,j,wf,dir,form)
  USE PARAM_MOD, ONLY: vi,uj
  !This subroutine is for finding the next element, following the bounds
  !  clockwise around water
  !This is all based on the element location, whether or not its on an 
  !  edge, the direction the path is coming from, and the element form
  IMPLICIT NONE

  INTEGER :: i,j,dir,form
  LOGICAL :: wf

  if(wf)then                      !If the element is on the edge
    if(dir == 2) then             !If path is going down,
      SELECT CASE(form)           !  then its on right edge

        CASE(0,1,4,5)             !If path stays on edge
          i = i - 1
          if(i == 1)then          !If path reaches bottom, #BEUG CL changed 2to1
            dir = 4               !  then go left on bottom edge
            j = j - 1
          endif

        CASE DEFAULT              !Path exits edge
          wf = .false.

      END SELECT
    elseif(dir == 4) then         !If path is going left,
      SELECT CASE(form)           !  then its on bottom edge

        CASE(0,1,2,3)             !If path stays on edge
          j = j - 1
          if(j == 1)then          !If path reaches left edge, #BEUG CL changed 2to1
            dir = 8               !  then go up on left edge
            i = i + 1
          endif

        CASE DEFAULT              !Path exits edge
          wf = .false.

      END SELECT
    elseif(dir == 6) then         !If path is going right,
      SELECT CASE(form)           !  then its on top edge

        CASE(0,4,8,12)            !If path stays on edge
          j = j + 1
          ! BEUG CL: next line uj-2 should be changed to uj-1
          if(j == uj-2)then       !If path reaches right edge,
            dir = 2               !  then go down on right edge
            i = i - 1
          endif

        CASE DEFAULT              !Path exits edge
          wf = .false.

      END SELECT
    elseif(dir == 8) then         !If path is going up,
      SELECT CASE(form)           !  then its on left edge

        CASE(0,2,8,10)            !If path stays on edge
          i = i + 1
          ! BEUG CL: next line vi-2 should be changed to vi-1
          if(i == vi-2)then       !If path reaches top edge,
            dir = 6               !  then go right on top edge
            j = j + 1
          endif

        CASE DEFAULT              !Path exits edge
          wf = .false.

      END SELECT
    else
      write(*,*) 'Error: wf direction not one of 2,4,6,8'
    endif





  else                            !If element is not on edge
    SELECT CASE(form)

    CASE(1,5,13)                  !These forms exit element bottom
      if(i == 2)then              !  if moving into bottom edge
        wf = .true.               !    wf is true
        i = i - 1                 !    move down & left
        j = j - 1
        if(j == 1)then            !  if moving to down left corner
          i = i + 1               !    move up
          dir = 8
        else
          dir = 4
        endif
      else                        !  if not moving to bottom edge
        i = i - 1                 !    just move down
        dir = 2
      endif

    CASE(2,3,7)                   !These forms exit element left
      if(j == 2)then              !  if moving into left edge
        wf = .true.               !    wf is true
        j = j - 1                 !    move left and up
        i = i + 1
        if(i == vi-1)then         !  if moving to top left corner
          j = j + 1               !    move right
          dir = 6
        else
          dir = 8
        endif
      else                        !  if not moving to left edge
        j = j - 1                 !    just move left
        dir = 4
      endif

    CASE(4,12,14)                 !These forms exit element right
      if(j == uj-2)then           !  if moving to right edge
        wf = .true.               !    wf is true
        j = j + 1                 !    move right and down
        i = i - 1
        if(i == 1)then            !  moving to down right corner
          j = j - 1               !    move left
          dir = 4
        else
          dir = 2
        endif
      else                        !  if not moving to right edge
        j = j + 1                 !    just move right
        dir = 6
      endif

    CASE(8,10,11)                 !These forms exit element top
      if(i == vi-2)then           !  if moving to top edge
        wf = .true.               !    wf is true
        i = i + 1                 !    move up and right
        j = j + 1
        if(j == uj-1)then         !  if moving to top right corner
          i = i - 1               !    move down
          dir = 2
        else
          dir = 6
        endif
      else                        !  if not moving to top edge
        i = i + 1                 !    just move up
        dir = 8
      endif

    CASE(16,17)                   !These forms exit right or left
                                  !16: u->r & d->l    17: d->r & u->l

                                  !if exit right
      if((form == 16 .AND. dir == 2).OR. (form == 17 .AND. dir == 8))then

        if(j == uj-2)then         !  if moving to right edge
          wf = .true.             !    wf is true
          j = j + 1               !    move right and down
          i = i - 1
          if(i == 1)then          !  moving to down right corner
            j = j - 1             !    move left
            dir = 4
          else
            dir = 2
          endif
        else                      !  if not moving to right edge
          j = j + 1               !    just move right
          dir = 6
        endif
      
      else                        !if exit left

        if(j == 2)then            !  if moving to left edge
          wf = .true.             !    wf is true
          j = j - 1               !    move left and up
          i = i + 1
          if(i == vi-1)then       !  if moving to top left corner
            j = j + 1             !    move right
            dir = 6
          else
            dir = 8
          endif
        else                      !  if not moving to left edge
          j = j - 1               !    just move left
          dir = 4
        endif
      
      endif

    CASE(18,19)                   !These forms exit up or down
                                  !18: r->u & l->d    19: l->u & r->d

                                  !if exit up
      if((form == 18 .AND. dir == 4).OR. (form == 19 .AND. dir == 6))then

        if(i == vi-2)then         !  if moving to top edge
          wf = .true.             !    wf is true
          i = i + 1               !    move up and right
          j = j + 1
          if(j == uj-1)then       !  if moving to top right corner
            i = i - 1             !    move down
            dir = 2
          else
            dir = 6
          endif
        else                      !  if not moving to top edge
          i = i + 1               !    just move up
          dir = 8
        endif

      else                        !if exit down

        if(i == 2)then            !  if moving to bottom edge
          wf = .true.             !    wf is true
          i = i - 1               !    move down and left
          j = j - 1
          if(j == 1)then          !  if moving to down left corner
            i = i + 1             !    move up
            dir = 8
          else
            dir = 4
          endif
        else                      !  if not moving to bottom edge
          i = i - 1               !    just move down
          dir = 2
        endif

      endif

    END SELECT
  endif

END SUBROUTINE getNext

  SUBROUTINE Get_coastdist(Ypos,Xpos,k_in,coastdist)

    USE PARAM_MOD, ONLY:us,Zgrid_depthinterp

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: k_in
    DOUBLE PRECISION, INTENT(IN) :: Ypos,Xpos
    DOUBLE PRECISION, INTENT(OUT) :: coastdist
    INTEGER i,j,start,bndpts,Mbnd,isle,islpts,klev
    DOUBLE PRECISION :: dXp1,dYp1,dX12,dY12,lenseg,param,Xseg,Yseg
    LOGICAL :: endMbnd,endIsle 
    
    !!!!_-----------------------------------------------------------------------

    i = 1                  !initialize i to 1
    klev=min(k_in,us_tridim)
    Mbnd = mid(i,klev)     !initialize isle to first boundary id
    start = 0              !initialize start to 0

    coastdist=9999999.
    do
      i = i + 1            !iterate through main boundary points
      endMbnd = .FALSE.
      if(i == tot_mbnd_pts(klev))then
        endMbnd = .TRUE.
      elseif(mid(i+1,klev) /= Mbnd) then
        endMbnd = .TRUE.
      endif
      If(endMbnd)then
        bndpts = i - start
        do j=2,bndpts
         dXp1 = Xpos - bx(start+j-1,klev)
         dYp1 = Ypos - by(start+j-1,klev)
         dX12 = bx(start+j,klev) - bx(start+j-1,klev)
         dY12 = by(start+j,klev) - by(start+j-1,klev)
         lenseg=( dX12 * dX12 + dY12 * dY12 )
         if(lenseg>0)then
          param =max(0.,min(1.,(dXp1*dX12+dYp1*dY12)/lenseg)) ! position (%) along segment of point closest to Xpos,Ypos
          Xseg=bx(start+j-1,klev)+ param*dX12                 ! point of boundary segment closest to (Xpos,Ypos)
          Yseg=by(start+j-1,klev)+ param*dY12                 ! point of boundary segment closest to (Xpos,Ypos)
          coastdist=min(coastdist,              &             ! closest distance between (Xpos,Ypos) and boundary segment
             sqrt((Xpos-Xseg)*(Xpos-Xseg)+(Ypos-Yseg)*(Ypos-Yseg))) 
         endif 
        enddo
        if(i == tot_mbnd_pts(klev) )exit
        start = i
        Mbnd = mid(i+1,klev)
      Endif
    enddo
    !!!!_-----------------------------------------------------------------------
    IF(num_ibnd(klev)>0)then
      i = 1                  !initialize i to 1
      isle = hid(i,klev)          !initialize isle to first island id
      start = 0              !initialize start to 0
      do
        i = i + 1            !iterate through island boundary points
        endIsle = .FALSE.
        if(i == tot_ibnd_pts(klev))then
          endIsle = .TRUE.
        elseif(hid(i+1,klev) /= isle) then
          endIsle = .TRUE.
        endif
        If(endIsle)then
          islpts = i - start
          do j=2,islpts
           dXp1 = Xpos - hx(start+j-1,klev)
           dYp1 = Ypos - hy(start+j-1,klev)
           dX12 = hx(start+j,klev) - hx(start+j-1,klev)
           dY12 = hy(start+j,klev) - hy(start+j-1,klev)
           lenseg=( dX12 * dX12 + dY12 * dY12 )
           if(lenseg>0)then
            param =max(0.,min(1.,(dXp1*dX12+dYp1*dY12)/lenseg)) ! position (%) along segment of point closest to Xpos,Ypos
            Xseg=hx(start+j-1,klev)+ param*dX12                 ! point of boundary segment closest to (Xpos,Ypos)
            Yseg=hy(start+j-1,klev)+ param*dY12                 ! point of boundary segment closest to (Xpos,Ypos)
            coastdist=min(coastdist,          &                 ! closest distance between (Xpos,Ypos) and boundary segment
               sqrt((Xpos-Xseg)*(Xpos-Xseg)+(Ypos-Yseg)*(Ypos-Yseg))) 
           endif 
          enddo 
          if(i == tot_ibnd_pts(klev) )exit
          start = i
          isle = hid(i+1,klev)
        Endif
      enddo
    ENDIF
    !!!!_-----------------------------------------------------------------------



  END SUBROUTINE Get_coastdist

  ! This subroutine uses the point-in-polygon approach to determine if a
  !   centroid is inside the model domain 
  ! It returns the value 1 in the variable inbounds if the centroid is in
  !   bounds and 0 if it is outside the boundaries
  SUBROUTINE mbounds(Ypos,Xpos,inbounds,mainboundnum,klev,nestdeg)

    USE PIP_MOD, ONLY: inpoly
    USE PARAM_MOD, ONLY:us,Zgrid_depthinterp
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: klev
    INTEGER, INTENT(OUT) :: inbounds,nestdeg
    DOUBLE PRECISION, INTENT(IN) :: Ypos,Xpos
    INTEGER, INTENT(OUT) :: mainboundnum


    LOGICAL :: endMbnd
    INTEGER :: i,j,start,bndpts
    INTEGER :: Mbnd

    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: blatlon
    mainboundnum=-1
    !read main boundaries into one variable for inpoly
    i = 1                  !initialize i to 1
    Mbnd = mid(i,klev)     !initialize isle to first boundary id
    start = 0              !initialize start to 0

    !write(*,*) 'looking for main boundary at level ',klev
    !write(*,*) 'num_mbnd(klev)=',num_mbnd(klev)
    !write(*,*) 'tot_mbnd_pts(klev)=',tot_mbnd_pts(klev)
    !if(klev<=us_tridim .and.tot_mbnd_pts(klev)>1)then
    do
      i = i + 1            !iterate through main boundary points
      !Determine if the end of a boundary has been reached by checking if it
      !  is either the very last point, or if not, if the next point is not
      !  part of the same island
      endMbnd = .FALSE.
      if(i == tot_mbnd_pts(klev))then
        endMbnd = .TRUE.
       !write(*,*) 'mid=',mid(i,klev)
      elseif(mid(i+1,klev) /= Mbnd) then
        endMbnd = .TRUE.
       !write(*,*) '    checking mid=',mid(i,klev)
      endif
      if(endMbnd)then
        !WRITE(*,'(5(a,i8))')'klev=',klev,' main bnd ',Mbnd,' from i=',start+1,' ->',start+i,' /',tot_mbnd_pts(klev)
        bndpts = i - start
        ALLOCATE(blatlon(2,bndpts))
        !write(*,*)
        do j=1,bndpts
          !write(*,'(2(a,i8,a,f20.12))')'bx(',start+j,')=',bx(start+j,klev),',  by(',start+j,')=',by(start+j,klev)
          blatlon(1,j) = bx(start+j,klev)
          blatlon(2,j) = by(start+j,klev)
        enddo
        if(inpoly(Xpos,Ypos,bndpts,blatlon)) then
          inbounds = 1
          !mainboundnum = dble(Mbnd)
          mainboundnum = max(mainboundnum,Mbnd)
          ! DEALLOCATE(blatlon)
          ! exit ! CL:OGS::: DO NOT EXIT AS WE WILL SEARCH THE HIGHER NESTED WATER
          ! BASIN DEGREE CONTAINING THE PARTICLE
        endif
        DEALLOCATE(blatlon)
        if(i == tot_mbnd_pts(klev) )exit
        start = i
        Mbnd = mid(i+1,klev)
      endif
    enddo
 
    if(mainboundnum>0)then
     nestdeg=polynestdeg(mainboundnum,klev)
    else
     nestdeg=-1
    endif

    !else
    !    write(*,*) 'ERROR mbounds klev=',klev,' mbound points=',tot_mbnd_pts(klev)
    !    stop
    !endif   
    !write(*,'(3(a,f9.1))') 'particle of pos',Xpos,' , ',Ypos,' found in main boundary number ',mainboundnum

  END SUBROUTINE mbounds



  ! This subroutine uses the point-in-polygon approach to determine if a
  ! centriod is inside one of the islands within the model domain 
  ! It returns the value 1 in the variable in_island if the centroid is in 
  !   island boundaries and 0 if it is outside the island boundaries
  SUBROUTINE ibounds(in_island,claty,clongx,island,klev,nestdeg)

    USE PIP_MOD, ONLY: inpoly
    USE PARAM_MOD, ONLY:us,Zgrid_depthinterp  
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: klev
    INTEGER, INTENT(OUT) :: in_island,nestdeg
    DOUBLE PRECISION, INTENT(IN) :: claty,clongx
    INTEGER, INTENT(OUT) :: island

    LOGICAL :: endIsle
    INTEGER :: i,j,start,islpts
    INTEGER :: isle
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: isbnds

    in_island = 0            !initialize in_island to 0
    island = 0             !initialize island to 0.0



    !write(*,*) 'looking for island boundary at level ',klev
!    write(*,*) 'tot_ibnd_pts(klev)=',tot_ibnd_pts(klev)
!    write(*,*) 'num_ibnd(klev)=',num_ibnd(klev)
 
    if(num_ibnd(klev)>0)then

      i = 1                  !initialize i to 1
      isle = hid(i,klev)          !initialize isle to first island id
      start = 0              !initialize start to 0
 !     write(*,*) 'hid(i+1,klev)=',hid(i+1,klev)

      do

        i = i + 1            !iterate through island boundary points

        !Determine if the end of an island has been reached by checking if it
        !  is either the very last point, or if not, if the next point is not
        !  part of the same island
        endIsle = .FALSE.
        if(i == tot_ibnd_pts(klev))then
          endIsle = .TRUE.
        ! write(*,*) 'hid=',hid(i,klev)
        elseif(hid(i+1,klev) /= isle) then
          endIsle = .TRUE.
        endif

        !if the end of an island has been reached:
        !  set islpts to the number of points around the island, allocate
        !  isbnds, then iterate through the island edge points storing
        !  the boundary locations into isbnds
        !  Next, call inpoly to see if the point is in the island
        !    If so, set island to the id of the island it is in, deallocate
        !      isbnds, and exit the subroutine
        !    If not, deallocate isbnds and repeat for all remaining islands
        if(endIsle)then
          islpts = i - start
          ALLOCATE(isbnds(2,islpts))
          do j=1,islpts
            isbnds(1,j) = hx(start+j,klev)
            isbnds(2,j) = hy(start+j,klev)
          enddo

          if(inpoly(clongx,claty,islpts,isbnds))then
            in_island = 1
            !island = isle 
            island = max(isle,island)
            !DEALLOCATE(isbnds)
            !write(*,*)'particle of pos (',clongx,claty,') at level ',klev,     &
             ! ' found to be in island number ',hid(i,klev)
            ! exit ! CL:OGS::: DO NOT EXIT AS WE WILL SEARCH THE HIGHER NESTED 
                   ! ISLAND  DEGREE CONTAINING THE PARTICLE
          endif
          DEALLOCATE(isbnds)
          if(i == tot_ibnd_pts(klev) )exit
          start = i
          isle = hid(i+1,klev)
        endif
      enddo

    endif

    if(island>0)then
     nestdeg=polynestdeg(island,klev)+num_mbnd(klev)-1000
    else
     nestdeg=-1
    endif
  END SUBROUTINE ibounds

  ! This subroutine calculates the intersection between the particle
  ! trajectory and the boundary line in a grid cell, and then calculates
  ! the reflection, returning the new particle location
  subroutine intersect_reflect(Xpos,Ypos,nXpos,nYpos,fintersectX,fintersectY,  &
    freflectX,freflectY,intersectf,skipbound,isWater,klev,bndintersect)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: klev
    INTEGER, INTENT(OUT) :: intersectf
    INTEGER, INTENT(INOUT) :: skipbound
    DOUBLE PRECISION, INTENT(IN) :: Xpos,Ypos,nXpos,nYpos
    DOUBLE PRECISION, INTENT(OUT) ::fintersectX,fintersectY,freflectX,freflectY
    LOGICAL, OPTIONAL, INTENT(OUT) :: isWater
    INTEGER :: i,intersect,skipboundi
    DOUBLE PRECISION :: crossk,dPBC,mBCperp,rx1,rx2,ry1,ry2,Bp,distBC,dist1,   &
      dist2,intersctx,interscty,rPxyzX,rPxyzY,Mbc,Bbc,Mp,bcx1,bcx2,bcy1,bcy2,  &
      bBCperp,xhigh,xlow,yhigh,ylow,d_Pinter,dtest,bxhigh,bxlow,byhigh,bylow
    DOUBLE PRECISION :: bndintersect(12),mindist     !--- CL:OGS
    INTEGER:: closerbound                                 !--- CL:OGS
    DOUBLE PRECISION :: bcx1_array(nbounds(klev)),bcx2_array(nbounds(klev)),bcy1_array(nbounds(klev)),bcy2_array(nbounds(klev)),dist(nbounds(klev)) 
     
    distBC=0.0
    Mbc = 0.0
    Bbc = 0.0
    Mp = 0.0
    Bp = 0.0
    intersect=0
    intersectf=0
    skipboundi = skipbound
    fintersectX = -999999.
    fintersectY = -999999.
    freflectX = -999999.
    freflectY = -999999.
    closerbound=-1
    dtest = 999999.
    mindist = 999999.
    isWater = .FALSE.

    if (Xpos.GE.nXpos) then
      xhigh = Xpos
      xlow = nXpos
    else
      xhigh = nXpos
      xlow = Xpos
    endif  

    if (Ypos.GE.nYpos) then
      yhigh = Ypos
      ylow = nYpos
    else
      yhigh = nYpos
      ylow = Ypos
    endif

    bcx1_array=bnd_x1(:nbounds(klev),klev)
    bcx2_array=bnd_x2(:nbounds(klev),klev)
    bcy1_array=bnd_y1(:nbounds(klev),klev)
    bcy2_array=bnd_y2(:nbounds(klev),klev)
    dist=    sqrt((bcx1_array- Xpos)**2+(bcy1_array- Ypos)**2) 
    dist=min(sqrt((bcx2_array- Xpos)**2+(bcy2_array- Ypos)**2),dist)
    dist=min(sqrt((bcx2_array-nXpos)**2+(bcy2_array-nYpos)**2),dist)

    do i=1,nbounds(klev)

      if (i /= skipbound)then
         continue
      else
         cycle
      endif

        intersect = 0
        bcx1=bcx1_array(i)
        bcx2=bcx2_array(i)
        bcy1=bcy1_array(i)
        bcy2=bcy2_array(i)
        if( (mindist>dist(i)))then
            closerbound=i
            mindist=dist(i)
       endif

        !If the boundary segment end points are both east, west, north, or 
        !  south of the particle's previous or new location, cycle to next 
        !  boundary
        if( ((bcx1 > xhigh) .AND. (bcx2 > xhigh)) .OR. &
            ((bcx1 < xlow ) .AND. (bcx2 < xlow )) .OR. &
            ((bcy1 > yhigh) .AND. (bcy2 > yhigh)) .OR. &
            ((bcy1 < ylow ) .AND. (bcy2 < ylow ))      ) cycle
        
        if (bcx1.GE.bcx2) then
          bxhigh = bcx1
          bxlow = bcx2
        else
          bxhigh = bcx2
          bxlow = bcx1
        endif  

        if (bcy1.GE.bcy2) then
          byhigh = bcy1
          bylow = bcy2
        else
          byhigh = bcy2
          bylow = bcy1
        endif

        !First determine if an undefined denominator is possible
        if (bcx1.EQ.bcx2 .OR. nXpos.EQ.Xpos ) then
          !test if they both vertical, if so cycle because they cannot intersect
          if (bcx1.EQ.bcx2 .AND. nXpos.EQ.Xpos ) cycle
          !test if perpendicular and parrallel to coordinate axes
          if (bcx1.EQ.bcx2 .AND. nYpos.EQ.Ypos ) then
            !undefined denominator, perp. & || to axes
            intersctx = bcx1
            interscty = nYpos
            if (intersctx.LE.xhigh  .AND. intersctx.GE.xlow .AND.              &
              interscty.LE.yhigh  .AND. interscty.GE.ylow .AND.                &
              intersctx.LE.bxhigh .AND. intersctx.GE.bxlow .AND.               &
              interscty.LE.byhigh .AND. interscty.GE.bylow  ) then
              dPBC=sqrt((intersctx-nXpos)**2+(interscty-nYpos)**2)
              rx1=nXpos+(DBLE(2.0)*dPBC)
              ry1=nYpos
              rx2=nXpos-(DBLE(2.0)*dPBC)
              ry2=nYpos
              dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
              dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 )
              if(dist1.LT.dist2) then
                rPxyzX= rx1
                rPxyzY= ry1
              else !elseif(dist1.GE.dist2) then ! BEUG corrected:GT->GE (CL)
                rPxyzX= rx2
                rPxyzY= ry2
              endif
              intersect=1
            endif
          elseif (nXpos.EQ.Xpos .AND. bcy1.EQ.bcy2 ) then
            !undefined denominator, perp. & || to axes
            intersctx = nXpos
            interscty = bcy1
            if (intersctx.LE.xhigh .AND.  intersctx.GE.xlow .AND.              &
              interscty.LE.yhigh  .AND. interscty.GE.ylow .AND.                &
              intersctx.LE.bxhigh .AND. intersctx.GE.bxlow .AND.               &
              interscty.LE.byhigh .AND. interscty.GE.bylow  ) then
              dPBC=sqrt((intersctx-nXpos)**2+(interscty-nYpos)**2)
              rx1=nXpos
              ry1=nYpos+(DBLE(2.0)*dPBC)
              rx2=nXpos
              ry2=nYpos-(DBLE(2.0)*dPBC)
              dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
              dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 )
              if(dist1.LT.dist2) then
                rPxyzX= rx1
                rPxyzY= ry1
              else !elseif(dist1.GE.dist2) then ! BEUG corrected:GT->GE (CL)
                rPxyzX= rx2
                rPxyzY= ry2
              endif
              intersect=1
            endif
          elseif (bcx1.EQ.bcx2 .AND. nYpos.NE.Ypos ) then
            !undefined denominator, not perpendicular
            Mp = (nYpos-Ypos)/(nXpos-Xpos)
            Bp = Ypos - Mp*Xpos
            intersctx = bcx1
            interscty = Mp*intersctx + Bp
            if (intersctx.LE.xhigh .AND.  intersctx.GE.xlow .AND.              &
                interscty.LE.yhigh  .AND. interscty.GE.ylow .AND.              &
                intersctx.LE.bxhigh .AND. intersctx.GE.bxlow .AND.             &
                interscty.LE.byhigh .AND. interscty.GE.bylow  ) then
              dPBC = nXpos-intersctx
              rx1=nXpos+(DBLE(2.0)*dPBC)
              ry1=nYpos
              rx2=nXpos-(DBLE(2.0)*dPBC)
              ry2=nYpos
              dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
              dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 )
              if(dist1.LT.dist2) then
                rPxyzX= rx1
                rPxyzY= ry1
              else !elseif(dist1.GE.dist2) then ! BEUG corrected:GT->GE (CL)
                rPxyzX= rx2
                rPxyzY= ry2
              endif
              intersect=1
            endif
          elseif (nXpos.EQ.Xpos .AND. bcy1.NE.bcy2  ) then
            !undefined denominator, not perpendicular
            Mbc = (bcy2-bcy1)/(bcx2-bcx1)
            Bbc = bcy2 - Mbc*bcx2
            intersctx = nXpos
            interscty = Mbc*intersctx + Bbc
            if (intersctx.LE.xhigh .AND.  intersctx.GE.xlow .AND.              &
                interscty.LE.yhigh  .AND. interscty.GE.ylow .AND.              &
                intersctx.LE.bxhigh .AND. intersctx.GE.bxlow .AND.             &
                interscty.LE.byhigh .AND. interscty.GE.bylow ) then
              !Now use cross product to determine the distance of the particle 
              !  from the boundary
              distBC = sqrt((bcx1-bcx2)**2+(bcy1-bcy2)**2)
              crossk= ((nXpos-bcx1)*(bcy2-bcy1)) - ((bcx2-bcx1)*(nYpos-bcy1))
              dPBC = sqrt(crossk**2)/distBC
              !find line perpendicular to boundary
              mBCperp = DBLE(-1.0)/Mbc
              bBCperp = nYpos - mBCperp*nXpos
              !find two potential reflection points
              rx1 = sqrt( ((DBLE(2.0)*dPBC)**2)/(DBLE(1.0)+mBCperp**2) ) +nXpos
              ry1 = mBCperp*rx1 + bBCperp
              rx2 = sqrt( ((DBLE(2.0)*dPBC)**2)/(DBLE(1.0)+mBCperp**2) )       &
                  * DBLE(-1.0) + nXpos
              ry2 = mBCperp*rx2 + bBCperp
              !point closest to intersection of boundary and particle trajectory
              !  is the right one
              dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
              dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 )
              if(dist1.LT.dist2) then
                rPxyzX= rx1
                rPxyzY= ry1
              else !elseif(dist1.GE.dist2) then ! BEUG corrected:GT->GE (CL)
                rPxyzX= rx2
                rPxyzY= ry2
              endif
              intersect=1
            endif
          endif
        else

          if(intersect == 0)then

            Mbc = (bcy2-bcy1)/(bcx2-bcx1)
            Bbc = bcy2 - Mbc*bcx2
            Mp = (nYpos-Ypos)/(nXpos-Xpos)
            Bp = Ypos - Mp*Xpos
            intersctx = (Bbc - Bp)/(Mp - Mbc)
            interscty = Mp*intersctx + Bp

            !when bc parallel with x-axis, byhigh=bylow=intersecty
            if (Mbc.EQ.0.0) interscty = byhigh
        
            if (intersctx.LE.xhigh .AND.  intersctx.GE.xlow .AND.              &
                interscty.LE.yhigh  .AND. interscty.GE.ylow .AND.              &
                intersctx.LE.bxhigh .AND. intersctx.GE.bxlow .AND.             &
                interscty.LE.byhigh .AND. interscty.GE.bylow  ) then

              if (Mbc.EQ.0.0) then  !inverse slope denominator not OK
                dPBC = nYpos-bcy1
                rx1=nXpos
                ry1=nYpos+(DBLE(2.0)*dPBC)
                rx2=nXpos
                ry2=nYpos-(DBLE(2.0)*dPBC)
                dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
                dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 ) 
                if(dist1.LT.dist2) then
                  rPxyzX= rx1
                  rPxyzY= ry1
                else !elseif(dist1.GE.dist2) then ! BEUG corrected:GT->GE (CL)
                  rPxyzX= rx2
                  rPxyzY= ry2
                endif
                  intersect=1
                endif  

              if(intersect == 0)then

                !Now use cross product to determine the distance of the
                !  particle from the boundary
                distBC = sqrt((bcx1-bcx2)**2+(bcy1-bcy2)**2)
                crossk= ((nXpos-bcx1)*(bcy2-bcy1)) - ((bcx2-bcx1)*(nYpos-bcy1))
                dPBC = sqrt(crossk**2)/distBC
                !find line perpendicular to boundary
                mBCperp = DBLE(-1.0)/Mbc
                bBCperp = nYpos - mBCperp*nXpos
                !find two potential reflection points
                rx1 = sqrt(((DBLE(2.0)*dPBC)**2)/(DBLE(1.0)+mBCperp**2)) +nXpos
                ry1 = mBCperp*rx1 + bBCperp
                rx2 = sqrt(((DBLE(2.0)*dPBC)**2)/(DBLE(1.0)+mBCperp**2))       &
                    * DBLE(-1.0) + nXpos
                ry2 = mBCperp*rx2 + bBCperp
                !point closest to intersection of boundary and particle 
                !  trajectory is the right one
                dist1 = sqrt( (intersctx-rx1)**2 + (interscty-ry1)**2 )
                dist2 = sqrt( (intersctx-rx2)**2 + (interscty-ry2)**2 )
                if(dist1.LT.dist2) then
                  rPxyzX= rx1
                  rPxyzY= ry1
                else !elseif(dist1.GE.dist2) then ! BEUG corrected:GT->GE (CL)
                  rPxyzX= rx2
                  rPxyzY= ry2
                endif
                intersect=1
                
              endif
            endif  
          endif
        endif

        
        d_Pinter = sqrt( (Xpos-intersctx)**2 + (Ypos-interscty)**2 )
        if( (intersect .EQ. 1) .AND. (d_Pinter .LT. dtest) ) then
          fintersectX = intersctx
          fintersectY = interscty
          freflectX = rPxyzX
          freflectY = rPxyzY
          intersectf = 1
          dtest = d_Pinter
          skipboundi = i
          isWater = .NOT. land(i,klev)
          bndintersect(1)=bcx1                      !--- CL:OGS
          bndintersect(3)=bcx2                      !--- CL:OGS
          bndintersect(2)=bcy1                      !--- CL:OGS
          bndintersect(4)=bcy2                      !--- CL:OGS

          bndintersect(5)=bnd_x1(max(1,i-1),klev)           !--- CL:OGS
          bndintersect(7)=bnd_y1(max(1,i-1),klev)           !--- CL:OGS
          bndintersect(6)=bnd_x2(max(1,i-1),klev)           !--- CL:OGS
          bndintersect(8)=bnd_y2(max(1,i-1),klev)           !--- CL:OGS

          bndintersect( 9)=bnd_x1(min(nbounds(klev),i+1),klev)           !--- CL:OGS
          bndintersect(11)=bnd_y1(min(nbounds(klev),i+1),klev)           !--- CL:OGS
          bndintersect(10)=bnd_x2(min(nbounds(klev),i+1),klev)           !--- CL:OGS
          bndintersect(12)=bnd_y2(min(nbounds(klev),i+1),klev)           !--- CL:OGS

        endif
       
    enddo
    if( (intersect .EQ. 0))then
          i=min(max(1,closerbound),nbounds(klev))
          bndintersect(1)=bnd_x1(i,klev)                                !--- CL:OGS
          bndintersect(3)=bnd_y1(i,klev)                                !--- CL:OGS
          bndintersect(2)=bnd_x2(i,klev)                                !--- CL:OGS
          bndintersect(4)=bnd_y2(i,klev)                                !--- CL:OGS

          bndintersect(5)=bnd_x1(max(1,i-1),klev)                       !--- CL:OGS
          bndintersect(7)=bnd_y1(max(1,i-1),klev)                       !--- CL:OGS
          bndintersect(6)=bnd_x2(max(1,i-1),klev)                       !--- CL:OGS
          bndintersect(8)=bnd_y2(max(1,i-1),klev)                       !--- CL:OGS

          bndintersect( 9)=bnd_x1(min(nbounds(klev),i+1),klev)           !--- CL:OGS
          bndintersect(11)=bnd_y1(min(nbounds(klev),i+1),klev)           !--- CL:OGS
          bndintersect(10)=bnd_x2(min(nbounds(klev),i+1),klev)           !--- CL:OGS
          bndintersect(12)=bnd_y2(min(nbounds(klev),i+1),klev)           !--- CL:OGS


    endif
    skipbound = skipboundi
  END SUBROUTINE intersect_reflect


SUBROUTINE output_xyBounds()
!This subroutine creates a blanking file for surfer of the model 
!  boundaries in meters
  USE PARAM_MOD, ONLY: ui,uj,vi,vj,us
  USE HYDRO_MOD, ONLY: getUVxy
  IMPLICIT NONE

  INTEGER :: i,j,klev,m,numpoly,bnd,STATUS,pstart,pend
  INTEGER, ALLOCATABLE, DIMENSION(:) :: polysizes
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: x_u,y_u,x_v,y_v

  ALLOCATE(x_u(ui,uj))
  ALLOCATE(y_u(ui,uj))
  ALLOCATE(x_v(vi,vj))
  ALLOCATE(y_v(vi,vj))

  write(*,*) 'output metric blanking file'

  !get the x/y locations of the U and V grid nodes
  CALL getUVxy(x_u,y_u,x_v,y_v)

  open(1,FILE='xybounds.bln',STATUS='REPLACE')
  DO klev=1,us_tridim
    if(bnum(klev)>0)then
      numpoly = bnds(bnum(klev),klev)%poly       !numpoly = # of polygons
     
      !Allocate polysizes to the number of polygons
      allocate(polysizes(numpoly),STAT=STATUS)
      IF(STATUS /= 0) write(*,*) 'Problem allocate polysizes'
     
      polysizes = 0       !initialize all polygon sizes to 0
     
      !iterate through boundaries incrementing polysize depending on
      !   polygon each boundary is part of
      do bnd=1,bnum(klev)
        polysizes(bnds(bnd,klev)%poly) = polysizes(bnds(bnd,klev)%poly) + 1
      enddo
     
      pstart = 1
     
      do m=1,numpoly                          !iterate through the polygons
        write(1,*) polysizes(m)+1,',',1       !write polysize,1 for .bln 1st row
        pend = pstart + polysizes(m) - 1      !calculate last bnds location
        do bnd=pstart,pend                      !iterate through polygon bnds
          i = bnds(bnd,klev)%ii                      !i position
          j = bnds(bnd,klev)%jj                      !j position
          if( bnds(bnd,klev)%onU ) then                !if on U Grid
            write(1,*) x_u(i,j),',',y_u(i,j),',',klev  !  write x & y of U node
          else                                      !else
            write(1,*) x_v(i,j),',',y_v(i,j),',',klev  !  write x & y of V node
          endif
        enddo
        bnd = pstart                            !go back to polygon's first bnd
        i = bnds(bnd,klev)%ii                        !  to close the polygon
        j = bnds(bnd,klev)%jj
        if( bnds(bnd,klev)%onU ) then
          write(1,*) x_u(i,j),',',y_u(i,j),',',klev
        else
          write(1,*) x_v(i,j),',',y_v(i,j),',',klev
        endif
     
        pstart = pstart + polysizes(m)        !calculate start of next polygon
     
      enddo
     
      deallocate(polysizes)
    else 
      write(*,*)'warning : output_xybounds: no bounds found for level',klev
    endif

  enddo ! klev=1,us_tridim

  close(1)

  DEALLOCATE(x_u,y_u,x_v,y_v)

END SUBROUTINE output_xyBounds


SUBROUTINE output_llBounds()
!This subroutine creates a blanking file for Surfer/Scripter of the model 
!  boundaries in longitude and latitude
  USE PARAM_MOD,   ONLY: ui,uj,vi,vj,us
  USE HYDRO_MOD,   ONLY: getUVxy
  USE CONVERT_MOD, ONLY: x2lon,y2lat
  IMPLICIT NONE

  INTEGER :: i,j,klev,bnd,m,numpoly,STATUS,pstart,pend
  INTEGER, ALLOCATABLE, DIMENSION(:) :: polysizes
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: x_u,y_u,x_v,y_v

  ALLOCATE(x_u(ui,uj))
  ALLOCATE(y_u(ui,uj))
  ALLOCATE(x_v(vi,vj))
  ALLOCATE(y_v(vi,vj))

  write(*,*) 'output lat/long blanking file'

  !get the x/y locations of the U and V grid nodes
  CALL getUVxy(x_u,y_u,x_v,y_v)

  open(1,FILE='llbounds.bln',STATUS='REPLACE')
  
  DO klev=1,us_tridim
    if(bnum(klev)>0)then
      numpoly = bnds(bnum(klev),klev)%poly       !numpoly = # of polygons
    
      !Allocate polysizes to the number of polygons
      allocate(polysizes(numpoly),STAT=STATUS)
      IF(STATUS /= 0) write(*,*) 'Problem allocate polysizes'
    
      polysizes = 0       !initialize all polygon sizes to 0
    
      !iterate through boundaries incrementing polysize depending on
      !   polygon each boundary is part of
      do bnd=1,bnum(klev)
        polysizes(bnds(bnd,klev)%poly) = polysizes(bnds(bnd,klev)%poly) + 1
      enddo
    
      pstart = 1
    
      do m=1,numpoly                      !iterate through the polygons
        if(m<=num_mbnd(klev))then
          write(1,*) polysizes(m)+1,',',1,',',polynestdeg(m,klev)   !write polysize,1,1 for first row of land bnd
        else                                  
          write(1,*) polysizes(m)+1,',',0,',',polynestdeg(m,klev)   !write polysize,1,1 for first row of island
        endif
        pend = pstart + polysizes(m) - 1  !calculate polygon last bnds location
        do bnd=pstart,pend                  !iterate through polygon bnds
          i = bnds(bnd,klev)%ii                  !i position
          j = bnds(bnd,klev)%jj                  !j position
          if( bnds(bnd,klev)%onU ) then          !if on U Grid write lon & lat of U node
            write(1,*) x2lon(x_u(i,j),y_u(i,j)),',',y2lat(y_u(i,j)),',',klev
          else                            !else write lon & lat of V node
            write(1,*) x2lon(x_v(i,j),y_v(i,j)),',',y2lat(y_v(i,j)),',',klev
          endif
        enddo
        bnd = pstart                        !go back to polygon's first bnd
        i = bnds(bnd,klev)%ii                    !  to close the polygon
        j = bnds(bnd,klev)%jj
        if( bnds(bnd,klev)%onU ) then
          write(1,*) x2lon(x_u(i,j),y_u(i,j)),',',y2lat(y_u(i,j)),',',klev
        else
          write(1,*) x2lon(x_v(i,j),y_v(i,j)),',',y2lat(y_v(i,j)),',',klev
        endif
    
        pstart = pstart + polysizes(m)    !calculate start of next polygon
    
      enddo
      deallocate(polysizes)
    else 
      write(*,*)'warning : output_llbounds: no bounds found for level',klev
    endif
  enddo ! klev=1,us_tridim
  close(1)

  DEALLOCATE(x_u,y_u,x_v,y_v)

END SUBROUTINE output_llBounds

!   !*********************************************************
!   !*        write nestingdegree in netcdf file             *
!   !*********************************************************

SUBROUTINE  createNestingDegreeNetCDF(nestingdegree)
    USE PARAM_MOD, ONLY: NCOutFile,outpath,outpathGiven,   &
        ExeDir,OutDir, ui,uj,vi,vj
    USE netcdf
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
    INTEGER, INTENT(IN) :: nestingdegree(:,:,:)
    INTEGER :: STATUS,NCID,NestDegID,       &
              niRID,njRID,nkRID,count,i,j
    CHARACTER(LEN=300) :: ncFile

    IF(outpathGiven)THEN
        ncFile = TRIM(outpath) // '/' // TRIM(NCOutFile) // '_NestingDegree.nc'
    ELSE
        ncFile = TRIM(NCOutFile) // '_NestingDegree.nc'
    ENDIF

    write(*,*)'Creating NetCDF Output File: ',TRIM(ncFile)

    STATUS = NF90_CREATE(TRIM(ncFile), NF90_CLOBBER, NCID)
    IF(STATUS /= NF90_NOERR) THEN
      WRITE(*,*) 'Problem creating NetCDF output file'
      WRITE(*,*) NF_STRERROR(STATUS)
      STOP
    ENDIF

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !NF90_DEF_DIM
        STATUS = NF90_DEF_DIM(NCID,'niR',vi,niRID)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: niR dim'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
        STATUS = NF90_DEF_DIM(NCID,'njR',uj,njRID)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: njR dim'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)
        STATUS = NF90_DEF_DIM(NCID,'nkR',us_tridim,nkRID)
        IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: nkR dim'
        IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !NF90_DEF_VAR

    STATUS = NF90_DEF_VAR(NCID,'NestDeg',NF_DOUBLE,(/niRID,njRID,nkRID/),NestDegID)
    IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: NestDeg var'
    IF(STATUS /= NF90_NOERR) WRITE(*,*) NF_STRERROR(STATUS)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    !NF_ENDDEF

      STATUS = NF90_ENDDEF(NCID)
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: EndDef'
      IF(status /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    !NF_PUT_VAR

      !NestDeg   
      STATUS = NF90_INQ_VARID(NCID, "NestDeg",          NestDegID)
      STATUS = NF90_PUT_VAR(NCID,    NestDegID,         nestingdegree)
      IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem put NestDeg '
      IF(STATUS /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    !NF_CLOSE

    STATUS = NF_CLOSE(NCID)
    IF(STATUS /= NF90_NOERR) WRITE(*,*) 'Problem createNetCDF: Close'
    IF(status /= NF90_NOERR) WRITE(*,*) NF90_STRERROR(STATUS)



END SUBROUTINE createNestingDegreeNetCDF


!   !*********************************************************
!   !*                     geteleform                        *
!   !*********************************************************
!SUBROUTINE geteleform(ni,nj,nk,eleform)      ! CL-OGS: created to get Depth of
!                                        ! particles in hydro module
!  USE PARAM_MOD, ONLY: ui,uj,vi,vj,us
!  IMPLICIT NONE
!  INTEGER,INTENT(IN):: ni,nj,nk
!  INTEGER, INTENT(INOUT):: eleform(:,:,:)
!     if((ni.ne.vi-1).or.(nj.ne.uj-1).or.(nk.ne.us))then
!        write(*,*)'ERROR DIM eleform requested to hydro mod'
!        stop
!     endif
!     eleform(:,:,:)=ele(:,:,:)%form
!END SUBROUTINE geteleform
!SUBROUTINE geteleform(i,j,eleform)      ! CL-OGS: created to get Depth of
!                                        ! particles in hydro module
!  USE PARAM_MOD, ONLY:us 
!  IMPLICIT NONE
!  INTEGER,INTENT(IN)::    i,j
!  INTEGER, INTENT(INOUT):: eleform(:)
!  INTEGER :: k
!  do k=1,us
!    eleform(k)=ele(i,j,k)%form
!END SUBROUTINE geteleform

SUBROUTINE finBoundary()
  IMPLICIT NONE

  DEALLOCATE(bnum)
  DEALLOCATE(hid,hx,hy)
  DEALLOCATE(mid,bx,by)
  DEALLOCATE(bnd_x1,bnd_x2,bnd_y1,bnd_y2)
  DEALLOCATE(land)
  DEALLOCATE(nbounds)
  DEALLOCATE(tot_mbnd_pts)
  DEALLOCATE(tot_ibnd_pts)
  DEALLOCATE(num_ibnd)
  DEALLOCATE(num_mbnd)
  DEALLOCATE(BND_SET)
  DEALLOCATE (ele)
  DEALLOCATE(polynestdeg)

END SUBROUTINE finBoundary


END MODULE
