!
!(c) Matthew Kennel, Institute for Nonlinear Science (2004)
!
! Licensed under the Academic Free License version 1.1 found in file LICENSE
! with additional provisions found in that same file.
!
module tools_kdtree2
  use tools_func, only: sphere_dist
  use tools_kdtree2_pq
  use tools_kinds, only: kind_real
  use tools_stripack, only: scoord
  ! K-D tree routines in Fortran 90 by Matt Kennel.
  ! Original program was written in Sather by Steve Omohundro and
  ! Matt Kennel.  Only the Euclidean metric is supported.
  !
  !
  ! This module is identical to 'kd_tree', except that the order
  ! of subscripts is reversed in the data file.
  ! In otherwords for an embedding of N D-dimensional vectors, the
  ! data file is here, in natural Fortran order  data(1:D, 1:N)
  ! because Fortran lays out columns first,
  !
  ! whereas conventionally (C-style) it is data(1:N,1:D)
  ! as in the original kd_tree module.
  !
  !-------------DATA TYPE, CREATION, DELETION---------------------
  public :: kdtree2, kdtree2_result, tree_node, kdtree2_create, kdtree2_destroy
  !---------------------------------------------------------------
  !-------------------SEARCH ROUTINES-----------------------------
  public :: kdtree2_n_nearest
  ! Return fixed number of nearest neighbors around arbitrary vector
  !----------------------------------------------------------------

  integer, parameter :: bucket_size = 12
  ! The maximum number of points to keep in a terminal node.

  type interval
      real(kind_real) :: lower,upper
  end type interval

  type :: tree_node
      ! an internal tree node
      private
      integer :: cut_dim
      ! the dimension to cut
      real(kind_real) :: cut_val
      ! where to cut the dimension
      real(kind_real) :: cut_val_left, cut_val_right
      ! improved cutoffs knowing the spread in child boxes.
      integer :: l, u
      type (tree_node), pointer :: left, right
      type(interval), pointer :: box(:) => null()
      ! child pointers
      ! Points included in this node are indexes[k] with k \in [l,u]
  end type tree_node

  type :: kdtree2
      ! Global information about the tree, one per tree
      integer :: n=0
      ! total # of points
      real(kind_real), pointer :: the_data(:,:) => null()
      ! pointer to the actual data array
      !
      !  IMPORTANT NOTE:  IT IS DIMENSIONED   the_data(1:d,1:N)
      !  which may be opposite of what may be conventional.
      !  This is, because in Fortran, the memory layout is such that
      !  the first dimension is in sequential order.  Hence, with
      !  (1:d,1:N), all components of the vector will be in consecutive
      !  memory locations.  The search time is dominated by the
      !  evaluation of distances in the terminal nodes.  Putting all
      !  vector components in consecutive memory location improves
      !  memory cache locality, and hence search speed, and may enable
      !  vectorization on some processors and compilers.

      integer, pointer :: ind(:) => null()
      ! permuted index into the data, so that indexes[l..u] of some
      ! bucket represent the indexes of the actual points in that
      ! bucket.
      logical       :: sort = .false.
      ! do we always sort output results?
      logical       :: rearrange = .false.
      real(kind_real), pointer :: rearranged_data(:,:) => null()
      ! if (rearrange .eqv. .true.) then rearranged_data has been
      ! created so that rearranged_data(:,i) = the_data(:,ind(i)),
      ! permitting search to use more cache-friendly rearranged_data, at
      ! some initial computation and storage cost.
      type (tree_node), pointer :: root => null()
      ! root pointer of the tree
  end type kdtree2


  type :: tree_search_record
      !
      ! One of these is created for each search.
      !
      private
      !
      ! Many fields are copied from the tree structure, in order to
      ! speed up the search.
      !
      integer           :: nn, nfound
      real(kind_real)      :: ballsize
      integer           :: centeridx=999, correltime=9999
      ! exclude points within 'correltime' of 'centeridx', iff centeridx >= 0
      integer           :: nalloc  ! how much allocated for results(:)?
      logical           :: rearrange  ! are the data rearranged or original?
      ! did the # of points found overflow the storage provided?
      logical           :: overflow
      real(kind_real), pointer :: qv(:)  ! query vector
      type(kdtree2_result), pointer :: results(:) ! results
      type(pq) :: pq
      real(kind_real), pointer :: data(:,:)  ! temp pointer to data
      integer, pointer      :: ind(:)     ! temp pointer to indexes
  end type tree_search_record

  real(kind_real),parameter :: rth = 1.0e-8 !< Reproducibility threshold

  private
  ! everything else is private.

  type(tree_search_record), save, target :: sr   ! A GLOBAL VARIABLE for search

contains

  function kdtree2_create(input_data,sort,rearrange) result (mr)
    !
    ! create the actual tree structure, given an input array of data.
    !
    ! Note, input data is input_data(1:d,1:N), NOT the other way around.
    ! THIS IS THE REVERSE OF THE PREVIOUS VERSION OF THIS MODULE.
    ! The reason for it is cache friendliness, improving performance.
    !
    ! Optional arguments:  if sort .eqv. .true. then output results
    !                      will be sorted by increasing distance.
    !                      default=.false., as it is faster to not sort.
    !
    !                      if rearrange .eqv. .true. then an internal
    !                      copy of the data, rearranged by terminal node,
    !                      will be made for cache friendliness.
    !                      default=.true., as it speeds searches, but
    !                      building takes longer, and extra memory is used.
    !
    ! .. Function Return Cut_value ..
    type (kdtree2), pointer :: mr
    logical, intent(in), optional      :: sort
    logical, intent(in), optional      :: rearrange
    ! ..
    ! .. Array Arguments ..
    real(kind_real), target :: input_data(:,:)
    !
    integer :: i
    ! ..
    allocate (mr)
    mr%the_data => input_data
    ! pointer assignment

    mr%n = size(input_data,2)

    call build_tree(mr)

    if (present(sort)) then
       mr%sort = sort
    else
       mr%sort = .false.
    endif

    if (present(rearrange)) then
       mr%rearrange = rearrange
    else
       mr%rearrange = .true.
    endif

    if (mr%rearrange) then
       allocate(mr%rearranged_data(3,mr%n))
       do i=1,mr%n
          mr%rearranged_data(:,i) = mr%the_data(:, &
           mr%ind(i))
       enddo
    else
       nullify(mr%rearranged_data)
    endif

  end function kdtree2_create

    subroutine build_tree(tp)
      type (kdtree2), pointer :: tp
      ! ..
      integer :: j
      type(tree_node), pointer :: dummy => null()
      ! ..
      allocate (tp%ind(tp%n))
      forall (j=1:tp%n)
         tp%ind(j) = j
      end forall
      tp%root => build_tree_for_range(tp,1,tp%n, dummy)
    end subroutine build_tree

    recursive function build_tree_for_range(tp,l,u,parent) result (res)
      ! .. Function Return Cut_value ..
      type (tree_node), pointer :: res
      ! ..
      ! .. Structure Arguments ..
      type (kdtree2), pointer :: tp
      type (tree_node),pointer           :: parent
      ! ..
      ! .. Scalar Arguments ..
      integer, intent (In) :: l, u
      ! ..
      ! .. Local Scalars ..
      integer :: i, c, m
      logical :: recompute
      real(kind_real)    :: average

!!$      If (.False.) Then
!!$         If ((l .Lt. 1) .Or. (l .Gt. tp%n)) Then
!!$            Stop 'illegal L value in build_tree_for_range'
!!$         End If
!!$         If ((u .Lt. 1) .Or. (u .Gt. tp%n)) Then
!!$            Stop 'illegal u value in build_tree_for_range'
!!$         End If
!!$         If (u .Lt. l) Then
!!$            Stop 'U is less than L, thats illegal.'
!!$         End If
!!$      Endif
!!$
      ! first compute min and max
      allocate (res)
      allocate(res%box(3))

      ! First, compute an APPROXIMATE bounding box of all points associated with this node.
      if ( u < l ) then
         ! no points in this box
         nullify(res)
         return
      end if

      if ((u-l)<=bucket_size) then
         !
         ! always compute true bounding box for terminal nodes.
         !
         do i=1,3
            call spread_in_coordinate(tp,i,l,u,res%box(i))
         end do
         res%cut_dim = 0
         res%cut_val = 0.0
         res%l = l
         res%u = u
         res%left =>null()
         res%right => null()
      else
         !
         ! modify approximate bounding box.  This will be an
         ! overestimate of the true bounding box, as we are only recomputing
         ! the bounding box for the dimension that the parent split on.
         !
         ! Going to a true bounding box computation would significantly
         ! increase the time necessary to build the tree, and usually
         ! has only a very small difference.  This box is not used
         ! for searching but only for deciding which coordinate to split on.
         !
         do i=1,3
            recompute=.true.
            if (associated(parent)) then
               if (i .ne. parent%cut_dim) then
                  recompute=.false.
               end if
            endif
            if (recompute) then
               call spread_in_coordinate(tp,i,l,u,res%box(i))
            else
               res%box(i) = parent%box(i)
            endif
         end do


         c = maxloc(res%box%upper-res%box%lower,1)
         !
         ! c is the identity of which coordinate has the greatest spread.
         !

         if (.false.) then
            ! select exact median to have fully balanced tree.
            m = (l+u)/2
            call select_on_coordinate(tp%the_data,tp%ind,c,m,l,u)
         else
            !
            ! select point halfway between min and max, as per A. Moore,
            ! who says this helps in some degenerate cases, or
            ! actual arithmetic average.
            !
            if (.true.) then
               ! actually compute average
               average = sum(tp%the_data(c,tp%ind(l:u))) / real(u-l+1,kind_real)
            else
               average = (res%box(c)%upper + res%box(c)%lower)/2.0
            endif

            res%cut_val = average
            m = select_on_coordinate_value(tp%the_data,tp%ind,c,average,l,u)
         endif

         ! moves indexes around
         res%cut_dim = c
         res%l = l
         res%u = u
!         res%cut_val = tp%the_data(c,tp%ind(m))

         res%left => build_tree_for_range(tp,l,m,res)
         res%right => build_tree_for_range(tp,m+1,u,res)

         if (associated(res%right) .eqv. .false.) then
            res%box = res%left%box
            res%cut_val_left = res%left%box(c)%upper
            res%cut_val = res%cut_val_left
         elseif (associated(res%left) .eqv. .false.) then
            res%box = res%right%box
            res%cut_val_right = res%right%box(c)%lower
            res%cut_val = res%cut_val_right
         else
            res%cut_val_right = res%right%box(c)%lower
            res%cut_val_left = res%left%box(c)%upper
            res%cut_val = (res%cut_val_left + res%cut_val_right)/2


            ! now remake the true bounding box for self.
            ! Since we are taking unions (in effect) of a tree structure,
            ! this is much faster than doing an exhaustive
            ! search over all points
            res%box%upper = max(res%left%box%upper,res%right%box%upper)
            res%box%lower = min(res%left%box%lower,res%right%box%lower)
         endif
      end if
    end function build_tree_for_range

    integer function select_on_coordinate_value(v,ind,c,alpha,li,ui) &
     result(res)
      ! Move elts of ind around between l and u, so that all points
      ! <= than alpha (in c cooordinate) are first, and then
      ! all points > alpha are second.

      !
      ! Algorithm (matt kennel).
      !
      ! Consider the list as having three parts: on the left,
      ! the points known to be <= alpha.  On the right, the points
      ! known to be > alpha, and in the middle, the currently unknown
      ! points.   The algorithm is to scan the unknown points, starting
      ! from the left, and swapping them so that they are added to
      ! the left stack or the right stack, as appropriate.
      !
      ! The algorithm finishes when the unknown stack is empty.
      !
      ! .. Scalar Arguments ..
      integer, intent (In) :: c, li, ui
      real(kind_real), intent(in) :: alpha
      ! ..
      real(kind_real) :: v(1:,1:)
      integer :: ind(1:)
      integer :: tmp
      ! ..
      integer :: lb, rb
      !
      ! The points known to be <= alpha are in
      ! [l,lb-1]
      !
      ! The points known to be > alpha are in
      ! [rb+1,u].
      !
      ! Therefore we add new points into lb or
      ! rb as appropriate.  When lb=rb
      ! we are done.  We return the location of the last point <= alpha.
      !
      !
      lb = li; rb = ui

      do while (lb < rb)
         if ( v(c,ind(lb)) <= alpha ) then
            ! it is good where it is.
            lb = lb+1
         else
            ! swap it with rb.
            tmp = ind(lb); ind(lb) = ind(rb); ind(rb) = tmp
            rb = rb-1
         endif
      end do

      ! now lb .eq. ub
      if (v(c,ind(lb)) <= alpha) then
         res = lb
      else
         res = lb-1
      endif

    end function select_on_coordinate_value

    subroutine select_on_coordinate(v,ind,c,k,li,ui)
      ! Move elts of ind around between l and u, so that the kth
      ! element
      ! is >= those below, <= those above, in the coordinate c.
      ! .. Scalar Arguments ..
      integer, intent (In) :: c, k, li, ui
      ! ..
      integer :: i, l, m, s, t, u
      ! ..
      real(kind_real) :: v(:,:)
      integer :: ind(:)
      ! ..
      l = li
      u = ui
      do while (l<u)
         t = ind(l)
         m = l
         do i = l + 1, u
            if (v(c,ind(i))<v(c,t)) then
               m = m + 1
               s = ind(m)
               ind(m) = ind(i)
               ind(i) = s
            end if
         end do
         s = ind(l)
         ind(l) = ind(m)
         ind(m) = s
         if (m<=k) l = m + 1
         if (m>=k) u = m - 1
      end do
    end subroutine select_on_coordinate

   subroutine spread_in_coordinate(tp,c,l,u,interv)
      ! the spread in coordinate 'c', between l and u.
      !
      ! Return lower bound in 'smin', and upper in 'smax',
      ! ..
      ! .. Structure Arguments ..
      type (kdtree2), pointer :: tp
      type(interval), intent(out) :: interv
      ! ..
      ! .. Scalar Arguments ..
      integer, intent (In) :: c, l, u
      ! ..
      ! .. Local Scalars ..
      real(kind_real) :: last, lmax, lmin, t, smin,smax
      integer :: i, ulocal
      ! ..
      ! .. Local Arrays ..
      real(kind_real), pointer :: v(:,:)
      integer, pointer :: ind(:)
      ! ..
      v => tp%the_data(1:,1:)
      ind => tp%ind(1:)
      smin = v(c,ind(l))
      smax = smin

      ulocal = u

      do i = l + 2, ulocal, 2
         lmin = v(c,ind(i-1))
         lmax = v(c,ind(i))
         if (lmin>lmax) then
            t = lmin
            lmin = lmax
            lmax = t
         end if
         if (smin>lmin) smin = lmin
         if (smax<lmax) smax = lmax
      end do
      if (i==ulocal+1) then
         last = v(c,ind(ulocal))
         if (smin>last) smin = last
         if (smax<last) smax = last
      end if

      interv%lower = smin
      interv%upper = smax

    end subroutine spread_in_coordinate


  subroutine kdtree2_destroy(tp)
    ! Deallocates all memory for the tree, except input data matrix
    ! .. Structure Arguments ..
    type (kdtree2), pointer :: tp
    ! ..
    call destroy_node(tp%root)

    deallocate (tp%ind)
    nullify (tp%ind)

    if (tp%rearrange) then
       deallocate(tp%rearranged_data)
       nullify(tp%rearranged_data)
    endif

    deallocate(tp)
    return

  contains
    recursive subroutine destroy_node(np)
      ! .. Structure Arguments ..
      type (tree_node), pointer :: np
      ! ..
      ! .. Intrinsic Functions ..
      intrinsic ASSOCIATED
      ! ..
      if (associated(np%left)) then
         call destroy_node(np%left)
         nullify (np%left)
      end if
      if (associated(np%right)) then
         call destroy_node(np%right)
         nullify (np%right)
      end if
      if (associated(np%box)) deallocate(np%box)
      deallocate(np)
      return

    end subroutine destroy_node

  end subroutine kdtree2_destroy

  subroutine kdtree2_n_nearest(tp,qv,nn,results)
    ! Find the 'nn' vectors in the tree nearest to 'qv' in euclidean norm
    ! returning their indexes and distances in 'indexes' and 'distances'
    ! arrays already allocated passed to this subroutine.
    type (kdtree2), pointer      :: tp
    real(kind_real), target, intent (In)    :: qv(:)
    integer, intent (In)         :: nn
    type(kdtree2_result), target :: results(:)


    sr%ballsize = huge(1.0)
    sr%qv => qv
    sr%nn = nn
    sr%nfound = 0
    sr%centeridx = -1
    sr%correltime = 0
    sr%overflow = .false.

    sr%results => results

    sr%nalloc = nn   ! will be checked

    sr%ind => tp%ind
    sr%rearrange = tp%rearrange
    if (tp%rearrange) then
       sr%Data => tp%rearranged_data
    else
       sr%Data => tp%the_data
    endif

    call validate_query_storage(nn)
    sr%pq = pq_create(results)

    call search(tp%root)

    if (tp%sort) then
       call kdtree2_sort_results(nn, results)
    endif
!    deallocate(sr%pqp)
    return
  end subroutine kdtree2_n_nearest

  subroutine validate_query_storage(n)
    !
    ! make sure we have enough storage for n
    !
    integer, intent(in) :: n

    if (size(sr%results,1) .lt. n) then
       write (*,*) 'KD_TREE_TRANS:  you did not provide enough storage for results(1:n)'
       stop
       return
    endif

    return
  end subroutine validate_query_storage

  function square_distance(iv, qv) result (res)
    ! distance between iv[1:n] and qv[1:n]
    ! .. Function Return Value ..
    ! re-implemented to improve vectorization.
    real(kind_real) :: res
    ! ..
    ! .. Array Arguments ..
    real(kind_real),intent(in) :: iv(3),qv(3)
    ! ..
    ! ..
    ! Quadratic norm
    res = sum( (iv-qv)**2 )

    ! Truncation for reproducibility
!    res = real(int(res/rth),kind_real)*rth
  end function square_distance

  function sdistance(iv, qv) result (res)
    ! spherical distance between iv[1:n] and qv[1:n]
    ! .. Function Return Value ..
    ! re-implemented to improve vectorization.
    real(kind_real) :: res, ilat, ilon, ir, qlat, qlon, qr
    ! ..
    ! .. Array Arguments ..
    real(kind_real),intent(in) :: iv(3),qv(3)
    ! ..
    ! ..
    ! .. Back to spherical coordinates
    call scoord(iv(1),iv(2),iv(3),ilat,ilon,ir)
    call scoord(qv(1),qv(2),qv(3),qlat,qlon,qr)

    ! Distance on sphere
    call sphere_dist(ilon,ilat,qlon,qlat,res)

  end function sdistance

  recursive subroutine search(node)
    !
    ! This is the innermost core routine of the kd-tree search.  Along
    ! with "process_terminal_node", it is the performance bottleneck.
    !
    ! This version uses a logically complete secondary search of
    ! "box in bounds", whether the sear
    !
    type (Tree_node), pointer          :: node
    ! ..
    type(tree_node),pointer            :: ncloser, nfarther
    !
    integer                            :: cut_dim, i
    ! ..
    real(kind_real)                               :: qval, dis
    real(kind_real)                               :: ballsize
    real(kind_real), pointer           :: qv(:)
    type(interval), pointer :: box(:)

    if ((associated(node%left) .and. associated(node%right)) .eqv. .false.) then
       ! we are on a terminal node
       if (sr%nn .eq. 0) then
          call process_terminal_node_fixedball(node)
       else
          call process_terminal_node(node)
       endif
    else
       ! we are not on a terminal node
       qv => sr%qv(1:)
       cut_dim = node%cut_dim
       qval = qv(cut_dim)

       if (qval < node%cut_val) then
          ncloser => node%left
          nfarther => node%right
          dis = (node%cut_val_right - qval)**2
!          extra = node%cut_val - qval
       else
          ncloser => node%right
          nfarther => node%left
          dis = (node%cut_val_left - qval)**2
!          extra = qval- node%cut_val_left
       endif

       if (associated(ncloser)) call search(ncloser)

       ! we may need to search the second node.
       if (associated(nfarther)) then
          ballsize = sr%ballsize
!          dis=extra**2
          if (dis <= ballsize) then
             !
             ! we do this separately as going on the first cut dimen is often
             ! a good idea.
             ! note that if extra**2 < sr%ballsize, then the next
             ! check will also be false.
             !
             box => node%box(1:)
             do i=1,3
                if (i .ne. cut_dim) then
                   dis = dis + dis2_from_bnd(qv(i),box(i)%lower,box(i)%upper)
                   if (dis > ballsize) then
                      return
                   endif
                endif
             end do

             !
             ! if we are still here then we need to search mroe.
             !
             call search(nfarther)
          endif
       endif
    end if
  end subroutine search

  real(kind_real) function dis2_from_bnd(x,amin,amax) result (res)
    real(kind_real), intent(in) :: x, amin,amax

    if (x > amax) then
       res = (x-amax)**2;
       return
    else
       if (x < amin) then
          res = (amin-x)**2;
          return
       else
          res = 0.0
          return
       endif
    endif
    return
  end function dis2_from_bnd

  subroutine process_terminal_node(node)
    !
    ! Look for actual near neighbors in 'node', and update
    ! the search results on the sr data structure.
    !
    type (tree_node), pointer          :: node
    !
    real(kind_real), pointer          :: qv(:)
    integer, pointer       :: ind(:)
    real(kind_real), pointer          :: data(:,:)
    !
    integer                :: i, indexofi, centeridx, correltime
    real(kind_real)        :: ballsize, sd, ssd, newpri
    logical                :: rearrange
    type(pq), pointer      :: pqp
    !
    ! copy values from sr to local variables
    !
    !
    ! Notice, making local pointers with an EXPLICIT lower bound
    ! seems to generate faster code.
    ! why?  I don't know.
    qv => sr%qv(1:)
    pqp => sr%pq
    ballsize = sr%ballsize
    rearrange = sr%rearrange
    ind => sr%ind(1:)
    data => sr%Data(1:,1:)
    centeridx = sr%centeridx
    correltime = sr%correltime

    !    doing_correl = (centeridx >= 0)  ! Do we have a decorrelation window?
    !    include_point = .true.    ! by default include all points
    ! search through terminal bucket.

    mainloop: do i = node%l, node%u
       if (rearrange) then
          sd = square_distance(data(:,i), qv)
          ssd = sdistance(data(:,i), qv)
          indexofi = ind(i)  ! only read it if we have not broken out
       else
          indexofi = ind(i)
          sd = square_distance(data(:,indexofi), qv)
          ssd = sdistance(data(:,indexofi), qv)
       endif
       if (sd>ballsize) cycle mainloop

       if (centeridx > 0) then ! doing correlation interval?
          if (abs(indexofi-centeridx) < correltime) cycle mainloop
       endif


       !
       ! two choices for any point.  The list so far is either undersized,
       ! or it is not.
       !
       ! If it is undersized, then add the point and its distance
       ! unconditionally.  If the point added fills up the working
       ! list then set the sr%ballsize, maximum distance bound (largest distance on
       ! list) to be that distance, instead of the initialized +infinity.
       !
       ! If the running list is full size, then compute the
       ! distance but break out immediately if it is larger
       ! than sr%ballsize, "best squared distance" (of the largest element),
       ! as it cannot be a good neighbor.
       !
       ! Once computed, compare to best_square distance.
       ! if it is smaller, then delete the previous largest
       ! element and add the new one.

       if (sr%nfound .lt. sr%nn) then
          !
          ! add this point unconditionally to fill list.
          !
          sr%nfound = sr%nfound +1
          newpri = pq_insert(pqp,sd,ssd,indexofi)
          if (sr%nfound .eq. sr%nn) ballsize = newpri
          ! we have just filled the working list.
          ! put the best square distance to the maximum value
          ! on the list, which is extractable from the PQ.

       else
          !
          ! now, if we get here,
          ! we know that the current node has a squared
          ! distance smaller than the largest one on the list, and
          ! belongs on the list.
          ! Hence we replace that with the current one.
          !
          ballsize = pq_replace_max(pqp,sd,ssd,indexofi)
       endif
    end do mainloop
    !
    ! Reset sr variables which may have changed during loop
    !
    sr%ballsize = ballsize

  end subroutine process_terminal_node

  subroutine process_terminal_node_fixedball(node)
    !
    ! Look for actual near neighbors in 'node', and update
    ! the search results on the sr data structure, i.e.
    ! save all within a fixed ball.
    !
    type (tree_node), pointer          :: node
    !
    real(kind_real), pointer          :: qv(:)
    integer, pointer       :: ind(:)
    real(kind_real), pointer          :: data(:,:)
    !
    integer                :: nfound
    integer                :: i, indexofi
    integer                :: centeridx, correltime, nn
    real(kind_real)                   :: ballsize, sd, ssd
    logical                :: rearrange

    !
    ! copy values from sr to local variables
    !
    qv => sr%qv(1:)
    ballsize = sr%ballsize
    rearrange = sr%rearrange
    ind => sr%ind(1:)
    data => sr%Data(1:,1:)
    centeridx = sr%centeridx
    correltime = sr%correltime
    nn = sr%nn ! number to search for
    nfound = sr%nfound

    ! search through terminal bucket.
    mainloop: do i = node%l, node%u

       !
       ! two choices for any point.  The list so far is either undersized,
       ! or it is not.
       !
       ! If it is undersized, then add the point and its distance
       ! unconditionally.  If the point added fills up the working
       ! list then set the sr%ballsize, maximum distance bound (largest distance on
       ! list) to be that distance, instead of the initialized +infinity.
       !
       ! If the running list is full size, then compute the
       ! distance but break out immediately if it is larger
       ! than sr%ballsize, "best squared distance" (of the largest element),
       ! as it cannot be a good neighbor.
       !
       ! Once computed, compare to best_square distance.
       ! if it is smaller, then delete the previous largest
       ! element and add the new one.

       ! which index to the point do we use?

       if (rearrange) then
          sd = square_distance(data(:,i), qv)
          ssd = sdistance(data(:,i), qv)
          indexofi = ind(i)  ! only read it if we have not broken out
       else
          indexofi = ind(i)
          sd = square_distance(data(:,indexofi), qv)
          ssd = sdistance(data(:,indexofi), qv)
       endif
       if (sd>ballsize) cycle mainloop

       if (centeridx > 0) then ! doing correlation interval?
          if (abs(indexofi-centeridx)<correltime) cycle mainloop
       endif

       nfound = nfound+1
       if (nfound .gt. sr%nalloc) then
          ! oh nuts, we have to add another one to the tree but
          ! there isn't enough room.
          sr%overflow = .true.
       else
          sr%results(nfound)%dis = sd
          sr%results(nfound)%sdis = ssd
          sr%results(nfound)%idx = indexofi
       endif
    end do mainloop
    !
    ! Reset sr variables which may have changed during loop
    !
    sr%nfound = nfound
  end subroutine process_terminal_node_fixedball

  subroutine kdtree2_sort_results(nfound,results)
    !  Use after search to sort results(1:nfound) in order of increasing
    !  distance.
    integer, intent(in)          :: nfound
    type(kdtree2_result), target :: results(:)
    !
    !

    !THIS IS BUGGY WITH INTEL FORTRAN
    !    If (nfound .Gt. 1) Call heapsort(results(1:nfound)%dis,results(1:nfound)%ind,nfound)
    !
    if (nfound .gt. 1) call heapsort_struct(results,nfound)

    return
  end subroutine kdtree2_sort_results

  subroutine heapsort_struct(a,n)
    !
    ! Sort a(1:n) in ascending order
    !
    !
    integer,intent(in)                 :: n
    type(kdtree2_result),intent(inout) :: a(:)

    !
    !
    type(kdtree2_result) :: value ! temporary value

    integer     :: i,j
    integer     :: ileft,iright

    ileft=n/2+1
    iright=n

    !    do i=1,n
    !       ind(i)=i
    ! Generate initial idum array
    !    end do

    if(n.eq.1) return

    do
       if(ileft > 1)then
          ileft=ileft-1
          value=a(ileft)
       else
          value=a(iright)
          a(iright)=a(1)
          iright=iright-1
          if (iright == 1) then
             a(1) = value
             return
          endif
       endif
       i=ileft
       j=2*ileft
       do while (j <= iright)
          if(j < iright) then
             if(a(j)%dis < a(j+1)%dis) j=j+1
          endif
          if(value%dis < a(j)%dis) then
             a(i)=a(j);
             i=j
             j=j+j
          else
             j=iright+1
          endif
       end do
       a(i)=value
    end do
  end subroutine heapsort_struct

end module tools_kdtree2
