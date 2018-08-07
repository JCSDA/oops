!
!(c) Matthew Kennel, Institute for Nonlinear Science (2004)
!
! Licensed under the Academic Free License version 1.1 found in file LICENSE
! with additional provisions found in that same file.
!
module tools_kdtree2_pq
  use tools_func, only: inf,sup,supeq
  use tools_kinds, only: kind_real
  !
  ! maintain a priority queue (PQ) of data, pairs of 'priority/payload',
  ! implemented with a binary heap.  This is the type, and the 'dis' field
  ! is the priority.
  !
  type kdtree2_result
      ! a pair of distances, indexes
      real(kind_real)    :: dis!=0.0
      real(kind_real)    :: sdis!=0.0
      integer :: idx!=-1   Initializers cause some bugs in compilers.
  end type kdtree2_result
  !
  ! A heap-based priority queue lets one efficiently implement the following
  ! operations, each in log(N) time, as opposed to linear time.
  !
  ! 1)  add a datum (push a datum onto the queue, increasing its length)
  ! 2)  return the priority value of the maximum priority element
  ! 3)  pop-off (and delete) the element with the maximum priority, decreasing
  !     the size of the queue.
  ! 4)  replace the datum with the maximum priority with a supplied datum
  !     (of either higher or lower priority), maintaining the size of the
  !     queue.
  !
  !
  ! In the k-d tree case, the 'priority' is the square distance of a point in
  ! the data set to a reference point.   The goal is to keep the smallest M
  ! distances to a reference point.  The tree algorithm searches terminal
  ! nodes to decide whether to add points under consideration.
  !
  ! A priority queue is useful here because it lets one quickly return the
  ! largest distance currently existing in the list.  If a new candidate
  ! distance is smaller than this, then the new candidate ought to replace
  ! the old candidate.  In priority queue terms, this means removing the
  ! highest priority element, and inserting the new one.
  !
  ! Algorithms based on Cormen, Leiserson, Rivest, _Introduction
  ! to Algorithms_, 1990, with further optimization by the author.
  !
  ! Originally informed by a C implementation by Sriranga Veeraraghavan.
  !
  ! This module is not written in the most clear way, but is implemented such
  ! for speed, as it its operations will be called many times during searches
  ! of large numbers of neighbors.
  !
  type pq
      !
      ! The priority queue consists of elements
      ! priority(1:heap_size), with associated payload(:).
      !
      ! There are heap_size active elements.
      ! Assumes the allocation is always sufficient.  Will NOT increase it
      ! to match.
      integer :: heap_size = 0
      type(kdtree2_result), pointer :: elems(:)
  end type pq

  public :: kdtree2_result

  public :: pq
  public :: pq_create
  public :: pq_delete, pq_insert
  public :: pq_extract_max, pq_max, pq_replace_max, pq_maxpri
  private

contains


  function pq_create(results_in) result(res)
    !
    ! Create a priority queue from ALREADY allocated
    ! array pointers for storage.  NOTE! It will NOT
    ! add any alements to the heap, i.e. any existing
    ! data in the input arrays will NOT be used and may
    ! be overwritten.
    !
    ! usage:
    !    real(kind_real), pointer :: x(:)
    !    integer, pointer :: k(:)
    !    allocate(x(1000),k(1000))
    !    pq => pq_create(x,k)
    !
    type(kdtree2_result), target:: results_in(:)
    type(pq) :: res
    !
    !
    integer :: nalloc

    nalloc = size(results_in,1)
    if (nalloc .lt. 1) then
       write (*,*) 'PQ_CREATE: error, input arrays must be allocated.'
    end if
    res%elems => results_in
    res%heap_size = 0
    return
  end function pq_create

  !
  ! operations for getting parents and left + right children
  ! of elements in a binary heap.
  !

!
! These are written inline for speed.
!
!  integer function parent(i)
!    integer, intent(in) :: i
!    parent = (i/2)
!    return
!  end function parent

!  integer function left(i)
!    integer, intent(in) ::i
!    left = (2*i)
!    return
!  end function left

!  integer function right(i)
!    integer, intent(in) :: i
!    right = (2*i)+1
!    return
!  end function right

!  logical function compare_priority(p1,p2)
!    real(kind_real), intent(in) :: p1, p2
!
!    compare_priority = (p1 .gt. p2)
!    return
!  end function compare_priority

  subroutine heapify(a,i_in)
    !
    ! take a heap rooted at 'i' and force it to be in the
    ! heap canonical form.   This is performance critical
    ! and has been tweaked a little to reflect this.
    !
    type(pq),pointer   :: a
    integer, intent(in) :: i_in
    !
    integer :: i, l, r, largest

    real(kind_real)    :: pri_i, pri_l, pri_r, pri_largest


    type(kdtree2_result) :: temp

    i = i_in

bigloop:  do
       l = 2*i ! left(i)
       r = l+1 ! right(i)
       !
       ! set 'largest' to the index of either i, l, r
       ! depending on whose priority is largest.
       !
       ! note that l or r can be larger than the heap size
       ! in which case they do not count.


       ! does left child have higher priority?
       if (l .gt. a%heap_size) then
          ! we know that i is the largest as both l and r are invalid.
          exit
       else
          pri_i = a%elems(i)%dis
          pri_l = a%elems(l)%dis
          if (sup(pri_l,pri_i)) then
             largest = l
             pri_largest = pri_l
          else
             largest = i
             pri_largest = pri_i
          endif

          !
          ! between i and l we have a winner
          ! now choose between that and r.
          !
          if (r .le. a%heap_size) then
             pri_r = a%elems(r)%dis
             if (sup(pri_r,pri_largest)) then
                largest = r
             endif
          endif
       endif

       if (largest .ne. i) then
          ! swap data in nodes largest and i, then heapify

          temp = a%elems(i)
          a%elems(i) = a%elems(largest)
          a%elems(largest) = temp
          !
          ! Canonical heapify() algorithm has tail-ecursive call:
          !
          !        call heapify(a,largest)
          ! we will simulate with cycle
          !
          i = largest
          cycle bigloop ! continue the loop
       else
          return   ! break from the loop
       end if
    enddo bigloop
    return
  end subroutine heapify

  subroutine pq_max(a,e)
    !
    ! return the priority and its payload of the maximum priority element
    ! on the queue, which should be the first one, if it is
    ! in heapified form.
    !
    type(pq),pointer :: a
    type(kdtree2_result),intent(out)  :: e

    if (a%heap_size .gt. 0) then
       e = a%elems(1)
    else
       write (*,*) 'PQ_MAX: ERROR, heap_size < 1'
       stop
    endif
    return
  end subroutine pq_max

  real(kind_real) function pq_maxpri(a)
    type(pq), pointer :: a

    if (a%heap_size .gt. 0) then
       pq_maxpri = a%elems(1)%dis
    else
       write (*,*) 'PQ_MAX_PRI: ERROR, heapsize < 1'
       stop
    endif
    return
  end function pq_maxpri

  subroutine pq_extract_max(a,e)
    !
    ! return the priority and payload of maximum priority
    ! element, and remove it from the queue.
    ! (equivalent to 'pop()' on a stack)
    !
    type(pq),pointer :: a
    type(kdtree2_result), intent(out) :: e

    if (a%heap_size .ge. 1) then
       !
       ! return max as first element
       !
       e = a%elems(1)

       !
       ! move last element to first
       !
       a%elems(1) = a%elems(a%heap_size)
       a%heap_size = a%heap_size-1
       call heapify(a,1)
       return
    else
       write (*,*) 'PQ_EXTRACT_MAX: error, attempted to pop non-positive PQ'
       stop
    end if

  end subroutine pq_extract_max


  real(kind_real) function pq_insert(a,dis,sdis,idx)
    !
    ! Insert a new element and return the new maximum priority,
    ! which may or may not be the same as the old maximum priority.
    !
    type(pq),pointer  :: a
    real(kind_real), intent(in) :: dis
    real(kind_real), intent(in) :: sdis
    integer, intent(in) :: idx
    !    type(kdtree2_result), intent(in) :: e
    !
    integer :: i, isparent
    real(kind_real)    :: parentdis,parentsdis
    !

    !    if (a%heap_size .ge. a%max_elems) then
    !       write (*,*) 'PQ_INSERT: error, attempt made to insert element on full PQ'
    !       stop
    !    else
    a%heap_size = a%heap_size + 1
    i = a%heap_size

    do while (i .gt. 1)
       isparent = int(i/2)
       parentdis = a%elems(isparent)%dis
       parentsdis = a%elems(isparent)%sdis
       if (sup(dis,parentdis)) then
          ! move what was in i's parent into i.
          a%elems(i)%dis = parentdis
          a%elems(i)%sdis = parentsdis
          a%elems(i)%idx = a%elems(isparent)%idx
          i = isparent
       else
          exit
       endif
    end do

    ! insert the element at the determined position
    a%elems(i)%dis = dis
    a%elems(i)%sdis = sdis
    a%elems(i)%idx = idx

    pq_insert = a%elems(1)%dis
    return
    !    end if

  end function pq_insert

  real(kind_real) function pq_replace_max(a,dis,sdis,idx)
    !
    ! Replace the extant maximum priority element
    ! in the PQ with (dis,sdis,idx).  Return
    ! the new maximum priority, which may be larger
    ! or smaller than the old one.
    !
    type(pq),pointer         :: a
    real(kind_real), intent(in) :: dis
    real(kind_real), intent(in) :: sdis
    integer, intent(in) :: idx
!    type(kdtree2_result), intent(in) :: e
    ! not tested as well!

    integer :: parent, child, N
    real(kind_real)    :: prichild, prichildp1

    type(kdtree2_result) :: etmp

    if (.true.) then
       N=a%heap_size
       if (N .ge. 1) then
          parent =1
          child=2

          loop: do while (child .le. N)
             prichild = a%elems(child)%dis

             !
             ! posibly child+1 has higher priority, and if
             ! so, get it, and increment child.
             !

             if (child .lt. N) then
                prichildp1 = a%elems(child+1)%dis
                if (inf(prichild,prichildp1)) then
                   child = child+1
                   prichild = prichildp1
                endif
             endif

             if (supeq(dis,prichild)) then
                exit loop
                ! we have a proper place for our new element,
                ! bigger than either children's priority.
             else
                ! move child into parent.
                a%elems(parent) = a%elems(child)
                parent = child
                child = 2*parent
             end if
          end do loop
          a%elems(parent)%dis = dis
          a%elems(parent)%sdis = sdis
          a%elems(parent)%idx = idx
          pq_replace_max = a%elems(1)%dis
       else
          a%elems(1)%dis = dis
          a%elems(1)%sdis = sdis
          a%elems(1)%idx = idx
          pq_replace_max = dis
       endif
    else
       !
       ! slower version using elementary pop and push operations.
       !
       call pq_extract_max(a,etmp)
       etmp%dis = dis
       etmp%sdis = sdis
       etmp%idx = idx
       pq_replace_max = pq_insert(a,dis,sdis,idx)
    endif
    return
  end function pq_replace_max

  subroutine pq_delete(a,i)
    !
    ! delete item with index 'i'
    !
    type(pq),pointer :: a
    integer           :: i

    if ((i .lt. 1) .or. (i .gt. a%heap_size)) then
       write (*,*) 'PQ_DELETE: error, attempt to remove out of bounds element.'
       stop
    endif

    ! swap the item to be deleted with the last element
    ! and shorten heap by one.
    a%elems(i) = a%elems(a%heap_size)
    a%heap_size = a%heap_size - 1

    call heapify(a,i)

  end subroutine pq_delete

end module tools_kdtree2_pq
