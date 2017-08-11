!----------------------------------------------------------------------
! qsort
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which gives the positions of the elements in the original order.
! Source: http://jblevins.org/mirror/amiller/qsort.f90
! Modified by Benjamin Menetrier for nicas
!----------------------------------------------------------------------
recursive subroutine qsort(n,list,order)

implicit none

! Passed variables
integer, intent(in) :: n          !< Input vector size
integer,intent(inout) :: list(n)  !< Vector to sort
integer,intent(inout) :: order(n) !< Positions of the elements in the original order

! Local variable
integer :: i

do i=1,n
  order(i) = i
end do

call quick_sort(n,1,n,list,order)

end subroutine qsort

recursive subroutine quick_sort(n,left_end,right_end,list,order)

implicit none

! Passed variables
integer,intent(in) :: n           !< Input vector size
integer,intent(in) :: left_end    !< Left end of the vector
integer,intent(in) :: right_end   !< Right end of the vector
integer,intent(inout) :: list(n)  !< Vector to sort
integer,intent(inout) :: order(n) !< Positions of the elements in the original order

! Local variables
integer,parameter :: max_simple_sort_size = 6
integer :: i,j,itemp
integer :: reference,temp

if (right_end<left_end+max_simple_sort_size) then
  ! Use interchange sort for small lists
  call interchange_sort(n,left_end,right_end,list,order)
else
  ! Use partition ("quick") sort
  reference = list((left_end+right_end)/2)
  i = left_end-1
  j = right_end+1
  do
    ! Scan list from left end until element >= reference is found
    do
      i = i+1
      if (list(i)>=reference) exit
    end do
    ! Scan list from right end until element <= reference is found
    do
      j = j-1
      if (list(j)<=reference) exit
    end do

    if (i<j) then
      ! Swap two out-of-order elements
      temp = list(i)
      list(i) = list(j)
      list(j) = temp
      itemp = order(i)
      order(i) = order(j)
      order(j) = itemp
    elseif (i==j) then
      i = i+1
      exit
    else
      exit
    end if
  end do

  if (left_end<j) call quick_sort(n,left_end,j,list,order)
  if (i<right_end) call quick_sort(n,i,right_end,list,order)
end if

end subroutine quick_sort

!----------------------------------------------------------------------
! Subroutine: quick_sort
!> Purpose: sort a subvector
!----------------------------------------------------------------------
subroutine interchange_sort(n,left_end,right_end,list,order)

implicit none

! Passed variables
integer,intent(in) :: n           !< Input vector size
integer,intent(in) :: left_end    !< Left end of the vector
integer,intent(in) :: right_end   !< Right end of the vector
integer,intent(inout) :: list(n)  !< Vector to sort
integer,intent(inout) :: order(n) !< Positions of the elements in the original order

! Local variables
integer :: i,j,itemp
integer :: temp

do i=left_end,right_end-1
  do j=i+1,right_end
    if (list(i)>list(j)) then
      temp = list(i)
      list(i) = list(j)
      list(j) = temp
      itemp = order(i)
      order(i) = order(j)
      order(j) = itemp
    end if
  end do
end do

end subroutine interchange_sort
