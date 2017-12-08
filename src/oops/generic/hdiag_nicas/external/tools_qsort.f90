!----------------------------------------------------------------------
! Module: tools_qsort.f90
!> Purpose: qsort routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module tools_qsort

use tools_kinds, only: kind_real

implicit none

interface qsort
  module procedure qsort_integer
  module procedure qsort_real
end interface

interface quick_sort
  module procedure quick_sort_integer
  module procedure quick_sort_real
end interface

interface interchange_sort
  module procedure interchange_sort_integer
  module procedure interchange_sort_real
end interface

private
public :: qsort

contains

!----------------------------------------------------------------------
! qsort_integer
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which gives the positions of the elements in the original order.
! Source: http://jblevins.org/mirror/amiller/qsort.f90
! Modified by Benjamin Menetrier for nicas
!----------------------------------------------------------------------
recursive subroutine qsort_integer(n,list,order)

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

end subroutine qsort_integer

!----------------------------------------------------------------------
! qsort_real
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which gives the positions of the elements in the original order.
! Source: http://jblevins.org/mirror/amiller/qsort.f90
! Modified by Benjamin Menetrier for hybrid_diag
!----------------------------------------------------------------------
recursive subroutine qsort_real(n,list,order)

implicit none

! Passed variables
integer, intent(in) :: n                 !< Input vector size
real(kind_real),intent(inout) :: list(n) !< Vector to sort
integer,intent(inout) :: order(n)        !< Positions of the elements in the original order

! Local variable
integer :: i

do i=1,n
  order(i) = i
end do

call quick_sort(n,1,n,list,order)

end subroutine qsort_real

!----------------------------------------------------------------------
! Subroutine: quick_sort_integer
!> Purpose: sort an integer subvector
!----------------------------------------------------------------------
recursive subroutine quick_sort_integer(n,left_end,right_end,list,order)

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

end subroutine quick_sort_integer

!----------------------------------------------------------------------
! Subroutine: quick_sort_real
!> Purpose: sort an real subvector
!----------------------------------------------------------------------
recursive subroutine quick_sort_real(n,left_end,right_end,list,order)

implicit none

! Passed variables
integer,intent(in) :: n                  !< Input vector size
integer,intent(in) :: left_end           !< Left end of the vector
integer,intent(in) :: right_end          !< Right end of the vector
real(kind_real),intent(inout) :: list(n) !< Vector to sort
integer,intent(inout) :: order(n)        !< Positions of the elements in the original order

! Local variables
integer,parameter :: max_simple_sort_size = 6
integer :: i,j,itemp
real(kind_real) :: reference,temp

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

end subroutine quick_sort_real

!----------------------------------------------------------------------
! Subroutine: interchange_sort_integer
!> Purpose: interchange integers
!----------------------------------------------------------------------
subroutine interchange_sort_integer(n,left_end,right_end,list,order)

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

end subroutine interchange_sort_integer

!----------------------------------------------------------------------
! Subroutine: interchange_sort_real
!> Purpose: interchange reals
!----------------------------------------------------------------------
subroutine interchange_sort_real(n,left_end,right_end,list,order)

implicit none

! Passed variables
integer,intent(in) :: n                  !< Input vector size
integer,intent(in) :: left_end           !< Left end of the vector
integer,intent(in) :: right_end          !< Right end of the vector
real(kind_real),intent(inout) :: list(n) !< Vector to sort
integer,intent(inout) :: order(n)        !< Positions of the elements in the original order

! Local variables
integer :: i,j,itemp
real(kind_real) :: temp

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

end subroutine interchange_sort_real

end module tools_qsort
