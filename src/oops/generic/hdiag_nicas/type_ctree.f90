!----------------------------------------------------------------------
! Module: type_ctree
!> Purpose: cover tree derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: ctree code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_ctree

use iso_c_binding, only: c_ptr,c_int,c_double
use tools_display, only: msgerror
use tools_kinds, only: kind_real
use tools_qsort, only: qsort

implicit none

type ctree_type
    type(c_ptr) :: ptr !< Pointer to the C++ class
contains
    procedure :: create => ctree_create
    procedure :: delete => ctree_delete
    procedure :: find_nearest_neighbors
end type ctree_type

real(kind_real),parameter :: rth = 1.0e-12 !< Reproducibility threshold

! C++ interface
interface
   function create_ctree_c(n,lon,lat,mask) bind(C,name="create_ctree")
   use iso_c_binding
   implicit none
   type(c_ptr) :: create_ctree_c
   integer(c_int),value :: n
   real(c_double) :: lon(*)
   real(c_double) :: lat(*)
   integer(c_int) :: mask(*)
   end function create_ctree_c
end interface
interface
   subroutine delete_ctree_c(ctree) bind(C,name="delete_ctree")
   use iso_c_binding
   implicit none
   type(c_ptr),value :: ctree
   end subroutine delete_ctree_c
end interface
interface
   subroutine find_nearest_neighbors_c(ctree,lon,lat,nn,nn_index,nn_dist) bind(C,name="find_nearest_neighbors")
   use iso_c_binding
   implicit none
   type(c_ptr),value :: ctree
   real(c_double),value :: lon
   real(c_double),value :: lat
   integer(c_int),value :: nn
   integer(c_int) :: nn_index(*)
   real(c_double) :: nn_dist(*)
   end subroutine find_nearest_neighbors_c
end interface

private
public :: ctree_type

contains

!----------------------------------------------------------------------
! Subroutine: ctree_create
!> Purpose: create a cover tree
!----------------------------------------------------------------------
subroutine ctree_create(ctree,n,lon,lat,mask)

implicit none

! Passed variables
class(ctree_type) :: ctree           !< Cover tree object
integer,intent(in) :: n              !< Number of points
real(kind_real),intent(in) :: lon(n) !< Points longitudes
real(kind_real),intent(in) :: lat(n) !< Points latitudes
logical,intent(in) :: mask(n)        !< Mask

! Local variable
integer :: i,imask(n)

! Check mask
if (count(mask)<1) call msgerror('mask should have at least one valid point to create a ctree')

! Convert mask
do i=1,n
   if (mask(i)) then
      imask(i) = 1
   else
      imask(i) = 0
   end if
end do

! Call C++ function
ctree%ptr = create_ctree_c(n,lon,lat,imask)

end subroutine ctree_create

!----------------------------------------------------------------------
! Subroutine: ctree_delete
!> Purpose: delete a cover tree
!----------------------------------------------------------------------
subroutine ctree_delete(ctree)

implicit none

! Passed variables
class(ctree_type),intent(inout) :: ctree !< Cover tree object

! Call C++ function
call delete_ctree_c(ctree%ptr)

end subroutine ctree_delete

!----------------------------------------------------------------------
! Subroutine: find_nearest_neighbors
!> Purpose: find nearest neighbors using a cover tree
!----------------------------------------------------------------------
subroutine find_nearest_neighbors(ctree,lon,lat,nn,nn_index,nn_dist)

implicit none

! Passed variables
class(ctree_type),intent(in) :: ctree      !< Cover tree object
real(kind_real),intent(in) :: lon          !< Point longitude
real(kind_real),intent(in) :: lat          !< Point latitude
integer,intent(in) :: nn                   !< Number of nearest neighbors to find
integer,intent(out) :: nn_index(nn)        !< Neareast neighbors index
real(kind_real),intent(out) :: nn_dist(nn) !< Neareast neighbors distance

! Local variables
integer :: i,j,nid
integer,allocatable :: order(:)

! Call C++ function
call find_nearest_neighbors_c(ctree%ptr,lon,lat,nn,nn_index,nn_dist)

! Indistinguishability threshold for cross-plateform reproducibility
i = 1
do while (i<nn)
   ! Count indistinguishable neighbors
   nid = 1
   do j=i+1,nn
      if (abs(nn_dist(i)-nn_dist(j))<rth*nn_dist(i)) nid = nid+1
   end do

   ! Reorder
   if (nid>1) then
      allocate(order(nid))
      call qsort(nid,nn_index(i:i+nid-1),order)
      do j=1,nid
         nn_dist(i+j-1) = nn_dist(i+order(j)-1)
      end do
      deallocate(order)
   end if

   ! Update
   i = i+nid
end do

end subroutine find_nearest_neighbors

end module type_ctree
