!----------------------------------------------------------------------
! Module: type_ctree
!> Purpose: cover tree derived type
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module type_ctree

use iso_c_binding, only: c_ptr,c_int,c_double
use tools_kinds, only: kind_real

implicit none

type ctreetype
    type(c_ptr) :: ptr !< Pointer to the C++ class
end type ctreetype

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
public :: ctreetype,create_ctree,delete_ctree,find_nearest_neighbors

contains

!----------------------------------------------------------------------
! Subroutine: create_ctree
!> Purpose: create a cover tree
!----------------------------------------------------------------------
function create_ctree(n,lon,lat,mask)

implicit none

! Passed variables
type(ctreetype) :: create_ctree      !< Cover tree object
integer,intent(in) :: n              !< Number of points
real(kind_real),intent(in) :: lon(n) !< Points longitudes
real(kind_real),intent(in) :: lat(n) !< Points latitudes
integer,intent(in) :: mask(n)        !< Mask

! Call C++ function
create_ctree%ptr = create_ctree_c(n,lon,lat,mask)

end function create_ctree

!----------------------------------------------------------------------
! Subroutine: delete_ctree
!> Purpose: delete a cover tree
!----------------------------------------------------------------------
subroutine delete_ctree(this)

implicit none

! Passed variables
type(ctreetype),intent(inout) :: this !< Cover tree object

! Call C++ function
call delete_ctree_c(this%ptr)

end subroutine delete_ctree

!----------------------------------------------------------------------
! Subroutine: find_nearest_neighbors
!> Purpose: find nearest neighbors using a cover tree
!----------------------------------------------------------------------
subroutine find_nearest_neighbors(this,lon,lat,nn,nn_index,nn_dist)

implicit none

! Passed variables
class(ctreetype),intent(in) :: this        !< Cover tree object
real(kind_real),intent(in) :: lon          !< Point longitude
real(kind_real),intent(in) :: lat          !< Point latitude
integer,intent(in) :: nn                   !< Number of nearest neighbors to find
integer,intent(out) :: nn_index(nn)        !< Neareast neighbors index
real(kind_real),intent(out) :: nn_dist(nn) !< Neareast neighbors distance

! Call C++ function
call find_nearest_neighbors_c(this%ptr,lon,lat,nn,nn_index,nn_dist)

end subroutine find_nearest_neighbors

end module type_ctree
