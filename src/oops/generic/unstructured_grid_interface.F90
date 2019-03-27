! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module unstructured_grid_interface

use iso_c_binding
use unstructured_grid_mod

implicit none

private
!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
!> Create unstructured grid
subroutine create_ug_c(key,colocated,nts) bind(c,name='create_ug_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: key    !< Unstructured grid
integer(c_int),intent(in) :: colocated !< Colocation flag
integer(c_int),intent(in) :: nts       !< Number of timeslots

! Local variables
type(unstructured_grid),pointer :: self

! Interface
call unstructured_grid_registry%init()
call unstructured_grid_registry%add(key)
call unstructured_grid_registry%get(key,self)

! Call Fortran
call create_ug(self,colocated,nts)

end subroutine create_ug_c
! ------------------------------------------------------------------------------
!> Delete unstructured grid
subroutine delete_ug_c(key) bind(c,name='delete_ug_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: key !< Unstructured grid

! Local variables
type(unstructured_grid),pointer :: self

! Interface
call unstructured_grid_registry%get(key,self)

! Call Fortran
call delete_ug(self)

! Clear interface
call unstructured_grid_registry%remove(key)

end subroutine delete_ug_c
! ------------------------------------------------------------------------------
end module unstructured_grid_interface
