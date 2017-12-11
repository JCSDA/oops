! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!>  Fortran module for handling generic column of data

module column_data_mod

use kinds

implicit none
private
public column_data, create_column_data, delete_column_data

!>  Derived type containing the data

type column_data
  real(kind=kind_real) :: lat                   !> Latitude of column
  real(kind=kind_real) :: lon                   !> Longitude of column
  real(kind=kind_real) :: area                  !> Area of column cell
  integer :: nlevs                              !> Number of levels
  integer :: nvar                               !> Number of  variables
  integer, allocatable :: mask(:)               !> Mask (size nlevs)
  real(kind=kind_real), allocatable :: fld(:,:) !> Fields values (size nlevs * nvar)
  integer :: glbind                             !> global index (optional)
end type column_data

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

subroutine create_column_data(self, plat, plon, parea, klevs, kvar, kmask, kglbind)
implicit none
type(column_data), intent(inout) :: self
real(kind=kind_real), intent(in) :: plat
real(kind=kind_real), intent(in) :: plon
real(kind=kind_real), intent(in) :: parea
integer, intent(in) :: klevs
integer, intent(in) :: kvar
integer, intent(in) :: kmask(klevs)
integer, intent(in) :: kglbind
self%lat = plat
self%lon = plon
self%area = parea
self%nlevs = klevs
self%nvar = kvar
allocate(self%mask(self%nlevs))
self%mask = kmask
allocate(self%fld(self%nlevs,self%nvar))
self%glbind = kglbind

end subroutine create_column_data

!-------------------------------------------------------------------------------

subroutine delete_column_data(self)
implicit none
type(column_data), intent(inout) :: self

deallocate(self%fld)

end subroutine delete_column_data

!-------------------------------------------------------------------------------

end module column_data_mod
