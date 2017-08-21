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
  integer :: nvar3d                             !> Number of 3d variables
  integer :: nvar2d                             !> Number of 2d variables
  integer, allocatable :: mask3d(:)             !> 3d mask (size nlevs)
  integer :: mask2d                             !> 2d mask
  real(kind=kind_real), allocatable :: fld3d(:) !> 3d fields values (size nlevs * nvar3d)
  real(kind=kind_real), allocatable :: fld2d(:) !> 2d fields values (size nfld2d)
  integer :: glbind                             !> global index (optional)
end type column_data

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

subroutine create_column_data(self, plat, plon, parea, klevs, kvar3d, kvar2d, kmask3d, kmask2d, kglbind)
implicit none
type(column_data), intent(inout) :: self
real(kind=kind_real), intent(in) :: plat
real(kind=kind_real), intent(in) :: plon
real(kind=kind_real), intent(in) :: parea
integer, intent(in) :: klevs
integer, intent(in) :: kvar3d
integer, intent(in) :: kvar2d
integer, intent(in) :: kmask3d(klevs)
integer, intent(in) :: kmask2d
integer, intent(in) :: kglbind
self%lat = plat
self%lon = plon
self%area = parea
self%nlevs = klevs
self%nvar3d = kvar3d
self%nvar2d = kvar2d
allocate(self%mask3d(self%nlevs))
self%mask3d = kmask3d
self%mask2d = kmask2d
allocate(self%fld3d(self%nlevs*self%nvar3d))
allocate(self%fld2d(self%nvar2d))
self%glbind = kglbind

end subroutine create_column_data

!-------------------------------------------------------------------------------

subroutine delete_column_data(self)
implicit none
type(column_data), intent(inout) :: self

deallocate(self%fld3d)
deallocate(self%fld2d)

end subroutine delete_column_data

!-------------------------------------------------------------------------------

end module column_data_mod
