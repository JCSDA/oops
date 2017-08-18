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
  real(kind=kind_real) :: lat  !> Latitude of column
  real(kind=kind_real) :: lon  !> Longitude of column
  real(kind=kind_real) :: area !> Area of column cell
  integer :: nlevs             !> Number of levels
  integer :: nvars             !> Number of variables in 3D fields
  integer :: nsurf             !> Number of variables in 2D fields
  integer, allocatable :: cmask(:)              !> levs mask (size nlevs)
  integer :: smask                              !> surface mask
  real(kind=kind_real), allocatable :: cols(:)  !> column of values (size nlevs * nvars)
  real(kind=kind_real), allocatable :: surf(:)  !> surface fields   (size nsurf)
end type column_data

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

subroutine create_column_data(self, plat, plon, parea, klevs, kvars, ksurf, kcmask, ksmask)
implicit none
type(column_data), intent(inout) :: self
real(kind=kind_real), intent(in) :: plat
real(kind=kind_real), intent(in) :: plon
real(kind=kind_real), intent(in) :: parea
integer, intent(in) :: klevs
integer, intent(in) :: kvars
integer, intent(in) :: ksurf
integer, intent(in) :: kcmask(klevs)
integer, intent(in) :: ksmask
self%lat = plat
self%lon = plon
self%area = parea
self%nlevs = klevs
self%nvars = kvars
self%nsurf = ksurf
allocate(self%cmask(self%nlevs))
self%cmask = kcmask
self%smask = ksmask
allocate(self%cols(self%nlevs*self%nvars))
allocate(self%surf(self%nsurf))

end subroutine create_column_data

!-------------------------------------------------------------------------------

subroutine delete_column_data(self)
implicit none
type(column_data), intent(inout) :: self

deallocate(self%cols)
deallocate(self%surf)

end subroutine delete_column_data

!-------------------------------------------------------------------------------

end module column_data_mod
