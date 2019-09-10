! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module oobump_interface

use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds
use missing_values_mod
use oobump_mod
use type_bump, only: bump_type
use type_nam, only: nvmax,nlmax,nc3max,nscalesmax,ndirmax,nldwvmax 
use unstructured_grid_mod

implicit none

private
! ------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
!> Create OOBUMP
subroutine create_oobump_c(c_key_oobump, c_key_ug, c_conf, ens1_ne, ens1_nsub, ens2_ne, ens2_nsub) bind(c, name='create_oobump_f90')

implicit none

! Passed variables
integer(c_int), intent(inout) :: c_key_oobump !< OOBUMP
integer(c_int), intent(in) :: c_key_ug        !< Unstructured grid
type(c_ptr), intent(in) :: c_conf             !< Configuration
integer(c_int), intent(in) :: ens1_ne         !< First ensemble size
integer(c_int), intent(in) :: ens1_nsub       !< Number of sub-ensembles in the first ensemble
integer(c_int), intent(in) :: ens2_ne         !< Second ensemble size
integer(c_int), intent(in) :: ens2_nsub       !< Number of sub-ensembles in the second ensemble

! Local variables
type(fckit_configuration) :: f_conf
type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug

! Interface
f_conf = fckit_configuration(c_conf)
call oobump_registry%init()
call oobump_registry%add(c_key_oobump)
call oobump_registry%get(c_key_oobump, self)
call unstructured_grid_registry%get(c_key_ug, ug)

! Call Fortran
call create_oobump(self, ug, f_conf, ens1_ne, ens1_nsub, ens2_ne, ens2_nsub)

end subroutine create_oobump_c
! ------------------------------------------------------------------------------
!> Delete OOBUMP
subroutine delete_oobump_c(c_key_oobump) bind(c, name='delete_oobump_f90')

implicit none

! Passed variables
integer(c_int), intent(inout) :: c_key_oobump !< OOBUMP

! Local variables
type(oobump_type), pointer :: self

! Interface
call oobump_registry%get(c_key_oobump, self)

! Delete OOBUMP
call delete_oobump(self)

! Clean interface
call oobump_registry%remove(c_key_oobump)

end subroutine delete_oobump_c
! ------------------------------------------------------------------------------
!> Get colocated flag
subroutine get_oobump_colocated_c(c_key_oobump, colocated) bind(c, name='get_oobump_colocated_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump !< OOBUMP
integer(c_int), intent(out) :: colocated   !< Colocated flag

! Local variables
type(oobump_type), pointer :: self

! Interface
call oobump_registry%get(c_key_oobump, self)

! Call Fortran
call get_oobump_colocated(self, colocated)

end subroutine get_oobump_colocated_c
! ------------------------------------------------------------------------------
!> Get number of timeslots
subroutine get_oobump_nts_c(c_key_oobump, nts) bind(c, name='get_oobump_nts_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump !< OOBUMP
integer(c_int), intent(out) :: nts         !< Number of timeslots

! Local variables
type(oobump_type), pointer :: self

! Interface
call oobump_registry%get(c_key_oobump, self)

! Call Fortran
call get_oobump_nts(self, nts)

end subroutine get_oobump_nts_c
! ------------------------------------------------------------------------------
!> Get control variable size
subroutine get_oobump_cv_size_c(c_key_oobump, n) bind(c, name='get_oobump_cv_size_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump !< OOBUMP
integer(c_int), intent(out) :: n           !< Control variable size

! Local variables
type(oobump_type), pointer :: self

! Interface
call oobump_registry%get(c_key_oobump, self)

! Call Fortran
call get_oobump_cv_size(self, n)

end subroutine get_oobump_cv_size_c
! ------------------------------------------------------------------------------
!> Add ensemble member
subroutine add_oobump_member_c(c_key_oobump, c_key_ug, ie, iens) bind(c, name='add_oobump_member_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump !< OOBUMP
integer(c_int), intent(in) :: c_key_ug     !< Unstructured grid
integer(c_int), intent(in) :: ie           !< Ensemble member index
integer(c_int), intent(in) :: iens         !< Ensemble index

! Local variables
type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug

! Interface
call oobump_registry%get(c_key_oobump, self)
call unstructured_grid_registry%get(c_key_ug, ug)

! Call Fortran
call add_oobump_member(self, ug, ie, iens)

end subroutine add_oobump_member_c
! ------------------------------------------------------------------------------
!> Remove ensemble member
subroutine remove_oobump_member_c(c_key_oobump, c_key_ug, ie, iens) bind(c, name='remove_oobump_member_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump !< OOBUMP
integer(c_int), intent(in) :: c_key_ug     !< Unstructured grid
integer(c_int), intent(in) :: ie           !< Ensemble member index
integer(c_int), intent(in) :: iens         !< Ensemble index

! Local variables
type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug

! Interface
call oobump_registry%get(c_key_oobump, self)
call unstructured_grid_registry%get(c_key_ug, ug)

! Call Fortran
call remove_oobump_member(self, ug, ie, iens)

end subroutine remove_oobump_member_c
! ------------------------------------------------------------------------------
!> Run BUMP drivers
subroutine run_oobump_drivers_c(c_key_oobump) bind(c, name='run_oobump_drivers_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump !< OOBUMP

! Local variables
type(oobump_type), pointer :: self

! Interface
call oobump_registry%get(c_key_oobump, self)

! Call Fortran
call run_oobump_drivers(self)

end subroutine run_oobump_drivers_c
! ------------------------------------------------------------------------------
!> Multiplication by BUMP vertical balance operator
subroutine multiply_oobump_vbal_c(c_key_oobump, c_key_ug) bind(c, name='multiply_oobump_vbal_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump !< OOBUMP
integer(c_int), intent(in) :: c_key_ug     !< Unstructured grid

! Local variables
type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug

! Interface
call oobump_registry%get(c_key_oobump, self)
call unstructured_grid_registry%get(c_key_ug, ug)

! Call Fortran
call multiply_oobump_vbal(self, ug)

end subroutine multiply_oobump_vbal_c
! ------------------------------------------------------------------------------
!> Multiplication by BUMP vertical balance operator inverse
subroutine multiply_oobump_vbal_inv_c(c_key_oobump, c_key_ug) bind(c, name='multiply_oobump_vbal_inv_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump !< OOBUMP
integer(c_int), intent(in) :: c_key_ug     !< Unstructured grid

! Local variables
type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug

! Interface
call oobump_registry%get(c_key_oobump, self)
call unstructured_grid_registry%get(c_key_ug, ug)

! Call Fortran
call multiply_oobump_vbal_inv(self, ug)

end subroutine multiply_oobump_vbal_inv_c
! ------------------------------------------------------------------------------
!> Multiplication by BUMP vertical balance operator adjoint
subroutine multiply_oobump_vbal_ad_c(c_key_oobump, c_key_ug) bind(c, name='multiply_oobump_vbal_ad_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump !< OOBUMP
integer(c_int), intent(in) :: c_key_ug     !< Unstructured grid

! Local variables
type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug

! Interface
call oobump_registry%get(c_key_oobump, self)
call unstructured_grid_registry%get(c_key_ug, ug)

! Call Fortran
call multiply_oobump_vbal_ad(self, ug)

end subroutine multiply_oobump_vbal_ad_c
! ------------------------------------------------------------------------------
!> Multiplication by BUMP vertical balance operator adjoint inverse
subroutine multiply_oobump_vbal_inv_ad_c(c_key_oobump, c_key_ug) bind(c, name='multiply_oobump_vbal_inv_ad_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump !< OOBUMP
integer(c_int), intent(in) :: c_key_ug     !< Unstructured grid

! Local variables
type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug

! Interface
call oobump_registry%get(c_key_oobump, self)
call unstructured_grid_registry%get(c_key_ug, ug)

! Call Fortran
call multiply_oobump_vbal_inv_ad(self, ug)

end subroutine multiply_oobump_vbal_inv_ad_c
! ------------------------------------------------------------------------------
!> Multiplication by BUMP NICAS operator
subroutine multiply_oobump_nicas_c(c_key_oobump, c_key_ug) bind(c, name='multiply_oobump_nicas_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump !< OOBUMP
integer(c_int), intent(in) :: c_key_ug     !< Unstructured grid

! Local variables
type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug

! Interface
call oobump_registry%get(c_key_oobump, self)
call unstructured_grid_registry%get(c_key_ug, ug)

! Call Fortran
call multiply_oobump_nicas(self, ug)

end subroutine multiply_oobump_nicas_c
! ------------------------------------------------------------------------------
!> Multiplication by BUMP NICAS operator square-root
subroutine multiply_oobump_nicas_sqrt_c(c_key_oobump, cv, c_key_ug) bind(c, name='multiply_oobump_nicas_sqrt_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump !< OOBUMP
real(c_double), intent(in) :: cv(:)        !< Control variable
integer(c_int), intent(in) :: c_key_ug     !< Unstructured grid

! Local variables
type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug

! Interface
call oobump_registry%get(c_key_oobump, self)
call unstructured_grid_registry%get(c_key_ug, ug)

! Call Fortran
call multiply_oobump_nicas_sqrt(self, cv, ug)

end subroutine multiply_oobump_nicas_sqrt_c
! ------------------------------------------------------------------------------
!> Multiplication by BUMP NICAS operator square-root adjoint
subroutine multiply_oobump_nicas_sqrt_ad_c(c_key_oobump, c_key_ug, cv) bind(c, name='multiply_oobump_nicas_sqrt_ad_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump !< OOBUMP
integer(c_int), intent(in) :: c_key_ug     !< Unstructured grid
real(c_double), intent(inout) :: cv(:)     !< Control variable

! Local variables
type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug

! Interface
call oobump_registry%get(c_key_oobump, self)
call unstructured_grid_registry%get(c_key_ug, ug)

! Call Fortran
call multiply_oobump_nicas_sqrt_ad(self, ug, cv)

end subroutine multiply_oobump_nicas_sqrt_ad_c
! ------------------------------------------------------------------------------
!> Randomize the BUMP NICAS operator
subroutine randomize_oobump_nicas_c(c_key_oobump, c_key_ug) bind(c, name='randomize_oobump_nicas_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump !< OOBUMP
integer(c_int), intent(in) :: c_key_ug     !< Unstructured grid

! Local variables
type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug

! Interface
call oobump_registry%get(c_key_oobump, self)
call unstructured_grid_registry%get(c_key_ug, ug)

! Call Fortran
call randomize_oobump_nicas(self, ug)

end subroutine randomize_oobump_nicas_c
! ------------------------------------------------------------------------------
!> Get BUMP parameter
subroutine get_oobump_param_c(c_key_oobump, nstr, cstr, c_key_ug) bind(c, name='get_oobump_param_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump       !< OOBUMP
integer(c_int), intent(in) :: nstr               !< Parameter name size
character(kind=c_char), intent(in) :: cstr(nstr) !< Parameter name
integer(c_int), intent(in) :: c_key_ug           !< Unstructured grid

! Local variables
type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug
integer :: istr
character(len=nstr) :: param

! Interface
call oobump_registry%get(c_key_oobump, self)
call unstructured_grid_registry%get(c_key_ug, ug)
param = ''
do istr=1,nstr
   param = trim(param)//cstr(istr)
end do

! Call Fortran
call get_oobump_param(self, param, ug)

end subroutine get_oobump_param_c
! ------------------------------------------------------------------------------
!> Set BUMP parameter
subroutine set_oobump_param_c(c_key_oobump, nstr, cstr, c_key_ug) bind(c, name='set_oobump_param_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump       !< OOBUMP
integer(c_int), intent(in) :: nstr               !< Parameter name size
character(kind=c_char), intent(in) :: cstr(nstr) !< Parameter name
integer(c_int), intent(in) :: c_key_ug           !< Unstructured grid

! Local variables
type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug
integer :: istr
character(len=nstr) :: param

! Interface
call oobump_registry%get(c_key_oobump, self)
call unstructured_grid_registry%get(c_key_ug, ug)
param = ''
do istr=1,nstr
   param = trim(param)//cstr(istr)
end do

! Call Fortran
call set_oobump_param(self, param, ug)

end subroutine set_oobump_param_c
! ------------------------------------------------------------------------------
end module oobump_interface
