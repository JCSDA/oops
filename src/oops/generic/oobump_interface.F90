! (C) Copyright 2017 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module oobump_interface

use atlas_module
use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm
use iso_c_binding
use kinds
use missing_values_mod
use oobump_mod
use type_bump, only: bump_type
use type_nam, only: nvmax,nlmax,nc3max,nscalesmax,ndirmax,nldwvmax

implicit none

private
! ------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
!> Create OOBUMP
subroutine oobump_create_c(c_key_oobump, c_comm, c_afunctionspace, c_afieldset, c_conf, c_grid, &
 & ens1_ne, ens1_nsub, ens2_ne, ens2_nsub) bind(c, name='oobump_create_f90')

implicit none

! Passed variables
integer(c_int), intent(inout) :: c_key_oobump    !< OOBUMP
type(c_ptr), value, intent(in) :: c_comm         !< Communicator
type(c_ptr),intent(in),value :: c_afunctionspace !< ATLAS function space pointer
type(c_ptr),intent(in),value :: c_afieldset      !< ATLAS fieldset pointer
type(c_ptr), intent(in) :: c_conf                !< Configuration
type(c_ptr), intent(in) :: c_grid                !< Grid configuration
integer(c_int), intent(in) :: ens1_ne            !< First ensemble size
integer(c_int), intent(in) :: ens1_nsub          !< Number of sub-ensembles in the first ensemble
integer(c_int), intent(in) :: ens2_ne            !< Second ensemble size
integer(c_int), intent(in) :: ens2_nsub          !< Number of sub-ensembles in the second ensemble

! Local variables
type(oobump_type), pointer :: self
type(fckit_mpi_comm) :: f_comm
type(atlas_functionspace) :: afunctionspace
type(atlas_fieldset) :: afieldset
type(fckit_configuration) :: f_conf
type(fckit_configuration) :: f_grid

! Interface
call oobump_registry%init()
call oobump_registry%add(c_key_oobump)
call oobump_registry%get(c_key_oobump, self)
f_comm = fckit_mpi_comm(c_comm)
afunctionspace = atlas_functionspace(c_afunctionspace)
afieldset = atlas_fieldset(c_afieldset)
f_conf = fckit_configuration(c_conf)
f_grid = fckit_configuration(c_grid)

! Call Fortran
call oobump_create(self, f_comm, afunctionspace, afieldset, f_conf, f_grid, &
 & ens1_ne, ens1_nsub, ens2_ne, ens2_nsub)

end subroutine oobump_create_c
! ------------------------------------------------------------------------------
!> Delete OOBUMP
subroutine oobump_delete_c(c_key_oobump) bind(c, name='oobump_delete_f90')

implicit none

! Passed variables
integer(c_int), intent(inout) :: c_key_oobump !< OOBUMP

! Local variables
type(oobump_type), pointer :: self

! Interface
call oobump_registry%get(c_key_oobump, self)

! Delete OOBUMP
call oobump_delete(self)

! Clean interface
call oobump_registry%remove(c_key_oobump)

end subroutine oobump_delete_c
! ------------------------------------------------------------------------------
!> Get control variable size
subroutine oobump_get_cv_size_c(c_key_oobump, n) bind(c, name='oobump_get_cv_size_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump !< OOBUMP
integer(c_int), intent(out) :: n           !< Control variable size

! Local variables
type(oobump_type), pointer :: self

! Interface
call oobump_registry%get(c_key_oobump, self)

! Call Fortran
call oobump_get_cv_size(self, n)

end subroutine oobump_get_cv_size_c
! ------------------------------------------------------------------------------
!> Add ensemble member
subroutine oobump_add_member_c(c_key_oobump, c_afieldset, ie, iens) bind(c, name='oobump_add_member_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump  !< OOBUMP
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer
integer(c_int), intent(in) :: ie            !< Ensemble member index
integer(c_int), intent(in) :: iens          !< Ensemble index

! Local variables
type(oobump_type), pointer :: self
type(atlas_fieldset) :: afieldset

! Interface
call oobump_registry%get(c_key_oobump, self)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call oobump_add_member(self, afieldset, ie, iens)

end subroutine oobump_add_member_c
! ------------------------------------------------------------------------------
!> Remove ensemble member
subroutine oobump_remove_member_c(c_key_oobump, c_afieldset, ie, iens) bind(c, name='oobump_remove_member_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump  !< OOBUMP
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer
integer(c_int), intent(in) :: ie            !< Ensemble member index
integer(c_int), intent(in) :: iens          !< Ensemble index

! Local variables
type(oobump_type), pointer :: self
type(atlas_fieldset) :: afieldset

! Interface
call oobump_registry%get(c_key_oobump, self)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call oobump_remove_member(self, afieldset, ie, iens)

end subroutine oobump_remove_member_c
! ------------------------------------------------------------------------------
!> Run BUMP drivers
subroutine oobump_run_drivers_c(c_key_oobump) bind(c, name='oobump_run_drivers_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump !< OOBUMP

! Local variables
type(oobump_type), pointer :: self

! Interface
call oobump_registry%get(c_key_oobump, self)

! Call Fortran
call oobump_run_drivers(self)

end subroutine oobump_run_drivers_c
! ------------------------------------------------------------------------------
!> Multiplication by BUMP vertical balance operator
subroutine oobump_multiply_vbal_c(c_key_oobump, c_afieldset) bind(c, name='oobump_multiply_vbal_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump  !< OOBUMP
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! Local variables
type(oobump_type), pointer :: self
type(atlas_fieldset) :: afieldset

! Interface
call oobump_registry%get(c_key_oobump, self)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call oobump_multiply_vbal(self, afieldset)

end subroutine oobump_multiply_vbal_c
! ------------------------------------------------------------------------------
!> Multiplication by BUMP vertical balance operator inverse
subroutine oobump_multiply_vbal_inv_c(c_key_oobump, c_afieldset) bind(c, name='oobump_multiply_vbal_inv_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump  !< OOBUMP
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! Local variables
type(oobump_type), pointer :: self
type(atlas_fieldset) :: afieldset

! Interface
call oobump_registry%get(c_key_oobump, self)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call oobump_multiply_vbal_inv(self, afieldset)

end subroutine oobump_multiply_vbal_inv_c
! ------------------------------------------------------------------------------
!> Multiplication by BUMP vertical balance operator adjoint
subroutine oobump_multiply_vbal_ad_c(c_key_oobump, c_afieldset) bind(c, name='oobump_multiply_vbal_ad_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump  !< OOBUMP
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! Local variables
type(oobump_type), pointer :: self
type(atlas_fieldset) :: afieldset

! Interface
call oobump_registry%get(c_key_oobump, self)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call oobump_multiply_vbal_ad(self, afieldset)

end subroutine oobump_multiply_vbal_ad_c
! ------------------------------------------------------------------------------
!> Multiplication by BUMP vertical balance operator adjoint inverse
subroutine oobump_multiply_vbal_inv_ad_c(c_key_oobump, c_afieldset) bind(c, name='oobump_multiply_vbal_inv_ad_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump  !< OOBUMP
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! Local variables
type(oobump_type), pointer :: self
type(atlas_fieldset) :: afieldset

! Interface
call oobump_registry%get(c_key_oobump, self)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call oobump_multiply_vbal_inv_ad(self, afieldset)

end subroutine oobump_multiply_vbal_inv_ad_c
! ------------------------------------------------------------------------------
!> Multiplication by BUMP NICAS operator
subroutine oobump_multiply_nicas_c(c_key_oobump, c_afieldset) bind(c, name='oobump_multiply_nicas_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump  !< OOBUMP
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! Local variables
type(oobump_type), pointer :: self
type(atlas_fieldset) :: afieldset

! Interface
call oobump_registry%get(c_key_oobump, self)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call oobump_multiply_nicas(self, afieldset)

end subroutine oobump_multiply_nicas_c
! ------------------------------------------------------------------------------
!> Multiplication by BUMP NICAS operator square-root
subroutine oobump_multiply_nicas_sqrt_c(c_key_oobump, cv, c_afieldset) bind(c, name='oobump_multiply_nicas_sqrt_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump  !< OOBUMP
real(c_double), intent(in) :: cv(:)         !< Control variable
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! Local variables
type(oobump_type), pointer :: self
type(atlas_fieldset) :: afieldset

! Interface
call oobump_registry%get(c_key_oobump, self)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call oobump_multiply_nicas_sqrt(self, cv, afieldset)

end subroutine oobump_multiply_nicas_sqrt_c
! ------------------------------------------------------------------------------
!> Multiplication by BUMP NICAS operator square-root adjoint
subroutine oobump_multiply_nicas_sqrt_ad_c(c_key_oobump, c_afieldset, cv) bind(c, name='oobump_multiply_nicas_sqrt_ad_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump  !< OOBUMP
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer
real(c_double), intent(inout) :: cv(:)      !< Control variable

! Local variables
type(oobump_type), pointer :: self
type(atlas_fieldset) :: afieldset

! Interface
call oobump_registry%get(c_key_oobump, self)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call oobump_multiply_nicas_sqrt_ad(self, afieldset, cv)

end subroutine oobump_multiply_nicas_sqrt_ad_c
! ------------------------------------------------------------------------------
!> Randomize the BUMP NICAS operator
subroutine oobump_randomize_nicas_c(c_key_oobump, c_afieldset) bind(c, name='oobump_randomize_nicas_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump  !< OOBUMP
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! Local variables
type(oobump_type), pointer :: self
type(atlas_fieldset) :: afieldset

! Interface
call oobump_registry%get(c_key_oobump, self)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call oobump_randomize_nicas(self, afieldset)

end subroutine oobump_randomize_nicas_c
! ------------------------------------------------------------------------------
!> Get BUMP parameter
subroutine oobump_get_param_c(c_key_oobump, nstr, cstr, c_afieldset) bind(c, name='oobump_get_param_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump       !< OOBUMP
integer(c_int), intent(in) :: nstr               !< Parameter name size
character(kind=c_char), intent(in) :: cstr(nstr) !< Parameter name
type(c_ptr),intent(in),value :: c_afieldset      !< ATLAS fieldset pointer

! Local variables
type(oobump_type), pointer :: self
integer :: istr
character(len=nstr) :: param
type(atlas_fieldset) :: afieldset

! Interface
call oobump_registry%get(c_key_oobump, self)
param = ''
do istr=1,nstr
  param = trim(param)//cstr(istr)
end do
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call oobump_get_param(self, param, afieldset)

end subroutine oobump_get_param_c
! ------------------------------------------------------------------------------
!> Set BUMP parameter
subroutine oobump_set_param_c(c_key_oobump, nstr, cstr, c_afieldset) bind(c, name='oobump_set_param_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump       !< OOBUMP
integer(c_int), intent(in) :: nstr               !< Parameter name size
character(kind=c_char), intent(in) :: cstr(nstr) !< Parameter name
type(c_ptr),intent(in),value :: c_afieldset      !< ATLAS fieldset pointer

! Local variables
type(oobump_type), pointer :: self
integer :: istr
character(len=nstr) :: param
type(atlas_fieldset) :: afieldset

! Interface
call oobump_registry%get(c_key_oobump, self)
param = ''
do istr=1,nstr
  param = trim(param)//cstr(istr)
end do
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call oobump_set_param(self, param, afieldset)

end subroutine oobump_set_param_c
! ------------------------------------------------------------------------------
end module oobump_interface
