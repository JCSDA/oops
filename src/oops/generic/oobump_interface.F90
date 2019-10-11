! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module oobump_interface

use fckit_configuration_module, only: fckit_configuration
use fckit_mpi_module, only: fckit_mpi_comm
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
subroutine oobump_create_c(c_key_oobump, c_key_ug, c_conf, ens1_ne, ens1_nsub, ens2_ne, ens2_nsub) bind(c, name='oobump_create_f90')

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
type(fckit_mpi_comm) :: f_comm

! Interface
f_conf = fckit_configuration(c_conf)
call oobump_registry%init()
call oobump_registry%add(c_key_oobump)
call oobump_registry%get(c_key_oobump, self)
call unstructured_grid_registry%get(c_key_ug, ug)

! Get MPI communicator (use default comm now, could use a specific name in the future) 
f_comm = fckit_mpi_comm()

! Call Fortran
call oobump_create(self, ug, f_conf, ens1_ne, ens1_nsub, ens2_ne, ens2_nsub, f_comm)

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
!> Get colocated flag
subroutine oobump_get_colocated_c(c_key_oobump, colocated) bind(c, name='oobump_get_colocated_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump !< OOBUMP
integer(c_int), intent(out) :: colocated   !< Colocated flag

! Local variables
type(oobump_type), pointer :: self

! Interface
call oobump_registry%get(c_key_oobump, self)

! Call Fortran
call oobump_get_colocated(self, colocated)

end subroutine oobump_get_colocated_c
! ------------------------------------------------------------------------------
!> Get number of timeslots
subroutine oobump_get_nts_c(c_key_oobump, nts) bind(c, name='oobump_get_nts_f90')

implicit none

! Passed variables
integer(c_int), intent(in) :: c_key_oobump !< OOBUMP
integer(c_int), intent(out) :: nts         !< Number of timeslots

! Local variables
type(oobump_type), pointer :: self

! Interface
call oobump_registry%get(c_key_oobump, self)

! Call Fortran
call oobump_get_nts(self, nts)

end subroutine oobump_get_nts_c
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
subroutine oobump_add_member_c(c_key_oobump, c_key_ug, ie, iens) bind(c, name='oobump_add_member_f90')

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
call oobump_add_member(self, ug, ie, iens)

end subroutine oobump_add_member_c
! ------------------------------------------------------------------------------
!> Remove ensemble member
subroutine oobump_remove_member_c(c_key_oobump, c_key_ug, ie, iens) bind(c, name='oobump_remove_member_f90')

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
call oobump_remove_member(self, ug, ie, iens)

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
subroutine oobump_multiply_vbal_c(c_key_oobump, c_key_ug) bind(c, name='oobump_multiply_vbal_f90')

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
call oobump_multiply_vbal(self, ug)

end subroutine oobump_multiply_vbal_c
! ------------------------------------------------------------------------------
!> Multiplication by BUMP vertical balance operator inverse
subroutine oobump_multiply_vbal_inv_c(c_key_oobump, c_key_ug) bind(c, name='oobump_multiply_vbal_inv_f90')

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
call oobump_multiply_vbal_inv(self, ug)

end subroutine oobump_multiply_vbal_inv_c
! ------------------------------------------------------------------------------
!> Multiplication by BUMP vertical balance operator adjoint
subroutine oobump_multiply_vbal_ad_c(c_key_oobump, c_key_ug) bind(c, name='oobump_multiply_vbal_ad_f90')

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
call oobump_multiply_vbal_ad(self, ug)

end subroutine oobump_multiply_vbal_ad_c
! ------------------------------------------------------------------------------
!> Multiplication by BUMP vertical balance operator adjoint inverse
subroutine oobump_multiply_vbal_inv_ad_c(c_key_oobump, c_key_ug) bind(c, name='oobump_multiply_vbal_inv_ad_f90')

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
call oobump_multiply_vbal_inv_ad(self, ug)

end subroutine oobump_multiply_vbal_inv_ad_c
! ------------------------------------------------------------------------------
!> Multiplication by BUMP NICAS operator
subroutine oobump_multiply_nicas_c(c_key_oobump, c_key_ug) bind(c, name='oobump_multiply_nicas_f90')

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
call oobump_multiply_nicas(self, ug)

end subroutine oobump_multiply_nicas_c
! ------------------------------------------------------------------------------
!> Multiplication by BUMP NICAS operator square-root
subroutine oobump_multiply_nicas_sqrt_c(c_key_oobump, cv, c_key_ug) bind(c, name='oobump_multiply_nicas_sqrt_f90')

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
call oobump_multiply_nicas_sqrt(self, cv, ug)

end subroutine oobump_multiply_nicas_sqrt_c
! ------------------------------------------------------------------------------
!> Multiplication by BUMP NICAS operator square-root adjoint
subroutine oobump_multiply_nicas_sqrt_ad_c(c_key_oobump, c_key_ug, cv) bind(c, name='oobump_multiply_nicas_sqrt_ad_f90')

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
call oobump_multiply_nicas_sqrt_ad(self, ug, cv)

end subroutine oobump_multiply_nicas_sqrt_ad_c
! ------------------------------------------------------------------------------
!> Randomize the BUMP NICAS operator
subroutine oobump_randomize_nicas_c(c_key_oobump, c_key_ug) bind(c, name='oobump_randomize_nicas_f90')

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
call oobump_randomize_nicas(self, ug)

end subroutine oobump_randomize_nicas_c
! ------------------------------------------------------------------------------
!> Get BUMP parameter
subroutine oobump_get_param_c(c_key_oobump, nstr, cstr, c_key_ug) bind(c, name='oobump_get_param_f90')

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
call oobump_get_param(self, param, ug)

end subroutine oobump_get_param_c
! ------------------------------------------------------------------------------
!> Set BUMP parameter
subroutine oobump_set_param_c(c_key_oobump, nstr, cstr, c_key_ug) bind(c, name='oobump_set_param_f90')

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
call oobump_set_param(self, param, ug)

end subroutine oobump_set_param_c
! ------------------------------------------------------------------------------
end module oobump_interface
