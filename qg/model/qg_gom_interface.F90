! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_gom_interface

use atlas_module, only: atlas_field
use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use datetime_mod
use qg_geom_mod
use qg_gom_mod
use qg_locs_mod
use oops_variables_mod

implicit none

private
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Setup GOM
subroutine qg_gom_setup_c(c_key_self,c_npaths,c_vars,nlevs) bind(c,name='qg_gom_setup_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< GOM
integer(c_int),intent(in)    :: c_npaths   !< Number of paths along which all variables
                                           !< will be interpolated
type(c_ptr),value,intent(in) :: c_vars     !< Variables
integer(c_int),intent(in)    :: nlevs

! Local variables
type(qg_gom),pointer :: self
type(oops_variables) :: vars

! Interface
call qg_gom_registry%init()
call qg_gom_registry%add(c_key_self)
call qg_gom_registry%get(c_key_self,self)
vars = oops_variables(c_vars)

! Call Fortran
call qg_gom_setup(self,c_npaths,vars,nlevs)

end subroutine qg_gom_setup_c
! ------------------------------------------------------------------------------
!> Create GOM and do nothing
subroutine qg_gom_create_c(c_key_self) bind(c,name='qg_gom_create_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< GOM

! Local variables
type(qg_gom),pointer :: self

! Interface
call qg_gom_registry%init()
call qg_gom_registry%add(c_key_self)
call qg_gom_registry%get(c_key_self,self)

end subroutine qg_gom_create_c
! ------------------------------------------------------------------------------
!> Delete GOM
subroutine qg_gom_delete_c(c_key_self) bind(c,name='qg_gom_delete_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< GOM

! Local variables
type(qg_gom),pointer :: self

! Interface
call qg_gom_registry%get(c_key_self,self)

! Call Fortran
call qg_gom_delete(self)

! Clear interface
call qg_gom_registry%remove(c_key_self)

end subroutine qg_gom_delete_c
! ------------------------------------------------------------------------------
!> Copy GOM
subroutine qg_gom_copy_c(c_key_self,c_key_other) bind(c,name='qg_gom_copy_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< GOM
integer(c_int),intent(in) :: c_key_other   !< Other GOM

! Local variables
type(qg_gom),pointer :: self
type(qg_gom),pointer :: other

! Interface
call qg_gom_registry%get(c_key_self,self)
call qg_gom_registry%get(c_key_other,other)

! Call Fortran
call qg_gom_copy(self,other)

end subroutine qg_gom_copy_c
! ------------------------------------------------------------------------------
subroutine qg_gom_fill_c(c_key, lvar, c_var, c_nloc, c_indx, c_nlev, c_vals) &
  bind(c, name="qg_gom_fill_f90")
implicit none
integer(c_int), intent(in) :: c_key
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(in) :: c_nloc
integer(c_int), intent(in) :: c_indx(c_nloc)
integer(c_int), intent(in) :: c_nlev
real(c_double), intent(in) :: c_vals(c_nlev, c_nloc)

type(qg_gom), pointer :: self

call qg_gom_registry%get(c_key,self)

call qg_gom_fill(self, lvar, c_var, c_nloc, c_indx, c_nlev, c_vals)

end subroutine qg_gom_fill_c

! ------------------------------------------------------------------------------
subroutine qg_gom_fillad_c(c_key, lvar, c_var, c_nloc, c_indx, c_nlev, c_vals) &
  bind(c, name="qg_gom_fillad_f90")
implicit none
integer(c_int), intent(in) :: c_key
integer(c_int), intent(in) :: lvar
character(kind=c_char, len=1), intent(in) :: c_var(lvar+1)
integer(c_int), intent(in) :: c_nloc
integer(c_int), intent(in) :: c_indx(c_nloc)
integer(c_int), intent(in) :: c_nlev
real(c_double), intent(inout) :: c_vals(c_nlev, c_nloc)

type(qg_gom),pointer :: self

call qg_gom_registry%get(c_key, self)

call qg_gom_fillad(self, lvar, c_var, c_nloc, c_indx, c_nlev, c_vals)

end subroutine qg_gom_fillad_c
! ------------------------------------------------------------------------------
!> Set GOM to zero
subroutine qg_gom_zero_c(c_key_self) bind(c,name='qg_gom_zero_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< GOM

! Local variables
type(qg_gom),pointer :: self

! Interface
call qg_gom_registry%get(c_key_self,self)

! Call Fortran
call qg_gom_zero(self)

end subroutine qg_gom_zero_c
! ------------------------------------------------------------------------------
!> Get GOM absolute value
subroutine qg_gom_abs_c(c_key_self) bind(c,name='qg_gom_abs_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< GOM

! Local variables
type(qg_gom),pointer :: self

! Interface
call qg_gom_registry%get(c_key_self,self)

! Call Fortran
call qg_gom_abs(self)

end subroutine qg_gom_abs_c
! ------------------------------------------------------------------------------
!> Generate random GOM
subroutine qg_gom_random_c(c_key_self) bind(c,name='qg_gom_random_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< GOM

! Local variables
type(qg_gom),pointer :: self

! Interface
call qg_gom_registry%get(c_key_self,self)

! Call Fortran
call qg_gom_random(self)

end subroutine qg_gom_random_c
! ------------------------------------------------------------------------------
!> Multiply GOM with a scalar
subroutine qg_gom_mult_c(c_key_self,c_zz) bind(c,name='qg_gom_mult_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< GOM
real(c_double),intent(in) :: c_zz       !< Multiplier

! Local variables
type(qg_gom),pointer :: self
integer :: jo,jv

! Interface
call qg_gom_registry%get(c_key_self,self)

! Call Fortran
call qg_gom_mult(self,c_zz)

end subroutine qg_gom_mult_c
! ------------------------------------------------------------------------------
!> Add GOM
subroutine qg_gom_add_c(c_key_self,c_key_other) bind(c,name='qg_gom_add_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< GOM
integer(c_int),intent(in) :: c_key_other !< Other GOM

! Local variables
type(qg_gom),pointer :: self
type(qg_gom),pointer :: other

! Interface
call qg_gom_registry%get(c_key_self,self)
call qg_gom_registry%get(c_key_other,other)

! Call Fortran
call qg_gom_add(self,other)

end subroutine qg_gom_add_c
! ------------------------------------------------------------------------------
!> Subtract GOM
subroutine qg_gom_diff_c(c_key_self,c_key_other) bind(c,name='qg_gom_diff_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< GOM
integer(c_int),intent(in) :: c_key_other !< Other GOM

! Local variables
type(qg_gom),pointer :: self
type(qg_gom),pointer :: other

! Interface
call qg_gom_registry%get(c_key_self,self)
call qg_gom_registry%get(c_key_other,other)

! Call Fortran
call qg_gom_diff(self,other)

end subroutine qg_gom_diff_c
! ------------------------------------------------------------------------------
!> Schur product for GOM
subroutine qg_gom_schurmult_c(c_key_self,c_key_other) bind(c,name='qg_gom_schurmult_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< GOM
integer(c_int),intent(in) :: c_key_other !< Other GOM

! Local variables
type(qg_gom),pointer :: self
type(qg_gom),pointer :: other

! Interface
call qg_gom_registry%get(c_key_self,self)
call qg_gom_registry%get(c_key_other,other)

! Call Fortran
call qg_gom_schurmult(self,other)

end subroutine qg_gom_schurmult_c
! ------------------------------------------------------------------------------
!> Schur division for GOM
subroutine qg_gom_divide_c(c_key_self,c_key_other) bind(c,name='qg_gom_divide_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< GOM
integer(c_int),intent(in) :: c_key_other !< Other GOM

! Local variables
type(qg_gom),pointer :: self
type(qg_gom),pointer :: other

! Interface
call qg_gom_registry%get(c_key_self,self)
call qg_gom_registry%get(c_key_other,other)

! Call Fortran
call qg_gom_divide(self,other)

end subroutine qg_gom_divide_c
! ------------------------------------------------------------------------------
!> Compute GOM RMS
subroutine qg_gom_rms_c(c_key_self,c_rms) bind(c,name='qg_gom_rms_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< GOM
real(c_double),intent(inout) :: c_rms   !< RMS

! Local variables
type(qg_gom),pointer :: self

! Interface
call qg_gom_registry%get(c_key_self,self)

! Call Fortran
call qg_gom_rms(self,c_rms)

end subroutine qg_gom_rms_c
! ------------------------------------------------------------------------------
!> GOM dot product
subroutine qg_gom_dotprod_c(c_key_gom1,c_key_gom2,c_prod) bind(c,name='qg_gom_dotprod_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_gom1 !< GOM 1
integer(c_int),intent(in) :: c_key_gom2 !< GOM 2
real(c_double),intent(inout) :: c_prod  !< Dot product

! Local variables
type(qg_gom),pointer :: gom1,gom2

! Interface
call qg_gom_registry%get(c_key_gom1,gom1)
call qg_gom_registry%get(c_key_gom2,gom2)

! Call Fortran
call qg_gom_dotprod(gom1,gom2,c_prod)

end subroutine qg_gom_dotprod_c
! ------------------------------------------------------------------------------
!> Compute GOM statistics
subroutine qg_gom_stats_c(c_key_self,c_kobs,c_pmin,c_pmax,c_prms) bind(c,name='qg_gom_stats_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< GOM
integer(c_int),intent(inout) :: c_kobs  !< Number of observations
real(c_double),intent(inout) :: c_pmin  !< Minimum value
real(c_double),intent(inout) :: c_pmax  !< Maximum value
real(c_double),intent(inout) :: c_prms  !< RMS

! Local variables
type(qg_gom),pointer :: self

! Interface
call qg_gom_registry%get(c_key_self,self)

! Call Fortran
call qg_gom_stats(self,c_kobs,c_pmin,c_pmax,c_prms)

end subroutine qg_gom_stats_c
! ------------------------------------------------------------------------------
!> Find and locate GOM max. value
subroutine qg_gom_maxloc_c(c_key_self,c_mxval,c_mxloc,c_mxvar) bind(c,name='qg_gom_maxloc_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< GOM
real(c_double),intent(inout) :: c_mxval !< Maximum value
integer(c_int),intent(inout) :: c_mxloc !< Location of maximum value
type(c_ptr),value,intent(in) :: c_mxvar !< Variable of maximum value

! Local variables
type(qg_gom),pointer :: self
type(oops_variables) :: mxvar

! Interface
call qg_gom_registry%get(c_key_self,self)
mxvar = oops_variables(c_mxvar)

! Call Fortran
call qg_gom_maxloc(self,c_mxval,c_mxloc,mxvar)

end subroutine qg_gom_maxloc_c
! ------------------------------------------------------------------------------
!> Read GOM from file
subroutine qg_gom_read_file_c(c_key_self,c_vars,c_conf) bind(c,name='qg_gom_read_file_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< GOM
type(c_ptr),value,intent(in) :: c_vars     !< Variables
type(c_ptr),value,intent(in) :: c_conf     !< Configuration

! Local variables
type(fckit_configuration) :: f_conf
type(qg_gom),pointer :: self
type(oops_variables) :: vars

! Interface
call qg_gom_registry%get(c_key_self,self)
f_conf = fckit_configuration(c_conf)
vars = oops_variables(c_vars)

! Call Fortran
call qg_gom_read_file(self,vars,f_conf)

end subroutine qg_gom_read_file_c
! ------------------------------------------------------------------------------
!> Write GOM to file
subroutine qg_gom_write_file_c(c_key_self,c_conf) bind(c,name='qg_gom_write_file_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< GOM
type(c_ptr),value,intent(in) :: c_conf  !< Configuration

! Local variables
type(fckit_configuration) :: f_conf
type(qg_gom),pointer :: self

! Interface
f_conf = fckit_configuration(c_conf)
call qg_gom_registry%get(c_key_self,self)

! Call Fortran
call qg_gom_write_file(self,f_conf)

end subroutine qg_gom_write_file_c
! ------------------------------------------------------------------------------
!> GOM analytic initialization
subroutine qg_gom_analytic_init_c(c_key_self,c_locs,c_conf) bind(c,name='qg_gom_analytic_init_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< GOM
type(c_ptr),value,intent(in) :: c_locs  !< Locations
type(c_ptr),value,intent(in) :: c_conf  !< Configuration

! Local variables
type(fckit_configuration) :: f_conf
type(qg_gom),pointer :: self
type(qg_locs) :: locs

! Interface
f_conf = fckit_configuration(c_conf)
call qg_gom_registry%get(c_key_self,self)
locs = qg_locs(c_locs)

! Call Fortran
call qg_gom_analytic_init(self,locs,f_conf)

end subroutine qg_gom_analytic_init_c
! ------------------------------------------------------------------------------
end module qg_gom_interface
