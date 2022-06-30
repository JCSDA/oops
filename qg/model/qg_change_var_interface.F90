! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_change_var_interface

use iso_c_binding
use qg_change_var_mod
use qg_fields_mod
use qg_geom_mod
use oops_variables_mod

implicit none

private
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Change of variable
subroutine qg_change_var_c(c_key_fld,c_vars) bind (c,name='qg_change_var_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld  !< Input field
type(c_ptr),value,intent(in) :: c_vars  !< Output variables

! Local variables
type(qg_fields),pointer :: fld
type(oops_variables) :: vars

! Interface
call qg_fields_registry%get(c_key_fld,fld)
vars = oops_variables(c_vars)

! Call Fortran
call qg_change_var(fld,vars)

end subroutine qg_change_var_c
! ------------------------------------------------------------------------------
!> Change of variable - tangent linear
subroutine qg_change_var_tl_c(c_key_fld,c_vars) bind (c,name='qg_change_var_tl_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld  !< Input field
type(c_ptr),value,intent(in) :: c_vars  !< Output variables

! Local variables
type(qg_fields),pointer :: fld
type(oops_variables) :: vars

! Interface
call qg_fields_registry%get(c_key_fld,fld)
vars = oops_variables(c_vars)

! Call Fortran
call qg_change_var_tl(fld,vars)

end subroutine qg_change_var_tl_c
! ------------------------------------------------------------------------------
!> Change of variable - adjoint
subroutine qg_change_var_ad_c(c_key_fld,c_vars) bind (c,name='qg_change_var_ad_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld  !< Input field
type(c_ptr),value,intent(in) :: c_vars  !< Output variables

! Local variables
type(qg_fields),pointer :: fld
type(oops_variables) :: vars

! Interface
call qg_fields_registry%get(c_key_fld,fld)
vars = oops_variables(c_vars)

! Call Fortran
call qg_change_var_ad(fld,vars)

end subroutine qg_change_var_ad_c
! ------------------------------------------------------------------------------
end module qg_change_var_interface
