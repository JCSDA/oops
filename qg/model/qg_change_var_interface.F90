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
!> Setup change of variable
subroutine qg_change_var_setup_c(c_key_self,c_vars_in,c_vars_out) bind (c,name='qg_change_var_setup_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self           !< Variable change
type(c_ptr),value,intent(in) :: c_vars_in  !< Input variables
type(c_ptr),value,intent(in) :: c_vars_out !< Output variables

! Local variable
type(qg_change_var_config),pointer :: self
type(oops_variables) :: vars_in,vars_out

! Interface
call qg_change_var_registry%init()
call qg_change_var_registry%add(c_key_self)
call qg_change_var_registry%get(c_key_self,self)
vars_in = oops_variables(c_vars_in)
vars_out = oops_variables(c_vars_out)

! Call Fortran
call qg_change_var_setup(self,vars_in,vars_out)

end subroutine qg_change_var_setup_c
! ------------------------------------------------------------------------------
!> Delete error covariance matrix
subroutine qg_change_var_delete_c(c_key_self) bind (c,name='qg_change_var_delete_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Error covariance configuration

! Clear interface
call qg_change_var_registry%remove(c_key_self)

end subroutine qg_change_var_delete_c
! ------------------------------------------------------------------------------
!> Change of variable
subroutine qg_change_var_c(c_key_conf,c_key_fld_in,c_key_fld_out) bind (c,name='qg_change_var_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_conf    !< Variable change
integer(c_int),intent(in) :: c_key_fld_in  !< Input field
integer(c_int),intent(in) :: c_key_fld_out !< Output field

! Local variables
type(qg_change_var_config),pointer :: conf
type(qg_fields),pointer :: fld_in,fld_out

! Interface
call qg_change_var_registry%get(c_key_conf,conf)
call qg_fields_registry%get(c_key_fld_in,fld_in)
call qg_fields_registry%get(c_key_fld_out,fld_out)

! Call Fortran
call qg_change_var(conf,fld_in,fld_out)

end subroutine qg_change_var_c
! ------------------------------------------------------------------------------
!> Change of variable - inverse
subroutine qg_change_var_inv_c(c_key_conf,c_key_fld_in,c_key_fld_out) bind (c,name='qg_change_var_inv_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_conf    !< Configuration
integer(c_int),intent(in) :: c_key_fld_in  !< Input field
integer(c_int),intent(in) :: c_key_fld_out !< Output field

! Local variables
type(qg_change_var_config),pointer :: conf
type(qg_fields),pointer :: fld_in,fld_out

! Interface
call qg_change_var_registry%get(c_key_conf,conf)
call qg_fields_registry%get(c_key_fld_in,fld_in)
call qg_fields_registry%get(c_key_fld_out,fld_out)

! Call Fortran
call qg_change_var_inv(conf,fld_in,fld_out)

end subroutine qg_change_var_inv_c
! ------------------------------------------------------------------------------
!> Change of variable - adjoint
subroutine qg_change_var_ad_c(c_key_conf,c_key_fld_in,c_key_fld_out) bind (c,name='qg_change_var_ad_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_conf    !< Configuration
integer(c_int),intent(in) :: c_key_fld_in  !< Input field
integer(c_int),intent(in) :: c_key_fld_out !< Output field

! Local variables
type(qg_change_var_config),pointer :: conf
type(qg_fields),pointer :: fld_in,fld_out

! Interface
call qg_change_var_registry%get(c_key_conf,conf)
call qg_fields_registry%get(c_key_fld_in,fld_in)
call qg_fields_registry%get(c_key_fld_out,fld_out)

! Call Fortran
call qg_change_var_ad(conf,fld_in,fld_out)

end subroutine qg_change_var_ad_c
! ------------------------------------------------------------------------------
!> Change of variable - inverse adjoint
subroutine qg_change_var_inv_ad_c(c_key_conf,c_key_fld_in,c_key_fld_out) bind (c,name='qg_change_var_inv_ad_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_conf    !< Configuration
integer(c_int),intent(in) :: c_key_fld_in  !< Input field
integer(c_int),intent(in) :: c_key_fld_out !< Output field

! Local variables
type(qg_change_var_config),pointer :: conf
type(qg_fields),pointer :: fld_in,fld_out

! Interface
call qg_change_var_registry%get(c_key_conf,conf)
call qg_fields_registry%get(c_key_fld_in,fld_in)
call qg_fields_registry%get(c_key_fld_out,fld_out)

! Call Fortran
call qg_change_var_inv_ad(conf,fld_in,fld_out)

end subroutine qg_change_var_inv_ad_c
! ------------------------------------------------------------------------------
end module qg_change_var_interface
