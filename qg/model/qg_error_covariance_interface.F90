! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_error_covariance_interface

use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log
use iso_c_binding
use kinds
use qg_error_covariance_mod
use qg_fields_mod
use qg_geom_mod
use oops_variables_mod

implicit none

public
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Setup error covariance matrix
subroutine qg_error_covariance_setup_c(c_key_self,c_conf,c_key_geom) bind (c,name='qg_error_covariance_setup_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Error covariance configuration
type(c_ptr),value,intent(in) :: c_conf     !< Configuration
integer(c_int),intent(in) :: c_key_geom    !< Geometry

! Local variables
type(fckit_configuration) :: f_conf
type(qg_error_covariance_config),pointer :: self
type(qg_geom),pointer :: geom

! Interface
f_conf = fckit_configuration(c_conf)
call qg_geom_registry%get(c_key_geom,geom)
call qg_error_covariance_registry%init()
call qg_error_covariance_registry%add(c_key_self)
call qg_error_covariance_registry%get(c_key_self,self)

! Call Fortran
call qg_error_covariance_setup(self,f_conf,geom)

end subroutine qg_error_covariance_setup_c
! ------------------------------------------------------------------------------
!> Delete error covariance matrix
subroutine qg_error_covariance_delete_c(c_key_self) bind (c,name='qg_error_covariance_delete_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Error covariance configuration

! Local variables
type(qg_error_covariance_config),pointer :: self

! Interface
call qg_error_covariance_registry%get(c_key_self,self)

! Call Fortran
call qg_error_covariance_delete(self)

! Clear interface
call qg_error_covariance_registry%remove(c_key_self)

end subroutine qg_error_covariance_delete_c
! ------------------------------------------------------------------------------
!> Multiply by error covariance matrix
subroutine qg_error_covariance_mult_c(c_key_self,c_key_in,c_key_out) bind(c,name='qg_error_covariance_mult_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< Error covariance configuration
integer(c_int),intent(in) :: c_key_in    !< Input field
integer(c_int),intent(in) :: c_key_out   !< Output field

! Local variables
type(qg_error_covariance_config),pointer :: self
type(qg_fields),pointer :: fld_in,fld_out

! Interface
call qg_error_covariance_registry%get(c_key_self,self)
call qg_fields_registry%get(c_key_in,fld_in)
call qg_fields_registry%get(c_key_out,fld_out)

! Call Fortran
call qg_error_covariance_mult(self,fld_in,fld_out)

end subroutine qg_error_covariance_mult_c
! ------------------------------------------------------------------------------
!> Randomize error covariance
subroutine qg_error_covariance_randomize_c(c_key_self,c_key_out) bind(c,name='qg_error_covariance_randomize_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Error covariance configuration
integer(c_int),intent(in) :: c_key_out  !< Output field

! Local variables
type(qg_error_covariance_config),pointer :: self
type(qg_fields),pointer :: fld_out

! Interface
call qg_error_covariance_registry%get(c_key_self,self)
call qg_fields_registry%get(c_key_out,fld_out)

! Call Fortran
call qg_error_covariance_randomize(self,fld_out)

end subroutine qg_error_covariance_randomize_c
! ------------------------------------------------------------------------------
end module qg_error_covariance_interface
