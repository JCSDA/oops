! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_error_stddev_interface

use iso_c_binding
use kinds
use qg_error_stddev_mod
use qg_fields_mod

implicit none

private
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Setup error standard deviation matrix
subroutine qg_error_stddev_setup_c(c_key_self,c_conf) bind (c,name='qg_error_stddev_setup_f90')

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Error standard deviation structure
type(c_ptr),intent(in) :: c_conf           !< Configuration

! Local variables
type(qg_error_stddev_config),pointer :: self

! Interface
call qg_error_stddev_registry%init()
call qg_error_stddev_registry%add(c_key_self)
call qg_error_stddev_registry%get(c_key_self,self)

! Call Fortran
call qg_error_stddev_setup(c_conf,self)

end subroutine qg_error_stddev_setup_c
! ------------------------------------------------------------------------------
!> Delete error standard deviation matrix
subroutine qg_error_stddev_delete_c(c_key_self) bind (c,name='qg_error_stddev_delete_f90')

! Passed variable
integer(c_int),intent(inout) :: c_key_self !< Error standard deviation

! Clear interface
call qg_error_stddev_registry%remove(c_key_self)

end subroutine qg_error_stddev_delete_c
! ------------------------------------------------------------------------------
!> Multiply by error standard deviation
subroutine qg_error_stddev_mult_c(c_key_conf,c_key_in,c_key_out) bind(c,name='qg_error_stddev_mult_f90')

! Passed variables
integer(c_int),intent(in) :: c_key_conf !< Error standard deviation
integer(c_int),intent(in) :: c_key_in   !< Input field
integer(c_int),intent(in) :: c_key_out  !< Output field

! Local variables
type(qg_error_stddev_config),pointer :: conf
type(qg_fields),pointer :: fld_in,fld_out

! Interface
call qg_error_stddev_registry%get(c_key_conf,conf)
call qg_fields_registry%get(c_key_in,fld_in)
call qg_fields_registry%get(c_key_out,fld_out)

! Call Fortran
call qg_error_stddev_mult(conf,fld_in,fld_out)

end subroutine qg_error_stddev_mult_c
! ------------------------------------------------------------------------------
!> Multiply by inverse of error standard deviation
subroutine qg_error_stddev_inv_mult_c(c_key_conf,c_key_in,c_key_out) bind(c,name='qg_error_stddev_inv_mult_f90')

! Passed variables
integer(c_int),intent(in) :: c_key_conf !< Error standard deviation
integer(c_int),intent(in) :: c_key_in   !< Input field
integer(c_int),intent(in) :: c_key_out  !< Output field

! Local variables
type(qg_error_stddev_config),pointer :: conf
type(qg_fields),pointer :: fld_in,fld_out

! Interface
call qg_error_stddev_registry%get(c_key_conf,conf)
call qg_fields_registry%get(c_key_in,fld_in)
call qg_fields_registry%get(c_key_out,fld_out)

! Call Fortran
call qg_error_stddev_inv_mult(conf,fld_in,fld_out)

end subroutine qg_error_stddev_inv_mult_c
! ------------------------------------------------------------------------------
end module qg_error_stddev_interface
