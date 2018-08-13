! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! ------------------------------------------------------------------------------

!> Setup for the QG model's background error standard deviation matrix

subroutine c_qg_bstddev_setup(c_key_self, c_conf) bind (c,name='qg_bstddev_setup_f90')

use iso_c_binding
use qg_bstddev_mod

implicit none
integer(c_int), intent(inout) :: c_key_self   !< The background std dev structure
type(c_ptr), intent(in)    :: c_conf     !< The configuration
type(qg_3d_bstddev_config), pointer :: self

call qg_3d_bstddev_registry%init()
call qg_3d_bstddev_registry%add(c_key_self)
call qg_3d_bstddev_registry%get(c_key_self, self)

call qg_3d_bstddev_setup(c_conf, self)

end subroutine c_qg_bstddev_setup

! ------------------------------------------------------------------------------
!> Delete for the QG model's background error std dev matrix

subroutine c_qg_bstddev_delete(c_key_self) bind (c,name='qg_bstddev_delete_f90')

use iso_c_binding
use qg_bstddev_mod

implicit none
integer(c_int), intent(inout) :: c_key_self  !< The background std dev structure
type(qg_3d_bstddev_config), pointer :: self

call qg_3d_bstddev_registry%get(c_key_self,self)
call qg_3d_bstddev_delete(self)
call qg_3d_bstddev_registry%remove(c_key_self)

end subroutine c_qg_bstddev_delete

! ------------------------------------------------------------------------------

!> Multiply streamfunction by inverse of std dev

subroutine c_qg_bstddev_inv_mult(c_key_conf, c_key_in, c_key_out) bind(c,name='qg_bstddev_invmult_f90')

use iso_c_binding
use qg_bstddev_mod
use qg_fields
use kinds

implicit none
integer(c_int), intent(in) :: c_key_conf  !< stddev config structure
integer(c_int), intent(in) :: c_key_in    !< Streamfunction: psi
integer(c_int), intent(in) :: c_key_out   !< Streamfunction: psi
type(qg_3d_bstddev_config), pointer :: conf
type(qg_field), pointer :: xin
type(qg_field), pointer :: xout

call qg_3d_bstddev_registry%get(c_key_conf,conf)
call qg_field_registry%get(c_key_in,xin)
call qg_field_registry%get(c_key_out,xout)

call zeros(xout)
call qg_3d_bstddev_inv_mult(xin,xout,conf)

end subroutine c_qg_bstddev_inv_mult

! ------------------------------------------------------------------------------

!> Multiply streamfunction by std dev

subroutine c_qg_bstddev_mult(c_key_conf, c_key_in, c_key_out) bind(c,name='qg_bstddev_mult_f90')

use iso_c_binding
use qg_bstddev_mod
use qg_fields
use kinds

implicit none
integer(c_int), intent(in) :: c_key_conf  !< stddev config structure
integer(c_int), intent(in) :: c_key_in    !< Streamfunction: psi
integer(c_int), intent(in) :: c_key_out   !< Streamfunction: psi
type(qg_3d_bstddev_config), pointer :: conf
type(qg_field), pointer :: xin
type(qg_field), pointer :: xout

call qg_3d_bstddev_registry%get(c_key_conf,conf)
call qg_field_registry%get(c_key_in,xin)
call qg_field_registry%get(c_key_out,xout)

call zeros(xout)
call qg_3d_bstddev_mult(xin,xout,conf)

end subroutine c_qg_bstddev_mult

! ------------------------------------------------------------------------------
