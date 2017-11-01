! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Fortran module for streamfunction observations for the QG model
module qg_stream_mod

use iso_c_binding
use config_mod
use duration_mod
use qg_obs_data
use qg_obs_vectors
use qg_obsoper_mod
use qg_vars_mod
use qg_locs_mod
use qg_goms_mod
use kinds

implicit none
private

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine c_qg_stream_setup(c_key_self, c_conf) bind(c,name='qg_stream_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self
type(c_ptr), intent(in)       :: c_conf

type(qg_obsoper), pointer :: self
character(len=1) :: svars(1) = (/"x"/)

call qg_obsoper_registry%init()
call qg_obsoper_registry%add(c_key_self)
call qg_obsoper_registry%get(c_key_self, self)

call qg_oper_setup(self, c_conf, svars, 1)

end subroutine c_qg_stream_setup

! ------------------------------------------------------------------------------

subroutine c_qg_stream_delete(c_key_self) bind(c,name='qg_stream_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_self

type(qg_obsoper), pointer :: self

call qg_obsoper_registry%get(c_key_self, self)
call qg_oper_delete(self)
call qg_obsoper_registry%remove(c_key_self)

end subroutine c_qg_stream_delete

! ------------------------------------------------------------------------------

subroutine qg_stream_equiv(c_key_gom, c_key_hofx, c_bias) &
 & bind(c,name='qg_stream_equiv_f90')
implicit none
integer(c_int), intent(in) :: c_key_gom
integer(c_int), intent(in) :: c_key_hofx
real(c_double), intent(in) :: c_bias
type(qg_goms), pointer  :: gom
type(obs_vect), pointer :: hofx
integer :: io, jo

call qg_goms_registry%get(c_key_gom, gom) 
call qg_obs_vect_registry%get(c_key_hofx,hofx)

do jo=1,gom%nobs
  io=gom%indx(jo)
  hofx%values(1,io)=gom%values(1,jo) + c_bias
enddo

end subroutine qg_stream_equiv

! ------------------------------------------------------------------------------

subroutine qg_stream_equiv_tl(c_key_gom, c_key_hofx, c_bias) &
 & bind(c,name='qg_stream_equiv_tl_f90')
implicit none
integer(c_int), intent(in) :: c_key_gom
integer(c_int), intent(in) :: c_key_hofx
real(c_double), intent(in) :: c_bias
type(qg_goms), pointer  :: gom
type(obs_vect), pointer :: hofx
integer :: io, jo

call qg_goms_registry%get(c_key_gom, gom)
call qg_obs_vect_registry%get(c_key_hofx,hofx)

do jo=1,gom%nobs
  io=gom%indx(jo)
  hofx%values(1,io)=gom%values(1,jo) + c_bias
enddo

end subroutine qg_stream_equiv_tl

! ------------------------------------------------------------------------------

subroutine qg_stream_equiv_ad(c_key_gom, c_key_hofx, c_bias) &
 & bind(c,name='qg_stream_equiv_ad_f90')
implicit none
integer(c_int), intent(in) :: c_key_gom
integer(c_int), intent(in) :: c_key_hofx
real(c_double), intent(inout) :: c_bias
type(qg_goms), pointer  :: gom
type(obs_vect), pointer :: hofx
integer :: io, jo

call qg_goms_registry%get(c_key_gom, gom)
call qg_obs_vect_registry%get(c_key_hofx,hofx)

do jo=1,gom%nobs
  io=gom%indx(jo)
  gom%values(1,jo)=hofx%values(1,io)
  c_bias = c_bias + hofx%values(1,io)
enddo

end subroutine qg_stream_equiv_ad

! ------------------------------------------------------------------------------

end module qg_stream_mod
