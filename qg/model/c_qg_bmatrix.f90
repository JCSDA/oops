! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

! ------------------------------------------------------------------------------

!> Setup for the QG model's background error covariance matrix

subroutine c_qg_b_setup(c_key_self, c_conf, c_key_geom) &
          & bind (c,name='qg_b_setup_f90')

use iso_c_binding
use qg_covariance_mod
use qg_geom_mod

implicit none
integer(c_int), intent(inout) :: c_key_self   !< The background covariance structure
type(c_ptr), intent(in)    :: c_conf     !< The configuration
integer(c_int), intent(in) :: c_key_geom !< Geometry
type(qg_3d_covar_config), pointer :: self
type(qg_geom),  pointer :: geom

call qg_geom_registry%get(c_key_geom, geom)
call qg_3d_cov_registry%init()
call qg_3d_cov_registry%add(c_key_self)
call qg_3d_cov_registry%get(c_key_self, self)

call qg_3d_covar_setup(c_conf, geom, self)

end subroutine c_qg_b_setup

! ------------------------------------------------------------------------------
!> Delete for the QG model's background error covariance matrix

subroutine c_qg_b_delete(c_key_self) bind (c,name='qg_b_delete_f90')

use iso_c_binding
use qg_covariance_mod

implicit none
integer(c_int), intent(inout) :: c_key_self  !< The background covariance structure
type(qg_3d_covar_config), pointer :: self

call qg_3d_cov_registry%get(c_key_self,self)
call qg_3d_covar_delete(self)
call qg_3d_cov_registry%remove(c_key_self)

end subroutine c_qg_b_delete

! ------------------------------------------------------------------------------

!> Multiply streamfunction by inverse of covariance

subroutine c_qg_b_inv_mult(c_key_conf, c_key_in, c_key_out) bind(c,name='qg_b_invmult_f90')

use iso_c_binding
use qg_covariance_mod
use qg_fields
use kinds

implicit none
integer(c_int), intent(in) :: c_key_conf  !< covar config structure
integer(c_int), intent(in) :: c_key_in    !< Streamfunction: psi
integer(c_int), intent(in) :: c_key_out   !< Streamfunction: psi
type(qg_3d_covar_config), pointer :: conf
type(qg_field), pointer :: xin
type(qg_field), pointer :: xout
real(kind=kind_real), allocatable :: xctl(:,:,:) ! Control vector

call qg_3d_cov_registry%get(c_key_conf,conf)
call qg_field_registry%get(c_key_in,xin)
call qg_field_registry%get(c_key_out,xout)

allocate(xctl(conf%nx, conf%ny, 2))
xctl(:,:,:)=0.0_kind_real

call qg_3d_covar_sqrt_inv_mult(conf%nx,conf%ny,xctl,xin,conf)
call zeros(xout)
call qg_3d_covar_sqrt_inv_mult_ad(conf%nx,conf%ny,xctl,xout,conf)

deallocate(xctl)

end subroutine c_qg_b_inv_mult

! ------------------------------------------------------------------------------

!> Multiply streamfunction by covariance

subroutine c_qg_b_mult(c_key_conf, c_key_in, c_key_out) bind(c,name='qg_b_mult_f90')

use iso_c_binding
use qg_covariance_mod
use qg_fields
use kinds

implicit none
integer(c_int), intent(in) :: c_key_conf  !< covar config structure
integer(c_int), intent(in) :: c_key_in    !< Streamfunction: psi
integer(c_int), intent(in) :: c_key_out   !< Streamfunction: psi
type(qg_3d_covar_config), pointer :: conf
type(qg_field), pointer :: xin
type(qg_field), pointer :: xout
real(kind=kind_real), allocatable :: xctl(:,:,:) ! Control vector

call qg_3d_cov_registry%get(c_key_conf,conf)
call qg_field_registry%get(c_key_in,xin)
call qg_field_registry%get(c_key_out,xout)

allocate(xctl(conf%nx, conf%ny, 2))

xctl(:,:,:)=0.0_kind_real
call qg_3d_covar_sqrt_mult_ad(conf%nx,conf%ny,xin,xctl,conf)
call zeros(xout)
call qg_3d_covar_sqrt_mult(conf%nx,conf%ny,xout,xctl,conf)

deallocate(xctl)

end subroutine c_qg_b_mult

! ------------------------------------------------------------------------------

!> Generate randomized increment

subroutine c_qg_b_randomize(c_key_conf, c_key_out) bind(c,name='qg_b_randomize_f90')

use iso_c_binding
use qg_covariance_mod
use qg_fields
use random_mod
use kinds

implicit none
integer(c_int), intent(in) :: c_key_conf  !< covar config structure
integer(c_int), intent(in) :: c_key_out   !< Randomized increment
type(qg_3d_covar_config), pointer :: conf
type(qg_field), pointer :: xout
real(kind=kind_real), allocatable :: xctl(:,:,:) ! Control vector
integer :: rseed = 7 ! for reproducible tests

call qg_3d_cov_registry%get(c_key_conf,conf)
call qg_field_registry%get(c_key_out,xout)

allocate(xctl(conf%nx, conf%ny, 2))
call normal_distribution(xctl,0.0_kind_real,1.0_kind_real,rseed)
call zeros(xout)
call qg_3d_covar_sqrt_mult(conf%nx,conf%ny,xout,xctl,conf)

deallocate(xctl)

end subroutine c_qg_b_randomize

! ------------------------------------------------------------------------------
