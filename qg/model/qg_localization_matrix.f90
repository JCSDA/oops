! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module QGLOCALIZATION
use iso_c_binding
use qg_fields
use kinds
use qg_constants
use qg_geom_mod
use config_mod
use qg_covariance_mod

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------

subroutine qg_localization_mult(c_key_conf, c_key_xincr) bind(c,name='qg_localization_mult_f90')
implicit none
integer(c_int), intent(in) :: c_key_conf
integer(c_int), intent(in) :: c_key_xincr

type(qg_3d_covar_config), pointer :: conf   !< Config structure
type(qg_field), pointer :: xincr
real(kind=kind_real), allocatable :: xctl(:,:,:) ! Control vector

call qg_3d_cov_registry%get(c_key_conf,conf)
call qg_field_registry%get(c_key_xincr,xincr)

allocate(xctl(conf%nx, conf%ny, 2))

xctl(:,:,:)=0.0_kind_real
call qg_3d_covar_sqrt_mult_ad(conf%nx,conf%ny,xincr,xctl,conf)
call zeros(xincr)
call qg_3d_covar_sqrt_mult(conf%nx,conf%ny,xincr,xctl,conf)

deallocate(xctl)
end subroutine qg_localization_mult

! ------------------------------------------------------------------------------

subroutine qg_localization_setup(c_key_conf, c_model, c_key_geom) bind(c,name='qg_localization_setup_f90')
implicit none
integer(c_int), intent(inout) :: c_key_conf
type(c_ptr), intent(in)    :: c_model  !< The configuration
integer(c_int), intent(in) :: c_key_geom !< Geometry
type(qg_3d_covar_config), pointer :: conf !< covar structure
type(qg_geom), pointer :: geom     !< Geometry

call qg_3d_cov_registry%init()
call qg_3d_cov_registry%add(c_key_conf)
call qg_3d_cov_registry%get(c_key_conf, conf)
call qg_geom_registry%get(c_key_geom, geom)
call qg_3d_covar_setup(c_model, geom, conf)

return
end subroutine qg_localization_setup

! ------------------------------------------------------------------------------

subroutine qg_localization_delete(c_key_conf) bind(c,name='qg_localization_delete_f90')
implicit none
integer(c_int), intent(inout) :: c_key_conf !< The model covariance structure

call qg_3d_covar_delete(c_key_conf)

end subroutine qg_localization_delete

! ------------------------------------------------------------------------------

end module QGLOCALIZATION
