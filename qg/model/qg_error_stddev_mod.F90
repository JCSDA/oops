! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_error_stddev_mod

use config_mod
use iso_c_binding
use kinds
!$ use omp_lib
use qg_constants_mod
use qg_geom_mod
use qg_fields_mod

implicit none

private
public :: qg_error_stddev_config
public :: qg_error_stddev_registry
public :: qg_error_stddev_setup,qg_error_stddev_mult,qg_error_stddev_inv_mult
! ------------------------------------------------------------------------------
type :: qg_error_stddev_config
  real(kind_real) :: sigma !< Error standard deviation
end type qg_error_stddev_config

#define LISTED_TYPE qg_error_stddev_config

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_error_stddev_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------
!> Setup error standard deviation
subroutine qg_error_stddev_setup(c_model,conf)

! Passed variables
type(c_ptr),intent(in) :: c_model                  !< Configuration
type(qg_error_stddev_config),intent(inout) :: conf !< Error standard deviation

! Get value
conf%sigma = config_get_real(c_model,'standard_deviation')

end subroutine qg_error_stddev_setup
! ------------------------------------------------------------------------------
!> Multiply by error standard deviation
subroutine qg_error_stddev_mult(conf,fld_in,fld_out)

! Passed variables
type(qg_error_stddev_config),intent(in) :: conf !< Error standard deviation
type(qg_fields),intent(in) :: fld_in            !< Input field
type(qg_fields),intent(inout) :: fld_out        !< Output field

! Local variables
integer :: ix,iy,iz

! Initialization
call qg_fields_zero(fld_out)

! Multiply by error standard deviation
!$omp parallel do schedule(static) private(iz,iy,ix)
do iz=1,fld_out%geom%nz
  do iy=1,fld_out%geom%ny
    do ix=1,fld_out%geom%nx
      fld_out%gfld3d(ix,iy,iz) = fld_in%gfld3d(ix,iy,iz)*conf%sigma
    enddo
  enddo
enddo
!$omp end parallel do

end subroutine qg_error_stddev_mult
! ------------------------------------------------------------------------------
!> Multiply by inverse of error standard deviation
subroutine qg_error_stddev_inv_mult(conf,fld_in,fld_out)

! Passed variables
type(qg_error_stddev_config),intent(in) :: conf !< Error standard deviation
type(qg_fields),intent(in) :: fld_in            !< Input field
type(qg_fields),intent(inout) :: fld_out        !< Output field

! Local variables
integer :: ix,iy,iz

! Initialization
call qg_fields_zero(fld_out)

! Multiply by inverse of error standard deviation
!$omp parallel do schedule(static) private(iz,iy,ix)
do iz=1,fld_out%geom%nz
  do iy=1,fld_out%geom%ny
    do ix=1,fld_out%geom%nx
      fld_out%gfld3d(ix,iy,iz) = fld_in%gfld3d(ix,iy,iz)/conf%sigma
    enddo
  enddo
enddo
!$omp end parallel do

end subroutine qg_error_stddev_inv_mult
! ------------------------------------------------------------------------------
end module qg_error_stddev_mod
