! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Structure holding configuration variables for the 3d error
!! std dev matrices of the QG analysis.

module qg_bstddev_mod

use kinds
implicit none

!> Fortran derived type to hold configuration data for the QG background/model std dev
type :: qg_3d_bstddev_config
  real(kind=kind_real)    :: sigma        !< Standard deviation
end type qg_3d_bstddev_config

#define LISTED_TYPE qg_3d_bstddev_config

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_3d_bstddev_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

!> Setup for the QG model's 3d error std dev matrices (B and Q_i)

!> This routine queries the configuration for the parameters that define the
!! std dev matrix, and stores the relevant values in the
!! error std dev structure.

subroutine qg_3d_bstddev_setup(c_model, config)

use qg_constants
use qg_geom_mod
use iso_c_binding
use config_mod
use kinds
use fckit_log_module, only : fckit_log

implicit none
type(c_ptr), intent(in)   :: c_model  !< The configuration
type(qg_3d_bstddev_config), intent(inout) :: config !< The std dev structure

config%sigma      = config_get_real(c_model,"standard_deviation")

return
end subroutine qg_3d_bstddev_setup

! ------------------------------------------------------------------------------

!> Delete for the QG model's 3d error std dev matrices

subroutine qg_3d_bstddev_delete(self)
implicit none
type(qg_3d_bstddev_config) :: self

end subroutine qg_3d_bstddev_delete

! ------------------------------------------------------------------------------

!> Multiply by inverse of std dev matrix

subroutine qg_3d_bstddev_inv_mult(xin,xout,config)
use iso_c_binding
use kinds
use qg_fields

implicit none
type(qg_field), intent(in)    :: xin
type(qg_field), intent(inout) :: xout
type(qg_3d_bstddev_config), intent(in) :: config !< bstddev config structure

integer :: i, j, k
real(kind=kind_real) :: zc

!--- multiply by inverse standard deviation

zc = 1.0_kind_real/config%sigma
do k=1,2
  do j=1,xout%geom%ny
    do i=1,xout%geom%nx
      xout%x(i,j,k) = zc * xin%x(i,j,k)
    enddo
  enddo
enddo

end subroutine qg_3d_bstddev_inv_mult

! ------------------------------------------------------------------------------

!> Multiply by std dev matrix

subroutine qg_3d_bstddev_mult(xin,xout,config)
use iso_c_binding
use kinds
use qg_fields

implicit none
type(qg_field), intent(in)    :: xin
type(qg_field), intent(inout) :: xout
type(qg_3d_bstddev_config), intent(in) :: config !< bstddev config structure

integer :: i, j, k

!--- multiply by standard deviation

do k=1,2
  do j=1,xout%geom%ny
    do i=1,xout%geom%nx
      xout%x(i,j,k) = config%sigma * xin%x(i,j,k)
    enddo
  enddo
enddo

end subroutine qg_3d_bstddev_mult

! ------------------------------------------------------------------------------

end module qg_bstddev_mod
