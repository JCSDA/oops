! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Structure holding configuration variables for the QG model

module qg_configs

use kinds
implicit none
private
public :: qg_config
public :: qg_config_registry

!> Fortran derived type to hold configuration data for the QG model
type :: qg_config
  integer :: nx     !< Zonal grid dimension
  integer :: ny     !< Meridional grid dimension
  ! dimensional parameters
  real(kind=kind_real) :: dt0       !< dimensional time (seconds)
  real(kind=kind_real) :: d1        !< depth of top layer (m)
  real(kind=kind_real) :: d2        !< depth of bottom layer (m)
  real(kind=kind_real) :: deltax0   !< zonal grid spacing (m)
  real(kind=kind_real) :: deltay0   !< meridional grid spacing (m)
  ! non-dimensional parameters
  real(kind=kind_real) :: dt        !< non-dimensional time
  real(kind=kind_real) :: f1        !< Coefficient of PV operator
  real(kind=kind_real) :: f2        !< Coefficient of PV operator
  real(kind=kind_real) :: deltax    !< non-dimensional zonal grid spacing
  real(kind=kind_real) :: deltay    !< non-dimensional meridional grid spacing
  real(kind=kind_real) :: rsmax     !< max non-dimensional height of orography
  real(kind=kind_real), allocatable :: rs(:,:) !< non-dimensional orography
end type qg_config

#define LISTED_TYPE qg_config

!> Linked list interface - defines registry_t type
#include "util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_config_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "util/linkedList_c.f"

end module qg_configs
