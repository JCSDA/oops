! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_locs_mod

use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log
use iso_c_binding
use kinds
use random_mod
use qg_constants_mod
use qg_geom_mod
use qg_obsvec_mod
use qg_projection_mod

implicit none
private
public :: qg_locs
public :: qg_locs_registry
public :: qg_locs_test,qg_locs_delete,qg_locs_nobs,qg_locs_element,qg_locs_from_obsvec
! ------------------------------------------------------------------------------
integer,parameter :: rseed = 1 !< Random seed (for reproducibility)

type :: qg_locs
  integer :: nlocs                      !< Number of locations
  real(kind_real),allocatable :: lon(:) !< Longitudes
  real(kind_real),allocatable :: lat(:) !< Latitudes
  real(kind_real),allocatable :: z(:)   !< Altitudes
  integer,allocatable :: indx(:)        !< Index
end type qg_locs

#define LISTED_TYPE qg_locs

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_locs_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------
!> Generate test locations
subroutine qg_locs_test(self,f_conf,klocs,klons,klats,kz)
  
implicit none

! Passed variables
type(qg_locs),intent(inout) :: self            !< Locations
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration
integer,intent(in) :: klocs                    !< Number of user-specified locations
real(kind_real),intent(in) :: klats(klocs)     !< User-specified latitudes (degrees)
real(kind_real),intent(in) :: klons(klocs)     !< User-specified longitudes (degrees)
real(kind_real),intent(in) :: kz(klocs)        !< User-specified altitudes (m)

! Local variables
integer :: jobs,Nrandom
real(kind_real),allocatable :: x(:),y(:)

! Get number of random locations
self%nlocs = klocs
if (f_conf%has("Nrandom")) then
  call f_conf%get_or_die("Nrandom",Nrandom)
  self%nlocs = self%nlocs+Nrandom
end if

! Check number of locations
if (self%nlocs<1) call abor1_ftn('qg_locs_test: no location')

! Allocation
allocate(self%indx(self%nlocs))
allocate(self%lon(self%nlocs))
allocate(self%lat(self%nlocs))
allocate(self%z(self%nlocs))

! Index
do jobs=1,self%nlocs
   self%indx(jobs)=jobs
enddo

! Copy specified locations
if (klocs>0) then
   self%lon(1:klocs) = klons
   self%lat(1:klocs) = klats
   self%z(1:klocs) = kz
endif

if (self%nlocs>klocs) then
  ! Allocation
  allocate(x(klocs+1:self%nlocs))
  allocate(y(klocs+1:self%nlocs))

  ! Generate random coordinates
  call uniform_distribution(x(klocs+1:self%nlocs),0.0_kind_real,domain_zonal,rseed)
  call uniform_distribution(y(klocs+1:self%nlocs),0.0_kind_real,domain_meridional,rseed)
  call uniform_distribution(self%z(klocs+1:self%nlocs),0.0_kind_real,domain_depth,rseed)

  ! Convert to lon/lat
  do jobs=klocs+1,self%nlocs
    call xy_to_lonlat(x(jobs),y(jobs),self%lon(jobs),self%lat(jobs))
  enddo
endif

end subroutine qg_locs_test
! ------------------------------------------------------------------------------
!> Delete locations
subroutine qg_locs_delete(self)

implicit none

! Passed variables
type(qg_locs),intent(inout) :: self !< Locations

! Release memory
if (allocated(self%lon)) deallocate(self%lon)
if (allocated(self%lat)) deallocate(self%lat)
if (allocated(self%z)) deallocate(self%z)

end subroutine qg_locs_delete
! ------------------------------------------------------------------------------
!> Get number of observations
subroutine qg_locs_nobs(self,kobs)

implicit none

! Passed variables
type(qg_locs),intent(inout) :: self  !< Locations
integer(c_int),intent(inout) :: kobs !< Number of observations

! Get number of observations
kobs = self%nlocs

end subroutine qg_locs_nobs
! ------------------------------------------------------------------------------
!> Get location element coordinates
subroutine qg_locs_element(self,iloc,lon,lat,z)

implicit none

! Passed variables
type(qg_locs),intent(in) :: self     !< Locations
integer,intent(in) :: iloc           !< Index
real(kind_real),intent(inout) :: lon !< Longitude
real(kind_real),intent(inout) :: lat !< Latitude
real(kind_real),intent(inout) :: z   !< Altitude

! Get location element coordinates
lon = self%lon(iloc+1)
lat = self%lat(iloc+1)
z = self%z(iloc+1)

end subroutine qg_locs_element
! ------------------------------------------------------------------------------
!> Setup locations from observation vector
subroutine qg_locs_from_obsvec(self,ovec,kobs)

implicit none

! Passed variables
type(qg_locs),intent(inout) :: self !< Locations
type(qg_obsvec),intent(in) :: ovec  !< Observation vector
integer,intent(in) :: kobs(:)       !< Observations index

! Local variables
integer :: jo

! Get number of observations
self%nlocs = ovec%nobs

! Allocation
allocate(self%indx(self%nlocs))
allocate(self%lon(self%nlocs))
allocate(self%lat(self%nlocs))
allocate(self%z(self%nlocs))

! Copy index
self%indx(:) = kobs(:)

! Copy coordinates
do jo=1,self%nlocs
  self%lon(jo) = ovec%values(1,jo)
  self%lat(jo) = ovec%values(2,jo)
  self%z(jo) = ovec%values(3,jo)
enddo

end subroutine qg_locs_from_obsvec
! ------------------------------------------------------------------------------
end module qg_locs_mod
