! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Fortran module handling observation locations

module qg_locs_mod

use iso_c_binding
use qg_obs_vectors
use kinds

implicit none
private
public :: qg_locs, qg_loc_setup
public :: qg_locs_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold observation locations
type :: qg_locs
  integer :: nloc
  real(kind=kind_real), allocatable :: xyz(:,:)
  integer, allocatable :: indx(:)
end type qg_locs

#define LISTED_TYPE qg_locs

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_locs_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine c_qg_loc_create(c_key_locs) bind(c,name='qg_loc_create_f90')

implicit none
integer(c_int), intent(inout) :: c_key_locs

call qg_locs_registry%init()
call qg_locs_registry%add(c_key_locs)

end subroutine c_qg_loc_create

! ------------------------------------------------------------------------------
!> Generate locations for interpolation test
!!
!! \details **c_qg_loc_test()** generates a list of user-specified and/or
!! randomly-generated locations.  It was originally intended to be used with the
!! test::testStateInterpolation() routine but it may also be adapted for other
!! applications.  The user can specify a list of specific locations to be tested
!! (optional) by setting the "lats" and "lons" items in the "StateTest" section
!! of the config file.  Alternatively or in addition to these specific locations,
!! the user may request that **Nrandom** random locations be generated.  If
!! neither lats/lons nor Nrandom is specified in the config file, then two
!! random locations are generated.
!!
!! \date April 2, 2018: Created (M. Miesch, JCSDA)
!!
!! \warning the **lats** and **lons** arrays in the input file are assumed to be
!! degrees (for uniformity relative to other models).  These are converted to normalized
!! x and y coordinates when they are stored in the LocationsQG object as defined here.
!! The height specified in the input file is assumed to be normalized already, so
!! it will take on a value between 0 and 1.
!!
!! \warning Since the interpolate() member function of State objects does not
!! interpolate in height (z). the z coordinate as recorded in the LocationsQG object
!! refers to an integer level of 1 or 2.  If the user does not explicitly specify
!! the (normalized) height for each latitude and longitude in the config file,
!! then the level is set to a default value of 1.
!!
!! \sa qg::LocationsQG
!!
subroutine c_qg_loc_test(c_key_locs,config,klocs,klats,klons,kz) bind(c,name='qg_loc_test_f90')
use config_mod
use fckit_log_module, only : fckit_log
  
implicit none
integer(c_int), intent(in) :: c_key_locs   !< key to F90 Locations object
type(c_ptr)   , intent(in) :: config       !< Configuration (typically State.StateGenerate)
integer(c_int), intent(in) :: klocs        !< Number of user-specified locations
real(c_double), intent(in) :: klats(klocs) !< user-specified latitudes (degrees)
real(c_double), intent(in) :: klons(klocs) !< user-specified longitudes (degrees)
real(c_double), intent(in) :: kz(klocs)    !< user-specified heights (normalized between 0-1)

type(qg_locs), pointer ::locs
real(kind_real), allocatable :: xx(:), yy(:), rnum(:)
integer :: nrand, nloc, i, jo, nseed
integer*4 :: rseed0
integer*4, allocatable :: rseed(:)

call fckit_log%warning("qg_locs_mod:qg_loc_test generating test locations")

if (config_element_exists(config,"Nrandom")) then
   nrand = config_get_int(config,"Nrandom")
else
   nrand = 0
endif

if (klocs > 0) then
   nloc = klocs + nrand
else
   nloc = nrand
endif

! pick 2 random locations if no locations are specified
if (nloc < 1) then
   nrand = 2
   nloc = nrand
endif

allocate(xx(nloc),yy(nloc))

! convert lat and lon to normalized x,y coordinate between 0 and 1   
! this is what interp_tl() wants
!!
!> \warning **qg_loc_test()** latitudes are converted to normalized locations
!! based on on an assumed latitudinal extent of 180 degrees 
!!

if (klocs > 0) then
   xx(1:klocs) = modulo((klons(:)+180.0_kind_real),360.0_kind_real)/360.0_kind_real
   yy(1:klocs) = (klats(:) + 90.0_kind_real)/180.0_kind_real
endif

if (config_element_exists(config,"random_seed")) then
   ! read in the (optional) seed as a real for higher precision
   ! included for reproducibility
   rseed0 = config_get_real(config,"random_seed")
else
   ! get the seed from the system clock
   call system_clock(count=rseed0)
endif

! define an optionally reproducible random number seed
call random_seed(size=nseed)
allocate(rseed(nseed))
do i=1,nseed
   rseed(i) = rseed0 + i**2
enddo
call random_seed(put=rseed)

allocate(rnum(3*nrand))
call random_number(rnum)

if (nrand > 0) then   
   xx(klocs+1:nloc) = rnum(1:nrand)
   yy(klocs+1:nloc) = rnum(nrand+1:2*nrand)
endif

! Now define the F90 locations object locs.  It's assumed that
! this object already exists in the registry.
call qg_locs_registry%get(c_key_locs,locs)

locs%nloc=nloc

! I think the index is simply the index of the locations array,
! ranging from 1 to nloc
allocate(locs%indx(nloc))

! Now initialize the lat and lon
allocate(locs%xyz(3,nloc))
do jo=1,nloc
   locs%xyz(1,jo)=xx(jo)
   locs%xyz(2,jo)=yy(jo)
   locs%indx(:)=jo
enddo

! The QG interpolation does not interpolate in height (z)
! So just define z as the midpoint of the level in question, in meters

do jo=1,klocs
   if (kz(jo) <= 0.5_kind_real) then
      locs%xyz(3,jo) = 1.0_kind_real
   else 
      locs%xyz(3,jo) = 2.0_kind_real
   endif
enddo

do jo=klocs+1,nloc
   i=2*nrand+jo-klocs
   locs%xyz(3,jo) = int(rnum(i)*2.0_kind_real)+1 
enddo

deallocate(xx,yy,rnum)

end subroutine c_qg_loc_test

! ------------------------------------------------------------------------------

subroutine qg_loc_setup(self, lvec, kobs)
implicit none
type(qg_locs), intent(inout) :: self
type(obs_vect), intent(in) :: lvec
integer, intent(in) :: kobs(:)
integer :: jc, jo

self%nloc=lvec%nobs
allocate(self%indx(self%nloc))
self%indx(:) = kobs(:)
allocate(self%xyz(3,self%nloc))
do jo=1,self%nloc
  do jc=1,3
    self%xyz(jc,jo)=lvec%values(jc,jo)
  enddo
enddo

end subroutine qg_loc_setup

! ------------------------------------------------------------------------------

subroutine c_qg_loc_delete(key) bind(c,name='qg_loc_delete_f90')

implicit none
integer(c_int), intent(inout) :: key
type(qg_locs), pointer :: self

call qg_locs_registry%get(key,self)
if (allocated(self%xyz)) deallocate(self%xyz)
call qg_locs_registry%remove(key)

end subroutine c_qg_loc_delete

! ------------------------------------------------------------------------------

subroutine c_qg_loc_nobs(key, kobs) bind(c,name='qg_loc_nobs_f90')

implicit none
integer(c_int), intent(in) :: key
integer(c_int), intent(inout) :: kobs
type(qg_locs), pointer :: self

call qg_locs_registry%get(key,self)
kobs = self%nloc

end subroutine c_qg_loc_nobs

! ------------------------------------------------------------------------------
subroutine c_qg_loc_element(key,iloc,xyz) bind(c,name='qg_loc_element_f90')

implicit none
integer(c_int), intent(in)    :: key
integer(c_int), intent(in)    :: iloc
real(c_double), intent(inout) :: xyz(3)
type(qg_locs), pointer :: self

call qg_locs_registry%get(key,self)

xyz(1) = self%xyz(2,iloc+1) ! latitude
xyz(2) = self%xyz(1,iloc+1) ! longitude
xyz(3) = self%xyz(3,iloc+1) ! height

end subroutine c_qg_loc_element

! ------------------------------------------------------------------------------

end module qg_locs_mod
