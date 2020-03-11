!
! (C) Copyright 2019-2020 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
!

#include "fckit/fctest.h"

!> Test interpolation interface in oops

!-------------------------------------------------------------------------------
! Manage test suite

TESTSUITE(test_util_interpinterface)

TESTSUITE_INIT
  use fckit_module
  use liboops_mod
  implicit none
  call liboops_initialise()
  call fckit_main%init()
END_TESTSUITE_INIT

TESTSUITE_FINALIZE
  use fckit_module
  use liboops_mod
  implicit none
  call fckit_main%final()
  call liboops_finalise()
END_TESTSUITE_FINALIZE

!-------------------------------------------------------------------------------
! define tests
!-------------------------------------------------------------------------------
!> Test application of interpolator
!

TEST(test_interpolator_apply)

use fckit_configuration_module, only : fckit_configuration
use fckit_module
use fckit_pathname_module, only : fckit_pathname
use kinds
use bump_interpolation_mod
use random_mod

type(bump_interpolator) :: interpolator
type(fckit_configuration) :: config
type(fckit_mpi_comm) :: comm
real(kind_real), parameter :: lat_range(2) = (/ -90, 90 /)
real(kind_real), parameter :: lon_range(2) = (/ 0, 360 /)
integer, parameter :: rseed = 8
integer, parameter :: nlat1 = 4, nlon1 = 8
integer, parameter :: nlat2 = 3, nlon2 = 6
real(kind_real) :: lat1(nlat1,nlon1), lon1(nlat1,nlon1)
real(kind_real) :: lat2(nlat2,nlon2), lon2(nlat2,nlon2)
real(kind_real) :: field1(nlat1,nlon1)
real(kind_real) :: field2_grid(nlat2,nlon2)
real(kind_real) :: field2_latlon(nlat2,nlon2)
real(kind_real), allocatable :: buffer(:)
character(len=:), allocatable :: filename
integer :: ii, jj
type(bump_grid) :: grid1, grid2
type(bump_field) :: infield, outfield

call fckit_log%info("Starting test_interpolator_apply")

comm = fckit_mpi_comm()

call fckit_resource("-config", "", filename)
config = fckit_YAMLConfiguration(fckit_pathname(filename))

! fill input arrays with some random lats, longitudes, and field values
allocate(buffer(nlat1))
do jj = 1, nlon1
   call uniform_distribution(buffer,lat_range(1),lat_range(2),seed=rseed,sort=.true.)
   lat1(:,jj) = buffer
enddo

deallocate(buffer)
allocate(buffer(nlon1))

do ii = 1, nlat1
   call uniform_distribution(buffer,lon_range(1),lon_range(2),seed=rseed,sort=.true.)
   lon1(ii,:) = buffer
   call normal_distribution(buffer, 0.0_kind_real, 1.0_kind_real, rseed)
   field1(ii,:) = buffer
enddo

! now define output grid
deallocate(buffer)
allocate(buffer(nlat2))
do jj = 1, nlon2
   call uniform_distribution(buffer,lat_range(1),lat_range(2),seed=rseed,sort=.true.)
   lat2(:,jj) = buffer
enddo

deallocate(buffer)
allocate(buffer(nlon2))

do ii = 1, nlat2
   call uniform_distribution(buffer,lon_range(1),lon_range(2),seed=rseed,sort=.true.)
   lon2(ii,:) = buffer
enddo

!------------------------
! grid init
call grid1%create(lat1, lon1)
call grid2%create(lat2, lon2)
call interpolator%init(comm, config, grid1, grid2)

! grid apply
call infield%create(field1)
call outfield%create(field2_grid, allocate_only = .true.)
call interpolator%apply(infield, outfield)
call outfield%result(field2_grid)

!------------------------
! latlon init and apply
call interpolator%init(comm, config, lat1, lon1, lat2, lon2)
call interpolator%apply(field1, field2_latlon)

!------------------------
! stop here for now.  When we have a real interpolator implemented, we can
! compare field2_grid to field2_latlon to verify that all the interfaces
! are functioning properly

FCTEST_CHECK_EQUAL(size(field2_grid), size(field2_latlon))

!------------------------
! clean up
call grid1%delete()
call grid2%delete()
call infield%delete()
call outfield%delete()
call interpolator%delete()
   
call fckit_log%info("Leaving test_interpolator_apply")

END_TEST

!-------------------------------------------------------------------------------

END_TESTSUITE
