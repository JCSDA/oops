! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_locs_interface

use config_mod
use fckit_log_module, only: fckit_log
use iso_c_binding
use kinds
use random_mod
use qg_constants_mod
use qg_geom_mod
use qg_locs_mod
use qg_projection_mod

implicit none
private
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Create locations
subroutine qg_locs_create_c(c_key_self) bind(c,name='qg_locs_create_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Locations

! Interface
call qg_locs_registry%init()
call qg_locs_registry%add(c_key_self)

end subroutine qg_locs_create_c
! ------------------------------------------------------------------------------
!> Generate test locations
subroutine qg_locs_test_c(c_key_self,c_conf,klocs,klons,klats,kz) bind(c,name='qg_locs_test_f90')
  
implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self   !< Locations
type(c_ptr),intent(in) :: c_conf          !< Configuration
integer(c_int),intent(in) :: klocs        !< Number of user-specified locations
real(c_double),intent(in) :: klats(klocs) !< User-specified latitudes (degrees)
real(c_double),intent(in) :: klons(klocs) !< User-specified longitudes (degrees)
real(c_double),intent(in) :: kz(klocs)    !< User-specified altitudes (m)

! Local variables
type(qg_locs),pointer :: self

! Interface
call qg_locs_registry%get(c_key_self,self)

! Call Fortran
call qg_locs_test(self,c_conf,klocs,klons,klats,kz)

end subroutine qg_locs_test_c
! ------------------------------------------------------------------------------
!> Delete locations
subroutine qg_locs_delete_c(c_key_self) bind(c,name='qg_locs_delete_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Locations

! Local variables
type(qg_locs),pointer :: self

! Interface
call qg_locs_registry%get(c_key_self,self)

! Call Fortran
call qg_locs_delete(self)

! Clear interface
call qg_locs_registry%remove(c_key_self)

end subroutine qg_locs_delete_c
! ------------------------------------------------------------------------------
!> Get number of observations
subroutine qg_locs_nobs_c(c_key_self,kobs) bind(c,name='qg_locs_nobs_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Locations
integer(c_int),intent(inout) :: kobs    !< Number of observations

! Local variables
type(qg_locs),pointer :: self

! Interface
call qg_locs_registry%get(c_key_self,self)

! Call Fortran
call qg_locs_nobs(self,kobs)

end subroutine qg_locs_nobs_c
! ------------------------------------------------------------------------------
!> Get location element coordinates
subroutine qg_locs_element_c(c_key_self,iloc,lon,lat,z) bind(c,name='qg_locs_element_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Locations
integer(c_int),intent(in) :: iloc       !< Index
real(c_double),intent(inout) :: lon     !< Longitude
real(c_double),intent(inout) :: lat     !< Latitude
real(c_double),intent(inout) :: z       !< Altitude

! Local variables
type(qg_locs),pointer :: self

! Interface
call qg_locs_registry%get(c_key_self,self)

! Call Fortran
call qg_locs_element(self,iloc,lon,lat,z)

end subroutine qg_locs_element_c
! ------------------------------------------------------------------------------
end module qg_locs_interface
