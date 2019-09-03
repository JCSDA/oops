! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_geom_interface

use fckit_configuration_module, only: fckit_configuration
use fckit_log_module,only: fckit_log
use kinds
use iso_c_binding
use qg_geom_mod

implicit none

private
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Setup geometry
subroutine qg_geom_setup_c(c_key_self,c_conf) bind(c,name='qg_geom_setup_f90')

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Geometry
type(c_ptr),intent(in) :: c_conf           !< Configuration

! Local variables
type(fckit_configuration) :: f_conf
type(qg_geom),pointer :: self

! Interface
f_conf = fckit_configuration(c_conf)
call qg_geom_registry%init()
call qg_geom_registry%add(c_key_self)
call qg_geom_registry%get(c_key_self,self)

! Call Fortran
call qg_geom_setup(self,f_conf)

end subroutine qg_geom_setup_c
! ------------------------------------------------------------------------------
!> Clone geometry
subroutine qg_geom_clone_c(c_key_self,c_key_other) bind(c,name='qg_geom_clone_f90')

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Geometry
integer(c_int),intent(in) :: c_key_other   !< Other geometry

! Local variables
type(qg_geom),pointer :: self,other

! Interface
call qg_geom_registry%get(c_key_other,other)
call qg_geom_registry%init()
call qg_geom_registry%add(c_key_self)
call qg_geom_registry%get(c_key_self ,self )

! Call Fortran
call qg_geom_clone(self,other)

end subroutine qg_geom_clone_c

! ------------------------------------------------------------------------------
!> Delete geometry
subroutine qg_geom_delete_c(c_key_self) bind(c,name='qg_geom_delete_f90')

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Geometry

! Local variables
type(qg_geom),pointer :: self

! Interface
call qg_geom_registry%get(c_key_self,self)

! Call Fortran
call qg_geom_delete(self)

! Clear interface
call qg_geom_registry%remove(c_key_self)

end subroutine qg_geom_delete_c
! ------------------------------------------------------------------------------
!> Get geometry info
subroutine qg_geom_info_c(c_key_self,c_nx,c_ny,c_nz,c_deltax,c_deltay) bind(c,name='qg_geom_info_f90')

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< Geometry
integer(c_int),intent(inout) :: c_nx     !< Number of points in the zonal direction
integer(c_int),intent(inout) :: c_ny     !< Number of points in the meridional direction
integer(c_int),intent(inout) :: c_nz     !< Number of vertical levels
real(c_double),intent(inout) :: c_deltax !< Zonal cell size
real(c_double),intent(inout) :: c_deltay !< Meridional cell size

! Local variables
type(qg_geom),pointer :: self

! Interface
call qg_geom_registry%get(c_key_self,self)

! Call Fortran
call qg_geom_info(self,c_nx,c_ny,c_nz,c_deltax,c_deltay)

end subroutine qg_geom_info_c
! ------------------------------------------------------------------------------
end module qg_geom_interface
