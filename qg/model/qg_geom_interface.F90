! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_geom_interface

use atlas_module, only: atlas_fieldset, atlas_functionspace_nodecolumns
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module,only: fckit_log
use kinds
use iso_c_binding
use qg_projection_mod
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
type(c_ptr),value,intent(in) :: c_conf     !< Configuration

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
!> Set lon/lat field
subroutine qg_geom_set_lonlat_c(c_key_self,c_afieldset) bind(c,name='qg_geom_set_lonlat_f90')

! Passed variables
integer(c_int),intent(in) :: c_key_self     !< Geometry
type(c_ptr),intent(in),value :: c_afieldset !< Fieldset pointer

! Local variables
type(qg_geom),pointer :: self
type(atlas_fieldset) :: afieldset

! Interface
call qg_geom_registry%get(c_key_self,self)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call qg_geom_set_lonlat(self,afieldset)

end subroutine qg_geom_set_lonlat_c
! ------------------------------------------------------------------------------
!> Set function space pointer
subroutine qg_geom_set_functionspace_pointer_c(c_key_self,c_afunctionspace) &
 & bind(c,name='qg_geom_set_functionspace_pointer_f90')

! Passed variables
integer(c_int),intent(in) :: c_key_self          !< Geometry
type(c_ptr),intent(in),value :: c_afunctionspace !< Function space pointer

! Local variables
type(qg_geom),pointer :: self

! Interface
call qg_geom_registry%get(c_key_self,self)
self%afunctionspace = atlas_functionspace_nodecolumns(c_afunctionspace)

end subroutine qg_geom_set_functionspace_pointer_c
! ------------------------------------------------------------------------------
!> Fill geometry fields
subroutine qg_geom_fill_geometry_fields_c(c_key_self,c_afieldset) &
 & bind(c,name='qg_geom_fill_geometry_fields_f90')

! Passed variables
integer(c_int),intent(in) :: c_key_self     !< Geometry
type(c_ptr),intent(in),value :: c_afieldset !< Fieldset pointer

! Local variables
type(qg_geom),pointer :: self
type(atlas_fieldset) :: afieldset

! Interface
call qg_geom_registry%get(c_key_self,self)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call qg_geom_fill_geometry_fields(self,afieldset)

end subroutine qg_geom_fill_geometry_fields_c
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
!> Get dimensions of computational domain
subroutine qg_geom_dimensions_c(lonmin, lonmax, latmin, latmax, zmax) bind(c,name='qg_geom_dimensions_f90')

    use qg_constants_mod
    implicit none

    ! Passed variables
    real(c_double), intent(inout) :: lonmin  !< longitude min
    real(c_double), intent(inout) :: lonmax  !< longitude max
    real(c_double), intent(inout) :: latmin  !< latitude min
    real(c_double), intent(inout) :: latmax  !< latitude max
    real(c_double), intent(inout) :: zmax    !< vertical max

    call xy_to_lonlat(0.0_kind_real, 0.0_kind_real, lonmin, latmin)
    call xy_to_lonlat(domain_zonal, domain_meridional, lonmax, latmax)

    zmax = domain_depth

end subroutine qg_geom_dimensions_c

! ------------------------------------------------------------------------------
end module qg_geom_interface
