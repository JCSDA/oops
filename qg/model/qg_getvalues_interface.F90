! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module qg_getvalues_interface

use datetime_mod
use iso_c_binding
use kinds
use qg_fields_mod
use qg_getvalues_mod
use qg_geom_mod
use qg_gom_mod
use qg_locs_mod

implicit none

private
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Create GetValues
subroutine qg_getvalues_create_c(c_key_self,c_key_geom,c_key_locs) bind(c,name='qg_getvalues_create_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self       !< GetValues
integer(c_int),intent(in) :: c_key_geom          !< Geometry
integer(c_int),intent(in) :: c_key_locs          !< Locations

! Local variables
type(qg_getvalues),pointer :: self
type(qg_geom),pointer :: geom
type(qg_locs),pointer :: locs

! Interface
call qg_getvalues_registry%init()
call qg_getvalues_registry%add(c_key_self)
call qg_getvalues_registry%get(c_key_self,self)
call qg_geom_registry%get(c_key_geom,geom)
call qg_locs_registry%get(c_key_locs,locs)
! Call Fortran
call qg_getvalues_create(self,geom,locs)

end subroutine qg_getvalues_create_c
! ------------------------------------------------------------------------------
!> Delete GetValues
subroutine qg_getvalues_delete_c(c_key_self) bind(c,name='qg_getvalues_delete_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Fields

! Local variables
type(qg_getvalues),pointer :: self

! Interface
call qg_getvalues_registry%get(c_key_self,self)

! Call Fortran
call qg_getvalues_delete(self)

! Clear interface
call qg_getvalues_registry%remove(c_key_self)

end subroutine qg_getvalues_delete_c
! ------------------------------------------------------------------------------
!> Interpolation from fields
subroutine qg_getvalues_interp_c(c_key_self,c_key_fld,c_t1,c_t2,c_key_gom) bind(c,name='qg_getvalues_interp_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self          !< GetValues
integer(c_int),intent(in) :: c_key_fld           !< Fields
type(c_ptr),value,intent(in) :: c_t1, c_t2       !< times
integer(c_int),intent(in) :: c_key_gom           !< Interpolated values

! Local variables
type(qg_getvalues),pointer :: self
type(qg_fields),pointer :: fld
type(qg_gom), pointer :: gom
type(datetime) :: t1, t2

! Interface
call qg_getvalues_registry%get(c_key_self,self)
call qg_fields_registry%get(c_key_fld,fld)
call qg_gom_registry%get(c_key_gom,gom)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)
! Call Fortran
call qg_getvalues_interp(self,fld,t1,t2,gom)

end subroutine qg_getvalues_interp_c
! ------------------------------------------------------------------------------
!> Interpolation from fields, tangent linear
subroutine qg_getvalues_interp_tl_c(c_key_self, c_key_fld,c_t1,c_t2,c_key_gom) bind(c,name='qg_getvalues_interp_tl_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self          !< GetValues
integer(c_int),intent(in) :: c_key_fld           !< Fields
type(c_ptr),value,intent(in) :: c_t1, c_t2       !< times
integer(c_int),intent(in) :: c_key_gom           !< Interpolated values

! Local variables
type(qg_getvalues),pointer :: self
type(qg_fields),pointer :: fld
type(datetime) :: t1, t2
type(qg_gom), pointer :: gom

! Interface
call qg_getvalues_registry%get(c_key_self, self)
call qg_fields_registry%get(c_key_fld,fld)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)
call qg_gom_registry%get(c_key_gom,gom)

! Call Fortran
call qg_getvalues_interp_tl(self,fld,t1,t2,gom)

end subroutine qg_getvalues_interp_tl_c
! ------------------------------------------------------------------------------
!> Interpolation from fields, adjoint
subroutine qg_getvalues_interp_ad_c(c_key_self,c_key_fld,c_t1,c_t2,c_key_gom) bind(c,name='qg_getvalues_interp_ad_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self          !< GetValues
integer(c_int),intent(in) :: c_key_fld           !< Fields
type(c_ptr),value,intent(in) :: c_t1, c_t2       !< times
integer(c_int),intent(in) :: c_key_gom           !< Interpolated values

! Local variables
type(qg_getvalues),pointer :: self
type(qg_fields),pointer :: fld
type(datetime) :: t1, t2
type(qg_gom), pointer :: gom

! Interface
call qg_getvalues_registry%get(c_key_self, self)
call qg_fields_registry%get(c_key_fld,fld)
call qg_gom_registry%get(c_key_gom,gom)
call c_f_datetime(c_t1, t1)
call c_f_datetime(c_t2, t2)

! Call Fortran
call qg_getvalues_interp_ad(self,fld,t1,t2,gom)

end subroutine qg_getvalues_interp_ad_c
! ------------------------------------------------------------------------------
end module qg_getvalues_interface
