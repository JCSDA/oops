! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

module unstructured_grid_interface

use atlas_module
use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use unstructured_grid_mod

implicit none

private
!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------
!> Create unstructured grid
subroutine create_ug_c(key,colocated,nts) bind(c,name='create_ug_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: key    !< Unstructured grid
integer(c_int),intent(in) :: colocated !< Colocated flag
integer(c_int),intent(in) :: nts       !< Number of timeslots
 
! Local variables
type(unstructured_grid),pointer :: self

! Interface
call unstructured_grid_registry%init()
call unstructured_grid_registry%add(key)
call unstructured_grid_registry%get(key,self)

! Call Fortran
call create_ug(self,colocated,nts)

end subroutine create_ug_c
! ------------------------------------------------------------------------------
!> Delete unstructured grid
subroutine delete_ug_c(key) bind(c,name='delete_ug_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: key !< Unstructured grid

! Local variables
type(unstructured_grid),pointer :: self

! Interface
call unstructured_grid_registry%get(key,self)

! Call Fortran
call delete_ug(self)

! Clear interface
call unstructured_grid_registry%remove(key)

end subroutine delete_ug_c
! ------------------------------------------------------------------------------
!> Get number of grids
subroutine ug_get_ngrid_c(key,ngrid) bind(c,name='ug_get_ngrid_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key    !< Unstructured grid
integer(c_int),intent(out) :: ngrid !< Number of grids

! Local variables
type(unstructured_grid),pointer :: self

! Interface
call unstructured_grid_registry%get(key,self)

! Call Fortran
ngrid = self%ngrid

end subroutine ug_get_ngrid_c
! ------------------------------------------------------------------------------
!> Get dimensions for a given grid
subroutine ug_get_dims_c(key,c_igrid,nl,nv,nts) bind(c,name='ug_get_dims_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key     !< Unstructured grid
integer(c_int),intent(in) :: c_igrid !< Grid index
integer(c_int),intent(out) :: nl     !< Number of levels
integer(c_int),intent(out) :: nv     !< Number of variables
integer(c_int),intent(out) :: nts    !< Number of timeslots

! Local variables
type(unstructured_grid),pointer :: self
integer :: igrid

! Interface
call unstructured_grid_registry%get(key,self)
igrid = c_igrid+1

! Call Fortran
nl = self%grid(igrid)%nl0
nv = self%grid(igrid)%nv
nts = self%nts

end subroutine ug_get_dims_c
! ------------------------------------------------------------------------------
!> Create ATLAS grid configuration from unstructured grid
subroutine ug_create_atlas_grid_conf_c(key,c_grid) bind(c,name='ug_create_atlas_grid_conf_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key !< Unstructured grid
type(c_ptr),intent(in) :: c_grid !< ATLAS grid configuration

! Local variables
type(unstructured_grid),pointer :: self
type(fckit_configuration) :: f_grid

! Interface
call unstructured_grid_registry%get(key,self)
f_grid = fckit_configuration(c_grid)

! Call Fortran
call ug_create_atlas_grid_conf(self,f_grid)

end subroutine ug_create_atlas_grid_conf_c
! ------------------------------------------------------------------------------
!> Set ATLAS function space pointer
subroutine ug_set_atlas_functionspace_pointer_c(key,c_afunctionspace) bind(c,name='ug_set_atlas_functionspace_pointer_f90')

! Passed variables
integer(c_int),intent(in) :: key                 !< Unstructured grid
type(c_ptr),intent(in),value :: c_afunctionspace !< ATLAS function space pointer

! Local variables
type(unstructured_grid),pointer :: self

! Interface
call unstructured_grid_registry%get(key,self)
self%afunctionspace = atlas_functionspace_nodecolumns(c_afunctionspace)

end subroutine ug_set_atlas_functionspace_pointer_c
! ------------------------------------------------------------------------------
!> Fill ATLAS fieldset from unstructured grid
subroutine ug_fill_atlas_fieldset_c(key,c_afieldset) bind(c,name='ug_fill_atlas_fieldset_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key            !< Unstructured grid
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! Local variables
type(unstructured_grid),pointer :: self
type(atlas_fieldset) :: afieldset

! Interface
call unstructured_grid_registry%get(key,self)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call ug_fill_atlas_fieldset(self,afieldset)

end subroutine ug_fill_atlas_fieldset_c
! ------------------------------------------------------------------------------
!> Set ATLAS fieldset pointer
subroutine ug_set_atlas_fieldset_pointer_c(key,c_afieldset) bind(c,name='ug_set_atlas_fieldset_pointer_f90')

! Passed variables
integer(c_int),intent(in) :: key            !< Unstructured grid
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS function space pointer

! Local variables
type(unstructured_grid),pointer :: self

! Interface
call unstructured_grid_registry%get(key,self)
self%afieldset = atlas_fieldset(c_afieldset)

end subroutine ug_set_atlas_fieldset_pointer_c
! ------------------------------------------------------------------------------
!> Set ATLAS fieldset
subroutine ug_set_atlas_c(key,c_afieldset) bind(c,name='ug_set_atlas_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key            !< Unstructured grid
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! Local variables
type(unstructured_grid),pointer :: self
type(atlas_fieldset) :: afieldset

! Interface
call unstructured_grid_registry%get(key,self)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call ug_set_atlas(self,afieldset)

end subroutine ug_set_atlas_c
! ------------------------------------------------------------------------------
!> Unstructured grid to ATLAS fieldset
subroutine ug_to_atlas_c(key,c_afieldset) bind(c,name='ug_to_atlas_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key                 !< Unstructured grid
type(c_ptr),intent(in),value :: c_afieldset      !< ATLAS fieldset pointer

! Local variables
type(unstructured_grid),pointer :: self
type(atlas_fieldset) :: afieldset

! Interface
call unstructured_grid_registry%get(key,self)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call ug_to_atlas(self,afieldset)

end subroutine ug_to_atlas_c
! ------------------------------------------------------------------------------
!> Unstructured grid from ATLAS fieldset
subroutine ug_from_atlas_c(key,c_afieldset) bind(c,name='ug_from_atlas_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: key            !< Unstructured grid
type(c_ptr),intent(in),value :: c_afieldset !< ATLAS fieldset pointer

! Local variables
type(unstructured_grid),pointer :: self
type(atlas_fieldset) :: afieldset

! Interface
call unstructured_grid_registry%get(key,self)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call ug_from_atlas(self,afieldset)

end subroutine ug_from_atlas_c
! ------------------------------------------------------------------------------
end module unstructured_grid_interface
