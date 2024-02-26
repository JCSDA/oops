! (C) Copyright 2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module qg_getvalues_interp_interface

use iso_c_binding

use atlas_module, only: atlas_fieldset
use kinds
use oops_variables_mod
use qg_geom_mod
use qg_getvalues_interp_mod

implicit none

private

contains

! ------------------------------------------------------------------------------
subroutine qg_getvalues_interp_c(c_key_geom, c_afieldset, c_vars, c_nlocs, c_locs, c_nvals, c_vals) &
    bind(c, name='qg_getvalues_interp_f90')

integer(c_int), intent(in)     :: c_key_geom
type(c_ptr), value, intent(in) :: c_vars
type(c_ptr), value, intent(in) :: c_afieldset
integer(c_int), intent(in)     :: c_nlocs
real(c_double), intent(in)     :: c_locs(2 * c_nlocs)
integer(c_int), intent(in)     :: c_nvals
real(c_double), intent(inout)  :: c_vals(c_nvals)

type(qg_geom), pointer :: geom
type(atlas_fieldset) :: afieldset
type(oops_variables) :: vars
real(kind_real) :: lats(c_nlocs), lons(c_nlocs)
integer :: ii, jj

call qg_geom_registry%get(c_key_geom,geom)
afieldset = atlas_fieldset(c_afieldset)
vars = oops_variables(c_vars)

ii = 0
do jj = 1, c_nlocs
  ii = ii + 1
  lats(jj) = c_locs(ii)
  ii = ii + 1
  lons(jj) = c_locs(ii)
enddo

call qg_getvalues_interp(geom, afieldset, vars, lats, lons, c_vals)

end subroutine qg_getvalues_interp_c

! ------------------------------------------------------------------------------
subroutine qg_getvalues_interp_ad_c(c_key_geom, c_afieldset, c_vars, c_nlocs, c_locs, c_nvals, c_vals) &
    bind(c, name='qg_getvalues_interp_ad_f90')

integer(c_int), intent(in)     :: c_key_geom
type(c_ptr), value, intent(in) :: c_afieldset
type(c_ptr), value, intent(in) :: c_vars
integer(c_int), intent(in)     :: c_nlocs
real(c_double), intent(in)     :: c_locs(2 * c_nlocs)
integer(c_int), intent(in)     :: c_nvals
real(c_double), intent(in)     :: c_vals(c_nvals)

type(qg_geom), pointer :: geom
type(atlas_fieldset) :: afieldset
type(oops_variables) :: vars
real(kind_real) :: lats(c_nlocs), lons(c_nlocs)
integer :: ii, jj

call qg_geom_registry%get(c_key_geom,geom)
afieldset = atlas_fieldset(c_afieldset)
vars = oops_variables(c_vars)

ii = 0
do jj = 1, c_nlocs
  ii = ii + 1
  lats(jj) = c_locs(ii)
  ii = ii + 1
  lons(jj) = c_locs(ii)
enddo

call qg_getvalues_interp_ad(geom, afieldset, vars, lats, lons, c_vals)

end subroutine qg_getvalues_interp_ad_c

! ------------------------------------------------------------------------------
end module qg_getvalues_interp_interface
