! (C) Copyright 2023 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module qg_getvalues_interp_mod

use iso_c_binding

use atlas_module, only: atlas_field, atlas_fieldset
use kinds
use oops_variables_mod
use qg_geom_mod
use qg_interp_mod

implicit none

private
public :: qg_getvalues_interp, qg_getvalues_interp_ad

contains

! ------------------------------------------------------------------------------
subroutine qg_getvalues_interp(geom, afieldset, vars, lats, lons, vals)

type(qg_geom), intent(in)        :: geom
type(atlas_fieldset), intent(in) :: afieldset
type(oops_variables), intent(in) :: vars
real(kind_real), intent(in)      :: lats(:)
real(kind_real), intent(in)      :: lons(:)
real(c_double), intent(inout)    :: vals(:)

integer :: nx, ny, nlev, nlocs
integer :: jvar, jloc, jlev, ii
character(len=1024) :: fname

type(atlas_field) :: afield
real(kind_real), pointer :: ptr(:,:)
real(kind_real), allocatable :: array(:,:,:)

nlocs = size(lats)
nx = geom%nx
ny = geom%ny
nlev = geom%nz
allocate(array(nx, ny, nlev))

ii = 0
do jvar=1, vars%nvars()
  fname = vars%variable(jvar)

  select case (trim(fname))
  case ('x', 'q', 'u', 'v')
    afield = afieldset%field(fname)
    call afield%data(ptr)
    do jlev=1,nlev
      array(:,:,jlev) = reshape(ptr(jlev,:), (/nx, ny/))
    enddo
    do jloc=1,nlocs
      call qg_interp_bilinear(geom, lons(jloc), lats(jloc), array(:,:,:), vals(ii+1:ii+nlev))
      ii = ii + nlev
    enddo
  case ('z')
    do jloc=1,nlocs
      vals(ii+1:ii+nlev) = geom%z(:)
      ii = ii + nlev
    enddo
  case default
    call abor1_ftn('qg_getvalues_interp: wrong input variable')
  endselect
enddo

if (size(vals) /= ii) call abor1_ftn('qg_getvalues_interp: error size')

end subroutine qg_getvalues_interp

! ------------------------------------------------------------------------------
subroutine qg_getvalues_interp_ad(geom, afieldset, vars, lats, lons, vals)

type(qg_geom), intent(in)           :: geom
type(atlas_fieldset), intent(inout) :: afieldset
type(oops_variables), intent(in)    :: vars
real(kind_real), intent(in)         :: lats(:)
real(kind_real), intent(in)         :: lons(:)
real(c_double), intent(in)          :: vals(:)

integer :: nx, ny, nlev, nlocs
integer :: jvar, jloc, jlev, ii
character(len=1024) :: fname

type(atlas_field) :: afield
real(kind_real), pointer :: ptr(:,:)
real(kind_real), allocatable :: array(:,:,:)

nlocs = size(lats)
nx = geom%nx
ny = geom%ny
nlev = geom%nz
allocate(array(nx, ny, nlev))

ii = 0
do jvar=1,vars%nvars()
  fname = vars%variable(jvar)

  select case (trim(fname))
  case ('x', 'q', 'u', 'v')
    afield = afieldset%field(fname)
    call afield%data(ptr)
    do jlev=1,nlev
      array(:,:,jlev) = reshape(ptr(jlev,:), (/nx, ny/))
    enddo
    do jloc=1,nlocs
      call qg_interp_bilinear_ad(geom, lons(jloc), lats(jloc), vals(ii+1:ii+nlev), array(:,:,:))
      ii = ii + nlev
    enddo
    do jlev=1,nlev
      ptr(jlev,:) = reshape(array(:,:,jlev), (/nx*ny/))
    enddo
  case ('z')
    ! do nothing
    ii = ii + nlocs * nlev
  case default
    call abor1_ftn('qg_getvalues_interp_ad: wrong input variable')
  endselect
enddo

if (size(vals) /= ii) call abor1_ftn('qg_getvalues_interp_ad: error size')

end subroutine qg_getvalues_interp_ad

! ------------------------------------------------------------------------------
end module qg_getvalues_interp_mod
