! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Interfaces to be called from C++ for Fortran handling of QG model fields

! ------------------------------------------------------------------------------

subroutine qg_field_create_c(c_key_self, c_key_geom, c_key_vars) bind(c,name='qg_field_create_f90')
use iso_c_binding
use qg_fields
use qg_geom_mod
use qg_vars_mod
implicit none
integer(c_int), intent(inout) :: c_key_self
integer(c_int), intent(in) :: c_key_geom !< Geometry
integer(c_int), intent(in) :: c_key_vars !< List of variables

type(qg_field), pointer :: self
type(qg_geom),  pointer :: geom
type(qg_vars),  pointer :: vars

call qg_geom_registry%get(c_key_geom, geom)
call qg_vars_registry%get(c_key_vars, vars)
call qg_field_registry%init()
call qg_field_registry%add(c_key_self)
call qg_field_registry%get(c_key_self,self)

call create(self, geom, vars)

end subroutine qg_field_create_c

! ------------------------------------------------------------------------------

subroutine qg_field_delete_c(c_key_self) bind(c,name='qg_field_delete_f90')
use iso_c_binding
use qg_fields
implicit none
integer(c_int), intent(inout) :: c_key_self
type(qg_field), pointer :: self

call qg_field_registry%get(c_key_self,self)

call delete(self)

call qg_field_registry%remove(c_key_self)

end subroutine qg_field_delete_c

! ------------------------------------------------------------------------------

subroutine qg_field_zero_c(c_key_self) bind(c,name='qg_field_zero_f90')
use iso_c_binding
use qg_fields
implicit none
integer(c_int), intent(in) :: c_key_self
type(qg_field), pointer :: self

call qg_field_registry%get(c_key_self,self)
call zeros(self)

end subroutine qg_field_zero_c

! ------------------------------------------------------------------------------

subroutine qg_field_random_c(c_key_self) bind(c,name='qg_field_random_f90')
use iso_c_binding
use qg_fields
implicit none
integer(c_int), intent(in) :: c_key_self
type(qg_field), pointer :: self

call qg_field_registry%get(c_key_self,self)
call random(self)

end subroutine qg_field_random_c

! ------------------------------------------------------------------------------

subroutine qg_field_copy_c(c_key_self,c_key_rhs) bind(c,name='qg_field_copy_f90')
use iso_c_binding
use qg_fields
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(qg_field), pointer :: self
type(qg_field), pointer :: rhs
call qg_field_registry%get(c_key_self,self)
call qg_field_registry%get(c_key_rhs,rhs)

call copy(self, rhs)

end subroutine qg_field_copy_c

! ------------------------------------------------------------------------------

subroutine qg_field_self_add_c(c_key_self,c_key_rhs) bind(c,name='qg_field_self_add_f90')
use iso_c_binding
use qg_fields
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(qg_field), pointer :: self
type(qg_field), pointer :: rhs
call qg_field_registry%get(c_key_self,self)
call qg_field_registry%get(c_key_rhs,rhs)

call self_add(self,rhs)

end subroutine qg_field_self_add_c

! ------------------------------------------------------------------------------

subroutine qg_field_self_schur_c(c_key_self,c_key_rhs) bind(c,name='qg_field_self_schur_f90')
use iso_c_binding
use qg_fields
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(qg_field), pointer :: self
type(qg_field), pointer :: rhs
call qg_field_registry%get(c_key_self,self)
call qg_field_registry%get(c_key_rhs,rhs)

call self_schur(self,rhs)

end subroutine qg_field_self_schur_c

! ------------------------------------------------------------------------------

subroutine qg_field_self_sub_c(c_key_self,c_key_rhs) bind(c,name='qg_field_self_sub_f90')
use iso_c_binding
use qg_fields
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs

type(qg_field), pointer :: self
type(qg_field), pointer :: rhs
call qg_field_registry%get(c_key_self,self)
call qg_field_registry%get(c_key_rhs,rhs)

call self_sub(self,rhs)

end subroutine qg_field_self_sub_c

! ------------------------------------------------------------------------------

subroutine qg_field_self_mul_c(c_key_self,c_zz) bind(c,name='qg_field_self_mul_f90')
use iso_c_binding
use qg_fields
use kinds
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: c_zz
type(qg_field), pointer :: self
real(kind=kind_real) :: zz

call qg_field_registry%get(c_key_self,self)
zz = c_zz

call self_mul(self,zz)

end subroutine qg_field_self_mul_c

! ------------------------------------------------------------------------------

subroutine qg_field_axpy_c(c_key_self,c_zz,c_key_rhs) bind(c,name='qg_field_axpy_f90')
use iso_c_binding
use qg_fields
use kinds
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(in) :: c_zz
integer(c_int), intent(in) :: c_key_rhs

type(qg_field), pointer :: self
type(qg_field), pointer :: rhs
real(kind=kind_real) :: zz

call qg_field_registry%get(c_key_self,self)
call qg_field_registry%get(c_key_rhs,rhs)
zz = c_zz

call axpy(self,zz,rhs)

end subroutine qg_field_axpy_c

! ------------------------------------------------------------------------------

subroutine qg_field_dot_prod_c(c_key_fld1,c_key_fld2,c_prod) bind(c,name='qg_field_dot_prod_f90')
use iso_c_binding
use qg_fields
use kinds
implicit none
integer(c_int), intent(in)    :: c_key_fld1, c_key_fld2
real(c_double), intent(inout) :: c_prod
real(kind=kind_real) :: zz
type(qg_field), pointer :: fld1, fld2

call qg_field_registry%get(c_key_fld1,fld1)
call qg_field_registry%get(c_key_fld2,fld2)

call dot_prod(fld1,fld2,zz)

c_prod = zz

end subroutine qg_field_dot_prod_c

! ------------------------------------------------------------------------------

subroutine qg_field_add_incr_c(c_key_self,c_key_rhs) bind(c,name='qg_field_add_incr_f90')
use iso_c_binding
use qg_fields
implicit none
integer(c_int), intent(in) :: c_key_self
integer(c_int), intent(in) :: c_key_rhs
type(qg_field), pointer :: self
type(qg_field), pointer :: rhs

call qg_field_registry%get(c_key_self,self)
call qg_field_registry%get(c_key_rhs,rhs)

call add_incr(self,rhs)

end subroutine qg_field_add_incr_c

! ------------------------------------------------------------------------------

subroutine qg_field_diff_incr_c(c_key_lhs,c_key_x1,c_key_x2) bind(c,name='qg_field_diff_incr_f90')
use iso_c_binding
use qg_fields
implicit none
integer(c_int), intent(in) :: c_key_lhs
integer(c_int), intent(in) :: c_key_x1
integer(c_int), intent(in) :: c_key_x2
type(qg_field), pointer :: lhs
type(qg_field), pointer :: x1
type(qg_field), pointer :: x2

call qg_field_registry%get(c_key_lhs,lhs)
call qg_field_registry%get(c_key_x1,x1)
call qg_field_registry%get(c_key_x2,x2)

call diff_incr(lhs,x1,x2)

end subroutine qg_field_diff_incr_c

! ------------------------------------------------------------------------------

subroutine qg_field_change_resol_c(c_key_fld,c_key_rhs) bind(c,name='qg_field_change_resol_f90')
use iso_c_binding
use qg_fields
implicit none
integer(c_int), intent(in) :: c_key_fld
integer(c_int), intent(in) :: c_key_rhs
type(qg_field), pointer :: fld, rhs

call qg_field_registry%get(c_key_fld,fld)
call qg_field_registry%get(c_key_rhs,rhs)

call change_resol(fld,rhs)

end subroutine qg_field_change_resol_c

! ------------------------------------------------------------------------------
subroutine qg_field_convert_to_c(c_key_fld, c_key_ug) bind (c,name='qg_field_convert_to_f90')
use iso_c_binding
use qg_fields
use unstructured_grid_mod
implicit none
integer(c_int), intent(in) :: c_key_fld
integer(c_int), intent(in) :: c_key_ug
type(qg_field), pointer :: fld
type(unstructured_grid), pointer :: ug

call qg_field_registry%get(c_key_fld,fld)
call unstructured_grid_registry%get(c_key_ug,ug)

call convert_to_ug(fld, ug)

end subroutine qg_field_convert_to_c
! ------------------------------------------------------------------------------
subroutine qg_field_convert_from_c(c_key_fld, c_key_ug) bind (c,name='qg_field_convert_from_f90')
use iso_c_binding
use qg_fields
use unstructured_grid_mod
implicit none
integer(c_int), intent(in) :: c_key_fld
integer(c_int), intent(in) :: c_key_ug
type(qg_field), pointer :: fld
type(unstructured_grid), pointer :: ug

call qg_field_registry%get(c_key_fld,fld)
call unstructured_grid_registry%get(c_key_ug,ug)

call convert_from_ug(fld, ug)

end subroutine qg_field_convert_from_c
! ------------------------------------------------------------------------------

subroutine qg_field_read_file_c(c_key_fld, c_conf, c_dt) bind(c,name='qg_field_read_file_f90')
use iso_c_binding
use qg_fields
use datetime_mod

implicit none
integer(c_int), intent(in) :: c_key_fld  !< Fields
type(c_ptr), intent(in)    :: c_conf !< Configuration
type(c_ptr), intent(inout) :: c_dt   !< DateTime

type(qg_field), pointer :: fld
type(datetime) :: fdate

call qg_field_registry%get(c_key_fld,fld)
call c_f_datetime(c_dt, fdate)
call read_file(fld, c_conf, fdate)

end subroutine qg_field_read_file_c

! ------------------------------------------------------------------------------

subroutine qg_field_write_file_c(c_key_fld, c_conf, c_dt) bind(c,name='qg_field_write_file_f90')
use iso_c_binding
use qg_fields
use datetime_mod

implicit none
integer(c_int), intent(in) :: c_key_fld  !< Fields
type(c_ptr), intent(in) :: c_conf !< Configuration
type(c_ptr), intent(in) :: c_dt   !< DateTime

type(qg_field), pointer :: fld
type(datetime) :: fdate

call qg_field_registry%get(c_key_fld,fld)
call c_f_datetime(c_dt, fdate)
call write_file(fld, c_conf, fdate)

end subroutine qg_field_write_file_c

! ------------------------------------------------------------------------------

subroutine qg_field_gpnorm_c(c_key_fld, kf, pstat) bind(c,name='qg_field_gpnorm_f90')
use iso_c_binding
use qg_fields
use kinds
implicit none
integer(c_int), intent(in) :: c_key_fld
integer(c_int), intent(in) :: kf
real(c_double), intent(inout) :: pstat(3*kf)

type(qg_field), pointer :: fld
real(kind=kind_real) :: zstat(3, kf)
integer :: jj, js, jf

call qg_field_registry%get(c_key_fld,fld)

call gpnorm(fld, kf, zstat)
jj=0
do jf = 1, kf
  do js = 1, 3
    jj=jj+1
    pstat(jj) = zstat(js,jf)
  enddo
enddo

end subroutine qg_field_gpnorm_c

! ------------------------------------------------------------------------------

subroutine qg_field_rms_c(c_key_fld, prms) bind(c,name='qg_field_rms_f90')
use iso_c_binding
use qg_fields
use kinds
implicit none
integer(c_int), intent(in) :: c_key_fld
real(c_double), intent(inout) :: prms

type(qg_field), pointer :: fld
real(kind=kind_real) :: zz

call qg_field_registry%get(c_key_fld,fld)

call fldrms(fld, zz)

prms = zz

end subroutine qg_field_rms_c

! ------------------------------------------------------------------------------

subroutine qg_field_interp_tl_c(c_key_fld,c_key_loc,c_key_gom) bind(c,name='qg_field_interp_tl_f90')
use iso_c_binding
use qg_fields
use qg_locs_mod
use qg_goms_mod
implicit none
integer(c_int), intent(in) :: c_key_fld
integer(c_int), intent(in) :: c_key_loc
integer(c_int), intent(in) :: c_key_gom
type(qg_field), pointer :: fld
type(qg_locs),  pointer :: locs
type(qg_goms),  pointer :: gom

call qg_field_registry%get(c_key_fld,fld)
call qg_locs_registry%get(c_key_loc,locs)
call qg_goms_registry%get(c_key_gom,gom)

call interp_tl(fld, locs, gom)

end subroutine qg_field_interp_tl_c

! ------------------------------------------------------------------------------

subroutine qg_field_interp_ad_c(c_key_fld,c_key_loc,c_key_gom) bind(c,name='qg_field_interp_ad_f90')
use iso_c_binding
use qg_fields
use qg_locs_mod
use qg_goms_mod
implicit none
integer(c_int), intent(in) :: c_key_fld
integer(c_int), intent(in) :: c_key_loc
integer(c_int), intent(in) :: c_key_gom
type(qg_field), pointer :: fld
type(qg_locs),  pointer :: locs
type(qg_goms),  pointer :: gom

call qg_field_registry%get(c_key_fld,fld)
call qg_locs_registry%get(c_key_loc,locs)
call qg_goms_registry%get(c_key_gom,gom)

call interp_ad(fld, locs, gom)

end subroutine qg_field_interp_ad_c

! ------------------------------------------------------------------------------

subroutine qg_fieldnum_c(c_key_fld, nx, ny, nf, nb) bind(c,name='qg_field_sizes_f90')
use iso_c_binding
use qg_fields
implicit none
integer(c_int), intent(in) :: c_key_fld
integer(kind=c_int), intent(inout) :: nx, ny, nf, nb
type(qg_field), pointer :: fld

call qg_field_registry%get(c_key_fld,fld)

nx = fld%nx
ny = fld%ny
nf = fld%nf
nb =0
if (fld%lbc) nb = 2

end subroutine qg_fieldnum_c

! ------------------------------------------------------------------------------
