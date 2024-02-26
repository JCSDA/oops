! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_fields_interface

use atlas_module, only: atlas_fieldset, atlas_field
use datetime_mod
use fckit_configuration_module, only: fckit_configuration
use iso_c_binding
use kinds
use oops_variables_mod
use qg_fields_mod
use qg_geom_mod
use qg_geom_iter_mod
use qg_locs_mod

implicit none

private
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Create fields from geometry and variables
subroutine qg_fields_create_c(c_key_self,c_key_geom,c_vars,c_lbc) bind(c,name='qg_fields_create_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self       !< Fields
integer(c_int),intent(in) :: c_key_geom          !< Geometry
type(c_ptr),value,intent(in) :: c_vars           !< List of variables
logical(c_bool),intent(in) :: c_lbc              !< Boundaries flag

! Local variables
type(qg_fields),pointer :: self
type(qg_geom),pointer :: geom
type(oops_variables) :: vars
logical :: lbc

! Interface
call qg_fields_registry%init()
call qg_fields_registry%add(c_key_self)
call qg_fields_registry%get(c_key_self,self)
call qg_geom_registry%get(c_key_geom,geom)
vars = oops_variables(c_vars)
lbc = c_lbc

! Call Fortran
call qg_fields_create(self,geom,vars,lbc)

end subroutine qg_fields_create_c
! ------------------------------------------------------------------------------
!> Create fields from another one
subroutine qg_fields_create_from_other_c(c_key_self,c_key_other,c_key_geom) bind(c,name='qg_fields_create_from_other_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self  !< Fields
integer(c_int),intent(in)    :: c_key_other !< Other fields
integer(c_int),intent(in) :: c_key_geom     !< Geometry

! Local variables
type(qg_fields),pointer :: self
type(qg_fields),pointer :: other
type(qg_geom),pointer :: geom

! Interface
call qg_fields_registry%get(c_key_other,other)
call qg_fields_registry%init()
call qg_fields_registry%add(c_key_self)
call qg_fields_registry%get(c_key_self,self)
call qg_geom_registry%get(c_key_geom,geom)

! Call Fortran
call qg_fields_create_from_other(self,other,geom)

end subroutine qg_fields_create_from_other_c
! ------------------------------------------------------------------------------
!> Delete fields
subroutine qg_fields_delete_c(c_key_self) bind(c,name='qg_fields_delete_f90')

implicit none

! Passed variables
integer(c_int),intent(inout) :: c_key_self !< Fields

! Local variables
type(qg_fields),pointer :: self

! Interface
call qg_fields_registry%get(c_key_self,self)

! Call Fortran
call qg_fields_delete(self)

! Clear interface
call qg_fields_registry%remove(c_key_self)

end subroutine qg_fields_delete_c
! ------------------------------------------------------------------------------
!> Set fields to zero
subroutine qg_fields_zero_c(c_key_self) bind(c,name='qg_fields_zero_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields

! Local variables
type(qg_fields),pointer :: self

! Interface
call qg_fields_registry%get(c_key_self,self)

! Call Fortran
call qg_fields_zero(self)

end subroutine qg_fields_zero_c
! ------------------------------------------------------------------------------
!> Set fields to ones
subroutine qg_fields_ones_c(c_key_self) bind(c,name='qg_fields_ones_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields

! Local variables
type(qg_fields),pointer :: self

! Interface
call qg_fields_registry%get(c_key_self,self)

! Call Fortran
call qg_fields_ones(self)

end subroutine qg_fields_ones_c
! ------------------------------------------------------------------------------
!> Set fields to Diracs
subroutine qg_fields_dirac_c(c_key_self,c_conf) bind(c,name='qg_fields_dirac_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields
type(c_ptr),value,intent(in) :: c_conf  !< Configuration

! Local variables
type(fckit_configuration) :: f_conf
type(qg_fields),pointer :: self

! Interface
f_conf = fckit_configuration(c_conf)
call qg_fields_registry%get(c_key_self,self)

! Call Fortran
call qg_fields_dirac(self,f_conf)

end subroutine qg_fields_dirac_c
! ------------------------------------------------------------------------------
!> Generate random fields
subroutine qg_fields_random_c(c_key_self) bind(c,name='qg_fields_random_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields

! Local variables
type(qg_fields),pointer :: self

! Interface
call qg_fields_registry%get(c_key_self,self)

! Call Fortran
call qg_fields_random(self)

end subroutine qg_fields_random_c
! ------------------------------------------------------------------------------
!> Copy fields
subroutine qg_fields_copy_c(c_key_self,c_key_other) bind(c,name='qg_fields_copy_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self  !< Fields
integer(c_int),intent(in) :: c_key_other !< Other fields

! Local variables
type(qg_fields),pointer :: self
type(qg_fields),pointer :: other

! Interface
call qg_fields_registry%get(c_key_self,self)
call qg_fields_registry%get(c_key_other,other)

! Call Fortran
call qg_fields_copy(self,other)

end subroutine qg_fields_copy_c
! ------------------------------------------------------------------------------
!> Add fields
subroutine qg_fields_self_add_c(c_key_self,c_key_rhs) bind(c,name='qg_fields_self_add_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields
integer(c_int),intent(in) :: c_key_rhs  !< Right-hand side

! Local variables
type(qg_fields),pointer :: self
type(qg_fields),pointer :: rhs

! Interface
call qg_fields_registry%get(c_key_self,self)
call qg_fields_registry%get(c_key_rhs,rhs)

! Call Fortran
call qg_fields_self_add(self,rhs)

end subroutine qg_fields_self_add_c
! ------------------------------------------------------------------------------
!> Subtract fields
subroutine qg_fields_self_sub_c(c_key_self,c_key_rhs) bind(c,name='qg_fields_self_sub_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields
integer(c_int),intent(in) :: c_key_rhs  !< Right-hand side

! Local variables
type(qg_fields),pointer :: self
type(qg_fields),pointer :: rhs

! Interface
call qg_fields_registry%get(c_key_self,self)
call qg_fields_registry%get(c_key_rhs,rhs)

! Call Fortran
call qg_fields_self_sub(self,rhs)

end subroutine qg_fields_self_sub_c
! ------------------------------------------------------------------------------
!> Multiply fields by a scalar
subroutine qg_fields_self_mul_c(c_key_self,c_zz) bind(c,name='qg_fields_self_mul_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields
real(c_double),intent(in) :: c_zz       !< Multiplier

! Local variables
type(qg_fields),pointer :: self

! Interface
call qg_fields_registry%get(c_key_self,self)

! Call Fortran
call qg_fields_self_mul(self,c_zz)

end subroutine qg_fields_self_mul_c
! ------------------------------------------------------------------------------
!> Apply axpy operator to fields
subroutine qg_fields_axpy_c(c_key_self,c_zz,c_key_rhs) bind(c,name='qg_fields_axpy_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields
real(c_double),intent(in) :: c_zz       !< Multiplier
integer(c_int),intent(in) :: c_key_rhs  !< Right-hand side

! Local variables
type(qg_fields),pointer :: self
type(qg_fields),pointer :: rhs

! Interface
call qg_fields_registry%get(c_key_self,self)
call qg_fields_registry%get(c_key_rhs,rhs)

! Call Fortran
call qg_fields_axpy(self,c_zz,rhs)

end subroutine qg_fields_axpy_c
! ------------------------------------------------------------------------------
!> Schur product of fields
subroutine qg_fields_self_schur_c(c_key_self,c_key_rhs) bind(c,name='qg_fields_self_schur_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields
integer(c_int),intent(in) :: c_key_rhs  !< Right-hand side

! Local variables
type(qg_fields),pointer :: self
type(qg_fields),pointer :: rhs

! Interface
call qg_fields_registry%get(c_key_self,self)
call qg_fields_registry%get(c_key_rhs,rhs)

! Call Fortran
call qg_fields_self_schur(self,rhs)

end subroutine qg_fields_self_schur_c
! ------------------------------------------------------------------------------
!> Compute dot product for fields
subroutine qg_fields_dot_prod_c(c_key_fld1,c_key_fld2,c_prod) bind(c,name='qg_fields_dot_prod_f90')

implicit none

! Passed variables
integer(c_int),intent(in)    :: c_key_fld1 !< First fields
integer(c_int),intent(in)    :: c_key_fld2 !< Second fields
real(c_double),intent(inout) :: c_prod     !< Dot product

! Local variables
type(qg_fields),pointer :: fld1,fld2

! Interface
call qg_fields_registry%get(c_key_fld1,fld1)
call qg_fields_registry%get(c_key_fld2,fld2)

! Call Fortran
call qg_fields_dot_prod(fld1,fld2,c_prod)

end subroutine qg_fields_dot_prod_c
! ------------------------------------------------------------------------------
!> Add increment to fields
subroutine qg_fields_add_incr_c(c_key_self,c_key_rhs) bind(c,name='qg_fields_add_incr_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self !< Fields
integer(c_int),intent(in) :: c_key_rhs  !< Right-hand side

! Local variables
type(qg_fields),pointer :: self
type(qg_fields),pointer :: rhs

! Interface
call qg_fields_registry%get(c_key_self,self)
call qg_fields_registry%get(c_key_rhs,rhs)

! Call Fortran
call qg_fields_add_incr(self,rhs)

end subroutine qg_fields_add_incr_c
! ------------------------------------------------------------------------------
!> Compute increment from the difference of two fields
subroutine qg_fields_diff_incr_c(c_key_lhs,c_key_fld1,c_key_fld2) bind(c,name='qg_fields_diff_incr_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_lhs  !< Left-hand side
integer(c_int),intent(in) :: c_key_fld1 !< First fields
integer(c_int),intent(in) :: c_key_fld2 !< Second fields

! Local variables
type(qg_fields),pointer :: lhs
type(qg_fields),pointer :: fld1
type(qg_fields),pointer :: fld2

! Interface
call qg_fields_registry%get(c_key_lhs,lhs)
call qg_fields_registry%get(c_key_fld1,fld1)
call qg_fields_registry%get(c_key_fld2,fld2)

! Call Fortran
call qg_fields_diff_incr(lhs,fld1,fld2)

end subroutine qg_fields_diff_incr_c
! ------------------------------------------------------------------------------
!> Change fields resolution
subroutine qg_fields_change_resol_c(c_key_fld,c_key_rhs) bind(c,name='qg_fields_change_resol_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld !< Fields
integer(c_int),intent(in) :: c_key_rhs !< Right-hand side

! Local variables
type(qg_fields),pointer :: fld,rhs

! Interface
call qg_fields_registry%get(c_key_fld,fld)
call qg_fields_registry%get(c_key_rhs,rhs)

! Call Fortran
call qg_fields_change_resol(fld,rhs)

end subroutine qg_fields_change_resol_c
! ------------------------------------------------------------------------------
!> Change fields resolution (adjoint)
subroutine qg_fields_change_resol_ad_c(c_key_fld,c_key_rhs) bind(c,name='qg_fields_change_resol_ad_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld !< Fields
integer(c_int),intent(in) :: c_key_rhs !< Right-hand side

! Local variables
type(qg_fields),pointer :: fld,rhs

! Interface
call qg_fields_registry%get(c_key_fld,fld)
call qg_fields_registry%get(c_key_rhs,rhs)

! Call Fortran
call qg_fields_change_resol_ad(fld,rhs)

end subroutine qg_fields_change_resol_ad_c
! ------------------------------------------------------------------------------
!> Read fields from file
subroutine qg_fields_read_file_c(c_key_fld,c_conf,c_dt) bind(c,name='qg_fields_read_file_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld  !< Fields
type(c_ptr),value,intent(in) :: c_conf  !< Configuration
type(c_ptr),value,intent(in) :: c_dt    !< Date and time

! Local variables
type(fckit_configuration) :: f_conf
type(qg_fields),pointer :: fld
type(datetime) :: fdate

! Interface
f_conf = fckit_configuration(c_conf)
call qg_fields_registry%get(c_key_fld,fld)
call c_f_datetime(c_dt,fdate)

! Call Fortran
call qg_fields_read_file(fld,f_conf,fdate)

end subroutine qg_fields_read_file_c
! ------------------------------------------------------------------------------
!> Write fields to file
subroutine qg_fields_write_file_c(c_key_fld,c_conf,c_dt) bind(c,name='qg_fields_write_file_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld !< Fields
type(c_ptr),value,intent(in) :: c_conf !< Configuration
type(c_ptr),value,intent(in) :: c_dt   !< Date and time

! Local variables
type(fckit_configuration) :: f_conf
type(qg_fields),pointer :: fld
type(datetime) :: fdate

! Interface
f_conf = fckit_configuration(c_conf)
call qg_fields_registry%get(c_key_fld,fld)
call c_f_datetime(c_dt,fdate)

! Call Fortran
call qg_fields_write_file(fld,f_conf,fdate)

end subroutine qg_fields_write_file_c
! ------------------------------------------------------------------------------
!> Analytic initialization of fields
subroutine qg_fields_analytic_init_c(c_key_fld,c_conf,c_dt) bind(c,name='qg_fields_analytic_init_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld  !< Fields
type(c_ptr),value,intent(in) :: c_conf  !< Configuration
type(c_ptr),value,intent(in) :: c_dt !< Date and time

! Local variables
type(fckit_configuration) :: f_conf
type(qg_fields),pointer :: fld
type(datetime) :: fdate

! Interface
f_conf = fckit_configuration(c_conf)
call qg_fields_registry%get(c_key_fld,fld)
call c_f_datetime(c_dt,fdate)

! Call Fortran
 call qg_fields_analytic_init(fld,f_conf,fdate)

end subroutine qg_fields_analytic_init_c
! ------------------------------------------------------------------------------
!> Fields statistics
subroutine qg_fields_gpnorm_c(c_key_fld,vpresent,vmin,vmax,vrms) bind(c,name='qg_fields_gpnorm_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld      !< Fields
integer(c_int),intent(inout) :: vpresent(6) !< Variables presence flag
real(c_double),intent(inout) :: vmin(6)     !< Variables minimum
real(c_double),intent(inout) :: vmax(6)     !< Variables maximum
real(c_double),intent(inout) :: vrms(6)     !< Variables RMS

! Local variables
type(qg_fields),pointer :: fld

! Interface
call qg_fields_registry%get(c_key_fld,fld)

! Call Fortran
call qg_fields_gpnorm(fld,vpresent,vmin,vmax,vrms)

end subroutine qg_fields_gpnorm_c
! ------------------------------------------------------------------------------
!> Fields RMS
subroutine qg_fields_rms_c(c_key_fld,prms) bind(c,name='qg_fields_rms_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld !< Fields
real(c_double),intent(inout) :: prms

! Local variables
type(qg_fields),pointer :: fld

! Interface
call qg_fields_registry%get(c_key_fld,fld)

! Call Fortran
call qg_fields_rms(fld,prms)

end subroutine qg_fields_rms_c
! ------------------------------------------------------------------------------
!> Get fields geometry
subroutine qg_fields_sizes_c(c_key_fld,c_nx,c_ny,c_nz) bind(c,name='qg_fields_sizes_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld !< Fields
integer(c_int),intent(inout) :: c_nx   !< X size
integer(c_int),intent(inout) :: c_ny   !< Y size
integer(c_int),intent(inout) :: c_nz   !< Z size

! Local variables
type(qg_fields),pointer :: fld

! Interface
call qg_fields_registry%get(c_key_fld,fld)

! Call Fortran
call qg_fields_sizes(fld,c_nx,c_ny,c_nz)

end subroutine qg_fields_sizes_c
! ------------------------------------------------------------------------------
!> Get fields geometry
subroutine qg_fields_lbc_c(c_key_fld,c_lbc) bind(c,name='qg_fields_lbc_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld !< Fields
integer(c_int),intent(inout) :: c_lbc  !< LBC presence

! Local variables
type(qg_fields),pointer :: fld

! Interface
call qg_fields_registry%get(c_key_fld,fld)

! Call Fortran
call qg_fields_lbc(fld,c_lbc)

end subroutine qg_fields_lbc_c
! ------------------------------------------------------------------------------
!> Convert fields to FieldSet
subroutine qg_fields_to_fieldset_c(c_key_fld,c_afieldset) bind (c,name='qg_fields_to_fieldset_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld      !< Fields
type(c_ptr),intent(in),value :: c_afieldset !< FieldSet pointer

! Local variables
type(qg_fields),pointer :: fld
type(atlas_fieldset) :: afieldset

! Interface
call qg_fields_registry%get(c_key_fld,fld)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call qg_fields_to_fieldset(fld,afieldset)

end subroutine qg_fields_to_fieldset_c
! ------------------------------------------------------------------------------
!> Convert Fieldset to fields
subroutine qg_fields_from_fieldset_c(c_key_fld,c_afieldset) bind (c,name='qg_fields_from_fieldset_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld      !< Fields
type(c_ptr),intent(in),value :: c_afieldset !< FieldSet pointer

! Local variables
type(qg_fields),pointer :: fld
type(atlas_fieldset) :: afieldset

! Interface
call qg_fields_registry%get(c_key_fld,fld)
afieldset = atlas_fieldset(c_afieldset)

! Call Fortran
call qg_fields_from_fieldset(fld,afieldset)

end subroutine qg_fields_from_fieldset_c
! ------------------------------------------------------------------------------
subroutine qg_fields_getvals_c(c_key, c_vars, c_nlocs, c_locs, c_nvals, c_vals) bind (c,name='qg_fields_getvals_f90')

implicit none
integer(c_int), intent(in)    :: c_key
type(c_ptr),value,intent(in)  :: c_vars
integer(c_int), intent(in)    :: c_nlocs
real(c_double), intent(in)    :: c_locs(2 * c_nlocs)
integer(c_int), intent(in)    :: c_nvals
real(c_double), intent(inout) :: c_vals(c_nvals)

type(qg_fields),pointer :: self
type(oops_variables) :: vars
real(kind_real) :: lats(c_nlocs), lons(c_nlocs)
integer :: ii, jj

call qg_fields_registry%get(c_key, self)
vars = oops_variables(c_vars)

ii = 0
do jj = 1, c_nlocs
  ii = ii + 1
  lats(jj) = c_locs(ii)
  ii = ii + 1
  lons(jj) = c_locs(ii)
enddo

call qg_fields_getvals(self, vars, lats, lons, c_vals)

end subroutine qg_fields_getvals_c
! ------------------------------------------------------------------------------
subroutine qg_fields_getvalsad_c(c_key, c_vars, c_nlocs, c_locs, c_nvals, c_vals) bind (c,name='qg_fields_getvalsad_f90')

implicit none
integer(c_int),intent(in)    :: c_key
type(c_ptr),value,intent(in) :: c_vars
integer(c_int), intent(in)   :: c_nlocs
real(c_double), intent(in)   :: c_locs(2 * c_nlocs)
integer(c_int), intent(in)   :: c_nvals
real(c_double), intent(in)   :: c_vals(c_nvals)

type(qg_fields),pointer :: self
type(oops_variables) :: vars
real(kind_real) :: lats(c_nlocs), lons(c_nlocs)
integer :: ii, jj

call qg_fields_registry%get(c_key, self)
vars = oops_variables(c_vars)

ii = 0
do jj = 1, c_nlocs
  ii = ii + 1
  lats(jj) = c_locs(ii)
  ii = ii + 1
  lons(jj) = c_locs(ii)
enddo

call qg_fields_getvalsad(self, vars, lats, lons, c_vals)

end subroutine qg_fields_getvalsad_c
! ------------------------------------------------------------------------------
!> Get points from fields
subroutine qg_fields_getpoint_c(c_key_fld,c_key_iter,c_nval,c_vals) bind(c,name='qg_fields_getpoint_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld         !< Fields
integer(c_int),intent(in) :: c_key_iter        !< Geometry iterator
integer(c_int),intent(in) :: c_nval            !< Number of values
real(c_double),intent(inout) :: c_vals(c_nval) !< Values

! Local variables
type(qg_fields),pointer :: fld
type(qg_geom_iter),pointer :: iter

! Interface
call qg_fields_registry%get(c_key_fld,fld)
call qg_geom_iter_registry%get(c_key_iter,iter)

! Call Fortran
call qg_fields_getpoint(fld,iter,c_nval,c_vals)

end subroutine qg_fields_getpoint_c
! ------------------------------------------------------------------------------
!> Set points for the fields
subroutine qg_fields_setpoint_c(c_key_fld,c_key_iter,c_nval,c_vals) bind(c,name='qg_fields_setpoint_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld         !< Fields
integer(c_int),intent(in) :: c_key_iter        !< Geometry iterator
integer(c_int),intent(in) :: c_nval            !< Number of values
real(c_double),intent(in) :: c_vals(c_nval)    !< Values

! Local variables
type(qg_fields),pointer :: fld
type(qg_geom_iter),pointer :: iter

! Interface
call qg_fields_registry%get(c_key_fld,fld)
call qg_geom_iter_registry%get(c_key_iter,iter)

! Call Fortran
call qg_fields_setpoint(fld,iter,c_nval,c_vals)

end subroutine qg_fields_setpoint_c
! ------------------------------------------------------------------------------
!> Serialize fields
subroutine qg_fields_serialize_c(c_key_fld,c_vsize,c_vect_fld) bind(c,name='qg_fields_serialize_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_fld            !< Fields
integer(c_int),intent(in) :: c_vsize              !< Size
real(c_double),intent(out) :: c_vect_fld(c_vsize) !< Vector

! Local variables
type(qg_fields),pointer :: fld

! Interface
call qg_fields_registry%get(c_key_fld,fld)

! Call Fortran
call qg_fields_serialize(fld,c_vsize,c_vect_fld)

end subroutine qg_fields_serialize_c
! ------------------------------------------------------------------------------
!> Deserialize fields
subroutine qg_fields_deserialize_c(c_key_self,c_vsize,c_vect_fld,c_index) bind(c,name='qg_fields_deserialize_f90')

implicit none

! Passed variables
integer(c_int),intent(in) :: c_key_self          !< Fields
integer(c_int),intent(in) :: c_vsize             !< Size
real(c_double),intent(in) :: c_vect_fld(c_vsize) !< Vector
integer(c_int), intent(inout):: c_index          !< Index

! Local variables
type(qg_fields),pointer :: self

! Interface
call qg_fields_registry%get(c_key_self,self)

! Call Fortran
call qg_fields_deserialize(self,c_vsize,c_vect_fld,c_index)

end subroutine qg_fields_deserialize_c
! ------------------------------------------------------------------------------
end module qg_fields_interface
