! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!>  Fortran module for handling generic BUMP

module oobump_mod

use config_mod
use iso_c_binding
use kinds
use missing_values_mod
use type_bump, only: bump_type
use type_nam, only: nvmax,nlmax,nc3max,nscalesmax,ndirmax,nldwvmax 
use unstructured_grid_mod

implicit none

type oobump_type
   integer :: ngrid                       !> Number of instances of BUMP
   type(bump_type),allocatable :: bump(:) !> Instances of BUMP
end type oobump_type

private
public oobump_type, create_oobump, delete_oobump, add_oobump_member, &
     & multiply_oobump_vbal, multiply_oobump_vbal_inv, multiply_oobump_vbal_ad, multiply_oobump_vbal_inv_ad, &
     & multiply_oobump_nicas, get_oobump_cv_size, multiply_oobump_nicas_sqrt, multiply_oobump_nicas_sqrt_ad, &
     & randomize_oobump_nicas, run_oobump_drivers, bump_read_conf

! ------------------------------------------------------------------------------

#define LISTED_TYPE oobump_type

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: oobump_registry

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

!> Linked list implementation
#include "oops/util/linkedList_c.f"

! ------------------------------------------------------------------------------
!  C++ interfaces
! ------------------------------------------------------------------------------

subroutine create_oobump_c(key, idx, c_conf, ens1_ne, ens1_nsub, ens2_ne, ens2_nsub) bind(c, name='create_oobump_f90')
implicit none
integer(c_int), intent(inout) :: key
integer(c_int), intent(in) :: idx
type(c_ptr), intent(in) :: c_conf
integer, intent(in) :: ens1_ne
integer, intent(in) :: ens1_nsub
integer, intent(in) :: ens2_ne
integer, intent(in) :: ens2_nsub

type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug

! Initialize BUMP registry
call oobump_registry%init()
call oobump_registry%add(key)
call oobump_registry%get(key, self)

! Get unstructured grid
call unstructured_grid_registry%get(idx, ug)

! Create BUMP
call create_oobump(self, ug, c_conf, ens1_ne, ens1_nsub, ens2_ne, ens2_nsub)

end subroutine create_oobump_c

! ------------------------------------------------------------------------------

subroutine delete_oobump_c(key) bind(c, name='delete_oobump_f90')
implicit none
integer(c_int), intent(inout) :: key

type(oobump_type), pointer :: self

! Get BUMP
call oobump_registry%get(key, self)

! Delete BUMP
call delete_oobump(self)

! Delete registry key
call oobump_registry%remove(key)

end subroutine delete_oobump_c

! ------------------------------------------------------------------------------

subroutine add_oobump_member_c(key, idx, ie, iens) bind(c, name='add_oobump_member_f90')
implicit none
integer(c_int), intent(in) :: key
integer(c_int), intent(in) :: idx
integer, intent(in) :: ie
integer, intent(in) :: iens

type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug

! Get BUMP
call oobump_registry%get(key, self)

! Get unstructured grid
call unstructured_grid_registry%get(idx, ug)

! Add BUMP member
call add_oobump_member(self, ug, ie, iens)

end subroutine add_oobump_member_c

! ------------------------------------------------------------------------------

subroutine run_oobump_drivers_c(key) bind(c, name='run_oobump_drivers_f90')
implicit none
integer(c_int), intent(in) :: key

type(oobump_type), pointer :: self

! Get BUMP
call oobump_registry%get(key, self)

! Run BUMP drivers
call run_oobump_drivers(self)

end subroutine run_oobump_drivers_c

! ------------------------------------------------------------------------------

subroutine multiply_oobump_vbal_c(key, idx) bind(c, name='multiply_oobump_vbal_f90')
implicit none
integer(c_int), intent(in) :: key
integer(c_int), intent(in) :: idx

type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug

! Get BUMP
call oobump_registry%get(key, self)

! Get unstructured grid
call unstructured_grid_registry%get(idx, ug)

! Multiply
call multiply_oobump_vbal(self, ug)

end subroutine multiply_oobump_vbal_c

! ------------------------------------------------------------------------------

subroutine multiply_oobump_vbal_inv_c(key, idx) bind(c, name='multiply_oobump_vbal_inv_f90')
implicit none
integer(c_int), intent(in) :: key
integer(c_int), intent(in) :: idx

type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug

! Get BUMP
call oobump_registry%get(key, self)

! Get unstructured grid
call unstructured_grid_registry%get(idx, ug)

! Multiply
call multiply_oobump_vbal_inv(self, ug)

end subroutine multiply_oobump_vbal_inv_c

! ------------------------------------------------------------------------------

subroutine multiply_oobump_vbal_ad_c(key, idx) bind(c, name='multiply_oobump_vbal_ad_f90')
implicit none
integer(c_int), intent(in) :: key
integer(c_int), intent(in) :: idx

type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug

! Get BUMP
call oobump_registry%get(key, self)

! Get unstructured grid
call unstructured_grid_registry%get(idx, ug)

! Multiply
call multiply_oobump_vbal_ad(self, ug)

end subroutine multiply_oobump_vbal_ad_c

! ------------------------------------------------------------------------------

subroutine multiply_oobump_vbal_inv_ad_c(key, idx) bind(c, name='multiply_oobump_vbal_inv_ad_f90')
implicit none
integer(c_int), intent(in) :: key
integer(c_int), intent(in) :: idx

type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug

! Get BUMP
call oobump_registry%get(key, self)

! Get unstructured grid
call unstructured_grid_registry%get(idx, ug)

! Multiply
call multiply_oobump_vbal_inv_ad(self, ug)

end subroutine multiply_oobump_vbal_inv_ad_c

! ------------------------------------------------------------------------------

subroutine multiply_oobump_nicas_c(key, idx) bind(c, name='multiply_oobump_nicas_f90')
implicit none
integer(c_int), intent(in) :: key
integer(c_int), intent(in) :: idx

type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug

! Get BUMP
call oobump_registry%get(key, self)

! Get unstructured grid
call unstructured_grid_registry%get(idx, ug)

! Multiply
call multiply_oobump_nicas(self, ug)

end subroutine multiply_oobump_nicas_c

! ------------------------------------------------------------------------------

subroutine get_oobump_cv_size_c(key, n) bind(c, name='get_oobump_cv_size_f90')
implicit none
integer(c_int), intent(in) :: key
integer(c_int), intent(out) :: n

type(oobump_type), pointer :: self

! Get BUMP
call oobump_registry%get(key, self)

! Get control variable size
call get_oobump_cv_size(self, n)

end subroutine get_oobump_cv_size_c

! ------------------------------------------------------------------------------

subroutine multiply_oobump_nicas_sqrt_c(key, cv, idx) bind(c, name='multiply_oobump_nicas_sqrt_f90')
implicit none
integer(c_int), intent(in) :: key
real(c_double), intent(in) :: cv(:)
integer(c_int), intent(in) :: idx

type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug

! Get BUMP
call oobump_registry%get(key, self)

! Get unstructured grid
call unstructured_grid_registry%get(idx, ug)

! Square-root multiply
call multiply_oobump_nicas_sqrt(self, cv, ug)

end subroutine multiply_oobump_nicas_sqrt_c

! ------------------------------------------------------------------------------

subroutine multiply_oobump_nicas_sqrt_ad_c(key, idx, cv) bind(c, name='multiply_oobump_nicas_sqrt_ad_f90')
implicit none
integer(c_int), intent(in) :: key
integer(c_int), intent(in) :: idx
real(c_double), intent(inout) :: cv(:)

type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug

! Get BUMP
call oobump_registry%get(key, self)

! Get unstructured grid
call unstructured_grid_registry%get(idx, ug)

! Square-root multiply
call multiply_oobump_nicas_sqrt_ad(self, ug, cv)

end subroutine multiply_oobump_nicas_sqrt_ad_c

! ------------------------------------------------------------------------------

subroutine randomize_oobump_nicas_c(key, idx) bind(c, name='randomize_oobump_nicas_f90')
implicit none
integer(c_int), intent(in) :: key
integer(c_int), intent(in) :: idx

type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug

! Get BUMP
call oobump_registry%get(key, self)

! Get unstructured grid
call unstructured_grid_registry%get(idx, ug)

! Multiply
call randomize_oobump_nicas(self, ug)

end subroutine randomize_oobump_nicas_c

! ------------------------------------------------------------------------------

subroutine get_oobump_param_c(key, nstr, cstr, idx) bind(c, name='get_oobump_param_f90')
implicit none
integer(c_int), intent(in) :: key
integer(c_int), intent(in) :: nstr
character(kind=c_char), intent(in) :: cstr(nstr)
integer(c_int), intent(in) :: idx

type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug
integer :: istr
character(len=nstr) :: param

! Get BUMP
call oobump_registry%get(key, self)

! Get unstructured grid
call unstructured_grid_registry%get(idx, ug)

! Copy string
param = ''
do istr=1,nstr
   param = trim(param)//cstr(istr)
end do

! Get parameter
call get_oobump_param(self, param, ug)

end subroutine get_oobump_param_c

! ------------------------------------------------------------------------------

subroutine set_oobump_param_c(key, nstr, cstr, idx) bind(c, name='set_oobump_param_f90')
implicit none
integer(c_int), intent(in) :: key
integer(c_int), intent(in) :: nstr
character(kind=c_char), intent(in) :: cstr(nstr)
integer(c_int), intent(in) :: idx

type(oobump_type), pointer :: self
type(unstructured_grid), pointer :: ug
integer :: istr
character(len=nstr) :: param

! Get BUMP
call oobump_registry%get(key, self)

! Get unstructured grid
call unstructured_grid_registry%get(idx, ug)

! Copy string
param = ''
do istr=1,nstr
   param = trim(param)//cstr(istr)
end do

! Set parameter
call set_oobump_param(self, param, ug)

end subroutine set_oobump_param_c

! ------------------------------------------------------------------------------
!  End C++ interfaces
! ------------------------------------------------------------------------------

subroutine create_oobump(self, ug, c_conf, ens1_ne, ens1_nsub, ens2_ne, ens2_nsub)

implicit none
type(oobump_type), intent(inout) :: self
type(unstructured_grid), intent(in) :: ug
type(c_ptr), intent(in) :: c_conf
integer, intent(in) :: ens1_ne
integer, intent(in) :: ens1_nsub
integer, intent(in) :: ens2_ne
integer, intent(in) :: ens2_nsub

integer :: igrid
real(kind_real) :: msvalr

! Allocation
self%ngrid = ug%ngrid
allocate(self%bump(self%ngrid))

do igrid=1,self%ngrid
   ! Initialize namelist
   call self%bump(igrid)%nam%init

   ! Read JSON
   call bump_read_conf(c_conf,self%bump(igrid))

   ! Get missing value
   msvalr  = missing_value(1.0_kind_real)

   ! Online setup
   call self%bump(igrid)%setup_online(ug%grid(igrid)%nmga,ug%grid(igrid)%nl0,ug%grid(igrid)%nv,ug%grid(igrid)%nts, &
 & ug%grid(igrid)%lon,ug%grid(igrid)%lat,ug%grid(igrid)%area,ug%grid(igrid)%vunit,ug%grid(igrid)%lmask,ens1_ne=ens1_ne, &
 & ens1_nsub=ens1_nsub, ens2_ne=ens2_ne, ens2_nsub=ens2_nsub,msvalr=msvalr)
end do

end subroutine create_oobump

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> Read BUMP configuration
subroutine bump_read_conf(conf,bump)

implicit none

! Passed variables
type(c_ptr), intent(in) :: conf        !< Configuration
type(bump_type), intent(inout) :: bump !< BUMP

! Local variables
integer :: n
integer,allocatable :: ivec(:)

! general_param
if (config_element_exists(conf,"datadir")) bump%nam%datadir = config_get_string(conf,1024,"datadir")
if (config_element_exists(conf,"prefix")) bump%nam%prefix = config_get_string(conf,1024,"prefix")
if (config_element_exists(conf,"verbosity")) bump%nam%verbosity = config_get_string(conf,1024,"verbosity")
if (config_element_exists(conf,"default_seed")) bump%nam%default_seed = integer_to_logical(config_get_int(conf,"default_seed"))

! driver_param
if (config_element_exists(conf,"method")) bump%nam%method = config_get_string(conf,1024,"method")
if (config_element_exists(conf,"strategy")) bump%nam%strategy = config_get_string(conf,1024,"strategy")
if (config_element_exists(conf,"new_cortrack")) bump%nam%new_cortrack = integer_to_logical(config_get_int(conf,"new_cortrack"))
if (config_element_exists(conf,"new_vbal")) bump%nam%new_vbal = integer_to_logical(config_get_int(conf,"new_vbal"))
if (config_element_exists(conf,"load_vbal")) bump%nam%load_vbal = integer_to_logical(config_get_int(conf,"load_vbal"))
if (config_element_exists(conf,"write_vbal")) bump%nam%write_vbal = integer_to_logical(config_get_int(conf,"write_vbal"))
if (config_element_exists(conf,"new_hdiag")) bump%nam%new_hdiag = integer_to_logical(config_get_int(conf,"new_hdiag"))
if (config_element_exists(conf,"write_hdiag")) bump%nam%write_hdiag = integer_to_logical(config_get_int(conf,"write_hdiag"))
if (config_element_exists(conf,"new_lct")) bump%nam%new_lct = integer_to_logical(config_get_int(conf,"new_lct"))
if (config_element_exists(conf,"write_lct")) bump%nam%write_lct = integer_to_logical(config_get_int(conf,"write_lct"))
if (config_element_exists(conf,"load_cmat")) bump%nam%load_cmat = integer_to_logical(config_get_int(conf,"load_cmat"))
if (config_element_exists(conf,"write_cmat")) bump%nam%write_cmat = integer_to_logical(config_get_int(conf,"write_cmat"))
if (config_element_exists(conf,"new_nicas")) bump%nam%new_nicas = integer_to_logical(config_get_int(conf,"new_nicas"))
if (config_element_exists(conf,"load_nicas")) bump%nam%load_nicas = integer_to_logical(config_get_int(conf,"load_nicas"))
if (config_element_exists(conf,"write_nicas")) bump%nam%write_nicas = integer_to_logical(config_get_int(conf,"write_nicas"))
if (config_element_exists(conf,"new_obsop")) bump%nam%new_obsop = integer_to_logical(config_get_int(conf,"new_obsop"))
if (config_element_exists(conf,"load_obsop")) bump%nam%load_obsop = integer_to_logical(config_get_int(conf,"load_obsop"))
if (config_element_exists(conf,"write_obsop")) bump%nam%write_obsop = integer_to_logical(config_get_int(conf,"write_obsop"))
if (config_element_exists(conf,"check_adjoints")) &
 & bump%nam%check_adjoints = integer_to_logical(config_get_int(conf,"check_adjoints"))
if (config_element_exists(conf,"check_pos_def")) &
 & bump%nam%check_pos_def = integer_to_logical(config_get_int(conf,"check_pos_def"))
if (config_element_exists(conf,"check_sqrt")) bump%nam%check_sqrt = integer_to_logical(config_get_int(conf,"check_sqrt"))
if (config_element_exists(conf,"check_dirac")) bump%nam%check_dirac = integer_to_logical(config_get_int(conf,"check_dirac"))
if (config_element_exists(conf,"check_randomization")) & 
 & bump%nam%check_randomization = integer_to_logical(config_get_int(conf,"check_randomization"))
if (config_element_exists(conf,"check_consistency")) &
 & bump%nam%check_consistency = integer_to_logical(config_get_int(conf,"check_consistency"))
if (config_element_exists(conf,"check_optimality")) &
 & bump%nam%check_optimality = integer_to_logical(config_get_int(conf,"check_optimality"))
if (config_element_exists(conf,"check_obsop")) &
 & bump%nam%check_obsop = integer_to_logical(config_get_int(conf,"check_obsop"))

! sampling_param
if (config_element_exists(conf,"sam_read")) bump%nam%sam_read = integer_to_logical(config_get_int(conf,"sam_read"))
if (config_element_exists(conf,"sam_write")) bump%nam%sam_write = integer_to_logical(config_get_int(conf,"sam_write"))
if (config_element_exists(conf,"mask_type")) bump%nam%mask_type = config_get_string(conf,1024,"mask_type")
if (config_element_exists(conf,"mask_th")) bump%nam%mask_th = config_get_real(conf,"mask_th")
if (config_element_exists(conf,"mask_check")) bump%nam%mask_check = integer_to_logical(config_get_int(conf,"mask_check"))
if (config_element_exists(conf,"draw_type")) bump%nam%draw_type = config_get_string(conf,1024,"draw_type")
if (config_element_exists(conf,"nc1")) bump%nam%nc1 = config_get_int(conf,"nc1")
if (config_element_exists(conf,"nc2")) bump%nam%nc2 = config_get_int(conf,"nc2")
if (config_element_exists(conf,"ntry")) bump%nam%ntry = config_get_int(conf,"ntry")
if (config_element_exists(conf,"nrep")) bump%nam%nrep = config_get_int(conf,"nrep")
if (config_element_exists(conf,"nc3")) bump%nam%nc3 = config_get_int(conf,"nc3")
if (config_element_exists(conf,"dc")) bump%nam%dc = config_get_real(conf,"dc")
if (config_element_exists(conf,"nl0r")) bump%nam%nl0r = config_get_int(conf,"nl0r")

! diag_param
if (config_element_exists(conf,"ne")) bump%nam%ne = config_get_int(conf,"ne")
if (config_element_exists(conf,"gau_approx")) bump%nam%gau_approx = integer_to_logical(config_get_int(conf,"gau_approx"))
if (config_element_exists(conf,"vbal_block")) then
  n = config_get_data_dimension(conf,"vbal_block")
  allocate(ivec(n))
  call config_get_int_vector(conf,"vbal_block",ivec)
  bump%nam%vbal_block(1:n) = integer_to_logical_vector(ivec)
  deallocate(ivec)
end if
if (config_element_exists(conf,"vbal_rad")) bump%nam%vbal_rad = config_get_real(conf,"vbal_rad")
if (config_element_exists(conf,"var_filter")) bump%nam%var_filter = integer_to_logical(config_get_int(conf,"var_filter"))
if (config_element_exists(conf,"var_full")) bump%nam%var_full = integer_to_logical(config_get_int(conf,"var_full"))
if (config_element_exists(conf,"var_niter")) bump%nam%var_niter = config_get_int(conf,"var_niter")
if (config_element_exists(conf,"var_rhflt")) bump%nam%var_rhflt = config_get_real(conf,"var_rhflt")
if (config_element_exists(conf,"local_diag")) bump%nam%local_diag = integer_to_logical(config_get_int(conf,"local_diag"))
if (config_element_exists(conf,"local_rad")) bump%nam%local_rad = config_get_real(conf,"local_rad")
if (config_element_exists(conf,"displ_diag")) bump%nam%displ_diag = integer_to_logical(config_get_int(conf,"displ_diag"))
if (config_element_exists(conf,"displ_rad")) bump%nam%displ_rad = config_get_real(conf,"displ_rad")
if (config_element_exists(conf,"displ_niter")) bump%nam%displ_niter = config_get_int(conf,"displ_niter")
if (config_element_exists(conf,"displ_rhflt")) bump%nam%displ_rhflt = config_get_real(conf,"displ_rhflt")

! fit_param
if (config_element_exists(conf,"minim_algo")) bump%nam%minim_algo = config_get_string(conf,1024,"minim_algo")
if (config_element_exists(conf,"double_fit")) then
  n = config_get_data_dimension(conf,"double_fit")
  allocate(ivec(n))
  call config_get_int_vector(conf,"double_fit",ivec)
  bump%nam%double_fit(1:n) = integer_to_logical_vector(ivec)
  deallocate(ivec)
end if
if (config_element_exists(conf,"lhomh")) bump%nam%lhomh = integer_to_logical(config_get_int(conf,"lhomh"))
if (config_element_exists(conf,"lhomv")) bump%nam%lhomv = integer_to_logical(config_get_int(conf,"lhomv"))
if (config_element_exists(conf,"rvflt")) bump%nam%rvflt = config_get_real(conf,"rvflt")
if (config_element_exists(conf,"lct_nscales")) bump%nam%lct_nscales = config_get_int(conf,"lct_nscales")
if (config_element_exists(conf,"lct_diag")) then
  allocate(ivec(config_get_data_dimension(conf,"lct_diag")))
  call config_get_int_vector(conf,"lct_diag",ivec)
  bump%nam%lct_diag = integer_to_logical_vector(ivec)
  deallocate(ivec)
end if

! nicas_param
if (config_element_exists(conf,"lsqrt")) bump%nam%lsqrt = integer_to_logical(config_get_int(conf,"lsqrt"))
if (config_element_exists(conf,"resol")) bump%nam%resol = config_get_real(conf,"resol")
if (config_element_exists(conf,"fast_sampling")) bump%nam%fast_sampling =  integer_to_logical(config_get_int(conf, &
 & "fast_sampling"))
if (config_element_exists(conf,"subsamp")) bump%nam%subsamp = config_get_string(conf,1024,"subsamp")
if (config_element_exists(conf,"nicas_interp")) bump%nam%nicas_interp = config_get_string(conf,1024,"nicas_interp")
if (config_element_exists(conf,"network")) bump%nam%network = integer_to_logical(config_get_int(conf,"network"))
if (config_element_exists(conf,"mpicom")) bump%nam%mpicom = config_get_int(conf,"mpicom")
if (config_element_exists(conf,"advmode")) bump%nam%advmode = config_get_int(conf,"advmode")
if (config_element_exists(conf,"forced_radii")) bump%nam%forced_radii = integer_to_logical(config_get_int(conf,"forced_radii"))
if (config_element_exists(conf,"rh")) bump%nam%rh = config_get_real(conf,"rh")
if (config_element_exists(conf,"rv")) bump%nam%rv = config_get_real(conf,"rv")
if (config_element_exists(conf,"write_grids")) bump%nam%write_grids = integer_to_logical(config_get_int(conf,"write_grids"))
if (config_element_exists(conf,"ndir")) bump%nam%ndir = config_get_int(conf,"ndir")
if (config_element_exists(conf,"londir")) then
  n = config_get_data_dimension(conf,"londir")
  call config_get_double_vector(conf,"londir",bump%nam%londir(1:n))
end if
if (config_element_exists(conf,"latdir")) then
  n = config_get_data_dimension(conf,"latdir")
  call config_get_double_vector(conf,"latdir",bump%nam%latdir(1:n))
end if
if (config_element_exists(conf,"levdir")) then
  n = config_get_data_dimension(conf,"levdir")
  call config_get_int_vector(conf,"levdir",bump%nam%levdir(1:n))
end if
if (config_element_exists(conf,"ivdir")) then
  n = config_get_data_dimension(conf,"ivdir")
  call config_get_int_vector(conf,"ivdir",bump%nam%ivdir(1:n))
end if
if (config_element_exists(conf,"itsdir")) then
  n = config_get_data_dimension(conf,"itsdir")
  call config_get_int_vector(conf,"itsdir",bump%nam%itsdir(1:n))
end if

! output_param
if (config_element_exists(conf,"nldwh")) bump%nam%nldwh = config_get_int(conf,"nldwh")
if (config_element_exists(conf,"il_ldwh")) then
  n = config_get_data_dimension(conf,"il_ldwh")
  call config_get_int_vector(conf,"il_ldwh",bump%nam%il_ldwh(1:n))
end if
if (config_element_exists(conf,"ic_ldwh")) then
  n = config_get_data_dimension(conf,"ic_ldwh")
  call config_get_int_vector(conf,"ic_ldwh",bump%nam%ic_ldwh(1:n))
end if
if (config_element_exists(conf,"nldwv")) bump%nam%nldwv = config_get_int(conf,"nldwv")
if (config_element_exists(conf,"lon_ldwv")) then
  n = config_get_data_dimension(conf,"lon_ldwv")
  call config_get_double_vector(conf,"lon_ldwv",bump%nam%lon_ldwv(1:n))
end if
if (config_element_exists(conf,"lat_ldwv")) then
  n = config_get_data_dimension(conf,"lat_ldwv")
  call config_get_double_vector(conf,"lat_ldwv",bump%nam%lat_ldwv(1:n))
end if
if (config_element_exists(conf,"diag_rhflt")) bump%nam%diag_rhflt = config_get_real(conf,"diag_rhflt")
if (config_element_exists(conf,"diag_interp")) bump%nam%diag_interp = config_get_string(conf,1024,"diag_interp")
if (config_element_exists(conf,"field_io")) bump%nam%field_io = integer_to_logical(config_get_int(conf,"field_io"))
if (config_element_exists(conf,"split_io")) bump%nam%split_io = integer_to_logical(config_get_int(conf,"split_io"))
if (config_element_exists(conf,"grid_output")) bump%nam%grid_output = integer_to_logical(config_get_int(conf,"grid_output"))
if (config_element_exists(conf,"grid_resol")) bump%nam%grid_resol = config_get_real(conf,"grid_resol")
if (config_element_exists(conf,"grid_interp")) bump%nam%grid_interp = config_get_string(conf,1024,"grid_interp")

end subroutine bump_read_conf
!-------------------------------------------------------------------------------
!> Convert integer to logical
function integer_to_logical(i)

implicit none

! Passed variables
integer,intent(in) :: i

! Returned argument
logical :: integer_to_logical

if (i==0) then
   integer_to_logical = .false.
elseif (i==1) then
   integer_to_logical = .true.
else
   call abor1_ftn('wrong integer in integer_to_logical')
end if

end function integer_to_logical
!-------------------------------------------------------------------------------
!> Convert integer to logical (vector
function integer_to_logical_vector(i)

implicit none

! Passed variables
integer,intent(in),dimension(:) :: i

! Returned argument
logical :: integer_to_logical_vector(size(i))

! Local variables
integer :: j

do j=1,size(i)
  if (i(j)==0) then
    integer_to_logical_vector(j) = .false.
  elseif (i(j)==1) then
    integer_to_logical_vector(j) = .true.
  else
     call abor1_ftn('wrong integer in integer_to_logical_vector')
  end if
end do

end function integer_to_logical_vector

!-------------------------------------------------------------------------------

subroutine delete_oobump(self)
implicit none
type(oobump_type), intent(inout) :: self
integer :: igrid

! Deallocate BUMP
if (allocated(self%bump)) then
   do igrid=1,self%ngrid  
      call self%bump(igrid)%dealloc
   end do
   deallocate(self%bump)
end if

end subroutine delete_oobump

!-------------------------------------------------------------------------------

subroutine add_oobump_member(self,ug,ie,iens)
implicit none
type(oobump_type), intent(inout) :: self
type(unstructured_grid), intent(inout) :: ug
integer, intent(in) :: ie
integer, intent(in) :: iens
integer :: igrid

! Add member
do igrid=1,self%ngrid
   call self%bump(igrid)%add_member(ug%grid(igrid)%fld,ie,iens)
end do

end subroutine add_oobump_member

!-------------------------------------------------------------------------------

subroutine run_oobump_drivers(self)
implicit none
type(oobump_type), intent(inout) :: self
integer :: igrid

! Run BUMP drivers
do igrid=1,self%ngrid
   call self%bump(igrid)%run_drivers
end do

end subroutine run_oobump_drivers

!-------------------------------------------------------------------------------

subroutine multiply_oobump_vbal(self,ug)
implicit none
type(oobump_type), intent(inout) :: self
type(unstructured_grid), intent(inout) :: ug
integer :: igrid

! Apply vertical balance
do igrid=1,self%ngrid
   call self%bump(igrid)%apply_vbal(ug%grid(igrid)%fld)
end do

end subroutine multiply_oobump_vbal

!-------------------------------------------------------------------------------

subroutine multiply_oobump_vbal_inv(self,ug)
implicit none
type(oobump_type), intent(inout) :: self
type(unstructured_grid), intent(inout) :: ug
integer :: igrid

! Apply vertical balance, inverse
do igrid=1,self%ngrid
   call self%bump(igrid)%apply_vbal_inv(ug%grid(igrid)%fld)
end do

end subroutine multiply_oobump_vbal_inv

!-------------------------------------------------------------------------------

subroutine multiply_oobump_vbal_ad(self,ug)
implicit none
type(oobump_type), intent(inout) :: self
type(unstructured_grid), intent(inout) :: ug
integer :: igrid

! Apply vertical balance, adjoint
do igrid=1,self%ngrid
   call self%bump(igrid)%apply_vbal_ad(ug%grid(igrid)%fld)
end do

end subroutine multiply_oobump_vbal_ad

!-------------------------------------------------------------------------------

subroutine multiply_oobump_vbal_inv_ad(self,ug)
implicit none
type(oobump_type), intent(inout) :: self
type(unstructured_grid), intent(inout) :: ug
integer :: igrid

! Apply vertical balance, inverse adjoint
do igrid=1,self%ngrid
   call self%bump(igrid)%apply_vbal_inv_ad(ug%grid(igrid)%fld)
end do

end subroutine multiply_oobump_vbal_inv_ad

!-------------------------------------------------------------------------------

subroutine multiply_oobump_nicas(self,ug)
implicit none
type(oobump_type), intent(inout) :: self
type(unstructured_grid), intent(inout) :: ug
integer :: igrid

! Apply NICAS
do igrid=1,self%ngrid
   call self%bump(igrid)%apply_nicas(ug%grid(igrid)%fld)
end do

end subroutine multiply_oobump_nicas

!-------------------------------------------------------------------------------

subroutine get_oobump_cv_size(self,n)
implicit none
type(oobump_type), intent(inout) :: self
integer, intent(out) :: n
integer :: igrid,nn

! Add control variable sizes for each grid
n = 0
do igrid=1,self%ngrid
   call self%bump(igrid)%get_cv_size(nn)
   n = n+nn
end do

end subroutine get_oobump_cv_size

!-------------------------------------------------------------------------------

subroutine multiply_oobump_nicas_sqrt(self,cv,ug)
implicit none
type(oobump_type), intent(inout) :: self
real(kind_real), intent(in) :: cv(:)
type(unstructured_grid), intent(inout) :: ug
integer :: offset,igrid,nn

! Initialization
offset = 0

do igrid=1,self%ngrid
   ! Get control variable size for this grid
   call self%bump(igrid)%get_cv_size(nn)

   ! Apply NICAS square-root
   call self%bump(igrid)%apply_nicas_sqrt(cv(offset+1:offset+nn), ug%grid(igrid)%fld)

   ! Update
   offset = offset+nn
end do

end subroutine multiply_oobump_nicas_sqrt

!-------------------------------------------------------------------------------

subroutine multiply_oobump_nicas_sqrt_ad(self,ug,cv)
implicit none
type(oobump_type), intent(inout) :: self
type(unstructured_grid), intent(in) :: ug
real(kind_real), intent(inout) :: cv(:)
integer :: offset,igrid,nn

! Initialization
offset = 0

do igrid=1,self%ngrid
   ! Get control variable size for this grid
   call self%bump(igrid)%get_cv_size(nn)

   ! Apply NICAS square-root
   call self%bump(igrid)%apply_nicas_sqrt_ad(ug%grid(igrid)%fld, cv(offset+1:offset+nn))

   ! Update
   offset = offset+nn
end do

end subroutine multiply_oobump_nicas_sqrt_ad

!-------------------------------------------------------------------------------

subroutine randomize_oobump_nicas(self,ug)
implicit none
type(oobump_type), intent(inout) :: self
type(unstructured_grid), intent(out) :: ug
integer :: igrid

! Randomize NICAS
do igrid=1,self%ngrid
   call self%bump(igrid)%randomize(ug%grid(igrid)%fld)
end do

end subroutine randomize_oobump_nicas

!-------------------------------------------------------------------------------

subroutine get_oobump_param(self,param,ug)
implicit none
type(oobump_type), intent(inout) :: self
character(len=*), intent(in) :: param
type(unstructured_grid), intent(inout) :: ug
integer :: igrid

! Get parameter
do igrid=1,self%ngrid
   call self%bump(igrid)%get_parameter(param,ug%grid(igrid)%fld)
end do

end subroutine get_oobump_param

!-------------------------------------------------------------------------------

subroutine set_oobump_param(self,param,ug)
implicit none
type(oobump_type), intent(inout) :: self
character(len=*),intent(in) :: param
type(unstructured_grid), intent(in) :: ug
integer :: igrid

! Set parameter
do igrid=1,self%ngrid
   call self%bump(igrid)%set_parameter(param,ug%grid(igrid)%fld)
end do

end subroutine set_oobump_param

!-------------------------------------------------------------------------------

end module oobump_mod
