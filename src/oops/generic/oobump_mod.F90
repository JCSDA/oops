! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!>  Fortran module for handling generic BUMP

module oobump_mod

use iso_c_binding
use kinds
use config_mod
use unstructured_grid_mod
#ifdef notdef_
use mpi
#endif
use type_bump, only: bump_type

implicit none

type oobump_type
   integer :: ngrid                       !> Number of instances of BUMP
   type(bump_type),allocatable :: bump(:) !> Instances of BUMP
end type oobump_type

private
public oobump_type, create_oobump, delete_oobump, add_oobump_member, &
     & multiply_oobump_vbal, multiply_oobump_vbal_inv, multiply_oobump_vbal_ad, multiply_oobump_vbal_inv_ad, &
     & multiply_oobump_nicas, run_oobump_drivers, bump_read_conf

#ifndef notdef_
  INCLUDE 'mpif.h'
#endif  
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
call oobump_registry%get(key,self)

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
call oobump_registry%get(key,self)

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
call oobump_registry%get(key,self)

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
call oobump_registry%get(key,self)

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
call oobump_registry%get(key,self)

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
call oobump_registry%get(key,self)

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
call oobump_registry%get(key,self)

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
call oobump_registry%get(key,self)

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
call oobump_registry%get(key,self)

! Get unstructured grid
call unstructured_grid_registry%get(idx, ug)

! Multiply
call multiply_oobump_nicas(self, ug)

end subroutine multiply_oobump_nicas_c

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
call oobump_registry%get(key,self)

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

! Allocation
self%ngrid = ug%ngrid
allocate(self%bump(self%ngrid))

do igrid=1,self%ngrid
   ! Initialize namelist
   call self%bump(igrid)%nam%init

   ! Read JSON
   call bump_read_conf(c_conf,self%bump(igrid))

   ! Online setup
   call self%bump(igrid)%setup_online(mpi_comm_world,ug%grid(igrid)%nmga,ug%grid(igrid)%nl0,ug%grid(igrid)%nv,ug%grid(igrid)%nts, &
 & ug%grid(igrid)%lon,ug%grid(igrid)%lat,ug%grid(igrid)%area,ug%grid(igrid)%vunit,ug%grid(igrid)%lmask, &
 & ens1_ne=ens1_ne,ens1_nsub=ens1_nsub,ens2_ne=ens2_ne,ens2_nsub=ens2_nsub)
end do

end subroutine create_oobump

!-------------------------------------------------------------------------------

subroutine bump_read_conf(c_conf,bump)
implicit none
type(c_ptr), intent(in) :: c_conf
type(bump_type), intent(inout) :: bump
integer :: il,its,iscales,ildwh,ildwv,idir
character(len=3) :: ilchar,itschar,ildwhchar,ildwvchar,idirchar

! Setup from configuration

! general_param
if (config_element_exists(c_conf,"prefix")) bump%nam%prefix = config_get_string(c_conf,1024,"prefix")
if (config_element_exists(c_conf,"default_seed")) bump%nam%default_seed = integer_to_logical(config_get_int(c_conf,"default_seed"))

! driver_param
if (config_element_exists(c_conf,"method")) bump%nam%method = config_get_string(c_conf,1024,"method")
if (config_element_exists(c_conf,"strategy")) bump%nam%strategy = config_get_string(c_conf,1024,"strategy")
if (config_element_exists(c_conf,"new_vbal")) bump%nam%new_vbal = integer_to_logical(config_get_int(c_conf,"new_vbal"))
if (config_element_exists(c_conf,"new_hdiag")) bump%nam%new_hdiag = integer_to_logical(config_get_int(c_conf,"new_hdiag"))
if (config_element_exists(c_conf,"new_lct")) bump%nam%new_lct = integer_to_logical(config_get_int(c_conf,"new_lct"))
if (config_element_exists(c_conf,"new_param")) bump%nam%new_param = integer_to_logical(config_get_int(c_conf,"new_param"))
if (config_element_exists(c_conf,"new_obsop")) bump%nam%new_obsop = integer_to_logical(config_get_int(c_conf,"new_obsop"))
if (config_element_exists(c_conf,"check_adjoints")) &
 & bump%nam%check_adjoints = integer_to_logical(config_get_int(c_conf,"check_adjoints"))
if (config_element_exists(c_conf,"check_pos_def")) &
 & bump%nam%check_pos_def = integer_to_logical(config_get_int(c_conf,"check_pos_def"))
if (config_element_exists(c_conf,"check_sqrt")) bump%nam%check_sqrt = integer_to_logical(config_get_int(c_conf,"check_sqrt"))
if (config_element_exists(c_conf,"check_dirac")) bump%nam%check_dirac = integer_to_logical(config_get_int(c_conf,"check_dirac"))
if (config_element_exists(c_conf,"check_randomization")) & 
 & bump%nam%check_randomization = integer_to_logical(config_get_int(c_conf,"check_randomization"))
if (config_element_exists(c_conf,"check_consistency")) &
 & bump%nam%check_consistency = integer_to_logical(config_get_int(c_conf,"check_consistency"))
if (config_element_exists(c_conf,"check_optimality")) &
 & bump%nam%check_optimality = integer_to_logical(config_get_int(c_conf,"check_optimality"))

! sampling_param
if (config_element_exists(c_conf,"sam_read")) bump%nam%sam_read = integer_to_logical(config_get_int(c_conf,"sam_read"))
if (config_element_exists(c_conf,"sam_write")) bump%nam%sam_write = integer_to_logical(config_get_int(c_conf,"sam_write"))
if (config_element_exists(c_conf,"mask_type")) bump%nam%mask_type = config_get_string(c_conf,1024,"mask_type")
if (config_element_exists(c_conf,"mask_th")) bump%nam%mask_th = config_get_real(c_conf,"mask_th")
if (config_element_exists(c_conf,"mask_check")) bump%nam%mask_check = integer_to_logical(config_get_int(c_conf,"mask_check"))
if (config_element_exists(c_conf,"draw_type")) bump%nam%draw_type = config_get_string(c_conf,1024,"draw_type")
if (config_element_exists(c_conf,"nc1")) bump%nam%nc1 = config_get_int(c_conf,"nc1")
if (config_element_exists(c_conf,"nc2")) bump%nam%nc2 = config_get_int(c_conf,"nc2")
if (config_element_exists(c_conf,"ntry")) bump%nam%ntry = config_get_int(c_conf,"ntry")
if (config_element_exists(c_conf,"nrep")) bump%nam%nrep = config_get_int(c_conf,"nrep")
if (config_element_exists(c_conf,"nc3")) bump%nam%nc3 = config_get_int(c_conf,"nc3")
if (config_element_exists(c_conf,"dc")) bump%nam%dc = config_get_real(c_conf,"dc")
if (config_element_exists(c_conf,"nl0r")) bump%nam%nl0r = config_get_int(c_conf,"nl0r")

! diag_param
if (config_element_exists(c_conf,"ne")) bump%nam%ne = config_get_int(c_conf,"ne")
if (config_element_exists(c_conf,"gau_approx")) bump%nam%gau_approx = integer_to_logical(config_get_int(c_conf,"gau_approx"))
if (config_element_exists(c_conf,"var_diag")) bump%nam%var_diag = integer_to_logical(config_get_int(c_conf,"var_diag"))
if (config_element_exists(c_conf,"var_filter")) bump%nam%var_filter = integer_to_logical(config_get_int(c_conf,"var_filter"))
if (config_element_exists(c_conf,"var_full")) bump%nam%var_full = integer_to_logical(config_get_int(c_conf,"var_full"))
if (config_element_exists(c_conf,"var_niter")) bump%nam%var_niter = config_get_int(c_conf,"var_niter")
   if (config_element_exists(c_conf,"var_rhflt")) bump%nam%var_rhflt = config_get_real(c_conf,"var_rhflt")
if (config_element_exists(c_conf,"local_diag")) bump%nam%local_diag = integer_to_logical(config_get_int(c_conf,"local_diag"))
if (bump%nam%local_diag) then
   if (config_element_exists(c_conf,"local_rad")) bump%nam%local_rad = config_get_real(c_conf,"local_rad")
end if
if (config_element_exists(c_conf,"displ_diag")) bump%nam%displ_diag = integer_to_logical(config_get_int(c_conf,"displ_diag"))
if (bump%nam%displ_diag) then
   if (config_element_exists(c_conf,"displ_rad")) bump%nam%displ_rad = config_get_real(c_conf,"displ_rad")
   if (config_element_exists(c_conf,"displ_niter")) bump%nam%displ_niter = config_get_int(c_conf,"displ_niter")
   if (config_element_exists(c_conf,"displ_rhflt")) bump%nam%displ_rhflt = config_get_real(c_conf,"displ_rhflt")
   if (config_element_exists(c_conf,"displ_tol")) bump%nam%displ_tol = config_get_real(c_conf,"displ_tol")
end if

! fit_param
if (config_element_exists(c_conf,"minim_algo")) bump%nam%minim_algo = config_get_string(c_conf,1024,"minim_algo")
if (config_element_exists(c_conf,"lhomh")) bump%nam%lhomh = integer_to_logical(config_get_int(c_conf,"lhomh"))
if (config_element_exists(c_conf,"lhomv")) bump%nam%lhomv = integer_to_logical(config_get_int(c_conf,"lhomv"))
if (config_element_exists(c_conf,"rvflt")) bump%nam%rvflt = config_get_real(c_conf,"rvflt")

! nicas_param
if (config_element_exists(c_conf,"lsqrt")) bump%nam%lsqrt = integer_to_logical(config_get_int(c_conf,"lsqrt"))
if (config_element_exists(c_conf,"resol")) bump%nam%resol = config_get_real(c_conf,"resol")
if (config_element_exists(c_conf,"nicas_interp")) bump%nam%nicas_interp = config_get_string(c_conf,1024,"nicas_interp")
if (config_element_exists(c_conf,"network")) bump%nam%network = integer_to_logical(config_get_int(c_conf,"network"))
if (config_element_exists(c_conf,"mpicom")) bump%nam%mpicom = config_get_int(c_conf,"mpicom")
if (config_element_exists(c_conf,"advmode")) bump%nam%advmode = config_get_int(c_conf,"advmode")
if (config_element_exists(c_conf,"forced_radii")) bump%nam%forced_radii = integer_to_logical(config_get_int(c_conf,"forced_radii"))
if (config_element_exists(c_conf,"rh")) bump%nam%rh = config_get_real(c_conf,"rh")
if (config_element_exists(c_conf,"rv")) bump%nam%rv = config_get_real(c_conf,"rv")
if (config_element_exists(c_conf,"ndir")) bump%nam%ndir = config_get_int(c_conf,"ndir")
do idir=1,bump%nam%ndir
   write(idirchar,'(i3)') idir
   if (config_element_exists(c_conf,"londir("//trim(adjustl(idirchar))//")")) &
 & bump%nam%londir(idir) = config_get_real(c_conf,"londir("//trim(adjustl(idirchar))//")")
   if (config_element_exists(c_conf,"latdir("//trim(adjustl(idirchar))//")")) & 
 & bump%nam%latdir(idir) = config_get_real(c_conf,"latdir("//trim(adjustl(idirchar))//")")
   if (config_element_exists(c_conf,"levdir("//trim(adjustl(idirchar))//")")) &
 & bump%nam%levdir(idir) = config_get_int(c_conf,"levdir("//trim(adjustl(idirchar))//")")
   if (config_element_exists(c_conf,"ivdir("//trim(adjustl(idirchar))//")")) &
 & bump%nam%ivdir(idir) = config_get_int(c_conf,"ivdir("//trim(adjustl(idirchar))//")")
   if (config_element_exists(c_conf,"itsdir("//trim(adjustl(idirchar))//")")) &
 & bump%nam%itsdir(idir) = config_get_int(c_conf,"itsdir("//trim(adjustl(idirchar))//")")
end do

! output_param
if (config_element_exists(c_conf,"nldwh")) bump%nam%nldwh = config_get_int(c_conf,"nldwh")
do ildwh=1,bump%nam%nldwh
   write(ildwvchar,'(i3)') ildwh
   if (config_element_exists(c_conf,"il_ldwh("//trim(adjustl(ildwhchar))//")")) &
 & bump%nam%il_ldwh(ildwh) = config_get_int(c_conf,"il_ldwh("//trim(adjustl(ildwhchar))//")")
   if (config_element_exists(c_conf,"ic_ldwh("//trim(adjustl(ildwhchar))//")")) &
 & bump%nam%ic_ldwh(ildwh) = config_get_int(c_conf,"ic_ldwh("//trim(adjustl(ildwhchar))//")")
end do
if (config_element_exists(c_conf,"nldwv")) bump%nam%nldwv = config_get_int(c_conf,"nldwv")
do ildwv=1,bump%nam%nldwv
   write(ildwvchar,'(i3)') ildwv
   if (config_element_exists(c_conf,"lon_ldwv("//trim(adjustl(ildwvchar))//")")) &
 & bump%nam%lon_ldwv(ildwv) = config_get_real(c_conf,"lon_ldwv("//trim(adjustl(ildwvchar))//")")
   if (config_element_exists(c_conf,"lat_ldwv("//trim(adjustl(ildwvchar))//")")) &
 & bump%nam%lat_ldwv(ildwv) = config_get_real(c_conf,"lat_ldwv("//trim(adjustl(ildwvchar))//")")
end do
if (config_element_exists(c_conf,"diag_rhflt")) bump%nam%diag_rhflt = config_get_real(c_conf,"diag_rhflt")
if (config_element_exists(c_conf,"diag_interp")) bump%nam%diag_interp = config_get_string(c_conf,1024,"diag_interp")
if (config_element_exists(c_conf,"field_io")) bump%nam%field_io = integer_to_logical(config_get_int(c_conf,"field_io"))
if (config_element_exists(c_conf,"split_io")) bump%nam%split_io = integer_to_logical(config_get_int(c_conf,"split_io"))
if (config_element_exists(c_conf,"grid_output")) bump%nam%grid_output = integer_to_logical(config_get_int(c_conf,"grid_output"))
if (bump%nam%grid_output) then
   if (config_element_exists(c_conf,"grid_resol")) bump%nam%grid_resol = config_get_real(c_conf,"grid_resol")
   if (config_element_exists(c_conf,"grid_interp")) bump%nam%grid_interp = config_get_string(c_conf,1024,"grid_interp")
end if

end subroutine bump_read_conf

!-------------------------------------------------------------------------------

logical function integer_to_logical(i)
implicit none
integer,intent(in) :: i

if (i==0) then
   integer_to_logical = .false.
elseif (i==1) then
   integer_to_logical = .true.
else
   call abor1_ftn('wrong integer in integer_to_logical')
end if

end function integer_to_logical

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
endif

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
type(oobump_type), intent(in) :: self
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
type(oobump_type), intent(in) :: self
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
type(oobump_type), intent(in) :: self
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
type(oobump_type), intent(in) :: self
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
type(oobump_type), intent(in) :: self
type(unstructured_grid), intent(inout) :: ug
integer :: igrid

! Apply NICAS
do igrid=1,self%ngrid
   call self%bump(igrid)%apply_nicas(ug%grid(igrid)%fld)
end do

end subroutine multiply_oobump_nicas

!-------------------------------------------------------------------------------

subroutine get_oobump_param(self,param,ug)
implicit none
type(oobump_type), intent(in) :: self
character(len=*),intent(in) :: param
type(unstructured_grid), intent(inout) :: ug
integer :: igrid

! Get parameter
do igrid=1,self%ngrid
   call self%bump(igrid)%get_parameter(param,ug%grid(igrid)%fld)
end do

end subroutine get_oobump_param

!-------------------------------------------------------------------------------

end module oobump_mod
