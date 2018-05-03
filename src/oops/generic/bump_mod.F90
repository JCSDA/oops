! (C) Copyright 2017 UCAR
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 

!>  Fortran module for handling generic BUMP localization

module bump_mod

use iso_c_binding
use kinds
use config_mod
use unstructured_grid_mod
use mpi
use tools_const, only: req
use tools_missing, only: msi,msr
use type_bump, only: bump_type

implicit none
private
public create_bump, delete_bump, bump_multiply

! ------------------------------------------------------------------------------

#define LISTED_TYPE bump_type

!> Linked list interface - defines registry_t type
#include "util/linkedList_i.f"
!#include "linkedList.intf.h"

!> Global registry
type(registry_t) :: bump_registry

!-------------------------------------------------------------------------------
contains
!-------------------------------------------------------------------------------

!> Linked list implementation
#include "util/linkedList_c.f"
!#include "linkedList.h"

! ------------------------------------------------------------------------------
!  C++ interfaces
! ------------------------------------------------------------------------------

subroutine create_bump_c(key, c_conf, nmga, nl0, nv, nts, lon, lat, area, vunit, imask, ens1_ne, ens1) &
 & bind(c, name='create_bump_f90')
implicit none
integer(c_int), intent(inout) :: key
type(c_ptr), intent(in) :: c_conf
integer(c_int), intent(in) :: nmga
integer(c_int), intent(in) :: nl0
integer(c_int), intent(in) :: nv
integer(c_int), intent(in) :: nts
real(c_double), intent(in) :: lon(nmga)
real(c_double), intent(in) :: lat(nmga)
real(c_double), intent(in) :: area(nmga)
real(c_double), intent(in) :: vunit(nmga*nl0)
integer(c_int), intent(in) :: imask(nmga*nl0)
integer(c_int), intent(in) :: ens1_ne
real(c_double), intent(in) :: ens1(nmga*nl0*nv*nts*ens1_ne)

type(bump_type), pointer :: self

! Initialize bump registry
call bump_registry%init()
call bump_registry%add(key)
call bump_registry%get(key,self)

! Create bump object
call create_bump(self, c_conf, nmga, nl0, nv, nts, lon, lat, area, vunit, imask, ens1_ne, ens1)

end subroutine create_bump_c

! ------------------------------------------------------------------------------

subroutine delete_bump_c(key) bind(c, name='delete_bump_f90')
implicit none
integer(c_int), intent(inout) :: key

type(bump_type), pointer :: self

call bump_registry%get(key,self)
call delete_bump(self)
call bump_registry%remove(key)

end subroutine delete_bump_c

! ------------------------------------------------------------------------------

subroutine bump_multiply_c(key, idx) bind(c, name='bump_multiply_f90')
implicit none
integer(c_int), intent(in) :: key
integer(c_int), intent(in) :: idx

type(bump_type), pointer :: self
type(unstructured_grid), pointer :: ug

call bump_registry%get(key,self)
call unstructured_grid_registry%get(idx, ug)
call bump_multiply(self, ug)

end subroutine bump_multiply_c

! ------------------------------------------------------------------------------
!  End C++ interfaces
! ------------------------------------------------------------------------------

subroutine create_bump(self, c_conf, nmga, nl0, nv, nts, lon, lat, area, vunit, imask, ens1_ne, ens1)

implicit none
type(bump_type), intent(inout) :: self
type(c_ptr), intent(in) :: c_conf
integer, intent(in) :: nmga
integer, intent(in) :: nl0
integer, intent(in) :: nv
integer, intent(in) :: nts
real(kind=kind_real), intent(in) :: lon(nmga)
real(kind=kind_real), intent(in) :: lat(nmga)
real(kind=kind_real), intent(in) :: area(nmga)
real(kind=kind_real), intent(in) :: vunit(nmga*nl0)
integer, intent(in) :: imask(nmga*nl0)
integer, intent(in) :: ens1_ne
real(kind=kind_real), intent(in) :: ens1(nmga*nl0*nv*nts*ens1_ne)

! Initialize namelist
call self%nam%init

! Read JSON
call bump_read_conf(c_conf,self)

! Online setup
call self%setup_online(mpi_comm_world,nmga,nl0,nv,nts,lon,lat,area,vunit,imask,ens1_ne,ens1)

end subroutine create_bump

!-------------------------------------------------------------------------------

subroutine bump_read_conf(c_conf,bump)
implicit none
type(c_ptr), intent(in) :: c_conf
type(bump_type), intent(inout) :: bump
integer :: il,its,iscales,ildwh,ildwv,idir
character(len=3) :: ilchar,itschar,ildwhchar,ildwvchar,idirchar

! Setup from configuration

! general_param
bump%nam%prefix = config_get_string(c_conf,1024,"prefix")
bump%nam%default_seed = integer_to_logical(config_get_int(c_conf,"default_seed"))

! driver_param
bump%nam%method = config_get_string(c_conf,1024,"method")
bump%nam%strategy = config_get_string(c_conf,1024,"strategy")
bump%nam%new_hdiag = integer_to_logical(config_get_int(c_conf,"new_hdiag"))
bump%nam%new_param = integer_to_logical(config_get_int(c_conf,"new_param"))
bump%nam%check_adjoints = integer_to_logical(config_get_int(c_conf,"check_adjoints"))
bump%nam%check_pos_def = integer_to_logical(config_get_int(c_conf,"check_pos_def"))
bump%nam%check_sqrt = integer_to_logical(config_get_int(c_conf,"check_sqrt"))
bump%nam%check_dirac = integer_to_logical(config_get_int(c_conf,"check_dirac"))
bump%nam%check_randomization = integer_to_logical(config_get_int(c_conf,"check_randomization"))
bump%nam%check_consistency = integer_to_logical(config_get_int(c_conf,"check_consistency"))
bump%nam%check_optimality = integer_to_logical(config_get_int(c_conf,"check_optimality"))
bump%nam%new_lct = integer_to_logical(config_get_int(c_conf,"new_lct"))
bump%nam%new_obsop = integer_to_logical(config_get_int(c_conf,"new_obsop"))

! sampling_param
bump%nam%sam_read = integer_to_logical(config_get_int(c_conf,"sam_read"))
bump%nam%sam_write = integer_to_logical(config_get_int(c_conf,"sam_write"))
bump%nam%mask_type = config_get_string(c_conf,1024,"mask_type")
bump%nam%mask_th = config_get_real(c_conf,"mask_th")
bump%nam%mask_check = integer_to_logical(config_get_int(c_conf,"mask_check"))
bump%nam%draw_type = config_get_string(c_conf,1024,"draw_type")
bump%nam%nc1 = config_get_int(c_conf,"nc1")
bump%nam%ntry = config_get_int(c_conf,"ntry")
bump%nam%nrep = config_get_int(c_conf,"nrep")
bump%nam%nc3 = config_get_int(c_conf,"nc3")
bump%nam%dc = config_get_real(c_conf,"dc")/req
bump%nam%nl0r = config_get_int(c_conf,"nl0r")

! diag_param
bump%nam%ne = config_get_int(c_conf,"ne")
bump%nam%gau_approx = integer_to_logical(config_get_int(c_conf,"gau_approx"))
bump%nam%full_var = integer_to_logical(config_get_int(c_conf,"full_var"))
bump%nam%local_diag = integer_to_logical(config_get_int(c_conf,"local_diag"))
if (bump%nam%local_diag) bump%nam%local_rad = config_get_real(c_conf,"local_rad")/req
bump%nam%displ_diag = integer_to_logical(config_get_int(c_conf,"displ_diag"))
if (bump%nam%local_diag) then
   bump%nam%displ_rad = config_get_real(c_conf,"displ_rad")/req
   bump%nam%displ_niter = config_get_int(c_conf,"displ_niter")
   bump%nam%displ_rhflt = config_get_real(c_conf,"displ_rhflt")/req
   bump%nam%displ_tol = config_get_real(c_conf,"displ_tol")
end if

! fit_param
bump%nam%minim_algo = config_get_string(c_conf,1024,"minim_algo")
bump%nam%lhomh = integer_to_logical(config_get_int(c_conf,"lhomh"))
bump%nam%lhomv = integer_to_logical(config_get_int(c_conf,"lhomv"))
bump%nam%rvflt = config_get_real(c_conf,"rvflt")

! nicas_param
bump%nam%lsqrt = integer_to_logical(config_get_int(c_conf,"lsqrt"))
bump%nam%resol = config_get_real(c_conf,"resol")
bump%nam%nicas_interp = config_get_string(c_conf,1024,"nicas_interp")
bump%nam%network = integer_to_logical(config_get_int(c_conf,"network"))
bump%nam%mpicom = config_get_int(c_conf,"mpicom")
bump%nam%advmode = config_get_int(c_conf,"advmode")
bump%nam%ndir = config_get_int(c_conf,"ndir")
do idir=1,bump%nam%ndir
   write(idirchar,'(i3)') idir
   bump%nam%londir(idir) = config_get_real(c_conf,"londir("//trim(adjustl(idirchar))//")")
   bump%nam%latdir(idir) = config_get_real(c_conf,"latdir("//trim(adjustl(idirchar))//")")
   bump%nam%levdir(idir) = config_get_int(c_conf,"levdir("//trim(adjustl(idirchar))//")")
   bump%nam%ivdir(idir) = config_get_int(c_conf,"ivdir("//trim(adjustl(idirchar))//")")
   bump%nam%itsdir(idir) = config_get_int(c_conf,"itsdir("//trim(adjustl(idirchar))//")")
end do

! output_param
bump%nam%nldwh = config_get_int(c_conf,"nldwh")
do ildwh=1,bump%nam%nldwh
   write(ildwvchar,'(i3)') ildwh
   bump%nam%il_ldwh(ildwh) = config_get_int(c_conf,"il_ldwh("//trim(adjustl(ildwhchar))//")")
   bump%nam%ic_ldwh(ildwh) = config_get_int(c_conf,"ic_ldwh("//trim(adjustl(ildwhchar))//")")
end do
bump%nam%nldwv = config_get_int(c_conf,"nldwv")
do ildwv=1,bump%nam%nldwv
   write(ildwvchar,'(i3)') ildwv
   bump%nam%lon_ldwv(ildwv) = config_get_real(c_conf,"lon_ldwv("//trim(adjustl(ildwvchar))//")")
   bump%nam%lat_ldwv(ildwv) = config_get_real(c_conf,"lat_ldwv("//trim(adjustl(ildwvchar))//")")
end do
bump%nam%diag_rhflt = config_get_real(c_conf,"diag_rhflt")/req
bump%nam%diag_interp = config_get_string(c_conf,1024,"diag_interp")
bump%nam%grid_output = integer_to_logical(config_get_int(c_conf,"grid_output"))
if (bump%nam%grid_output) then
   bump%nam%grid_resol = config_get_real(c_conf,"grid_resol")
   bump%nam%grid_interp = config_get_string(c_conf,1024,"grid_interp")
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

subroutine delete_bump(self)
implicit none
type(bump_type), intent(inout) :: self

end subroutine delete_bump

!-------------------------------------------------------------------------------

subroutine bump_multiply(self,ug)
implicit none
type(bump_type), intent(inout) :: self
type(unstructured_grid), intent(inout) :: ug

! Apply localization
call self%apply_nicas(ug%fld)

end subroutine bump_multiply

!-------------------------------------------------------------------------------

end module bump_mod
