! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_change_var_mod

use kinds
use qg_convert_q_to_x_mod
use qg_convert_x_to_q_mod
use qg_convert_x_to_u_mod
use qg_convert_x_to_v_mod
use qg_fields_mod
use oops_variables_mod

implicit none

private
public :: qg_change_var_registry
public :: qg_change_var
public :: qg_change_var_tl,qg_change_var_ad
! ------------------------------------------------------------------------------
type :: qg_change_var_config
  ! Input variables flags
  logical :: x_in
  logical :: q_in
  logical :: u_in
  logical :: v_in

  ! Output variables flags
  logical :: x_out
  logical :: q_out
  logical :: u_out
  logical :: v_out

  ! Conversions
  logical :: x_to_q
  logical :: q_to_x
  logical :: x_to_u
  logical :: x_to_v
end type qg_change_var_config

#define LISTED_TYPE qg_change_var_config

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_change_var_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------
!> Setup change of variables setup
subroutine qg_change_var_setup(conf,fld_in,fld_out,ad)

implicit none

! Passed variables
type(qg_change_var_config),intent(inout) :: conf !< Variable change
type(qg_fields),intent(in) :: fld_in             !< Input field
type(qg_fields),intent(in) :: fld_out            !< Output field
logical,intent(in),optional :: ad                !< Adjoint flag

! Local variables
logical :: lad

! Local flag
lad = .false.
if (present(ad)) lad = ad

! Check what is allocated in fields
if (lad) then
  ! Adjoint case
  conf%x_in = allocated(fld_out%x)
  conf%q_in = allocated(fld_out%q)
  conf%u_in = allocated(fld_out%u)
  conf%v_in = allocated(fld_out%v)
  conf%x_out = allocated(fld_in%x)
  conf%q_out = allocated(fld_in%q)
  conf%u_out = allocated(fld_in%u)
  conf%v_out = allocated(fld_in%v)
else
  ! Normal case
  conf%x_in = allocated(fld_in%x)
  conf%q_in = allocated(fld_in%q)
  conf%u_in = allocated(fld_in%u)
  conf%v_in = allocated(fld_in%v)
  conf%x_out = allocated(fld_out%x)
  conf%q_out = allocated(fld_out%q)
  conf%u_out = allocated(fld_out%u)
  conf%v_out = allocated(fld_out%v)
endif

! Check input/output consistency
if (conf%x_out.and.(.not.(conf%x_in.or.conf%q_in))) call abor1_ftn('qg_change_var_setup: x or q required to compute x')
if (conf%q_out.and.(.not.(conf%x_in.or.conf%q_in))) call abor1_ftn('qg_change_var_setup: x or q required to compute q')
if (conf%u_out.and.(.not.(conf%x_in.or.conf%q_in.or.conf%u_in))) &
 & call abor1_ftn('qg_change_var_setup: x, q or u required to compute u')
if (conf%v_out.and.(.not.(conf%x_in.or.conf%q_in.or.conf%v_in))) &
 & call abor1_ftn('qg_change_var_setup: x, q or v required to compute v')

! Initialize required conversions
conf%q_to_x = .false.
conf%x_to_q = .false.
conf%x_to_u = .false.
conf%x_to_v = .false.

! Define required conversions
if (conf%x_out) conf%q_to_x = (.not.conf%x_in)
if (conf%q_out) conf%x_to_q = (.not.conf%q_in)
if (conf%u_out) then
  if (.not.conf%u_in) then
    conf%x_to_u = .true.
    conf%q_to_x = (.not.conf%x_in)
  endif
endif
if (conf%v_out) then
  if (.not.conf%v_in) then
    conf%x_to_v = .true.
    conf%q_to_x = (.not.conf%x_in)
  endif
endif

end subroutine qg_change_var_setup
! ------------------------------------------------------------------------------
!> Get variables
subroutine qg_change_var_get(conf,fld,x,q,u,v)

implicit none

! Passed variables
type(qg_change_var_config),intent(in) :: conf                         !< Variable change
type(qg_fields),intent(in) :: fld                                     !< Fields
real(kind_real),intent(out) :: x(fld%geom%nx,fld%geom%ny,fld%geom%nz) !< Streamfunction
real(kind_real),intent(out) :: q(fld%geom%nx,fld%geom%ny,fld%geom%nz) !< Potential vorticity
real(kind_real),intent(out) :: u(fld%geom%nx,fld%geom%ny,fld%geom%nz) !< Zonal wind
real(kind_real),intent(out) :: v(fld%geom%nx,fld%geom%ny,fld%geom%nz) !< Meridional wind

if (allocated(fld%x)) then
  x = fld%x
else
  x = 0.0_kind_real
endif
if (allocated(fld%q)) then
  q = fld%q
else
  q = 0.0_kind_real
endif
if (allocated(fld%u)) then
  u = fld%u
else
  u = 0.0_kind_real
endif
if (allocated(fld%v)) then
  v = fld%v
else
  v = 0.0_kind_real
endif

end subroutine qg_change_var_get
! ------------------------------------------------------------------------------
!> Set variables
subroutine qg_change_var_set(conf,fld,x,q,u,v)

implicit none

! Passed variables
type(qg_change_var_config),intent(in) :: conf                        !< Variable change
type(qg_fields),intent(inout) :: fld                                 !< Fields
real(kind_real),intent(in) :: x(fld%geom%nx,fld%geom%ny,fld%geom%nz) !< Streamfunction
real(kind_real),intent(in) :: q(fld%geom%nx,fld%geom%ny,fld%geom%nz) !< Potential vorticity
real(kind_real),intent(in) :: u(fld%geom%nx,fld%geom%ny,fld%geom%nz) !< Zonal wind
real(kind_real),intent(in) :: v(fld%geom%nx,fld%geom%ny,fld%geom%nz) !< Meridional wind

if (allocated(fld%x)) fld%x = x
if (allocated(fld%q)) fld%q = q
if (allocated(fld%u)) fld%u = u
if (allocated(fld%v)) fld%v = v

end subroutine qg_change_var_set
! ------------------------------------------------------------------------------
!> Change of variable
subroutine qg_change_var(fld_in,fld_out)

implicit none

! Passed variables
type(qg_fields),intent(in) :: fld_in     !< Input fields
type(qg_fields),intent(inout) :: fld_out !< Output fields

! Local variables
real(kind_real) :: x(fld_in%geom%nx,fld_in%geom%ny,fld_in%geom%nz),q(fld_in%geom%nx,fld_in%geom%ny,fld_in%geom%nz)
real(kind_real) :: u(fld_in%geom%nx,fld_in%geom%ny,fld_in%geom%nz),v(fld_in%geom%nx,fld_in%geom%ny,fld_in%geom%nz)
type(qg_change_var_config) :: conf

! Check resolution
call qg_fields_check_resolution(fld_in,fld_out)

! Copy boundary conditions
call qg_fields_copy_lbc(fld_out,fld_in)

! Define change of variable configuration
call qg_change_var_setup(conf,fld_in,fld_out)

! Get x, q, u and v
call qg_change_var_get(conf,fld_in,x,q,u,v)

! Conversions
if (conf%x_to_q) call convert_x_to_q(fld_in%geom,x,fld_in%x_north,fld_in%x_south,q)
if (conf%q_to_x) call convert_q_to_x(fld_in%geom,q,fld_in%x_north,fld_in%x_south,x)
if (conf%x_to_u) call convert_x_to_u(fld_in%geom,x,fld_in%x_north,fld_in%x_south,u)
if (conf%x_to_v) call convert_x_to_v(fld_in%geom,x,v)

! Set x, q, u and v
call qg_change_var_set(conf,fld_out,x,q,u,v)

end subroutine qg_change_var
! ------------------------------------------------------------------------------
!> Change of variable
subroutine qg_change_var_tl(fld_in,fld_out)

implicit none

! Passed variables
type(qg_fields),intent(in) :: fld_in     !< Input fields
type(qg_fields),intent(inout) :: fld_out !< Output fields

! Local variables
real(kind_real) :: x(fld_in%geom%nx,fld_in%geom%ny,fld_in%geom%nz),q(fld_in%geom%nx,fld_in%geom%ny,fld_in%geom%nz)
real(kind_real) :: u(fld_in%geom%nx,fld_in%geom%ny,fld_in%geom%nz),v(fld_in%geom%nx,fld_in%geom%ny,fld_in%geom%nz)
type(qg_change_var_config) :: conf

! Check resolution
call qg_fields_check_resolution(fld_in,fld_out)

! Copy boundary conditions
call qg_fields_copy_lbc(fld_out,fld_in)

! Define change of variable configuration
call qg_change_var_setup(conf,fld_in,fld_out)

! Get x, q, u and v
call qg_change_var_get(conf,fld_in,x,q,u,v)

! Conversions
if (conf%x_to_q) call convert_x_to_q_tl(fld_in%geom,x,q)
if (conf%q_to_x) call convert_q_to_x_tl(fld_in%geom,q,x)
if (conf%x_to_u) call convert_x_to_u_tl(fld_in%geom,x,u)
if (conf%x_to_v) call convert_x_to_v_tl(fld_in%geom,x,v)

! Set x, q, u and v
call qg_change_var_set(conf,fld_out,x,q,u,v)

end subroutine qg_change_var_tl
! ------------------------------------------------------------------------------
!> Change of variable - adjoint
subroutine qg_change_var_ad(fld_in,fld_out)

implicit none

! Passed variables
type(qg_fields),intent(in) :: fld_in     !< Input fields
type(qg_fields),intent(inout) :: fld_out !< Output fields

! Local variables
real(kind_real) :: x(fld_in%geom%nx,fld_in%geom%ny,fld_in%geom%nz),q(fld_in%geom%nx,fld_in%geom%ny,fld_in%geom%nz)
real(kind_real) :: u(fld_in%geom%nx,fld_in%geom%ny,fld_in%geom%nz),v(fld_in%geom%nx,fld_in%geom%ny,fld_in%geom%nz)
type(qg_change_var_config) :: conf

! Checks
call qg_fields_check_resolution(fld_in,fld_out)

! Copy boundary conditions
call qg_fields_copy_lbc(fld_out,fld_in)

! Define change of variable configuration
call qg_change_var_setup(conf,fld_in,fld_out,.true.)

! Get x, q, u and v
call qg_change_var_get(conf,fld_in,x,q,u,v)

! Conversions
if (conf%x_to_v) call convert_x_to_v_ad(fld_in%geom,v,x)
if (conf%x_to_u) call convert_x_to_u_ad(fld_in%geom,u,x)
if (conf%q_to_x) call convert_q_to_x_ad(fld_in%geom,x,q)
if (conf%x_to_q) call convert_x_to_q_ad(fld_in%geom,q,x)

! Set x, q, u and v
call qg_change_var_set(conf,fld_out,x,q,u,v)

end subroutine qg_change_var_ad
! ------------------------------------------------------------------------------
end module qg_change_var_mod
