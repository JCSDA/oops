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
!> Setup change of variables
subroutine qg_change_var_setup(conf,fld,vars,ad)

implicit none

! Passed variables
type(qg_change_var_config),intent(inout) :: conf !< Variable change
type(qg_fields),intent(inout) :: fld             !< field
type(oops_variables),intent(in) :: vars          !< Output variables
logical,intent(in),optional :: ad                !< Adjoint flag

! Local variables
logical :: lad

! Local flag
lad = .false.
if (present(ad)) lad = ad

if (lad) then
  ! Adjoint case
  conf%x_in = vars%has('x')
  conf%q_in = vars%has('q')
  conf%u_in = vars%has('u')
  conf%v_in = vars%has('v')
  conf%x_out = fld%vars%has('x')
  conf%q_out = fld%vars%has('q')
  conf%u_out = fld%vars%has('u')
  conf%v_out = fld%vars%has('v')
else
  ! Normal case
  conf%x_in = fld%vars%has('x')
  conf%q_in = fld%vars%has('q')
  conf%u_in = fld%vars%has('u')
  conf%v_in = fld%vars%has('v')
  conf%x_out = vars%has('x')
  conf%q_out = vars%has('q')
  conf%u_out = vars%has('u')
  conf%v_out = vars%has('v')
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

! Allocation and initialization
if (lad) then
  ! Adjoint case
  if (conf%q_to_x) then
    if (.not. allocated(fld%q)) then
        allocate(fld%q(fld%geom%nx,fld%geom%ny,fld%geom%nz))
        fld%q = 0.0_kind_real
    endif
  endif
  if (conf%x_to_q.or.conf%x_to_v.or.conf%x_to_u) then
    if (.not. allocated(fld%x)) then
        allocate(fld%x(fld%geom%nx,fld%geom%ny,fld%geom%nz))
        fld%x = 0.0_kind_real
    endif
  endif
else
  ! Normal case
  if (conf%q_to_x) then
    if (.not. allocated(fld%x)) allocate(fld%x(fld%geom%nx,fld%geom%ny,fld%geom%nz))
    fld%x = 0.0_kind_real
  endif
  if (conf%x_to_q) then
    if (.not. allocated(fld%q)) allocate(fld%q(fld%geom%nx,fld%geom%ny,fld%geom%nz))
    fld%q = 0.0_kind_real
  endif
  if (conf%x_to_v) then
    if (.not. allocated(fld%v)) allocate(fld%v(fld%geom%nx,fld%geom%ny,fld%geom%nz))
    fld%v = 0.0_kind_real
  endif
  if (conf%x_to_u) then
    if (.not. allocated(fld%u)) allocate(fld%u(fld%geom%nx,fld%geom%ny,fld%geom%nz))
    fld%u = 0.0_kind_real
  endif
endif

end subroutine qg_change_var_setup

! ------------------------------------------------------------------------------
!> Free memory and reset fld%vars
subroutine qg_change_var_cleanup(fld,vars)

implicit none

! Passed variables
type(qg_fields),intent(inout) :: fld             !< field
type(oops_variables),intent(in) :: vars          !< Output variables

! Release memory
if (fld%vars%has('x') .and. .not.vars%has('x')) deallocate(fld%x)
if (fld%vars%has('q') .and. .not.vars%has('q')) deallocate(fld%q)
if (fld%vars%has('u') .and. .not.vars%has('u')) deallocate(fld%u)
if (fld%vars%has('v') .and. .not.vars%has('v')) deallocate(fld%v)

! Reset variables
call fld%vars%destruct()
fld%vars = oops_variables(vars)

end subroutine qg_change_var_cleanup

! ------------------------------------------------------------------------------
!> Change of variable
subroutine qg_change_var(fld,vars)
implicit none
type(qg_fields),intent(inout) :: fld     !< Fields to be transformed
type(oops_variables),intent(in) :: vars  !< Output variables

type(qg_change_var_config) :: conf

! Define change of variable configuration
call qg_change_var_setup(conf, fld, vars)

! Conversions
if (conf%x_to_q) call convert_x_to_q(fld%geom, fld%x, fld%x_north, fld%x_south, fld%q)
if (conf%q_to_x) call convert_q_to_x(fld%geom, fld%q, fld%x_north, fld%x_south, fld%x)
if (conf%x_to_u) call convert_x_to_u(fld%geom, fld%x, fld%x_north, fld%x_south, fld%u)
if (conf%x_to_v) call convert_x_to_v(fld%geom, fld%x, fld%v)

call qg_change_var_cleanup(fld, vars)

end subroutine qg_change_var

! ------------------------------------------------------------------------------
!> Change of variable
subroutine qg_change_var_tl(fld,vars)
implicit none
type(qg_fields),intent(inout) :: fld     !< Fields to be transformed
type(oops_variables),intent(in) :: vars  !< Output variables

type(qg_change_var_config) :: conf

! Define change of variable configuration
call qg_change_var_setup(conf, fld, vars)

! Conversions
if (conf%x_to_q) call convert_x_to_q_tl(fld%geom, fld%x, fld%q)
if (conf%q_to_x) call convert_q_to_x_tl(fld%geom, fld%q, fld%x)
if (conf%x_to_u) call convert_x_to_u_tl(fld%geom, fld%x, fld%u)
if (conf%x_to_v) call convert_x_to_v_tl(fld%geom, fld%x, fld%v)

call qg_change_var_cleanup(fld, vars)

end subroutine qg_change_var_tl
! ------------------------------------------------------------------------------
!> Change of variable - adjoint
subroutine qg_change_var_ad(fld,vars)
implicit none
type(qg_fields),intent(inout) :: fld     !< Fields to be transformed
type(oops_variables),intent(in) :: vars  !< Output variables

type(qg_change_var_config) :: conf

! Define change of variable configuration
call qg_change_var_setup(conf, fld, vars, .true.)

! Conversions
if (conf%x_to_v) call convert_x_to_v_ad(fld%geom, fld%v, fld%x)
if (conf%x_to_u) call convert_x_to_u_ad(fld%geom, fld%u, fld%x)
if (conf%q_to_x) call convert_q_to_x_ad(fld%geom, fld%x, fld%q)
if (conf%x_to_q) call convert_x_to_q_ad(fld%geom, fld%q, fld%x)

call qg_change_var_cleanup(fld, vars)

end subroutine qg_change_var_ad
! ------------------------------------------------------------------------------
end module qg_change_var_mod