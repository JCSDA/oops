! (C) Copyright 2009-2016 ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_change_var_mod

use qg_convert_q_to_x_mod
use qg_convert_x_to_q_mod
use qg_fields_mod
use oops_variables_mod

implicit none

private
public :: qg_change_var_config
public :: qg_change_var_registry
public :: qg_change_var_setup,qg_change_var,qg_change_var_inv,qg_change_var_ad,qg_change_var_inv_ad
! ------------------------------------------------------------------------------
type :: qg_change_var_config
  character(len=1024) :: varchange !< Variable change name
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
!> Setup change of variable
subroutine qg_change_var_setup(self,vars_in,vars_out)

implicit none

! Passed variables
type(qg_change_var_config),intent(inout) :: self !< Variable change
type(oops_variables),intent(in) :: vars_in       !< Input variables
type(oops_variables),intent(in) :: vars_out      !< Output variables

! Check
if ((vars_in%nvars() /= 1) .or. (vars_out%nvars() /= 1)) then
  call abor1_ftn('qg_change_var_setup: wrong change of variable')
endif
if (vars_in%variable(1) == vars_out%variable(1)) then
  self%varchange = 'identity'
elseif ((vars_in%variable(1) == 'x') .and. (vars_out%variable(1) == 'q')) then
  self%varchange = 'x_to_q'
elseif ((vars_in%variable(1) == 'q') .and. (vars_out%variable(1) == 'x')) then
  self%varchange = 'q_to_x'
else
  call abor1_ftn('qg_change_var_setup: wrong change of variable')
endif

end subroutine qg_change_var_setup
! ------------------------------------------------------------------------------
!> Change of variable
subroutine qg_change_var(conf,fld_in,fld_out)

implicit none

! Passed variables
type(qg_change_var_config),intent(in) :: conf !< Variable change
type(qg_fields),intent(in) :: fld_in          !< Input fields
type(qg_fields),intent(inout) :: fld_out      !< Output fields

! Check fields resolution
call qg_fields_check_resolution(fld_in,fld_out)

select case (trim(conf%varchange))
case ('identity')
  ! Copy fields
  call qg_fields_copy(fld_out,fld_in)
  fld_out%lq = fld_in%lq
case ('x_to_q')
  ! Check fields variables
  if (fld_in%lq) call abor1_ftn('qg_change_var: wrong input fields variables for '//trim(conf%varchange))

  ! Conversion
  call convert_x_to_q_tl(fld_in%geom,fld_in%gfld3d,fld_out%gfld3d)
  fld_out%lq = .true.

  ! Copy boundary conditions
  call qg_fields_copy(fld_out,fld_in,.true.)
case ('q_to_x')
  ! Check fields variables
  if (.not.fld_in%lq) call abor1_ftn('qg_change_var: wrong input fields variables for '//trim(conf%varchange))

  ! Conversion
  call convert_q_to_x_tl(fld_in%geom,fld_in%gfld3d,fld_out%gfld3d)
  fld_out%lq = .false.

  ! Copy boundary conditions
  call qg_fields_copy(fld_out,fld_in,.true.)
case default
  call abor1_ftn('qg_change_var: wrong variable change')
end select

end subroutine qg_change_var
! ------------------------------------------------------------------------------
!> Change of variable - inverse
subroutine qg_change_var_inv(conf,fld_in,fld_out)

implicit none

! Passed variables
type(qg_change_var_config),intent(in) :: conf !< Variable change
type(qg_fields),intent(in) :: fld_in          !< Input fields
type(qg_fields),intent(inout) :: fld_out      !< Output fields

! Check fields resolution
call qg_fields_check_resolution(fld_in,fld_out)

select case (trim(conf%varchange))
case ('identity')
  ! Copy fields
  call qg_fields_copy(fld_out,fld_in)
  fld_out%lq = fld_in%lq
case ('x_to_q')
  ! Check fields variables
  if (.not.fld_in%lq) call abor1_ftn('qg_change_var_inv: wrong input fields variables for '//trim(conf%varchange))

  ! Conversion
  call convert_q_to_x_tl(fld_in%geom,fld_in%gfld3d,fld_out%gfld3d)
  fld_out%lq = .false.

  ! Copy boundary conditions
  call qg_fields_copy(fld_out,fld_in,.true.)
case ('q_to_x')
  ! Check fields variables
  if (fld_in%lq) call abor1_ftn('qg_change_var_inv: wrong input fields variables for '//trim(conf%varchange))

  ! Conversion
  call convert_x_to_q_tl(fld_in%geom,fld_in%gfld3d,fld_out%gfld3d)
  fld_out%lq = .true.

  ! Copy boundary conditions
  call qg_fields_copy(fld_out,fld_in,.true.)
case default
  call abor1_ftn('qg_change_var_inv: wrong variable change')
end select

end subroutine qg_change_var_inv
! ------------------------------------------------------------------------------
!> Change of variable - adjoint
subroutine qg_change_var_ad(conf,fld_in,fld_out)

implicit none

! Passed variables
type(qg_change_var_config),intent(in) :: conf !< Variable change
type(qg_fields),intent(in) :: fld_in          !< Input fields
type(qg_fields),intent(inout) :: fld_out      !< Output fields

! Check fields resolution
call qg_fields_check_resolution(fld_in,fld_out)

select case (trim(conf%varchange))
case ('identity')
  ! Copy fields
  call qg_fields_copy(fld_out,fld_in)
  fld_out%lq = fld_in%lq
case ('x_to_q')
  ! Check fields variables
  if (.not.fld_in%lq) call abor1_ftn('qg_change_var_ad: wrong input fields variables for '//trim(conf%varchange))

  ! Conversion
  call convert_x_to_q_ad(fld_in%geom,fld_in%gfld3d,fld_out%gfld3d)
  fld_out%lq = .false.

  ! Copy boundary conditions
  call qg_fields_copy(fld_out,fld_in,.true.)
case ('q_to_x')
  ! Check fields variables
  if (fld_in%lq) call abor1_ftn('qg_change_var_ad: wrong input fields variables for '//trim(conf%varchange))

  ! Conversion
  call convert_q_to_x_ad(fld_in%geom,fld_in%gfld3d,fld_out%gfld3d)
  fld_out%lq = .true.

  ! Copy boundary conditions
  call qg_fields_copy(fld_out,fld_in,.true.)
case default
  call abor1_ftn('qg_change_var_ad: wrong variable change')
end select

end subroutine qg_change_var_ad
! ------------------------------------------------------------------------------
!> Change of variable - inverse adjoint
subroutine qg_change_var_inv_ad(conf,fld_in,fld_out)

implicit none

! Passed variables
type(qg_change_var_config),intent(in) :: conf !< Variable change
type(qg_fields),intent(in) :: fld_in          !< Input fields
type(qg_fields),intent(inout) :: fld_out      !< Output fields

! Check fields resolution
call qg_fields_check_resolution(fld_in,fld_out)

select case (trim(conf%varchange))
case ('identity')
  ! Copy fields
  call qg_fields_copy(fld_out,fld_in)
  fld_out%lq = fld_in%lq
case ('x_to_q')
  ! Check fields variables
  if (fld_in%lq) call abor1_ftn('qg_change_var_inv_ad: wrong input fields variables for '//trim(conf%varchange))

  ! Conversion
  call convert_q_to_x_ad(fld_in%geom,fld_in%gfld3d,fld_out%gfld3d)
  fld_out%lq = .true.

  ! Copy boundary conditions
  call qg_fields_copy(fld_out,fld_in,.true.)
case ('q_to_x')
  ! Check fields variables
  if (.not.fld_in%lq) call abor1_ftn('qg_change_var_inv_ad: wrong input fields variables for '//trim(conf%varchange))

  ! Conversion
  call convert_x_to_q_ad(fld_in%geom,fld_in%gfld3d,fld_out%gfld3d)
  fld_out%lq = .false.

  ! Copy boundary conditions
  call qg_fields_copy(fld_out,fld_in,.true.)
case default
  call abor1_ftn('qg_change_var_inv_ad: wrong variable change')
end select

end subroutine qg_change_var_inv_ad
! ------------------------------------------------------------------------------
end module qg_change_var_mod
