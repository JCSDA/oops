! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_model_mod

use duration_mod
use fckit_configuration_module, only: fckit_configuration
use fckit_log_module,only: fckit_log
use iso_c_binding
use kinds
use qg_advect_q_mod
use qg_constants_mod
use qg_convert_q_to_x_mod
use qg_convert_x_to_q_mod
use qg_convert_x_to_u_mod
use qg_convert_x_to_v_mod
use qg_fields_mod
use random_mod

implicit none

private
public :: qg_model_config
public :: qg_model_registry
public :: qg_model_setup,qg_model_propagate,qg_model_propagate_tl,qg_model_propagate_ad
! ------------------------------------------------------------------------------
type :: qg_model_config
  real(kind_real) :: dt !< Time step (seconds)
end type qg_model_config

#define LISTED_TYPE qg_model_config

!> Linked list interface - defines registry_t type
#include "oops/util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_model_registry
! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
! Public
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "oops/util/linkedList_c.f"
! ------------------------------------------------------------------------------
!> Setup model
subroutine qg_model_setup(self,f_conf)

implicit none

! Passed variables
type(qg_model_config),intent(inout) :: self    !< Model configuration
type(fckit_configuration),intent(in) :: f_conf !< FCKIT configuration

! Local variables
type(duration) :: dtstep
character(len=20) :: ststep
character(len=160) :: record
character(len=:),allocatable :: str

! Define time step
call f_conf%get_or_die("tstep",str)
ststep = str
dtstep = trim(ststep)
self%dt = duration_seconds(dtstep)
write(record,*) 'qg_model_setup: dt = ',self%dt
call fckit_log%info(record)

end subroutine qg_model_setup
! ------------------------------------------------------------------------------
!> Perform a timestep of the QG model
subroutine qg_model_propagate(conf,fld)

implicit none

! Passed variables
type(qg_model_config),intent(in) :: conf !< Model configuration
type(qg_fields),intent(inout) :: fld     !< State fields

! Local variables
real(kind_real) :: q(fld%geom%nx,fld%geom%ny,fld%geom%nz),qnew(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: u(fld%geom%nx,fld%geom%ny,fld%geom%nz),v(fld%geom%nx,fld%geom%ny,fld%geom%nz)

! Check input
if (.not.allocated(fld%x)) call abor1_ftn('qg_model_propagate: x required')
if (.not.allocated(fld%x_north)) call abor1_ftn('qg_model_propagate: x_north required')
if (.not.allocated(fld%x_south)) call abor1_ftn('qg_model_propagate: x_south required')
if (.not.allocated(fld%q_north)) call abor1_ftn('qg_model_propagate: q_north required')
if (.not.allocated(fld%q_south)) call abor1_ftn('qg_model_propagate: q_south required')

! Compute potential vorticity
call convert_x_to_q(fld%geom,fld%x,fld%x_north,fld%x_south,q)

! Compute wind
call convert_x_to_u(fld%geom,fld%x,fld%x_north,fld%x_south,u)
call convert_x_to_v(fld%geom,fld%x,v)

! Advect potential vorticity
call advect_q(fld%geom,conf%dt,u,v,q,fld%q_north,fld%q_south,qnew)

! Compute streamfunction
call convert_q_to_x(fld%geom,qnew,fld%x_north,fld%x_south,fld%x)

! Complete other fields
call qg_fields_complete(fld,'x')

end subroutine qg_model_propagate
! ------------------------------------------------------------------------------
!> Perform a timestep of the QG model - tangent linear
subroutine qg_model_propagate_tl(conf,traj,fld)

implicit none

! Passed variables
type(qg_model_config),intent(in)  :: conf !< Model configuration
type(qg_fields),intent(in) :: traj        !< Trajectory fields
type(qg_fields),intent(inout) :: fld      !< Increment fields

! Local variables
real(kind_real) :: q_traj(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: u_traj(fld%geom%nx,fld%geom%ny,fld%geom%nz),v_traj(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: q(fld%geom%nx,fld%geom%ny,fld%geom%nz),qnew(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: u(fld%geom%nx,fld%geom%ny,fld%geom%nz),v(fld%geom%nx,fld%geom%ny,fld%geom%nz)

! Trajectory

! Check input
if (.not.allocated(traj%x)) call abor1_ftn('qg_model_propagate_tl: x trajectory required')
if (.not.allocated(traj%x_north)) call abor1_ftn('qg_model_propagate_tl: x_north required')
if (.not.allocated(traj%x_south)) call abor1_ftn('qg_model_propagate_tl: x_south required')
if (.not.allocated(traj%q_north)) call abor1_ftn('qg_model_propagate_tl: q_north required')
if (.not.allocated(traj%q_south)) call abor1_ftn('qg_model_propagate_tl: q_south required')

! Compute potential vorticity
call convert_x_to_q(traj%geom,traj%x,traj%x_north,traj%x_south,q_traj)

! Compute wind
call convert_x_to_u(traj%geom,traj%x,traj%x_north,traj%x_south,u_traj)
call convert_x_to_v(traj%geom,traj%x,v_traj)

! Perturbation

! Check input
if (.not.allocated(fld%x)) call abor1_ftn('qg_model_propagate_tl: x perturbation required')

! Compute potential vorticity
call convert_x_to_q_tl(fld%geom,fld%x,q)

! Compute wind
call convert_x_to_u_tl(fld%geom,fld%x,u)
call convert_x_to_v_tl(fld%geom,fld%x,v)

! Advect potential vorticity
call advect_q_tl(fld%geom,conf%dt,u_traj,v_traj,q_traj,traj%q_north,traj%q_south,u,v,q,qnew)

! Compute streamfunction
call convert_q_to_x_tl(fld%geom,qnew,fld%x)

! Complete other fields
call qg_fields_complete(fld,'x')

end subroutine qg_model_propagate_tl
! ------------------------------------------------------------------------------
!> Perform a timestep of the QG model - adjoint
subroutine qg_model_propagate_ad(conf,traj,fld)

implicit none

! Passed variables
type(qg_model_config),intent(in)  :: conf !< Model configuration
type(qg_fields),intent(in) :: traj        !< Trajectory fields
type(qg_fields),intent(inout) :: fld      !< Increment fields

! Local variables
real(kind_real) :: q_traj(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: u_traj(fld%geom%nx,fld%geom%ny,fld%geom%nz),v_traj(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: q(fld%geom%nx,fld%geom%ny,fld%geom%nz),qnew(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: u(fld%geom%nx,fld%geom%ny,fld%geom%nz),v(fld%geom%nx,fld%geom%ny,fld%geom%nz)

! Trajectory

! Check input
if (.not.allocated(traj%x)) call abor1_ftn('qg_model_propagate_tl: x trajectory required')
if (.not.allocated(traj%x_north)) call abor1_ftn('qg_model_propagate_tl: x_north required')
if (.not.allocated(traj%x_south)) call abor1_ftn('qg_model_propagate_tl: x_south required')
if (.not.allocated(traj%q_north)) call abor1_ftn('qg_model_propagate_tl: q_north required')
if (.not.allocated(traj%q_south)) call abor1_ftn('qg_model_propagate_tl: q_south required')

! Compute potential vorticity
call convert_x_to_q(traj%geom,traj%x,traj%x_north,traj%x_south,q_traj)

! Compute wind
call convert_x_to_u(traj%geom,traj%x,traj%x_north,traj%x_south,u_traj)
call convert_x_to_v(traj%geom,traj%x,v_traj)

! Perturbation

! Check input
if (.not.allocated(fld%x)) call abor1_ftn('qg_model_propagate_tl: x perturbation required')

! Initialization
u = 0.0_kind_real
v = 0.0_kind_real
q = 0.0_kind_real
qnew = 0.0_kind_real

! Compute streamfunction
call convert_q_to_x_ad(fld%geom,fld%x,qnew)

! Advect potential vorticity
call advect_q_ad(fld%geom,conf%dt,u_traj,v_traj,q_traj,traj%q_north,traj%q_south,qnew,u,v,q)

! Initialize x
fld%x = 0.0_kind_real

! Compute wind
call convert_x_to_v_ad(fld%geom,v,fld%x)
call convert_x_to_u_ad(fld%geom,u,fld%x)

! Compute potential vorticity
call convert_x_to_q_ad(fld%geom,q,fld%x)

! Complete other fields
call qg_fields_complete(fld,'x')

end subroutine qg_model_propagate_ad
! ------------------------------------------------------------------------------
end module qg_model_mod

