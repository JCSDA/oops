! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

module qg_model_mod

use config_mod
use duration_mod
use fckit_log_module,only: fckit_log
use iso_c_binding
use kinds
use qg_advect_q_mod
use qg_constants_mod
use qg_convert_q_to_x_mod
use qg_convert_x_to_q_mod
use qg_convert_x_to_uv_mod
use qg_fields_mod
use random_mod

implicit none

private
public :: qg_model_config
public :: qg_model_registry
public :: qg_model_setup,qg_model_propagate,qg_model_propagate_tl,qg_model_propagate_ad,qg_model_propagate_test
! ------------------------------------------------------------------------------
type :: qg_model_config
  real(kind_real) :: dt !< Time step (seconds)
end type qg_model_config

real(kind_real),parameter :: eps = 1.0e-10 !< Epsilon value for adjoint tests

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
subroutine qg_model_setup(self,conf)

implicit none

! Passed variables
type(qg_model_config),intent(inout) :: self !< Model configuration
type(c_ptr),intent(in) :: conf              !< Configuration

! Local variables
type(duration) :: dtstep
character(len=20) :: ststep
character(len=160) :: record

! Define time step
ststep = config_get_string(conf,len(ststep),'tstep')
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
real(kind_real) :: x(fld%geom%nx,fld%geom%ny,fld%geom%nz),q(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: u(fld%geom%nx,fld%geom%ny,fld%geom%nz),v(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: qnew(fld%geom%nx,fld%geom%ny,fld%geom%nz)

! Initialize streamfunction and potential vorticity
if (fld%lq) then
  q = fld%gfld3d
  call convert_q_to_x(fld%geom,q,fld%x_north,fld%x_south,x)
else
  x = fld%gfld3d
  call convert_x_to_q(fld%geom,x,fld%x_north,fld%x_south,q)
endif

! Initialize wind
call convert_x_to_uv(fld%geom,x,fld%x_north,fld%x_south,u,v)

! Advect potential vorticity
call advect_q(fld%geom,conf%dt,u,v,q,fld%q_north,fld%q_south,qnew)

! Save streamfunction or potential vorticity
if (fld%lq) then
  fld%gfld3d = qnew
else
  call convert_q_to_x(fld%geom,qnew,fld%x_north,fld%x_south,x)
  fld%gfld3d = x
endif

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
real(kind_real) :: x_traj(fld%geom%nx,fld%geom%ny,fld%geom%nz),q_traj(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: u_traj(fld%geom%nx,fld%geom%ny,fld%geom%nz),v_traj(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: x(fld%geom%nx,fld%geom%ny,fld%geom%nz),q(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: u(fld%geom%nx,fld%geom%ny,fld%geom%nz),v(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: qnew(fld%geom%nx,fld%geom%ny,fld%geom%nz)

! Trajectory

! Initialize streamfunction and potential vorticity
if (traj%lq) then
  q_traj = traj%gfld3d
  call convert_q_to_x(fld%geom,q_traj,traj%x_north,traj%x_south,x_traj)
else
  x_traj = traj%gfld3d
  call convert_x_to_q(fld%geom,x_traj,traj%x_north,traj%x_south,q_traj)
endif

! Compute wind
call convert_x_to_uv(fld%geom,x_traj,traj%x_north,traj%x_south,u_traj,v_traj)

! Perturbation

! Initialize streamfunction and potential vorticity
if (fld%lq) then
  q = fld%gfld3d
  call convert_q_to_x_tl(fld%geom,q,x)
else
  x = fld%gfld3d
  call convert_x_to_q_tl(fld%geom,x,q)
endif

! Compute wind
call convert_x_to_uv_tl(fld%geom,x,u,v)

! Advect PV
call advect_q_tl(fld%geom,conf%dt,u_traj,v_traj,q_traj,traj%q_north,traj%q_south,u,v,q,qnew)

! Save streamfunction or potential vorticity
if (fld%lq) then
  fld%gfld3d = qnew
else
  call convert_q_to_x_tl(fld%geom,qnew,x)
  fld%gfld3d = x
endif

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
real(kind_real) :: x_traj(fld%geom%nx,fld%geom%ny,fld%geom%nz),q_traj(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: u_traj(fld%geom%nx,fld%geom%ny,fld%geom%nz),v_traj(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: x(fld%geom%nx,fld%geom%ny,fld%geom%nz),q(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: u(fld%geom%nx,fld%geom%ny,fld%geom%nz),v(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: qnew(fld%geom%nx,fld%geom%ny,fld%geom%nz)

! Trajectory

! Initialize streamfunction and potential vorticity
if (traj%lq) then
  q_traj = traj%gfld3d
  call convert_q_to_x(fld%geom,q_traj,traj%x_north,traj%x_south,x_traj)
else
  x_traj = traj%gfld3d
  call convert_x_to_q(fld%geom,x_traj,traj%x_north,traj%x_south,q_traj)
endif

! Compute wind
call convert_x_to_uv(fld%geom,x_traj,traj%x_north,traj%x_south,u_traj,v_traj)

! Perturbation

! Initialization
x = 0.0
q = 0.0
u = 0.0
v = 0.0
qnew = 0.0

! Save streamfunction or potential vorticity
if (fld%lq) then
  qnew = fld%gfld3d
else
  x = fld%gfld3d
  call convert_q_to_x_ad(fld%geom,x,qnew)
endif

! Advect PV
call advect_q_ad(fld%geom,conf%dt,u_traj,v_traj,q_traj,traj%q_north,traj%q_south,qnew,u,v,q)

! Compute wind
x = 0.0
call convert_x_to_uv_ad(fld%geom,u,v,x)

! Initialize streamfunction and potential vorticity
if (fld%lq) then
  call convert_q_to_x_ad(fld%geom,x,q)
  fld%gfld3d = q
else
  call convert_x_to_q_ad(fld%geom,q,x)
  fld%gfld3d = x
endif

end subroutine qg_model_propagate_ad
! ------------------------------------------------------------------------------
!> Test propagation adjoint
subroutine qg_model_propagate_test(conf,traj,fld)

implicit none

! Passed variables
type(qg_model_config),intent(in)  :: conf !< Model configuration
type(qg_fields),intent(in) :: traj        !< Trajectory fields
type(qg_fields),intent(inout) :: fld      !< Increment fields

! Local variables
real(kind_real) :: x_traj(fld%geom%nx,fld%geom%ny,fld%geom%nz),q_traj(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: u_traj(fld%geom%nx,fld%geom%ny,fld%geom%nz),v_traj(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: x(fld%geom%nx,fld%geom%ny,fld%geom%nz),q(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: u(fld%geom%nx,fld%geom%ny,fld%geom%nz),v(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: qnew(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: x_save(fld%geom%nx,fld%geom%ny,fld%geom%nz),q_save(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: u_save(fld%geom%nx,fld%geom%ny,fld%geom%nz),v_save(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: qnew_save(fld%geom%nx,fld%geom%ny,fld%geom%nz)
real(kind_real) :: dot1,dot2
type(qg_fields) :: fld1,fld2,fld1_save,fld2_save
character(len=1024) :: record

! Trajectory

! Initialize streamfunction and potential vorticity
if (traj%lq) then
  q_traj = traj%gfld3d
  call convert_q_to_x(fld%geom,q_traj,traj%x_north,traj%x_south,x_traj)
else
  x_traj = traj%gfld3d
  call convert_x_to_q(fld%geom,x_traj,traj%x_north,traj%x_south,q_traj)
endif

! Compute wind
call convert_x_to_uv(fld%geom,x_traj,traj%x_north,traj%x_south,u_traj,v_traj)

! Perturbation

! Test q_to_x
call normal_distribution(x_save,0.0_kind_real,1.0_kind_real)
call normal_distribution(q_save,0.0_kind_real,1.0_kind_real)
x = x_save
q = 0.0
call convert_q_to_x_ad(fld%geom,x,q)
call convert_q_to_x_tl(fld%geom,q_save,x)
dot1 = sum(x*x_save)
dot2 = sum(q*q_save)
if (abs(dot1-dot2)/abs(dot1)>eps) then
  write(record,*) 'qg_model_propagate_test: convert_q_to_x adjoint test failed:',dot1,dot2,abs(dot1-dot2)/abs(dot1)
  call abor1_ftn(record)
endif

! Test x_to_q
call normal_distribution(x_save,0.0_kind_real,1.0_kind_real)
call normal_distribution(q_save,0.0_kind_real,1.0_kind_real)
q = q_save
x = 0.0
call convert_x_to_q_ad(fld%geom,q,x)
call convert_x_to_q_tl(fld%geom,x_save,q)
dot1 = sum(x*x_save)
dot2 = sum(q*q_save)
if (abs(dot1-dot2)/abs(dot1)>eps) then
  write(record,*) 'qg_model_propagate_test: convert_x_to_q adjoint test failed:',dot1,dot2,abs(dot1-dot2)/abs(dot1)
  call abor1_ftn(record)
endif

! Test x_to_uv
call normal_distribution(x_save,0.0_kind_real,1.0_kind_real)
call normal_distribution(u_save,0.0_kind_real,1.0_kind_real)
call normal_distribution(v_save,0.0_kind_real,1.0_kind_real)
u = u_save
v = v_save
x = 0.0
call convert_x_to_uv_ad(fld%geom,u,v,x)
call convert_x_to_uv_tl(fld%geom,x_save,u,v)
dot1 = sum(x*x_save)
dot2 = sum(u*u_save)+sum(v*v_save)
if (abs(dot1-dot2)/abs(dot1)>eps) then
  write(record,*) 'qg_model_propagate_test: convert_x_to_uv adjoint test failed:',dot1,dot2,abs(dot1-dot2)/abs(dot1)
  call abor1_ftn(record)
endif

! Test advect_q
call normal_distribution(u_save,0.0_kind_real,1.0_kind_real)
call normal_distribution(v_save,0.0_kind_real,1.0_kind_real)
call normal_distribution(q_save,0.0_kind_real,1.0_kind_real)
call normal_distribution(qnew_save,0.0_kind_real,1.0_kind_real)
qnew = qnew_save
u = 0.0
v = 0.0
q = 0.0
call advect_q_ad(fld%geom,conf%dt,u_traj,v_traj,q_traj,traj%q_north,traj%q_south,qnew,u,v,q)
call advect_q_tl(fld%geom,conf%dt,u_traj,v_traj,q_traj,traj%q_north,traj%q_south,u_save,v_save,q_save,qnew)
dot1 = sum(u*u_save)+sum(v*v_save)+sum(q*q_save)
dot2 = sum(qnew*qnew_save)
if (abs(dot1-dot2)/abs(dot1)>eps) then
    write(record,*) 'qg_model_propagate_test: advect_q adjoint test failed:',dot1,dot2,abs(dot1-dot2)/abs(dot1)
    call abor1_ftn(record)
endif

! Test propagate
call qg_fields_create_from_other(fld1,fld)
call qg_fields_create_from_other(fld2,fld)
call qg_fields_create_from_other(fld1_save,fld)
call qg_fields_create_from_other(fld2_save,fld)
call qg_fields_random(fld1_save)
call qg_fields_random(fld2_save)
call qg_fields_copy(fld1,fld1_save)
call qg_fields_copy(fld2,fld2_save)
call qg_model_propagate_ad(conf,traj,fld2)
call qg_model_propagate_tl(conf,traj,fld1)
call qg_fields_dot_prod(fld1,fld2_save,dot1)
call qg_fields_dot_prod(fld2,fld1_save,dot2)
if (abs(dot1-dot2)/abs(dot1)>eps) then
    write(record,*) 'qg_model_propagate_test:propagate adjoint test failed:',dot1,dot2,abs(dot1-dot2)/abs(dot1)
    call abor1_ftn(record)
endif

end subroutine qg_model_propagate_test
! ------------------------------------------------------------------------------
end module qg_model_mod

