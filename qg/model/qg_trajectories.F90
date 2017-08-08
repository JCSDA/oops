! (C) Copyright 2009-2016 ECMWF.
! 
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
! In applying this licence, ECMWF does not waive the privileges and immunities 
! granted to it by virtue of its status as an intergovernmental organisation nor
! does it submit to any jurisdiction.

!> Handle the trajectory for the QG model

module qg_trajectories

use kinds
implicit none
private

public :: qg_trajectory, set_traj, get_traj, delete_traj
public :: qg_traj_registry

! ------------------------------------------------------------------------------

!> Fortran derived type to hold the QG model trajectory
type :: qg_trajectory
  real(kind=kind_real), allocatable :: x(:,:,:)
  real(kind=kind_real), allocatable :: xn(:)
  real(kind=kind_real), allocatable :: xs(:)
  real(kind=kind_real), allocatable :: qn(:,:)
  real(kind=kind_real), allocatable :: qs(:,:)
end type qg_trajectory

#define LISTED_TYPE qg_trajectory

!> Linked list interface - defines registry_t type
#include "util/linkedList_i.f"

!> Global registry
type(registry_t) :: qg_traj_registry

! ------------------------------------------------------------------------------
contains
! ------------------------------------------------------------------------------
!> Linked list implementation
#include "util/linkedList_c.f"

! ------------------------------------------------------------------------------

subroutine set_traj(self,kx,ky,px,pxn,pxs,pqn,pqs)
implicit none
type(qg_trajectory), intent(inout) :: self
integer, intent(in) :: kx,ky
real(kind=kind_real), intent(in)    :: px(kx,ky,2)
real(kind=kind_real), intent(in)    :: pxn(2),pxs(2)
real(kind=kind_real), intent(in)    :: pqn(kx,2),pqs(kx,2)

allocate(self%x(kx,ky,2))
allocate(self%xn(2),self%xs(2))
allocate(self%qn(kx,2),self%qs(kx,2))

self%x(:,:,:)=px(:,:,:)
self%xn(:)   =pxn(:)
self%xs(:)   =pxs(:)
self%qn(:,:) =pqn(:,:)
self%qs(:,:) =pqs(:,:)

return
end subroutine set_traj

! ------------------------------------------------------------------------------

subroutine get_traj(self,kx,ky,px,pxn,pxs,pqn,pqs)
implicit none
type(qg_trajectory), intent(in) :: self
integer, intent(in) :: kx,ky
real(kind=kind_real), intent(inout) :: px(kx,ky,2)
real(kind=kind_real), intent(inout) :: pxn(2),pxs(2)
real(kind=kind_real), intent(inout) :: pqn(kx,2),pqs(kx,2)

px(:,:,:)=self%x(:,:,:)
pxn(:)   =self%xn(:)
pxs(:)   =self%xs(:)
pqn(:,:) =self%qn(:,:)
pqs(:,:) =self%qs(:,:)

return
end subroutine get_traj

! ------------------------------------------------------------------------------

subroutine delete_traj(self)
implicit none
type(qg_trajectory), intent(inout) :: self

deallocate(self%x)
deallocate(self%xn)
deallocate(self%xs)
deallocate(self%qn)
deallocate(self%qs)

return
end subroutine delete_traj

! ------------------------------------------------------------------------------

subroutine c_minmax_traj(c_key_self, pminmax) bind(c,name='qg_traj_minmaxrms_f90')
use iso_c_binding
implicit none
integer(c_int), intent(in) :: c_key_self
real(c_double), intent(inout) :: pminmax(3,5)
type(qg_trajectory), pointer :: self
real(kind=kind_real) :: zz

call qg_traj_registry%get(c_key_self,self)

zz=real(size(self%x),kind_real)
pminmax(1,1)=minval(self%x(:,:,:))
pminmax(2,1)=maxval(self%x(:,:,:))
pminmax(3,1)=sqrt(sum(self%x(:,:,:)**2)/zz)

zz=2.0_kind_real
pminmax(1,2)=minval(self%xn(:))
pminmax(2,2)=maxval(self%xn(:))
pminmax(3,2)=sqrt(sum(self%xn(:)**2)/zz)

pminmax(1,3)=minval(self%xs(:))
pminmax(2,3)=maxval(self%xs(:))
pminmax(3,3)=sqrt(sum(self%xs(:)**2)/zz)

zz=real(size(self%qn),kind_real)
pminmax(1,4)=minval(self%qn(:,:))
pminmax(2,4)=maxval(self%qn(:,:))
pminmax(3,4)=sqrt(sum(self%qn(:,:)**2)/zz)

pminmax(1,5)=minval(self%qs(:,:))
pminmax(2,5)=maxval(self%qs(:,:))
pminmax(3,5)=sqrt(sum(self%qs(:,:)**2)/zz)

end subroutine c_minmax_traj

! ------------------------------------------------------------------------------

end module qg_trajectories
