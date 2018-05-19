!----------------------------------------------------------------------
! Module: nicas_apply_interp
!> Purpose: interpolation routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module nicas_apply_interp

use omp_lib
use tools_kinds, only: kind_real
use tools_display, only: msgerror
use tools_missing, only: msr
use type_geom, only: geomtype
use type_linop, only: apply_linop,apply_linop_ad
use type_mpl, only: mpl
use type_nam, only: namtype
use type_ndata, only: ndatatype
use yomhook, only: lhook,dr_hook

implicit none

private
public :: apply_interp,apply_interp_s,apply_interp_v,apply_interp_h
public :: apply_interp_ad,apply_interp_h_ad,apply_interp_v_ad,apply_interp_s_ad

contains

!----------------------------------------------------------------------
! Subroutine: apply_interp
!> Purpose: apply interpolation
!----------------------------------------------------------------------
subroutine apply_interp(geom,ndata,alpha,fld)

implicit none

! Passed variables
type(geomtype),intent(in) :: geom                      !< Geometry
type(ndatatype),intent(in) :: ndata                    !< NICAS data
real(kind_real),intent(in) :: alpha(ndata%nsb)         !< Subgrid field
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0) !< Field

! Local variables
real(kind_real) :: gamma(ndata%nc1b,ndata%nl1),delta(ndata%nc1b,geom%nl0)
real(kind_real) :: zhook_handle

if (lhook) call dr_hook('apply_interp',0,zhook_handle)

! Subsampling horizontal interpolation
call apply_interp_s(ndata,alpha,gamma)

! Vertical interpolation
call apply_interp_v(geom,ndata,gamma,delta)

! Horizontal interpolation
call apply_interp_h(geom,ndata,delta,fld)

! Normalization
fld = fld*ndata%norm

if (lhook) call dr_hook('apply_interp',1,zhook_handle)

end subroutine apply_interp

!----------------------------------------------------------------------
! Subroutine: apply_interp_s
!> Purpose: apply subsampling interpolation
!----------------------------------------------------------------------
subroutine apply_interp_s(ndata,alpha,gamma)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata                        !< NICAS data
real(kind_real),intent(in) :: alpha(ndata%nsb)             !< Subgrid field
real(kind_real),intent(out) :: gamma(ndata%nc1b,ndata%nl1) !< Subset Sc1 field, limited levels

! Local variables
integer :: isb,il1
real(kind_real) :: beta(ndata%nc1b,ndata%nl1)

! Initialization
beta = 0.0

! Copy
!$omp parallel do schedule(static) private(isb)
do isb=1,ndata%nsb
   beta(ndata%sb_to_c1b(isb),ndata%sb_to_l1(isb)) = alpha(isb)
end do
!$omp end parallel do

! Subsampling horizontal interpolation
!$omp parallel do schedule(static) private(il1)
do il1=1,ndata%nl1
   call apply_linop(ndata%s(il1),beta(:,il1),gamma(:,il1))
end do
!$omp end parallel do

end subroutine apply_interp_s

!----------------------------------------------------------------------
! Subroutine: apply_interp_v
!> Purpose: apply vertical interpolation
!----------------------------------------------------------------------
subroutine apply_interp_v(geom,ndata,gamma,delta)

implicit none

! Passed variables
type(geomtype),intent(in) :: geom                         !< Geometry
type(ndatatype),intent(in) :: ndata                       !< NICAS data
real(kind_real),intent(in) :: gamma(ndata%nc1b,ndata%nl1) !< Subset Sc1 field, limited levels
real(kind_real),intent(out) :: delta(ndata%nc1b,geom%nl0) !< Subset Sc1 field, full levels

! Local variables
integer :: ic1b
real(kind_real),allocatable :: gamma_tmp(:),delta_tmp(:)

! Vertical interpolation
!$omp parallel do schedule(static) private(ic1b,gamma_tmp,delta_tmp)
do ic1b=1,ndata%nc1b
   ! Allocation
   allocate(gamma_tmp(ndata%nl1))
   allocate(delta_tmp(geom%nl0))

   ! Copy data
   gamma_tmp = gamma(ic1b,:)

   ! Apply interpolation
   call apply_linop(ndata%v,gamma_tmp,delta_tmp)

   ! Copy data
   delta(ic1b,:) = delta_tmp

   ! Release memory
   deallocate(gamma_tmp)
   deallocate(delta_tmp)
end do
!$omp end parallel do

end subroutine apply_interp_v

!----------------------------------------------------------------------
! Subroutine: apply_interp_h
!> Purpose: apply horizontal interpolation
!----------------------------------------------------------------------
subroutine apply_interp_h(geom,ndata,delta,fld)

implicit none

! Passed variables
type(geomtype),intent(in) :: geom                        !< Geometry
type(ndatatype),intent(in) :: ndata                      !< NICAS data
real(kind_real),intent(in) :: delta(ndata%nc1b,geom%nl0) !< Subset Sc1 field, full levels
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0)   !< Field

! Local variables
integer :: il0

! Horizontal interpolation
!$omp parallel do schedule(static) private(il0)
do il0=1,geom%nl0
   call apply_linop(ndata%h(min(il0,geom%nl0i)),delta(:,il0),fld(:,il0))
end do
!$omp end parallel do

end subroutine apply_interp_h

!----------------------------------------------------------------------
! Subroutine: apply_interp_ad
!> Purpose: apply interpolation adjoint
!----------------------------------------------------------------------
subroutine apply_interp_ad(geom,ndata,fld,alpha)

implicit none

! Passed variables
type(geomtype),intent(in) :: geom                     !< Geometry
type(ndatatype),intent(in) :: ndata                   !< NICAS data
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0) !< Field
real(kind_real),intent(out) :: alpha(ndata%nsb)       !< Subgrid field

! Local variables
real(kind_real) :: gamma(ndata%nc1b,ndata%nl1),delta(ndata%nc1b,geom%nl0)
real(kind_real) :: fld_tmp(geom%nc0a,geom%nl0)
real(kind_real) :: zhook_handle

if (lhook) call dr_hook('apply_interp_ad',0,zhook_handle)

! Normalization
fld_tmp = fld*ndata%norm

! Horizontal interpolation
call apply_interp_h_ad(geom,ndata,fld_tmp,delta)

! Vertical interpolation
call apply_interp_v_ad(geom,ndata,delta,gamma)

! Subsampling horizontal interpolation
call apply_interp_s_ad(ndata,gamma,alpha)

if (lhook) call dr_hook('apply_interp_ad',1,zhook_handle)

end subroutine apply_interp_ad

!----------------------------------------------------------------------
! Subroutine: apply_interp_h_ad
!> Purpose: apply horizontal interpolation adjoint
!----------------------------------------------------------------------
subroutine apply_interp_h_ad(geom,ndata,fld,delta)

implicit none

! Passed variables
type(geomtype),intent(in) :: geom                         !< Geometry
type(ndatatype),intent(in) :: ndata                       !< NICAS data
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0)     !< Field
real(kind_real),intent(out) :: delta(ndata%nc1b,geom%nl0) !< Subset Sc1 field, full levels

! Local variables
integer :: il0

!$omp parallel do schedule(static) private(il0)
do il0=1,geom%nl0
   call apply_linop_ad(ndata%h(min(il0,geom%nl0i)),fld(:,il0),delta(:,il0))
end do
!$omp end parallel do

end subroutine apply_interp_h_ad

!----------------------------------------------------------------------
! Subroutine: apply_interp_v_ad
!> Purpose: apply vertical interpolation adjoint
!----------------------------------------------------------------------
subroutine apply_interp_v_ad(geom,ndata,delta,gamma)

implicit none

! Passed variables
type(geomtype),intent(in) :: geom                          !< Geometry
type(ndatatype),intent(in) :: ndata                        !< NICAS data
real(kind_real),intent(in) :: delta(ndata%nc1b,geom%nl0)   !< Subset Sc1 field, full levels
real(kind_real),intent(out) :: gamma(ndata%nc1b,ndata%nl1) !< Subset Sc1 field, limited levels

! Local variables
integer :: ic1b
real(kind_real),allocatable :: gamma_tmp(:),delta_tmp(:)

! Vertical interpolation
!$omp parallel do schedule(static) private(ic1b,gamma_tmp,delta_tmp)
do ic1b=1,ndata%nc1b
   ! Allocation
   allocate(gamma_tmp(ndata%nl1))
   allocate(delta_tmp(geom%nl0))

   ! Copy data
   delta_tmp = delta(ic1b,:)

   ! Apply interpolation
   call apply_linop_ad(ndata%v,delta_tmp,gamma_tmp)

   ! Copy data
   gamma(ic1b,:) = gamma_tmp

   ! Release memory
   deallocate(gamma_tmp)
   deallocate(delta_tmp)
end do
!$omp end parallel do

end subroutine apply_interp_v_ad

!----------------------------------------------------------------------
! Subroutine: apply_interp_s_ad
!> Purpose: apply subsampling interpolation adjoint
!----------------------------------------------------------------------
subroutine apply_interp_s_ad(ndata,gamma,alpha)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata                       !< NICAS data
real(kind_real),intent(in) :: gamma(ndata%nc1b,ndata%nl1) !< Subset Sc1 field, limited levels
real(kind_real),intent(out) :: alpha(ndata%nsb)           !< Subgrid field

! Local variables
integer :: il1,isb
real(kind_real) :: beta(ndata%nc1b,ndata%nl1)

! Subsampling horizontal interpolation
!$omp parallel do schedule(static) private(il1)
do il1=1,ndata%nl1
   call apply_linop_ad(ndata%s(il1),gamma(:,il1),beta(:,il1))
end do
!$omp end parallel do

! Copy
!$omp parallel do schedule(static) private(isb)
do isb=1,ndata%nsb
   alpha(isb) = beta(ndata%sb_to_c1b(isb),ndata%sb_to_l1(isb))
end do
!$omp end parallel do

end subroutine apply_interp_s_ad

end module nicas_apply_interp
