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
use tools_kinds,only: kind_real
use tools_missing, only: msr
use type_geom, only: geomtype
use type_linop, only: apply_linop,apply_linop_ad
use type_mpl, only: mpl
use type_nam, only: namtype
use type_ndata, only: ndatatype

implicit none

private
public :: apply_interp,apply_interp_ad

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
integer :: isb,il1,ic1b,il0
real(kind_real) :: beta(ndata%nc1b,ndata%nl1),gamma(ndata%nc1b,ndata%nl1),delta(ndata%nc1b,geom%nl0)
real(kind_real) :: gammaT(ndata%nl1,ndata%nc1b),deltaT(geom%nl0,ndata%nc1b)

!$omp parallel do schedule(static) private(isb)
do isb=1,ndata%nsb
   beta(ndata%sb_to_c2b(isb),ndata%sb_to_l1(isb)) = alpha(isb)
end do
!$omp end parallel do

! Subsampling horizontal interpolation
!$omp parallel do schedule(static) private(il1)
do il1=1,ndata%nl1
   call apply_linop(ndata%s(il1),beta(1:ndata%nc2b(il1),il1),gamma(:,il1))
end do
!$omp end parallel do

! Transpose data
gammaT = transpose(gamma)

! Vertical interpolation
!$omp parallel do schedule(static) private(ic1b)
do ic1b=1,ndata%nc1b
   call apply_linop(ndata%v,gammaT(:,ic1b),deltaT(:,ic1b))
end do
!$omp end parallel do

! Transpose data
delta = transpose(deltaT)

! Horizontal interpolation
!$omp parallel do schedule(static) private(il0)
do il0=1,geom%nl0
   call apply_linop(ndata%h(min(il0,geom%nl0i)),delta(:,il0),fld(:,il0))
end do
!$omp end parallel do

! Normalization
fld = fld*ndata%norm

end subroutine apply_interp

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
integer :: isb,il1,ic1b,il0
real(kind_real) :: beta(ndata%nc1b,ndata%nl1),gamma(ndata%nc1b,ndata%nl1),delta(ndata%nc1b,geom%nl0)
real(kind_real) :: gammaT(ndata%nl1,ndata%nc1b),deltaT(geom%nl0,ndata%nc1b)
real(kind_real) :: fld_tmp(geom%nc0a,geom%nl0)

! Normalization
fld_tmp = fld*ndata%norm

! Horizontal interpolation
!$omp parallel do schedule(static) private(il0)
do il0=1,geom%nl0
   call apply_linop_ad(ndata%h(min(il0,geom%nl0i)),fld_tmp(:,il0),delta(:,il0))
end do
!$omp end parallel do

! Transpose data
deltaT = transpose(delta)

! Vertical interpolation
!$omp parallel do schedule(static) private(ic1b)
do ic1b=1,ndata%nc1b
   call apply_linop_ad(ndata%v,deltaT(:,ic1b),gammaT(:,ic1b))
end do
!$omp end parallel do

! Transpose data
gamma = transpose(gammaT)

! Subsampling horizontal interpolation
!$omp parallel do schedule(static) private(il1)
do il1=1,ndata%nl1
   call apply_linop_ad(ndata%s(il1),gamma(:,il1),beta(1:ndata%nc2b(il1),il1))
end do
!$omp end parallel do

!$omp parallel do schedule(static) private(isb)
do isb=1,ndata%nsb
   alpha(isb) = beta(ndata%sb_to_c2b(isb),ndata%sb_to_l1(isb))
end do
!$omp end parallel do

end subroutine apply_interp_ad

end module nicas_apply_interp
