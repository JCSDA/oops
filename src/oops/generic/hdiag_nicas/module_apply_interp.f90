!----------------------------------------------------------------------
! Module: module_apply_interp
!> Purpose: interpolation routines
!> <br>
!> Author: Benjamin Menetrier
!> <br>
!> Licensing: this code is distributed under the CeCILL-C license
!> <br>
!> Copyright © 2017 METEO-FRANCE
!----------------------------------------------------------------------
module module_apply_interp

use omp_lib
use tools_kinds,only: kind_real
use tools_missing, only: msr
use type_geom, only: geomtype
use type_linop, only: apply_linop,apply_linop_ad
use type_mpl, only: mpl
use type_nam, only: namtype
use type_ndata, only: ndatatype,ndataloctype

implicit none

interface interp
   module procedure interp_global
   module procedure interp_local
end interface

interface interp_ad
   module procedure interp_ad_global
   module procedure interp_ad_local
end interface

private
public :: interp,interp_ad

contains

!----------------------------------------------------------------------
! Subroutine: interp_global
!> Purpose: interpolation, global
!----------------------------------------------------------------------
subroutine interp_global(ndata,alpha,fld)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata                               !< NICAS data
real(kind_real),intent(in) :: alpha(ndata%ns)                     !< Subgrid field
real(kind_real),intent(out) :: fld(ndata%geom%nc0,ndata%geom%nl0) !< Field

! Local variables
integer :: is,il1,ic1,il0
real(kind_real) :: beta(ndata%nc1,ndata%nl1),gamma(ndata%nc1,ndata%nl1),delta(ndata%nc1,ndata%geom%nl0)
real(kind_real) :: gammaT(ndata%nl1,ndata%nc1),deltaT(ndata%geom%nl0,ndata%nc1)

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

!$omp parallel do schedule(static) private(is)
do is=1,ndata%ns
   beta(ndata%is_to_ic2(is),ndata%is_to_il1(is)) = alpha(is)
   if (nam%lsqrt) beta(ndata%is_to_ic2(is),ndata%is_to_il1(is)) = &
 & beta(ndata%is_to_ic2(is),ndata%is_to_il1(is))*ndata%norm_sqrt(is)
end do
!$omp end parallel do

! Subsampling horizontal interpolation
!$omp parallel do schedule(static) private(il1)
do il1=1,ndata%nl1
   call apply_linop(ndata%s(il1),beta(1:ndata%nc2(il1),il1),gamma(:,il1))
end do
!$omp end parallel do

! Transpose data
gammaT = transpose(gamma)

! Vertical interpolation
!$omp parallel do schedule(static) private(ic1)
do ic1=1,ndata%nc1
   call apply_linop(ndata%v,gammaT(:,ic1),deltaT(:,ic1))
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

! End associate
end associate

end subroutine interp_global

!----------------------------------------------------------------------
! Subroutine: interp_local
!> Purpose: interpolation, local
!----------------------------------------------------------------------
subroutine interp_local(nam,geom,ndataloc,alpha,fld)

implicit none

! Passed variables
type(namtype),intent(in) :: nam                        !< Namelist
type(geomtype),intent(in) :: geom                      !< Geometry
type(ndataloctype),intent(in) :: ndataloc              !< NICAS data, local
real(kind_real),intent(in) :: alpha(ndataloc%nsb)      !< Subgrid field, local
real(kind_real),intent(out) :: fld(geom%nc0a,geom%nl0) !< Field, local

! Local variables
integer :: isb,il1,ic1b,il0
real(kind_real) :: beta(ndataloc%nc1b,ndataloc%nl1),gamma(ndataloc%nc1b,ndataloc%nl1),delta(ndataloc%nc1b,geom%nl0)
real(kind_real) :: gammaT(ndataloc%nl1,ndataloc%nc1b),deltaT(geom%nl0,ndataloc%nc1b)

!$omp parallel do schedule(static) private(isb)
do isb=1,ndataloc%nsb
   beta(ndataloc%isb_to_ic2b(isb),ndataloc%isb_to_il1(isb)) = alpha(isb)
   if (nam%lsqrt) beta(ndataloc%isb_to_ic2b(isb),ndataloc%isb_to_il1(isb)) = &
 & beta(ndataloc%isb_to_ic2b(isb),ndataloc%isb_to_il1(isb))*ndataloc%norm_sqrt(isb)
end do
!$omp end parallel do

! Subsampling horizontal interpolation
!$omp parallel do schedule(static) private(il1)
do il1=1,ndataloc%nl1
   call apply_linop(ndataloc%s(il1),beta(1:ndataloc%nc2b(il1),il1),gamma(:,il1))
end do
!$omp end parallel do

! Transpose data
gammaT = transpose(gamma)

! Vertical interpolation
!$omp parallel do schedule(static) private(ic1b)
do ic1b=1,ndataloc%nc1b
   call apply_linop(ndataloc%v,gammaT(:,ic1b),deltaT(:,ic1b))
end do
!$omp end parallel do

! Transpose data
delta = transpose(deltaT)

! Horizontal interpolation
!$omp parallel do schedule(static) private(il0)
do il0=1,geom%nl0
   call apply_linop(ndataloc%h(min(il0,geom%nl0i)),delta(:,il0),fld(:,il0))
end do
!$omp end parallel do

! Normalization
fld = fld*ndataloc%norm

end subroutine interp_local

!----------------------------------------------------------------------
! Subroutine: interp_ad_global
!> Purpose: interpolation adjoint, global
!----------------------------------------------------------------------
subroutine interp_ad_global(ndata,fld,alpha)

implicit none

! Passed variables
type(ndatatype),intent(in) :: ndata                              !< NICAS data
real(kind_real),intent(in) :: fld(ndata%geom%nc0,ndata%geom%nl0) !< Field
real(kind_real),intent(out) :: alpha(ndata%ns)                   !< Subgrid field

! Local variables
integer :: is,il1,ic1,il0
real(kind_real) :: beta(ndata%nc1,ndata%nl1),gamma(ndata%nc1,ndata%nl1),delta(ndata%nc1,ndata%geom%nl0)
real(kind_real) :: gammaT(ndata%nl1,ndata%nc1),deltaT(ndata%geom%nl0,ndata%nc1)
real(kind_real) :: fld_tmp(ndata%geom%nc0,ndata%geom%nl0)

! Associate
associate(nam=>ndata%nam,geom=>ndata%geom)

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
!$omp parallel do schedule(static) private(ic1)
do ic1=1,ndata%nc1
   call apply_linop_ad(ndata%v,deltaT(:,ic1),gammaT(:,ic1))
end do
!$omp end parallel do

! Transpose data
gamma = transpose(gammaT)

! Subsampling horizontal interpolation
!$omp parallel do schedule(static) private(il1)
do il1=1,ndata%nl1
   call apply_linop_ad(ndata%s(il1),gamma(:,il1),beta(1:ndata%nc2(il1),il1))
end do
!$omp end parallel do

!$omp parallel do schedule(static) private(is)
do is=1,ndata%ns
   alpha(is) = beta(ndata%is_to_ic2(is),ndata%is_to_il1(is))
   if (nam%lsqrt) alpha(is) = alpha(is)*ndata%norm_sqrt(is)
end do
!$omp end parallel do

! End associate
end associate

end subroutine interp_ad_global

!----------------------------------------------------------------------
! Subroutine: interp_ad_local
!> Purpose: interpolation adjoint, local
!----------------------------------------------------------------------
subroutine interp_ad_local(nam,geom,ndataloc,fld,alpha)

implicit none

! Passed variables
type(namtype),intent(in) :: nam                       !< Namelist
type(geomtype),intent(in) :: geom                     !< Geometry
type(ndataloctype),intent(in) :: ndataloc             !< NICAS data, local
real(kind_real),intent(in) :: fld(geom%nc0a,geom%nl0) !< Field, local
real(kind_real),intent(out) :: alpha(ndataloc%nsb)    !< Subgrid field, local

! Local variables
integer :: isb,il1,ic1b,il0
real(kind_real) :: beta(ndataloc%nc1b,ndataloc%nl1),gamma(ndataloc%nc1b,ndataloc%nl1),delta(ndataloc%nc1b,geom%nl0)
real(kind_real) :: gammaT(ndataloc%nl1,ndataloc%nc1b),deltaT(geom%nl0,ndataloc%nc1b)
real(kind_real) :: fld_tmp(geom%nc0a,geom%nl0)

! Normalization
fld_tmp = fld*ndataloc%norm

! Horizontal interpolation
!$omp parallel do schedule(static) private(il0)
do il0=1,geom%nl0
   call apply_linop_ad(ndataloc%h(min(il0,geom%nl0i)),fld_tmp(:,il0),delta(:,il0))
end do
!$omp end parallel do

! Transpose data
deltaT = transpose(delta)

! Vertical interpolation
!$omp parallel do schedule(static) private(ic1b)
do ic1b=1,ndataloc%nc1b
   call apply_linop_ad(ndataloc%v,deltaT(:,ic1b),gammaT(:,ic1b))
end do
!$omp end parallel do

! Transpose data
gamma = transpose(gammaT)

! Subsampling horizontal interpolation
!$omp parallel do schedule(static) private(il1)
do il1=1,ndataloc%nl1
   call apply_linop_ad(ndataloc%s(il1),gamma(:,il1),beta(1:ndataloc%nc2b(il1),il1))
end do
!$omp end parallel do

!$omp parallel do schedule(static) private(isb)
do isb=1,ndataloc%nsb
   alpha(isb) = beta(ndataloc%isb_to_ic2b(isb),ndataloc%isb_to_il1(isb))
   if (nam%lsqrt) alpha(isb) = alpha(isb)*ndataloc%norm_sqrt(isb)
end do
!$omp end parallel do

end subroutine interp_ad_local

end module module_apply_interp
