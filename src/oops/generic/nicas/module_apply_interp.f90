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

use module_apply_com, only: alpha_com_AB
use module_namelist, only: nam
use omp_lib
use tools_kinds,only: kind_real
use tools_missing, only: msr
use type_fields, only: fldtype,alphatype
use type_linop, only: apply_linop,apply_linop_ad
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
type(ndatatype),intent(in) :: ndata !< Sampling data
real(kind_real),intent(in) :: alpha(ndata%ns) !< Subgrid variable
real(kind_real),intent(out) :: fld(ndata%nc0,ndata%nl0)  !< Field

! Local variables
integer :: is,il1,ic1,il0
real(kind_real) :: beta(ndata%nc1,ndata%nl1),gamma(ndata%nc1,ndata%nl1),delta(ndata%nc1,ndata%nl0)

!$omp parallel do private(is)
do is=1,ndata%ns
   if (nam%lsqrt) then
      ! Internal normalization
      beta(ndata%is_to_ic2(is),ndata%is_to_il1(is)) = alpha(is)*ndata%norm_sqrt(is)
   else
      ! Copy
      beta(ndata%is_to_ic2(is),ndata%is_to_il1(is)) = alpha(is)
   end if
end do
!$omp end parallel do

! Subsampling horizontal interpolation
!$omp parallel do private(il1)
do il1=1,ndata%nl1
   call apply_linop(ndata%s(il1),beta(1:ndata%nc2(il1),il1),gamma(:,il1))
end do
!$omp end parallel do

! Vertical interpolation
!$omp parallel do private(ic1)
do ic1=1,ndata%nc1
   call apply_linop(ndata%v(ndata%vbot(ic1)),gamma(ic1,:),delta(ic1,:))
end do
!$omp end parallel do

! Horizontal interpolation
!$omp parallel do private(il0)
do il0=1,ndata%nl0
   call apply_linop(ndata%h(min(il0,ndata%nl0i)),delta(:,il0),fld(:,il0))
end do
!$omp end parallel do

! Normalization
fld = fld*ndata%norm

end subroutine interp_global

!----------------------------------------------------------------------
! Subroutine: interp_local
!> Purpose: interpolation, local
!----------------------------------------------------------------------
subroutine interp_local(ndataloc,alpha,fld)

implicit none

! Passed variables
type(ndataloctype),intent(in) :: ndataloc !< Sampling data
real(kind_real),intent(in) :: alpha(ndataloc%nsb) !< Subgrid variable
real(kind_real),intent(out) :: fld(ndataloc%nc0a,ndataloc%nl0)  !< Field

! Local variables
integer :: isb,il1,ic1b,il0
real(kind_real) :: beta(ndataloc%nc1b,ndataloc%nl1),gamma(ndataloc%nc1b,ndataloc%nl1),delta(ndataloc%nc1b,ndataloc%nl0)

!$omp parallel do private(isb)
do isb=1,ndataloc%nsb
   if (nam%lsqrt) then
      ! Internal normalization
      beta(ndataloc%isb_to_ic2b(isb),ndataloc%isb_to_il1(isb)) = alpha(isb)*ndataloc%norm_sqrt(isb)
   else
      ! Copy
      beta(ndataloc%isb_to_ic2b(isb),ndataloc%isb_to_il1(isb)) = alpha(isb)
   end if
end do
!$omp end parallel do

! Subsampling horizontal interpolation
!$omp parallel do private(il1)
do il1=1,ndataloc%nl1
   call apply_linop(ndataloc%s(il1),beta(1:ndataloc%nc2b(il1),il1),gamma(:,il1))
end do
!$omp end parallel do

! Vertical interpolation
!$omp parallel do private(ic1b)
do ic1b=1,ndataloc%nc1b
   call apply_linop(ndataloc%v(ndataloc%vbot(ic1b)),gamma(ic1b,:),delta(ic1b,:))
end do
!$omp end parallel do

! Horizontal interpolation
!$omp parallel do private(il0)
do il0=1,ndataloc%nl0
   call apply_linop(ndataloc%h(min(il0,ndataloc%nl0i)),delta(:,il0),fld(:,il0))
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
type(ndatatype),intent(in) :: ndata    !< Sampling data
real(kind_real),intent(in) :: fld(ndata%nc0,ndata%nl0)  !< Field
real(kind_real),intent(out) :: alpha(ndata%ns) !< Subgrid variable

! Local variables
integer :: is,il1,ic1,il0
real(kind_real) :: beta(ndata%nc1,ndata%nl1),gamma(ndata%nc1,ndata%nl1),delta(ndata%nc1,ndata%nl0)
real(kind_real) :: fld_tmp(ndata%nc0,ndata%nl0)

! Normalization
fld_tmp = fld*ndata%norm

! Horizontal interpolation
!$omp parallel do private(il0)
do il0=1,ndata%nl0
   call apply_linop_ad(ndata%h(min(il0,ndata%nl0i)),fld_tmp(:,il0),delta(:,il0))
end do
!$omp end parallel do

! Vertical interpolation
!$omp parallel do private(ic1)
do ic1=1,ndata%nc1
   call apply_linop_ad(ndata%v(ndata%vbot(ic1)),delta(ic1,:),gamma(ic1,:))
end do
!$omp end parallel do

! Subsampling horizontal interpolation
!$omp parallel do private(il1)
do il1=1,ndata%nl1
   call apply_linop_ad(ndata%s(il1),gamma(:,il1),beta(1:ndata%nc2(il1),il1))
end do
!$omp end parallel do

!$omp parallel do private(is)
do is=1,ndata%ns
   if (nam%lsqrt) then
      ! Internal normalization
      alpha(is) = beta(ndata%is_to_ic2(is),ndata%is_to_il1(is))*ndata%norm_sqrt(is)
   else
      ! Copy
      alpha(is) = beta(ndata%is_to_ic2(is),ndata%is_to_il1(is))
   end if
end do
!$omp end parallel do

end subroutine interp_ad_global

!----------------------------------------------------------------------
! Subroutine: interp_ad_local
!> Purpose: interpolation adjoint, local
!----------------------------------------------------------------------
subroutine interp_ad_local(ndataloc,fld,alpha)

implicit none

! Passed variables
type(ndataloctype),intent(in) :: ndataloc  !< Sampling data
real(kind_real),intent(in) :: fld(ndataloc%nc0a,ndataloc%nl0)  !< Field
real(kind_real),intent(out) :: alpha(ndataloc%nsb) !< Subgrid variable

! Local variables
integer :: isb,il1,ic1b,il0
real(kind_real) :: beta(ndataloc%nc1b,ndataloc%nl1),gamma(ndataloc%nc1b,ndataloc%nl1),delta(ndataloc%nc1b,ndataloc%nl0)
real(kind_real) :: fld_tmp(ndataloc%nc0a,ndataloc%nl0)

! Normalization
fld_tmp = fld*ndataloc%norm

! Horizontal interpolation
!$omp parallel do private(il0)
do il0=1,ndataloc%nl0
   call apply_linop_ad(ndataloc%h(min(il0,ndataloc%nl0i)),fld_tmp(:,il0),delta(:,il0))
end do
!$omp end parallel do

! Vertical interpolation
!$omp parallel do private(ic1b)
do ic1b=1,ndataloc%nc1b
   call apply_linop_ad(ndataloc%v(ndataloc%vbot(ic1b)),delta(ic1b,:),gamma(ic1b,:))
end do
!$omp end parallel do

! Subsampling horizontal interpolation
!$omp parallel do private(il1)
do il1=1,ndataloc%nl1
   call apply_linop_ad(ndataloc%s(il1),gamma(:,il1),beta(1:ndataloc%nc2b(il1),il1))
end do
!$omp end parallel do

!$omp parallel do private(isb)
do isb=1,ndataloc%nsb
   if (nam%lsqrt) then
      ! Internal normalization
      alpha(isb) = beta(ndataloc%isb_to_ic2b(isb),ndataloc%isb_to_il1(isb))*ndataloc%norm_sqrt(isb)
   else
      ! Copy
      alpha(isb) = beta(ndataloc%isb_to_ic2b(isb),ndataloc%isb_to_il1(isb))
   end if
end do
!$omp end parallel do

end subroutine interp_ad_local

end module module_apply_interp
